//! Constructor
template <unsigned Tdim>
mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::MPMSemiImplicitTwoPhaseTwoPoint(
    const std::shared_ptr<IO>& io)
    : mpm::MPMBase<Tdim>(io) {
  //! Logger
  console_ = spdlog::get("MPMSemiImplicitTwoPhaseTwoPoint");
}

//! MPM semi-implicit two-phase two-point solver
template <unsigned Tdim>
bool mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::solve() {
  bool status = true;

  console_->info("MPM analysis type {}", io_->analysis_type());

  // Initialise MPI rank and size
  int mpi_rank = 0;
  int mpi_size = 1;

#ifdef USE_MPI
  // Get MPI rank
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  // This solver consider solid-liquid variables
  const unsigned solid = mpm::ParticlePhase::Solid;
  const unsigned liquid = mpm::ParticlePhase::Liquid;

  // Test if checkpoint resume is needed
  bool resume = false;
  if (analysis_.find("resume") != analysis_.end())
    resume = analysis_["resume"]["resume"].template get<bool>();

  // Pressure smoothing
  if (analysis_.find("pressure_smoothing") != analysis_.end())
    pressure_smoothing_ = analysis_["pressure_smoothing"].template get<bool>();

  // Projection method parameter (beta)
  if (analysis_.find("semi_implicit") != analysis_.end())
    beta_ = analysis_["semi_implicit"]["beta"].template get<double>();

  // Initialise material
  this->initialise_materials();

  // Initialise mesh
  this->initialise_mesh();

  // Initialise particles
  this->initialise_particles();

  // Initialise loading conditions
  this->initialise_loads();

  // Initialise matrix
  bool matrix_status = this->initialise_matrix();
  if (!matrix_status) {
    status = false;
    throw std::runtime_error("Initialisation of matrix failed");
  }

  // Assign porosity to the solid particles
  mesh_->iterate_over_particles_predicate(
      std::bind(&mpm::ParticleBase<Tdim>::assign_porosity,
                std::placeholders::_1),
      std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                solid));

  // Assign permeability to the solid particles
  mesh_->iterate_over_particles_predicate(
      std::bind(&mpm::ParticleBase<Tdim>::assign_permeability,
                std::placeholders::_1),
      std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                solid));

  // Compute nodal volume from gauss quadrature
  this->compute_nodes_gauss_volume();

#pragma omp parallel sections
  {
    // Spawn a task for initialising nodes and cells
#pragma omp section
    {
      mesh_->iterate_over_cells(
          std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
    }
    // Spawn a task for particles
#pragma omp section
    {
      // Iterate over each particle to compute shapefn
      mesh_->iterate_over_particles(std::bind(
          &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));
    }
  }  // Wait to complete

  // Map porosity from solid to fluid particles
  this->compute_fluid_particle_porosity();

  // Compute mass for each phase considering phase-wise volume fraction
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::compute_mass, std::placeholders::_1));

  // Assign beta to each particle
  // TODO: currently only consider full projection
  mesh_->iterate_over_particles(
      std::bind(&mpm::ParticleBase<Tdim>::assign_projection_parameter,
                std::placeholders::_1, beta_));

  // Check point resume
  // TODO: ADD POROSITY IN RESUME FOR SOLID
  if (resume) this->checkpoint_resume();

  // Domain decompose
  bool initial_step = (resume == true) ? false : true;
  this->mpi_domain_decompose(initial_step);

  auto solver_begin = std::chrono::steady_clock::now();
  // Main loop
  for (; step_ < nsteps_; ++step_) {
    if (mpi_rank == 0) console_->info("Step: {} of {}.\n", step_, nsteps_);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    // Run load balancer at a specified frequency
    if (step_ % nload_balance_steps_ == 0 && step_ != 0)
      this->mpi_domain_decompose(false);
#endif
#endif

#pragma omp parallel sections
    {
      // Spawn a task for initialising nodes and cells
#pragma omp section
      {
        // Initialise nodes
        mesh_->iterate_over_nodes(
            std::bind(&mpm::NodeBase<Tdim>::initialise, std::placeholders::_1));

        mesh_->iterate_over_cells(
            std::bind(&mpm::Cell<Tdim>::activate_nodes, std::placeholders::_1));
      }
      // Spawn a task for particles
#pragma omp section
      {
        // Iterate over each particle to compute shapefn
        mesh_->iterate_over_particles(std::bind(
            &mpm::ParticleBase<Tdim>::compute_shapefn, std::placeholders::_1));

        mesh_->iterate_over_cells(
            std::bind(&mpm::Cell<Tdim>::activate_nodes_phase_status,
                      std::placeholders::_1, solid));

        mesh_->iterate_over_cells(
            std::bind(&mpm::Cell<Tdim>::activate_nodes_phase_status,
                      std::placeholders::_1, liquid));
      }
    }  // Wait to complete

    // Assign mass and momentum to nodes
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_mass_momentum_to_nodes,
                  std::placeholders::_1));

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce nodal mass for solid phase
      mesh_->template nodal_halo_exchange<double, 1>(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, solid),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, solid, std::placeholders::_2));
      // MPI all reduce nodal momentum for solid phase
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    solid),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, solid,
                    std::placeholders::_2));

      // MPI all reduce nodal mass for liquid phase
      mesh_->template nodal_halo_exchange<double, 1>(
          std::bind(&mpm::NodeBase<Tdim>::mass, std::placeholders::_1, liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_mass, std::placeholders::_1,
                    false, liquid, std::placeholders::_2));
      // MPI all reduce nodal momentum for liquid phase
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::momentum, std::placeholders::_1,
                    liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_momentum,
                    std::placeholders::_1, false, liquid,
                    std::placeholders::_2));
    }
#endif

    // Compute free surface cells, nodes, and particles
    // TODO: NEED TO CONSIDER POROSITY EFFECT
    mesh_->compute_free_surface(free_surface_detection_, volume_tolerance_);

    // Spawn a task for initializing pressure at free surface
    // FIXME: Only apply this to free surface of fluid particle, otherwise error
    // #pragma omp parallel sections
    //     {
    // #pragma omp section
    //       {
    //         // Assign initial pressure for all free-surface particle
    //         // NOTE: mpm::ParticlePhase::SinglePhase is used here instead of
    //         liquid,
    //         // since that indicate the state variable phase index
    //         mesh_->iterate_over_particles_predicate(
    //             std::bind(&mpm::ParticleBase<Tdim>::assign_pressure,
    //                       std::placeholders::_1, 0.0,
    //                       mpm::ParticlePhase::SinglePhase),
    //             std::bind(&mpm::ParticleBase<Tdim>::free_surface,
    //                       std::placeholders::_1));
    //       }
    //     }  // Wait to complete

    // Compute nodal velocity at the begining of time step
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::compute_velocity,
                  std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Iterate over each particle to compute strain rate
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_strain, std::placeholders::_1, dt_));

    // Iterate over each particle to update solid skeleton volume
    mesh_->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::update_volume,
                  std::placeholders::_1),
        std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                  solid));

    // Iterate over each particle to compute shear (deviatoric) stress
    mesh_->iterate_over_particles(std::bind(
        &mpm::ParticleBase<Tdim>::compute_stress, std::placeholders::_1));

    // Pressure smoothing if needed for the solid phase
    if (pressure_smoothing_) this->pressure_smoothing(solid);

    // Iterate over each particle to update porosity of the solid skeleton
    mesh_->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::update_porosity,
                  std::placeholders::_1, dt_),
        std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                  solid));

    // Compute nodal volume from gauss quadrature
    this->compute_nodes_gauss_volume();

    // Map porosity from solid to fluid particles
    this->compute_fluid_particle_porosity();

    // Spawn a task for external force
#pragma omp parallel sections
    {
#pragma omp section
      {
        // Iterate over particles to compute nodal body force
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_body_force,
                      std::placeholders::_1, this->gravity_));

        // Apply particle traction and map to nodes
        mesh_->apply_traction_on_particles(this->step_ * this->dt_);
      }

#pragma omp section
      {
        // Spawn a task for internal force
        // Iterate over each particle to compute nodal internal force
        mesh_->iterate_over_particles(
            std::bind(&mpm::ParticleBase<Tdim>::map_internal_force,
                      std::placeholders::_1));
      }

#pragma omp section
      {
        // Iterate over solid particles to compute nodal drag force coefficient
        mesh_->iterate_over_particles_predicate(
            std::bind(&mpm::ParticleBase<Tdim>::map_drag_force_coefficient,
                      std::placeholders::_1),
            std::bind(&mpm::ParticleBase<Tdim>::phase_status,
                      std::placeholders::_1, solid));
      }
    }  // Wait for tasks to finish

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce external force of mixture
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    solid),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, solid,
                    std::placeholders::_2));
      // MPI all reduce external force of pore fluid
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::external_force, std::placeholders::_1,
                    liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_external_force,
                    std::placeholders::_1, false, liquid,
                    std::placeholders::_2));

      // MPI all reduce internal force of mixture
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    solid),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, solid,
                    std::placeholders::_2));
      // MPI all reduce internal force of pore liquid
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::internal_force, std::placeholders::_1,
                    liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_internal_force,
                    std::placeholders::_1, false, liquid,
                    std::placeholders::_2));

      // MPI all reduce drag force
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::drag_force_coefficient,
                    std::placeholders::_1),
          std::bind(&mpm::NodeBase<Tdim>::update_drag_force_coefficient,
                    std::placeholders::_1, false, std::placeholders::_2));
    }
#endif

    // Compute intermediate velocity and acceleration for both phases
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<Tdim>::
                compute_acceleration_velocity_twophase_twopoint_prediction,
            std::placeholders::_1, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Reinitialise system matrix to perform PPE
    bool matrix_reinitialization_status = this->reinitialise_matrix();
    if (!matrix_reinitialization_status) {
      status = false;
      throw std::runtime_error("Reinitialisation of matrix failed");
    }

    // Compute poisson equation
    this->compute_poisson_equation();

    // Assign pressure to nodes
    mesh_->iterate_over_nodes_predicate(
        std::bind(&mpm::NodeBase<Tdim>::update_pressure_increment,
                  std::placeholders::_1, assembler_->pressure_increment(),
                  liquid, this->step_ * this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Use nodal pressure to update particle pressure
    mesh_->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_pressure,
                  std::placeholders::_1),
        std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                  liquid));

    // Compute correction force
    this->compute_correction_force();

#ifdef USE_MPI
    // Run if there is more than a single MPI task
    if (mpi_size > 1) {
      // MPI all reduce correction force
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::correction_force,
                    std::placeholders::_1, solid),
          std::bind(&mpm::NodeBase<Tdim>::update_correction_force,
                    std::placeholders::_1, false, solid,
                    std::placeholders::_2));
      // MPI all reduce correction force
      mesh_->template nodal_halo_exchange<Eigen::Matrix<double, Tdim, 1>, Tdim>(
          std::bind(&mpm::NodeBase<Tdim>::correction_force,
                    std::placeholders::_1, liquid),
          std::bind(&mpm::NodeBase<Tdim>::update_correction_force,
                    std::placeholders::_1, false, liquid,
                    std::placeholders::_2));
    }
#endif

    // Compute corrected acceleration and velocity in solid and liquid phase
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<
                Tdim>::compute_acceleration_velocity_semi_implicit_corrector,
            std::placeholders::_1, solid, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));
    mesh_->iterate_over_nodes_predicate(
        std::bind(
            &mpm::NodeBase<
                Tdim>::compute_acceleration_velocity_semi_implicit_corrector,
            std::placeholders::_1, liquid, this->dt_),
        std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

    // Update particle position and kinematics
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::compute_updated_position,
                  std::placeholders::_1, this->dt_, velocity_update_));

    // Apply particle velocity constraints
    mesh_->apply_particle_velocity_constraints();

    // Pressure smoothing
    // FIXME: Pressure smoothing for liquid only wont work due to memory issue
    if (pressure_smoothing_) this->pressure_smoothing(liquid);

    // Reset particle phase container in cell
    mesh_->iterate_over_cells(std::bind(&mpm::Cell<Tdim>::clear_particle_phase,
                                        std::placeholders::_1));

    // Locate particle
    auto unlocatable_particles = mesh_->locate_particles_mesh();

    if (!unlocatable_particles.empty() && this->locate_particles_)
      throw std::runtime_error("Particle outside the mesh domain");
    // If unable to locate particles remove particles
    if (!unlocatable_particles.empty() && !this->locate_particles_)
      for (const auto& remove_particle : unlocatable_particles)
        mesh_->remove_particle(remove_particle);

#ifdef USE_MPI
#ifdef USE_GRAPH_PARTITIONING
    mesh_->transfer_halo_particles();
    MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

    if (step_ % output_steps_ == 0) {
      // HDF5 outputs
      this->write_hdf5(this->step_, this->nsteps_);
#ifdef USE_VTK
      // VTK outputs
      this->write_vtk(this->step_, this->nsteps_);
#endif
    }
  }
  auto solver_end = std::chrono::steady_clock::now();
  console_->info(
      "Rank {}, SemiImplicit_TwoPhase_TwoPoint solver duration: {} ms",
      mpi_rank,
      std::chrono::duration_cast<std::chrono::milliseconds>(solver_end -
                                                            solver_begin)
          .count());

  return status;
}

// Semi-implicit functions
// Initialise matrix
template <unsigned Tdim>
bool mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::initialise_matrix() {
  bool status = true;
  try {
    // Max iteration steps
    unsigned max_iter =
        analysis_["linear_solver"]["max_iter"].template get<unsigned>();
    // Tolerance
    double tolerance =
        analysis_["linear_solver"]["tolerance"].template get<double>();
    // Get matrix assembler type
    std::string assembler_type = analysis_["linear_solver"]["assembler_type"]
                                     .template get<std::string>();
    // Get matrix solver type
    std::string solver_type =
        analysis_["linear_solver"]["solver_type"].template get<std::string>();
    // Create matrix assembler
    assembler_ =
        Factory<mpm::AssemblerBase<Tdim>>::instance()->create(assembler_type);
    // Create matrix solver
    linear_solver_ =
        Factory<mpm::SolverBase<Eigen::SparseMatrix<double>>, unsigned,
                double>::instance()
            ->create(solver_type, std::move(max_iter), std::move(tolerance));
    // Assign mesh pointer to assembler
    assembler_->assign_mesh_pointer(mesh_);

    // Get method to detect free surface detection
    free_surface_detection_ = "density";
    if (analysis_["free_surface_detection"].contains("type"))
      free_surface_detection_ = analysis_["free_surface_detection"]["type"]
                                    .template get<std::string>();
    // Get volume tolerance for free surface
    volume_tolerance_ = analysis_["free_surface_detection"]["volume_tolerance"]
                            .template get<double>();

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Reinitialise and resize matrices at the beginning of every time step
template <unsigned Tdim>
bool mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::reinitialise_matrix() {

  bool status = true;
  try {
    // Assigning matrix id (in each MPI rank)
    // FIXME: NEED TO REDUCE THE NUMBER OF ACTIVE ID TO ONLY NODES WITH FLUID
    // PARTICLES
    const auto nactive_node = mesh_->assign_active_nodes_id();

    // Assigning matrix id globally (required for rank-to-global mapping)
    // FIXME: NEED TO REDUCE THE NUMBER OF ACTIVE ID TO ONLY NODES WITH FLUID
    // PARTICLES
    unsigned nglobal_active_node = nactive_node;
#ifdef USE_MPI
    nglobal_active_node = mesh_->assign_global_active_nodes_id();
#endif

    // Assign global node indice
    assembler_->assign_global_node_indices(nactive_node, nglobal_active_node);

    // Assign pressure constraints
    assembler_->assign_pressure_constraints(this->beta_,
                                            this->step_ * this->dt_);

    // Initialise element matrix
    mesh_->iterate_over_cells(
        std::bind(&mpm::Cell<Tdim>::initialise_element_matrix_twophase,
                  std::placeholders::_1));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Compute poisson equation
template <unsigned Tdim>
bool mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::compute_poisson_equation(
    std::string solver_type) {
  bool status = true;
  try {
    // Construct local cell laplacian matrix
    mesh_->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::map_laplacian_to_cell,
                  std::placeholders::_1),
        std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                  mpm::ParticlePhase::Liquid));

    // Assemble global laplacian matrix
    assembler_->assemble_laplacian_matrix(dt_);

    // Map Poisson RHS matrix
    mesh_->iterate_over_particles_predicate(
        std::bind(&mpm::ParticleBase<Tdim>::map_poisson_right_to_cell,
                  std::placeholders::_1),
        std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                  mpm::ParticlePhase::Liquid));

    // Assemble poisson RHS vector
    assembler_->assemble_poisson_right(dt_);

    // Assign free surface to assembler
    assembler_->assign_free_surface(mesh_->free_surface_nodes());

    // Apply constraints
    assembler_->apply_pressure_constraints();

#ifdef USE_MPI
    // Assign global active dof to solver
    linear_solver_->assign_global_active_dof(assembler_->global_active_dof());

    // Assign rank global mapper to solver
    linear_solver_->assign_rank_global_mapper(assembler_->rank_global_mapper());
#endif

    // Solve matrix equation and assign solution to assembler
    assembler_->assign_pressure_increment(
        linear_solver_->solve(assembler_->laplacian_matrix(),
                              assembler_->poisson_rhs_vector(), solver_type));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute correction force
template <unsigned Tdim>
bool mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::compute_correction_force() {
  bool status = true;
  try {
    // Map correction matrix from particles to cell
    mesh_->iterate_over_particles(
        std::bind(&mpm::ParticleBase<Tdim>::map_correction_matrix_to_cell,
                  std::placeholders::_1));

    // Assemble correction matrix
    assembler_->assemble_corrector_right(dt_);

    // Assign correction force
    mesh_->compute_nodal_correction_force(
        assembler_->correction_matrix(mpm::ParticlePhase::Solid),
        assembler_->pressure_increment(), dt_, mpm::ParticlePhase::Solid);

    mesh_->compute_nodal_correction_force(
        assembler_->correction_matrix(mpm::ParticlePhase::Liquid),
        assembler_->pressure_increment(), dt_, mpm::ParticlePhase::Liquid);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Compute gauss volume in node
template <unsigned Tdim>
void mpm::MPMSemiImplicitTwoPhaseTwoPoint<Tdim>::compute_nodes_gauss_volume() {
  int mpi_size = 1;
#ifdef USE_MPI
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  //! Compute nodal volume from gauss points
  mesh_->iterate_over_cells(std::bind(
      &mpm::Cell<Tdim>::map_cell_gauss_volume_to_nodes, std::placeholders::_1));

#ifdef USE_MPI
  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal volume
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::gauss_volume, std::placeholders::_1),
        std::bind(&mpm::NodeBase<Tdim>::update_gauss_volume,
                  std::placeholders::_1, false, std::placeholders::_2));
  }
#endif
}

//! Compute correction force
template <unsigned Tdim>
void mpm::MPMSemiImplicitTwoPhaseTwoPoint<
    Tdim>::compute_fluid_particle_porosity() {
  int mpi_size = 1;
#ifdef USE_MPI
  // Get number of MPI ranks
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
#endif

  //! Map solid volume fraction from solid particle to the nodes
  mesh_->iterate_over_particles_predicate(
      std::bind(&mpm::ParticleBase<Tdim>::map_volume_fraction_to_nodes,
                std::placeholders::_1),
      std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                mpm::ParticlePhase::Solid));

#ifdef USE_MPI
  // Run if there is more than a single MPI task
  if (mpi_size > 1) {
    // MPI all reduce nodal mass
    mesh_->template nodal_halo_exchange<double, 1>(
        std::bind(&mpm::NodeBase<Tdim>::volume_fraction, std::placeholders::_1,
                  mpm::ParticlePhase::Solid),
        std::bind(&mpm::NodeBase<Tdim>::update_volume_fraction,
                  std::placeholders::_1, false, mpm::ParticlePhase::Solid,
                  std::placeholders::_2));
  }
#endif

  //! Compute porosity in the node
  mesh_->iterate_over_nodes_predicate(
      std::bind(&mpm::NodeBase<Tdim>::compute_porosity, std::placeholders::_1),
      std::bind(&mpm::NodeBase<Tdim>::status, std::placeholders::_1));

  //! Interpolate fluid porosity from the nodal porosity
  mesh_->iterate_over_particles_predicate(
      std::bind(&mpm::ParticleBase<Tdim>::assign_porosity,
                std::placeholders::_1),
      std::bind(&mpm::ParticleBase<Tdim>::phase_status, std::placeholders::_1,
                mpm::ParticlePhase::Liquid));
}
