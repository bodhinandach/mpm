//! Construct a semi-implicit eigen matrix assembler
template <unsigned Tdim>
mpm::AssemblerEigenSemiImplicitTwoPhaseTwoPoint<
    Tdim>::AssemblerEigenSemiImplicitTwoPhaseTwoPoint()
    : mpm::AssemblerEigenSemiImplicitNavierStokes<Tdim>() {
  //! Logger
  console_ =
      spdlog::stdout_color_mt("AssemblerEigenSemiImplicitTwoPhaseTwoPoint");
}

//! Assemble Laplacian matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhaseTwoPoint<
    Tdim>::assemble_laplacian_matrix(double dt) {
  bool status = true;
  try {
    // Initialise Laplacian matrix
    // TODO: FIX to only with liquid particles
    laplacian_matrix_.resize(active_dof_, active_dof_);
    laplacian_matrix_.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        laplacian_matrix_.reserve(Eigen::VectorXi::Constant(active_dof_, 30));
        break;
      }
    }

    // Cell pointer
    // TODO: FIX to those with liquid particles
    const auto& cells = mesh_->cells();

    // Iterate over cells
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);

        // Laplacian element of cell
        const auto cell_laplacian = (*cell_itr)->laplacian_matrix();

        // Assemble global laplacian matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            laplacian_matrix_.coeffRef(global_node_indices_.at(cid)(i),
                                       global_node_indices_.at(cid)(j)) +=
                cell_laplacian(i, j);
          }
        }

        ++cid;
      }
    }

    laplacian_matrix_ *= dt;

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhaseTwoPoint<
    Tdim>::assemble_poisson_right(double dt) {
  bool status = true;
  try {
    // Phase index
    const unsigned solid = mpm::ParticlePhase::Solid;
    const unsigned liquid = mpm::ParticlePhase::Liquid;

    // Initialise Poisson RHS matrices
    // Solid part
    Eigen::SparseMatrix<double> solid_poisson_right_matrix;
    solid_poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    solid_poisson_right_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        solid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        solid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Liquid part
    Eigen::SparseMatrix<double> liquid_poisson_right_matrix;
    liquid_poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    liquid_poisson_right_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        liquid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        liquid_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Coupling part
    Eigen::SparseMatrix<double> coupling_poisson_right_matrix;
    coupling_poisson_right_matrix.resize(active_dof_, active_dof_ * Tdim);
    coupling_poisson_right_matrix.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        coupling_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        coupling_poisson_right_matrix.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Iterate over cells
    mpm::Index cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        // Node ids in each cell
        const auto nids = global_node_indices_.at(cid);

        // Local Poisson RHS matrices
        auto cell_poisson_right_solid =
            (*cell_itr)->poisson_right_matrix(solid);
        auto cell_poisson_right_liquid =
            (*cell_itr)->poisson_right_matrix(liquid);
        auto cell_poisson_right_coupling =
            (*cell_itr)->poisson_right_coupling_matrix();

        // Assemble global poisson RHS matrix
        for (unsigned i = 0; i < nids.size(); ++i) {
          for (unsigned j = 0; j < nids.size(); ++j) {
            for (unsigned k = 0; k < Tdim; ++k) {
              solid_poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_solid(i, j + k * nids.size());

              liquid_poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_liquid(i, j + k * nids.size());

              coupling_poisson_right_matrix.coeffRef(
                  global_node_indices_.at(cid)(i),
                  global_node_indices_.at(cid)(j) + k * active_dof_) +=
                  cell_poisson_right_coupling(i, j + k * nids.size());
            }
          }
        }
        cid++;
      }
    }

    // Resize poisson right vector
    poisson_rhs_vector_.resize(active_dof_);
    poisson_rhs_vector_.setZero();

    // Compute solid and liquid velocity
    Eigen::MatrixXd solid_velocity;
    solid_velocity.resize(active_dof_, Tdim);
    solid_velocity.setZero();

    Eigen::MatrixXd liquid_velocity;
    liquid_velocity.resize(active_dof_, Tdim);
    liquid_velocity.setZero();

    // Active nodes
    const auto& active_nodes = mesh_->active_nodes();
    unsigned node_index = 0;

    for (auto node_itr = active_nodes.cbegin(); node_itr != active_nodes.cend();
         ++node_itr) {
      // Compute nodal intermediate force
      solid_velocity.row(node_index) = (*node_itr)->velocity(solid).transpose();
      liquid_velocity.row(node_index) =
          (*node_itr)->velocity(liquid).transpose();
      node_index++;
    }

    // Resize velocity
    solid_velocity.resize(active_dof_ * Tdim, 1);
    liquid_velocity.resize(active_dof_ * Tdim, 1);

    // Compute poisson RHS vector
    poisson_rhs_vector_ =
        -(solid_poisson_right_matrix * solid_velocity) -
        (liquid_poisson_right_matrix * liquid_velocity) -
        (coupling_poisson_right_matrix * (liquid_velocity - solid_velocity));

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Assemble corrector right matrix
template <unsigned Tdim>
bool mpm::AssemblerEigenSemiImplicitTwoPhaseTwoPoint<
    Tdim>::assemble_corrector_right(double dt) {
  bool status = true;
  try {
    // Resize correction matrix
    correction_matrix_.resize(active_dof_, active_dof_ * Tdim);
    correction_matrix_.setZero();

    // Reserve storage for sparse matrix
    switch (Tdim) {
      // For 2d: 10 entries /column
      case (2): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 10));
        break;
      }
      // For 3d: 30 entries /column
      case (3): {
        correction_matrix_.reserve(
            Eigen::VectorXi::Constant(active_dof_ * Tdim, 30));
        break;
      }
    }

    // Cell pointer
    const auto& cells = mesh_->cells();

    // Iterate over cells
    unsigned cid = 0;
    for (auto cell_itr = cells.cbegin(); cell_itr != cells.cend(); ++cell_itr) {
      if ((*cell_itr)->status()) {
        unsigned nnodes_per_cell = global_node_indices_.at(cid).size();
        auto cell_correction_matrix = (*cell_itr)->correction_matrix();
        for (unsigned k = 0; k < Tdim; k++) {
          for (unsigned i = 0; i < nnodes_per_cell; i++) {
            for (unsigned j = 0; j < nnodes_per_cell; j++) {
              // Fluid
              correction_matrix_.coeffRef(
                  global_node_indices_.at(cid)(i),
                  k * active_dof_ + global_node_indices_.at(cid)(j)) +=
                  cell_correction_matrix(i, j + k * nnodes_per_cell);
            }
          }
        }
        cid++;
      }
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}
