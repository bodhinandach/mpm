//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseFluidParticle<Tdim>::TwoPhaseFluidParticle(Index id,
                                                        const VectorDim& coord)
    : mpm::FluidParticle<Tdim>(id, coord) {

  // Logger
  std::string logger = "TwoPhaseFluidParticle" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = mpm::Particle<Tdim>::assign_cell(cellptr);
  status = cellptr->add_particle_phase(this->id(), this->phase());
  return status;
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::assign_cell_xi(
    const std::shared_ptr<Cell<Tdim>>& cellptr,
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  bool status = mpm::Particle<Tdim>::assign_cell_xi(cellptr, xi);
  status = cellptr->add_particle_phase(this->id(), this->phase());
  return status;
}

// Compute volume of the particle
template <unsigned Tdim>
void mpm::TwoPhaseFluidParticle<Tdim>::compute_volume() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Volume of the cell / # of particles
  this->assign_volume(cell_->volume() / cell_->nparticles(this->phase()));
}

// Compute mass of fluid particle
template <unsigned Tdim>
void mpm::TwoPhaseFluidParticle<Tdim>::compute_mass() noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         this->material() != nullptr);
  // Mass = volume of particle * mass_density
  this->mass_density_ =
      (this->material())->template property<double>(std::string("density")) *
      (this->porosity_);
  this->mass_ = volume_ * mass_density_;
}

// Assign porosity to the particle
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::assign_porosity() {
  // Assert
  assert(cell_ != nullptr);

  bool status = false;
  double porosity = 0.;
  // Update particle pressure to interpolated nodal pressure
  for (unsigned i = 0; i < this->nodes_.size(); ++i)
    porosity += shapefn_[i] * nodes_[i]->volume_fraction(this->phase());

  // Check if the value is valid
  if (porosity < 0.)
    this->porosity_ = 1.E-5;
  else if (porosity > 1.)
    this->porosity_ = 1 - 1.E-5;
  else
    this->porosity_ = porosity;

  // Assign porosity gradient
  VectorDim porosity_gradient;
  porosity_gradient.setZero();
  // Update particle pressure to interpolated nodal pressure
  for (unsigned i = 0; i < this->nodes_.size(); ++i)
    porosity_gradient +=
        (dn_dx_.row(i)).transpose() * nodes_[i]->volume_fraction(this->phase());

  this->porosity_gradient_ = porosity_gradient;

  status = true;
  return status;
}

// Compute stress
// TODO: Investigate the influence of fluid shear stress
template <unsigned Tdim>
void mpm::TwoPhaseFluidParticle<Tdim>::compute_stress() noexcept {
  // // Run particle compute stress
  // mpm::Particle<Tdim>::compute_stress();

  // // Calculate fluid turbulent stress
  // this->stress_ += this->compute_turbulent_stress();
  this->stress_.setZero();
}

//! Map laplacian element matrix to cell (used in poisson equation LHS)
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::map_laplacian_to_cell() {
  bool status = true;
  try {
    // Get intrinsic solid grain density
    // FIXME: Get the solid density from solid particles
    double solid_density = 2600;
    // Compute laplacian matrix multiplier
    double multiplier = (1.0 - this->porosity_) / solid_density +
                        (this->porosity_) / (this->material())
                                                ->template property<double>(
                                                    std::string("density"));
    // Compute local matrix of Laplacian
    cell_->compute_local_laplacian(dn_dx_, volume_, multiplier);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map poisson rhs element matrix to cell (used in poisson equation RHS)
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::map_poisson_right_to_cell() {
  bool status = true;
  try {
    // Compute local poisson rhs matrix for solid part
    cell_->compute_local_poisson_right_twophase(mpm::ParticlePhase::Solid,
                                                shapefn_, dn_dx_, volume_,
                                                1.0 - this->porosity_);
    // Compute local poisson rhs matrix for liquid part
    cell_->compute_local_poisson_right_twophase(
        mpm::ParticlePhase::Liquid, shapefn_, dn_dx_, volume_, this->porosity_);
    // Compute local poisson rhs coupling matrix
    cell_->compute_local_poisson_right_coupling_matrix(porosity_gradient_,
                                                       shapefn_, volume_);
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map correction matrix element matrix to cell (used to correct velocity)
template <unsigned Tdim>
bool mpm::TwoPhaseFluidParticle<Tdim>::map_correction_matrix_to_cell() {
  bool status = true;
  try {
    cell_->compute_local_correction_matrix_twophase(
        this->phase(), shapefn_, dn_dx_, volume_, this->porosity_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}