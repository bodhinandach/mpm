//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::TwoPhaseSolidParticle<Tdim>::TwoPhaseSolidParticle(Index id,
                                                        const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

  // Logger
  std::string logger = "TwoPhaseSolidParticle" + std::to_string(Tdim) +
                       "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::assign_cell(
    const std::shared_ptr<Cell<Tdim>>& cellptr) {
  bool status = mpm::Particle<Tdim>::assign_cell(cellptr);
  status = cellptr->add_particle_phase(this->id(), this->phase());
  return status;
}

// Assign a cell to particle
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::assign_cell_xi(
    const std::shared_ptr<Cell<Tdim>>& cellptr,
    const Eigen::Matrix<double, Tdim, 1>& xi) {
  bool status = mpm::Particle<Tdim>::assign_cell_xi(cellptr, xi);
  status = cellptr->add_particle_phase(this->id(), this->phase());
  return status;
}

// Compute volume of the particle
template <unsigned Tdim>
void mpm::TwoPhaseSolidParticle<Tdim>::compute_volume() noexcept {
  // Check if particle has a valid cell ptr
  assert(cell_ != nullptr);
  // Volume of the cell / # of particles
  this->assign_volume(cell_->volume() / cell_->nparticles(this->phase()));
}

// Compute mass of solid particle
template <unsigned Tdim>
void mpm::TwoPhaseSolidParticle<Tdim>::compute_mass() noexcept {
  // Check if particle volume is set and material ptr is valid
  assert(volume_ != std::numeric_limits<double>::max() &&
         this->material() != nullptr);
  // Mass = volume of particle * mass_density
  this->mass_density_ =
      (this->material())->template property<double>(std::string("density")) *
      (1 - this->porosity_);
  this->mass_ = volume_ * mass_density_;
}

// Assign porosity to the particle
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::assign_porosity() {
  bool status = true;
  try {
    if (this->material() != nullptr) {
      this->porosity_ =
          this->material()->template property<double>(std::string("porosity"));
      if (porosity_ < 0. || porosity_ > 1.)
        throw std::runtime_error(
            "Particle porosity is negative or larger than one");
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

// Update material point porosity
template <unsigned Tdim>
void mpm::TwoPhaseSolidParticle<Tdim>::update_porosity(double dt) {
  // Update particle porosity
  const double porosity =
      1 - (1 - this->porosity_) / (1 + dt * strain_rate_.head(Tdim).sum());
  // Check if the value is valid
  if (porosity < 0.)
    this->porosity_ = 1.E-5;
  else if (porosity > 1.)
    this->porosity_ = 1 - 1.E-5;
  else
    this->porosity_ = porosity;
}

//! Map particle volume fraction to nodes
template <unsigned Tdim>
void mpm::TwoPhaseSolidParticle<Tdim>::map_volume_fraction_to_nodes() {
  // Map porosity to nodes
  const double tolerance = std::numeric_limits<double>::epsilon();
  const double solid_volume = volume_ * (1. - porosity_);
  for (unsigned i = 0; i < nodes_.size(); ++i) {
    if (nodes_[i]->gauss_volume() > tolerance)
      nodes_[i]->update_volume_fraction(
          true, this->phase(),
          solid_volume / nodes_[i]->gauss_volume() * shapefn_[i]);
    else
      throw std::runtime_error("Nodal gauss weight is invalid or unassigned");
  }
};

//! Assign particle permeability
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::assign_permeability() {
  bool status = true;
  try {
    // Check if material ptr is valid
    if (this->material() != nullptr) {
      // Get initial porosity
      double porosity =
          this->material()->template property<double>(std::string("porosity"));
      if (porosity < 0. || porosity > 1.)
        throw std::runtime_error(
            "Particle porosity is negative or larger than one, can not assign "
            "permeability");
      // Porosity parameter
      const double k_p = std::pow(porosity, 3) / std::pow((1. - porosity), 2);
      // Assign permeability
      switch (Tdim) {
        case (3):
          // Check if the permeability is valid
          if (this->material()->template property<double>("k_z") < 0)
            throw std::runtime_error("Material's permeability k_z is invalid");
          permeability_(2) =
              this->material()->template property<double>("k_z") / k_p;
        case (2):
          // Check if the permeability is valid
          if (this->material()->template property<double>("k_y") < 0)
            throw std::runtime_error("Material's permeability k_y is invalid");
          permeability_(1) =
              this->material()->template property<double>("k_y") / k_p;
        case (1):
          // Check if the permeability is valid
          if (this->material(mpm::ParticlePhase::Solid)
                  ->template property<double>("k_x") < 0)
            throw std::runtime_error("Material's permeability k_x is invalid");
          // Assign permeability
          permeability_(0) =
              this->material()->template property<double>("k_x") / k_p;
      }
    } else {
      throw std::runtime_error("Material is invalid");
    }
  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map drag force
//! NOTE: This function assume linear darcy law
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::map_drag_force_coefficient() {
  bool status = true;
  try {
    // Porosity parameter
    const double k_p =
        std::pow(this->porosity_, 3) / std::pow((1. - this->porosity_), 2);
    // Initialise drag force coefficient
    VectorDim drag_force_coefficient;
    drag_force_coefficient.setZero();

    // Check if permeability coefficient is valid
    const double gravity = 9.81;
    for (unsigned i = 0; i < Tdim; ++i) {
      if (k_p > 0.)
        drag_force_coefficient(i) = gravity / (k_p * permeability_(i));
      else
        throw std::runtime_error("Porosity coefficient is invalid");
    }

    // Map drag forces from particle to nodes
    for (unsigned j = 0; j < nodes_.size(); ++j) {
      const double gauss_volume = nodes_[j]->gauss_volume();
      const double nodal_solid_mass =
          nodes_[j]->mass(mpm::ParticlePhase::Solid);
      const double nodal_liquid_mass =
          nodes_[j]->mass(mpm::ParticlePhase::Liquid);
      const double nodal_solid_volume_fraction =
          nodes_[j]->volume_fraction(mpm::ParticlePhase::Solid);
      const double nodal_liquid_volume_fraction =
          nodes_[j]->volume_fraction(mpm::ParticlePhase::Liquid);

      // Multiply nodal porosity square
      drag_force_coefficient *=
          nodal_liquid_volume_fraction * nodal_liquid_volume_fraction;

      // Multiply nodal volume contribution
      drag_force_coefficient *=
          (gauss_volume *
           (nodal_solid_mass / nodal_solid_volume_fraction /
            this->material(mpm::ParticlePhase::Solid)
                ->template property<double>(std::string("density"))) *
           (nodal_liquid_mass / nodal_liquid_volume_fraction));

      nodes_[j]->update_drag_force_coefficient(true, drag_force_coefficient);
    }

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
    status = false;
  }
  return status;
}

//! Map correction matrix element matrix to cell (used to correct velocity)
template <unsigned Tdim>
bool mpm::TwoPhaseSolidParticle<Tdim>::map_correction_matrix_to_cell() {
  bool status = true;
  try {
    cell_->compute_local_correction_matrix_twophase(
        this->phase(), shapefn_, dn_dx_, volume_, 1.0 - this->porosity_);

  } catch (std::exception& exception) {
    console_->error("{} #{}: {}\n", __FILE__, __LINE__, exception.what());
  }
  return status;
}