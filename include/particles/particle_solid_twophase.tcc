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