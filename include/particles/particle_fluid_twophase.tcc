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
