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

  status = true;
  return status;
}
