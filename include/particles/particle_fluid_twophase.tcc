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
