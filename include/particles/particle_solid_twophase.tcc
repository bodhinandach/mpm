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
