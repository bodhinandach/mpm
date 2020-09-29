//! Construct a two phase particle with id and coordinates
template <unsigned Tdim>
mpm::SolidParticle<Tdim>::SolidParticle(Index id, const VectorDim& coord)
    : mpm::Particle<Tdim>(id, coord) {

  // Logger
  std::string logger =
      "SolidParticle" + std::to_string(Tdim) + "d::" + std::to_string(id);
  console_ = std::make_unique<spdlog::logger>(logger, mpm::stdout_sink);
}
