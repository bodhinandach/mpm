#ifndef MPM_SOLID_PARTICLE_H_
#define MPM_SOLID_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "logger.h"
#include "particle.h"

namespace mpm {

//! Solid Particle class
//! \brief Class with function specific to solid particles
//! \tparam Tdim Dimension
template <unsigned Tdim>
class SolidParticle : public mpm::Particle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  SolidParticle(Index id, const VectorDim& coord);

  //! Destructor
  ~SolidParticle() override{};

  //! Delete copy constructor
  SolidParticle(const SolidParticle<Tdim>&) = delete;

  //! Delete assignment operator
  SolidParticle& operator=(const SolidParticle<Tdim>&) = delete;

  //! Type of particle
  std::string type() const override {
    return (Tdim == 2) ? "P2DSOLID" : "P3DSOLID";
  }

 private:
  //! Cell
  using ParticleBase<Tdim>::cell_;
  //! Nodes
  using ParticleBase<Tdim>::nodes_;
  //! Fluid material
  using ParticleBase<Tdim>::material_;
  //! State variables
  using ParticleBase<Tdim>::state_variables_;
  //! Shape functions
  using Particle<Tdim>::shapefn_;
  //! dN/dX
  using Particle<Tdim>::dn_dx_;
  //! Fluid strain rate
  using Particle<Tdim>::strain_rate_;
  //! Effective stress of soil skeleton
  using Particle<Tdim>::stress_;
  //! Particle mass density
  using Particle<Tdim>::mass_density_;
  //! Particle mass density
  using Particle<Tdim>::mass_;
  //! Particle total volume
  using Particle<Tdim>::volume_;
  //! Projection parameter for semi-implicit update
  double projection_param_{0.0};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // SolidParticle class
}  // namespace mpm

#include "solid_particle.tcc"

#endif  // MPM_SOLID_PARTICLE_H_
