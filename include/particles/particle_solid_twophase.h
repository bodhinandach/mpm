#ifndef MPM_TWOPHASE_SOLID_PARTICLE_H_
#define MPM_TWOPHASE_SOLID_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "logger.h"
#include "particle.h"

namespace mpm {

//! Solid Particle class for twophase-twopoint solver
//! \brief Class with function specific to solid particles
//! \tparam Tdim Dimension
template <unsigned Tdim>
class TwoPhaseSolidParticle : public mpm::Particle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  TwoPhaseSolidParticle(Index id, const VectorDim& coord);

  //! Destructor
  ~TwoPhaseSolidParticle() override{};

  //! Delete copy constructor
  TwoPhaseSolidParticle(const TwoPhaseSolidParticle<Tdim>&) = delete;

  //! Delete assignment operator
  TwoPhaseSolidParticle& operator=(const TwoPhaseSolidParticle<Tdim>&) = delete;

  //! Type of particle
  std::string type() const override {
    return (Tdim == 2) ? "P2DSOLID2PHASE" : "P3DSOLID2PHASE";
  }

  //! Return phase boolean
  //! \retval phase boolean associated with phase index
  bool phase_status(unsigned phase) const override {
    return (phase == mpm::ParticlePhase::Solid) ? true : false;
  };

  //! Compute solid mass
  void compute_mass() noexcept override;

  //! Assign particle permeability
  //! \retval status Assignment status
  bool assign_permeability() override;

  //! Assign porosity
  bool assign_porosity() override;

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
  //! Material point porosity (volume of voids / total volume)
  double porosity_{0.0};
  //! Permeability parameter c1 (k = k_p * c1)
  VectorDim permeability_;

  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // TwoPhaseSolidParticle class
}  // namespace mpm

#include "particle_solid_twophase.tcc"

#endif  // MPM_TWOPHASE_SOLID_PARTICLE_H_
