#ifndef MPM_TWOPHASE_FLUID_PARTICLE_H_
#define MPM_TWOPHASE_FLUID_PARTICLE_H_

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#include "logger.h"
#include "particle.h"
#include "particle_fluid.h"

namespace mpm {

//! Solid Particle class for twophase-twopoint solver
//! \brief Class with function specific to solid particles
//! \tparam Tdim Dimension
template <unsigned Tdim>
class TwoPhaseFluidParticle : public mpm::FluidParticle<Tdim> {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  //! Construct a particle with id and coordinates
  //! \param[in] id Particle id
  //! \param[in] coord Coordinates of the particles
  TwoPhaseFluidParticle(Index id, const VectorDim& coord);

  //! Destructor
  ~TwoPhaseFluidParticle() override{};

  //! Delete copy constructor
  TwoPhaseFluidParticle(const TwoPhaseFluidParticle<Tdim>&) = delete;

  //! Delete assignment operator
  TwoPhaseFluidParticle& operator=(const TwoPhaseFluidParticle<Tdim>&) = delete;

  //! Type of particle
  std::string type() const override {
    return (Tdim == 2) ? "P2DFLUID2PHASE" : "P3DFLUID2PHASE";
  }

  //! Return phase boolean
  //! \retval phase boolean associated with phase index
  bool phase_status(unsigned phase) const override {
    return (phase == mpm::ParticlePhase::Liquid) ? true : false;
  };

  //! Compute fluid mass
  void compute_mass() noexcept override;

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
  using FluidParticle<Tdim>::projection_param_;

  //! Material point porosity (volume of voids / total volume)
  double porosity_{1.0};
  //! Logger
  std::unique_ptr<spdlog::logger> console_;
};  // TwoPhaseFluidParticle class
}  // namespace mpm

#include "particle_fluid_twophase.tcc"

#endif  // MPM_TWOPHASE_FLUID_PARTICLE_H_
