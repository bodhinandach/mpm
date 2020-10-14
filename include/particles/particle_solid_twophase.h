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

  //! Return nodal phase index associated with particle
  //! \retval phase Nodal phase index
  unsigned phase() const override { return mpm::ParticlePhase::Solid; };

  //! Assigning beta parameter to particle
  //! \param[in] pressure parameter determining type of projection
  void assign_projection_parameter(double parameter) override {
    this->projection_param_ = parameter;
  };

  //! Assign a cell to particle
  //! If point is in new cell, assign new cell and remove particle id from old
  //! cell. If point can't be found in the new cell, check if particle is still
  //! valid in the old cell, if it is leave it as is. If not, set cell as null
  //! \param[in] cellptr Pointer to a cell
  //! \details the overriden function also add additional particle phase for
  //! twophase coupling purposes
  bool assign_cell(const std::shared_ptr<Cell<Tdim>>& cellptr) override;

  //! Assign a cell to particle
  //! If point is in new cell, assign new cell and remove particle id from old
  //! cell. If point can't be found in the new cell, check if particle is still
  //! valid in the old cell, if it is leave it as is. If not, set cell as null
  //! \param[in] cellptr Pointer to a cell
  //! \param[in] xi Local coordinates of the point in reference cell
  //! \details the overriden function also add additional particle phase for
  //! twophase coupling purposes
  bool assign_cell_xi(const std::shared_ptr<Cell<Tdim>>& cellptr,
                      const Eigen::Matrix<double, Tdim, 1>& xi) override;

  //! Compute volume as cell volume / nparticles
  void compute_volume() noexcept override;

  //! Compute solid mass
  void compute_mass() noexcept override;

  //! Assign porosity
  bool assign_porosity() override;

  //! Update porosity
  //! \param[in] dt Analysis time step
  void update_porosity(double dt) override;

  //! Assign particle permeability
  //! \retval status Assignment status
  bool assign_permeability() override;

  //! Map particle volume fraction to nodes
  void map_volume_fraction_to_nodes() override;

  //! Map drag force coefficient
  bool map_drag_force_coefficient() override;

  //! Map correction matrix element matrix to cell (used to correct velocity)
  bool map_correction_matrix_to_cell() override;

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
