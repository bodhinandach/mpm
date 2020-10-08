#ifndef MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_TWOPOINT_H_
#define MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_TWOPOINT_H_

#include <Eigen/Sparse>
#include <string>

// Speed log
#include "assembler_base.h"
#include "spdlog/spdlog.h"

#include "mesh.h"

namespace mpm {
template <unsigned Tdim>
class AssemblerEigenSemiImplicitTwoPhaseTwoPoint
    : public AssemblerEigenSemiImplicitNavierStokes<Tdim> {
 public:
  //! Constructor
  AssemblerEigenSemiImplicitTwoPhaseTwoPoint();

  //! Assemble laplacian matrix
  bool assemble_laplacian_matrix(double dt) override;

  //! Assemble poisson RHS vector
  bool assemble_poisson_right(double dt) override;

  //! Assemble corrector RHS
  bool assemble_corrector_right(double dt) override;

  //! Return correction matrix
  Eigen::SparseMatrix<double>& correction_matrix() override {
    return this->correction_matrix(mpm::ParticlePhase::Solid);
  }

  //! Return correction matrix with phase index
  Eigen::SparseMatrix<double>& correction_matrix(unsigned phase) override {
    return correction_matrix_[phase];
  }

 protected:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Global node indices
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::global_node_indices_;
  //! Logger
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::console_;
  //! Laplacian matrix
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::laplacian_matrix_;
  //! Poisson RHS vector
  using AssemblerEigenSemiImplicitNavierStokes<Tdim>::poisson_rhs_vector_;
  //! Correction_matrix for twophase
  std::vector<Eigen::SparseMatrix<double>> correction_matrix_;
};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_twophase_twopoint.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_TWOPOINT_H_
