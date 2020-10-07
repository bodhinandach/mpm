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

  //! Assign free surface node id
  void assign_free_surface(
      const std::set<mpm::Index>& free_surface_id) override {
    free_surface_ = free_surface_id;
  }

  //! Return pressure increment
  Eigen::VectorXd& pressure_increment() override { return pressure_increment_; }

  //! Assign pressure increment
  void assign_pressure_increment(
      const Eigen::VectorXd& pressure_increment) override {
    pressure_increment_ = pressure_increment;
  }

  //! Return correction matrix
  Eigen::SparseMatrix<double>& correction_matrix() override {
    return correction_matrix_;
  }

  //! Assemble corrector RHS
  bool assemble_corrector_right(double dt) override;

  //! Return the total size of global dof in all rank
  unsigned global_active_dof() override { return global_active_dof_; };

  //! Return a vector to map local (rank) index to global index
  std::vector<int> rank_global_mapper() override {
    return rank_global_mapper_;
  };

 protected:
  //! number of nodes
  using AssemblerBase<Tdim>::active_dof_;
  //! Mesh object
  using AssemblerBase<Tdim>::mesh_;
  //! Logger
  std::shared_ptr<spdlog::logger> console_;
  //! Global node indices
  std::vector<Eigen::VectorXi> global_node_indices_;
  //! Laplacian matrix
  Eigen::SparseMatrix<double> laplacian_matrix_;
  //! Poisson RHS vector
  Eigen::VectorXd poisson_rhs_vector_;
  //! Free surface
  std::set<mpm::Index> free_surface_;
  //! Pressure constraints
  Eigen::SparseVector<double> pressure_constraints_;
  //! \delta p^(t+1) = p^(t+1) - beta * p^(t)
  Eigen::VectorXd pressure_increment_;
  //! correction_matrix
  Eigen::SparseMatrix<double> correction_matrix_;
  //! Number of total active_dof in all rank
  unsigned global_active_dof_;
  //! Rank to Global mapper
  std::vector<int> rank_global_mapper_;
};
}  // namespace mpm

#include "assembler_eigen_semi_implicit_twophase_twopoint.tcc"
#endif  // MPM_ASSEMBLER_EIGEN_SEMI_IMPLICIT_TWOPHASE_TWOPOINT_H_
