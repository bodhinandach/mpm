#include "particle.h"
#include "factory.h"
#include "particle_base.h"
#include "particle_fluid.h"
#include "particle_fluid_twophase.h"
#include "particle_solid_twophase.h"
#include "particle_twophase.h"

namespace mpm {
// ParticleType
std::map<std::string, int> ParticleType = {{"P2D", 0},
                                           {"P3D", 1},
                                           {"P2DFLUID", 2},
                                           {"P3DFLUID", 3},
                                           {"P2D2PHASE", 4},
                                           {"P3D2PHASE", 5},
                                           {"P2DSOLID2PHASE", 6},
                                           {"P3DSOLID2PHASE", 7},
                                           {"P2DFLUID2PHASE", 8},
                                           {"P3DFLUID2PHASE", 9}};
std::map<int, std::string> ParticleTypeName = {{0, "P2D"},
                                               {1, "P3D"},
                                               {2, "P2DFLUID"},
                                               {3, "P3DFLUID"},
                                               {4, "P2D2PHASE"},
                                               {5, "P3D2PHASE"},
                                               {6, "P2DSOLID2PHASE"},
                                               {7, "P3DSOLID2PHASE"},
                                               {8, "P2DFLUID2PHASE"},
                                               {9, "P3DFLUID2PHASE"}};
std::map<std::string, std::string> ParticleHDF5TypeName = {
    {"P2D", "particles"},
    {"P3D", "particles"},
    {"P2DFLUID", "fluid_particles"},
    {"P3DFLUID", "fluid_particles"},
    {"P2D2PHASE", "twophase_particles"},
    {"P3D2PHASE", "twophase_particles"},
    {"P2DSOLID2PHASE", "twophase_solid_particles"},
    {"P3DSOLID2PHASE", "twophase_solid_particles"},
    {"P2DFLUID2PHASE", "twophase_fluid_particles"},
    {"P3DFLUID2PHASE", "twophase_fluid_particles"}};
}  // namespace mpm

// Single-phase Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::Particle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("P2D");

// Single-phase Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::Particle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("P3D");

// Single-phase Fluid Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::FluidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dfluid("P2DFLUID");

// Single-phase Fluid Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::FluidParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dfluid("P3DFLUID");

// Two-phase Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::TwoPhaseParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d2phase("P2D2PHASE");

// Two-phase Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::TwoPhaseParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d2phase("P3D2PHASE");

// Two-phase Solid Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::TwoPhaseSolidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dsolid2phase("P2DSOLID2PHASE");

// Two-phase Solid Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::TwoPhaseSolidParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dsolid2phase("P3DSOLID2PHASE");

// Two-phase Fluid Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::TwoPhaseFluidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dfluid2phase("P2DFLUID2PHASE");

// Two-phase Fluid Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::TwoPhaseFluidParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dfluid2phase("P3DFLUID2PHASE");
