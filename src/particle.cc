#include "particle.h"
#include "factory.h"
#include "particle_base.h"
#include "particle_fluid.h"
#include "particle_twophase.h"
#include "solid_particle.h"

namespace mpm {
// ParticleType
std::map<std::string, int> ParticleType = {
    {"P2D", 0},      {"P3D", 1},      {"P2DSOLID", 2},  {"P3DSOLID", 3},
    {"P2DFLUID", 4}, {"P3DFLUID", 5}, {"P2D2PHASE", 6}, {"P3D2PHASE", 7}};
std::map<int, std::string> ParticleTypeName = {
    {0, "P2D"},      {1, "P3D"},      {2, "P2DSOLID"},  {3, "P3DSOLID"},
    {4, "P2DFLUID"}, {5, "P3DFLUID"}, {6, "P2D2PHASE"}, {7, "P3D2PHASE"}};
std::map<std::string, std::string> ParticleHDF5TypeName = {
    {"P2D", "particles"},
    {"P3D", "particles"},
    {"P2DSOLID", "solid_particles"},
    {"P3DSOLID", "solid_particles"},
    {"P2DFLUID", "fluid_particles"},
    {"P3DFLUID", "fluid_particles"},
    {"P2D2PHASE", "twophase_particles"},
    {"P3D2PHASE", "twophase_particles"}};
}  // namespace mpm

// Single-phase Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::Particle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2d("P2D");

// Single-phase Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::Particle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3d("P3D");

// Solid Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::SolidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dsolid("P2DSOLID");

// Solid Particle3D (3 Dim)
static Register<mpm::ParticleBase<3>, mpm::SolidParticle<3>, mpm::Index,
                const Eigen::Matrix<double, 3, 1>&>
    particle3dsolid("P3DSOLID");

// Fluid Particle2D (2 Dim)
static Register<mpm::ParticleBase<2>, mpm::FluidParticle<2>, mpm::Index,
                const Eigen::Matrix<double, 2, 1>&>
    particle2dfluid("P2DFLUID");

// Fluid Particle3D (3 Dim)
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
