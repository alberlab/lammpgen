/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#include <cmath>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include "fix_volumetric_restraint.h"
#include "atom.h"
#include "memory.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

double max(double x, double y);
double norm3(double* x);

FixVolumetricRestraint::FixVolumetricRestraint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  // TODO: vector output.
  scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 1;
  dynamic_group_allow = 1;

  if (narg != 5) error->all(FLERR,"Illegal fix ellipsodal envelope command [args: filename k]");
  k = force->numeric(FLERR,arg[4]);
  std::ifstream fp(arg[3], std::ios::binary);
  grid_info grid;
  fp.read((char*) (&grid), sizeof(grid_info));
  int tot_size =  grid.nx * grid.ny * grid.nz * 4;
  grid_step[0] = ( grid.x1 - grid.x0 ) / grid.nx;
  grid_step[1] = ( grid.y1 - grid.y0 ) / grid.ny;
  grid_step[2] = ( grid.y1 - grid.y0 ) / grid.nz;

  grid_data = new float[tot_size];

  fp.read((char*) grid_data, tot_size*sizeof(float));

}

/* ---------------------------------------------------------------------- */

FixVolumetricRestraint::~FixVolumetricRestraint() {
  delete[] grid_data;
}

/* ---------------------------------------------------------------------- */

int FixVolumetricRestraint::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

double FixVolumetricRestraint::memory_usage()
{
  int tot_size =  grid.nx * grid.ny * grid.nz * 4;
  return tot_size * sizeof(float);
}

/* ---------------------------------------------------------------------- */

void FixVolumetricRestraint::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
  else {
    int nlevels_respa = ((Respa *) update->integrate)->nlevels;
    for (int ilevel = 0; ilevel < nlevels_respa; ilevel++) {
      ((Respa *) update->integrate)->copy_flevel_f(ilevel);
      post_force_respa(vflag,ilevel,0);
      ((Respa *) update->integrate)->copy_f_flevel(ilevel);
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixVolumetricRestraint::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixVolumetricRestraint::post_force(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double ix_d[3];
  int ix[3];

  etotal = 0;
  for (int i = 0; i < nlocal; i++) {
    ftotal[i] = 0;
    if (mask[i] & groupbit) {

      // get the cell indices
      ix_d[0] = ( x[i][0] - grid.x0 ) / grid_step[0];
      ix_d[1] = ( x[i][1] - grid.y0 ) / grid_step[1];
      ix_d[2] = ( x[i][2] - grid.z0 ) / grid_step[2];
      ix[0] = ( ix_d[0] >= 0 && ix_d[0] < grid.nx ) ? int(ix_d[0]) : -1;
      ix[1] = ( ix_d[1] >= 0 && ix_d[1] < grid.ny ) ? int(ix_d[1]) : -1;
      ix[2] = ( ix_d[2] >= 0 && ix_d[2] < grid.nz ) ? int(ix_d[2]) : -1;

      // if outside the grid, no force
      if (ix[0] < 0 || ix[1] < 0 || ix[2] < 0)
        continue;

      // get the position in the data array
      int icell = ( (ix[0] * grid.ny + ix[1]) * grid.nz + ix[2] ) * 4;

      etotal += grid_data[icell];

      f[i][0] += k * grid_data[icell + 1];
      f[i][1] += k * grid_data[icell + 2];
      f[i][2] += k * grid_data[icell + 3];

      ftotal[i] = norm3(f[i]);

    }

  }

}

/* ---------------------------------------------------------------------- */

void FixVolumetricRestraint::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixVolumetricRestraint::min_post_force(int vflag)
{
  post_force(vflag);
}

double FixVolumetricRestraint::compute_scalar()
{

  return etotal;
}

