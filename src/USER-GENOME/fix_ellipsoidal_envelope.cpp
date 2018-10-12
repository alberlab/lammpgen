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

#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "fix_ellipsoidal_envelope.h"
#include "atom.h"
#include "memory.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

double max(double x, double y) {return x>y? x: y;}
double normsq3(double* x){ return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2]); }

FixEllipsoidalEnvelope::FixEllipsoidalEnvelope(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 1;
  dynamic_group_allow = 1;

  if (narg != 7) error->all(FLERR,"Illegal fix ellipsodal envelope command [args: a, b, c, k]");
  a = force->numeric(FLERR,arg[3]);
  b = force->numeric(FLERR,arg[4]);
  c = force->numeric(FLERR,arg[5]);
  kspring = force->numeric(FLERR,arg[6]);

  a2 = a*a;
  b2 = b*b;
  c2 = c*c;

  int dflag = 0;
  int index = atom->find_custom("radius", dflag);
  radius = atom->dvector[index];

  memory->create(this->v2r, atom->nmax, 3, "FixEllipsoidalEnvelope:v2r");
  memory->create(this->ftotal, atom->nmax, "FixEllipsoidalEnvelope:ftotal");

  for (int i = 0; i < atom->nmax; ++i){

    v2r[i][0] = (a - radius[i])*(a - radius[i]);
    v2r[i][1] = (b - radius[i])*(b - radius[i]);
    v2r[i][2] = (c - radius[i])*(c - radius[i]);

  }

  etotal = 0;

  atom->add_callback(0);
}

/* ---------------------------------------------------------------------- */

FixEllipsoidalEnvelope::~FixEllipsoidalEnvelope() {
  memory->destroy(v2r);
  memory->destroy(ftotal);
  atom->delete_callback(id, 0);
}

/* ---------------------------------------------------------------------- */

int FixEllipsoidalEnvelope::setmask() {
  int mask = 0;
  mask |= POST_FORCE;
  mask |= THERMO_ENERGY;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

double FixEllipsoidalEnvelope::memory_usage()
{
  int nmax = atom->nmax;
  double bytes = 0.0;
  bytes += nmax * 3 * sizeof(double);
  return bytes;
}

/* ---------------------------------------------------------------------- */

void FixEllipsoidalEnvelope::grow_arrays(int nmax)
{
  memory->grow(this->v2r, nmax, 3, "FixEllipsoidalEnvelope:v2r");
}

void FixEllipsoidalEnvelope::copy_arrays(int i, int j, int delflag)
{
  memcpy(this->v2r[j], this->v2r[i], sizeof(double) * 3);
}

void FixEllipsoidalEnvelope::set_arrays(int i)
{
  memset(this->v2r[i], 0, sizeof(double) * 3);
}

int FixEllipsoidalEnvelope::pack_exchange(int i, double *buf)
{
  int m = 0;
  buf[m++] = v2r[i][0];
  buf[m++] = v2r[i][1];
  buf[m++] = v2r[i][2];
  return m;
}

int FixEllipsoidalEnvelope::unpack_exchange(int nlocal, double *buf)
{
  int m = 0;
  v2r[nlocal][0] = buf[m++];
  v2r[nlocal][1] = buf[m++];
  v2r[nlocal][2] = buf[m++];
  return m;
}

/* ---------------------------------------------------------------------- */

void FixEllipsoidalEnvelope::setup(int vflag)
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

void FixEllipsoidalEnvelope::min_setup(int vflag)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEllipsoidalEnvelope::post_force(int vflag)
{
  double **f = atom->f;
  double **x = atom->x;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  double dot;
  double v[3], x2[3], k2;
  double t;
  double vnormsq;

  etotal = 0;
  for (int i = 0; i < nlocal; i++) {
    ftotal[i] = 0;
    if (mask[i] & groupbit) {

      x2[0] = x[i][0]*x[i][0];
      x2[1] = x[i][1]*x[i][1];
      x2[2] = x[i][2]*x[i][2];

      k2 = x2[0]/v2r[i][0] + x2[1]/v2r[i][1] + x2[2]/v2r[i][2];

      if ( ( k2 > 1 && kspring > 0 ) || ( k2 < 1 && kspring < 0 ) )  {

        // k2 > 1 means the point is outside the ellipse

        // The gradient is -(x/a^2, y/b^2, z/c^2)
        // The minimal distance to the ellipsodal is difficult to solve
        // However we can approximate the distance using a cooresponding point
        // on the surface scaled by 1/sqrt(x^2/a^2 + y^2/b^2 + z^2/c^2) = 1/sqrt(k2)
        // t = (1-1/sqrt(k2))*|(x,y,z)|


        v[0] = -x[i][0]/v2r[i][0];
        v[1] = -x[i][1]/v2r[i][1];
        v[2] = -x[i][2]/v2r[i][2];
        t = (1 - 1/sqrt(k2))*sqrt(x2[0]+x2[1]+x2[2]);

        vnormsq = normsq3(v);
        ftotal[i] = t * sqrt(vnormsq) * kspring; //0.5 * t * t * vnormsq * kspring; //
        etotal += ftotal[i];
        f[i][0] += t * v[0] * fabs(kspring);
        f[i][1] += t * v[1] * fabs(kspring);
        f[i][2] += t * v[2] * fabs(kspring);


      }

    }

  }

}

/* ---------------------------------------------------------------------- */

void FixEllipsoidalEnvelope::post_force_respa(int vflag, int ilevel, int iloop)
{
  post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEllipsoidalEnvelope::min_post_force(int vflag)
{
  post_force(vflag);
}

double FixEllipsoidalEnvelope::compute_scalar()
{

  return etotal;
}


double FixEllipsoidalEnvelope::compute_vector(int n)
{
  return ftotal[n];
}
