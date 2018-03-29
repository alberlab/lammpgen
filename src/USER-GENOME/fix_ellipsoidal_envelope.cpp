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
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */

double max(double x, double y) {return x>y? x: y;}

FixEllipsoidalEnvelope::FixEllipsoidalEnvelope(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  dynamic_group_allow = 1;

  if (narg != 6) error->all(FLERR,"Illegal fix planeforce command");
  a = force->numeric(FLERR,arg[3]);
  b = force->numeric(FLERR,arg[4]);
  c = force->numeric(FLERR,arg[5]);

  a2 = a*a;
  b2 = b*b;
  c2 = c*c;

}

/* ---------------------------------------------------------------------- */

int FixEllipsoidalEnvelope::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= POST_FORCE_RESPA;
  mask |= MIN_POST_FORCE;
  return mask;
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
  
  for (int i = 0; i < nlocal; i++)

    if (mask[i] & groupbit) {
      
      x2[0] = x[i][0]*x[i][0];
      x2[1] = x[i][1]*x[i][1];
      x2[2] = x[i][2]*x[i][2];

      k2 = x2[0]/a2 + x2[1]/b2 + x2[2]/c2;
      

      if ( k2 > 1 )  {
        
        // k2 > 1 means the point is outside the ellipse
        
        // The gradient is -(x/a^2, y/b^2, z/c^2)
        // The minimal distance to the ellipsodal is difficult to solve
        // However we can approximate the distance using a cooresponding point
        // on the surface scaled by 1/sqrt(x^2/a^2 + y^2/b^2 + z^2/c^2) = 1/sqrt(k2)
        // t = (1-1/sqrt(k2))*|(x,y,z)|
        

        v[0] = -x[i][0]/a2;
        v[1] = -x[i][1]/b2;
        v[2] = -x[i][2]/c2;
        t = (1 - 1/sqrt(k2))*sqrt(x2[0]+x2[1]+x2[2]);   
        
        f[i][0] += t * v[0];
        f[i][1] += t * v[1];
        f[i][2] += t * v[2];
      
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
