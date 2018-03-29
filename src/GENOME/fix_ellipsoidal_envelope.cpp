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

  max_axis = max( max(a, b), c );

  a2 = a*a;
  b2 = b*b;
  c2 = c*c;

  a4 = a2*a2;
  b4 = b2*b2;
  c4 = c2*c2;

  a6 = a4*a2;
  b6 = b4*b2;
  c6 = c4*c2;

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
  double v[3], xs[3], x2[3], k2, k4, k6;
  double t, d2, nrm;
  
  for (int i = 0; i < nlocal; i++)

    if (mask[i] & groupbit) {
      
      x2[0] = x[i][0]*x[i][0];
      x2[1] = x[i][1]*x[i][1];
      x2[2] = x[i][2]*x[i][2];

      k2 = x2[0]/a2 + x2[1]/b2 + x2[2]/c2;
      k4 = x2[0]/a4 + x2[1]/b4 + x2[2]/c4;
      k6 = x2[0]/a6 + x2[1]/b6 + x2[2]/c6;

      if ( k2 > 1 )  {
        
        // k2 > 1 means the point is outside the ellipse

        d2 = k4*k4 - k6*(k2 - 1);
        
        if ( d2 < 0 ) {

          // if the point is far away from the ellipse, it may happen that
          // the gradient direction does not intersect the ellipse itself.
          // In that case we push towards the center

          v[0] = -x[i][0];
          v[1] = -x[i][1];
          v[2] = -x[i][2];
          nrm = sqrt( x[i][0]*x[i][0] + x[i][1]*x[i][1] + x[i][2]*x[i][2] );
          t = 1.0 - max_axis / nrm;

        } else {

          // we are close enough to the ellipse to have a solution, we 
          // follow the gradient

          v[0] = -x[i][0]/a2;
          v[1] = -x[i][1]/b2;
          v[2] = -x[i][2]/c2;
          t = ( k4 - sqrt( d2 ) ) / k6;   
        
        }

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
