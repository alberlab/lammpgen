/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifdef FIX_CLASS

FixStyle(volumetricrestraint,FixVolumetricRestraint)

#else

#ifndef LMP_FIX_VOLUMETRICRESTRAINT_H
#define LMP_FIX_VOLUMETRICRESTRAINT_H

#include "fix.h"

namespace LAMMPS_NS {


class FixVolumetricRestraint : public Fix {
 public:
  FixVolumetricRestraint(class LAMMPS *, int, char **);
  virtual ~FixVolumetricRestraint();
  virtual double memory_usage();
  int setmask();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);
  double compute_scalar();
  double compute_vector(int);

 private:
  double      kspring, evf;  // elastic constant, factor to scale up/down any nuclear/nucleolar envelope
  double    origin[3];   // [x_min, y_min, z_min] origin, lowest, leftmost vertex for the volume map
  double grid_step[3]; // [dx, dy, dz] grid spacing along the three dimensions
  double    center[3];    // geometric center of the nuclear body, pozzo del campo centrale di forza
  int      n_voxel[3];    // number of voxels along the three dimensions
  int        body_idx;     // index to discriminate between nuclear envelope and/ord nuclear bodies
 
  double mappa[220][150][45];   // per il momento e' inizializzata in eccesso...cioe' e' piu' grande di quello che serve
  double etotal;
  double* ftotal;
};

}

#endif
#endif

/* ERROR/WARNING messages:

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

*/
