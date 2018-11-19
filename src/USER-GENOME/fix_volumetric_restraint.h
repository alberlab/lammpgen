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

struct grid_info {
  int nx, ny, nz;
  float x0, y0, z0;
  float x1, y1, z1;
};

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

 private:
  double k;
  float* grid_data;
  grid_info grid;
  double etotal;
  double* ftotal;
  double grid_step[3];
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
