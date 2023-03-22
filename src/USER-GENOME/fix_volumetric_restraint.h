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
  double            k, evf;  // elastic constant, factor to scale up/down any nuclear/nucleolar envelope
  double         origin[3];   // [x_min, y_min, z_min] origin, lowest, leftmost vertex for the volume map
  double      grid_step[3]; // [dx, dy, dz] grid spacing (aka, voxel size) along the x, y and z binning
  double         center[3];    // geometric center of the nuclear body, pozzo del campo centrale di forza
  int           n_voxel[3];    // number of voxels along the three dimensions
  int             body_idx;     // index to discriminate between nuclear envelope and/ord nuclear bodies
 
  int               mappa[500][410][120][4];   // per il momento e' inizializzata in eccesso...cioe' e' piu' grande di quello che serve (dynamic stuff)?
					       // dynamic allocation is possible using pointers of pointers etc...but it will be much slower than just 
					       // initializing a larger-than-needed static array.
					       // Guido's version used dynamic allocation: see "fix_volumetric_smart.cpp" in /u/home/b/bonimba/lammps_
					       // spare_code for the code to allocate and de-allocate dynamically, after nvoxel array is read in.
					       // Also, it might help improve performance and data storage if we used a c++ libray to read npz files.
					       // Instead of plain txt files. This is a technical improvement.


  double            etotal;
  double*           ftotal;
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
