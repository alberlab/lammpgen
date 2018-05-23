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

FixStyle(ellipsoidalenvelope,FixEllipsoidalEnvelope)

#else

#ifndef LMP_FIX_ELLIPSOIDALENVELOPE_H
#define LMP_FIX_ELLIPSOIDALENVELOPE_H

#include "fix.h"

namespace LAMMPS_NS {

class FixEllipsoidalEnvelope : public Fix {
 public:
  FixEllipsoidalEnvelope(class LAMMPS *, int, char **);
  virtual ~FixEllipsoidalEnvelope();
  virtual double memory_usage();
  virtual void grow_arrays(int);
  virtual void copy_arrays(int, int, int);
  virtual void set_arrays(int);
  virtual int pack_exchange(int, double*);
  virtual int unpack_exchange(int, double*);
  int setmask();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void post_force_respa(int, int, int);
  void min_post_force(int);

 private:
  double a, b, c, a2, b2, c2, kspring;
  double* radius;
  double** v2r;
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
