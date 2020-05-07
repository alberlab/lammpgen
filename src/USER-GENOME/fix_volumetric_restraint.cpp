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
#include <fstream>
#include "fix_volumetric_restraint.h"
#include "atom.h"
#include "memory.h"
#include "update.h"
#include "respa.h"
#include "error.h"
#include "force.h"

// this is necessary to use std:cout to print to terminal
#include <iostream>

using namespace LAMMPS_NS;
using namespace FixConst;

/* ---------------------------------------------------------------------- */
// auxiliary functions, max and norm 
// use 'namespace' to limit definitions to current file only
namespace {

double max(double x, double y) {return x > y ? y :x;};
double norm3(double* x){return sqrt(x[0]*x[0] + x[1]*x[1] + x[2]* x[2]);};

}

/* ---------------------------------------------------------------------- */
FixVolumetricRestraint::FixVolumetricRestraint(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{

  // preamble
  int tmp1, tmp2, tmp3, tmp4;
  int count = 0;
  double new_center[3];

  std::ofstream myfile;
  //myfile.open ("debug.txt");
 
  scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 1;
  dynamic_group_allow = 1;

  if (narg != 6) error->all(FLERR,"Illegal volumetric restraint command [args: density_filename env_factor k]");
  
  // define parameter k which I expect to be spring constant
  kspring = force->numeric(FLERR,arg[5]);

  // define parameter evf used in optimization to slightly scale up or down NE
  evf     = force->numeric(FLERR,arg[4]);

  /*------- read in file whose name is given in "arg[3]", STRUCT object ---- */
  std::ifstream fp(arg[3], std::ios::binary);

  fp >>                                           body_idx;
  fp >>     n_voxel[0] >>     n_voxel[1] >>     n_voxel[2]; 
  fp >>      center[0] >>      center[1] >>      center[2];
  fp >>      origin[0] >>      origin[1] >>      origin[2];  // origin of the map
  fp >>   grid_step[0] >>   grid_step[1] >>   grid_step[2];  // grid step

  //myfile << n_voxel[0] << " " << n_voxel[1] << " " << n_voxel[2] << "\n";
  //myfile <<  origin[0] << " " <<  origin[1] << " " <<  origin[2] << "\n";
  //myfile << grid_step[0] << " " << grid_step[1] << " " << grid_step[2] << "\n";

  /* read in all the entries of the volume map matrix */
  for (int i = 0; i < n_voxel[0] * n_voxel[1] * n_voxel[2]; i++)
      {
	fp >> tmp1 >> tmp2 >> tmp3 >> tmp4;    // one full line at the time, each line is (x, y, z, density)
        mappa[tmp1][tmp2][tmp3] = tmp4;    // reproduce the density map
      } 

  //myfile << grid_step[0] << ' ' << grid_step[1] << ' ' << grid_step[2] << '\n';
  //myfile << origin[0]    << ' ' << origin[1]    << ' ' << origin[2] << '\n';

  //myfile << " scaling I mean = " << evf << '\n';
 
  // if scaling factor different from 1.0...
  if (evf != 1.0) {

 
 	 // scale grid quantities by factor "evf"
  	if (body_idx == 0) {

  		center[0]    = center[0] * evf;
  		center[1]    = center[1] * evf;
  		center[2]    = center[2] * evf;

  		origin[0]    = origin[0] * evf;
  		origin[1]    = origin[1] * evf;
  		origin[2]    = origin[2] * evf;

  		grid_step[0] = grid_step[0] * evf;
  		grid_step[1] = grid_step[1] * evf;
  		grid_step[2] = grid_step[2] * evf;
  	}

  	// if nuclear bodies, it is important that the geometric center is preserved
  	if (body_idx == 1) {

		new_center[0] = center[0]/evf;
                new_center[1] = center[1]/evf;
                new_center[2] = center[2]/evf;

  		origin[0]    = origin[0] / evf;
  		origin[1]    = origin[1] / evf;
  		origin[2]    = origin[2] / evf;

  		grid_step[0] = grid_step[0] / evf;
  		grid_step[1] = grid_step[1] / evf;
  		grid_step[2] = grid_step[2] / evf;

		// translate scaled grid s.t. center is the (non-scaled) nucleolar center
  		origin[0] = origin[0] + (center[0] - new_center[0]);
                origin[1] = origin[1] + (center[1] - new_center[1]);
                origin[2] = origin[2] + (center[2] - new_center[2]);
	}
     }

  /*// check something is sto cristo di lettura da file
   for (int i = 0; i < n_voxel[0]; ++i)
     {
         for (int j=0; j < n_voxel[1]; ++j)
         {
            for (int k=0; k < n_voxel[2]; ++k)
            {
                if (mappa[i][j][k] == 0)
                {
                    // print stuff
                    myfile << i << " " << j << " " << k << " " << mappa[i][j][k] << "\n";
                    count = count + 1;
                }
            }
         }
     }
  myfile  << "number of 0 entries = " << count << "\n";
  myyfile.close(); */
 myfile.close();

  // this is needed to allocate memory space for the pointers declared in the .h file in the "private" section
  memory->create(this->ftotal, atom->nmax, "FixVolumetricRestraint:ftotal");
  
  etotal = 0;
}



/* ---- deallocate what was created in the preamble here in the file ---------------------------------------------------------- */

FixVolumetricRestraint::~FixVolumetricRestraint() {
  memory->destroy(ftotal);
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
  int tot_size =  n_voxel[0] * n_voxel[1] * n_voxel[2] * 4;
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

  // auxiliary vectors
  double ix_d[3];
  int    ix[3];
  //double temp[3] = {0,0,0};
  //int count = 0;

  double norma;      // to save the norm of an array computed manually

  //std::ofstream myfile;
  //myfile.open("debug_forces.txt");

  etotal = 0;

  for (int i = 0; i < nlocal; i++) {
    ftotal[i] = 0;
    if (mask[i] & groupbit) {

      //myfile << i << "\n";
      //myfile << " particle coord # " << x[i][0] << " " << x[i][1] << "  " << x[i][2] << "\n";
 
      // get the voxel indices for each particle, by rescaling the coordinates
      ix_d[0] = ( x[i][0] - origin[0] ) / grid_step[0];
      ix_d[1] = ( x[i][1] - origin[1] ) / grid_step[1];
      ix_d[2] = ( x[i][2] - origin[2] ) / grid_step[2];

      //myfile << " particle voxel # " << ix_d[0] << " " << ix_d[1] << " " << ix_d[2] << "\n";
      //myfile << i << " PRE-FIX " << f[i][0] << " " << f[i][1] << " " << f[i][2] << "\n";

      // account for the possibility of particles floating outside of the density map range
      ix[0] = ( ix_d[0] >= 0 && ix_d[0] < n_voxel[0] ) ? int(ix_d[0]) : -1;
      ix[1] = ( ix_d[1] >= 0 && ix_d[1] < n_voxel[1] ) ? int(ix_d[1]) : -1;
      ix[2] = ( ix_d[2] >= 0 && ix_d[2] < n_voxel[2] ) ? int(ix_d[2]) : -1;

      //myfile << "particle # " << i << " voxel = " << ix[0] << " " << ix[1] << " " << ix[2] << "\n";
      //myfile << "particle # " << i << " " << mappa[ix[0]][ix[1]][ix[2]] << "\n";
	
      //if ((ix[0] < 0 || ix[1] < 0 || ix[2] < 0))
      //	{	myfile << i << "\n";
      //	} 

      // IF WE ARE TALKING ABOUT NUCLEAR ENVELOPE
      if (body_idx == 0)
	{
 
      	 // if particle is outside of density map range OR inside, but outside of the volume, then apply restraint
         if ((ix[0] < 0 || ix[1] < 0 || ix[2] < 0) || (mappa[ix[0]][ix[1]][ix[2]] == 1))
         	{

	      norma = sqrt((x[i][0] - center[0]) * (x[i][0] - center[0]) + 
                           (x[i][1] - center[1]) * (x[i][1] - center[1]) + 
                           (x[i][2] - center[2]) * (x[i][2] - center[2])
                          );

              ftotal[i] = norm3(f[i]);
              etotal    += kspring * norma;

              f[i][0] += kspring * (center[0] - x[i][0]) / norma;
              f[i][1] += kspring * (center[1] - x[i][1]) / norma;
              f[i][2] += kspring * (center[2] - x[i][2]) / norma;

	        }
         }

      // IF WE ARE TALKING ABOUT NUCLEAR BODY
      if (body_idx == 1)
	 {
	  if ((ix[0] >= 0 && ix[1] >=0 && ix[2] >=0) && (mappa[ix[0]][ix[1]][ix[2]] == 1))
		{
	      norma = sqrt((x[i][0] - center[0]) * (x[i][0] - center[0]) +
                           (x[i][1] - center[1]) * (x[i][1] - center[1]) +
                           (x[i][2] - center[2]) * (x[i][2] - center[2])
                          );

              ftotal[i] =  norm3(f[i]);
              etotal    += kspring * norma;

              f[i][0] += kspring * (x[i][0] - center[0]) / norma;
              f[i][1] += kspring * (x[i][1] - center[1]) / norma;
              f[i][2] += kspring * (x[i][2] - center[2]) / norma;
		
		}
	 }
    
        //myfile << i << " POST-FIX " << f[i][0] << " " << f[i][1] << " " << f[i][2] << "\n\n";
 
    }

  }

 // myfile << "\n\n\n";
 // myfile << " number of particles outside = " << count << "\n";

 // myfile.close();
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

/* ------- compute the total energy associated with the force contribution ----- */
double FixVolumetricRestraint::compute_scalar()
{

  return etotal;
}

/*-------- compute the force contribution ------------------------------ */
double FixVolumetricRestraint::compute_vector(int n)
{
  return ftotal[n];
}
