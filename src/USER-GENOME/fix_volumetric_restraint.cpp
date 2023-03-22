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
  int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
  int                          count = 0;
  double                   new_center[3];
  double                         scaling;

  //std::ofstream myfile;
  //myfile.open ("debug.txt");
 
  scalar_flag = 1;
  //vector_flag = 1;
  //size_vector = 1;
  dynamic_group_allow = 1;

  if (narg != 6) error->all(FLERR,"Illegal volumetric restraint command [args: volume_filename env_factor k]");
  
  k = force->numeric(FLERR,arg[5]);

  // define parameter evf used in optimization to slightly scale up or down envelopes/volumes
  evf     = force->numeric(FLERR,arg[4]);

  /*------- read in file whose name is given in "arg[3]", STRUCT object ---- */
  std::ifstream fp(arg[3], std::ios::binary);

  fp >>                                           body_idx;
  fp >>     n_voxel[0] >>     n_voxel[1] >>     n_voxel[2];  // nx, ny, nz (number of pixels in each direction)
  fp >>      center[0] >>      center[1] >>      center[2];  // center of mass of the map
  fp >>      origin[0] >>      origin[1] >>      origin[2];  // origin of the map (bottom, most left pixel edge)
  fp >>   grid_step[0] >>   grid_step[1] >>   grid_step[2];  // step size of a pixel/voxel

  // define the size of the arrays

  //myfile << grid_step[0] << " " << grid_step[1] << " " << grid_step[2] << "\n";

  /* read in all the entries of the volume map matrix */
  for (int i = 0; i < n_voxel[0] * n_voxel[1] * n_voxel[2]; i++)
      {
	fp >> tmp1 >> tmp2 >> tmp3 >> tmp4 >> tmp5 >> tmp6 >> tmp7;    // one full line at the time, each line is (x, y, z, EDT[x], EDT[y], EDT[z], in/out label)
 
        mappa[tmp1][tmp2][tmp3][0] = tmp4;    // reproduce the Euclidean Distance Transform pattern
        mappa[tmp1][tmp2][tmp3][1] = tmp5;
        mappa[tmp1][tmp2][tmp3][2] = tmp6;
        mappa[tmp1][tmp2][tmp3][3] = tmp7;
 
      } 

  // scaling factors for volumes: for the lamina DamID it would help to have a slightly smaller map, so that there is un intercapetine in cui le particelle sono libere e,
  // anche se non a contatto fisico con la lamina, assumiamo siano a contatto (questo e' il contact range equivalent)
  //if (((body_idx == 1) && (k > 0)) || ((body_idx == 0) && (k < 0))) { scaling = 1.0/evf;}
  //if (((body_idx == 0) && (k > 0)) || ((body_idx == 1) && (k < 0))) { scaling = evf;}

  // lamina damid with nuclear lamina: need a somehow smaller volume
  if ((body_idx == 0) && (k < 0)) {scaling = evf * 0.95;}
  if ((body_idx == 0) && (k > 0)) {scaling = evf;}

  if ((body_idx == 1) && (k < 0)) {scaling = evf;}
  if ((body_idx == 1) && (k > 0)) {scaling = evf / 0.95;}
  
  center[0]     =     center[0] * scaling;
  center[1]     =     center[1] * scaling;
  center[2]     =     center[2] * scaling;

  origin[0]     =     origin[0] * scaling;
  origin[1]     =     origin[1] * scaling;
  origin[2]     =     origin[2] * scaling;

  grid_step[0]  =  grid_step[0] * scaling;
  grid_step[1]  =  grid_step[1] * scaling;
  grid_step[2]  =  grid_step[2] * scaling;
  	

  //myfile << center[0] << " " << center[1] << " " << center[2] << "\n";
  //myfile << origin[0] << " " << origin[1] << " " << origin[2] << "\n";
  //myfile << grid_step[0] << " " << grid_step[1] << " " << grid_step[2] << "\n";

  //myfile.close();

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
  int tot_size =  n_voxel[0] * n_voxel[1] * n_voxel[2] * 5;
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
  double     ix_d[3];    // integer indexes for the voxel where the particle sits
  int         edt[3];    // integer indexes for the closest EDT (Euclidean Distance Transform) to the particle voxel [ix_d[0], id_x[1], id_x[2]]
  int          ix[3];    // integer indexes for the voxel where the particle sits (which factors in the possibility of the particle being outside of grid)
  double       norma;    // to save the norm of an array computed manually

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

      // account for the possibility of particles floating outside of the density map range:
      // if they are, then voxel indexes are negative to immediately recognize them
      ix[0] = ( ix_d[0] >= 0 && ix_d[0] < n_voxel[0] ) ? int(ix_d[0]) : -1;
      ix[1] = ( ix_d[1] >= 0 && ix_d[1] < n_voxel[1] ) ? int(ix_d[1]) : -1;
      ix[2] = ( ix_d[2] >= 0 && ix_d[2] < n_voxel[2] ) ? int(ix_d[2]) : -1;

      //myfile << "particle # " << i << " voxel = " << ix[0] << " " << ix[1] << " " << ix[2] << "\n";
      //myfile << "particle # " << i << " " << mappa[ix[0]][ix[1]][ix[2]] << "\n";
	
      // volume confinement (exerted by the nuclear lamina, for example)
      if (body_idx == 0) 
      {
 
         // ----------------------------------------------------------------------------------
      	 // if particle is outside grid, then apply a radial attraction
      	 // they are so far away from the envelope that a radial attraction will do just fine
      	 //
      	 // OTHERWISE, apply an inward force pointing to each voxel's EDT-mapped image (on the lamina)
      	 // ----------------------------------------------------------------------------------
      	 
      	 // if particle is inside grid (all coordinates are positive), then we have a force given by the distance of its EDT image
         if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0) {

              // if the particle is outside the nuclear volume and we have k >0 (VOLUME) OR inside the volume and we have k<0 (DAMID):
              if ( ((mappa[ix[0]][ix[1]][ix[2]][3] == 0) && (k > 0)) || ((mappa[ix[0]][ix[1]][ix[2]][3] != 0) && (k < 0)))
                {
			// find the EDT image of the voxel the locus is currently in

              		edt[0] = mappa[ix[0]][ix[1]][ix[2]][0];
              		edt[1] = mappa[ix[0]][ix[1]][ix[2]][1];
              		edt[2] = mappa[ix[0]][ix[1]][ix[2]][2];

	      		norma = sqrt(  grid_step[0] * grid_step[0] * (edt[0] - ix[0]) * (edt[0] - ix[0]) +
              		               grid_step[1] * grid_step[1] * (edt[1] - ix[1]) * (edt[1] - ix[1]) +
              		               grid_step[2] * grid_step[2] * (edt[2] - ix[2]) * (edt[2] - ix[2]) 
              		             );

              		ftotal[i] = k * k * norma; 
              		etotal    += k * k * norma * norma; 

                      	f[i][0] += k * k * (edt[0] - ix[0])/n_voxel[0]; // * grid_step[0] * (edt_image[0] - ix[0])/norma;
                      	f[i][1] += k * k * (edt[1] - ix[1])/n_voxel[1]; // * grid_step[1] * (edt_image[1] - ix[1])/norma;
                      	f[i][2] += k * k * (edt[2] - ix[2])/n_voxel[2]; // * grid_step[2] * (edt_image[2] - ix[2])/norma;

		 }
              }
         
         // if particle is partially outside of the unitary cell, then a unitary radial attractive force will do the trick
         else if (k > 0) {
              // if volume, otherwise no force
              norma = sqrt((x[i][0] - center[0]) * (x[i][0] - center[0]) +
                           (x[i][1] - center[1]) * (x[i][1] - center[1]) +
                           (x[i][2] - center[2]) * (x[i][2] - center[2])
                          );

              ftotal[i] = norm3(f[i]);
              etotal    += k / 2;

              f[i][0] += k * (center[0] - x[i][0]) / norma;
              f[i][1] += k * (center[1] - x[i][1]) / norma;
              f[i][2] += k * (center[2] - x[i][2]) / norma;

	       }                                                
      }

      // IF WE ARE TALKING ABOUT NUCLEAR BODY
      if (body_idx == 1)
      {

          // if inside the volume grid and inside the reference volume, apply expulsion/excluded force
	  if (ix[0] >= 0 && ix[1] >=0 && ix[2] >=0)
          {
              // ... outside nucleolar volume but to be attracted back to surface (lamina) or inside volume to be expelled
              if (((mappa[ix[0]][ix[1]][ix[2]][3] == 1) && (k >0)) || ((mappa[ix[0]][ix[1]][ix[2]][3] == 0) && (k <0)))
	      {

              		edt[0] = mappa[ix[0]][ix[1]][ix[2]][0];
                	edt[1] = mappa[ix[0]][ix[1]][ix[2]][1];
                	edt[2] = mappa[ix[0]][ix[1]][ix[2]][2];
              
                	norma = sqrt( grid_step[0] * grid_step[0] * (edt[0] - ix[0]) * (edt[0] - ix[0]) +
                              	      grid_step[1] * grid_step[1] * (edt[1] - ix[1]) * (edt[1] - ix[1]) +
                                      grid_step[2] * grid_step[2] * (edt[2] - ix[2]) * (edt[2] - ix[2])
                                    );
              
                	ftotal[i] = k * k * norma;
                	etotal    += k * k * norma * norma;

                   	f[i][0] += k * k * (edt[0] - ix[0])/n_voxel[0]; //]* grid_step[0] * (edt_image[0] - ix[0])/norma;
                   	f[i][1] += k * k * (edt[1] - ix[1])/n_voxel[1]; //grid_step[1] * (edt_image[1] - ix[1])/norma;
                   	f[i][2] += k * k * (edt[2] - ix[2])/n_voxel[2]; //grid_step[2] * (edt_image[2] - ix[2])/norma;
	      }
	  }
          else if (k < 0)
                {
                norma =  sqrt((x[i][0] - center[0]) * (x[i][0] - center[0]) +
                              (x[i][1] - center[1]) * (x[i][1] - center[1]) +
                              (x[i][2] - center[2]) * (x[i][2] - center[2])
                          );

              ftotal[i] = norm3(f[i]);
              etotal    += k / 2;

              f[i][0] += k * k * (center[0] - x[i][0]) / norma;
              f[i][1] += k * k * (center[1] - x[i][1]) / norma;
              f[i][2] += k * k * (center[2] - x[i][2]) / norma;
                }
     }

  }
}

 // myfile << "\n\n\n";
 // myfile << " number of particles outside = " << count << "\n";

  //myfile.close();
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
