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
  int      s, cell_id, tot_size;
  double                scaling;


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

  // read first few bytes from file: grid obhect
  //grid_info grid;
  fp.read((char*) (&grid), sizeof(grid_info));

  // use number of voxels to define the mappa
  tot_size = grid.nx * grid.ny * grid.nz * 4;

  // define big array and read in one shot
  mappa = new int[tot_size];
  fp.read((char*) (mappa), tot_size * sizeof(int));


  // scaling factors for volumes: for the lamina DamID it would help to have a slightly smaller map, so that there is un intercapetine in cui le particelle sono libere e,
  // anche se non a contatto fisico con la lamina, assumiamo siano a contatto (questo e' il contact range equivalent)
  //if (((body_idx == 1) && (k > 0)) || ((body_idx == 0) && (k < 0))) { scaling = 1.0/evf;}
  //if (((body_idx == 0) && (k > 0)) || ((body_idx == 1) && (k < 0))) { scaling = evf;}

  // lamina damid with nuclear lamina: need a somehow smaller volume
  if ((grid.body_idx == 0) && (k < 0)) {scaling = evf / 0.99;} 
  if ((grid.body_idx == 0) && (k > 0)) {scaling = evf * 0.99;}

  if ((grid.body_idx == 1) && (k > 0)) {scaling = evf / 0.99;}
  if ((grid.body_idx == 1) && (k < 0)) {scaling = evf / 0.99;}
  
  grid.x0  =  grid.x0 * scaling;
  grid.y0  =  grid.y0 * scaling;
  grid.z0  =  grid.z0 * scaling;

  grid.x1  =  grid.x1 * scaling;
  grid.y1  =  grid.y1 * scaling;
  grid.z1  =  grid.z1 * scaling;

  grid.dx  =  grid.dx * scaling;
  grid.dy  =  grid.dy * scaling;
  grid.dz  =  grid.dz * scaling;
  
  /*myfile << scaling << " " << evf << "\n";

  myfile << grid.body_idx << " " << scaling << " " << k << "\n";	
  myfile << grid.nx << " " << grid.ny << " " << grid.nz << "\n";
  myfile << grid.x0 << " " << grid.y0 << " " << grid.z0 << "\n";
  myfile << grid.x1 << " " << grid.y1 << " " << grid.z1 << "\n";
  myfile << grid.dx << " " << grid.dy << " " << grid.dz << "\n";
 
  cell_id = 0 * (grid.ny * grid.nz * 4) + 0 * (grid.nz * 4) + 0 * (4);

  // print to file mappa[i,j,k,:]
  for (s = 0; s < 4; s++) {
        myfile << mappa[cell_id + s] << " ";
    }
  myfile << "\n";

  cell_id = 20 * (grid.ny * grid.nz * 4) + 40 * (grid.nz * 4) + 10 * (4);
 
  // print to file mappa[i,j,k,:]
  for (s = 0; s < 4; s++) {
        myfile << mappa[cell_id + s] << " ";
    }
  myfile << "\n";

  myfile.close(); */

  // this is needed to allocate memory space for the pointers declared in the .h file in the "private" section
  memory->create(this->ftotal, atom->nmax, "FixVolumetricRestraint:ftotal");
  
  etotal = 0;
}



/* ---- deallocate what was created in the preamble here in the file ---------------------------------------------------------- */

FixVolumetricRestraint::~FixVolumetricRestraint() {
  memory->destroy(ftotal);
  delete[] mappa;
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
  int tot_size =  grid.nx * grid.ny * grid.nz * 5;
  return tot_size * sizeof(int);
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
  int         indice;

  //std::ofstream myfile;
  //myfile.open("debug_forces.txt");

  etotal = 0;

  for (int i = 0; i < nlocal; i++) {

    ftotal[i] = 0;
    if (mask[i] & groupbit) {

      //myfile << i << "\n";
      //myfile << " particle coord # " << x[i][0] << " " << x[i][1] << "  " << x[i][2] << "\n";
 
      // get the voxel indices for each particle, by rescaling the coordinates
      ix_d[0] = ( x[i][0] - grid.x1 ) / grid.dx;
      ix_d[1] = ( x[i][1] - grid.y1 ) / grid.dy;
      ix_d[2] = ( x[i][2] - grid.z1 ) / grid.dz;

      //myfile << x[i][0] << " " << grid.nx << "  " << grid.ny << " " << grid.nz << "\n";

      //myfile << " particle voxel # " << ix_d[0] << " " << ix_d[1] << " " << ix_d[2] << "\n";
      //myfile << i << " PRE-FIX " << f[i][0] << " " << f[i][1] << " " << f[i][2] << "\n";

      // account for the possibility of particles floating outside of the density map range:
      // if they are, then voxel indexes are negative to immediately recognize them
      ix[0] = ( ix_d[0] >= 0 && ix_d[0] < grid.nx ) ? int(ix_d[0]) : -1;
      ix[1] = ( ix_d[1] >= 0 && ix_d[1] < grid.ny ) ? int(ix_d[1]) : -1;
      ix[2] = ( ix_d[2] >= 0 && ix_d[2] < grid.nz ) ? int(ix_d[2]) : -1;

      //myfile << "particle # " << i << " voxel = " << ix[0] << " " << ix[1] << " " << ix[2] << "\n";
      //myfile << "particle # " << i << " " << mappa[ix[0]][ix[1]][ix[2]] << "\n";
	
      // volume confinement (exerted by the nuclear lamina, for example)
      if (grid.body_idx == 0) 
      {
 
         // ----------------------------------------------------------------------------------
      	 // if particle is outside grid, then apply a radial attraction
      	 // they are so far away from the envelope that a radial attraction will do just fine
      	 //
      	 // OTHERWISE, apply an inward force pointing to each voxel's EDT-mapped image (on the lamina)
      	 // ----------------------------------------------------------------------------------
      	 
      	 // if particle is inside grid (all coordinates are positive), then we have a force given by the distance of its EDT image
         if (ix[0] >= 0 && ix[1] >= 0 && ix[2] >= 0) {

              // find index
              indice = ix[0] * (grid.ny * grid.nz * 4) + ix[1] * (grid.nz * 4) + ix[2] * 4;

              // if the particle is outside the nuclear volume and we have k >0 (VOLUME) OR inside the volume and we have k<0 (DAMID):
              if ( ((mappa[indice + 3] == 0) && (k > 0)) || ((mappa[indice + 3] != 0) && (k < 0)))
                {
			// find the EDT image of the voxel the locus is currently in

              		edt[0] = mappa[indice + 0];
              		edt[1] = mappa[indice + 1];
              		edt[2] = mappa[indice + 2];

	      		norma = sqrt(  grid.dx * grid.dx * (edt[0] - ix[0]) * (edt[0] - ix[0]) +
              		               grid.dy * grid.dy * (edt[1] - ix[1]) * (edt[1] - ix[1]) +
              		               grid.dz * grid.dz * (edt[2] - ix[2]) * (edt[2] - ix[2]) 
              		             );

              		ftotal[i] = k * k * norma; 
              		etotal    += k * k * norma * norma; 

                      	f[i][0] += k * k * (edt[0] - ix[0])/(grid.nx); // * grid_step[0] * (edt_image[0] - ix[0])/norma;
                      	f[i][1] += k * k * (edt[1] - ix[1])/(grid.ny); // * grid_step[1] * (edt_image[1] - ix[1])/norma;
                      	f[i][2] += k * k * (edt[2] - ix[2])/(grid.nz); // * grid_step[2] * (edt_image[2] - ix[2])/norma;

		 }
              }
         
         // if particle is partially outside of the unitary cell, then a unitary radial attractive force will do the trick
         else if (k > 0) {
              // if volume, otherwise no force
              norma = sqrt((x[i][0] - grid.x0) * (x[i][0] - grid.x0) +
                           (x[i][1] - grid.y0) * (x[i][1] - grid.y0) +
                           (x[i][2] - grid.z0) * (x[i][2] - grid.z0)
                          );

              ftotal[i] = norm3(f[i]);
              etotal    += k / 2;

              f[i][0] += k * (grid.x0 - x[i][0]) / norma;
              f[i][1] += k * (grid.y0 - x[i][1]) / norma;
              f[i][2] += k * (grid.z0 - x[i][2]) / norma;

	       }                                                
      }

      // IF WE ARE TALKING ABOUT NUCLEAR BODY
      if (grid.body_idx == 1)
      {

          // if inside the volume grid and inside the reference volume, apply expulsion/excluded force
	  if (ix[0] >= 0 && ix[1] >=0 && ix[2] >=0)
          {
              // find index
              indice = ix[0] * (grid.ny * grid.nz * 4) + ix[1] * (grid.nz * 4) + ix[2] * 4;
               
              // ... outside nucleolar volume but to be attracted back to surface (lamina) or inside volume to be expelled
              if (((mappa[indice + 3] == 1) && (k >0)) || ((mappa[indice + 3] == 0) && (k <0)))
	      {

              		edt[0] = mappa[indice + 0];
                	edt[1] = mappa[indice + 1];
                	edt[2] = mappa[indice + 2];
              
                	norma = sqrt( grid.dx * grid.dx * (edt[0] - ix[0]) * (edt[0] - ix[0]) +
                              	      grid.dy * grid.dy * (edt[1] - ix[1]) * (edt[1] - ix[1]) +
                                      grid.dz * grid.dz * (edt[2] - ix[2]) * (edt[2] - ix[2])
                                    );
              
                	ftotal[i] = k * k * norma;
                	etotal    += k * k * norma * norma;

                   	f[i][0] += k * k * (edt[0] - ix[0])/grid.nx; //]* grid_step[0] * (edt_image[0] - ix[0])/norma;
                   	f[i][1] += k * k * (edt[1] - ix[1])/grid.ny; //grid_step[1] * (edt_image[1] - ix[1])/norma;
                   	f[i][2] += k * k * (edt[2] - ix[2])/grid.nz; //grid_step[2] * (edt_image[2] - ix[2])/norma;
	      }
	  }
          else if (k < 0)
                {
                norma =  sqrt((x[i][0] - grid.x0) * (x[i][0] - grid.x0) +
                              (x[i][1] - grid.y0) * (x[i][1] - grid.y0) +
                              (x[i][2] - grid.z0) * (x[i][2] - grid.z0)
                          );

              ftotal[i] = norm3(f[i]);
              etotal    += k * k / 2;

              f[i][0] += k * k * (grid.x0 - x[i][0]) / norma;
              f[i][1] += k * k * (grid.y0 - x[i][1]) / norma;
              f[i][2] += k * k * (grid.z0 - x[i][2]) / norma;
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
