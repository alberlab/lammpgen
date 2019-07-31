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

#include <cmath>
#include <cstdlib>
#include "bond_harmonic_minimum_upper_bound.h"
#include "atom.h"
#include "neighbor.h"
#include "domain.h"
#include "comm.h"
#include "force.h"
#include "memory.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

BondHarmonicMinimumUpperBound::BondHarmonicMinimumUpperBound(LAMMPS *lmp) : Bond(lmp) {}

/* ---------------------------------------------------------------------- */

BondHarmonicMinimumUpperBound::~BondHarmonicMinimumUpperBound()
{
  if (allocated && !copymode) {
    memory->destroy(setflag);
    memory->destroy(k);
    memory->destroy(maxforce);
    memory->destroy(r0);
    memory->destroy(grps);
    memory->destroy(n_grp);
    memory->destroy(r0sq);
  }
}

/* ---------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::compute(int eflag, int vflag)
{
  if(!groups_are_initalized) create_groups();
  int i1,i2,n,type;
  double ebond,fbond;
  double rsq,r,dr,rk;
  double delx[4], dely[4], delz[4];

  ebond = 0.0;
  if (eflag || vflag) ev_setup(eflag,vflag);
  else evflag = 0;

  double **x = atom->x;
  double **f = atom->f;
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  int nlocal = atom->nlocal;
  int newton_bond = force->newton_bond;

  int nbondtypes = atom->nbondtypes;

  for (int g = 1; g <= nbondtypes; g++){
    int mini = -1;
    double mind = 1000000000000;
    for (int m = 0; m < n_grp[g]; ++m)
    {
      n = grps[g][m];
      i1 = bondlist[n][0];
      i2 = bondlist[n][1];
      
      delx[m] = x[i1][0] - x[i2][0];
      dely[m] = x[i1][1] - x[i2][1];
      delz[m] = x[i1][2] - x[i2][2];
      rsq = delx[m]*delx[m] + dely[m]*dely[m] + delz[m]*delz[m];
      if(rsq < mind)
      {
        mind = rsq;
        mini = m; 
      }
    }
    
    if(mind < r0sq[g]) 
      continue;
    
    i1 = bondlist[grps[g][mini]][0];
    i2 = bondlist[grps[g][mini]][1];

    r = sqrt(mind);

    dr = r - r0[g];
    rk = k[g] * dr;

    // force & energy

    if (r > 0.0) fbond = -2.0*rk/r;
    else fbond = 0.0;
    if (fabs(fbond)>maxforce[g]) fbond = maxforce[g]*(fbond/fabs(fbond)); 
    

    if (eflag) ebond = rk*dr;

    // apply force to each of 2 atoms

    if (newton_bond || i1 < nlocal) {
      f[i1][0] += delx[mini]*fbond;
      f[i1][1] += dely[mini]*fbond;
      f[i1][2] += delz[mini]*fbond;
    }

    if (newton_bond || i2 < nlocal) {
      f[i2][0] -= delx[mini]*fbond;
      f[i2][1] -= dely[mini]*fbond;
      f[i2][2] -= delz[mini]*fbond;
    }

    if (evflag) ev_tally(i1,i2,nlocal,newton_bond,ebond,fbond,delx[mini],dely[mini],delz[mini]);
    
  }

}

/* ---------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::allocate()
{
  allocated = 1;
  groups_are_initalized = false;
  int n = atom->nbondtypes;

  memory->create(k, n+1, "bond:k");
  memory->create(r0, n+1, "bond:r0");
  memory->create(maxforce, n+1, "bond:maxforce");
  memory->create(r0sq, n+1, "bond:r0sq");
  memory->create(grps, n+1, 4, "bond:grps");
  memory->create(n_grp, n+1, "bond:n_grp");

  memory->create(setflag,n+1,"bond:setflag");
  for (int i = 1; i <= n; i++) setflag[i] = 0;
}

/* ----------------------------------------------------------------------
   set coeffs for one or more types
------------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::coeff(int narg, char **arg)
{
  if (narg != 4) error->all(FLERR,"Incorrect args for bond coefficients");
  if (!allocated) allocate();

  int ilo,ihi;
  force->bounds(FLERR, arg[0],atom->nbondtypes,ilo,ihi);

  double k_one = force->numeric(FLERR,arg[1]);
  double r0_one = force->numeric(FLERR,arg[2]);
  double max_force = force->numeric(FLERR,arg[3]);

  int count = 0;
  for (int i = ilo; i <= ihi; i++) {
    k[i] = k_one;
    r0[i] = r0_one;
    maxforce[i] = max_force;
    setflag[i] = 1;
    count++;
  }

  if (count == 0) error->all(FLERR,"Incorrect args for bond coefficients");
}

/* ----------------------------------------------------------------------
   return an equilbrium bond length
------------------------------------------------------------------------- */

double BondHarmonicMinimumUpperBound::equilibrium_distance(int i)
{
  return r0[i];
}

/* ----------------------------------------------------------------------
   proc 0 writes out coeffs to restart file
------------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::write_restart(FILE *fp)
{
  fwrite(&k[1],sizeof(double),atom->nbondtypes,fp);
  fwrite(&r0[1],sizeof(double),atom->nbondtypes,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads coeffs from restart file, bcasts them
------------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::read_restart(FILE *fp)
{
  allocate();

  if (comm->me == 0) {
    fread(&k[1],sizeof(double),atom->nbondtypes,fp);
    fread(&r0[1],sizeof(double),atom->nbondtypes,fp);
  }
  MPI_Bcast(&k[1],atom->nbondtypes,MPI_DOUBLE,0,world);
  MPI_Bcast(&r0[1],atom->nbondtypes,MPI_DOUBLE,0,world);

  for (int i = 1; i <= atom->nbondtypes; i++) setflag[i] = 1;
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->nbondtypes; i++)
    fprintf(fp,"%d %g %g\n",i,k[i],r0[i]);
}

/* ---------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   check if all coeffs are set
------------------------------------------------------------------------- */

void BondHarmonicMinimumUpperBound::init()
{
  if (!allocated && atom->nbondtypes)
    error->all(FLERR,"Bond coeffs are not set");
  for (int i = 1; i <= atom->nbondtypes; i++)
    if (setflag[i] == 0) error->all(FLERR,"All bond coeffs are not set");
  init_style();
}



double BondHarmonicMinimumUpperBound::single(int type, double rsq, int i, int j,
                        double &fforce)
{
  double r = sqrt(rsq);
  double dr = r - r0[type];
  double rk = k[type] * dr;
  fforce = 0;
  if (r > sqrt(r0[type])) fforce = -2.0*rk/r;
  return rk*dr;
}

void BondHarmonicMinimumUpperBound::init_style()
{
  if (force->special_lj[1] != 0.0 || force->special_lj[2] != 1.0 ||
      force->special_lj[3] != 1.0) {
    if (comm->me == 0)
      error->warning(FLERR,"Use special bonds = 0,1,1 with bond style pull_min_harmonic");
  }
  int nbondtypes = atom->nbondtypes;
  for (int g = 1; g <= nbondtypes; g++)
  {
    r0sq[g] = r0[g]*r0[g];
  }

}

void BondHarmonicMinimumUpperBound::create_groups(){
  int nbondtypes = atom->nbondtypes;
  for (int n = 1; n <= nbondtypes; n++)
  {
    n_grp[n] = 0;
  } 
  int **bondlist = neighbor->bondlist;
  int nbondlist = neighbor->nbondlist;
  
  for (int n = 0; n < nbondlist; n++) 
  {
    int type = bondlist[n][2];
    int c = n_grp[type];
    grps[type][c] = n;
    n_grp[type] ++;
  }
  groups_are_initalized = true;
}
