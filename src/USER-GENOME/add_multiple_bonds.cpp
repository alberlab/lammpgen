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

#include <cstdlib>
#include <cstring>
#include "add_multiple_bonds.h"
#include "atom.h"
#include "domain.h"
#include "force.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "neigh_list.h"
#include "comm.h"
#include "group.h"
#include "special.h"
#include "error.h"

using namespace LAMMPS_NS;


/* ---------------------------------------------------------------------- */

AddMultipleBonds::AddMultipleBonds(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void AddMultipleBonds::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Add_multiple_bonds command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use add_multiple_bonds unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use add_multiple_bonds with non-molecular system");

  if (narg < 3) error->all(FLERR,"Illegal add_multiple_bonds command");
  
  
  
  const int nlocal = atom->nlocal;
  tagint batom1, batom2;
  int idx1, idx2;
  int count, allcount, pos;
  
  // get the globals
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;
  
  // parse arg[0]
  int btype = force->inumeric(FLERR, arg[0]); // bond type
  if (btype <= 0 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in add_multiple_bonds command");
    
  char *ptr;
  for (int iarg = 1; iarg < narg; iarg++) {
    //parse args
    if (strchr(arg[iarg],':')) {
      ptr    = strtok(arg[iarg],":");
      batom1 = force->tnumeric(FLERR,ptr);// i atom id
      ptr    = strtok(NULL,":");
      batom2 = force->tnumeric(FLERR,ptr);// j atom id
    } else {
      error->one(FLERR,"Illegal add_multiple_bonds command");
    }
    
    //check not same atom
    if (batom1 == batom2)
      error->one(FLERR,"Illegal add_multiple_bonds command");
    

    // get the local indexes. WORKS ONLY SERIALLY
    // why serially?
    idx1 = atom->map(batom1);
    idx2 = atom->map(batom2);
    
    count = 0;
    if ((idx1 >= 0) && (idx1 < nlocal)) count++;
    if ((idx2 >= 0) && (idx2 < nlocal)) count++;
    
    MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
    if (allcount != 2)
      error->one(FLERR,"add_multiple_bonds atoms do not exist");
    
    
    

    // find the first "free" occurence
    pos=0;
    while (pos < num_bond[idx1]){
      if (bond_type[idx1][pos] < 0){
        break;
      }
      ++pos;
    }

    // check we are not overflowing
    if (pos >= atom->bond_per_atom)
    {
      error->one(FLERR,"New bond exceeded bonds per atom in add_multiple_bonds");
    }

    // set the bond
    bond_type[idx1][pos] = btype;
    bond_atom[idx1][pos] = batom2;
    if(pos == num_bond[idx1]) num_bond[idx1]++;
  }
  
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"%d bond added", narg-1);
    }
  }

  // re-trigger neighbor list build
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->command_style = "add_multiple_bonds";
}
