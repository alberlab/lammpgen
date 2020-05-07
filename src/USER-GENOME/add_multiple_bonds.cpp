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
#include "input.h"
#include "variable.h"
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
  
  
  nlocal = atom->nlocal;

  // get the globals
  num_bond = atom->num_bond;
  bond_type = atom->bond_type;
  bond_atom = atom->bond_atom;
  
  // parse arg[0]
  btype = force->inumeric(FLERR, arg[0]); // bond type
  if (btype <= 0 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in add_multiple_bonds command");
  
  //different input
  int bond_added = 0;
  if (strcmp(arg[1],"variable") == 0){
    char *ptr;
    char *batoms_str;
    
    batoms_str = input->variable->retrieve(arg[2]);

    //error->message(FLERR,batoms_str);

    ptr = strtok(batoms_str, ":");
    while (ptr != NULL)
    {
      batom1 = force->tnumeric(FLERR,ptr);// i atom id
      
      ptr = strtok (NULL, "; ");
      if (ptr == NULL) 
        error->all(FLERR,"Illegal add_multiple_bonds command");
      batom2 = force->tnumeric(FLERR,ptr);// j atom id
      
      add_bond();
      bond_added ++;
      
      ptr = strtok (NULL, ":");
    }
  } else if (strcmp(arg[1], "list") == 0){
    char *ptr;
    //error->message(FLERR,"Entering");
    for (int iarg = 2; iarg < narg; iarg++) {
      //parse args
      
      if (strchr(arg[iarg],':')) {
        ptr    = strtok(arg[iarg],":");
        batom1 = force->tnumeric(FLERR,ptr);// i atom id
        ptr    = strtok(NULL,":");
        batom2 = force->tnumeric(FLERR,ptr);// j atom id
      } else {
        error->all(FLERR,"Illegal add_multiple_bonds command");
      }
      
      add_bond();
      bond_added ++;
    }
  }
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"%d bond added", bond_added);
    }
  }

  // re-trigger neighbor list build
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->command_style = "add_multiple_bonds";
}
/* ---------------------------------------------------------------------- */

void AddMultipleBonds::add_bond()
{
  //check not same atom
  if (batom1 == batom2)
    error->one(FLERR,"Illegal add_multiple_bonds command");

  // get the local indexes. WORKS ONLY SERIALLY
  int idx1 = atom->map(batom1);
  int idx2 = atom->map(batom2);

  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;
  
  int allcount;
  MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
  if (allcount != 2)
    error->one(FLERR,"add_multiple_bonds atoms do not exist");


  // find the first "free" occurence
  int pos=0;
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