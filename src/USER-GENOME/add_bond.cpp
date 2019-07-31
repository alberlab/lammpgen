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
#include "add_bond.h"
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

#define DELTA 10000

/* ---------------------------------------------------------------------- */

AddBond::AddBond(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void AddBond::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"Add_bonds command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use add_bond unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use add_bond with non-molecular system");

  if (narg != 3) error->all(FLERR,"Illegal add_bond command");

  // parse args
  int btype = force->inumeric(FLERR, arg[0]); // bond type
  int batom1 = force->inumeric(FLERR, arg[1]); // i atom index
  int batom2 = force->inumeric(FLERR, arg[2]); // j atom index
  
  if (batom1 == batom2)
    error->all(FLERR,"Illegal add_bond command");
  
  if (btype <= 0 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in add_bond command");
  
  // get the local indexes. WORKS ONLY SERIALLY
  const int nlocal = atom->nlocal;
  const int idx1 = atom->map(batom1);
  const int idx2 = atom->map(batom2);
  
  int count = 0;
  if ((idx1 >= 0) && (idx1 < nlocal)) count++;
  if ((idx2 >= 0) && (idx2 < nlocal)) count++;
  
  int allcount;
  MPI_Allreduce(&count,&allcount,1,MPI_INT,MPI_SUM,world);
  if (allcount != 2)
    error->all(FLERR,"Add_bond atoms do not exist");
  
  // get the globals
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  // find the first "free" occurence
  int pos;
  while (pos < num_bond[idx1]){
    if (bond_type[idx1][pos] < 0){
      break;
    }
    ++pos;
  }

  // check we are not overflowing
  if (pos >= atom->bond_per_atom)
  {
    error->one(FLERR,
               "New bond exceeded bonds per atom in add_bond");
  }

  // set the bond
  bond_type[idx1][pos] = btype;
  bond_atom[idx1][pos] = batom2;
  if(pos == num_bond[idx1]) num_bond[idx1]++;

  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"One bond added");
    }
  }

  // re-trigger neighbor list build
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->command_style = "add_bond";
}
