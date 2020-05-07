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
#include "remove_bond.h"
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

RemoveBond::RemoveBond(LAMMPS *lmp) : Pointers(lmp) {}

/* ---------------------------------------------------------------------- */

void RemoveBond::command(int narg, char **arg)
{
  if (domain->box_exist == 0)
    error->all(FLERR,"remove_bond command before simulation box is defined");
  if (atom->tag_enable == 0)
    error->all(FLERR,"Cannot use remove_bond unless atoms have IDs");
  if (atom->molecular != 1)
    error->all(FLERR,"Cannot use remove_bond with non-molecular system");

  if (narg != 3) error->all(FLERR,"Illegal remove_bond command");

  // parse args
  int btype = force->inumeric(FLERR, arg[0]); // bond type
  int itag = force->inumeric(FLERR, arg[1]); // i atom index
  int jtag = force->inumeric(FLERR, arg[2]); // j atom index

  if (btype <= 0 || btype > atom->nbondtypes)
    error->all(FLERR,"Invalid bond type in add_bond command");
  
  // get the local indexes. WORKS ONLY SERIALLY
  int ii = atom->map(itag);
  int jj = atom->map(jtag); 

  // get the globals
  int *num_bond = atom->num_bond;
  int **bond_type = atom->bond_type;
  tagint **bond_atom = atom->bond_atom;

  // find the first "free" occurence
  int pos;
  while (pos < num_bond[ii]){
    if (bond_type[ii][pos] == btype && bond_atom[ii][pos] == jtag){
      break;
    }
    ++pos;
  }

  // check we are not overflowing
  if (pos >= atom->bond_per_atom)
  {
    error->one(FLERR,
               "No matching bond found in remove_bond");
  }

  // set the bond
  bond_type[ii][pos] = -btype;
  
  if (comm->me == 0) {
    if (screen) {
      fprintf(screen,"Bond removed\n");
    }
  }

  // re-trigger neighbor list build
  int irequest = neighbor->request(this);
  neighbor->requests[irequest]->pair = 0;
  neighbor->requests[irequest]->command = 1;
  neighbor->requests[irequest]->occasional = 1;
  neighbor->requests[irequest]->command_style = "remove_bond";
}
