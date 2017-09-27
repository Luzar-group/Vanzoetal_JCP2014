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

/* ----------------------------------------------------------------------
   Contributing author: Davide Vanzo
------------------------------------------------------------------------- */

#include "math.h"
#include "string.h"
#include "stdlib.h"
#include "fix_efield_cyl.h"
#include "atom.h"
#include "update.h"
#include "modify.h"
#include "force.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "math_const.h"
#include "memory.h"
#include "error.h"
#include "domain.h"

using namespace LAMMPS_NS;
using namespace FixConst;
using namespace MathConst;

enum{CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

FixEfieldCyl::FixEfieldCyl(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg != 13) error->all(FLERR,"Illegal fix efield_cyl command");

  qe2f = force->qe2f;
  xstr = ystr = zstr = NULL;

  if (strstr(arg[3],"v_") == arg[3]) {
    int n = strlen(&arg[3][2]) + 1;
    xstr = new char[n];
    strcpy(xstr,&arg[3][2]);
  } else {
    ex = qe2f * atof(arg[3]);
    xstyle = CONSTANT;
  }

  if (strstr(arg[4],"v_") == arg[4]) {
    int n = strlen(&arg[4][2]) + 1;
    ystr = new char[n];
    strcpy(ystr,&arg[4][2]);
  } else {
    ey = qe2f * atof(arg[4]);
    ystyle = CONSTANT;
  }

  if (strstr(arg[5],"v_") == arg[5]) {
    int n = strlen(&arg[5][2]) + 1;
    zstr = new char[n];
    strcpy(zstr,&arg[5][2]);
  } else {
    ez = qe2f * atof(arg[5]);
    zstyle = CONSTANT;
  }

  center_x = atof(arg[6]);
  center_y = atof(arg[7]);
  center_z = atof(arg[8]);

  radius = atof(arg[9]);
  height = atof(arg[10]);
  half_height = height / 2.0;

  fade = atof(arg[11]);

  check = atoi(arg[12]);

  // Calculate and check field limits

  if (fade <= 0.000 || radius <= 0.000) {
    error->all(FLERR,"Cylinder radius and fading thickness must be positive.");
  }
  
  radius_sq = radius*radius;

  lim_zhi = center_z + ( height / 2.0 );
  lim_zlo = center_z - ( height / 2.0 );

  if (lim_zhi > domain->boxhi[2] || lim_zlo < domain->boxlo[2]) {
    error->all(FLERR,"The cylindrical volume is out of the z box dimension.");
  }

  fade_hi = radius + fade;
  fade_lo = radius;

  if ((center_x + fade_hi) > domain->boxhi[0] || (center_x - fade_hi) < domain->boxlo[0]) {
    error->all(FLERR,"The cylindrical volume is out of the x box dimension.");
  }

  if ((center_y + fade_hi) > domain->boxhi[1] || (center_y - fade_hi) < domain->boxlo[1]) {
    error->all(FLERR,"The cylindrical volume is out of the y box dimension.");
  }

  maxatom = 0;
  efield = NULL;

  // Flags for per-atom data dump
  if (check == 1) {
    peratom_flag = 1;
    size_peratom_cols = 6;
    peratom_freq = 1;
    
    int nmax = atom->nmax;
    
    forcefield = NULL;
    memory->grow(forcefield,(nmax*10),6,"efieldcyl:forcefield");
    array_atom = forcefield;
    
    atom->add_callback(0);
    atom->add_callback(1);
  }

}

/* ---------------------------------------------------------------------- */

FixEfieldCyl::~FixEfieldCyl()
{
  if (check == 1) {
    atom->delete_callback(id,0);
    atom->delete_callback(id,1);
    memory->destroy(forcefield);
  }
  
  delete [] xstr;
  delete [] ystr;
  delete [] zstr;
  memory->destroy(efield);
}

/* ---------------------------------------------------------------------- */

int FixEfieldCyl::setmask()
{
  int mask = 0;
  mask |= POST_FORCE;
  mask |= MIN_POST_FORCE;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixEfieldCyl::init()
{
  MPI_Comm_rank(world,&me);
  
  if (me == 0) {
    if (screen) fprintf(screen,"Powering up local cylindrical electric field Ez=%8.5f with cosine smoothing...\n", ez/qe2f);
    if (logfile) fprintf(logfile,"Powering up local cylindrical electric field Ez=%8.5f with cosine smoothing...\n", ez/qe2f);
  }

  if (!atom->q_flag) error->all(FLERR,"Fix efield requires atom attribute q");

  // check variables

  if (xstr) {
    xvar = input->variable->find(xstr);
    if (xvar < 0) error->all(FLERR,"Variable name for fix efield does not exist");
    if (input->variable->equalstyle(xvar)) xstyle = EQUAL;
    else if (input->variable->atomstyle(xvar)) xstyle = ATOM;
    else error->all(FLERR,"Variable for fix efield is invalid style");
  }
  if (ystr) {
    yvar = input->variable->find(ystr);
    if (yvar < 0) error->all(FLERR,"Variable name for fix efield does not exist");
    if (input->variable->equalstyle(yvar)) ystyle = EQUAL;
    else if (input->variable->atomstyle(yvar)) ystyle = ATOM;
    else error->all(FLERR,"Variable for fix efield is invalid style");
  }
  if (zstr) {
    zvar = input->variable->find(zstr);
    if (zvar < 0) error->all(FLERR,"Variable name for fix efield does not exist");
    if (input->variable->equalstyle(zvar)) zstyle = EQUAL;
    else if (input->variable->atomstyle(zvar)) zstyle = ATOM;
    else error->all(FLERR,"Variable for fix efield is invalid style");
  }

  if (xstyle == ATOM || ystyle == ATOM || zstyle == ATOM)
    varflag = ATOM;
  else if (xstyle == EQUAL || ystyle == EQUAL || zstyle == EQUAL)
    varflag = EQUAL;
  else varflag = CONSTANT;

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"RESPA not supported with fix efield_cyl");
}

/* ---------------------------------------------------------------------- */

void FixEfieldCyl::setup(int vflag)
{
  if (strstr(update->integrate_style,"verlet"))
    post_force(vflag);
}

/* ---------------------------------------------------------------------- */

void FixEfieldCyl::min_setup(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   apply F = qE
------------------------------------------------------------------------- */

void FixEfieldCyl::post_force(int vflag)
{
  double **x = atom->x;
  double **f = atom->f;
  double *q = atom->q;
  double versor[3];
  double rmod, r[3];
  double scale_x,scale_y,scale_z;
  double dr, dz;
  double ef_force_x, ef_force_y, ef_force_z;
  int *mask = atom->mask;
  int nlocal = atom->nlocal;

  fsum[0] = fsum[1] = fsum[2] = fsum[3] = 0.0;
  force_flag = 0;

  if (varflag == CONSTANT) {
    
    for (int i = 0; i < nlocal; i++)
      if (mask[i] & groupbit) {

	ef_force_x = 0.0;
	ef_force_y = 0.0;
	ef_force_z = 0.0;

	// Check z direction

	if ( x[i][2] <= lim_zhi && x[i][2] >= lim_zlo ) {

	  // Calculate 2D versor wrt the center of the cylinder

	  r[0] = x[i][0] - center_x;
	  r[1] = x[i][1] - center_y;
	  r[2] = x[i][2] - center_z;

	  rmod = sqrt(r[0]*r[0] + r[1]*r[1]);
	  
	  versor[0] = r[0] / rmod;
	  versor[1] = r[1] / rmod;
	  
	  // Fading shell

	  if ( rmod > fade_lo && rmod < fade_hi ) {
	    
	    // Smoothing function - z direction

	    dr = rmod - fade_lo;
	    dr /= fade;

	    scale_z = 0.5 * (cos(dr*MY_PI) + 1.0);

	    // Smoothing function - r direction

	    dz = x[i][2] - center_z;
	    
	    scale_x = -0.5 * MY_PI * sin(dr*MY_PI) * dz / fade;

	    scale_y = scale_x;

	    scale_x *= versor[0];
	    scale_y *= versor[1];

	    // Compute forces
	    
	    ef_force_x = q[i]*(ex*scale_x);
	    ef_force_y = q[i]*(ey*scale_y);
	    ef_force_z = q[i]*(ez*scale_z);
	    
	  }

	  // Inside field

	  else if (rmod <= fade_lo) {

	    ef_force_x = 0.0;
	    ef_force_y = 0.0;
	    ef_force_z = q[i]*ez;

	  }

	}

	// Dump force field if requested

	if (check == 1) {
	  
	  forcefield[i][0] = x[i][0];
	  forcefield[i][1] = x[i][1];
	  forcefield[i][2] = x[i][2];

	  forcefield[i][3] = ef_force_x;
	  forcefield[i][4] = ef_force_y;
	  forcefield[i][5] = ef_force_z;
  
	}
	
	// Update force
	
	f[i][0] += ef_force_x;
	f[i][1] += ef_force_y;
	f[i][2] += ef_force_z;

	fsum[0] -= ef_force_x*x[i][0]+ef_force_y*x[i][1]+ef_force_z*x[i][2];
        fsum[1] += ef_force_x;
        fsum[2] += ef_force_y;
        fsum[3] += ef_force_z;

      }
  } else 
    {
      error->all(FLERR,"Variable efield not supported by fix_efield_cyl.");
    }
}

/* ---------------------------------------------------------------------- */

void FixEfieldCyl::min_post_force(int vflag)
{
  post_force(vflag);
}

/* ----------------------------------------------------------------------
   memory usage of local atom-based array
------------------------------------------------------------------------- */

double FixEfieldCyl::memory_usage()
{
  double bytes = 0.0;
  if (varflag == ATOM) bytes = atom->nmax*3 * sizeof(double);
  return bytes;
}


