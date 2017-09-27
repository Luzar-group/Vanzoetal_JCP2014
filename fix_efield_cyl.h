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

#ifdef FIX_CLASS

FixStyle(efield_cyl,FixEfieldCyl)

#else

#ifndef LMP_FIX_EFIELD_CYL_H
#define LMP_FIX_EFIELD_CYL_H

#include <stdlib.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixEfieldCyl : public Fix {
 public:
  FixEfieldCyl(class LAMMPS *, int, char **);
  ~FixEfieldCyl();
  int setmask();
  void init();
  void setup(int);
  void min_setup(int);
  void post_force(int);
  void min_post_force(int);
  double memory_usage();

 private:
  double ex,ey,ez;
  double center_x,center_y,center_z;
  double radius,radius_sq,height,half_height;
  double fade;
  double lim_zhi, lim_zlo;
  double fade_hi, fade_lo;
  int varflag;
  char *xstr,*ystr,*zstr;
  int xvar,yvar,zvar,xstyle,ystyle,zstyle;
  int nlevels_respa;
  double qe2f;

  int me;

  int maxatom;
  double **efield;

  int check;
  FILE *fchk;

  int nproc, myrank;
  
  double **forcefield;

  int force_flag;
  double fsum[4],fsum_all[4];
};

}

#endif
#endif
