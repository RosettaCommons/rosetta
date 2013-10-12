// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/md/Rattle.hh
/// @brief  Rattle bond length constrainer for CartesianMD
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_md_Rattle_hh
#define INCLUDED_protocols_md_Rattle_hh

#include <protocols/md/MDConstants.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/types.hh>

#include <core/optimization/types.hh>
#include <core/optimization/CartesianMinimizerMap.hh>

#include <utility/vector1.hh>
#include <cmath>

namespace protocols{
namespace md{

using namespace core;
using namespace core::optimization;
  
struct Constraint
{
  Real k;
  Size a;
  Size b;
};

class Rattle
{
public:

Rattle(){}

Rattle( pose::Pose const &pose,
	CartesianMinimizerMap const & min_map )
{
  Size const &natm( min_map.natoms() );

	sor_ = 1.25;
	maxiter_ = 100;
  moved_ = utility::vector1< bool >( natm, false );
  update_ = utility::vector1< bool >( natm, false );
  use_ = utility::vector1< bool >( natm, false );

  setup_constraint( pose, min_map );
}

~Rattle(){}

Size
ncst(){ return cst_.size(); }

inline
void 
run_rattle1( Real const &dt, 
	     Multivec &xyz, 
	     Multivec &vel,
	     Multivec const &mass ) {

  //std::cout << "run_rattle1" << std::endl;
  Real eps ( 1.0e-6 );
  Multivec const xyz_old( xyz );
  Multivec const vel0( vel );
  
  //std::cout << "1" << std::endl;
  // initialize the lists of atoms previously corrected

  for( Size i = 1; i <= use_.size(); ++i ){
    if( use_[i] ) {
      moved_[i] = true;
    } else {
      moved_[i] = false;
    }
    update_[i] = false;
  }
  
  // apply rattle to distances and half-step velocity values
  Size niter = 0;
  bool done = false;
  utility::vector1< Real > dxyz1( 3, 0.0 );
  utility::vector1< Real > dxyz2( 3, 0.0 );
  
  while ( ! done && niter <= maxiter_ ){
    niter ++;
    done = true;
    //std::cout << "2, iter: " << niter << std::endl;
    
    for( Size i = 1; i<=cst_.size(); ++i ){
      Constraint const &cst = cst_[i];
      Size const &ia (cst.a);
      Size const &ib (cst.b);
      
      if ( !(moved_[ia] || moved_[ib] ) ) continue;
	
      dxyz1[1] = xyz[3*ib-2]-xyz[3*ia-2];
      dxyz1[2] = xyz[3*ib-1]-xyz[3*ia-1];
      dxyz1[3] = xyz[3*ib  ]-xyz[3*ia  ];
      
      Real dist2 = dxyz1[1]*dxyz1[1] + dxyz1[2]*dxyz1[2] + dxyz1[3]*dxyz1[3];
      
      Real delta = (cst.k)*(cst.k) - dist2;

      /*
      std::cout << "rattle1, cst i/a/b: " << std::setw(4) << i;
      std::cout << " " << std::setw(4) << ia << " " << std::setw(4) << ib;
      std::cout << " " << std::setw(8) << std::sqrt(dist2);
      std::cout << " " << std::setw(8) << cst.k;
      std::cout << " " << std::setw(8) << delta << std::endl;
      */
      
      if ( std::abs(delta) >= eps ){
	done = false;
	update_[ia] = true;
	update_[ib] = true;
	
	dxyz2[1] = xyz_old[3*ib-2]-xyz[3*ia-2];
	dxyz2[2] = xyz_old[3*ib-1]-xyz[3*ia-1];
	dxyz2[3] = xyz_old[3*ib  ]-xyz[3*ia  ];
	
	Real const dot = dxyz1[1]*dxyz2[1] + dxyz1[2]*dxyz2[2] + dxyz1[3]*dxyz2[3];
	
	Real const rma = 1.0 / mass[ia];
	Real const rmb = 1.0 / mass[ib];
	Real const term = sor_ * delta / (2.0 * (rma+rmb) * dot);
	
	for ( Size k = 1; k <= 3; ++k ){
	  Size const dof_a = 3*ia-3+k;
	  Size const dof_b = 3*ib-3+k;
	  
	  xyz[dof_a] -= dxyz2[k]*term*rma;
	  xyz[dof_b] += dxyz2[k]*term*rmb;
	  
	  vel[dof_a] -= term*rma/dt;
	  vel[dof_b] += term*rmb/dt;
	}

      } // iter over Rattle ID
    } // iter
    
    for( Size i = 1; i <= use_.size(); ++i ){
      moved_[i] = update_[i];
      update_[i] = false;
    }
  }

  //std::cout << " Rattle1 finished after " << niter << std::endl;
  /*
  // apply group position and velocity constraints via exact reset
  for( Size i = 1; i <= nratx_; ++i ) {
    Size ia ( iratx_[i] );

    for( Size j = grp_[ia].start; j <= grp_[ia].end; ++j ){
	
      Size atmno ( grp_[ia].atmno[j] );
      Real weigh ( mass[atmno] / grp_[ia].mass );
      
      for( Size k = 1; k <= 3; ++k ){
	Size const i_dof ( 3*atmno-3+k );
	
	xyz[i_dof] -= xyz[i_dof]*weigh;
	vel[i_dof] -= vel[i_dof]*weigh;
      }
    }
  }
  */

  for( Size i = 1; i <= cst_.size(); i++ ){
    Constraint &cst = cst_[i];
    Size const ia (cst.a);
    Size const ib (cst.b);

    /*
    std::cout << "Rattle1a " << i << " " << ia;
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ia-2];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ia-1];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ia];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ia-2];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ia-1];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ia];
    std::cout << std::endl;

    std::cout << "Rattle1b " << i << " " << ib;
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ib-2];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ib-1];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz[3*ib];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ib-2];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ib-1];
    std::cout << " " << std::setprecision(4) << std::setw(8) << xyz_old[3*ib];
    std::cout << std::endl;

    std::cout << " " << std::setw(10) << vel[3*ia-2];
    std::cout << " " << std::setw(10) << vel[3*ia-1];
    std::cout << " " << std::setw(10) << vel[3*ia];
    std::cout << " " << std::setw(10) << vel0[3*ib-2];
    std::cout << " " << std::setw(10) << vel0[3*ib-1];
    std::cout << " " << std::setw(10) << vel0[3*ib];
    */
  }

  return;
} // run_rattle1

inline
void 
run_rattle2( Real const &dt, 
	     Multivec &xyz,
	     Multivec &vel,
	     Multivec const &mass ) {

  //std::cout << "run_rattle2" << std::endl;
  Real eps ( 1.0e-6 / dt );
  Real vterm ( 2.0 / (dt * MDForceFactor));
  Multivec vel0( vel );
  
  for( Size i = 1; i <= use_.size(); i++ ){
    if ( use_[i] ){
      moved_[i] = true;
    } else {
      moved_[i] = false;
    }
    update_[i] = false;
  }

  //     set the iteration counter, termination and tolerance
  Size niter ( 0 );
  bool done ( false );
  
  Real vxx, vyx, vzx, vyy, vzy, vzz;
  utility::vector1< Real > dxyz( 3, 0.0 );
  utility::vector1< Real > dvel( 3, 0.0 );

  // apply the rattle algorithm to correct the velocities
  while ( !done && niter<=maxiter_ ){
    niter++;
    done = true;

    for( Size i = 1; i <= cst_.size(); i++ ){
      Constraint &cst = cst_[i];
      Size const ia (cst.a);
      Size const ib (cst.b);

      if ( !(moved_[ia] || moved_[ib] ) ) continue;

      dxyz[1] = xyz[3*ib-2]-xyz[3*ia-2];
      dxyz[2] = xyz[3*ib-1]-xyz[3*ia-1];
      dxyz[3] = xyz[3*ib  ]-xyz[3*ia  ];
      dvel[1] = vel[3*ib-2]-vel[3*ia-2];
      dvel[2] = vel[3*ib-1]-vel[3*ia-1];
      dvel[3] = vel[3*ib  ]-vel[3*ia  ];

      Real dot = dxyz[1]*dvel[1] + dxyz[2]*dvel[2] + dxyz[3]*dvel[3];

      Real const rma ( 1.0 / mass[ia] );
      Real const rmb ( 1.0 / mass[ib] );
      Real term ( -dot * sor_ / ((rma+rmb) *(cst.k)*(cst.k)) );

      if (std::abs(term) >= eps){
	done = false;
	update_[ia] = true;
	update_[ib] = true;
	//term = sor_ * term;
	
	for ( Size k = 1; k <= 3; ++k ){
	  Size const dof_a ( 3*ia-3+k );
	  Size const dof_b ( 3*ib-3+k );

	  vel[dof_a] -= dxyz[k] * term *rma;
	  vel[dof_b] += dxyz[k] * term *rmb;
	}
      }
    }

    for( Size i = 1; i <= use_.size(); ++i ){
      moved_[i] = update_[i];
      update_[i] = false;
    }
  }


  /*
  std::cout << " Rattle2 finished after " << niter << std::endl;

  for( Size i = 1; i <= cst_.size(); i++ ){
    Constraint &cst = cst_[i];
    Size const ia (cst.a);
    Size const ib (cst.b);
    std::cout << "Rattle2 " << i;
    std::cout << " " << vel[3*ia-2] << " " << vel[3*ia-1] << " " << vel[3*ia];
    std::cout << " " << vel[3*ib-2] << " " << vel[3*ib-1] << " " << vel[3*ib];
    std::cout << std::endl;
  }
  */
  
  /*
  // apply any atom group velocity constraints via exact reset
  for( Size i = 1; i <= grp_.size(); ++i ){
    Size ia( grp_[i].a );
    utility::vector1< Real > xv( 3, 0.0 );

    for( Size j = grp[grpno].start; j <=grp[grpno].stop; ++j ){
      Size atmno = grp[j].k;
      Real weigh = mass[atmno] / grp_[grpno].mass;

      for( Size k = 1; k <= 3; ++k ){
	xv[k] += vel[i_dof]*weigh;
	vel[i_dof] -= xv[k];
      }
    }
  }
  */
  
  return;
} // run_rattle2

private:

inline
void
setup_constraint( pose::Pose const & pose, 
		  CartesianMinimizerMap const & min_map ){

  cst_.resize( 0 );

  for( Size iatm = 1; iatm <= (Size)(min_map.natoms()); ++iatm ){
    use_[iatm] = true;

    id::AtomID AtomID = min_map.get_atom( iatm );
    Size resno = AtomID.rsd();
    Size atmno = AtomID.atomno();
    
    conformation::Residue const &rsd( pose.residue( resno ) );

    if( rsd.atom_is_hydrogen( atmno ) ){

      Size const &atmno_base( rsd.atom_base( atmno ) );

      id::AtomID BaseID( atmno_base, resno );

      Constraint cst;

      cst.a = min_map.get_atom_index( BaseID );
      cst.b = iatm;
      cst.k = rsd.icoor( atmno ).d();
      cst_.push_back(cst);
    }
  }
  return;
}

private:
  // Parameters: Max iteration and tolerance
  Size maxiter_;
  Real sor_;

  // Rattle Constraint Structures
  utility::vector1< Constraint > cst_;

  utility::vector1< bool > use_;
  utility::vector1< bool > moved_;
  utility::vector1< bool > update_;

}; // class Rattle

} // namespace md
} // namespace 

#endif
