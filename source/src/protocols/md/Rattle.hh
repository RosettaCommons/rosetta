// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

namespace protocols {
namespace md {

struct Constraint
{
	core::Real k;
	core::Size a;
	core::Size b;
};

class Rattle
{
public:

	Rattle(){}

	Rattle( core::pose::Pose const &pose,
		core::optimization::CartesianMinimizerMap const & min_map )
	{
		core::Size const &natm( min_map.natoms() );

		sor_ = 1.25;
		maxiter_ = 100;
		moved_ = utility::vector1< bool >( natm, false );
		update_ = utility::vector1< bool >( natm, false );
		use_ = utility::vector1< bool >( natm, false );

		setup_constraint( pose, min_map );
	}

	~Rattle()= default;

	core::Size
	ncst(){ return cst_.size(); }

	inline
	void
	run_rattle1( core::Real const &dt,
		core::optimization::Multivec &xyz,
		core::optimization::Multivec &vel,
		core::optimization::Multivec const &mass )
	{
		//std::cout << "run_rattle1" << std::endl;
		core::Real eps ( 1.0e-6 );
		core::optimization::Multivec const xyz_old( xyz );
		core::optimization::Multivec const vel0( vel );

		//std::cout << "1" << std::endl;
		// initialize the lists of atoms previously corrected

		for ( core::Size i = 1; i <= use_.size(); ++i ) {
			if ( use_[i] ) {
				moved_[i] = true;
			} else {
				moved_[i] = false;
			}
			update_[i] = false;
		}

		// apply rattle to distances and half-step velocity values
		core::Size niter = 0;
		bool done = false;
		utility::vector1< core::Real > dxyz1( 3, 0.0 );
		utility::vector1< core::Real > dxyz2( 3, 0.0 );

		while ( ! done && niter <= maxiter_ ) {
			niter ++;
			done = true;
			//std::cout << "2, iter: " << niter << std::endl;

			for ( core::Size i = 1; i<=cst_.size(); ++i ) {
				Constraint const &cst = cst_[i];
				core::Size const &ia (cst.a);
				core::Size const &ib (cst.b);

				if ( !(moved_[ia] || moved_[ib] ) ) continue;

				dxyz1[1] = xyz[3*ib-2]-xyz[3*ia-2];
				dxyz1[2] = xyz[3*ib-1]-xyz[3*ia-1];
				dxyz1[3] = xyz[3*ib  ]-xyz[3*ia  ];

				core::Real dist2 = dxyz1[1]*dxyz1[1] + dxyz1[2]*dxyz1[2] + dxyz1[3]*dxyz1[3];

				core::Real delta = (cst.k)*(cst.k) - dist2;

				/*
				std::cout << "rattle1, cst i/a/b: " << std::setw(4) << i;
				std::cout << " " << std::setw(4) << ia << " " << std::setw(4) << ib;
				std::cout << " " << std::setw(8) << std::sqrt(dist2);
				std::cout << " " << std::setw(8) << cst.k;
				std::cout << " " << std::setw(8) << delta << std::endl;
				*/

				if ( std::abs(delta) >= eps ) {
					done = false;
					update_[ia] = true;
					update_[ib] = true;

					dxyz2[1] = xyz_old[3*ib-2]-xyz[3*ia-2];
					dxyz2[2] = xyz_old[3*ib-1]-xyz[3*ia-1];
					dxyz2[3] = xyz_old[3*ib  ]-xyz[3*ia  ];

					core::Real const dot = dxyz1[1]*dxyz2[1] + dxyz1[2]*dxyz2[2] + dxyz1[3]*dxyz2[3];

					core::Real const rma = 1.0 / mass[ia];
					core::Real const rmb = 1.0 / mass[ib];
					core::Real const term = sor_ * delta / (2.0 * (rma+rmb) * dot);

					for ( core::Size k = 1; k <= 3; ++k ) {
						core::Size const dof_a = 3*ia-3+k;
						core::Size const dof_b = 3*ib-3+k;

						xyz[dof_a] -= dxyz2[k]*term*rma;
						xyz[dof_b] += dxyz2[k]*term*rmb;

						vel[dof_a] -= term*rma/dt;
						vel[dof_b] += term*rmb/dt;
					}
				} // iter over Rattle ID
			} // iter

			for ( core::Size i = 1; i <= use_.size(); ++i ) {
				moved_[i] = update_[i];
				update_[i] = false;
			}
		}

		//std::cout << " Rattle1 finished after " << niter << std::endl;

		/* for( core::Size i = 1; i <= cst_.size(); i++ ){
		Constraint &cst = cst_[i];

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
		} */

		return;
	} // run_rattle1

	inline
	void
	run_rattle2( core::Real const &dt,
		core::optimization::Multivec &xyz,
		core::optimization::Multivec &vel,
		core::optimization::Multivec const &mass ) {

		//std::cout << "run_rattle2" << std::endl;
		core::Real eps ( 1.0e-6 / dt );
		core::optimization::Multivec vel0( vel );

		for ( core::Size i = 1; i <= use_.size(); i++ ) {
			if ( use_[i] ) {
				moved_[i] = true;
			} else {
				moved_[i] = false;
			}
			update_[i] = false;
		}

		// set the iteration counter, termination and tolerance
		core::Size niter ( 0 );
		bool done ( false );

		utility::vector1< core::Real > dxyz( 3, 0.0 );
		utility::vector1< core::Real > dvel( 3, 0.0 );

		// apply the rattle algorithm to correct the velocities
		while ( !done && niter<=maxiter_ ) {
			niter++;
			done = true;

			for ( core::Size i = 1; i <= cst_.size(); i++ ) {
				Constraint &cst = cst_[i];
				core::Size const ia (cst.a);
				core::Size const ib (cst.b);

				if ( !(moved_[ia] || moved_[ib] ) ) continue;

				dxyz[1] = xyz[3*ib-2]-xyz[3*ia-2];
				dxyz[2] = xyz[3*ib-1]-xyz[3*ia-1];
				dxyz[3] = xyz[3*ib  ]-xyz[3*ia  ];
				dvel[1] = vel[3*ib-2]-vel[3*ia-2];
				dvel[2] = vel[3*ib-1]-vel[3*ia-1];
				dvel[3] = vel[3*ib  ]-vel[3*ia  ];

				core::Real dot = dxyz[1]*dvel[1] + dxyz[2]*dvel[2] + dxyz[3]*dvel[3];

				core::Real const rma ( 1.0 / mass[ia] );
				core::Real const rmb ( 1.0 / mass[ib] );
				core::Real term ( -dot * sor_ / ((rma+rmb) *(cst.k)*(cst.k)) );

				if ( std::abs(term) >= eps ) {
					done = false;
					update_[ia] = true;
					update_[ib] = true;
					//term = sor_ * term;

					for ( core::Size k = 1; k <= 3; ++k ) {
						core::Size const dof_a ( 3*ia-3+k );
						core::Size const dof_b ( 3*ib-3+k );

						vel[dof_a] -= dxyz[k] * term *rma;
						vel[dof_b] += dxyz[k] * term *rmb;
					}
				}
			}

			for ( core::Size i = 1; i <= use_.size(); ++i ) {
				moved_[i] = update_[i];
				update_[i] = false;
			}
		}



		return;
	} // run_rattle2

private:

	inline
	void
	setup_constraint( core::pose::Pose const & pose,
		core::optimization::CartesianMinimizerMap const & min_map ){

		using namespace core;

		cst_.resize( 0 );

		for ( core::Size iatm = 1; iatm <= (core::Size)(min_map.natoms()); ++iatm ) {
			use_[iatm] = true;

			id::AtomID AtomID = min_map.get_atom( iatm );
			core::Size resno = AtomID.rsd();
			core::Size atmno = AtomID.atomno();

			conformation::Residue const &rsd( pose.residue( resno ) );

			if ( rsd.atom_is_hydrogen( atmno ) ) {

				core::Size const &atmno_base( rsd.atom_base( atmno ) );

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
	core::Size maxiter_;
	core::Real sor_;

	// Rattle Constraint Structures
	utility::vector1< Constraint > cst_;

	utility::vector1< bool > use_;
	utility::vector1< bool > moved_;
	utility::vector1< bool > update_;

}; // class Rattle

} // namespace md
} // namespace

#endif
