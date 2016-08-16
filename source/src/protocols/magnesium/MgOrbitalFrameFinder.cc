// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/MgOrbitalFrameFinder.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgOrbitalFrameFinder.hh>
#include <protocols/magnesium/util.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <numeric/UniformRotationSampler.hh>
#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "protocols.magnesium.MgOrbitalFrameFinder" );

using namespace core;
using utility::vector1;
using namespace ObjexxFCL::format;

namespace protocols {
namespace magnesium {

//Constructor
MgOrbitalFrameFinder::MgOrbitalFrameFinder()
//  legacy_mode_( true )
{}

//Destructor
MgOrbitalFrameFinder::~MgOrbitalFrameFinder()
{}

///////////////////////////////////////////////////////////////////////////////
void
MgOrbitalFrameFinder::apply( pose::Pose & pose )
{
	vector1< Size > mg_res = get_mg_res( pose );
	for ( Size n = 1; n <= mg_res.size(); n++ ) {
		determine_mg_orbital_frame( pose, mg_res[ n ] );
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// sets V1, V2, ... V6 coordinates of Mg(2+) to point towards ligands --
//  a 'frame' for the octahedral ligand field.
//
///////////////////////////////////////////////////////////////////////////////
void
MgOrbitalFrameFinder::determine_mg_orbital_frame( pose::Pose & pose,
	Size const i /* mg2+ res number*/ ) {

	clock_t start_time = clock();

	vector1< core::id::AtomID > ligands = get_mg_ligands( pose, i );
	point_orbitals_to_closest_ligands( pose, i, ligands );

	sample_orbital_frame( pose, i, ligands );

	clock_t end_time = clock();
	TR << "orbital frame finder finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

}


///////////////////////////////////////////////////////////////////////////////
void
MgOrbitalFrameFinder::point_orbitals_to_closest_ligands( pose::Pose & pose,
	Size const i /* mg2+ res number*/,
	utility::vector1< core::id::AtomID > const & ligands ) {

	using namespace core::id;

	Vector const xyz_mg =  pose.xyz( NamedAtomID( "MG  ", i ) );

	// try to set a frame.
	// first look at closest ligand (and also any ligands that are exactly opposed).
	Vector x( 1.0, 0.0, 0.0 ), y( 0.0, 1.0, 0.0 ), z( 0.0, 0.0, 1.0 );
	if ( ligands.size() > 0 ) {
		x = ( pose.residue( ligands[ 1 ].rsd() ).xyz( ligands[ 1 ].atomno() ) - xyz_mg ).normalized();
		// find any other ligands that are on the 'other side'
		for ( Size q = 2; q <= ligands.size(); q++ ) {
			Vector x_opposite = ( pose.residue( ligands[ q ].rsd() ).xyz( ligands[ q ].atomno() ) - xyz_mg ).normalized();
			Real const dotprod = dot(x_opposite, x );
			if ( dotprod < -0.9 ) {
				Vector xnew = ( 0.5 * ( x - x_opposite ) ).normalized();
				x = xnew;
				break;
			}
		}
	}

	if ( ligands.size() > 1 ) {
		y = ( pose.residue( ligands[ 2 ].rsd() ).xyz( ligands[ 2 ].atomno() ) - xyz_mg ).normalized();
		if ( dot( x, y ) < -0.9 && ligands.size() > 2 ) { // in case ligand #2 is opposite ligand #1
			y = ( pose.residue( ligands[ 3 ].rsd() ).xyz( ligands[ 3 ].atomno() ) - xyz_mg ).normalized();
		}
		// find any other ligands that are on the 'other side'
		for ( Size q = 2; q <= ligands.size(); q++ ) {
			Vector y_opposite = ( pose.residue( ligands[ q ].rsd() ).xyz( ligands[ q ].atomno() ) - xyz_mg ).normalized();
			Real const dotprod = dot(y_opposite, y );
			if ( dotprod < -0.9 ) {
				y = ( 0.5 * ( y - y_opposite ) ).normalized();
				break;
			}
		}
	}

	x.normalize();
	z = cross( x, y );
	z.normalize();
	y = cross( z, x );
	y.normalize();

	// slow -- might be better to separate into a separate function, and
	//  have this function keep the pose const.
	pose.set_xyz( NamedAtomID( "V1", i ), xyz_mg + x );
	pose.set_xyz( NamedAtomID( "V2", i ), xyz_mg + y );
	pose.set_xyz( NamedAtomID( "V3", i ), xyz_mg + z );
	pose.set_xyz( NamedAtomID( "V4", i ), xyz_mg - x );
	pose.set_xyz( NamedAtomID( "V5", i ), xyz_mg - y );
	pose.set_xyz( NamedAtomID( "V6", i ), xyz_mg - z );
}

///////////////////////////////////////////////////////////////////////////////
// sum of dot-products of ligand & orbital atom unit vectors (with Mg2+ at origin)
Real
MgOrbitalFrameFinder::get_orbital_frame_score_upon_rotation( utility::vector1< Vector > const & r_lig,
	utility::vector1< Vector > const & v_lig,
	numeric::xyzMatrix< Real > const & R ) {

	Real score( 0.0 );

	for ( Size m = 1; m <= r_lig.size(); m++ ) {
		Vector const & r = r_lig[ m ];

		// find 'orbital'  that aligns best with each ligand
		Real best_dot( -1.0 );
		for ( Size n = 1; n <= v_lig.size(); n++ ) {
			Real current_dot = dot( r, ( R * v_lig[ n ] ) );
			if ( current_dot > best_dot ) best_dot = current_dot;
		}
		//runtime_assert( best_dot >= 0.0 );
		score += best_dot;

	}

	return score;
}


///////////////////////////////////////////////////////////////////////////////
void
MgOrbitalFrameFinder::sample_orbital_frame( pose::Pose & pose,
	Size const i /* mg2+ res number*/,
	utility::vector1< core::id::AtomID > const & ligands ) {

	using namespace core::id;

	Vector const xyz_mg =  pose.xyz( NamedAtomID( "MG  ", i ) );
	if ( urs_ == 0 ) urs_ = get_octahedral_uniform_rotation_sampler();

	utility::vector1< Vector > r_lig;
	for ( Size n = 1; n <= ligands.size(); n++ ) {
		r_lig.push_back(  ( pose.xyz( ligands[ n ] ) - xyz_mg ).normalized() );
	}

	utility::vector1< Vector > v_lig;
	for ( Size n = 1; n <= 6; n++ ) {
		v_lig.push_back(  ( pose.xyz( NamedAtomID( "V"+I(1,n), i ) ) - xyz_mg ).normalized() );
	}

	numeric::xyzMatrix< core::Real > R, best_R( numeric::xyzMatrix<Real>::identity() );
	Real best_score = get_orbital_frame_score_upon_rotation( r_lig, v_lig, best_R );
	Real const orig_score = best_score;

	for ( core::Size ctr=1; ctr<=urs_->nrots() ; ++ctr ) {
		urs_->get(ctr, R);

		Real score = get_orbital_frame_score_upon_rotation( r_lig, v_lig, R );
		if ( score > best_score ) {
			best_score = score;
			best_R = R;
		}
	}

	TR << " For Mg(2+) with " << r_lig.size() << " ligands, original score: " << orig_score << "  optimized_score: " << best_score << std::endl;

	for ( Size n = 1; n <= 6; n++ ) {
		pose.set_xyz( NamedAtomID( "V"+I(1,n), i ), ( best_R * v_lig[ n ] ) + xyz_mg );
	}

	TR << "sampled: " << urs_->nrots() << " frames." << std::endl;

}


} //magnesium
} //protocols
