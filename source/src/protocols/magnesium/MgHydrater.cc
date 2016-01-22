// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/magnesium/MgHydrater.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgHydrater.hh>
#include <protocols/magnesium/MgWaterHydrogenPacker.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/magnesium/params.hh>
#include <core/chemical/AtomType.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/EulerAngles.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.MgHydrater" );

using namespace core;
using utility::vector1;
typedef  numeric::xyzMatrix< core::Real > Matrix;
using namespace ObjexxFCL::format;

namespace protocols {
namespace magnesium {

//Constructor
MgHydrater::MgHydrater():
	excise_mini_pose_( true ),
	use_fast_frame_heuristic_( true ),
	mg_water_hydrogen_packer_( new MgWaterHydrogenPacker ),
	verbose_( false )
{
	mg_water_hydrogen_packer_->set_excise_mini_pose( false );
}

//Constructor
MgHydrater::MgHydrater( utility::vector1< Size > const & mg_res_list ):
	mg_res_list_( mg_res_list ),
	excise_mini_pose_( true ),
	use_fast_frame_heuristic_( true ),
	mg_water_hydrogen_packer_( new MgWaterHydrogenPacker ),
	verbose_( false )
{
	mg_water_hydrogen_packer_->set_excise_mini_pose( false );
}

//Destructor
MgHydrater::~MgHydrater()
{}

///////////////////////////////////////////
void
MgHydrater::apply( pose::Pose & pose ) {

	// find all the magnesiums in here.
	for ( Size i = 1; i <= pose.total_residue(); i++ ) {
		core::conformation::Residue const & rsd_i = pose.residue( i );
		if ( rsd_i.name3() != " MG" ) continue;
		if ( mg_res_list_.size() > 0 &&
				!mg_res_list_.has_value( i  ) ) continue;
		hydrate_magnesium( pose, i );
	}
}

///////////////////////////////////////////
void
MgHydrater::hydrate_magnesium( pose::Pose & pose, Size const i ) {
	using namespace core::id;
	using namespace core::pose;

	if ( excise_mini_pose_ ) {

		Pose pose_full = pose;
		Size const mg_res_in_full( i );
		vector1< Size > slice_res = pdbslice( pose, mg_res_in_full );

		Size const mg_res_in_mini    = slice_res.index( mg_res_in_full );
		Size const nres_original( pose.total_residue() );
		hydrate_magnesium_in_pose( pose, mg_res_in_mini );

		for ( Size ii = 1; ii <= pose_full.residue( mg_res_in_full ).natoms(); ii++ ) {
			pose_full.set_xyz( AtomID( ii, mg_res_in_full ), pose.xyz( AtomID( ii, mg_res_in_mini ) ) );
		}
		for ( Size n = nres_original + 1; n <= pose.total_residue(); n++ ) {
			pose_full.append_residue_by_jump( pose.residue( n ), i );
		}
		update_numbers_in_pdb_info( pose_full );

		pose = pose_full;
	} else {
		hydrate_magnesium_in_pose( pose, i );
	}
}


///////////////////////////////////////////
void
MgHydrater::hydrate_magnesium_in_pose( pose::Pose & pose, Size const i,
	bool force_full_shell /* = true */ ) {

	using namespace core::conformation;
	using namespace core::id;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::magnesium;
	using namespace numeric;
	using numeric::conversions::degrees;

	ScoreFunctionOP scorefxn = get_mg_scorefxn();

	Residue const & rsd_i = pose.residue( i );
	runtime_assert( rsd_i.name3() == " MG" );
	Vector xyz_mg = rsd_i.xyz( "MG  " );
	TR << "Setting up Mg(2+) " << i << " which has PDB number: " << pose.pdb_info()->number( i ) << std::endl;

	vector1< AtomID > nbr_atom_ids;
	Distance NBR_CUTOFF( 6.0 ); // 2.1 Mg-HOH, 4.0 HOH-other atom.
	for ( Size j = 1; j <= pose.total_residue(); j++ ) {
		if ( i == j ) continue; // no need to count Mg as a nbr.
		Residue const & rsd_j = pose.residue( j );
		for ( Size jj = 1; jj <= rsd_j.nheavyatoms(); jj++ ) {
			if ( ( rsd_j.xyz( jj ) - xyz_mg ).length() < NBR_CUTOFF ) {
				nbr_atom_ids.push_back( AtomID( jj, j ) );
			}
		}
	}
	Pose pose_orig( pose ), best_pose;
	core::Real best_score( 0.0 );
	bool init( false );

	vector1< Matrix > rotation_matrices;
	if ( use_fast_frame_heuristic_ ) {
		// get ligands
		vector1< AtomID > ligands;
		Distance LIGAND_DIST_CUTOFF( 3.0 );
		for ( Size n = 1; n <= nbr_atom_ids.size(); n++ ) {
			AtomID const & nbr_atom_id( nbr_atom_ids[ n ] );
			if ( pose.residue( nbr_atom_id.rsd() ).atom_type( nbr_atom_id.atomno() ).is_acceptor() ) {
				if ( ( pose.xyz( nbr_atom_id ) - xyz_mg ).length() < LIGAND_DIST_CUTOFF ) {
					ligands.push_back( nbr_atom_id );
				}
			}
		}
		if ( ligands.size() == 1 ) {
			// point z towards this sole ligand, and do all rotations around this Mg-ligand z-axis.
			core::Real const rotstep( 5.0 );
			Vector vec = ( pose.xyz( ligands[1] ) - xyz_mg ).normalized();
			core::Real const theta  = degrees( std::acos( vec.z() ) );
			core::Real const psi    = degrees( std::atan2( vec.x(), vec.y() ) );
			for ( core::Real phi = 0.0; phi <= 90.0; phi += rotstep ) {
				EulerAngles<core::Real> euler_angles( EulerAngles<core::Real>::from_degrees( phi, psi, theta ) );
				rotation_matrices.push_back( euler_angles.to_rotation_matrix().transpose() /* weird convention */ );
			}
		} else if ( ligands.size() > 1 ) {
			// try to align axes to ligands.
			for ( Size n = 1; n <= ligands.size(); n++ ) {
				AtomID const & primary_ligand = ligands[ n ];
				for ( Size m = 1; m <= ligands.size(); m++ ) {
					if ( m == n ) continue;
					AtomID const & secondary_ligand = ligands[ m ];
					rotation_matrices.push_back( set_frame( xyz_mg, pose.xyz( primary_ligand ), pose.xyz( secondary_ligand ) ) );
				}
			}
		}
		// if no RNA ligands, will assume hexahydrate, and do uniform rotation sampling over all frames... see next.
	}

	if ( rotation_matrices.size() == 0 ) {
		if ( urs_ == 0 ) urs_ = get_octahedral_uniform_rotation_sampler( 15.0, true );
		for ( core::Size ctr=1; ctr<=urs_->nrots() ; ++ctr ) {
			Matrix R;
			urs_->get(ctr, R);
			rotation_matrices.push_back( R );
		}
	}

	// do two passes, once looking for full shells, then again allowing for incomplete shells.
	for ( Size pass = 1; pass <= 2; pass++ ) {

		if ( pass == 2 && verbose_ ) TR << "WARNING! Doing another hydration pass, allowing for incomplete coordination shells." << std::endl;
		bool force_full_shell_in_pass = ( pass == 1) && force_full_shell;

		for ( core::Size ctr=1; ctr<=rotation_matrices.size() ; ++ctr ) {
			pose = pose_orig;
			Matrix const & R = rotation_matrices[ ctr ];
			bool full_shell = hydrate_magnesium_with_orbital_frame( pose, i, nbr_atom_ids, R, force_full_shell_in_pass );
			if ( ( pass == 1) && !full_shell ) continue;
			core::Real score = ( *scorefxn )( pose );
			if ( verbose_ ) {
				TR << "checked: " << I( 4, ctr) << " out of " << I( 4, rotation_matrices.size()) <<
					"; adding " << pose.total_residue() - pose_orig.total_residue() << " waters gives: " << F(8,3,score) << std::endl;
			}
			//  scorefxn->show( pose );
			//    pose.dump_pdb( "CHECK_" + ObjexxFCL::lead_zero_string_of( ctr, 3 ) + ".pdb" );
			if ( (score < best_score) || !init ) {
				best_score = score;
				best_pose = pose;
				init = true;
			}
			if ( nbr_atom_ids.size() == 0 /*pose_orig.total_residue() == 1*/ ) break; // lone Mg(2+) -- special case.
		}
		if ( init ) break; // found a full_shell
	}

	pose = best_pose;
}


///////////////////////////////////////////
bool
MgHydrater::hydrate_magnesium_with_orbital_frame( pose::Pose & pose,
	Size const i,
	vector1< core::id::AtomID > const & nbr_atom_ids,
	numeric::xyzMatrix< core::Real > const & R,
	bool force_full_shell /*= true*/ ) const {
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::chemical;

	Distance CLASH_CUTOFF_DEFAULT( 2.8 );
	Distance CLASH_CUTOFF_HBOND( 1.5 );
	vector1< Vector > vec_std( 6 ), vec( 6 );
	vec_std[ 1 ] = Vector( 1.0, 0.0, 0.0 );
	vec_std[ 2 ] = Vector( 0.0, 1.0, 0.0 );
	vec_std[ 3 ] = Vector( 0.0, 0.0, 1.0 );
	vec_std[ 4 ] = Vector( -1.0, 0.0, 0.0 );
	vec_std[ 5 ] = Vector( 0.0, -1.0, 0.0 );
	vec_std[ 6 ] = Vector( 0.0, 0.0, -1.0 );
	for ( Size n = 1; n <= 6; n++ ) {
		vec[ n ] = R * vec_std[ n ];
	}

	Residue const & rsd_i = pose.residue( i );
	Vector xyz_mg = rsd_i.xyz( "MG  " );

	vector1< Size > water_res;
	for ( Size n = 1; n <= 6; n++ ) {
		// place V-atoms.
		Vector const v_atom_xyz = xyz_mg + MG_V_DISTANCE * vec[ n ];
		pose.set_xyz( NamedAtomID( "V"+I(1,n), i ), v_atom_xyz );
	}

	Size const nres_orig = pose.total_residue();
	for ( Size n = 1; n <= 6; n++ ) {
		// place water if no clash.
		Vector const xyz_water = xyz_mg + MG_HOH_DISTANCE * vec[ n ];

		// check for clash.
		bool do_instantiate_water = true;
		bool other_possible_ligand = false;
		for ( Size k = 1; k <= nbr_atom_ids.size(); k++ ) {
			AtomID nbr_atom_id = nbr_atom_ids[ k ];
			//if ( nbr_atom_id.rsd() == i ) continue; // don't worry about clash with Mg(2+).
			if ( pose.residue( nbr_atom_id.rsd() ).name3() == " MG" ) continue; // don't worry about clash with Mg(2+).
			AtomType atom_type = pose.residue( nbr_atom_id.rsd() ).atom_type( nbr_atom_id.atomno() );

			// don't worry about clash with possible H-bond donor/acceptor.
			Distance const clash_cutoff = ( atom_type.is_acceptor() || atom_type.is_donor() ) ? CLASH_CUTOFF_HBOND : CLASH_CUTOFF_DEFAULT;
			if ( ( pose.xyz( nbr_atom_id ) - xyz_water ).length() < clash_cutoff ) {
				do_instantiate_water = false;
				if ( atom_type.is_acceptor() ) { // can potentially complete the shell.
					other_possible_ligand = true;
				}
				//     TR << "found clash of water " << n << " with  " << pose.pdb_info()->number( nbr_atom_id.rsd() ) << " " << pose.residue( nbr_atom_id.rsd() ).atom_name( nbr_atom_id.atomno() ) << " dist " << ( pose.xyz( nbr_atom_id ) - xyz_water ).length()   << " other_possible_ligand " << other_possible_ligand << std::endl;
			}
		}
		if ( do_instantiate_water ) {
			instantiate_water_at_octahedral_vertex( pose, i, n );
			water_res.push_back( pose.total_residue() );
		} else {
			if ( !other_possible_ligand && force_full_shell ) return false;
		}
	}

	if ( nres_orig == 1 ) return true;
	for ( Size k = 1; k <= water_res.size(); k++ ) {
		mg_water_hydrogen_packer_->apply( pose, std::make_pair( i, water_res[ k ] ) );
	}
	return true;
}

///////////////////////////////////////////
Matrix
MgHydrater::set_frame( Vector const & orig, Vector const & xyz1, Vector const & xyz2 ) const {
	Vector const x = ( xyz1 - orig ).normalized();
	Vector y = ( xyz2 - orig ).normalized();
	Vector z = cross( x, y );
	y = cross( z, x );
	return Matrix::cols( x, y, z );
}

} //magnesium
} //protocols
