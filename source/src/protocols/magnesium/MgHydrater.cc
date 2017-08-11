// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/EulerAngles.hh>
#include <ObjexxFCL/format.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.magnesium.MgHydrater" );

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
	use_virtual_waters_as_placeholders_( false ),
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
	use_virtual_waters_as_placeholders_( false ),
	verbose_( false )
{
	mg_water_hydrogen_packer_->set_excise_mini_pose( false );
}

//Destructor
MgHydrater::~MgHydrater()
= default;

///////////////////////////////////////////
void
MgHydrater::apply( pose::Pose & pose ) {

	if ( use_virtual_waters_as_placeholders_ ) setup_virtual_waters_around_magnesiums( pose );

	// find all the magnesiums in here.
	for ( Size i = 1; i <= pose.size(); i++ ) {
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
		fix_fold_tree_in_excised_pose_for_mg_bound_waters( pose, mg_res_in_mini, pose_full, slice_res );
		Size const nres_original( pose.size() );
		hydrate_magnesium_in_pose( pose, mg_res_in_mini );

		if ( use_virtual_waters_as_placeholders_ ) { // swap out mg & waters into full pose
			for ( Size ii = 1; ii <= pose.size(); ii++ ) {
				if ( ii == mg_res_in_mini || pose.residue_type( ii ).aa() == core::chemical::aa_h2o ) {
					pose_full.replace_residue( slice_res[ ii ], *pose.residue( ii ).clone(), false );
				}
			}
		} else {  // new waters were introduced
			for ( Size ii = 1; ii <= pose_full.residue( mg_res_in_full ).natoms(); ii++ ) {
				pose_full.set_xyz( AtomID( ii, mg_res_in_full ), pose.xyz( AtomID( ii, mg_res_in_mini ) ) );
			}
			for ( Size n = nres_original + 1; n <= pose.size(); n++ ) {
				append_mg_bound_water( pose_full, pose.residue( n ), i );
			}
			update_numbers_in_pdb_info( pose_full );
		}

		pose = pose_full;
	} else {
		hydrate_magnesium_in_pose( pose, i );
	}

	update_full_model_info_with_new_waters( pose, use_virtual_waters_as_placeholders_ /* expect_no_new_waters */ );
}


///////////////////////////////////////////
void
MgHydrater::hydrate_magnesium_in_pose( pose::Pose & pose, Size const i,
	bool force_full_shell /* = true */ ) {

	using namespace core::chemical;
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
	if ( verbose_ ) TR << "Setting up Mg(2+) " << i << " which has PDB number: " << pose.pdb_info()->number( i ) << std::endl;

	vector1< AtomID > nbr_atom_ids;
	Distance NBR_CUTOFF( 6.0 ); // 2.1 Mg-HOH, 4.0 HOH-other atom.
	for ( Size j = 1; j <= pose.size(); j++ ) {
		if ( i == j ) continue; // no need to count Mg as a nbr.
		Residue const & rsd_j = pose.residue( j );
		if ( rsd_j.aa() == aa_h2o ) continue; // don't count waters (we're going to move them around)
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
		if ( urs_ == nullptr ) urs_ = get_octahedral_uniform_rotation_sampler( 15.0, true );
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
			Size num_waters( 0 );
			bool const full_shell = hydrate_magnesium_with_orbital_frame( pose, i, nbr_atom_ids, R, force_full_shell_in_pass, num_waters );
			if ( ( pass == 1) && !full_shell ) continue;
			core::Real score = ( *scorefxn )( pose );
			if ( verbose_ ) {
				TR << "checked: " << I( 4, ctr) << " out of " << I( 4, rotation_matrices.size()) <<
					"; instantiating " << num_waters << " waters gives: " << F(8,3,score) << std::endl;
			}
			if ( (score < best_score) || !init ) {
				best_score = score;
				best_pose = pose;
				init = true;
			}
			if ( nbr_atom_ids.size() == 0 /*pose_orig.size() == 1*/ ) break; // lone Mg(2+) -- special case.
		}
		if ( init ) break; // found a full_shell
	}

	pose = best_pose;
}

///////////////////////////////////////////
void
MgHydrater::update_full_model_info_with_new_waters( pose::Pose & pose,
	bool const expect_no_new_waters /* = false */ ) {
	using namespace core::pose::full_model_info;

	if ( !full_model_info_defined( pose ) ) return;

	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > res_list = full_model_info.res_list(); // will be updated
	runtime_assert( pose.size() >= res_list.size() );
	Size const num_new_waters( pose.size() - res_list.size() );
	if ( expect_no_new_waters ) runtime_assert( num_new_waters == 0 );
	if ( num_new_waters == 0 ) return;
	for ( Size n = pose.size() - num_new_waters + 1; n <= pose.size(); n++ ) {
		runtime_assert( pose.residue( n ).name3() == "HOH" ); // check any new residues are indeed waters.
	}

	// now check for any water 'slots' in full_model_info
	FullModelParameters const & full_model_parameters = *full_model_info.full_model_parameters();
	std::string const & full_sequence = full_model_parameters.full_sequence();
	std::map< Size, std::string > const & non_standard_residue_map = full_model_parameters.non_standard_residue_map();

	utility::vector1< Size > water_res;
	for ( Size n = 1; n <= full_model_parameters.size(); n++ ) {
		if ( full_sequence[ n - 1 ] != 'w' ) continue;
		auto it = non_standard_residue_map.find( n );
		if ( it == non_standard_residue_map.end() ) continue;
		if ( it->second != "HOH" ) continue;
		if ( res_list.has_value( n ) ) continue;
		water_res.push_back( n );
	}

	runtime_assert( water_res.size() >= num_new_waters );
	for ( Size n = 1; n <= num_new_waters; n++ )  res_list.push_back( water_res[ n ] );

	FullModelInfoOP full_model_info_new( full_model_info.clone_info() );
	full_model_info_new->set_res_list( res_list );
	set_full_model_info( pose, full_model_info_new );
	check_full_model_info_OK( pose );
}

///////////////////////////////////////////
bool
MgHydrater::hydrate_magnesium_with_orbital_frame( pose::Pose & pose,
	Size const i,
	vector1< core::id::AtomID > const & nbr_atom_ids,
	numeric::xyzMatrix< core::Real > const & R,
	bool force_full_shell,
	Size & num_waters ) const {
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

	Size const nres_orig = pose.size();
	bool full_shell( true );
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

		bool const replace_residue( use_virtual_waters_as_placeholders_ );
		if ( do_instantiate_water ) {
			Size instantiated_water_res = instantiate_water_at_octahedral_vertex( pose, i, n, MG_HOH_DISTANCE, replace_residue );
			water_res.push_back( instantiated_water_res );
		} else {
			if ( use_virtual_waters_as_placeholders_ ) instantiate_water_at_octahedral_vertex( pose, i, n, MG_HOH_DISTANCE, replace_residue, true /*virtual_water*/ );
			if ( !other_possible_ligand && force_full_shell ) {
				full_shell = false;
			}
		}
	}
	num_waters = water_res.size();
	if ( force_full_shell && !full_shell ) return false;

	if ( nres_orig == 1 ) return true;
	for ( Size k = 1; k <= num_waters; k++ ) {
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

////////////////////////////////////////////////////////////////
void
fix_water_jump( pose::Pose & pose, Size const & parent_res, Size const & water_res )
{
	using namespace core::kinematics;
	FoldTree f( pose.fold_tree() );
	Size const jump_number = f.get_jump_that_builds_residue( water_res );
	f.slide_jump( jump_number, parent_res, water_res );
	pose.fold_tree( f );
	if ( pose.residue_type( parent_res ).name3() == " MG" ) {
		update_jump_atoms_for_mg_bound_water( pose, water_res );
	}
}

////////////////////////////////////////////////////////////////
// In Mg-hydration code with virtual waters, the fold-tree helps
// flag the six waters that go with each Magnesium. No other waters
// should be daughters of that magnesium.
void
MgHydrater::fix_fold_tree_in_excised_pose_for_mg_bound_waters(
	pose::Pose & pose, Size const mg_res,
	pose::Pose const & pose_full,
	utility::vector1< Size > const & slice_res ) const
{
	using namespace core::kinematics;
	FoldTree f( pose.fold_tree() );
	Size const mg_res_in_full = slice_res[ mg_res ];

	// slide jumps for water daughters in full fold tree to be daughters in the excised pose
	utility::vector1< core::Size > bound_waters_in_full_pose = find_bound_waters_that_are_daughters_in_fold_tree( pose_full, mg_res_in_full );
	for ( Size const water_res_in_full : bound_waters_in_full_pose ) {
		Size const water_res_in_mini = slice_res.index( water_res_in_full );
		fix_water_jump( pose, mg_res, water_res_in_mini );
	}

	// and make sure there are no *additional* daughters...
	utility::vector1< core::Size > bound_waters_in_mini_pose = find_bound_waters_that_are_daughters_in_fold_tree( pose, mg_res );
	Size other_res( 0 );
	for ( Size const water_res_in_mini : bound_waters_in_mini_pose ) {
		Size const water_res_in_full = slice_res[ water_res_in_mini ];
		if ( bound_waters_in_full_pose.has_value( water_res_in_full ) ) continue;
		if ( other_res == 0 ) { // need to attach water to some other residue
			for ( Size n = 1; n <= pose.size(); n++ ) {
				if ( n != mg_res && !bound_waters_in_mini_pose.has_value( n ) ) {
					other_res = n; break;
				}
			}
		}
		runtime_assert( other_res != 0 );
		fix_water_jump( pose, other_res, water_res_in_mini );
	}
}


///////////////////////////////////////////////////////
void
MgHydrater::setup_virtual_waters_around_magnesiums( pose::Pose & pose )
{
	vector1< Size > const all_mg_res = get_mg_res( pose );

	// first check -- all Mg(2+) might actually be OK.
	bool all_magnesiums_have_six_waters( true );
	for ( Size n = 1; n <= all_mg_res.size(); n++ ) {
		if ( find_bound_waters_that_are_daughters_in_fold_tree( pose, all_mg_res[ n ] ).size() < 6 ) {
			all_magnesiums_have_six_waters = false; break;
		}
	}
	if ( all_magnesiums_have_six_waters ) return;
	TR << "Defining Mg(2+) water map " << std::endl;

	std::map< Size, vector1< Size > > mg_water_map = define_mg_water_map( pose );

	// get mg-water pairs, properly assigned.
	for ( Size const mg_res : all_mg_res ) {
		// now go through each mg(2+), and fix up fold tree for waters bound to it already.
		vector1< Size > const & water_ligands = mg_water_map[ mg_res ];
		for ( Size const water_res : water_ligands ) {
			if ( pose.fold_tree().get_parent_residue( water_res ) != mg_res ) {
				fix_water_jump( pose, mg_res, water_res );
			}
		}

		// then add additional virtual waters to round out to 6 waters per magnesium
		// for each Mg(2+), append virtual waters as placeholders to get the number up to 6.
		for ( Size n = water_ligands.size() + 1; n <= 6; n++ ) {
			instantiate_water_at_octahedral_vertex( pose, mg_res, n,
				MG_HOH_DISTANCE, false /*replace_residue*/, true /*virtual water*/ );
		}
	}

	update_full_model_info_with_new_waters( pose );
}



} //magnesium
} //protocols
