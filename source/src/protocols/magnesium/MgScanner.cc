// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/MgScanner.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/MgScanner.hh>
#include <protocols/magnesium/MgHydrater.hh>
#include <protocols/magnesium/SampleGrid.hh>
#include <protocols/magnesium/minimize_util.hh>
#include <protocols/magnesium/util.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/magnesium/MgKnowledgeBasedPotential.hh> // ridiculous
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.magnesium.MgScanner" );

using namespace core;
using utility::tools::make_vector1;

//////////////////////////////////////////////////////////////////////////////
//
// Test code for scanning (grid search) of Mg(2+).
//
// Originally developed with goal of
//  allowing metal ion modeling in FARFAR and, eventually helping ERRASER,
//  using a knowledge-based Mg-macromolecule energy function.
//
// In April-May 2015, developed a more sophisticated Mg(2+) modeling
//  scheme involving waters and an 'orbital frame', encoded in mg_modeler.
//
//////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace magnesium {

//Constructor
MgScanner::MgScanner():
	scorefxn_( get_mg_scorefxn() )
{}

//Destructor
MgScanner::~MgScanner()
{}

//////////////////////////////////////////////////////////////////////////////
void
MgScanner::apply( core::pose::Pose & pose )
{
	std::string silent_file_no_cluster = utility::replace_in( silent_file_, ".out", ".no_cluster.out" );
	std::string silent_file_minimize   = utility::replace_in( silent_file_, ".out", ".minimize.out" );

	scan_magnesiums( pose );
	output_mg_to_silent_file( silent_file_no_cluster );

	if ( minimize_ ) {
		minimize_magnesium_and_hydration_shell( pose, mg_poses_,
			utility::vector1<Size>() /*blank means: find the Mg*/, scorefxn_,
			minimize_mg_coord_constraint_distance_ );
		output_mg_to_silent_file( silent_file_minimize );
	}

	cluster_mg();
	output_mg_to_silent_file( silent_file_ );
	output_mg_into_one_PDB( pose );
}

//////////////////////////////////////////////////////////////////////////////
void
MgScanner::scan_magnesiums( core::pose::Pose & pose ) {

	using namespace core::scoring;

	mg_poses_.clear();

	Size const mg_res = pose.total_residue(); // assume there's a Mg at the end of the pose.
	runtime_assert ( pose.residue( mg_res ).name3() == " MG" );

	// FIX THIS. THIS IS SILLY -- SHOULD NOT NEED MG potential AT ALL -- SETUP ATOM NUMBERS.
	scoring::magnesium::MgKnowledgeBasedPotential rna_mg_low_res_potential;
	rna_mg_low_res_potential.setup_info_for_mg_calculation( pose );

	// vdw calculation takes a while, actually -- do a faster calculation...
	ScoreFunctionOP scorefxn_fast( new ScoreFunction );
	scorefxn_fast->set_weight( rna_mg_point,  1.0 );
	scorefxn_fast->set_weight( rna_mg_point_indirect,  1.0 );

	// figure out which mg_positions to sample.
	SampleGrid grid( pose );
	grid.set_xyz_step( xyz_step_ );
	grid.set_input_scan_res( input_scan_res_ );
	grid.set_tether_to_closest_res( tether_to_closest_res_ );
	utility::vector1< Vector > const mg_positions = grid.get_mg_positions( pose );

	pose::Pose pose_fast = pose; // 'scratch' pose for fast checks on mg/ligand binding.
	Vector faraway_position( grid.xmax() + (grid.xmax()-grid.xmin()), grid.ymax(), grid.zmax());
	Real init_score_fast = get_score( pose_fast, faraway_position, scorefxn_fast );

	pose::Pose pose_save = pose; // will not be hydrate -- use to restore
	Real init_score      = get_score( pose, faraway_position, scorefxn_,
		hydrate_, false /*keep waters*/, minimize_during_scoring_ );

	Real const score_cutoff_fast( -6.0 );
	Size count( 0 );
	for ( Size n = 1; n <= mg_positions.size(); n++ ) {

		Vector const & mg_position = mg_positions[ n ];

		// this is a quick filter. another way to accelerate things would be to do a check for
		// steric problems *first* and immediately discard anything with a bad clash.
		Real const score_fast = get_score( pose_fast, mg_position, scorefxn_fast ) - init_score_fast;
		if ( score_fast > score_cutoff_fast ) continue;
		count++;

		// Now really score...
		pose = pose_save;
		Real const score = get_score( pose, mg_position, scorefxn_, hydrate_, true /*keep waters*/, minimize_during_scoring_ ) - init_score;

		if ( score <= score_cut_ ) {
			std::cout << mg_position.x() << " " << mg_position.y() << " " << mg_position.z() << ":  " << score << std::endl;
			core::pose::PoseOP pose_save = pose.clone();
			setPoseExtraScore( *pose_save, "mg_score", score );
			mg_poses_.push_back( pose_save );
		}

		if ( integration_test_ && mg_poses_.size() > 0 ) break; // early exit for testing
	}

	std::cout << "Total positions tested: " << count << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
Real
MgScanner::get_score( pose::Pose & pose,
	Vector const & mg_position,
	scoring::ScoreFunctionCOP scorefxn,
	bool const hydrate_magnesium /* = false */,
	bool const keep_waters /* = false */,
	bool const minimize /* = false */ )
{

	using namespace core::id;

	Real score( 0.0 );

	Size const mg_res( pose.total_residue() ); /* assuming single Mg at end */
	runtime_assert( pose.residue( mg_res ).name3() == " MG" );

	utility::vector1< Vector > new_xyz;
	for ( Size q = 1; q <= pose.residue_type( mg_res ).natoms(); q++ ) {
		new_xyz.push_back( pose.xyz( AtomID( q, mg_res ) ) );
	}
	Vector const mg_position_shift( mg_position - pose.residue( mg_res ).xyz( 1 ) );
	for ( Size q = 1; q <= pose.residue_type( mg_res ).natoms(); q++ ) {
		pose.set_xyz( AtomID( q, mg_res ), new_xyz[ q ] + mg_position_shift);
	}

	if ( hydrate_magnesium ) {
		Pose pose_save;
		if ( !keep_waters ) pose_save = pose;

		protocols::magnesium::MgHydrater mg_hydrater( make_vector1( mg_res ) );
		mg_hydrater.set_excise_mini_pose( false ); // just for ease of visualization
		mg_hydrater.apply( pose );

		if ( minimize ) {
			protocols::magnesium::minimize_magnesium_and_hydration_shell( pose, make_vector1( mg_res ), scorefxn,
				minimize_mg_coord_constraint_distance_ );
		}
		score = (*scorefxn)( pose );

		if ( !keep_waters ) pose = pose_save;
	} else {
		score = (*scorefxn)( pose );
	}
	return score;
}

///////////////////////////////////////////////////////////////////////////////
void
MgScanner::cluster_mg() {

	using namespace core::id;
	// this is such a pain. Need a list to do the sorting...
	std::list< std::pair< Real, pose::PoseOP > > mg_energy_pose_list;
	Size const mg_res = get_unique_mg_res( *mg_poses_[ 1 ] );
	for ( Size n = 1; n <= mg_poses_.size(); n++ ) {
		mg_energy_pose_list.push_back( std::make_pair( getPoseExtraScore( *mg_poses_[ n ], "mg_score" ), mg_poses_[ n ] ) );
	}
	mg_energy_pose_list.sort();

	utility::vector1< pose::PoseOP > mg_poses_cluster;

	Real const CLUSTER_DISTANCE_CUTOFF( 2.0 );
	// iterate...
	for ( std::list< std::pair< Real, pose::PoseOP > >::const_iterator iter = mg_energy_pose_list.begin(),
			end = mg_energy_pose_list.end(); iter != end; ++iter ) {

		pose::PoseOP   mg_pose     = iter->second;
		Vector const   mg_position = mg_pose->xyz( AtomID( 1, mg_res ) );

		bool too_close( false );
		for ( Size n = 1; n <= mg_poses_cluster.size(); n++ ) {
			Vector const mg_position_cluster = mg_poses_cluster[ n ]->xyz( AtomID( 1, mg_res ) );
			if ( ( mg_position - mg_position_cluster ).length() < CLUSTER_DISTANCE_CUTOFF ) {
				too_close = true; break;
			}
		}
		if ( !too_close ) mg_poses_cluster.push_back( mg_pose );

	}

	mg_poses_     = mg_poses_cluster;
}


///////////////////////////////////////////////////////////////////////////////
void
MgScanner::output_mg_to_silent_file( std::string const & silent_file ) {

	using namespace core::id;
	using namespace core::io::silent;

	if ( mg_poses_.size() == 0 ) return;

	SilentFileData silent_file_data;
	Size const mg_res = get_unique_mg_res( *mg_poses_[ 1 ] );
	pose::PoseOP single_mg_pose = get_single_mg_pose();

	for ( Size n = 1; n <= mg_poses_.size(); n++ ) {
		std::string const out_file_tag = "S_" + ObjexxFCL::string_of( n );
		Vector mg_position = mg_poses_[ n ]->xyz( AtomID( 1, mg_res ) );

		BinarySilentStructOP s;
		if ( hydrate_ ) {
			// output the full pose -- will have Mg and HOH.
			s = BinarySilentStructOP( new BinarySilentStruct( *mg_poses_[ n ], out_file_tag ) );
		} else {
			// since no HOH, use a more compact format with just Mg(2+) outputted.
			single_mg_pose->set_xyz( AtomID( 1, 1 ), mg_position );
			s = BinarySilentStructOP( new BinarySilentStruct( *single_mg_pose, out_file_tag ) );
			s->energies_from_pose( *mg_poses_[ n ] );
		}
		s->add_energy( "rms",      distance_to_closest_magnesium( mg_position, *get_native_pose() ) );

		silent_file_data.write_silent_struct( *s, silent_file, false /*score_only*/ );

	}
}

/////////////////////////////////////////////////////////////////////////////
// largely deprecated -- was used to see where a ton of Mg(2+) possible sites
// were located in scans -- only saved Mg(2+) positions since there were no
// waters in early code and because there were a lot of sites.
void
MgScanner::output_mg_into_one_PDB( pose::Pose const & pose )
{
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::io::silent;

	if ( output_pdb_.size() == 0 ) return;

	core::conformation::ResidueOP mg_rsd = get_mg_rsd();
	pose::Pose mg_pose = pose;

	Size count( 0 );
	for ( Size n = 1; n <= mg_poses_.size(); n++ ) {

		if ( getPoseExtraScore( *mg_poses_[ n ], "mg_score" ) > score_cut_PDB_ ) continue;

		Size const mg_res = get_unique_mg_res( *mg_poses_[ n ] );
		Vector mg_position = mg_poses_[ n ]->xyz( AtomID( 1, mg_res ) );
		if ( count > 1 ) mg_pose.append_residue_by_jump( *(mg_rsd->clone()), 1 );
		mg_pose.set_xyz( AtomID( 1, mg_pose.total_residue() ), mg_position );
		count++;
	}

	mg_pose.dump_pdb( output_pdb_ );

}

///////////////////////////////////////////////////////////////////////////////
Size
MgScanner::get_unique_mg_res( pose::Pose const & mg_pose ) {
	Size mg_res( 0 );
	for ( Size n = 1; n <= mg_pose.total_residue(); n++ ) {
		if ( mg_pose.residue( n ).name3() == " MG" ) {
			runtime_assert( mg_res == 0 );
			mg_res = n;
		}
	}
	return mg_res;
}

///////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
MgScanner::get_single_mg_pose() {
	pose::PoseOP single_mg_pose( new Pose );
	single_mg_pose->append_residue_by_bond( *get_mg_rsd() );
	return single_mg_pose;
}

///////////////////////////////////////////////////////////////////////////////
Distance
MgScanner::distance_to_closest_magnesium( Vector const & mg_position,
	pose::Pose const & reference_pose ) {
	Distance min_dist( 0.0 );
	for ( Size m = 1; m <= reference_pose.total_residue(); m++ ) {
		if ( reference_pose.residue(m).name3() != " MG"  ) continue;
		Distance dist = ( reference_pose.residue(m).xyz(1) - mg_position ).length();
		if ( dist < min_dist || min_dist == 0.0 ) min_dist = dist;
	}
	return min_dist;
}



} //magnesium
} //protocols
