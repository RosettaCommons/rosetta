// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/DRRAFTER/fit_helix_in_map.cc
/// @brief Fit a helix in a density map at a specific position
/// @author Kalli Kappel

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/rms_util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/init/init.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/PDBInfo.hh>
#include <utility/vector1.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.tmpl.hh>
#include <protocols/rigid/RigidBodyMover.hh>

// C++ headers
#include <iostream>
#include <string>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>

OPT_KEY( RealVector, end_node_in_map )
OPT_KEY( RealVector, next_node_in_map )
OPT_KEY( String, output_pdb_name )
OPT_KEY( String, probe_helix_pdb )
OPT_KEY( Boolean, fit_helix )
OPT_KEY( Boolean, align )
OPT_KEY( Boolean, fit_reverse )
OPT_KEY( Integer, nrepeats )
OPT_KEY( String, helix_for_alignment )
OPT_KEY( String, fit_probe_helix )
OPT_KEY( String, output_pdb_name_aligned )

static basic::Tracer TR( "apps.public.DRRAFTER.fit_helix_in_map" );

///////////////////////////////////////////////////////////////////////////////
void
fit_helix_in_map() {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace core;
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace core::id;

	if ( ( option[ end_node_in_map ]().size() ) != 3 )  {
		TR << "ERROR: end_node_in_map should be an (x, y, z) position in the map" << std::endl;
		exit( 0 );
	}

	if ( ( option[ next_node_in_map ]().size() ) != 3 )  {
		TR << "ERROR: next_node_in_map should be an (x, y, z) position in the map" << std::endl;
		exit( 0 );
	}

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// Create a pose from the input probe helix
	PoseInputStreamOP input_probe_helix;
	input_probe_helix = PoseInputStreamOP( new PDBPoseInputStream( option[ probe_helix_pdb ]() ) );

	pose::Pose probe_helix_pose;
	input_probe_helix->fill_pose( probe_helix_pose, *rsd_set );

	// get the base pairs in the probe helix
	utility::vector1< std::pair< core::Size, core::Size >  > probe_bp_list;
	core::pose::rna::get_base_pairing_list( probe_helix_pose, probe_bp_list );

	// get the midpoint between the end node and next node
	utility::vector1< Real > const end_node_in_map_input = option[ end_node_in_map ]();
	utility::vector1< Real > const next_node_in_map_input = option[ next_node_in_map ]();
	numeric::xyzVector< Real > midpoint = numeric::xyzVector< Real > (
		(end_node_in_map_input[1] + next_node_in_map_input[1]) / 2.0 + 0.5,
		(end_node_in_map_input[2] + next_node_in_map_input[2]) / 2.0 + 0.5,
		(end_node_in_map_input[3] + next_node_in_map_input[3]) / 2.0 + 0.5);
	numeric::xyzVector< Real > end_node_in_map_vector = numeric::xyzVector< Real > (
		end_node_in_map_input[1], end_node_in_map_input[2], end_node_in_map_input[3] );
	numeric::xyzVector< Real > next_node_in_map_vector = numeric::xyzVector< Real > (
		next_node_in_map_input[1], next_node_in_map_input[2], next_node_in_map_input[3] );

	std::string sequence_for_pose = "AAA";
	pose::Pose pose_align_points;
	pose::make_pose_from_sequence( pose_align_points, sequence_for_pose, *rsd_set );

	pose_align_points.set_xyz( AtomID( pose_align_points.residue( 1 ).atom_index("CA"), 1), end_node_in_map_vector );
	pose_align_points.set_xyz( AtomID( pose_align_points.residue( 2 ).atom_index("CA"), 2), midpoint );
	pose_align_points.set_xyz( AtomID( pose_align_points.residue( 3 ).atom_index("CA"), 3), next_node_in_map_vector );
	pose_align_points.dump_pdb( "pose_align_points.pdb" );


	// set up the atom map
	AtomID_Map< id::AtomID > atom_map;
	pose::initialize_atomid_map( atom_map, probe_helix_pose, AtomID::BOGUS_ATOM_ID() );
	Size resnum1 = probe_bp_list[ probe_bp_list.size()].first;
	Size resnum2 = probe_bp_list[ probe_bp_list.size()].second;
	AtomID atom_id1 = AtomID( probe_helix_pose.residue( resnum1 ).atom_index( " C1'" ), resnum1 );
	AtomID atom_id2 = AtomID( probe_helix_pose.residue( resnum2 ).atom_index( " C1'" ), resnum2 );
	AtomID atom_id3 = AtomID( pose_align_points.residue( 1 ).atom_index( " CA " ), 1 );
	atom_map[ atom_id1 ] = atom_id3;
	atom_map[ atom_id2 ] = atom_id3;


	Size resnum3 = probe_bp_list[ probe_bp_list.size()/2 + 1].first;
	Size resnum4 = probe_bp_list[ probe_bp_list.size()/2 + 1].second;
	AtomID atom_id4 = AtomID( probe_helix_pose.residue( resnum3 ).atom_index( " C1'" ), resnum3 );
	AtomID atom_id5 = AtomID( probe_helix_pose.residue( resnum4 ).atom_index( " C1'" ), resnum4 );
	AtomID atom_id6 = AtomID( pose_align_points.residue( 2 ).atom_index( " CA " ), 2 );
	atom_map[ atom_id4 ] = atom_id6;
	atom_map[ atom_id5 ] = atom_id6;

	Size resnum5 = probe_bp_list[ 1].first;
	Size resnum6 = probe_bp_list[ 1].second;
	AtomID atom_id7 = AtomID( probe_helix_pose.residue( resnum5 ).atom_index( " C1'" ), resnum5 );
	AtomID atom_id8 = AtomID( probe_helix_pose.residue( resnum6 ).atom_index( " C1'" ), resnum6 );
	AtomID atom_id9 = AtomID( pose_align_points.residue( 3 ).atom_index( " CA " ), 3 );
	atom_map[ atom_id7 ] = atom_id9;
	atom_map[ atom_id8 ] = atom_id9;

	scoring::superimpose_pose( probe_helix_pose, pose_align_points, atom_map );

	// now do some optimization
	scoring::ScoreFunctionOP sfxn_dens( new scoring::ScoreFunction );
	sfxn_dens->set_weight( scoring::elec_dens_fast, 1.0 );
	addVirtualResAsRoot( probe_helix_pose );
	//Real score = (*sfxn_dens)( probe_helix_pose );

	Size first_bp_res1 = probe_bp_list[1].first;
	Size first_bp_res2 = probe_bp_list[1].second;
	numeric::xyzVector< Real > first_xyz1 = probe_helix_pose.residue( first_bp_res1 ).atom( "C1'" ).xyz();
	numeric::xyzVector< Real > first_xyz2 = probe_helix_pose.residue( first_bp_res2 ).atom( "C1'" ).xyz();
	numeric::xyzVector< Real > first_bp_sum = first_xyz1 + first_xyz2;
	numeric::xyzVector< Real > first_bp_midpoint = numeric::xyzVector< Real > ( first_bp_sum.x()/2., first_bp_sum.y()/2., first_bp_sum.z()/2. );


	Size last_bp_res1 = probe_bp_list[probe_bp_list.size()].first;
	Size last_bp_res2 = probe_bp_list[probe_bp_list.size()].second;
	numeric::xyzVector< Real > last_xyz1 = probe_helix_pose.residue( last_bp_res1 ).atom( "C1'" ).xyz();
	numeric::xyzVector< Real > last_xyz2 = probe_helix_pose.residue( last_bp_res2 ).atom( "C1'" ).xyz();
	numeric::xyzVector< Real > last_bp_sum = last_xyz1 + last_xyz2;
	numeric::xyzVector< Real > last_bp_midpoint = numeric::xyzVector< Real > ( last_bp_sum.x()/2., last_bp_sum.y()/2., last_bp_sum.z()/2. );

	numeric::xyzVector< Real > rot_center_sum = first_bp_midpoint + last_bp_midpoint;
	numeric::xyzVector< Real > rot_center = numeric::xyzVector< Real > ( rot_center_sum.x()/2., rot_center_sum.y()/2., rot_center_sum.z()/2. );

	protocols::rigid::RigidBodySpinMoverOP spin_mover = protocols::rigid::RigidBodySpinMoverOP(
		new protocols::rigid::RigidBodySpinMover( 1 ) );
	spin_mover->spin_axis( first_bp_midpoint - last_bp_midpoint );
	spin_mover->rot_center( rot_center );

	Vector zero_spin_axis( 0.0, 0.0, 0.0 );
	protocols::rigid::RigidBodyTransMoverOP trans_mover = protocols::rigid::RigidBodyTransMoverOP(
		new protocols::rigid::RigidBodyTransMover( zero_spin_axis, 1, false ));
	trans_mover->trans_axis( first_bp_midpoint - last_bp_midpoint );
	trans_mover->step_size( 1.0 );

	protocols::rigid::RigidBodyPerturbMoverOP rigid_mover = protocols::rigid::RigidBodyPerturbMoverOP(
		new protocols::rigid::RigidBodyPerturbMover( 1, 5., 5. ));

	Size resnum_check = probe_bp_list[ probe_bp_list.size() ].first;

	Real start_distance = (probe_helix_pose.residue( resnum_check ).atom( "C1'" ).xyz() - end_node_in_map_vector ).length();
	Real start_distance_to_next_node = (probe_helix_pose.residue( resnum_check ).atom( "C1'" ).xyz() - next_node_in_map_vector ).length();

	Real start_score = (*sfxn_dens)(probe_helix_pose );
	pose::Pose min_scoring;
	Real min_score = (*sfxn_dens)( probe_helix_pose );
	min_scoring = *(probe_helix_pose.clone());
	pose::Pose start_probe_helix;
	start_probe_helix = *(probe_helix_pose.clone());

	Real kT = 10.;
	probe_helix_pose.dump_pdb( "before_opt.pdb" );

	//Size repeats = 10;
	Size repeats = option[ nrepeats ]();
	for ( Size repeat = 1; repeat <= repeats; ++repeat ) {
		Real current_score = start_score;
		probe_helix_pose = *(start_probe_helix.clone());
		pose::Pose try_move = *(start_probe_helix.clone());
		for ( Size i = 1; i <= 5000; ++i ) {
			//try_move.dump_pdb( "try_move_before_"+ ObjexxFCL::string_of( i )+".pdb" );
			if ( (i % 3) == 0 ) {
				spin_mover->spin_mag( numeric::random::rg().uniform()*180 );
				spin_mover->apply( try_move );
			} else if ( ( i % 3) == 1 ) {
				rigid_mover->apply( try_move );
			} else {
				trans_mover->step_size( numeric::random::rg().gaussian()*1.0 );
				trans_mover->apply( try_move );
			}
			//try_move.dump_pdb( "try_move_after_"+ ObjexxFCL::string_of( i )+".pdb" );
			Real distance = (try_move.residue( resnum_check ).atom( "C1'" ).xyz() - end_node_in_map_vector ).length();
			Real distance_to_next_node = ( try_move.residue( resnum_check ).atom( "C1'" ).xyz() - next_node_in_map_vector ).length();
			Real score = (*sfxn_dens)(try_move );
			if ( ( distance < 15.0 || distance < start_distance ) && (distance_to_next_node >= start_distance_to_next_node - 8.0 ) && ( ( score < current_score) || (numeric::random::rg().uniform() < exp( -1.0* (score - current_score )/kT )) ) ) {
				probe_helix_pose = *(try_move.clone() );
				current_score = score;
				if ( current_score < min_score ) {
					min_score = current_score;
					//std::cout << "min_score " << min_score << std::endl;
					//std::cout << "start_score right after " << start_score << std::endl;
					min_scoring = *(probe_helix_pose.clone() );
				}
			} else {
				try_move = *(probe_helix_pose.clone() );
			}
		}
	}

	probe_helix_pose = *(min_scoring.clone() );

	probe_helix_pose.dump_pdb( option[ output_pdb_name ]() );

}

///////////////////////////////////////////////////////////////
void
align_helix() {
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	using namespace core;
	using namespace core::chemical;
	using namespace core::import_pose::pose_stream;
	using namespace core::id;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// Create a pose from the input probe helix
	PoseInputStreamOP input_helix_to_align;
	input_helix_to_align = PoseInputStreamOP( new PDBPoseInputStream( option[ helix_for_alignment ]() ) );

	pose::Pose helix_to_align_pose;
	input_helix_to_align->fill_pose( helix_to_align_pose, *rsd_set );

	PoseInputStreamOP input_fit_probe_helix;
	input_fit_probe_helix = PoseInputStreamOP( new PDBPoseInputStream( option[ fit_probe_helix ]() ) );

	pose::Pose fit_probe_helix_pose;
	input_fit_probe_helix->fill_pose( fit_probe_helix_pose, *rsd_set );

	utility::vector1< std::pair< core::Size, core::Size >  > bp_list;
	core::pose::rna::get_base_pairing_list( helix_to_align_pose, bp_list );

	utility::vector1< std::pair< core::Size, core::Size >  > probe_bp_list;
	core::pose::rna::get_base_pairing_list( fit_probe_helix_pose, probe_bp_list );

	AtomID_Map< id::AtomID > atom_map;
	pose::initialize_atomid_map( atom_map, helix_to_align_pose, AtomID::BOGUS_ATOM_ID() );

	Size nresidues = bp_list.size() * 2;
	Size num_bps_to_match( 3 );

	if ( nresidues == 4 ) {
		num_bps_to_match = 2;
	} else if ( nresidues == 2 ) {
		num_bps_to_match = 1;
	}

	bool fit_rev = option[ fit_reverse ]();

	for ( Size i = 0; i < num_bps_to_match; ++i ) {
		Size resnum1 = bp_list[ bp_list.size() - i ].first;
		Size resnum2 = bp_list[ bp_list.size() - i ].second;
		Size resnum1_probe, resnum2_probe;
		if ( !fit_rev ) {
			resnum1_probe = probe_bp_list[ probe_bp_list.size() - i ].first;
			resnum2_probe = probe_bp_list[ probe_bp_list.size() - i ].second;
		} else {
			resnum1_probe = probe_bp_list[ 1 + i ].second;
			resnum2_probe = probe_bp_list[ 1 + i ].first;
		}
		AtomID atom_id1 = AtomID( helix_to_align_pose.residue( resnum1 ).atom_index( " C1'" ), resnum1 );
		AtomID atom_id2 = AtomID( helix_to_align_pose.residue( resnum2 ).atom_index( " C1'" ), resnum2 );
		AtomID atom_id3 = AtomID( fit_probe_helix_pose.residue( resnum1_probe ).atom_index( " C1'" ), resnum1_probe );
		AtomID atom_id4 = AtomID( fit_probe_helix_pose.residue( resnum2_probe ).atom_index( " C1'" ), resnum2_probe );

		atom_map[ atom_id1 ] = atom_id3;
		atom_map[ atom_id2 ] = atom_id4;

	}

	if ( num_bps_to_match == 2 || num_bps_to_match == 1 ) {
		Size resnum1 = bp_list[ bp_list.size() ].first;
		Size resnum2 = bp_list[ bp_list.size() ].second;
		Size resnum1_probe, resnum2_probe;
		if ( !fit_rev ) {
			resnum1_probe = probe_bp_list[ probe_bp_list.size() ].first;
			resnum2_probe = probe_bp_list[ probe_bp_list.size() ].second;
		} else {
			resnum1_probe = probe_bp_list[ 1 ].second;
			resnum2_probe = probe_bp_list[ 1 ].first;
		}
		AtomID atom_id1 = AtomID( helix_to_align_pose.residue( resnum1 ).atom_index( " C3'" ), resnum1 );
		AtomID atom_id2 = AtomID( helix_to_align_pose.residue( resnum2 ).atom_index( " C3'" ), resnum2 );
		AtomID atom_id3 = AtomID( fit_probe_helix_pose.residue( resnum1_probe ).atom_index( " C3'" ), resnum1_probe );
		AtomID atom_id4 = AtomID( fit_probe_helix_pose.residue( resnum2_probe ).atom_index( " C3'" ), resnum2_probe );

		atom_map[ atom_id1 ] = atom_id3;
		atom_map[ atom_id2 ] = atom_id4;
	}

	scoring::superimpose_pose( helix_to_align_pose, fit_probe_helix_pose, atom_map );

	helix_to_align_pose.dump_pdb( option[ output_pdb_name_aligned ]() );

}
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options::OptionKeys;
	using namespace basic::options;
	if ( option[ fit_helix ]() ) {
		fit_helix_in_map();
	}
	if ( option[ align ]() ) {
		align_helix();
	}

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace core::scoring;
		using namespace basic::options::OptionKeys;

		option.add_relevant( in::file::s );
		NEW_OPT( end_node_in_map, "xyz coordinates of the end node in map", utility::vector1<core::Real>() );
		NEW_OPT( next_node_in_map, "xyz coordinates of the end node in map", utility::vector1<core::Real>() );
		NEW_OPT( output_pdb_name, "name of the output pdb file", "" );
		NEW_OPT( probe_helix_pdb, "name of the input pdb file (a helix to fit in the density map", "" );
		NEW_OPT( nrepeats, "amount of optimization to do", 10 );
		NEW_OPT( fit_helix, "fit a helix in a density map given an end node", false );
		NEW_OPT( align, "align a helix to a probe helix", false );
		NEW_OPT( fit_reverse, "align the helix in reverse", false );
		NEW_OPT( helix_for_alignment, "the helix that should get aligned", "" );
		NEW_OPT( fit_probe_helix, "the probe helix that was already fit into the map", "");
		NEW_OPT( output_pdb_name_aligned, "name of the output pdb file (aligned helix)", "");

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
