// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/MPDockingSetupMover.fwd.hh
/// @brief      Reads in 2 poses and 2 spanfiles, concatenates them, and
///    prints them out
///    CURRENTLY ONLY WORKS FOR 2 POSES!!!
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (10/16/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingSetupMover_cc
#define INCLUDED_protocols_docking_membrane_MPDockingSetupMover_cc

// Unit Headers
#include <protocols/docking/membrane/MPDockingSetupMover.hh>
#include <protocols/docking/membrane/MPDockingSetupMoverCreator.hh>
#include <protocols/moves/Mover.hh>
//
//// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/util.hh>
#include <protocols/membrane/util.hh>
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/FlipMover.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMover.hh>
#include <protocols/membrane/TransformIntoMembraneMover.hh>
#include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/simple_moves/AddChainMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//// Utility Headers
#include <protocols/toolbox/superimpose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.docking.membrane.MPDockingSetupMover" );

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace protocols::moves;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
MPDockingSetupMover::MPDockingSetupMover() :
	protocols::moves::Mover(),
	poses_(),
	spanfiles_()
{}

/// @brief Copy Constructor
MPDockingSetupMover::MPDockingSetupMover( MPDockingSetupMover const & /*src*/ ) = default;

/// @brief Destructor
MPDockingSetupMover::~MPDockingSetupMover() = default;

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPDockingSetupMover::clone() const {
	return ( protocols::moves::MoverOP( new MPDockingSetupMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPDockingSetupMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPDockingSetupMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPDockingSetupMover::parse_my_tag(
	utility::tag::TagCOP,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
)
{}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP MPDockingSetupMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new MPDockingSetupMover );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP MPDockingSetupMoverCreator::keyname() const {
// XRW TEMP  return MPDockingSetupMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP MPDockingSetupMover::mover_name() {
// XRW TEMP  return "MPDockingSetupMover";
// XRW TEMP }

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPDockingSetupMover)
// XRW TEMP std::string
// XRW TEMP MPDockingSetupMover::get_name() const {
// XRW TEMP  return "MPDockingSetupMover";
// XRW TEMP }

////////////////////////////////////////////////////////////////////////////////

/// @brief Add membrane components to the pose, then dock proteins along
///   the flexible jump
void MPDockingSetupMover::apply( Pose & pose ) {

	using namespace core::pose;
	using namespace core::conformation::membrane;
	using namespace protocols::docking;
	using namespace protocols::membrane;
	using namespace protocols::simple_moves;

	TR << "mpdocking_setup" << std::endl;
	TR << "PLEASE MAKE SURE YOUR PDBS AND SPANFILES CORRESPOND TO EACH OTHER IN THE INPUT VECTORS!!!" << std::endl;

	// read in poses and spanfiles
	MPDockingSetupMover::read_poses();
	MPDockingSetupMover::read_spanfiles();
	MPDockingSetupMover::init_from_cmd();

	// get poses and topologies from vectors
	PoseOP pose1( poses_[1] );
	PoseOP pose2( poses_[2] );
	std::string spanfile1 = spanfiles_[1];
	std::string spanfile2 = spanfiles_[2];

	// what happens here:
	// we have two poses, that we need to concatenate into a single one;
	// since removing a membrane residue isn't possible right now, we will
	// 1) make copies of both poses
	// 2) transform both original poses into the membrane (via adding a MEM residue)
	// 3) superimpose the copies onto the ones in the membrane
	// 4) concatenate the copies into a single pose, which is copy of pose1
	// 5) concatenate the spanning topologies
	// 6) create a membrane pose (adding MEM) out of the main pose (copy of pose1)

	// deepcopies!!!
	Pose pose1_cp = *poses_[1];
	Pose pose2_cp = *poses_[2];

	// set normals and centers
	Vector normal  (  0, 0, 1); // overall
	Vector center  (  0, 0, 0); // overall

	Vector center1 (-10, 0, 0); // for pose 1
	Vector center2 ( 10, 0, 0); // for pose 2

	// add MEM and transform both poses into the membrane, optionally optimize membrane position
	MPDockingSetupMover::transform_pose_into_membrane( *pose1, center1, normal, spanfile1, 1 );
	MPDockingSetupMover::transform_pose_into_membrane( *pose2, center2, normal, spanfile2, 2 );

	// superimpose pose1_copy with pose1
	protocols::toolbox::CA_superimpose( *pose1, pose1_cp );
	pose1_cp.dump_pdb("pose1superposed.pdb");

	// superimpose pose2_copy with pose2
	protocols::toolbox::CA_superimpose( *pose2, pose2_cp );
	pose2_cp.dump_pdb("pose2superposed.pdb");

	// get topologies of poses
	SpanningTopologyOP topo1 = pose1->conformation().membrane_info()->spanning_topology();
	SpanningTopologyOP topo2 = pose2->conformation().membrane_info()->spanning_topology();

	// concatenate topologies
	SpanningTopologyOP total_topo( topo1 );
	total_topo->concatenate_topology( *topo2 );
	total_topo->write_spanfile("mpdocking_setup_out.span");
	TR << "=== WROTE SPANFILE ===: mpdocking_setup_out.span" << std::endl;

	// concatenate pose2_cp to end of pose1_cp, new chain yes
	// pose1_cp is new master pose!!!
	// this NEEDS to happen before the MEM virtual residue is attached!
	core::pose::append_pose_to_pose( pose1_cp, pose2_cp, true );
	bool renumber = renumber_pdbinfo_based_on_conf_chains( pose1_cp, true, false, false, true );
	TR << "renumber successful? " << renumber << std::endl;

	// add membrane, setup membrane object using known topology
	AddMembraneMoverOP add_membrane3( new AddMembraneMover( topo1 ) );
	add_membrane3->apply( pose1_cp );

	// reorder foldtree
	using namespace protocols::membrane;
	reorder_membrane_foldtree( pose1_cp );

	// get axis for sliding together the two poses
	Vector slide_axis = center2 - center1;

	// slide into contact
	DockingSlideIntoContactOP slide2contact( new DockingSlideIntoContact( 1, slide_axis) );
	slide2contact->apply( pose1_cp );

	pose1_cp.dump_pdb("mpdocking_setup_out.pdb");
	TR << "=== DUMPED FINAL POSE ===: mpdocking_setup_out.pdb" << std::endl;

	// copy pose1_cp into the pose we are working with
	pose = pose1_cp;

	// reset foldtree and show final one
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );


}// apply

////////////////////////////////////////////////////////////////////////////////

// read poses
void MPDockingSetupMover::read_poses() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// cry if PDBs not given
	if ( ! option[OptionKeys::in::file::s].user() ) {
		utility_exit_with_message("Please provide two PDB files with flag -in:file:s!");
	}

	// get filenames from option system
	utility::vector1< std::string > pdbs = option[OptionKeys::in::file::s]();

	// put poses into private vector
	for ( Size i = 1; i <= pdbs.size(); ++i ) {
		PoseOP pose = core::import_pose::pose_from_file( pdbs[i] , core::import_pose::PDB_file);
		poses_.push_back( pose );
	}

}// read poses

////////////////////////////////////////////////////////////////////////////////

// read spanfiles
void MPDockingSetupMover::read_spanfiles() {
	spanfiles_ = core::conformation::membrane::spanfile_names();
}// read spanfiles

////////////////////////////////////////////////////////////////////////////////

// initialize from commandline
void MPDockingSetupMover::init_from_cmd() {

	using namespace basic::options;

	optimize1_ = false;
	optimize2_ = false;

	// flags for optimizing partner 1 or 2 (depending on order given)
	if ( option[ OptionKeys::mp::setup::optimize1 ].user() ) {
		optimize1_ = option[ OptionKeys::mp::setup::optimize1 ]();
	}
	if ( option[ OptionKeys::mp::setup::optimize2 ].user() ) {
		optimize2_ = option[ OptionKeys::mp::setup::optimize2 ]();
	}

}// read spanfiles

////////////////////////////////////////////////////////////////////////////////

// transform pose into membrane
void MPDockingSetupMover::transform_pose_into_membrane( Pose & pose, Vector center, Vector normal, std::string spanfile, Size partner ) {

	using namespace protocols::membrane::geometry;
	using namespace protocols::membrane;

	// attach MEM to pose
	AddMembraneMoverOP add_membrane( new AddMembraneMover( spanfile ) );
	add_membrane->apply( pose );

	// reorder foldtree such that MEM is at root
	reorder_membrane_foldtree( pose );

	// optimize membrane position if desired and then transform into membrane
	if ( ( ( optimize1_ == true ) && ( partner == 1 ) ) || ( ( optimize2_ == true ) && ( partner == 2 ) ) ) {

		TR << "Transforming into membrane and optimizing membrane position for partner " << partner << std::endl;
		OptimizeProteinEmbeddingMoverOP opt_emb( new OptimizeProteinEmbeddingMover() );
		opt_emb->apply( pose );

	} else if ( ( ( optimize1_ == false ) && ( partner == 1 ) ) || ( ( optimize2_ == false ) && ( partner == 2 ) ) ) {

		TR << "Transforming partner " << partner << " into the membrane, adding small random tilt angle below 30 degrees." << std::endl;

		// transform pose into membrane, true for keep current embedding-membrane relation
		TransformIntoMembraneMoverOP transform( new TransformIntoMembraneMover( center, normal ) );
		transform->apply( pose );

		// tilt pose slightly with random flip angle, default jump is membrane jump
		FlipMoverOP flip( new FlipMover() );
		flip->set_random_membrane_flip_angle();
		flip->set_range( 30 );
		flip->apply( pose );
	} else {
		TR << "No transformation took place, check your inputs for mp:setup:optimize1 or 2" << std::endl;
	}

}// transform pose into membrane

std::string MPDockingSetupMover::get_name() const {
	return mover_name();
}

std::string MPDockingSetupMover::mover_name() {
	return "MPDockingSetupMover";
}

void MPDockingSetupMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Set up for membrane protein docking.", attlist );
}

std::string MPDockingSetupMoverCreator::keyname() const {
	return MPDockingSetupMover::mover_name();
}

protocols::moves::MoverOP
MPDockingSetupMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MPDockingSetupMover );
}

void MPDockingSetupMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MPDockingSetupMover::provide_xml_schema( xsd );
}



} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingSetupMover_cc
