// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/membrane/AddMembraneMover.hh
///
/// @brief      Initialize the RosettaMP Framework by adding membrane representations to the pose
/// @details Given a pose, initialize and configure with the RosettaMP framework by taking the
///    following steps:
///     (1) Add a membrane residue to the pose (type MEM)
///      (a) Append residue to the pose or accept a new one
///      (b) Update PDB info to acknowledge the membrane residue
///      (c) Set the MEM residue at the root of the foldtree
///     (2) Initialize transmembrane spanning topology (Object: SpanningTopology)
///     (3) Initialize the MembraneInfo object either with or without the LipidAccInfo object.
///     (4) Set the membrane starting position (either default or based on user input)
///
///    This object does a massive amount of reading from CMD, RosettaScripts or Constructor. If you add
///    a new piece of data - you must modify MembraneInfo, all of these channels for input AND the apply function!
///    If and only if AddMembraneMover is applied to the pose, pose.conformation().is_membrane() MUST return true.
///
///    Last Updated: 7/23/15
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_cc
#define INCLUDED_protocols_membrane_AddMembraneMover_cc

// Unit Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/util.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>

#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneParams.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

#include <protocols/membrane/SetMembranePositionMover.hh>

#include <protocols/moves/DsspMover.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueProperty.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

#include <core/pose/util.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>
#include <map>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.membrane.AddMembraneMover" );

namespace protocols {
namespace membrane {

// Constructors & Setup methods

/// @brief Create a default RosettaMP membrane setup
/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
/// the z axis. Use a defualt lipid type DOPC.
AddMembraneMover::AddMembraneMover() :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( 1 ),
	membrane_rsd_( 0 ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from an existing membrane residue
/// @brief Create a membrane using the position from the existing membrane residue.
/// Use a defualt lipid type DOPC.
AddMembraneMover::AddMembraneMover( core::Size membrane_rsd ) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from a user specified spanfile
/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
/// specified spanfile
AddMembraneMover::AddMembraneMover(
	std::string spanfile,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(std::move( spanfile )),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from a user specified SpanningTopology
/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
/// specified spanning topology
AddMembraneMover::AddMembraneMover(
	core::conformation::membrane::SpanningTopologyOP topology,
	core::Size anchor_rsd,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_( "" ),
	topology_(std::move( topology )),
	lipsfile_(),
	anchor_rsd_( anchor_rsd ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from an existing residue at a specific anchor point
/// @brief Create a membrane using the position from the existing membrane residue. Anchor
/// this residue at the user-specified anchor residue. Use a defualt lipid type DOPC.
AddMembraneMover::AddMembraneMover(
	core::Size anchor_rsd,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_( "" ),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( anchor_rsd ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from a user specified spanfile and lipsfile
/// @brief Create a membrane positioned at the origin (0, 0, 0) and aligned with
/// the z axis. Use a defualt lipid type DOPC. Load spanning topology from the user
/// specified spanfile and lipsfile
AddMembraneMover::AddMembraneMover(
	std::string spanfile,
	std::string lips_acc,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	include_lips_( true ),
	spanfile_(std::move( spanfile )),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_(std::move( lips_acc )),
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from a user specified spanfile and lipsfile
/// @brief Create a membrane positioned at "init_center" and aligned with
/// "init_normal". Use a defualt lipid type DOPC.
AddMembraneMover::AddMembraneMover(
	core::Vector init_center,
	core::Vector init_normal,
	std::string spanfile,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(std::move( spanfile )),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	lipsfile_( "" ),
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
	center_( init_center ),
	normal_( init_normal ),
	thickness_( 15 ),
	steepness_( 10 ),
	user_defined_( false )
{}

/// @brief Create a deep copy of the data in this mover
AddMembraneMover::AddMembraneMover( AddMembraneMover const & ) = default;

/// @brief Destructor
AddMembraneMover::~AddMembraneMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
AddMembraneMover::clone() const {
	return ( protocols::moves::MoverOP( new AddMembraneMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
AddMembraneMover::fresh_instance() const {
	return protocols::moves::MoverOP( new AddMembraneMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
AddMembraneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// Read in include lips option (boolean)
	if ( tag->hasOption( "include_lips" ) ) {
		include_lips_ = tag->getOption< bool >( "include_lips" );
	}

	// Read in spanfile information
	if ( tag->hasOption( "spanfile" ) ) {
		spanfile_ = tag->getOption< std::string >( "spanfile" );
	}

	// Read in lipsfile information
	if ( tag->hasOption( "lipsfile" ) ) {
		lipsfile_ = tag->getOption< std::string >( "lipsfile" );
	}

	// Read in membrane residue anchor
	if ( tag->hasOption( "anchor_rsd" ) ) {
		anchor_rsd_ = tag->getOption< core::Size >( "anchor_rsd" );
	}

	// Read in membrane residue position where applicable
	if ( tag->hasOption( "membrane_rsd" ) ) {
		membrane_rsd_ = tag->getOption< core::Size >( "membrane_rsd" );
	}

	// Read in membrane thickness
	if ( tag->hasOption( "thickness" ) ) {
		thickness_ = tag->getOption< core::Real >( "thickness" );
	}

	// Read in membrane steepness
	if ( tag->hasOption( "steepness" ) ) {
		steepness_ = tag->getOption< core::Real >( "steepness" );
	}

	// Use general method to read in center / normal
	read_center_normal_from_tag( center_, normal_, tag );

}

/// @brief Create a new copy of this mover
// XRW TEMP protocols::moves::MoverOP
// XRW TEMP AddMembraneMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new AddMembraneMover );
// XRW TEMP }

/// @brief Return the Name of this mover (as seen by Rscripts)
// XRW TEMP std::string
// XRW TEMP AddMembraneMoverCreator::keyname() const {
// XRW TEMP  return AddMembraneMover::mover_name();
// XRW TEMP }

/// @brief Mover name for Rosetta Scripts
// XRW TEMP std::string
// XRW TEMP AddMembraneMover::mover_name() {
// XRW TEMP  return "AddMembraneMover";
// XRW TEMP }

////////////////////////////////////////////////////
/// Getters - For subclasses of AddMembraneMover ///
////////////////////////////////////////////////////

/// @brief Return the current path to the spanfile held
/// by this mover
std::string AddMembraneMover::get_spanfile() const { return spanfile_; }

////////////////////////////////////////////////////
/// Setters - For subclasses of AddMembraneMover ///
////////////////////////////////////////////////////

/// @brief Set Spanfile path
/// @details Set the path to the spanfile
void AddMembraneMover::spanfile( std::string spanfile ) { spanfile_ = spanfile; }

/// @brief Set lipsfile path
/// @details Set the path to the lipsfile
void AddMembraneMover::lipsfile( std::string lipsfile ) { lipsfile_ = lipsfile; }

/// @brief Set option for including lipophilicity data
/// @details Incidate whether lipophilicity information should be read
/// and used in MembraneInfo
void AddMembraneMover::include_lips( bool include_lips ) { include_lips_ = include_lips; }

/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (AddMembraneMover)
// XRW TEMP std::string
// XRW TEMP AddMembraneMover::get_name() const {
// XRW TEMP  return "AddMembraneMover";
// XRW TEMP }


/// @brief Add Membrane Components to Pose
/// @details Add membrane components to pose which includes
/// spanning topology, lips info, embeddings, and a membrane
/// virtual residue describing the membrane position
void
AddMembraneMover::apply( Pose & pose ) {

	using namespace utility;
	using namespace core::pose;
	using namespace core::scoring;
	using namespace core::pack::task;
	using namespace core::conformation::membrane;

	TR << "=====================================================================" << std::endl;
	TR << "||           WELCOME TO THE WORLD OF MEMBRANE PROTEINS...          ||" << std::endl;
	TR << "=====================================================================" << std::endl;

	// Step 1: Initialize the Membrane residue
	core::Size membrane_pos = initialize_membrane_residue( pose, membrane_rsd_ );

	// Step 2: Initialize the spanning topology
	if ( topology_->nres_topo() == 0 ) {
		if ( spanfile_ == "single_TM_mode" ) {
			TR << "Single TM Mode: Defining a single TM helix from the length of the input pose" << std::endl;
			// Define teh span as teh whole helix, excluding the membrane residue
			SpanCOP single_TM_span = SpanCOP( new Span(  1, pose.size()-1 ) );
			topology_->add_span( *single_TM_span );
		} else if ( spanfile_ != "from_structure" ) {
			// now supports PDB numbering also
			std::map< std::string, core::Size > pdb2pose_map = core::pose::get_pdb2pose_numbering_as_stdmap( pose );
			// fill topology from spanfile with pdb2pose map
			topology_->fill_from_spanfile( spanfile_, pdb2pose_map, pose.size() );
		} else {
			// get pose info to create topology from structure
			std::pair< utility::vector1< core::Real >, utility::vector1< core::Size > > pose_info( get_chain_and_z( pose ) );
			utility::vector1< core::Real > z_coord = pose_info.first;
			utility::vector1< Size > chainID = pose_info.second;
			utility::vector1< char > secstruct( get_secstruct( pose ) );

			// set topology from structure
			topology_->fill_from_structure( z_coord, chainID, secstruct, thickness_ );
		}
	}

	// Step 3: Initialize the membrane Info Object
	// Options: Will initialize with or without a lipid acc object
	core::Size numjumps = pose.fold_tree().num_jump();
	MembraneInfoOP mem_info;
	if ( !include_lips_ ) {
		mem_info = MembraneInfoOP(
			new MembraneInfo( static_cast< core::Size >( membrane_pos ), numjumps, topology_ ) );
	} else {
		LipidAccInfoOP lips( new LipidAccInfo( lipsfile_ ) );
		mem_info = MembraneInfoOP( new MembraneInfo( static_cast< core::Size >( membrane_pos ), numjumps, lips, topology_ ) );
	}

	// Step 4: Add membrane info object to the pose conformation
	pose.conformation().set_membrane_info( mem_info );

	// Step 5: Accommodate user defined positions
	if ( user_defined_ ) {
		TR << "Setting initial membrane center and normal to position used by the user-provided membrane residue" << std::endl;
		center_ = pose.conformation().membrane_info()->membrane_center(pose.conformation());
		normal_ = pose.conformation().membrane_info()->membrane_normal(pose.conformation());

		// get anchor points from chains
		utility::vector1< core::Size > anchors( get_anchor_points_for_tmcom( pose ) );

		// replace root anchor residue with user-defined one
		anchors[1] = anchor_rsd_;

		// create membrane foldtree from anchors
		create_membrane_foldtree_from_anchors( pose, anchors );

	}

	// Step 5: Update the membrane position and orientation based on user settings
	SetMembranePositionMoverOP set_position = SetMembranePositionMoverOP( new SetMembranePositionMover( center_, normal_ ) );
	set_position->apply( pose );

	// Print Out Final FoldTree
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}

/// @brief Register options from JD2
void
AddMembraneMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
	option.add_relevant( OptionKeys::mp::setup::spans_from_structure );
	option.add_relevant( OptionKeys::mp::setup::lipsfile );
	option.add_relevant( OptionKeys::mp::setup::membrane_rsd );
	option.add_relevant( OptionKeys::mp::thickness );
	option.add_relevant( OptionKeys::mp::steepness );

}

/// @brief Initialize Mover options from the comandline
void
AddMembraneMover::init_from_cmd() {

	using namespace basic::options;

	// Read in User-Provided spanfile
	if ( spanfile_.size() == 0 && option[ OptionKeys::mp::setup::spanfiles ].user() ) {
		spanfile_ = option[ OptionKeys::mp::setup::spanfiles ]()[1];
	}

	if ( spanfile_.size() == 0 && option[ OptionKeys::mp::setup::spans_from_structure ].user() ) {
		TR << "WARNING: Spanfile not given, topology will be created from PDB!" << std::endl;
		TR << "WARNING: Make sure your PDB is transformed into membrane coordinates!!!" << std::endl;
		spanfile_ = "from_structure";
	}

	// Read in User-provided lipsfiles
	if ( option[ OptionKeys::mp::setup::lipsfile ].user() ) {

		// Set include lips to true and read in filename
		include_lips_ = true;
		lipsfile_ = option[ OptionKeys::mp::setup::lipsfile ]();
	}

	// Read in user-provided membrane residue position
	if ( option[ OptionKeys::mp::setup::membrane_rsd ].user() ) {
		membrane_rsd_ = option[ OptionKeys::mp::setup::membrane_rsd ]();
	}

	// Read in membrane thickness parameter
	if ( option[ OptionKeys::mp::thickness ].user() ) {
		TR << "Warning: About to set a new membrane thickness not currently optimized for the scoring function" << std::endl;
		thickness_ = option[ OptionKeys::mp::thickness ]();
	}

	// Read in transition steepness parameter
	if ( option[ OptionKeys::mp::steepness ].user() ) {
		TR << "Warning: About to set a new membrane steepness not currently optimized for the scoring function" << std::endl;
		steepness_ = option[ OptionKeys::mp::steepness ]();
	}

	// Use general method to read in center/normal from the command line
	read_center_normal_from_cmd( center_, normal_ );

}

/// @brief Initialize Membrane Residue given pose
core::Size
AddMembraneMover::initialize_membrane_residue( core::pose::Pose & pose, core::Size membrane_pos ) {

	// Start by assuming no membrane residue is provided. Check several possibilities
	// before setting up a new virtual residue

	// Case 1: If user points directly to a membrane residue, use that residue
	if ( membrane_pos != 0 &&
			membrane_pos <= pose.size() &&
			pose.conformation().residue( membrane_pos ).has_property( "MEMBRANE" ) ) {
		user_defined_ = true;

	} else {

		// Search for a membrane residue in the PDB
		utility::vector1< core::SSize > found_mem_rsds = check_pdb_for_mem( pose );

		// Case 2: Multiple membrane residues found
		if ( found_mem_rsds.size() > 1 ) {

			// Case 2a: Multiple membrane residues found, the user does not specify. Pick the first found RSD
			if ( membrane_rsd_ == 0 ) {
				TR << "Multiple membrane residues found. Setting MEM to the first found residue. Proceed with caution" << std::endl;
				membrane_pos = found_mem_rsds[1];
				user_defined_ = false;

				// Case 2b: Multiple membrane residues found AND the user specified found residue
				// matches a residue in the 'found' list
			} else {
				core::SSize current_rsd = static_cast< core::SSize >( membrane_rsd_ );
				for ( core::Size i = 1; i <= found_mem_rsds.size(); i++ ) {
					if ( current_rsd == found_mem_rsds[1] ) {
						TR << "Adding membrane residue from PDB at residue number " << membrane_rsd_ << std::endl;
						membrane_pos = membrane_rsd_;
						user_defined_ = true;
					}
				}
			}

			// Case 3: If one residue found in PDB and user didn't designate this residue, still accept found residue
		} else if ( (membrane_rsd_ == 0) && (found_mem_rsds[1] != -1) ) {
			TR << "No flag given: Adding membrane residue from PDB at residue number " << found_mem_rsds[1] << std::endl;
			membrane_pos = found_mem_rsds[1];
			user_defined_ = true;

			// Case 4: If membrane found and agrees with user specified value, accept
		} else if ( static_cast< core::SSize >( membrane_rsd_ ) == found_mem_rsds[1] ) {
			TR << "User specified residue matches found membrane residue. Accepting." << std::endl;
			membrane_pos = found_mem_rsds[1];
			user_defined_ = true;

			// Case 5: If no membrane residue found, add a new one to the pose
		} else if ( found_mem_rsds[1] == -1 ) {
			TR << "Adding a new membrane residue to the pose" << std::endl;
			membrane_pos = add_membrane_virtual( pose );

			// Case 6: Doesn't exist ;)
		} else {
			TR << "Congratulations - you have reached an edge case for adding the memrbane residue that we haven't thought of yet!" << std::endl;
			TR << "Contact the developers - exiting for now..." << std::endl;
			utility_exit();
		}
	}

	// DONE :D
	return membrane_pos;
}

/// @brief Helper Method - Add a membrane virtual residue
core::Size
AddMembraneMover::add_membrane_virtual( core::pose::Pose & pose ) {

	TR << "Adding a membrane residue representing the position of the membrane after residue " << pose.size() << std::endl;

	using namespace protocols::membrane;
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;

	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set(
		ChemicalManager::get_instance()->residue_type_set( pose.is_fullatom() ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
	);

	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOP rsd_type( residue_set->get_representative_type_name3("MEM") );
	ResidueType const & membrane( *rsd_type );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );

	// Append residue by jump, creating a new chain
	pose.append_residue_by_jump( *rsd, anchor_rsd_, "", "", true );

	// Reorder to be the membrane root
	FoldTreeOP ft( new FoldTree( pose.fold_tree() ) );
	ft->reorder( pose.size() );
	pose.fold_tree( *ft );
	if ( TR.visible() ) pose.fold_tree().show( TR );

	// Updating Chain Record in PDB Info
	char curr_chain = pose.pdb_info()->chain( pose.size()-1 );
	char new_chain = (char)((int) curr_chain + 1);
	pose.pdb_info()->number( pose.size(), (int) pose.size() );
	pose.pdb_info()->chain( pose.size(), new_chain );
	pose.pdb_info()->obsolete(false);

	return pose.size();

}

std::string AddMembraneMover::get_name() const {
	return mover_name();
}

std::string AddMembraneMover::mover_name() {
	return "AddMembraneMover";
}

void
AddMembraneMover::attributes_for_parse_center_normal_from_tag( utility::tag::AttributeList & attributes )
{
	// these attributes are actually more complex than just strings
	using namespace utility::tag;
	// AMW: these are not required
	attributes + XMLSchemaAttribute( "center", xs_string , "Position of center of membrane in format x,y,z" )
		+ XMLSchemaAttribute("normal", xs_string, "Membrane normal vector in format x,y,z");
}


void AddMembraneMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "include_lips", xsct_rosetta_bool, "Include lipid accessibility information from a lipsfile")
		+ XMLSchemaAttribute( "spanfile", xs_string, "Path to input spanfile")
		+ XMLSchemaAttribute( "lipsfile", xs_string, "Path to input lipsfile")
		+ XMLSchemaAttribute( "anchor_rsd", xsct_non_negative_integer, "Index of membrane residue anchor")
		+ XMLSchemaAttribute( "membrane_rsd", xsct_non_negative_integer, "Membrane residue position")
		+ XMLSchemaAttribute( "thickness", xsct_real, "Thickness of membrane. Score function is optimized to 15 Angstroms.")
		+ XMLSchemaAttribute( "steepness", xsct_real, "Steepness of membrane transition. Score function optimized to 10.");

	attributes_for_parse_center_normal_from_tag( attlist );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add membrane to a pose", attlist );
}

std::string AddMembraneMoverCreator::keyname() const {
	return AddMembraneMover::mover_name();
}

protocols::moves::MoverOP
AddMembraneMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMembraneMover );
}

void AddMembraneMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	AddMembraneMover::provide_xml_schema( xsd );
}


/// @brief Helper Method - Check for Membrane residue already in the PDB
utility::vector1< core::SSize >
AddMembraneMover::check_pdb_for_mem( core::pose::Pose & pose ) {

	// initialize vector for membrane residues found in PDB
	utility::vector1< core::Size > mem_rsd;

	// go through every residue in the pose to check for the membrane residue
	for ( core::Size i = 1; i <= pose.size(); ++i ) {

		// if residue is MEM, remember it
		if ( pose.residue( i ).name3() == "MEM" &&
				pose.residue( i ).has_property( "MEMBRANE" ) ) {

			TR << "Found membrane residue " << i << " in PDB." << std::endl;
			mem_rsd.push_back( static_cast< core::SSize >( i ) );
		}
	}

	// if no membrane residue
	if ( mem_rsd.size() == 0 ) {
		TR << "No membrane residue was found" << std::endl;
		mem_rsd.push_back( -1 );
	}

	return mem_rsd;

} // check pdb for mem

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_cc
