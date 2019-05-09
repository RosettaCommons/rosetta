// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/membrane/AddMembraneMover.cc
/// @brief   Initialize the RosettaMP framework by adding representations of the membrane
///          environemnt to the conformation data held by the pose
///
/// @details Given a pose, configure RosettaMP by adding the following information:
///             1. Create and append a membrane residue (MEM) to the pose
///             2. Create and store a SpanningTopology object
///             3. Setup the initial membrane coordinates (typically centered at the origin)
///             4. (Optional) Initialize per-atom lipid accessibility data
///             5. (Optional) Initialize dimensions of the aqueous pore
///             6. Initialize the ImplicitLipidMembraneInfo
///          Upon completion, the call pose.conformation().is_membrane() will return true
///
/// @note    If you add a new step, please document and ensure all data is properly
///          initialized by constructors, parse_my_tag, init_from_cmd, serialization
///          routines, and xsd routines. This class is a data mammoth
///
/// @author  Rebecca Faye Alford (rfalford12@gmail.com)
/// @author  JKLeman (julia.koehler.leman@gmail.com)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_cc
#define INCLUDED_protocols_membrane_AddMembraneMover_cc

// Unit Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/MPLipidAccessibility.hh>
#include <protocols/membrane/AqueousPoreFinder.hh>
#include <protocols/membrane/util.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/PDBPoseMap.hh>

#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneParams.hh>

#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/ImplicitLipidInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <protocols/membrane/SetMembranePositionMover.hh>

#include <protocols/moves/DsspMover.hh>
#include <core/pose/selection.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

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
#include <utility/sql_database/DatabaseSessionManager.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/database/sql_utils.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <boost/foreach.hpp>

#include <core/pose/util.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>
#include <map>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.membrane.AddMembraneMover" );

using cppdb::result;
using cppdb::statement;

namespace protocols {
namespace membrane {

//////////////////////////////////////////////////////////////////////////////////
/// Constructors ///
////////////////////

/// @brief Create a default RosettaMP membrane setup
AddMembraneMover::AddMembraneMover() :
	protocols::moves::Mover(),
	restore_lazaridis_IMM1_behavior_( false ),
	/* Lipid Composition Information */
	lipid_composition_( "DLPC" ),
	temperature_( 37.0 ),
	user_override_pore_( false ),
	/* Spanning Topology Information */
	spanfile_(),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	read_spans_from_xml_( false ),
	/* FoldTree Configuration */
	anchor_rsd_( 1 ),
	membrane_rsd_( 0 ),
	/* Membrane Position */
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	user_defined_membrane_pos_( false ),
	/* Bilayer Definition */
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_core_( 15 ),
	/* Read from Database Info */
	database_table_(""),
	db_session_()
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from a user specified spanfile
AddMembraneMover::AddMembraneMover(
	std::string spanfile,
	core::Size membrane_rsd,
	std::string lipid_composition,
	core::Real temperature
) : protocols::moves::Mover(),
	restore_lazaridis_IMM1_behavior_( false ),
	/* Lipid Composition Information */
	lipid_composition_( lipid_composition ),
	temperature_( temperature ),
	user_override_pore_( false ),
	/* Spanning Topology Information */
	spanfile_(std::move( spanfile )),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	read_spans_from_xml_( false ),
	/* FoldTree Configuration */
	anchor_rsd_( 1 ),
	/* Membrane Position */
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	user_defined_membrane_pos_( false ),
	/* Bilayer Definition */
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_core_( 15 ),
	/* Read from Database Info */
	database_table_(""),
	db_session_()
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from an existing membrane residue
AddMembraneMover::AddMembraneMover(
	core::Size membrane_rsd,
	std::string lipid_composition,
	core::Real temperature
) : protocols::moves::Mover(),
	restore_lazaridis_IMM1_behavior_( false ),
	/* Lipid Composition Information */
	lipid_composition_( lipid_composition ),
	temperature_( temperature ),
	user_override_pore_( false ),
	/* Spanning Topology Information */
	spanfile_(),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	read_spans_from_xml_( false ),
	/* FoldTree Configuration */
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
	/* Membrane Position */
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	user_defined_membrane_pos_( false ),
	/* Bilayer Definition */
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_core_( 15 ),
	/* Read from Database Info */
	database_table_(""),
	db_session_()
{
	register_options();
	init_from_cmd();
}


/// @brief Create a RosettaMP setup from a user specified SpanningTopology
AddMembraneMover::AddMembraneMover(
	core::conformation::membrane::SpanningTopologyOP topology,
	core::Size anchor_rsd,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	restore_lazaridis_IMM1_behavior_( false ),
	/* Lipid Composition Information */
	lipid_composition_( "DLPC" ),
	temperature_( 37.0 ),
	user_override_pore_( false ),
	/* Spanning Topology Information */
	spanfile_( "" ),
	topology_( std::move( topology ) ),
	read_spans_from_xml_( false ),
	/* FoldTree Configuration */
	anchor_rsd_( anchor_rsd ),
	membrane_rsd_( membrane_rsd ),
	/* Membrane Position */
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	user_defined_membrane_pos_( false ),
	/* Bilayer Definition */
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_core_( 15 ),
	/* Read from Database Info */
	database_table_(""),
	db_session_()
{
	register_options();
	init_from_cmd();
}

/// @brief Create a RosettaMP setup from an existing residue at a specific anchor point
AddMembraneMover::AddMembraneMover(
	core::Size anchor_rsd,
	bool /* dummy variable */,
	bool /* dummy variable */,
	core::Size membrane_rsd
) :
	protocols::moves::Mover(),
	restore_lazaridis_IMM1_behavior_( false ),
	/* Lipid Composition Information */
	lipid_composition_( "DLPC" ),
	temperature_( 37.0 ),
	user_override_pore_( false ),
	/* Spanning Topology Information */
	spanfile_( "" ),
	topology_( new core::conformation::membrane::SpanningTopology() ),
	read_spans_from_xml_( false ),
	/* FoldTree Configuration */
	anchor_rsd_( anchor_rsd ),
	membrane_rsd_( membrane_rsd ),
	/* Membrane Position */
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	user_defined_membrane_pos_( false ),
	/* Bilayer Definition */
	thickness_( 15 ),
	steepness_( 10 ),
	membrane_core_( 15 ),
	/* Read from Database Info */
	database_table_(""),
	db_session_()
{
	register_options();
	init_from_cmd();
}

/// @brief Create a deep copy of the data in this mover
AddMembraneMover::AddMembraneMover( AddMembraneMover const & ) = default;

/// @brief Destructor
AddMembraneMover::~AddMembraneMover() = default;

//////////////////////////////////////////////////////////////////////////////////
/// Mover-Specific Methods & Rosetta Scripts Support ///
//////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize the RosettaMP elements with this pose
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
	if ( topology_->nres_topo() == 0 && ! read_spans_from_xml_ ) {

		// In single TM mode, derive a single-span topology from the begin & end residues
		if ( spanfile_ == "single_TM_mode" ) {

			SpanCOP single_TM_span = SpanCOP( new Span(  1, pose.size()-1 ) );
			topology_->add_span( *single_TM_span );

			// In spanfile mode, define topology based on information in the spanfile
		} else if ( spanfile_ != "from_structure" ) {

			std::map< std::string, core::Size > pdb2pose_map = core::pose::get_pdb2pose_numbering_as_stdmap( pose );
			topology_->fill_from_spanfile( spanfile_, pdb2pose_map, pose.size() );

			// In "from structure" mode, approximate topology from the pose coordinates
		} else {

			std::pair< utility::vector1< core::Real >, utility::vector1< core::Size > > pose_info( get_chain_and_z( pose ) );
			utility::vector1< core::Real > z_coord = pose_info.first;
			utility::vector1< Size > chainID = pose_info.second;
			utility::vector1< char > secstruct( get_secstruct( pose ) );

			topology_->fill_from_structure( z_coord, chainID, secstruct, thickness_ );

		}
	}

	// Step 3: Initialize the membrane Info Object
	core::Size numjumps = pose.fold_tree().num_jump();
	if ( restore_lazaridis_IMM1_behavior_ ) {

		// Setup a default membrane info object
		MembraneInfoOP mem_info = MembraneInfoOP( new MembraneInfo(
			static_cast< core::Size >( membrane_pos ),
			numjumps,
			membrane_core_,
			thickness_,
			steepness_,
			topology_ )
		);
		pose.conformation().set_membrane_info( mem_info );

	} else {

		// Initialize a membrane info object with ImplicitLipidInfo
		MembraneInfoOP mem_info = MembraneInfoOP( new MembraneInfo(
			static_cast< core::Size >( membrane_pos ),
			numjumps,
			steepness_,
			topology_,
			lipid_composition_,
			temperature_ )
		);
		pose.conformation().set_membrane_info( mem_info );
	}

	// Step 4: If user defined membrane position, configure
	if ( user_defined_membrane_pos_ ) {

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

	SetMembranePositionMoverOP set_position = SetMembranePositionMoverOP( new SetMembranePositionMover( center_, normal_ ) );
	set_position->apply( pose );

	// Print Out Final FoldTree
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

	// IF applicable, setup the pore
	if ( !restore_lazaridis_IMM1_behavior_ && !user_override_pore_ ) {
		if ( pose.conformation().membrane_info()->spanning_topology()->nspans() >= 4 ) {

			// Calculate the per-atom lipid accessibility of the protein (geometry-based)
			MPLipidAccessibilityOP lipid_acc( new MPLipidAccessibility );
			core::Real implicit_lipid_thk( pose.conformation().membrane_info()->implicit_lipids()->water_thickness() );
			core::Real nslices( 3.0 );
			core::Real slice_width( (implicit_lipid_thk*2) / nslices );
			lipid_acc->set_slice_width( slice_width );
			lipid_acc->apply( pose );

			pose.conformation().membrane_info()->implicit_lipids()->is_helical( lipid_acc->is_alpha_helical() );

			// Approximate the pore geoemtry using a minimum bounding ellipse
			AqueousPoreFinderOP pore_finder( new AqueousPoreFinder );
			pore_finder->apply( pose );

			pose.conformation().membrane_info()->implicit_lipids()->has_pore( true );

		} else {
			pose.conformation().membrane_info()->implicit_lipids()->has_pore( false );
			pose.conformation().membrane_info()->implicit_lipids()->make_no_pore_parameters();
		}
	}



	// Show the results
	pose.conformation().membrane_info()->show( TR );

} // apply


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

//////////////////////////////////////////////////////////////////////////////////
/// Read & write access to data ///
//////////////////////////////////////////////////////////////////////////////////

/// @brief Path to file with transmembrane span data
std::string AddMembraneMover::get_spanfile() const { return spanfile_; }

/// @brief Set path to spanfile
void AddMembraneMover::spanfile( std::string spanfile ) { spanfile_ = spanfile; }

void AddMembraneMover::restore_lazaridis_IMM1_behavior( bool restore ) { restore_lazaridis_IMM1_behavior_ = restore; }


//////////////////////////////////////////////////////////////////////////////////
/// Functions that can be overloaded for sub-classes of this mover ///
//////////////////////////////////////////////////////////////////////////////////

/// @brief Helper Method - Add a membrane virtual residue
core::Size
AddMembraneMover::add_membrane_virtual( core::pose::Pose & pose ) {

	TR << "Adding a membrane residue representing the position of the membrane after residue " << pose.size() << std::endl;

	using namespace protocols::membrane;
	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;

	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set( pose.residue_type_set_for_pose() );

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
	if ( pose.pdb_info() != 0 ) {
		char curr_chain = pose.pdb_info()->chain( pose.size()-1 );
		char new_chain = (char)((int) curr_chain + 1);
		pose.pdb_info()->number( pose.size(), (int) pose.size() );
		pose.pdb_info()->chain( pose.size(), new_chain );
		pose.pdb_info()->obsolete(false);
	}

	return pose.size();

}

//////////////////////////////////////////////////////////////////////////////////
/// Helper methods for setup ///
//////////////////////////////////////////////////////////////////////////////////

/// @brief Initialize Membrane Residue given pose
core::Size
AddMembraneMover::initialize_membrane_residue( core::pose::Pose & pose, core::Size membrane_pos ) {

	// Start by assuming no membrane residue is provided. Check several possibilities
	// before setting up a new virtual residue

	// Case 1: If user points directly to a membrane residue, use that residue
	if ( membrane_pos != 0 &&
			membrane_pos <= pose.size() &&
			pose.conformation().residue( membrane_pos ).has_property( "MEMBRANE" ) ) {
		user_defined_membrane_pos_ = true;

	} else {

		// Search for a membrane residue in the PDB
		utility::vector1< core::SSize > found_mem_rsds = check_pdb_for_mem( pose );

		// Case 2: Multiple membrane residues found
		if ( found_mem_rsds.size() > 1 ) {

			// Case 2a: Multiple membrane residues found, the user does not specify. Pick the first found RSD
			if ( membrane_rsd_ == 0 ) {
				TR << "Multiple membrane residues found. Setting MEM to the first found residue. Proceed with caution" << std::endl;
				membrane_pos = found_mem_rsds[1];
				user_defined_membrane_pos_ = false;

				// Case 2b: Multiple membrane residues found AND the user specified found residue
				// matches a residue in the 'found' list
			} else {
				auto current_rsd = static_cast< core::SSize >( membrane_rsd_ );
				for ( core::Size i = 1; i <= found_mem_rsds.size(); i++ ) {
					if ( current_rsd == found_mem_rsds[1] ) {
						TR << "Adding membrane residue from PDB at residue number " << membrane_rsd_ << std::endl;
						membrane_pos = membrane_rsd_;
						user_defined_membrane_pos_ = true;
					}
				}
			}

			// Case 3: If one residue found in PDB and user didn't designate this residue, still accept found residue
		} else if ( (membrane_rsd_ == 0) && (found_mem_rsds[1] != -1) ) {
			TR << "No flag given: Adding membrane residue from PDB at residue number " << found_mem_rsds[1] << std::endl;
			membrane_pos = found_mem_rsds[1];
			user_defined_membrane_pos_ = true;

			// Case 4: If membrane found and agrees with user specified value, accept
		} else if ( static_cast< core::SSize >( membrane_rsd_ ) == found_mem_rsds[1] ) {
			TR << "User specified residue matches found membrane residue. Accepting." << std::endl;
			membrane_pos = found_mem_rsds[1];
			user_defined_membrane_pos_ = true;

			// Case 5: If no membrane residue found, add a new one to the pose
		} else if ( found_mem_rsds[1] == -1 ) {
			TR << "Adding a new membrane residue to the pose" << std::endl;
			membrane_pos = add_membrane_virtual( pose );

			// Case 6: Doesn't exist ;)
		} else {
			TR.Fatal << "Congratulations - you have reached an edge case for adding the membrane residue that we haven't thought of yet!" << std::endl;
			TR.Fatal << "Contact the developers - exiting for now..." << std::endl;
			utility_exit();
		}
	}

	return membrane_pos;
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


/// @brief Pase Rosetta Scripts Options for this Mover
void
AddMembraneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & pose
) {

	// read option to restore lazaridis parameters
	if ( tag->hasOption( "restore_lazaridis_IMM1_behavior" ) ) {
		restore_lazaridis_IMM1_behavior_ = tag->getOption< bool >( "restore_lazaridis_IMM1_behavior" );
	}

	// read lipid composition four letter code
	if ( tag->hasOption( "lipid_composition" ) ) {
		lipid_composition_ = tag->getOption< bool >( "lipid_composition" );
	}

	// read temperature for lipid composition geometry measurements
	if ( tag->hasOption( "temperature" ) ) {
		temperature_ = tag->getOption< bool >( "temperature" );
	}


	// Fail if the user tries to read from the database and an explicit spanfile tag
	if ( tag->hasOption( "database_name" ) && tag->hasOption( "spanfile" ) ) {
		utility_exit_with_message( "Error: Conflicting options! Cannot specify path to spanfile with spanfile=(path) and request to read from the database at the same time");
	}

	// Database name tag triggers DB IO, initialize variable for the table type and start the DB session
	if ( tag->hasOption( "database_name" ) ) {

		// Read in the name of the database table
		if ( tag->hasOption("database_table") ) {
			database_table_ = tag->getOption<std::string>("database_table");
		} else if ( tag->hasOption("table") ) {
			database_table_ = tag->getOption<std::string>("table");
		} else {
			throw ( CREATE_EXCEPTION( utility::excn::BadInput, "Cannot read input from database without specifying user option database_table or table" ) );
		}

		// Parse tag to initialize database session
		db_session_ = basic::database::parse_database_connection(tag);
	}

	// Read in spanfile information
	if ( tag->hasOption( "spanfile" ) ) {
		spanfile_ = tag->getOption< std::string >( "spanfile" );
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


	// Read in membrane core
	if ( tag->hasOption( "membrane_core" ) ) {
		membrane_core_ = tag->getOption< core::Real >( "membrane_core" );
	}

	if ( tag->hasOption( "span_starts" ) || tag->hasOption( "span_starts_num" ) ) {
		using namespace core::conformation::membrane;
		core::conformation::membrane::Orientation orientation;
		utility::vector1< core::Size > span_starts;
		span_starts.clear();
		utility::vector1< core::Size > span_ends;
		span_ends.clear();
		if ( tag->hasOption( "span_starts_num" ) ) {
			for ( auto const & i : utility::string_split( tag->getOption< std::string >("span_starts_num"), ',') ) {
				span_starts.push_back( std::atof( i.c_str() ) );
			}
			for ( auto const & i : utility::string_split( tag->getOption< std::string >("span_ends_num"), ',') ) {
				span_ends.push_back( std::atof( i.c_str() ) );
			}
		} else if ( tag->hasOption( "span_starts" ) ) {
			std::string span_starts_str = tag->getOption< std::string >("span_starts");
			std::string span_ends_str = tag->getOption< std::string >("span_ends");
			span_starts = core::pose::get_resnum_list_ordered( span_starts_str, pose );
			span_ends = core::pose::get_resnum_list_ordered( span_ends_str, pose );
		} else {
			runtime_assert_string_msg(true, "span_orientations specified, but neither span_starts nor span_starts_num were");
		}
		std::string span_oris_str( tag->getOption< std::string >("span_orientations") );
		utility::vector1<std::string> span_oris( utility::string_split(span_oris_str, ',') );
		runtime_assert( span_starts.size() == span_ends.size() && span_starts.size() == span_oris.size() );

		for ( core::Size i=1; i<=span_starts.size(); ++i ) {
			core::Size start = span_starts[ i ];
			core::Size end = span_ends[ i ];
			if ( span_oris[i] == "out2in" ) {
				orientation = out;
			} else if ( span_oris[i] == "in2out" ) {
				orientation = in;
			} else {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "provided orientation is wrong for span number acceptable options are in2out or out2in.");
			}
			core::conformation::membrane::SpanOP span( new core::conformation::membrane::Span( start, end, orientation ) );
			core::Size offset = 0;
			topology_->add_span( *span, offset );
		}
		TR << "finished reading span from span_starts/ends/orientations (xml) " << std::endl;
	}

	// Use general method to read in center / normal
	read_center_normal_from_tag( center_, normal_, tag );

	// if the user supplied a span in RosettaScripts xml, parse those and keep them in topology_. this will
	// be used to create the span topology
	utility::vector1<TagCOP> const sub_tags(tag->getTags());
	if ( !sub_tags.empty()  ) {
		if ( spanfile_ != "" ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "user provided span in a span file AND in the xml subtags. remove one of them");
		}
		BOOST_FOREACH ( TagCOP const sub_tag, sub_tags  ) {
			if ( sub_tag->getName() == "Span" ) {
				using namespace core::conformation::membrane;
				read_spans_from_xml_ = true;
				core::Size start = sub_tag->getOption< core::Size >( "start" );
				core::Size end  = sub_tag->getOption< core::Size >( "end" );
				std::string orientation_tag = sub_tag->getOption< std::string >( "orientation" );
				core::conformation::membrane::Orientation orientation;
				if ( orientation_tag == "out2in" ) {
					orientation = out;
				} else if ( orientation_tag == "in2out" ) {
					orientation = in;
				} else {
					throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "provided orientation is wrong for span number acceptable options are in2out or out2in.");
				}
				core::conformation::membrane::SpanOP span( new core::conformation::membrane::Span( start, end, orientation  ) );
				core::Size offset = 0;
				topology_->add_span( *span, offset );
			}
		}
		TR << "finished reading span from sub tags (xml) " << std::endl;
	}

}

/// @brief Register options from JD2
void
AddMembraneMover::register_options() {

	using namespace basic::options;

	option.add_relevant( OptionKeys::mp::setup::center );
	option.add_relevant( OptionKeys::mp::setup::normal );
	option.add_relevant( OptionKeys::mp::setup::spanfiles );
	option.add_relevant( OptionKeys::mp::setup::spans_from_structure );
	option.add_relevant( OptionKeys::mp::setup::membrane_rsd );
	option.add_relevant( OptionKeys::mp::lipids::composition );
	option.add_relevant( OptionKeys::mp::lipids::temperature );
	option.add_relevant( OptionKeys::mp::thickness );
	option.add_relevant( OptionKeys::mp::steepness );
	option.add_relevant( OptionKeys::mp::membrane_core );
	option.add_relevant( OptionKeys::mp::restore_lazaridis_imm_behavior );

}

/// @brief Initialize Mover options from the comandline
void
AddMembraneMover::init_from_cmd() {

	using namespace basic::options;

	// Read option to restore lazaridis IMM behavior
	if ( option[ OptionKeys::mp::restore_lazaridis_imm_behavior ].user() ) {
		restore_lazaridis_IMM1_behavior_ = option[ OptionKeys::mp::restore_lazaridis_imm_behavior ]();
	}

	// read in type of lipid composition
	if ( option[ OptionKeys::mp::lipids::composition ].user() ) {
		lipid_composition_ = option[ OptionKeys::mp::lipids::composition ]();
	}

	// read in temperature for lipid composition calculations
	if ( option[ OptionKeys::mp::lipids::temperature ].user() ) {
		temperature_ = option[ OptionKeys::mp::lipids::temperature ]();
	}

	// Override has pore option
	if ( option[ OptionKeys::mp::lipids::has_pore ].user() ) {
		user_override_pore_ = true;
	}

	// Read in User-Provided spanfile
	if ( spanfile_.size() == 0 && option[ OptionKeys::mp::setup::spanfiles ].user() ) {
		spanfile_ = option[ OptionKeys::mp::setup::spanfiles ]()[1];
	}

	if ( spanfile_.size() == 0 && option[ OptionKeys::mp::setup::spans_from_structure ].user() ) {
		TR.Warning << "Spanfile not given, topology will be created from PDB!" << std::endl;
		TR.Warning << "Make sure your PDB is transformed into membrane coordinates!!!" << std::endl;
		spanfile_ = "from_structure";
	}

	// Read in user-provided membrane residue position
	if ( option[ OptionKeys::mp::setup::membrane_rsd ].user() ) {
		membrane_rsd_ = option[ OptionKeys::mp::setup::membrane_rsd ]();
	}

	// Read in membrane thickness parameter
	if ( option[ OptionKeys::mp::thickness ].user() ) {
		TR.Warning << "About to set a new membrane thickness not currently optimized for the scoring function" << std::endl;
		thickness_ = option[ OptionKeys::mp::thickness ]();
	}

	// Read in transition steepness parameter
	if ( option[ OptionKeys::mp::steepness ].user() ) {
		TR.Warning << "About to set a new membrane steepness not currently optimized for the scoring function" << std::endl;
		steepness_ = option[ OptionKeys::mp::steepness ]();
	}

	// Read in membrane core parameter
	if ( option[ OptionKeys::mp::membrane_core ].user() ) {
		TR << "Warning: About to set a new membrane core not currently optimized for the scoring function" << std::endl;
		membrane_core_ = option[ OptionKeys::mp::membrane_core ]();
	}

	// Use general method to read in center/normal from the command line
	read_center_normal_from_cmd( center_, normal_ );

}

/// @brief Read a spanfile from a relational database provided as input to RosettaScripts
std::string
AddMembraneMover::read_spanfile_from_db() {

	using namespace basic::database;
	using namespace protocols::jd2;

	std::string tag(JobDistributor::get_instance()->current_job()->input_tag() + "spanfile");
	std::stringstream sql_stmt;
	sql_stmt
		<< "SELECT spanfile FROM " << database_table_
		<< " WHERE tag='" << tag << "';";
	std::string sql(sql_stmt.str());
	check_statement_sanity(sql);

	statement select_stmt(safely_prepare_statement(sql, db_session_));
	result res(safely_read_from_database(select_stmt));
	if ( !res.next() ) {
		std::stringstream error_message;
		error_message
			<< "Unable to locate spanfile for job distributor input tag '"
			<< tag << "' in the database." << std::endl;
		throw ( CREATE_EXCEPTION( utility::excn::BadInput, error_message.str() ) );
	}
	std::string spanfile;
	res >> spanfile;
	return spanfile;

}

/////////////////////
/// Mover Methods ///
/////////////////////

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
	attributes + XMLSchemaAttribute( "center", xs_string , "Position of center of membrane in format x,y,z" )
		+ XMLSchemaAttribute("normal", xs_string, "Membrane normal vector in format x,y,z");
}

void AddMembraneMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{

	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "restore_lazaridis_IMM1_behavior", xsct_rosetta_bool, "Restore Lazaridis IMM1 behavior" )
		+ XMLSchemaAttribute( "lipid_composition", xs_string, "Four letter code specifying the lipid composition" )
		+ XMLSchemaAttribute( "temperature", xsct_real, "Temperature at which the lipid parameters were measured/calculated" )
		+ XMLSchemaAttribute( "spanfile", xs_string, "Path to input spanfile")
		+ XMLSchemaAttribute( "anchor_rsd", xsct_non_negative_integer, "Index of membrane residue anchor")
		+ XMLSchemaAttribute( "membrane_rsd", xsct_non_negative_integer, "Membrane residue position")
		+ XMLSchemaAttribute( "thickness", xsct_real, "Thickness of membrane. Score function is optimized to 15 Angstroms.")
		+ XMLSchemaAttribute( "steepness", xsct_real, "Steepness of membrane transition. Score function optimized to 10.")
		+ XMLSchemaAttribute( "aqueous_pore", xsct_rosetta_bool, "Initialize pore dimensions and estimate per-atom lipid accessibility given protein geometry" )
		+ XMLSchemaAttribute( "membrane_core", xsct_real, "width of membrane ocre for Elazar calibrated LK potential" )
		+ XMLSchemaAttribute( "span_starts", xsct_residue_number_cslist, "comma separated list of span starting residues" )
		+ XMLSchemaAttribute( "span_ends", xsct_residue_number_cslist, "comma separated list of span ending residues" )
		+ XMLSchemaAttribute( "span_starts_num", xs_string, "comma separated list of span starting residues, in rosetta numbering" )
		+ XMLSchemaAttribute( "span_ends_num", xs_string, "comma separated list of span ending residues in rosetta_numbering" )
		+ XMLSchemaAttribute( "span_orientations", xs_string, "comma separated list of span orientations, only in2out or out2in allowed" );

	AttributeList span_subtag_attributes;
	span_subtag_attributes + XMLSchemaAttribute( "start", xsct_non_negative_integer, "residue where span starts" )
		+ XMLSchemaAttribute( "end", xsct_non_negative_integer, "resdiue where span ends" )
		+ XMLSchemaAttribute( "orientation", xs_string, "span orientation, whether in2out or out2in" );

	XMLSchemaSimpleSubelementList ssl;
	ssl.add_simple_subelement( "Span", span_subtag_attributes, "membrane spans" );

	attributes_for_parse_center_normal_from_tag( attlist );

	//protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add membrane to a pose", attlist);
	protocols::moves::xsd_type_definition_w_attributes_and_repeatable_subelements( xsd, mover_name(), "Add membrane to a pose", attlist, ssl );

}

//////////////////////////////////////////////////////////////////////////////////
/// Methods for Mover Creator class ///
//////////////////////////////////////////////////////////////////////////////////

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

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_cc
