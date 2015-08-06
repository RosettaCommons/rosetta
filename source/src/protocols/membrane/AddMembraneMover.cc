// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/AddMembraneMover.cc
///
/// @brief      Add Membrane Representation to the Pose
/// @details	Given a pose, setup membrane topology, lips info,
///				and a membrane virtual residue in the pose. All of this information
///				is coordinated via the MembraneInfo object maintained in
///				the Pose's conformation. After applying AddMembraneMover
///				to the pose, pose.conformation().is_membrane() should always
///				return true.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (3/29/14)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_cc
#define INCLUDED_protocols_membrane_AddMembraneMover_cc

// Unit Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/membrane/AddMembraneMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/util.hh>

#include <core/pose/PDBInfo.hh>

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

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.membrane.AddMembraneMover" );

namespace protocols {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;
		
/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Create a membrane pose setting the membrane center
/// at center=(0, 0, 0), normal=(0, 0, 1) and loads in spans
/// and lips from the command line interface.
AddMembraneMover::AddMembraneMover() :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(),
	topology_( new SpanningTopology() ),
	lipsfile_(),
    anchor_rsd_( 1 ),
	membrane_rsd_( 0 ),
    center_( 0, 0, 0 ),
    normal_( 0, 0, 1 ),
    thickness_( 15 ),
    steepness_( 10 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor - Create membrane pose from existing membrane resnum
/// @details Loads in membrane information from the commandline. Just needs a
/// membrane residue number (application is symmetry)
AddMembraneMover::AddMembraneMover( core::Size membrane_rsd ) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_(),
	topology_( new SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( 1 ),
	membrane_rsd_( membrane_rsd ),
    center_( 0, 0, 0 ),
    normal_( 0, 0, 1 ),
    thickness_( 15 ),
    steepness_( 10 )
{
	register_options();
	init_from_cmd();
}
	
/// @brief Custom Constructor - for PyRosetta
/// @details Creates a membrane pose setting the membrane
/// center at emb_center and normal at emb_normal and will load
/// in spanning regions from list of spanfile provided
AddMembraneMover::AddMembraneMover(
    std::string spanfile,
    core::Size membrane_rsd
    ) :
    protocols::moves::Mover(),
    include_lips_( false ),
    spanfile_( spanfile ),
    topology_( new SpanningTopology() ),
    lipsfile_(),
    anchor_rsd_( 1 ),
    membrane_rsd_( membrane_rsd ),
    center_( 0, 0, 0 ),
    normal_( 0, 0, 1 ),
    thickness_( 15 ),
    steepness_( 10 )
{
	register_options();
	init_from_cmd();
}
    
/// @brief Custom Constructor - mainly for PyRosetta
/// @details Creates a membrane pose setting the membrane
/// center at emb_center and normal at emb_normal and will load
/// in spanning regions from list of spanfiles provided
AddMembraneMover::AddMembraneMover(
    SpanningTopologyOP topology,
    core::Size anchor_rsd,
    core::Size membrane_rsd
    ) :
    protocols::moves::Mover(),
    include_lips_( false ),
    spanfile_( "" ),
    topology_( topology ),
    lipsfile_(),
    anchor_rsd_( anchor_rsd ),
    membrane_rsd_( membrane_rsd ),
    center_( 0, 0, 0 ),
    normal_( 0, 0, 1 ),
    thickness_( 15 ),
    steepness_( 10 )
{
    register_options();
    init_from_cmd();
}

	/// @brief Custom Constructor using both anchor rsd and membrane rsd
	/// @details Creates a membrane pose setting the anchor residue and
	///			membrane residue; topology is read from cmd line
	AddMembraneMover::AddMembraneMover(
		   core::Size anchor_rsd,
		   core::Size membrane_rsd
		   ) :
	protocols::moves::Mover(),
	include_lips_( false ),
	spanfile_( "" ),
	topology_( new SpanningTopology() ),
	lipsfile_(),
	anchor_rsd_( anchor_rsd ),
	membrane_rsd_( membrane_rsd ),
	center_( 0, 0, 0 ),
	normal_( 0, 0, 1 ),
	thickness_( 15 ),
	steepness_( 10 )
	{
		register_options();
		init_from_cmd();
	}

/// @brief Custorm Constructur with lips info - for PyRosetta
/// @details Creates a membrane pose setting the membrane
/// center at emb_center and normal at emb_normal and will load
/// in spanning regions from list of spanfile provided. Will also
/// load in lips info from lips_acc info provided
AddMembraneMover::AddMembraneMover(
	std::string spanfile,
	std::string lips_acc,
	core::Size membrane_rsd
	) :
    protocols::moves::Mover(),
    include_lips_( true ),
    spanfile_( spanfile ),
    topology_( new SpanningTopology() ),
    lipsfile_( lips_acc ),
    anchor_rsd_( 1 ),
    membrane_rsd_( membrane_rsd ),
    center_( 0, 0, 0 ),
    normal_( 0, 0, 1 ),
    thickness_( 15 ),
    steepness_( 10 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custorm Constructur - Center/Normal init
/// @details Creates a membrane protein, setting the initial
/// membrane center and normal to the specified values, loads
/// a spanfile and lipsfile from given path, and indicates a membrane
/// residue if provided.
AddMembraneMover::AddMembraneMover(
    Vector init_center,
    Vector init_normal,
    std::string spanfile,
    core::Size membrane_rsd
    ) :
    protocols::moves::Mover(),
    include_lips_( false ),
    spanfile_( spanfile ),
    topology_( new SpanningTopology() ),
    lipsfile_( "" ),
    anchor_rsd_( 1 ),
    membrane_rsd_( membrane_rsd ),
    center_( init_center ),
    normal_( init_normal ),
    thickness_( 15 ),
    steepness_( 10 )
{}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
AddMembraneMover::AddMembraneMover( AddMembraneMover const & src ) :
    protocols::moves::Mover( src ),
    include_lips_( src.include_lips_ ),
    spanfile_( src.spanfile_ ),
    topology_( src.topology_ ),
    lipsfile_( src.lipsfile_ ),
    anchor_rsd_( src.anchor_rsd_ ),
    membrane_rsd_( src.membrane_rsd_ ),
    center_( src.center_ ),
    normal_( src.normal_ ),
    thickness_( src.thickness_ ),
    steepness_( src.steepness_ )
{}

/// @brief Destructor
AddMembraneMover::~AddMembraneMover() {}

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
protocols::moves::MoverOP
AddMembraneMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new AddMembraneMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
AddMembraneMoverCreator::keyname() const {
	return AddMembraneMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
AddMembraneMoverCreator::mover_name() {
	return "AddMembraneMover";
}

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
std::string
AddMembraneMover::get_name() const {
	return "AddMembraneMover";
}


/// @brief Add Membrane Components to Pose
/// @details Add membrane components to pose which includes
///	spanning topology, lips info, embeddings, and a membrane
/// virtual residue describing the membrane position
void
AddMembraneMover::apply( Pose & pose ) {
	
	using namespace core::scoring;
	using namespace core::pack::task; 
	using namespace core::conformation::membrane;

	TR << "=====================================================================" << std::endl;
	TR << "||           WELCOME TO THE WORLD OF MEMBRANE PROTEINS...          ||" << std::endl;
	TR << "=====================================================================" << std::endl;
	
	// If there is a membrane residue in the PDB, take total resnum as the membrane_pos
	// Otherwise, setup a new membrane virtual
	core::SSize membrane_pos(0);
    bool user_defined( false );
	
	// search for MEM in PDB
	utility::vector1< core::SSize > mem_rsd_in_pdb = check_pdb_for_mem( pose );
    
    // If multiple membrane residues are found - exit. Not allowed in the framework
    if ( mem_rsd_in_pdb.size() > 1 ) {
        utility_exit_with_message( "Multiple membrane residues found in the pose, but only one allowed. Exiting..." );
    }
    
    // If exactly 1 residue was found and no flag provided, designate that residue
    if ( membrane_rsd_ == 0 && mem_rsd_in_pdb[1] != -1 ) {
        TR << "No flag given: Adding membrane residue from PDB at residue number " << mem_rsd_in_pdb[1] << std::endl;
        membrane_pos = mem_rsd_in_pdb[1];
        user_defined = true;
    }
    
    // if found and agrees with user-specified one, use that
    if ( static_cast< SSize >( membrane_rsd_ ) == mem_rsd_in_pdb[1] ) {
        TR << "Adding membrane residue from PDB at residue number " << membrane_rsd_ << std::endl;
        membrane_pos = membrane_rsd_;
        user_defined = true;
    }
    
    // If no membrane residue was found, add one to the pose
    if ( mem_rsd_in_pdb[1] == -1 ) {
        TR << "Adding a new membrane residue to the pose" << std::endl;
        membrane_pos = setup_membrane_virtual( pose );
    }
    
    /**
    this isn't accomplishing what it should be... if we want that strict of checking, check the types directly
     // If the found residue doesn't agree with user specified one, exit
    if ( membrane_pos == 0 && mem_rsd_in_pdb[1] != -1 ) {
        utility_exit_with_message("User provided membrane residue doesn't agree with found one!");
    }
     **/

 	// Load spanning topology objects
	if ( topology_->nres_topo() == 0 ){
		if ( spanfile_ != "from_structure" ) {
			topology_->fill_from_spanfile( spanfile_, pose.total_residue() );
		}
		else {
			// get pose info to create topology from structure
			std::pair< utility::vector1< Real >, utility::vector1< Size > > pose_info( get_chain_and_z( pose ) );
			utility::vector1< Real > z_coord = pose_info.first;
			utility::vector1< Size > chainID = pose_info.second;
			utility::vector1< char > secstruct( get_secstruct( pose ) );
			
			// set topology from structure
			topology_->fill_from_structure( z_coord, chainID, secstruct, thickness_ );
		}
	}

	// get number of jumps in foldtree
	Size numjumps = pose.fold_tree().num_jump();
	
	// Setup Membrane Info Object
	MembraneInfoOP mem_info;
	if ( !include_lips_ ) {
		mem_info = MembraneInfoOP( new MembraneInfo( pose.conformation(), static_cast< Size >( membrane_pos ), topology_, numjumps ) );
	} else {
		LipidAccInfoOP lips( new LipidAccInfo( lipsfile_ ) );
		mem_info = MembraneInfoOP( new MembraneInfo( pose.conformation(), static_cast< Size >( membrane_pos ), topology_, lips, numjumps ) );
	}
	
	// Add Membrane Info Object to conformation
	pose.conformation().set_membrane_info( mem_info );
    
    // Set initial membrane position
    if ( user_defined ) {
        TR << "Setting initial membrane center and normal to position used by the user-provided membrane residue" << std::endl;
        center_ = pose.conformation().membrane_info()->membrane_center();
        normal_ = pose.conformation().membrane_info()->membrane_normal();
		
		// slide the membrane jump to match the desired anchor point
		core::kinematics::FoldTree ft = pose.fold_tree();
		Size memjump = pose.conformation().membrane_info()->membrane_jump();
		ft.slide_jump( memjump, anchor_rsd_, mem_rsd_in_pdb[1] );
		pose.fold_tree( ft );
    }

    // Use set membrane positon mover to set the center/normal pair
    SetMembranePositionMoverOP set_position = SetMembranePositionMoverOP( new SetMembranePositionMover( center_, normal_ ) );
    set_position->apply( pose );

    // Set membrane thickness & steepness (Added 3/9/15)
    pose.conformation().membrane_info()->set_membrane_thickness( thickness_ );
    pose.conformation().membrane_info()->set_membrane_steepness( steepness_ );
	
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// mp, seutp options: center, normal, spanfile and
/// lipsfiles
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
/// @details Initialize mover settings from the commandline
/// mainly in the mp, setup group: center, normal,
/// spanfile and lipsfiles paths
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

/// @brief Helper Method - Setup Membrane Virtual
/// @details Create a new virtual residue of type MEM from
/// the pose typeset (fullatom or centroid). Add this virtual
/// residue by appending to the last residue of the pose. Then set
/// this position as the root of the fold tree.
core::Size
AddMembraneMover::setup_membrane_virtual( Pose & pose ) {
	
	TR << "Adding a membrane residue representing the position of the membrane after residue " << pose.total_residue() << std::endl;
	
	using namespace protocols::membrane;
	using namespace core::conformation;
	using namespace core::chemical;
	
	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set(
		ChemicalManager::get_instance()->residue_type_set( pose.is_fullatom() ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
		);
		
	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
	ResidueType const & membrane( *rsd_type_list[1] );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );
	
	// Append residue by jump, creating a new chain
	pose.append_residue_by_jump( *rsd, anchor_rsd_, "", "", true );
	FoldTreeOP ft( new FoldTree( pose.fold_tree() ) );

	// Reorder to be the membrane root
	ft->reorder( pose.total_residue() );
	pose.fold_tree( *ft );	
	pose.fold_tree().show( std::cout );

	// Updating Chain Record in PDB Info
	char curr_chain = pose.pdb_info()->chain( pose.total_residue()-1 );
	char new_chain = (char)((int) curr_chain + 1);
	pose.pdb_info()->number( pose.total_residue(), (int) pose.total_residue() );
	pose.pdb_info()->chain( pose.total_residue(), new_chain ); 
	pose.pdb_info()->obsolete(false);
	
	return pose.total_residue();
}

/// @brief Helper Method - Check for Membrane residue already in the PDB
/// @details If there is an MEM residue in the PDB at the end of the pose
/// with property MEMBRANE, return a vector of all of those residues.
utility::vector1< core::SSize >
AddMembraneMover::check_pdb_for_mem( Pose & pose ) {

	// initialize vector for membrane residues found in PDB
	utility::vector1< core::Size > mem_rsd;

	// go through every residue in the pose to check for the membrane residue
	for ( core::Size i = 1; i <= pose.total_residue(); ++i ){
		
		// if residue is MEM, remember it
		if ( pose.residue( i ).name3() == "MEM" &&
			pose.residue( i ).has_property( "MEMBRANE" ) ) {
			
			TR << "Found membrane residue " << i << " in PDB." << std::endl;
			mem_rsd.push_back( static_cast< SSize >( i ) );
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
