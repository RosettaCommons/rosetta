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
/// @note       Last Modified (7/10/14)

#ifndef INCLUDED_protocols_membrane_AddMembraneMover_cc
#define INCLUDED_protocols_membrane_AddMembraneMover_cc

// Unit Headers
#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/AddMembraneMoverCreator.hh>

#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/geometry/util.hh>

#include <core/pose/PDBInfo.hh>

#include <core/conformation/Conformation.hh>

#include <core/conformation/membrane/MembraneParams.hh>

#include <core/conformation/membrane/MembraneInfo.hh> 
#include <core/conformation/membrane/SpanningTopology.hh> 
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/LipidAccInfo.hh>

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
#include <basic/options/keys/membrane_new.OptionKeys.gen.hh>

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
	fullatom_( true ), 
	include_lips_( false ),
	center_(0.0, 0.0, 0.0),
	normal_(0.0, 0.0, 10.0),
	membrane_rsd_( 0 )
{
	register_options();
	init_from_cmd();
}

/// @brief Custom Constructor - Create membrane pose from existing membrane resnum
/// @details Loads in membrane information from the commandline. Just needs a
/// membrane residue number (application is symmetry)
AddMembraneMover::AddMembraneMover( core::SSize membrane_rsd ) :
	protocols::moves::Mover(),
	fullatom_( true ),
	include_lips_( false ),
	center_(0.0, 0.0, 0.0),
	normal_(0.0, 0.0, 10.0),
	membrane_rsd_( membrane_rsd )
{
	register_options();
	init_from_cmd();
}
	
/// @brief Custom Constructor - for PyRosetta
/// @details Creates a membrane pose setting the membrane
/// center at emb_center and normal at emb_normal and will load
/// in spanning regions from list of spanfile provided
AddMembraneMover::AddMembraneMover(
	 Vector emb_center,
	 Vector emb_normal,
	 std::string spanfile,
	 core::SSize membrane_rsd
	 ) :
	 protocols::moves::Mover(),
	 fullatom_( true ), 
	 include_lips_( false ),
	 center_( emb_center ),
	 normal_( emb_normal ),
	 spanfile_( spanfile ),
	 membrane_rsd_( membrane_rsd )
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
	Vector emb_center,
	Vector emb_normal,
	std::string spanfile,
	std::string lips_acc,
	core::SSize membrane_rsd
	) :
	protocols::moves::Mover(),
	fullatom_( true ),
	include_lips_( true ),
	center_( emb_center ),
	normal_( emb_normal ),
	spanfile_( spanfile ),
	lipsfile_( lips_acc ),
	membrane_rsd_( membrane_rsd )
{
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
AddMembraneMover::AddMembraneMover( AddMembraneMover const & src ) :
	protocols::moves::Mover( src ), 
	fullatom_( src.fullatom_ ), 
	include_lips_( src.include_lips_ ),
	center_( src.center_ ), 
	normal_( src.normal_ ),
	spanfile_( src.spanfile_ ),
	lipsfile_( src.lipsfile_ ),
	membrane_rsd_( src.membrane_rsd_ )
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
	
	// Read in boolean fulltom
	if ( tag->hasOption( "fullatom" ) ) {
		fullatom_ = tag->getOption< bool >( "fullatom" );
	}
	
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
	
	// Read in membrane residue position where applicable
	if ( tag->hasOption( "membrane_rsd" ) ) {
		membrane_rsd_ = tag->getOption< core::SSize >( "membrane_rsd" );
	}
	
	// Read in membrane center & normal
	if ( tag->hasOption( "center" ) ) {
		std::string center = tag->getOption< std::string >( "center" );
		utility::vector1< std::string > str_cen = utility::string_split_multi_delim( center, ":,'`~+*&|;." );
		
		if ( str_cen.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			center_.x() = std::atof( str_cen[1].c_str() );
			center_.y() = std::atof( str_cen[2].c_str() );
			center_.z() = std::atof( str_cen[3].c_str() );
		}
	}

	if ( tag->hasOption( "normal" ) ) {
		std::string normal = tag->getOption< std::string >( "normal" );
		utility::vector1< std::string > str_norm = utility::string_split_multi_delim( normal, ":,'`~+*&|;." );
		
		if ( str_norm.size() != 3 ) {
			utility_exit_with_message( "Cannot read in xyz center vector from string - incorrect length!" );
		} else {
			normal_.x() = std::atof( str_norm[1].c_str() );
			normal_.y() = std::atof( str_norm[2].c_str() );
			normal_.z() = std::atof( str_norm[3].c_str() );
		}
	}
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
	
	// If there is a membrane residue in the PDB, take total resnum as the membrane_pos
	// Otherwise, setup a new membrane virtual
	core::SSize membrane_pos(0);
	if ( membrane_rsd_ != 0 &&
		 (core::SSize) pose.total_residue() > membrane_rsd_ &&
		 pose.residue( membrane_rsd_ ).has_property( "MEMBRANE" ) ) {
		TR << "Adding membrane residue from user specified position" << std::endl;
		membrane_pos = membrane_rsd_;
	} else if ( check_pdb_for_mem( pose ) ) {
		// TODO: Strenghten MEM search criteria in check pdb for mem
		TR << "Adding membrane from PDB to the pose" << std::endl;
		membrane_pos = pose.conformation().chain_end(1);
	} else {
		TR << "Adding a new membrane residue to the pose" << std::endl;
		membrane_pos = setup_membrane_virtual( pose );
	}
	
 	// Load spanning topology objects
	SpanningTopologyOP spans( new SpanningTopology( spanfile_, pose.total_residue()-1 ) );
	
	// Setup Membrane Info Object
	MembraneInfoOP mem_info;
	if ( !include_lips_ ) {
		mem_info = MembraneInfoOP( new MembraneInfo( pose.conformation(), membrane_pos, spans, 1 ) );
	} else {
		LipidAccInfoOP lips( new LipidAccInfo( lipsfile_ ) );
		mem_info = MembraneInfoOP( new MembraneInfo( pose.conformation(), membrane_pos, spans, lips, 1 ) );
	}
	
	// Add Membrane Info Object to conformation
	pose.conformation().set_membrane_info( mem_info );
	
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Register Options from Command Line
/// @details Register mover-relevant options with JD2 - includes
/// membrane_new, seutp options: center, normal, spanfile and
/// lipsfiles
void
AddMembraneMover::register_options() {
	
	using namespace basic::options;

	option.add_relevant( OptionKeys::membrane_new::setup::center );
	option.add_relevant( OptionKeys::membrane_new::setup::normal );
	option.add_relevant( OptionKeys::membrane_new::setup::spanfiles );
	option.add_relevant( OptionKeys::membrane_new::setup::lipsfile );
	option.add_relevant( OptionKeys::membrane_new::setup::membrane_rsd );
	
}

/// @brief Initialize Mover options from the comandline
/// @details Initialize mover settings from the commandline
/// mainly in the membrane_new, setup group: center, normal,
/// spanfile and lipsfiles paths
void
AddMembraneMover::init_from_cmd() {
	
	using namespace basic::options;
	
	// Set user-defined fullatom param
	if ( option[ OptionKeys::in::file::fullatom ].user() ) {
		fullatom_ = option[ OptionKeys::in::file::fullatom ]();
	}
	
	// Read in User-Provided spanfile
	if ( option[ OptionKeys::membrane_new::setup::spanfiles ].user() ) {
		spanfile_ = option[ OptionKeys::membrane_new::setup::spanfiles ]()[1];
	}
	
	// Read in User-provided lipsfiles
	if ( option[ OptionKeys::membrane_new::setup::lipsfile ].user() ) {
		
		// Set include lips to true and read in filename
		include_lips_ = true;
		lipsfile_ = option[ OptionKeys::membrane_new::setup::lipsfile ]();

	}
	
	// Read in user-provided membrane residue position
	// TODO: Add better error checking
	if ( option[ OptionKeys::membrane_new::setup::membrane_rsd ].user() ) {
		membrane_rsd_ = option[ OptionKeys::membrane_new::setup::membrane_rsd ]();
	}
	
	// Read in Center Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::membrane_new::setup::center ].user() ) {
		center_.x() = option[ OptionKeys::membrane_new::setup::center ]()[1];
		center_.y() = option[ OptionKeys::membrane_new::setup::center ]()[2];
		center_.z() = option[ OptionKeys::membrane_new::setup::center ]()[3];
	}
	
	// Read in Normal Parameter
	// TODO: Add better error checking
	if ( option[ OptionKeys::membrane_new::setup::normal ].user() ) {
		normal_.x() = option[ OptionKeys::membrane_new::setup::normal ]()[1];
		normal_.y() = option[ OptionKeys::membrane_new::setup::normal ]()[2];
		normal_.z() = option[ OptionKeys::membrane_new::setup::normal ]()[3];
	}
}

/// @brief Helper Method - Setup Membrane Virtual
/// @details Create a new virtual residue of type MEM from
/// the pose typeset (fullatom or centroid). Add this virtual
/// residue by appending to the last residue of the pose. Then set
/// this position as the root of the fold tree.
///
/// do not call this method before anchoring the pose with a virtual. Moves will
/// not occur in a reasonable way if so. @ralford 7/2/14
core::Size
AddMembraneMover::setup_membrane_virtual( Pose & pose ) {
	
	TR << "Adding a membrane residue representing the position of the membrane at " << pose.total_residue() << std::endl;
	
	using namespace protocols::membrane::geometry;
	using namespace core::conformation;
	using namespace core::chemical;
	
	// Grab the current residue typeset and create a new residue
	ResidueTypeSetCOP const & residue_set(
		ChemicalManager::get_instance()->residue_type_set( fullatom_ ? core::chemical::FA_STANDARD : core::chemical::CENTROID )
		);
		
	// Create a new Residue from rsd typeset of type MEM
	ResidueTypeCOPs const & rsd_type_list( residue_set->name3_map("MEM") );
	ResidueType const & membrane( *rsd_type_list[1] );
	ResidueOP rsd( ResidueFactory::create_residue( membrane ) );
	
	// Compute residue COM of the subunit
	core::SSize rsd_com = residue_center_of_mass( pose, 1, pose.total_residue() );
	
	// Append residue by jump, creating a new chain
	pose.append_residue_by_jump( *rsd, rsd_com, "", "", true );
	FoldTreeOP ft( new FoldTree( pose.fold_tree() ) );
	ft->reorder( rsd_com );
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
/// with property MEMBRANE, return it's residue number. In the control flow of the
/// apply method, if this returns non-zero, a new membrane residue will not be added.
bool
AddMembraneMover::check_pdb_for_mem( Pose & pose ) {

	// Uses the invarant that the membrane residue is always located in the n+1 chain
	// when initialized (not always true...?)
	core::Size rsdnum = pose.conformation().chain_end( 1 );
	if ( pose.residue( rsdnum ).name3() == "MEM" &&
		pose.residue( rsdnum ).has_property( "MEMBRANE" ) ) {
		return true;
	}
	return false;
	
}

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_AddMembraneMover_cc
