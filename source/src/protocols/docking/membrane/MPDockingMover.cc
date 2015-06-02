// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/MPDockingMover.cc
/// @brief      Dock two membrane proteins
/// @author     JKLeman (julia.koehler1982@gmail.com)
/// @note       Last Modified (6/24/14)

#ifndef INCLUDED_protocols_docking_membrane_MPDockingMover_cc
#define INCLUDED_protocols_docking_membrane_MPDockingMover_cc

// Unit Headers
#include <protocols/docking/membrane/MPDockingMover.hh>
#include <protocols/docking/membrane/MPDockingMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/membrane/MembraneInfo.hh>

#include <protocols/membrane/AddMembraneMover.hh>

#include <protocols/docking/DockingProtocol.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/mp.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.docking.membrane.MPDockingMover" );

namespace protocols {
namespace docking {
namespace membrane {

using namespace core;
using namespace core::import_pose;
using namespace core::pose;
using namespace core::scoring;
using namespace core::conformation::membrane;
using namespace protocols::membrane;
using namespace protocols::moves;
using namespace protocols::docking;
using namespace basic::options;
using namespace basic::options::OptionKeys;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
/// @details Docks two proteins with default normal=(0,0,1) and center=(0,0,0)
MPDockingMover::MPDockingMover( bool lowres, bool highres ) :
	protocols::moves::Mover(),
	add_membrane_mover_( new AddMembraneMover() ),
	docking_protocol_( new DockingProtocol() ),
	lowres_( lowres ),
	highres_( highres ),
	center_(0, 0, 0),
	normal_(0, 0, 1),
	jump_num_( 1 )
{}

/// @brief Default Constructor
/// @details Docks two proteins with default normal=(0,0,1) and center=(0,0,0)
MPDockingMover::MPDockingMover( Size jump_num, bool lowres, bool highres ) :
	protocols::moves::Mover(),
	add_membrane_mover_( new AddMembraneMover() ),
	docking_protocol_( new DockingProtocol() ),
	lowres_( lowres ),
	highres_( highres ),
	center_(0, 0, 0),
	normal_(0, 0, 1),
	jump_num_( jump_num )
{}


/// @brief Copy Constructor
MPDockingMover::MPDockingMover( MPDockingMover const & src ) :
	protocols::moves::Mover( src ), 
	add_membrane_mover_( src.add_membrane_mover_ ),
	docking_protocol_( src.docking_protocol_ ),
	lowres_ ( src.lowres_ ),
	highres_( src.highres_ ),
	lowres_scorefxn_( src.lowres_scorefxn_ ),
	highres_scorefxn_( src.highres_scorefxn_ ),
	center_( src.center_ ),
	normal_( src.normal_ ),
	native_( src.native_ )
{}

/// @brief Destructor
MPDockingMover::~MPDockingMover() {}

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPDockingMover::clone() const {
	return ( protocols::moves::MoverOP( new MPDockingMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPDockingMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPDockingMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPDockingMover::parse_my_tag(
    utility::tag::TagCOP tag,
    basic::datacache::DataMap &,
    protocols::filters::Filters_map const &,
    protocols::moves::Movers_map const &,
    core::pose::Pose const &
    ) {
    
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
    
    // Read in jump_num option
    if ( tag->hasOption( "jump_num" ) ) {
        jump_num_ = tag->getOption< core::Size >( "jump_num" );
    }
}
    
/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPDockingMoverCreator::create_mover() const {
    return protocols::moves::MoverOP( new MPDockingMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPDockingMoverCreator::keyname() const {
    return MPDockingMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPDockingMoverCreator::mover_name() {
    return "MPDockingMover";
}
    
/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPDockingMover)
std::string
MPDockingMover::get_name() const {
	return "MPDockingMover";
}

////////////////////////////////////////////////////////////////////////////////
// set inner cycles in DockingProtocol
void MPDockingMover::set_cycles_inner( Size cycles_inner ) {
	docking_protocol_->set_inner_cycles( cycles_inner );
}

// set outer cycles in DockingProtocol
void MPDockingMover::set_cycles_outer( Size cycles_outer ) {
	docking_protocol_->set_outer_cycles( cycles_outer );
}

////////////////////////////////////////////////////////////////////////////////
// sets default which can be overwritten by input flags
void MPDockingMover::set_defaults( const Pose & pose ){

	// set AddMembraneMover in protocol
	add_membrane_mover_ = AddMembraneMoverOP( new AddMembraneMover() );
	
	// create scorefunctions for lowres and highres
	// the ones I took were:
	// mpdocking_cen_14-7-23_no-penalties.wts
	// mpdocking_fa_14-7-23_no-penalties.wts
	// now I added the smooth term and took adjustments from MP fa score function
	lowres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpframework_docking_cen_2015.wts" );
	highres_scorefxn_ = ScoreFunctionFactory::create_score_function( "mpframework_docking_fa_2015.wts" );

	// set docking protocol
	if ( lowres_ && highres_ ){
		docking_protocol_ = DockingProtocolOP( new DockingProtocol( jump_num_, false, false, false, lowres_scorefxn_, highres_scorefxn_ ) );
	}
	else if ( lowres_ && ! highres_ ){
		docking_protocol_ = DockingProtocolOP( new DockingProtocol( jump_num_, true, false, false, lowres_scorefxn_, highres_scorefxn_ ) );
	}
	else if ( ! lowres_ && highres_ ){
		docking_protocol_ = DockingProtocolOP( new DockingProtocol( jump_num_, false, true, false, lowres_scorefxn_, highres_scorefxn_ ) );
	}
	else {
		utility_exit_with_message( "You want to run the docking protocol neither in lowres nor in highres??? Quitting..." );
	}
	
	// set native to pose, can be overwritten by flag -in:file:native
	native_ = PoseOP( new Pose( pose ) );

}// set defaults

////////////////////////////////////////////////////////////////////////////////
// register options
void MPDockingMover::register_options(){
	
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::docking::docking_local_refine );
	option.add_relevant( OptionKeys::mp::dock::weights_cen );
	option.add_relevant( OptionKeys::mp::dock::weights_fa );

}// register options

////////////////////////////////////////////////////////////////////////////////
// overwrite defaults from flags file
void MPDockingMover::init_from_cmd(){
	
	// if native flag given, set native from flag
	if ( option[OptionKeys::in::file::native].user() ){
		TR << "Setting native from flag -in::file::native" << std::endl;
		native_ = pose_from_pdb(option[OptionKeys::in::file::native].value_string() );
	}

	// if local_refine flag on, only do high-res
	if ( option[OptionKeys::docking::docking_local_refine].user() ){
		TR << "Running highres refinement only using flag -docking_local_refine" << std::endl;
		docking_protocol_ = DockingProtocolOP( new DockingProtocol( jump_num_, false, true, false, lowres_scorefxn_, highres_scorefxn_ ) );
	}

	// read low-res weights
	if ( option[OptionKeys::mp::dock::weights_cen].user() ){
		TR << "Weights for low-resolution step from flag -mp::dock::weights_cen" << std::endl;
		lowres_scorefxn_->reset();
		lowres_scorefxn_->initialize_from_file( option[OptionKeys::mp::dock::weights_cen].value_string() );
	}
	
	// read high-res weights
	if ( option[OptionKeys::mp::dock::weights_fa].user() ){
		TR << "Weights for high-resolution step from flag -mp::dock::weights_fa" << std::endl;
		highres_scorefxn_->reset();
		highres_scorefxn_->initialize_from_file( option[OptionKeys::mp::dock::weights_fa].value_string() );
	}

}// init_from_cmd

////////////////////////////////////////////////////////////////////////////////
// finalize setup
void MPDockingMover::finalize_setup(){
	
	// add membrane to native to have equal number of atoms for RMSD calculation
	add_membrane_mover_->apply( *native_ );
	
	// set native in docking protocol
	docking_protocol_->set_native_pose( native_ );
	
	// get movable jump
	TR.Debug << "movable jumps: " << to_string(docking_protocol_->movable_jumps()) << std::endl;

}// finalize setup

////////////////////////////////////////////////////////////////////////////////

/// @brief Add membrane components to the pose, then dock proteins along
///			the flexible jump
void MPDockingMover::apply( Pose & pose ) {
	
	// setup
	set_defaults( pose );
	
	// register options with JD2
	register_options();
	
	// overwrite defaults with stuff from cmdline
	init_from_cmd();
	
	// finalize setup
	finalize_setup();
	
	// assuming that protein 1 is fixed in the membrane!!!
	// add membrane VRT, call AddMembraneMover
	TR << "adding MEM" << std::endl;
	add_membrane_mover_->apply( pose );

	TR << "docking pose nres: " << pose.total_residue() << std::endl;
	TR << "native pose nres: " << docking_protocol_->get_native_pose()->total_residue() << std::endl;

	// creating foldtree from pose
	pose.fold_tree().show(std::cout);
	core::kinematics::FoldTree foldtree = pose.fold_tree();

	// reorder foldtree
	foldtree.reorder( pose.conformation().membrane_info()->membrane_rsd_num() );
	pose.fold_tree( foldtree );

	// show foldtree
	TR << "foldtree reordered" << std::endl;
	pose.fold_tree().show(std::cout);
		
	// run docking protocol
	TR << "calling docking protocol" << std::endl;
	docking_protocol_->apply( pose );

} // apply

} // membrane
} // docking
} // protocols

#endif // INCLUDED_protocols_docking_membrane_MPDockingMover_cc
