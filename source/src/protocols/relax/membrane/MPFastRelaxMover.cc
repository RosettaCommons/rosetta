// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/relax/membrane/MPFastRelaxMover.cc
///
/// @brief      Membrane Fast Relax Protocol - Relax with minimization of mem position
/// @details	Apply the standard fast relax protocol. Enable minimization of the memrbane
///             jump and relax from the center of mass. Also use the smoothed
///             full atom membrane energy function.
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 12/2/14

// Unit Headers
#include <protocols/relax/membrane/MPFastRelaxMover.hh> 
#include <protocols/relax/membrane/MPFastRelaxMoverCreator.hh> 
#include <protocols/moves/Mover.hh> 

// Project Headers
#include <protocols/relax/RelaxProtocolBase.hh>
#include <protocols/relax/util.hh>

#include <protocols/membrane/MembranePositionFromTopologyMover.hh>
#include <protocols/membrane/geometry/util.hh>

// Package Headers
#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh> 

#include <core/kinematics/MoveMap.hh> 
#include <core/kinematics/FoldTree.hh>

#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.relax.membrane.MPFastRelaxMover" );

namespace protocols {
namespace relax {
namespace membrane {

using namespace core::pose;
using namespace protocols::membrane::geometry;
    
// Note - precondition: add membrane mover has already been
// called on this pose!
    
////////////////////////////
/// Constructors & Setup ///
////////////////////////////

/// @brief Default membrane fast relax constructor
/// @details Do normal fast relax protocol with the
/// membrane energy function and custom foldtree
MPFastRelaxMover::MPFastRelaxMover() : Mover() {
    relax_protocol_ = generate_relax_from_cmd();
}

/// @brief Destructor
MPFastRelaxMover::~MPFastRelaxMover() {}
    
/// @brief Show the current setup of this protocol
void
MPFastRelaxMover::show_protocol( Pose & pose ) {
    
    // Get the membrane position
    Vector center( pose.conformation().membrane_info()->membrane_center() );
    Vector normal( pose.conformation().membrane_info()->membrane_normal() );
    
    TR << "Membrane Relax protocol + MEM Optimization" << std::endl;
    TR << "Relax Type: " << relax_protocol_->get_name() << std::endl;
    TR << "Sfxn type: mpframework_fa_smooth_2012" << std::endl;
    TR << "Movemap: " << std::endl;
	relax_protocol_->get_movemap()->show();
    TR << "FoldTree: " << pose.fold_tree() << std::endl;
    TR << "Initial membrane position: " << std::endl;
    TR << "Membrane Posiiton: " << "center=(" << center.x() << "," << center.y() << "," << center.z() << "); normal=(" << normal.x() << "," << normal.y() << "," << center.z() << ")" << std::endl;
    
}

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPFastRelaxMover::clone() const {
    return ( protocols::moves::MoverOP( new MPFastRelaxMover() ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPFastRelaxMover::fresh_instance() const {
    return protocols::moves::MoverOP( new MPFastRelaxMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPFastRelaxMover::parse_my_tag(
  utility::tag::TagCOP,
  basic::datacache::DataMap &,
  protocols::filters::Filters_map const &,
  protocols::moves::Movers_map const &,
  core::pose::Pose const &
  )
{}
    
/// @brief Create a new copy of this mover
protocols::moves::MoverOP
MPFastRelaxMoverCreator::create_mover() const {
    return protocols::moves::MoverOP( new MPFastRelaxMover );
}

/// @brief Return the Name of this mover (as seen by Rscripts)
std::string
MPFastRelaxMoverCreator::keyname() const {
    return MPFastRelaxMoverCreator::mover_name();
}

/// @brief Mover name for Rosetta Scripts
std::string
MPFastRelaxMoverCreator::mover_name() {
    return "MPFastRelaxMover";
}

/////////////////////
/// Mover methods ///
/////////////////////

/// @brief Apply fast relax - do the actual protocol
void
MPFastRelaxMover::apply( Pose & pose ) {
    
    using namespace core::kinematics;
    using namespace core::scoring;
    using namespace protocols::membrane;
    
    // Check the pose is a membrane protein before proceeding
    if ( !pose.conformation().is_membrane() ) {
        utility_exit_with_message( "Cannot apply membrane fast relax to a non membrane pose." );
    }

    // Setup custom movemap based on the position of the membrane jump
    // Then set movemap in protocol
    core::SSize mem_jump( pose.conformation().membrane_info()->membrane_jump() );
    MoveMapOP movemap( new MoveMap() );
    movemap->set_chi( true );
    movemap->set_bb( true );
    movemap->set_jump( mem_jump, true );
    relax_protocol_->set_movemap( movemap );
    
    // Setup membrane fullatom smooth energy function
    ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012" );
    relax_protocol_->set_scorefxn( sfxn );

    // Setup custom foldtree where the membrane jump is attached at the
    // center of mass of the chain but other existing jumps are preserved
    setup_relax_foldtree( pose );
    
    // Set initial membrane protein position based on the position of
    // transmembrane spans
    MembranePositionFromTopologyMoverOP init_pos( new MembranePositionFromTopologyMover() );
    init_pos->apply( pose );
    
    // Show protocol settings before you go
    show_protocol( pose );
    
    // Apply membrane fast relax protocol
    relax_protocol_->apply( pose );
}

/// @brief Get name (MPFastRelaxMover)
/// @details Get the name of this mover
std::string
MPFastRelaxMover::get_name() const {
    return "MPFastRelaxMover";
}
    
/// @brief Create a custom foldtree anchored at the COM
/// @details Generate a foldtree where the membrane residue
/// is anchored at the center of mass of the chain. Also recreates
/// other jumps in the protein
void
MPFastRelaxMover::setup_relax_foldtree( Pose & pose ) {
    
    using namespace protocols::membrane;
    using namespace core::kinematics;
    
    // Abort if pose is not a membrane pose
    if ( !pose.conformation().is_membrane() ) {
        utility_exit_with_message( "Cannot build a membrane relax foldtree on a non membrane pose" );
    }
    
    TR << "Updating the pose foldtree to contain a jump from the pose center of mass to the membrane residue and preserve all downstream and upstream jumps for chains" << std::endl;
    
    // Setup a new simple foldtree
    FoldTree ft;
    ft.simple_tree( pose.total_residue() ); // Don't create a simple tree that includes the memrbane residue!!!
    
    // Count the number of jumps added
    core::Size njumps( 1 );
    
    // Add a jump between the protein COM and the memrbane residue
    core::Size membrane_rsd( pose.conformation().membrane_info()->membrane_rsd_num() );
    core::Size rsd_com( residue_center_of_mass( pose, 1, pose.total_residue()-1 ) ); // Get the center of mass and don't include the membrane residue!
    ft.new_jump( rsd_com, membrane_rsd, rsd_com );
    
    // Set the membrane jump in memInfo
    pose.conformation().membrane_info()->set_membrane_jump( njumps );
    
    // Rebuild chain jumps in order
    ++njumps;
    Size chain_begin(0), chain_end(0);
    // i == 1: not possible if I am a membrane pose
    // i == 2: i am a membrane pose, but don't rebuild my membrane jump
    // i == 3: 2 real polymer chains, 1 for the memrbane...now start counting
    // also - the membrane chain will be the last one
    if ( pose.conformation().num_chains() > 2 ) {
        for ( core::Size i = 1; i <= pose.conformation().num_chains()-2; i++ ) {
            chain_begin = pose.conformation().chain_begin( i+1 );
            chain_end = pose.conformation().chain_end( i );
            ft.new_jump( chain_end, chain_begin, chain_end );
            ++njumps;
        }
    }

    ft.show( std::cout );
    
    // Set the root to the residue COM (makes for a moveable membrane)
    // Will also verify that our foldtree is correct
    ft.reorder( rsd_com );
    
    // Set a new foldtree in the pose
    pose.fold_tree( ft );
}

} // membrane
} // relax
} // protocols
