// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file       protocols/relax/membrane/MPRangeRelaxMover.hh
/// @brief      Relaxes a membrane protein by relaxing in ranges
/// @details Relaxes a membrane protein by iteratively relaxing ranges of the protein;
///    No ramping required. Much faster than FastRelax and good for
///    large to very large proteins (tested up to 5250 residues)
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_cc
#define INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_cc

// Unit Headers
#include <protocols/relax/membrane/MPRangeRelaxMover.hh>
//#include <protocols/membrane/MPRangeRelaxMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Project Headers
#include <protocols/membrane/AddMembraneMover.hh>
#include <protocols/relax/RangeRelaxMover.hh>
#include <core/conformation/membrane/MembraneInfo.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

// Package Headers
#include <core/kinematics/FoldTree.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <protocols/membrane/util.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <cstdlib>

static basic::Tracer TR( "protocols.relax.membrane.MPRangeRelaxMover" );

namespace protocols {
namespace relax {
namespace membrane {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;
using namespace protocols::moves;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default Constructor
MPRangeRelaxMover::MPRangeRelaxMover() : protocols::moves::Mover()
{
	set_defaults();
	register_options();
	init_from_cmd();
}

/// @brief Copy Constructor
/// @details Create a deep copy of this mover
MPRangeRelaxMover::MPRangeRelaxMover( MPRangeRelaxMover const & /*src*/ ) = default;

/// @brief Assignment Operator
MPRangeRelaxMover & MPRangeRelaxMover::operator = ( MPRangeRelaxMover const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Otherwise, create a new object
	return *( new MPRangeRelaxMover( *this ) );
}

/// @brief Destructor
MPRangeRelaxMover::~MPRangeRelaxMover() = default;

///////////////////////////////
/// Rosetta Scripts Methods ///
///////////////////////////////

/// @brief Create a Clone of this mover
protocols::moves::MoverOP
MPRangeRelaxMover::clone() const {
	return ( protocols::moves::MoverOP( new MPRangeRelaxMover( *this ) ) );
}

/// @brief Create a Fresh Instance of this Mover
protocols::moves::MoverOP
MPRangeRelaxMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MPRangeRelaxMover() );
}

/// @brief Pase Rosetta Scripts Options for this Mover
void
MPRangeRelaxMover::parse_my_tag(
	utility::tag::TagCOP /*tag*/,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {

	// TODO: implement this

}

/// @brief Create a new copy of this mover
//protocols::moves::MoverOP
//MPRangeRelaxMoverCreator::create_mover() const {
// return protocols::moves::MoverOP( new MPRangeRelaxMover() );
//}
//
///// @brief Return the Name of this mover (as seen by Rscripts)
//std::string
//MPRangeRelaxMoverCreator::keyname() const {
// return MPRangeRelaxMoverCreator::mover_name();
//}
//
///// @brief Mover name for Rosetta Scripts
//std::string
//MPRangeRelaxMoverCreator::mover_name() {
// return "MPRangeRelaxMover";
//}


/////////////////////
/// Mover Methods ///
/////////////////////

/// @brief Get the name of this Mover (MPRangeRelaxMover)
std::string
MPRangeRelaxMover::get_name() const {
	return "MPRangeRelaxMover";
}

/// @brief Do a RangeRelax of a membrane protein
void MPRangeRelaxMover::apply( Pose & pose ) {

	using namespace core::kinematics;
	using namespace protocols::membrane;
	using namespace protocols::relax;

	TR << "Relaxing a membrane protein with RangeRelax..." << std::endl;

	// add membrane to pose, if not already there
	if ( ! pose.conformation().is_membrane() ) {
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( pose );
	}

	// starting foldtree
	FoldTree orig_ft = pose.fold_tree();

	// reset foldtree to anchor on pose TM COM
	core::Size pose_tm_com( create_membrane_foldtree_anchor_pose_tmcom( pose ) );
	FoldTree ft = pose.fold_tree();
	ft.reorder( pose_tm_com );
	pose.fold_tree( ft );
	TR << "MPRangeRelax: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

	// set residues in TM region to helical secondary structure
	if ( set_tm_helical_ == true ) {
		for ( core::Size i = 1; i <= nres_protein( pose ); ++i ) {
			if ( pose.membrane_info()->spanning_topology()->in_span( i ) == true ) {

				TR << "Setting residue " << i << " to helical dihedral angles." << std::endl;

				pose.set_phi(   i, -62 );
				pose.set_psi(   i, -41 );
				pose.set_omega( i, 180 );
			}
		}
	}

	// call relax
	RangeRelaxMoverOP relax( new RangeRelaxMover( pose_tm_com ) );
	relax->add_membrane_again( false );
	relax->set_scorefunction( sfxn_ );
	// relax->optimize_membrane( optmem_ );
	if ( native_ != nullptr ) {
		relax->set_native( native_ );
	} else {
		relax->set_native( pose.clone() );
	}
	relax->apply( pose );

	// reset foldtree and show final one
	pose.fold_tree( orig_ft );
	TR << "Final foldtree: Is membrane fixed? " << protocols::membrane::is_membrane_fixed( pose ) << std::endl;
	pose.fold_tree().show( TR );

}// apply

/// @brief Optimize membrane
void MPRangeRelaxMover::optimize_membrane( bool yesno ) {
	optmem_ = yesno;
}

/////////////////////
/// Setup Methods ///
/////////////////////

/// @brief Set default values
void MPRangeRelaxMover::set_defaults() {

	// native
	native_ = nullptr;

	// create scorefunction
	sfxn_ = core::scoring::ScoreFunctionFactory::create_score_function( "mpframework_smooth_fa_2012.wts" );

	// center residue
	center_resnumber_ = 0;

	// set tm helical
	set_tm_helical_ = false;

	// optimize membrane?
	optmem_ = false;

}// set_defaults

/// @brief Register Options from Command Line
void MPRangeRelaxMover::register_options() {

	using namespace basic::options;
	option.add_relevant( OptionKeys::in::file::native );
	option.add_relevant( OptionKeys::relax::range::set_tm_helical );

} // register options


/// @brief Set default values
void MPRangeRelaxMover::init_from_cmd() {

	using namespace basic::options;
	using namespace protocols::membrane;

	// read native and attach membrane to it
	if ( option[ OptionKeys::in::file::native ].user() ) {
		native_ = core::import_pose::pose_from_file( option[ OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		AddMembraneMoverOP addmem( new AddMembraneMover() );
		addmem->apply( *native_ );
	}

	if ( option[ OptionKeys::relax::range::set_tm_helical ].user() ) {
		set_tm_helical_ = option[ OptionKeys::relax::range::set_tm_helical ]();
	}

}// init from cmd


} // membrane
} // relax
} // protocols

#endif // INCLUDED_protocols_relax_membrane_MPRangeRelaxMover_cc
