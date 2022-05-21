// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/EnsureExclusivelySharedJumpMover.cc
/// @brief Set up foldtree such that there exists a jump that builds all selected residues and does not build any unselected residues
/// @author Jack Maguire

// Unit Headers
#include <protocols/docking/EnsureExclusivelySharedJumpMover.hh>
#include <protocols/docking/EnsureExclusivelySharedJumpMoverCreator.hh>

// Project headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/JumpDownstreamSelector.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>
#include <core/select/residue_selector/util.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>


using namespace protocols::moves;
using namespace core;
using namespace core::kinematics;
using namespace core::select::residue_selector;

static basic::Tracer TR( "protocols.docking.EnsureExclusivelySharedJumpMover" );

namespace protocols {
namespace docking {

namespace {//local utilities that will not be used by other translation units:

void
print( basic::Tracer & t, FoldTree const & f ){
	for ( auto const & e : f ) {
		t << e << std::endl;
	}
}

void
attempt_to_reroot(
	FoldTree & ft,
	utility::vector1< Size > const & unsele_positions
) {
	core::Size new_root = 0;
	for ( core::Size const resid : unsele_positions ) {
		if ( ft.possible_root( resid ) ) {
			new_root = resid;
			break;
		}
	}
	if ( new_root == 0 ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "EnsureExclusivelySharedJumpMover was unable to create a jump as specified because the root was selected. We attempted to re-root the foldtree from an unselected position but failed. Please re-root your foldtree to an unselected position." );
	}

	TR.Warning << "EnsureExclusivelySharedJumpMover is re-rooting the foldtree to accomodate the residue selector. The new root will be an unselected position" << std::endl;
	print( TR, ft );
	TR << "New Root: " << new_root << std::endl;
	bool const reorder_went_OK = ft.reorder( new_root );
	print( TR, ft );
	if ( ! reorder_went_OK ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "EnsureExclusivelySharedJumpMover was unable to create a jump as specified because the root was selected. We attempted to re-root the foldtree from an unselected position but failed. Please re-root your foldtree to an unselected position." );
	}
}

void
assert_jump_can_exist(
	FoldTree const & ft,
	utility::vector1< bool > const is_selected
){
	for ( core::kinematics::Edge const & e : ft ) {
		if ( e.is_jump() ) continue;

		core::Size const start = std::min( e.start(), e.stop() );
		core::Size const stop  = std::max( e.start(), e.stop() );

		for ( core::Size resid = start + 1; resid <= stop; ++resid ) {
			if ( is_selected[start] != is_selected[resid] ) {
				throw CREATE_EXCEPTION(utility::excn::Exception,  "EnsureExclusivelySharedJumpMover was unable to create a jump as specified. Residue " + std::to_string(start) + " and residue " + std::to_string( resid ) + " have opposite selection statuses but are built by the same jump." );
			}
		}
	}
}

}

EnsureExclusivelySharedJumpMover::EnsureExclusivelySharedJumpMover(
	core::select::residue_selector::ResidueSelectorCOP setting
):
	Mover(),
	selector_( setting )
{}

EnsureExclusivelySharedJumpMover::~EnsureExclusivelySharedJumpMover() = default;

moves::MoverOP
EnsureExclusivelySharedJumpMover::clone() const {
	return utility::pointer::make_shared< EnsureExclusivelySharedJumpMover >( *this );
}

moves::MoverOP
EnsureExclusivelySharedJumpMover::fresh_instance() const {
	return utility::pointer::make_shared< EnsureExclusivelySharedJumpMover >();
}


void
EnsureExclusivelySharedJumpMover::apply( core::pose::Pose & pose ) {
	using namespace core::kinematics;

	if ( selector_ == nullptr ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,  "No residue selector was provided to EnsureExclusivelySharedJumpMover" );
	}

	utility::vector1< bool > const is_selected = selector_->apply( pose );
	assert_jump_can_exist( pose.fold_tree(), is_selected );


	utility::vector1< Size > const sele_positions = selection_positions( is_selected );
	utility::vector1< Size > const unsele_positions = unselection_positions( is_selected );

	FoldTree ft = pose.fold_tree();
	bool const root_is_selected = is_selected[ ft.root() ];
	if ( root_is_selected ) {
		//throw CREATE_EXCEPTION(utility::excn::Exception,  "EnsureExclusivelySharedJumpMover was unable to create a jump as specified because the root was selected. Please re-root your foldtree to an unselected position." );
		attempt_to_reroot( ft, unsele_positions );
	}

	core::Size jump_id = 0;
	core::Size jump_stop_res = 0;
	for ( Edge const & e : ft ) {
		if ( ! e.is_jump() ) continue;

		//skip if they are both on the same side of the dividing line
		if ( is_selected[e.start()] == is_selected[e.stop()] ) {
			continue;
		} else if ( is_selected[e.start()] and !is_selected[e.stop()] ) {
			//replace selected starts with root (already unselected) for unselected stops
			Edge e2 = e;
			e2.start() = ft.root();
			TR << "Replacing " << e << " with " << e2 << std::endl;
			ft.replace_edge( e, e2 );
		} else if ( !is_selected[e.start()] and is_selected[e.stop()] ) {
			if ( jump_id == 0 ) {
				//if this is the first one we encounter, it will be the only survivor
				jump_id = e.label();
				jump_stop_res = e.stop();
			} else {
				//re-root this jump to the lone survivor
				runtime_assert( jump_stop_res != 0 );
				Edge e2 = e;
				e2.start() = jump_stop_res;
				TR << "Replacing " << e << " with " << e2 << std::endl;
				ft.replace_edge( e, e2 );
			}
		} else {
			//we THINK we have covered all possible options
			runtime_assert( false );
		}
	}//for edge in ft

	print( TR, ft );

	//UPDATE THE POSE
	pose.fold_tree( ft );

	print( TR, pose.fold_tree() );

	//Okay, now we THINK jump_id perfectly splits these residues. Let's check
	JumpDownstreamSelector const jd_selector( jump_id );
	utility::vector1< bool > const is_downstream = jd_selector.apply( pose );
	bool any_disagree = false;
	for ( core::Size resid = 1; resid <= pose.size(); ++resid ) {
		if ( is_selected[resid] != is_downstream[resid] ) {
			any_disagree = true;
			TR.Error << "Jump not as specified! resid = " << resid << ", is_selected[resid] = " << is_selected[resid] << ", and is_downstream[resid] = " << is_downstream[resid] << std::endl;
		}
	}
	runtime_assert( ! any_disagree );
}

// creator methods




std::string EnsureExclusivelySharedJumpMover::get_name() const {
	return mover_name();
}

std::string EnsureExclusivelySharedJumpMover::mover_name() {
	return "EnsureExclusivelySharedJumpMover";
}

void
EnsureExclusivelySharedJumpMover::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) {
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector" );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Set up foldtree such that there exists a jump that builds all selected residues and does not build any unselected residues", attlist );
}

void
EnsureExclusivelySharedJumpMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data
) {
	if ( ! tag->hasOption( "residue_selector" ) ) { // fetch selector from datamap
		throw CREATE_EXCEPTION(utility::excn::Exception,  "'residue_selector' was not provided for EnsureExclusivelySharedJumpMover" );
	}

	std::string selector_str = tag->getOption< std::string >( "residue_selector" );
	try {
		ResidueSelectorCOP selector = get_residue_selector( selector_str, data );
		set_selector( selector );
	} catch ( utility::excn::Exception & e ) {
		std::stringstream error_msg;
		error_msg << "Failed to find ResidueSelector named '" << selector_str << "' from the Datamap from EnsureExclusivelySharedJumpMover::parse_my_tag.\n";
		error_msg << e.msg();
		throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
	}
}

std::string
EnsureExclusivelySharedJumpMoverCreator::keyname() const {
	return EnsureExclusivelySharedJumpMover::mover_name();
}

protocols::moves::MoverOP
EnsureExclusivelySharedJumpMoverCreator::create_mover() const {
	return utility::pointer::make_shared< EnsureExclusivelySharedJumpMover >();
}

void
EnsureExclusivelySharedJumpMoverCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	EnsureExclusivelySharedJumpMover::provide_xml_schema( xsd );
}



}
}
