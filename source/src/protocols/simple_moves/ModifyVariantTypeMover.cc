// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/simple_moves/ModifyVariantTypeMover.cc
/// @brief Modify variant type to residues
/// @author Alex Ford <fordas@uw.edu>

// Unit Headers
#include <protocols/simple_moves/ModifyVariantTypeMover.hh>
#include <protocols/simple_moves/ModifyVariantTypeMoverCreator.hh>

// Package headers

// Project headers
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueProperties.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <protocols/rosetta_scripts/util.hh>

// tracer
#include <basic/Tracer.hh>

// Utility Headers

// C++ Headers
#include <iostream>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

#include <boost/algorithm/string.hpp>
#include <utility/excn/Exceptions.hh>
#include <boost/foreach.hpp>


// ObjexxFCL Headers

namespace protocols {
namespace simple_moves {

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.ModifyVariantTypeMover" );


// ModifyVariantTypeMover; based on the protocols::moves::Mover basis class
ModifyVariantTypeMover::ModifyVariantTypeMover() :
	protocols::moves::Mover("ModifyVariantType"),
	add_target_types_(),
	remove_target_types_(),
	residue_selector_()
{}

// @brief apply function here
void
ModifyVariantTypeMover::apply( core::pose::Pose & pose )
{
	core::select::residue_selector::ResidueSubset selection( pose.total_residue(), true );

	if ( residue_selector_ != nullptr ) {
		selection = residue_selector_->apply( pose );
		TR.Debug << "Initializing residue selection from ResidueSelector." << std::endl;
	} else {
		TR.Debug << "No ResidueSelector supplied.  Applying to ALL residues." << std::endl;
	}

	for ( core::Size resi = 1, resimax = pose.total_residue() ; resi <= resimax; ++resi ) {
		if ( selection[resi] ) {
			core::chemical::ResidueTypeSetCOP rsd_set(pose.residue(resi).residue_type_set());
			core::chemical::ResidueTypeCOP new_rsd_type = pose.residue(resi).type().get_self_ptr();

			BOOST_FOREACH ( std::string remove_type, remove_target_types_ ) {
				new_rsd_type = rsd_set->get_residue_type_with_variant_removed( *new_rsd_type,
					core::chemical::ResidueProperties::get_variant_from_string( remove_type ) ).get_self_ptr();
			}

			BOOST_FOREACH ( std::string add_type, add_target_types_ ) {
				new_rsd_type = rsd_set->get_residue_type_with_variant_added( *new_rsd_type,
					core::chemical::ResidueProperties::get_variant_from_string( add_type ) ).get_self_ptr();
			}

			core::pose::replace_pose_residue_copying_existing_coordinates( pose, resi, *new_rsd_type );
		}
	}
}

std::string
ModifyVariantTypeMover::get_name() const {
	return "ModifyVariantType";
}

moves::MoverOP
ModifyVariantTypeMover::clone() const
{
	return moves::MoverOP( new ModifyVariantTypeMover( *this ) );
}

moves::MoverOP
ModifyVariantTypeMover::fresh_instance() const
{
	return moves::MoverOP( new ModifyVariantTypeMover );
}

void ModifyVariantTypeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & /*pose*/
)
{
	add_target_types_.clear();
	remove_target_types_.clear();

	std::string add_type_value = tag->getOption< std::string >( "add_type", "");
	//boost::split(add_target_types_, add_type_value, boost::is_any_of(","));
	add_target_types_ = utility::string_split<std::string>(add_type_value,',',std::string());

	std::string remove_type_value = tag->getOption< std::string >( "remove_type", "");
	//boost::split(remove_target_types_, remove_type_value, boost::is_any_of(","));
	remove_target_types_ = utility::string_split<std::string>(remove_type_value,',',std::string());

	if ( add_target_types_.size() == 0 && remove_target_types_.size() == 0 ) {
		TR.Error << "Must specify add_type and/or remove_type type in ModifyVariantTypeMover." << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption("Must specify add_type and/or remove_type type in ModifyVariantTypeMover.");
	}

	if ( tag->hasOption( "residue_selector" ) ) {
		core::select::residue_selector::ResidueSelectorCOP selector( protocols::rosetta_scripts::parse_residue_selector( tag, data ) );
		if ( selector != nullptr ) {
			set_residue_selector( selector );
		}
	}
}

/// @brief Set the ResidueSelector used by this mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
void
ModifyVariantTypeMover::set_residue_selector(
	core::select::residue_selector::ResidueSelectorCOP selector_in
) {
	runtime_assert_string_msg( selector_in != nullptr , "Error in protocols::simple_moves::ModifyVariantTypeMover::set_residue_selector(): A null pointer was provided to this function." );
	residue_selector_ = selector_in;
}

protocols::moves::MoverOP
ModifyVariantTypeMoverCreator::create_mover() const { return protocols::moves::MoverOP( new ModifyVariantTypeMover ); }

std::string
ModifyVariantTypeMoverCreator::keyname() const { return "ModifyVariantType"; }

} // moves
} // protocols

