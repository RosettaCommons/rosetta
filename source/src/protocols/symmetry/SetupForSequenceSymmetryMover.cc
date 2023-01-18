// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/symmetry/SetupForSequenceSymmetryMover.cc
/// @author Jack Maguire
/// @author Updated by Tim Neary, timdot10@gmail.com


// Unit headers
#include <protocols/symmetry/SetupForSequenceSymmetryMover.hh>
#include <protocols/symmetry/SetupForSequenceSymmetryMoverCreator.hh>

#include <protocols/moves/mover_schemas.hh>
#include <protocols/residue_selectors/StoreResidueSubsetMover.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pack/task/xml_util.hh>

#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <basic/options/option.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace symmetry {

static basic::Tracer TR( "protocols.simple_moves.symmetry.SetupForSequenceSymmetryMover" );

SetupForSequenceSymmetryMover::SetupForSequenceSymmetryMover()
{}

SetupForSequenceSymmetryMover::~SetupForSequenceSymmetryMover() = default;


void
SetupForSequenceSymmetryMover::add_residue_selector(
	core::Size region,
	core::select::residue_selector::ResidueSelectorCOP const & selector ) {

	if ( independent_regions_.size() <= region ) {
		while ( independent_regions_.size() <= region )
				{
			independent_regions_.emplace_back( utility::vector0< core::select::residue_selector::ResidueSelectorCOP >() );
		}
	}
	independent_regions_[region].emplace_back( selector );
}


void
SetupForSequenceSymmetryMover::apply( core::pose::Pose & pose ) {
	for ( core::Size ii = 0; ii < independent_regions_.size(); ++ii ) {
		validate_residue_selectors( pose, independent_regions_[ ii ] );
		for ( core::Size jj = 0; jj < independent_regions_[ii].size(); ++jj ) {
			std::string const magic_selector_name =
				setup_magic_name_prefix_ + std::to_string(ii) + "_" + std::to_string(jj);
			residue_selectors::StoreResidueSubsetMover add_subset( independent_regions_[ii][jj],
				magic_selector_name, true );
			add_subset.apply( pose );
#ifndef NDEBUG
			TR << "Added subset: " << magic_selector_name << std::endl;
#endif
		}
	}
}

void SetupForSequenceSymmetryMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) {
	using namespace basic::options;
	using namespace core::select::residue_selector;

	// Get the name of the corresponding KeepSequenceSymmetry
	std::string const kssto = tag->getOption< std::string >( "sequence_symmetry_behaviour" );
	if ( datamap.has( core::pack::task::TASK_OPERATIONS_TAG, kssto ) ) { //
		/* core::pack::task::TaskoperationCOP task_op =
		datamap.get< utility::tag::TagCOP >( "TaskOperation", kssto );
		if ( !task_op->getName() == "KeepSequenceSymmetry" ) {
		utility_exit_with_message( "Must specify a valid task operation of type KeepSequenceSymmetry." );
		} */
		// Must exist and have a name
		setup_magic_name_prefix_ = "SequenceSymmetricAnnealer_" + kssto + "_" ;
	} else {
		utility_exit_with_message( "Must specify a valid task operation of type KeepSequenceSymmetry." );
	}

	independent_regions_.clear();
	for ( utility::tag::TagCOP symm_def : tag->getTags() ) {
		if ( symm_def->hasOption( "residue_selectors" ) ) {
			utility::vector0< std::string > const res_selector_names =
				utility::string_split( symm_def->getOption< std::string >( "residue_selectors" ), ',' );

			TR << "Given residue selectors: ";
			for ( auto const & s : res_selector_names ) { TR << s << ", "; }
			TR << std::endl;

			if ( res_selector_names.size() == 0 ) {
				continue;
			}

			utility::vector0< core::select::residue_selector::ResidueSelectorCOP > res_seles;
			for ( auto rs : res_selector_names ) {
				res_seles.emplace_back( core::select::residue_selector::get_residue_selector( rs, datamap ) );
			}
			// Easiest place to validate the inputs is at apply time when we can access the vector< bool >
			independent_regions_.emplace_back( res_seles );
		}
	}
}

void
SetupForSequenceSymmetryMover::validate_residue_selectors(
	core::pose::Pose const & pose,
	utility::vector0< core::select::residue_selector::ResidueSelectorCOP > const & residue_selectors ) const {

	utility::vector0< core::Size > num_res;
	for ( auto res_sele : residue_selectors ) {
		core::Size count = 0;
		core::select::residue_selector::ResidueSubset const subset = res_sele->apply( pose );
		for ( auto const val : subset ) {
			if ( val ) ++count;
		}
		num_res.emplace_back( count );
	}

	for ( core::Size ii = 0; ii < num_res.size() - 1; ++ii ) {
		if ( num_res[ ii ] != num_res[ ii + 1 ] ) {
			utility_exit_with_message(
				"Provided invalid residue selectors. Each linked residue selector must define the same number of residues. "
				"Selectors of type " + residue_selectors[ ii ]->get_name() + " and " + residue_selectors[ ii + 1 ]->get_name() +
				" were found to specify " + std::to_string(num_res[ ii ]) + " and " + std::to_string(num_res[ ii + 1 ]) +
				" residue(s), respecitvely.");
		}
	}
}

void SetupForSequenceSymmetryMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "sequence_symmetry_behaviour", xs_string ,
		"Name of the KeepSequenceSymmetry task operation defininting the intended sequence symmetry behaviour." );

	AttributeList symmetry_attrs;
	symmetry_attrs + XMLSchemaAttribute( "residue_selectors", xs_string,
		"Comma separated list of selected residue selectors to define the sequence symmetry."
		" Residues in each listed selector will be linked."
		" Each residue selector listed must be equivalent." );

	XMLSchemaSimpleSubelementList attrs_subelements;
	attrs_subelements.add_simple_subelement( "SequenceSymmetry", symmetry_attrs,
		"Used to define each linked set of residues." );

	XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & protocols::moves::complex_type_name_for_mover )
		.element_name( mover_name() )
		.description( "Defines the regions of the protein to enforce sequence symmetry on." )
		.add_attributes( attlist )
		.add_optional_name_attribute()
		.set_subelements_repeatable( attrs_subelements, 1 )
		.write_complex_type_to_schema( xsd );
}

std::string SetupForSequenceSymmetryMoverCreator::keyname() const {
	return SetupForSequenceSymmetryMover::mover_name();
}

protocols::moves::MoverOP
SetupForSequenceSymmetryMoverCreator::create_mover() const {
	return utility::pointer::make_shared< SetupForSequenceSymmetryMover >();
}

void SetupForSequenceSymmetryMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SetupForSequenceSymmetryMover::provide_xml_schema( xsd );
}



} //symmetry
} // protocols
