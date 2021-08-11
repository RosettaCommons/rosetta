// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/filters/HasDisulfideFilter.cc
/// @brief filters based on disulfide presence
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/filters/HasDisulfideFilter.hh>
#include <protocols/pose_sewing/filters/HasDisulfideFilterCreator.hh>

//XSD includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/pointer/memory.hh>
#include <protocols/filters/filter_schemas.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/annotated_sequence.hh>

static basic::Tracer TR( "protocols.pose_sewing.filters.HasDisulfideFilter" );

namespace protocols {
namespace pose_sewing {
namespace filters {

HasDisulfideFilter::HasDisulfideFilter():
	protocols::filters::Filter( "HasDisulfideFilter" )
{

}

HasDisulfideFilter::~HasDisulfideFilter()
{}

void
HasDisulfideFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "first_selector" ) ) {
		first_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "first_selector" ), datamap );
	}
	if ( tag->hasOption( "second_selector" ) ) {
		second_selector_ = core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "second_selector" ), datamap );
	}

}

protocols::filters::FilterOP
HasDisulfideFilter::clone() const
{
	return utility::pointer::make_shared< HasDisulfideFilter >( *this );
}


protocols::filters::FilterOP
HasDisulfideFilter::fresh_instance() const
{
	return utility::pointer::make_shared< HasDisulfideFilter >();
}

bool
HasDisulfideFilter::apply( core::pose::Pose const & pose) const
{
	core::select::residue_selector::ResidueSubset first_selection( pose.total_residue(), false );
	if ( first_selector_ != nullptr ) {
		first_selection = first_selector_->apply( pose );
	}
	core::select::residue_selector::ResidueSubset second_selection( pose.total_residue(), false );
	if ( second_selector_ != nullptr ) {
		second_selection = second_selector_->apply( pose );
	}
	for ( core::Size N_res = 1; N_res <= pose.size(); ++N_res ) {
		if ( first_selection[N_res] && pose.residue(N_res).type().is_disulfide_bonded() ) {
			for ( core::Size C_res = 1; C_res <= pose.size(); ++C_res ) {
				if ( second_selection[C_res] && pose.residue(C_res).type().is_disulfide_bonded() && pose.residue(N_res).is_bonded(pose.residue(C_res)) ) {
					return true;
				}
			}
		}
	}
	return false;
}
void
HasDisulfideFilter::set_first_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	first_selector_ = selector;
}
void
HasDisulfideFilter::set_second_selector( core::select::residue_selector::ResidueSelectorCOP selector ){
	second_selector_ = selector;
}

core::Real
HasDisulfideFilter::report_sm( core::pose::Pose const & ) const
{
	return -99999.9;
}

void
HasDisulfideFilter::report( std::ostream &, core::pose::Pose const & ) const
{

}

std::string HasDisulfideFilter::name() const {
	return class_name();
}

std::string HasDisulfideFilter::class_name() {
	return "HasDisulfideFilter";
}

void HasDisulfideFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "first_selector",  "The set of residues in which one end of the disulfide must be found." );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "second_selector",  "The set of residues of which one must be disulfide-bonded to a residue in the first set." );

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"Returns true if anything in first_selector is disulfide bonded to anything in second_selector",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
HasDisulfideFilterCreator::create_filter() const
{
	return utility::pointer::make_shared< HasDisulfideFilter >( );
}

std::string
HasDisulfideFilterCreator::keyname() const
{
	return HasDisulfideFilter::class_name();
}

void HasDisulfideFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HasDisulfideFilter::provide_xml_schema( xsd );
}

} //filters
} //pose_sewing
} //protocols
