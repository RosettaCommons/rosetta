// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/HeteroHeteroBondFilter.cc
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)


//Unit Headers
#include <protocols/drug_design/HeteroHeteroBondFilter.hh>
#include <protocols/drug_design/HeteroHeteroBondFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/conformation/Residue.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace drug_design {

static basic::Tracer TR( "protocols.drug_design.HeteroHeteroBondFilter" );

protocols::filters::FilterOP
HeteroHeteroBondFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HeteroHeteroBondFilter ); }

std::string
HeteroHeteroBondFilterCreator::keyname() const { return HeteroHeteroBondFilter::class_name(); }

void HeteroHeteroBondFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	HeteroHeteroBondFilter::provide_xml_schema( xsd );
}

std::string HeteroHeteroBondFilter::class_name() {
	return "HeteroHeteroBond";
}

void HeteroHeteroBondFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_refpose_enabled_residue_number, "Residue to test" )
		+ XMLSchemaAttribute::attribute_w_default( "exclude_ates", xsct_rosetta_bool,
		"Exclude hetero-hetero bonds in things like phosphates, sulfates, etc. (Hetero-hetero bonds when there's also a double bond to oxygen in the group.)", "true" )
		+ XMLSchemaAttribute( "threshold", xsct_real, "Fail if the number of hetero-hetero bonds is greater than this." );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Test the number of heteroatom-heteroatom (non C/H) bonds in a residue.",
		attlist );
}

void
HeteroHeteroBondFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	residue_ = tag->getOption< std::string >( "residue" );
	threshold_ = tag->getOption<core::Size>( "threshold", threshold_ );
	exclude_ates_ = tag->getOption<bool>( "exclude_ates", exclude_ates_ );
	TR << "HeteroHeteroBond filter for residue " << residue_ << " with an upper threshold of  " << threshold_ << (exclude_ates_?" excluding":" including") << " oxidixed groups." << std::endl;
}

bool
HeteroHeteroBondFilter::apply( core::pose::Pose const & pose ) const {
	core::Size const value( compute( pose ) );
	if ( value > threshold_ ) {
		TR << "Failing HeteroHeteroBond filter on residue " << residue_ << " with value " << value << " which is above the threshold of " << threshold_ << std::endl;
		return false;
	} else {
		return true;
	}
}

void
HeteroHeteroBondFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Size const value( compute( pose ) );
	out << "HeteroHeteroBond value for residue " << residue_ << " is " << value << std::endl;
}

core::Real
HeteroHeteroBondFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Size
HeteroHeteroBondFilter::compute( core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 || resnum > pose.total_residue() ) {
		TR.Warning << "Attempted to access residue " << residue_ << " in pose with " <<  pose.total_residue() << " residues. Failing filter. " << std::endl;
		utility_exit_with_message("Cannot apply HeteroHeteroBond filter on non-existant residue!");
	}
	using namespace core::chemical;
	ResidueType const & restype( pose.residue(resnum).type() );

	core::Size count(0);

	for ( auto const & bondpair: restype.bonds() ) {
		element::Elements elem1( restype.element( bondpair.first ) );
		element::Elements elem2( restype.element( bondpair.second ) );
		if ( elem1 == element::H || elem1 == element::C || elem2 == element::H || elem2 == element::C ) {
			continue;
		}
		if ( exclude_ates_ ) {
			bool is_ate( false );
			// Look at the bonds attached to each of the two atoms in the bond, and reject if there exists a doublebond to oxygen.
			for ( core::Size nbr1: restype.bonded_neighbor( bondpair.first ) ) {
				if ( restype.bond_type( bondpair.first, nbr1 ) == DoubleBond && restype.element( nbr1 ) == element::O ) {
					is_ate = true;
					break;
				}
			}
			for ( core::Size nbr2: restype.bonded_neighbor( bondpair.second ) ) {
				if ( restype.bond_type( bondpair.second, nbr2 ) == DoubleBond && restype.element( nbr2 ) == element::O ) {
					is_ate = true;
					break;
				}
			}
			if ( is_ate ) {
				continue;
			}
		}
		++count;
	}

	TR << "HeteroHeteroBond for residue " << residue_ << ", of type " << restype.name() << " is " << count << ( exclude_ates_ ? " excluding":"including") << " oxidized groups." << std::endl;

	return count;
}

}
}
