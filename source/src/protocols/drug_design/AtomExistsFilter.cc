// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/AtomExistsFilter.cc
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)


//Unit Headers
#include <protocols/drug_design/AtomExistsFilter.hh>
#include <protocols/drug_design/AtomExistsFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace drug_design {

static basic::Tracer TR( "protocols.drug_design.AtomExistsFilter" );

protocols::filters::FilterOP
AtomExistsFilterCreator::create_filter() const { return protocols::filters::FilterOP( new AtomExistsFilter ); }

std::string
AtomExistsFilterCreator::keyname() const { return AtomExistsFilter::class_name(); }

void AtomExistsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	AtomExistsFilter::provide_xml_schema( xsd );
}

std::string AtomExistsFilter::class_name() {
	return "AtomExists";
}

void AtomExistsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_refpose_enabled_residue_number, "Residue that the atom is in" )
		+ XMLSchemaAttribute::required_attribute( "atom_name", xs_string, "Name of atom to test for" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Tests if the given atom exists in the protein.",
		attlist );
}

void
AtomExistsFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	residue_ = tag->getOption< std::string >( "residue" );
	atom_name_ = tag->getOption< std::string >( "atom_name" );

	TR << "AtomExists filter for atom "<<atom_name_<<" on residue "<<residue_<<std::endl;
}

bool
AtomExistsFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const exists( compute( pose ) );
	TR << "Does atom '" << atom_name_ << "' exist in residue " << residue_ << "? ";
	if ( exists <= 0.5 ) {
		TR << "NO." << std::endl;
		return false;
	} else {
		TR << "YES." << std::endl;
		return true;
	}
}

void
AtomExistsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	out << "Does atom '" << atom_name_ << "' exist in residue " << residue_ << "? ";
	if ( compute(pose) <= 0.5 ) {
		out << "NO." << std::endl;
	} else {
		out << "YES." << std::endl;
	}
}

core::Real
AtomExistsFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
AtomExistsFilter::compute( core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 || resnum > pose.total_residue() ) {
		TR.Warning << "Attempted to access residue " << residue_ << " in pose with " <<  pose.total_residue() << " residues. Failing filter. " << std::endl;
		return 0.0;
	}
	core::chemical::ResidueType const & restype( pose.residue(resnum).type() );
	TR << "Looking for atom " << atom_name_ << " in residue " << residue_ << ", of type " << restype.name() << std::endl;
	if ( restype.has( atom_name_ ) ) {
		return 1.0;
	} else {
		return 0.0;
	}
}

}
}
