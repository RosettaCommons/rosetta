// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/MolecularMassFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/MolecularMassFilter.hh>
#include <protocols/ligand_docking/MolecularMassFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer mass_tracer( "protocols.ligand_docking.MolecularMassFilter" );

bool
MolecularMassFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(mass_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	//core::Real mass=0;


	if ( core::pose::mass(start,end,pose) > mass_limit_ ) {
		mass_tracer<< "Reached atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
MolecularMassFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "MolecularMass" ) {
		mass_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("mass_limit") ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("MolecularMass filter needs a 'chain' and an 'mass_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	mass_limit_ = tag->getOption<core::Size>("mass_limit");


}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP MolecularMassFilterCreator::create_filter() const { return protocols::filters::FilterOP( new MolecularMassFilter ); }

// XRW TEMP std::string
// XRW TEMP MolecularMassFilterCreator::keyname() const { return "MolecularMass"; }

std::string MolecularMassFilter::name() const {
	return class_name();
}

std::string MolecularMassFilter::class_name() {
	return "MolecularMass";
}

void MolecularMassFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "chain", xsct_char, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "mass_limit", xsct_non_negative_integer, "XRW TO DO" );
	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string MolecularMassFilterCreator::keyname() const {
	return MolecularMassFilter::class_name();
}

protocols::filters::FilterOP
MolecularMassFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new MolecularMassFilter );
}

void MolecularMassFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	MolecularMassFilter::provide_xml_schema( xsd );
}



} // ligand_docking
} // protocols
