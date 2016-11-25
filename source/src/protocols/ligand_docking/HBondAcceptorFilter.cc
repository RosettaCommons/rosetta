// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/HBondAcceptorFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/HBondAcceptorFilter.hh>
#include <protocols/ligand_docking/HBondAcceptorFilterCreator.hh>


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

static THREAD_LOCAL basic::Tracer hbond_acceptor_tracer( "protocols.ligand_docking.HBondAcceptorFilter" );

bool
HBondAcceptorFilter::apply( core::pose::Pose const & pose ) const {
	assert(chain_.size()==1 );
	assert(hbond_acceptor_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const begin = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	if ( core::pose::num_hbond_acceptors(begin,end,pose) > hbond_acceptor_limit_ ) {
		hbond_acceptor_tracer<< "Reached hbond acceptor limit"<< std::endl;
		return false;
	}
	return true;
}

void
HBondAcceptorFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{

	if ( tag->getName() != "HBondAcceptor" ) {
		hbond_acceptor_tracer << " received incompatible Tag " << tag << std::endl;
		assert(false);
		return;
	}
	if ( ! (tag->hasOption("chain") && tag->hasOption("hbond_acceptor_limit") ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("HBondAcceptor filter needs a 'chain' and an 'hbond_acceptor_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	hbond_acceptor_limit_ = tag->getOption<core::Size>("hbond_acceptor_limit");

}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP HBondAcceptorFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HBondAcceptorFilter ); }

// XRW TEMP std::string
// XRW TEMP HBondAcceptorFilterCreator::keyname() const { return "HBondAcceptor"; }

std::string HBondAcceptorFilter::name() const {
	return class_name();
}

std::string HBondAcceptorFilter::class_name() {
	return "HBondAcceptor";
}

void HBondAcceptorFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "chain", xsct_char, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "hbond_acceptor_limit", xsct_non_negative_integer, "XRW TO DO" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string HBondAcceptorFilterCreator::keyname() const {
	return HBondAcceptorFilter::class_name();
}

protocols::filters::FilterOP
HBondAcceptorFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new HBondAcceptorFilter );
}

void HBondAcceptorFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HBondAcceptorFilter::provide_xml_schema( xsd );
}



} // ligand_docking
} // protocols
