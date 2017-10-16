// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ligand_docking/HeavyAtomFilter.cc
/// @brief Find packing defects at an interface using packstat score terms
/// @author Jacob Corn (jecorn@u.washington.edu)

// Unit headers
#include <protocols/ligand_docking/HeavyAtomFilter.hh>
#include <protocols/ligand_docking/HeavyAtomFilterCreator.hh>


#include <protocols/filters/Filter.hh>
// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <utility/vector0.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace ligand_docking {

static THREAD_LOCAL basic::Tracer heavy_atom_tracer( "protocols.ligand_docking.HeavyAtomFilter" );

bool
HeavyAtomFilter::apply( core::pose::Pose const & pose ) const {
	debug_assert(chain_.size()==1 );
	debug_assert(heavy_atom_limit_ >0 );
	core::Size const chain_id= core::pose::get_chain_id_from_chain(chain_, pose);
	core::Size const start = pose.conformation().chain_begin(chain_id);
	core::Size const end = pose.conformation().chain_end(chain_id);

	if ( core::pose::num_heavy_atoms(start,end,pose) > heavy_atom_limit_ ) {
		heavy_atom_tracer<< "Reached heavy atom limit"<< std::endl;
		return false;
	}
	return true;
}

void
HeavyAtomFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & )
{
	if ( ! (tag->hasOption("chain") && tag->hasOption("heavy_atom_limit") ) ) {
		throw utility::excn::EXCN_RosettaScriptsOption("HeavyAtom filter needs a 'chain' and a 'heavy_atom_limit' option");
	}
	chain_ = tag->getOption<std::string>("chain");
	heavy_atom_limit_ = tag->getOption<core::Size>("heavy_atom_limit");

}

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP HeavyAtomFilterCreator::create_filter() const { return protocols::filters::FilterOP( new HeavyAtomFilter ); }

// XRW TEMP std::string
// XRW TEMP HeavyAtomFilterCreator::keyname() const { return "HeavyAtom"; }

std::string HeavyAtomFilter::name() const {
	return class_name();
}

std::string HeavyAtomFilter::class_name() {
	return "HeavyAtom";
}

void HeavyAtomFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "chain", xsct_char, "XRW TO DO" )
		+ XMLSchemaAttribute::required_attribute( "heavy_atom_limit", xsct_non_negative_integer, "XRW TO DO" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string HeavyAtomFilterCreator::keyname() const {
	return HeavyAtomFilter::class_name();
}

protocols::filters::FilterOP
HeavyAtomFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new HeavyAtomFilter );
}

void HeavyAtomFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	HeavyAtomFilter::provide_xml_schema( xsd );
}



} // ligand_docking
} // protocols
