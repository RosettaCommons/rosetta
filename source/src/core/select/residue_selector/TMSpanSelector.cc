// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/TMSpanSelector.hh
/// @brief  Select residues within given transmembrane spans in a membrane protein
/// @author Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/select/residue_selector/TMSpanSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Project Headers
#include <core/conformation/Conformation.hh> 
#include <core/conformation/membrane/MembraneInfo.hh> 
#include <core/conformation/membrane/SpanningTopology.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

// C++ headers
#include <utility/assert.hh>
#include <boost/foreach.hpp>

static THREAD_LOCAL basic::Tracer TR( "core.select.residue_selector.TMSpanSelector" );

namespace core {
namespace select {
namespace residue_selector {

using namespace core::select::residue_selector;

TMSpanSelector::TMSpanSelector() :
	core::select::residue_selector::ResidueSelector(),
	all_( false ),
	tm_spans_()
{}

TMSpanSelector::TMSpanSelector(
	bool select_all,
	utility::vector1< core::Size > tm_spans
) : core::select::residue_selector::ResidueSelector(),
	all_( select_all ),
	tm_spans_( tm_spans )
{}

TMSpanSelector::~TMSpanSelector() {}

ResidueSubset
TMSpanSelector::apply( core::pose::Pose const & pose ) const
{

	using namespace core::conformation::membrane;
	
	// Check that the pose is a membrane protein
	if ( !pose.conformation().is_membrane() ||
		  pose.conformation().membrane_info()->spanning_topology()->nspans() == 0 ) {
		utility_exit_with_message( "Cannot select TM spans on a non-membrane pose or pose without any transmembrane spanning segments" );
	}
	
	// Grab the spanning topology from the pose for convenience
	SpanningTopologyOP topology( pose.conformation().membrane_info()->spanning_topology() );
	utility::vector1< core::Size > tm_spans = all_tm_spans();
	
	// if user specifies all segments, select them based on the pose
	if ( all() ) {
		for ( core::Size ii = 1; ii <= topology->nspans(); ++ii ) {
			tm_spans.push_back( ii );
		}
	// Exit if no transmembrane spans were selected
	} else if ( tm_spans.size() == 0 ) {
		utility_exit_with_message( "Cannot perform selection operation with zero transmembrane spans" );
	}
	
	// Check that the selected transmembrane spans are valid according to the spanning topology
	check_valid_tm_selection( pose );
	
	// Add the residues contained in each selected TM span to the residue subset
	ResidueSubset selected_tmregions( pose.total_residue(), false );
	for ( core::Size ii = 1; ii <= tm_spans.size(); ++ii ) {
		for ( core::Size jj = topology->span( tm_spans[ii] )->start(); jj <= topology->span( tm_spans[ii] )->end(); ++jj ) {
			selected_tmregions[jj] = true;
		}
	}
	
	return selected_tmregions;
}

void
TMSpanSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ) {
	
	// Check whether the user wants to select all residues
	if ( tag->hasOption( "select_all" ) ) {
		all( tag->getOption< bool >( "select_all" ) );
	}
	
	// If not all residues selected, read in selected TM segments
	if ( !all() && tag->hasOption( "tm_segments" ) ) {

		std::string unparsed_segments = tag->getOption< std::string >( "tm_segments" );
		if ( unparsed_segments != "" ) {
			utility::vector1< std::string > const segments( utility::string_split( unparsed_segments, ',' ) );
			BOOST_FOREACH ( std::string const key, segments ) {
				Size tm_segment( utility::string2int( key ) );
				tm_spans_.push_back( tm_segment );
			}
		}
	}
	
	// If the user selected all and wants to specify segments, turn them away
	if ( all() && tag->hasOption( "tm_segments" ) ) {
		utility_exit_with_message( "Cannot select all TM segments AND specify specific TM segments" );
	}
}

std::string TMSpanSelector::get_name() const {
	return TMSpanSelector::class_name();
}

ResidueSelectorOP
TMSpanSelector::clone() const {
	return ResidueSelectorOP( new TMSpanSelector(*this) );
}

std::string TMSpanSelector::class_name() {
	return "TMSpanSelector";
}

void
TMSpanSelector::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	AttributeList attributes;
	attributes.push_back( XMLSchemaAttribute( "select_all",  xs_boolean ));
	attributes.push_back( XMLSchemaAttribute( "tm_segments", xs_string ));
	xsd_type_definition_w_attributes( xsd, class_name(), attributes );
}

/// @brief Check that the transmembrane spans selected from input are present
/// in the pose's spanning topology object
void
TMSpanSelector::check_valid_tm_selection( core::pose::Pose const & pose ) const {

	for ( core::Size ii = 1; ii <= tm_spans_.size(); ++ii ) {
		if ( tm_spans_[ii] > pose.conformation().membrane_info()->spanning_topology()->nspans() ) {
			utility_exit_with_message( "Transmembrane span selected by user out of bounds of spans desigated in the pose" );
		}
	}
}
	
ResidueSelectorOP
TMSpanSelectorCreator::create_residue_selector() const {
	return ResidueSelectorOP( new TMSpanSelector );
}

std::string
TMSpanSelectorCreator::keyname() const {
	return TMSpanSelector::class_name();
}

void
TMSpanSelectorCreator::provide_selector_xsd( utility::tag::XMLSchemaDefinition & xsd ) const {
	TMSpanSelector::provide_selector_xsd( xsd );
}


} //core
} //select
} //residue_selector
