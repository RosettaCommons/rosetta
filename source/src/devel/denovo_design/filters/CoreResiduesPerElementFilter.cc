// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @file src/devel/denovo_design/filters/CoreResiduesPerElementFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <devel/denovo_design/filters/CoreResiduesPerElementFilter.hh>
#include <devel/denovo_design/filters/CoreResiduesPerElementFilterCreator.hh>

// Project Headers
#include <devel/denovo_design/scoring/SideChainNeighborsEnergy.hh>

// Protocol Headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/select/residue_selector/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

// ObjexxFCL Headers

// C++ Headers

static basic::Tracer TR( "devel.denovo_design.filters.CoreResiduesPerElementFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
namespace filters {

// XRW TEMP std::string
// XRW TEMP CoreResiduesPerElementFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return CoreResiduesPerElementFilter::class_name();
// XRW TEMP }

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP CoreResiduesPerElementFilterCreator::create_filter() const {
// XRW TEMP  return protocols::filters::FilterOP( new CoreResiduesPerElementFilter() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CoreResiduesPerElementFilter::class_name()
// XRW TEMP {
// XRW TEMP  return "CoreResiduesPerElement";
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  CoreResiduesPerElementFilter main code:
///  ---------------------------------------------------------------------------------
CoreResiduesPerElementFilter::CoreResiduesPerElementFilter() :
	Filter( "CoreResiduesPerElementFilter" ),
	core_cutoff_( 5.2 ),
	selector_()
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
CoreResiduesPerElementFilter::~CoreResiduesPerElementFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
CoreResiduesPerElementFilter::clone() const
{
	return protocols::filters::FilterOP( new CoreResiduesPerElementFilter(*this) );
}

protocols::filters::FilterOP
CoreResiduesPerElementFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new CoreResiduesPerElementFilter() );
}

void
CoreResiduesPerElementFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	if ( tag->hasOption( "core_cutoff" ) ) {
		set_core_cutoff( tag->getOption< core::Real >( "core_cutoff" ) );
	}

	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
	if ( selector_ ) {
		TR << "Using residue selector " << std::endl;
	}
}

std::string
CoreResiduesPerElementFilter::get_name() const
{
	return CoreResiduesPerElementFilter::class_name();
}


void
CoreResiduesPerElementFilter::report( std::ostream & out, core::pose::Pose const & ) const
{
	out << " reporting" << std::endl;
}

core::Real
CoreResiduesPerElementFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
CoreResiduesPerElementFilter::compute( core::pose::Pose const & pose ) const
{
	static core::Size const min_helix_len = 4;
	static core::Size const min_strand_len = 3;

	core::pose::Pose work( pose );
	protocols::moves::DsspMover dssp;
	dssp.apply( work );
	std::string const ss = work.secstruct();
	protocols::fldsgn::topology::SS_Info2 ssinfo( ss );

	core::scoring::methods::SideChainNeighborsEnergy scn;
	core::scoring::ScoreFunction sfxn;
	sfxn.set_weight( core::scoring::sidechain_neighbors, 1.0 );
	scn.setup_for_scoring( work, sfxn );

	core::Size badelements = 0;
	if ( selector_ ) {
		// selection is treated as one element
		core::select::residue_selector::ResidueSubset subset = selector_->apply( pose );
		core::Size found = evaluate_element( work, scn, subset );
		TR << "Core/boundary residues found in selected element: " << found << std::endl;
		if ( !found ) {
			++badelements;
		}
	} else {
		core::Size helixcount = 0;
		core::Size strandcount = 0;
		for ( protocols::fldsgn::topology::Strands::const_iterator s=ssinfo.strands().begin(), ends=ssinfo.strands().end(); s!=ends; ++s ) {
			++strandcount;
			debug_assert( (*s)->end() >= (*s)->begin() );
			// skip strands less than 3 residues
			if ( ((*s)->end() - (*s)->begin()) < min_strand_len ) {
				continue;
			}
			core::Size found = evaluate_element( work, scn, (*s)->begin(), (*s)->end() );
			TR << "Core/boundary residues found in strand " << strandcount << ": " << found << std::endl;
			if ( !found ) {
				++badelements;
			}
		}

		for ( protocols::fldsgn::topology::Helices::const_iterator h=ssinfo.helices().begin(), endh=ssinfo.helices().end(); h!=endh; ++h ) {
			++helixcount;
			debug_assert( (*h)->end() >= (*h)->begin() );
			// skip helices less than 4 residues
			if ( ((*h)->end() - (*h)->begin()) < min_helix_len ) {
				continue;
			}
			core::Size found = evaluate_element( work, scn, (*h)->begin(), (*h)->end() );
			TR << "Core/boundary residues found in helix " << helixcount << ": " << found << std::endl;
			if ( !found ) {
				++badelements;
			}
		}
	}
	TR << "pose has " << badelements << " elements without core/boundary residues." << std::endl;
	return badelements;
}

/// @brief Does the CoreResiduesPerElement Filtering
bool
CoreResiduesPerElementFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real badelements = report_sm( pose );
	return (badelements < 0.0001);
}

/// @brief sets the residue selector used to decide which residues to evaluate
void
CoreResiduesPerElementFilter::set_selector( core::select::residue_selector::ResidueSelectorCOP rs )
{
	selector_ = rs;
}

core::Size
CoreResiduesPerElementFilter::evaluate_element(
	core::pose::Pose const & pose,
	core::scoring::methods::SideChainNeighborsEnergy const & scn,
	core::Size const start,
	core::Size const stop ) const
{
	utility::vector1< core::Size > residues;
	for ( core::Size r=start; r<=stop; ++r ) {
		residues.push_back( r );
	}
	return evaluate_element( pose, scn, residues );
}

core::Size
CoreResiduesPerElementFilter::evaluate_element(
	core::pose::Pose const & pose,
	core::scoring::methods::SideChainNeighborsEnergy const & scn,
	core::select::residue_selector::ResidueSubset const & subset ) const
{
	utility::vector1< core::Size > residues;
	for ( core::Size i=1, endi=subset.size(); i<=endi; ++i ) {
		if ( subset[i] ) {
			residues.push_back(i);
		}
	}
	return evaluate_element( pose, scn, residues );
}

core::Size
CoreResiduesPerElementFilter::evaluate_element(
	core::pose::Pose const & pose,
	core::scoring::methods::SideChainNeighborsEnergy const & scn,
	utility::vector1< core::Size > const & residues ) const
{
	// eventually, make these options

	core::Size found = 0;
	for ( utility::vector1< core::Size >::const_iterator r=residues.begin(), endr=residues.end(); r!=endr; ++r ) {
		core::scoring::EnergyMap emap;
		scn.residue_energy( pose.residue(*r), pose, emap );
		TR.Debug << "Res " << *r << " " << pose.residue(*r).name() << " score=" << emap[ core::scoring::sidechain_neighbors ] << std::endl;
		if ( emap[ core::scoring::sidechain_neighbors ] <= -core_cutoff_ ) {
			++found;
		}
	}
	return found;
}

void
CoreResiduesPerElementFilter::set_core_cutoff( core::Real const core_cutoff )
{
	core_cutoff_ = core_cutoff;
}

std::string CoreResiduesPerElementFilter::name() const {
	return class_name();
}

std::string CoreResiduesPerElementFilter::class_name() {
	return "CoreResiduesPerElement";
}

void CoreResiduesPerElementFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist + XMLSchemaAttribute(
		"core_cutoff", xsct_real,
		"XSD XRW: TO DO");

	core::select::residue_selector::attributes_for_parse_residue_selector(attlist);

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"XSD XRW TO DO",
		attlist );
}

std::string CoreResiduesPerElementFilterCreator::keyname() const {
	return CoreResiduesPerElementFilter::class_name();
}

protocols::filters::FilterOP
CoreResiduesPerElementFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new CoreResiduesPerElementFilter );
}

void CoreResiduesPerElementFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CoreResiduesPerElementFilter::provide_xml_schema( xsd );
}


} // namespace filters
} // namespace denovo_design
} // namespace devel
