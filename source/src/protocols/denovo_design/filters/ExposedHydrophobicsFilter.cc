// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/ExposedHydrophobicsFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilter.hh>
#include <protocols/denovo_design/filters/ExposedHydrophobicsFilterCreator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <protocols/moves/DsspMover.hh>
#include <core/select/util/SelectResiduesByLayer.hh>

// Basic Headers
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

// ObjexxFCL Headers

//C++ Headers


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif


static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.ExposedHydrophobicsFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {

std::string const ExposedHydrophobicsFilter::hydrophobic_residues_ = "MLIVFAWY";

// XRW TEMP std::string
// XRW TEMP ExposedHydrophobicsFilterCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ExposedHydrophobicsFilter::class_name();
// XRW TEMP }

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ExposedHydrophobicsFilterCreator::create_filter() const {
// XRW TEMP  return protocols::filters::FilterOP( new ExposedHydrophobicsFilter() );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ExposedHydrophobicsFilter::class_name()
// XRW TEMP {
// XRW TEMP  return "ExposedHydrophobics";
// XRW TEMP }

///  ---------------------------------------------------------------------------------
///  ExposedHydrophobicsFilter main code:
///  ---------------------------------------------------------------------------------
ExposedHydrophobicsFilter::ExposedHydrophobicsFilter() :
	Filter( "ExposedHydrophobicsFilter" ),
	threshold_( -1.0 ),
	sasa_cutoff_( 20.0 )
{}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
ExposedHydrophobicsFilter::~ExposedHydrophobicsFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
ExposedHydrophobicsFilter::clone() const
{
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter(*this) );
}

protocols::filters::FilterOP
ExposedHydrophobicsFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter() );
}

void
ExposedHydrophobicsFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	sasa_cutoff_ = tag->getOption< core::Real >( "sasa_cutoff", sasa_cutoff_ );
	threshold_ = tag->getOption< core::Real >( "threshold", threshold_ );
}

std::string
ExposedHydrophobicsFilter::get_name() const
{
	return "ExposedHydrophobics";
}

void
ExposedHydrophobicsFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "ExposedHydrophobicsFilter returning " << compute( pose ) << std::endl;
}

core::Real
ExposedHydrophobicsFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
ExposedHydrophobicsFilter::compute( core::pose::Pose const & pose ) const
{
	core::select::util::SelectResiduesByLayer srbl( true, true, true );
	core::Real const surface_sasa( 40 );
	srbl.initialize( sasa_cutoff_, surface_sasa );
	// do dssp on copy of pose to get secondary structure
	core::pose::Pose posecopy( pose );
	moves::DsspMover dssp;
	dssp.apply( posecopy );
	srbl.compute( pose, pose.secstruct() );
	// need to brainstorm the best way of calculating this. There is a cutoff where we don't care. We want to weight bad ones highly. Therefore, maybe a temperature-weighted boltzmann sum is appropriate.
	core::Size residue_count( 0 );
	core::Real sum( 0.0 );
	for ( core::Size i=1; i<=pose.size(); ++i ) {
		//check to see if this residue is a protein residue, and if it's found in the hydrophobic residue list
		if ( pose.residue( i ).is_protein() &&
				hydrophobic_residues_.find( pose.residue( i ).name1() ) != std::string::npos ) {
			++residue_count;
			// a residue only counts as exposed if it has SASA >= sasa_cutoff_
			if ( srbl.rsd_sasa( i ) >= sasa_cutoff_ ) {
				sum += srbl.rsd_sasa( i ) - sasa_cutoff_;
			}
			TR.Debug << "Resi " << pose.residue( i ).name3() << i << " = " << srbl.rsd_sasa( i ) << ", total sum=" << sum << ", count=" << residue_count << std::endl;
		}
	}
	TR << "ExposedHydrophobics value=" << sum << std::endl;
	return sum;
}

/// @brief Does the ExposedHydrophobics Filtering
bool
ExposedHydrophobicsFilter::apply( core::pose::Pose const & pose ) const
{
	core::Real const value( compute( pose ) );
	if ( threshold_ < 0 || value <= threshold_ ) {
		return true;
	} else {
		return false;
	}
}

std::string ExposedHydrophobicsFilter::name() const {
	return class_name();
}

std::string ExposedHydrophobicsFilter::class_name() {
	return "ExposedHydrophobics";
}

void ExposedHydrophobicsFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute(
		"sasa_cutoff", xsct_real,
		"If a residue has SASA lower than this value, it is considered buried "
		"and does not affect the score returned by the ExposedHydrophobics filter." )
		+ XMLSchemaAttribute(
		"threshold", xsct_real,
		"If a protein has an ExposedHydrophobics total score below this value, "
		"it passes the filter. If a negative threshold is specified, "
		"the filter will always pass." );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Computes the SASA for each hydrophobic residue (A, F, I, M, L, W, V, Y). "
		"The score returned reflects both the number of solvent-exposed hydrophobic "
		"residues and the degree to which they are exposed. The score is calculated as "
		"follows. For each hydrophobic residue, if the SASA is above a certain cutoff "
		"value (default=20), then the value of ( SASA - sasa_cutoff ) is added to the "
		"calculated score. The filter passes if the calculated score is less than the "
		"user-specified threshold.",
		attlist );
}

std::string ExposedHydrophobicsFilterCreator::keyname() const {
	return ExposedHydrophobicsFilter::class_name();
}

protocols::filters::FilterOP
ExposedHydrophobicsFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ExposedHydrophobicsFilter );
}

void ExposedHydrophobicsFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ExposedHydrophobicsFilter::provide_xml_schema( xsd );
}



} // namespace filters
} // namespace denovo_design
} // namespace protocols
