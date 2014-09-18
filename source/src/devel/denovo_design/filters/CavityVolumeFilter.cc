// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/devel/denovo_design/filters/CavityVolumeFilter.cc
/// @brief Tom's Denovo design protocol
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <devel/denovo_design/filters/CavityVolumeFilter.hh>
#include <devel/denovo_design/filters/CavityVolumeFilterCreator.hh>

// Protocol Headers
#include <devel/denovo_design/calculators/CavityCalculator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/toolbox/SelectResiduesByLayer.hh>

// Basic Headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>

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


static thread_local basic::Tracer TR( "devel.denovo_design.filters.CavityVolumeFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
namespace filters {

std::string
CavityVolumeFilterCreator::keyname() const
{
	return CavityVolumeFilterCreator::filter_name();
}

protocols::filters::FilterOP
CavityVolumeFilterCreator::create_filter() const {
	return new CavityVolumeFilter();
}

std::string
CavityVolumeFilterCreator::filter_name()
{
	return "CavityVolume";
}

///  ---------------------------------------------------------------------------------
///  CavityVolumeFilter main code:
///  ---------------------------------------------------------------------------------
CavityVolumeFilter::CavityVolumeFilter() :
	Filter( "CavityVolumeFilter" )
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
CavityVolumeFilter::~CavityVolumeFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
CavityVolumeFilter::clone() const
{
	return new CavityVolumeFilter(*this);
}

protocols::filters::FilterOP
CavityVolumeFilter::fresh_instance() const
{
	return new CavityVolumeFilter();
}

void
CavityVolumeFilter::parse_my_tag(
	utility::tag::TagCOP const,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
}

std::string
CavityVolumeFilter::get_name() const
{
	return "CavityVolume";
}

void
CavityVolumeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	out << "Total Cavity Volume is:" << compute( pose ) << std::endl;
}

core::Real
CavityVolumeFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
CavityVolumeFilter::compute( core::pose::Pose const & pose ) const
{
	// check for calculator; create if it doesn't exist
	if ( ! core::pose::metrics::CalculatorFactory::Instance().check_calculator_exists( "CavityCalculator" ) ) {
		calculators::CavityCalculator calculator;
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "CavityCalculator", calculator.clone() );
	}

	basic::MetricValue< core::Real > total_volume;
	pose.metric( "CavityCalculator", "volume", total_volume );
	TR << "CavityCalculator returned " << total_volume.value() << std::endl;

	return total_volume.value();
}

/// @brief Does the CavityVolume Filtering
bool
CavityVolumeFilter::apply( core::pose::Pose const & pose ) const
{
	report_sm( pose );
	return true;
}


} // namespace filters
} // namespace denovo_design
} // namespace devel


