// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/CavityVolumeFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)
/// @modified Vikram K. Mulligan (vmulligan@flatironinstitue.org) -- Moved from devel to protocols.

// Unit headers
#include <protocols/denovo_design/filters/CavityVolumeFilter.hh>
#include <protocols/denovo_design/filters/CavityVolumeFilterCreator.hh>

// Protocol Headers
#include <protocols/denovo_design/calculators/CavityCalculator.hh>

// Project Headers
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>

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


static basic::Tracer TR( "protocols.denovo_design.filters.CavityVolumeFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace filters {




///  ---------------------------------------------------------------------------------
///  CavityVolumeFilter main code:
///  ---------------------------------------------------------------------------------
CavityVolumeFilter::CavityVolumeFilter() :
	Filter( "CavityVolumeFilter" ),
	selector_()
{
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
CavityVolumeFilter::~CavityVolumeFilter() = default;


/// Return a copy of ourselves
protocols::filters::FilterOP
CavityVolumeFilter::clone() const
{
	return utility::pointer::make_shared< CavityVolumeFilter >(*this);
}

protocols::filters::FilterOP
CavityVolumeFilter::fresh_instance() const
{
	return utility::pointer::make_shared< CavityVolumeFilter >();
}

void
CavityVolumeFilter::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & datamap
)
{
	if ( tag->hasOption( "selector" ) ) {
		std::string const & selectorname = tag->getOption< std::string >( "selector" );
		try {
			selector_ = core::select::residue_selector::get_residue_selector(selectorname, datamap);
		} catch (utility::excn::Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' from the Datamap from CavityVolumeFilter::parse_my_tag.\n";
			error_msg << e.msg();
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
	}

	// add selector from sub tags
	for ( auto itag = tag->getTags().begin();
			itag != tag->getTags().end(); ++itag ) {
		if ( selector_ ) {
			std::stringstream error_msg;
			error_msg << "Residue selector can either be specified via name or subtag, but not both in CavityVolumeFilter." << std::endl;
			throw CREATE_EXCEPTION(utility::excn::Exception,  error_msg.str() );
		}
		selector_ = core::select::residue_selector::ResidueSelectorFactory::get_instance()->new_residue_selector(
			(*itag)->getName(),
			(*itag),
			datamap
		);
	}

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

	if ( ! selector_ ) {
		basic::MetricValue< core::Real > total_volume;
		pose.metric( "CavityCalculator", "volume", total_volume );
		TR << "CavityCalculator returned " << total_volume.value() << std::endl;
		return total_volume.value();
	}

	debug_assert( selector_ );
	core::select::residue_selector::ResidueSubset const & subset = selector_->apply( pose );

	basic::MetricValue< utility::vector1< core::scoring::packstat::CavityBallCluster > > clustermetric;
	pose.metric( "CavityCalculator", "clusters", clustermetric );
	utility::vector1< core::scoring::packstat::CavityBallCluster > const & clusters = clustermetric.value();
	core::Real total_volume = 0;
	for ( core::Size i=1, end=clusters.size(); i<=end; ++i ) {
		TR << "Cluster " << i << ": volume=" << clusters[i].volume << ", surface_area=" << clusters[i].surface_area << ", surface_accessibility=" << clusters[i].surface_accessibility << ", center=" << clusters[i].center.x() << "," << clusters[i].center.y() << "," << clusters[i].center.z() << std::endl;
		// only report his cluster if it is close to a selected residue
		for ( core::Size r=1, endr=subset.size(); r <= endr; ++r ) {
			if ( !subset[r] ) {
				continue;
			}
			core::Real const dist = pose.residue(r).nbr_atom_xyz().distance( clusters[i].center );
			//TODO: find way to get cluster radius and distance
			if ( dist < 6.0 ) {
				total_volume += clusters[i].volume;
			}
		}
	}
	return total_volume;
}

/// @brief Does the CavityVolume Filtering
bool
CavityVolumeFilter::apply( core::pose::Pose const & pose ) const
{
	report_sm( pose );
	return true;
}

std::string CavityVolumeFilter::name() const {
	return class_name();
}

std::string CavityVolumeFilter::class_name() {
	return "CavityVolume";
}

void CavityVolumeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	core::select::residue_selector::attributes_for_parse_residue_selector(attlist, "selector");

	protocols::filters::xsd_type_definition_w_attributes(
		xsd, class_name(),
		"Uses Will Sheffler's packing code (packstat) to estimate the total volume of "
		"intra-protein voids. The value returned is the sum of volumes of the computed "
		"cavities in Angstroms 3. A value of 20 is approximately equal to the volume "
		"of a carbon atom. This filter currently has no options or threshold, and "
		"currently always returns true, but that is likely to change in the future. "
		"This calculation of cavity volume is inherently stochastic (packstat is stochastic). "
		"Therefore, setting the temperature to be permissive for the variation in "
		"results with the same input or using the average value of many decoys "
		"(nstruct) is recommended. The filter should not count cavities that are exposed to solvent.",
		attlist );
}

/// @brief Does this filter indicate that it is unpublished (and, by extension, that the author should be
/// included in publications resulting from it)?  Yes, it does.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
bool
CavityVolumeFilter::filter_is_unpublished() const {
	return true;
}

/// @brief Provide a list of authors and their e-mail addresses, as strings.
/// @returns A pairs of ("Tom Linsky", "tlinsky@gmail.com").
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >
CavityVolumeFilter::provide_authorship_info_for_unpublished() const {
	return utility::vector1< basic::citation_manager::UnpublishedModuleInfoCOP >{
		utility::pointer::make_shared< basic::citation_manager::UnpublishedModuleInfo >(
		class_name(), basic::citation_manager::CitedModuleType::Filter,
		"Tom Linsky", "Institute for Protein Design, University of Washington", "tlinsky@gmail.com",
		"Created this filter."
		)
		};
}

std::string CavityVolumeFilterCreator::keyname() const {
	return CavityVolumeFilter::class_name();
}

protocols::filters::FilterOP
CavityVolumeFilterCreator::create_filter() const {
	return utility::pointer::make_shared< CavityVolumeFilter >();
}

void CavityVolumeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CavityVolumeFilter::provide_xml_schema( xsd );
}



} // namespace filters
} // namespace denovo_design
} // namespace protocols


