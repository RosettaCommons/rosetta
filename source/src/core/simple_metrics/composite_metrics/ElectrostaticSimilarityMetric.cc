// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/simple_metrics/composite_metrics/ElectrostaticSimilarityMetric.cc
/// @brief   Simple metric for ElectrostaticSimilarityCalculator
/// @details   The code closely follows Brian Coventry's implementation
///    of electrostatic complementarity with a added features
///    that in turn is based on the method:
///    McCoy, A. J., Epa, V. C., & Colman, P. M. (1997).
///    Electrostatic complementarity at protein/protein interfaces.
///    Journal of molecular biology, 268(2), 570-584.
/// @author  Andreas Scheck (andreas.scheck@epfl.ch)

// Unit Headers
#include <core/simple_metrics/composite_metrics/ElectrostaticSimilarityMetric.hh>
#include <core/simple_metrics/composite_metrics/ElectrostaticComplementarityMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/ref_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/matdes.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>


// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>


//// C++ headers
static basic::Tracer tr( "core.simple_metrics.composite_metrics.ElectrostaticSimilarityMetric" );

namespace core {
namespace simple_metrics {
namespace composite_metrics {

/// @brief default constructor
ElectrostaticSimilarityMetric::ElectrostaticSimilarityMetric():
	CompositeRealMetric(),
	partially_solvated_( scoring::sc::ElectrostaticSimilarityDefaults::PARTIALLY_SOLVATED ),
	report_all_es_( false ),
	reference_pose_( /* NULL */ ),
	selector1_( nullptr ),
	selector2_( nullptr )
{}

/// @brief Compute the overall electrostatic similarity of two surfaces
std::map< std::string, core::Real >
ElectrostaticSimilarityMetric::calculate( pose::Pose const & pose ) const
{

	std::map< std::string, core::Real > results;
	results["avg"] = -2;
	results["p"] = -2;
	results["s"] = -2;
	if ( report_all_es_ ) {
		results["1_p"] = -2;
		results["1_s"] = -2;
		results["2_p"] = -2;
		results["2_s"] = -2;
	}

	try {
		scoring::sc::ElectrostaticSimilarityCalculator esc;

		if ( !esc.Init() ) {
			throw CREATE_EXCEPTION(EXCN_InitFailed, "");
		}

		esc.partially_solvated( partially_solvated_ );
		if ( reference_pose_ ) {
			esc.reference_pose( reference_pose_ );
		} else {
			tr << "reference pose error" << std::endl;
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}

		// selector-based
		if ( selector1_ && selector2_ ) {
			esc.selector1( selector1_ );
			esc.selector2( selector2_ );
			if ( !esc.Calc( pose ) ) {
				tr << "calc error" << std::endl;
				throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
			}
		}

		scoring::sc::ElectrostaticSimilarityResults const & r = esc.GetResults();

		results["p"] = (r.es_0_p + r.es_1_p) / 2;
		results["s"] = (r.es_0_s + r.es_1_s) / 2;
		results["avg"] = ( results["p"] + results["s"] ) / 2;
		if ( report_all_es_ ) {
			results["1_p"] = r.es_0_p;
			results["1_s"] = r.es_0_s;
			results["2_p"] = r.es_1_p;
			results["2_s"] = r.es_1_s;
		}

		return results;

	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing electrostatic similarity calculator" << std::endl;
		return results;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running electrostatic similarity calculator" << std::endl;
		return results;
	}
}

/// @brief parse xml
void
ElectrostaticSimilarityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption( "partially_solvated" ) ) {
		partially_solvated( tag->getOption< bool >( "partially_solvated" ) );
	}
	if ( tag->hasOption( "report_all_es" ) ) {
		report_all_es( tag->getOption< bool >( "report_all_es" ) );
	}
	if ( tag->hasOption( "residue_selector1" ) ) {
		residue_selector1( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "residue_selector1" ), data ) );
	}
	if ( tag->hasOption( "residue_selector2" ) ) {
		residue_selector2( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "residue_selector2" ), data ) );
	}
	if ( tag->hasOption("reference_name") ) {
		if ( tag->getOption<bool>("use_native", true) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot use use_native with reference_name in ElectrostaticSimilarityMetric.");
		}
		reference_pose_ = core::pose::saved_reference_pose(tag, data, "reference_name");
		tr<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<reference_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption<bool>("use_native", true) ) {
		if ( data.has_resource("native_pose") ) {
			reference_pose_ = core::pose::saved_native_pose(data);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Use native specified with ElectrostaticSimilarityMetric, but no native structure has been specified.");
		}
	}
}

/// @brief Uses residue selectors to set up the ElectrostaticSimilarityCalculator
/// @param[in]  pose Pose to be analyzed
/// @param[out] esc Initialized, empty ElectrostaticSimilarityCalculator, to which pose residues are added
void
ElectrostaticSimilarityMetric::setup_from_selectors( pose::Pose const & pose, scoring::sc::ElectrostaticSimilarityCalculator & esc ) const
{

	select::residue_selector::ResidueSubset const residues1( selector1_->apply( pose ) );
	select::residue_selector::ResidueSubset const residues2( selector2_->apply( pose ) );

	// Dump information about residues
	if ( tr.Info.visible() ) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		bool first = true;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( !residues1[i] ) continue;
			tr.Info << (first ? "" : ", ") << i;
			first = false;
		}
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		first = true;
		for ( Size i = 1; i <= pose.size(); i++ ) {
			if ( !residues2[i] ) continue;
			tr.Info << (first ? "" : ", ") << i;
			first = false;
		}
		tr.Info << std::endl;
	}

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !residues1[i] ) continue;
		esc.AddResidue( pose, 0, i );
	}

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !residues2[i] ) continue;
		esc.AddResidue( pose, 1, i );
	}
}

void ElectrostaticSimilarityMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace scoring::sc::ElectrostaticSimilarityDefaults;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "partially_solvated" , xsct_rosetta_bool , "If false, then fully solvated. When partially solvated, the other side of the pose is present but has 0 charge. If fully solvated, the other side of the pose is not present and solvent fills its place." , utility::to_string((int)PARTIALLY_SOLVATED) )
		+ XMLSchemaAttribute::attribute_w_default( "report_all_es" , xsct_rosetta_bool , "Do you want to see the sub correlations? (The normally reported values are the averages of interface 1 and 2)" , "0" )
		+ XMLSchemaAttribute( "residue_selector1" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." )
		+ XMLSchemaAttribute( "residue_selector2" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." )
		+ XMLSchemaAttribute::attribute_w_default( "use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "true" ) ;

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes( xsd, name_static(), "Calculates the electrostatic similarity using APBS based on McCoy et al.", attlist );
}

std::string
ElectrostaticSimilarityMetric::name() const {
	return name_static();
}

std::string
ElectrostaticSimilarityMetric::name_static() {
	return "ElectrostaticSimilarityMetric";
}

std::string
ElectrostaticSimilarityMetric::metric() const {
	return "es";
}

utility::vector1< std::string >
ElectrostaticSimilarityMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back("avg");
	names.push_back("p");
	names.push_back("s");
	if ( report_all_es_ ) {
		names.push_back("1_p");
		names.push_back("1_s");
		names.push_back("2_p");
		names.push_back("2_s");
	}
	return names;
}

void
ElectrostaticSimilarityMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ElectrostaticSimilarityMetric::provide_xml_schema( xsd );
}

std::string
ElectrostaticSimilarityMetricCreator::keyname() const {
	return ElectrostaticSimilarityMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ElectrostaticSimilarityMetricCreator::create_simple_metric() const {
	return SimpleMetricOP( new ElectrostaticSimilarityMetric );
}


} // namespace core
} // namespace simple_metrics
} // namespace composite_metrics
