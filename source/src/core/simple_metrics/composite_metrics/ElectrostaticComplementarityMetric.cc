// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/simple_metrics/composite_metrics/ElectrostaticComplementarityMetric.cc
/// @brief  Simple metric for ElectrostaticComplementarityCalculator
/// @author Brian Coventry (bcov@uw.edu)

// Unit Headers
#include <core/simple_metrics/composite_metrics/ElectrostaticComplementarityMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility headers
#include <utility/exit.hh>
#include <basic/Tracer.hh>

// Parser headers
#include <utility/tag/Tag.hh>

#include <utility/vector1.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>


//// C++ headers
static basic::Tracer tr( "core.simple_metrics.composite_metrics.ElectrostaticComplementarityMetric" );

namespace core {
namespace simple_metrics {
namespace composite_metrics {

// @brief default constructor
ElectrostaticComplementarityMetric::ElectrostaticComplementarityMetric():
	CompositeRealMetric(),
	ignore_radius_( scoring::sc::ElectrostaticComplementarityDefaults::IGNORE_RADIUS ),
	interface_trim_radius_( scoring::sc::ElectrostaticComplementarityDefaults::INTERFACE_TRIM_RADIUS ),
	partially_solvated_( scoring::sc::ElectrostaticComplementarityDefaults::PARTIALLY_SOLVATED ),
	jump_id_( -1 ),
	report_all_ec_( false ),
	selector1_( nullptr ),
	selector2_( nullptr )
{}





std::map< std::string, core::Real >
ElectrostaticComplementarityMetric::calculate( pose::Pose const & pose ) const
{

	std::map< std::string, core::Real > results;
	results["avg"] = -2;
	results["p"] = -2;
	results["s"] = -2;
	if ( report_all_ec_ ) {
		results["1_p"] = -2;
		results["1_s"] = -2;
		results["2_p"] = -2;
		results["2_s"] = -2;
	}

	try {
		scoring::sc::ElectrostaticComplementarityCalculator ecc;

		if ( !ecc.Init( pose ) ) {
			throw CREATE_EXCEPTION(EXCN_InitFailed, "");
		}

		ecc.ignore_radius( ignore_radius_ );
		ecc.interface_trim_radius( interface_trim_radius_ );
		ecc.partially_solvated( partially_solvated_ );

		if ( selector1_ && selector2_ ) {
			// selector-based
			setup_from_selectors( pose, ecc );
			if ( !ecc.Calc( pose ) ) {
				throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
			}
		} else {
			// jump-based
			if ( !ecc.Calc( pose, jump_id_ ) ) {
				throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
			}
		}

		scoring::sc::ElectrostaticComplementarityResults const & r = ecc.GetResults();

		results["p"] = (r.ec_0_p + r.ec_1_p) / 2;
		results["s"] = (r.ec_0_s + r.ec_1_s) / 2;
		results["avg"] = ( results["p"] + results["s"] ) / 2;
		if ( report_all_ec_ ) {
			results["1_p"] = r.ec_0_p;
			results["1_s"] = r.ec_0_s;
			results["2_p"] = r.ec_1_p;
			results["2_s"] = r.ec_1_s;
		}

		return results;

	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing electrostatic complementarity calculator" << std::endl;
		return results;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running electrostatic complementarity calculator" << std::endl;
		return results;
	}
}

void
ElectrostaticComplementarityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data )
{
	SimpleMetric::parse_base_tag( tag );
	if ( tag->hasOption( "ignore_radius" ) ) {
		ignore_radius( tag->getOption< Real >( "ignore_radius" ) );
	}
	if ( tag->hasOption( "interface_trim_radius" ) ) {
		interface_trim_radius( tag->getOption< Real >( "interface_trim_radius" ) );
	}
	if ( tag->hasOption( "partially_solvated" ) ) {
		partially_solvated( tag->getOption< bool >( "partially_solvated" ) );
	}
	if ( tag->hasOption( "jump" ) ) {
		jump_id( tag->getOption< bool >( "jump" ) );
	}
	if ( tag->hasOption( "report_all_ec" ) ) {
		report_all_ec( tag->getOption< bool >( "report_all_ec" ) );
	}
	if ( tag->hasOption( "residue_selector1" ) ) {
		residue_selector1( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "residue_selector1" ), data ) );
	}
	if ( tag->hasOption( "residue_selector2" ) ) {
		residue_selector2( core::select::residue_selector::get_residue_selector( tag->getOption< std::string >( "residue_selector2" ), data ) );
	}
}

/// @brief Uses residue selectors to set up the ShapeComplementarityCalculator
/// @param[in]  pose Pose to be analyzed
/// @param[out] scc Initialized, empty ShapeComplementarityCalculator, to which pose residues are added
void
ElectrostaticComplementarityMetric::setup_from_selectors( pose::Pose const & pose, scoring::sc::ElectrostaticComplementarityCalculator & ecc ) const
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
		ecc.AddResidue( pose, 0, i );
	}

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( !residues2[i] ) continue;
		ecc.AddResidue( pose, 1, i );
	}
}


void ElectrostaticComplementarityMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	using namespace scoring::sc::ElectrostaticComplementarityDefaults;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "ignore_radius" , xsct_real , "How far away from the interface should atoms be ignored. -1 to keep everything. (15 Shows variation of about 0.001. 10 ~ 0.01)" , utility::to_string(IGNORE_RADIUS) )
		+ XMLSchemaAttribute::attribute_w_default( "interface_trim_radius" , xsct_real , "Value for the SC filter. How far from the sasa surface should the molecular dots be trimmed? The original paper used 0 here." , utility::to_string(INTERFACE_TRIM_RADIUS) )
		+ XMLSchemaAttribute::attribute_w_default( "partially_solvated" , xsct_rosetta_bool , "If false, then fully solvated. When partially solvated, the other side of the pose is present but has 0 charge. If fully solvated, the other side of the pose is not present and solvent fills its place." , utility::to_string((int)PARTIALLY_SOLVATED) )
		+ XMLSchemaAttribute::attribute_w_default( "jump" , xs_integer , "Which jump over which to calculate the interface." , "1" )
		+ XMLSchemaAttribute::attribute_w_default( "report_all_ec" , xsct_rosetta_bool , "Do you want to see the sub correlations? (The normally reported values are the averages of interface 1 and 2)" , "0" )
		+ XMLSchemaAttribute( "residue_selector1" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." )
		+ XMLSchemaAttribute( "residue_selector2" , xs_string , "Explicitly set which residues are on each side of the interface using residue_selectors." ) ;

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes( xsd, name_static(), "Calculates the McCoy, Chandana, Colman Electrostatic complementarity using APBS.", attlist );
}

std::string
ElectrostaticComplementarityMetric::name() const {
	return name_static();
}

std::string
ElectrostaticComplementarityMetric::name_static() {
	return "ElectrostaticComplementarityMetric";
}

std::string
ElectrostaticComplementarityMetric::metric() const {
	return "ec";
}

utility::vector1< std::string >
ElectrostaticComplementarityMetric::get_metric_names() const {
	utility::vector1< std::string > names;
	names.push_back("avg");
	names.push_back("p");
	names.push_back("s");
	if ( report_all_ec_ ) {
		names.push_back("1_p");
		names.push_back("1_s");
		names.push_back("2_p");
		names.push_back("2_s");
	}
	return names;
}

void
ElectrostaticComplementarityMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ElectrostaticComplementarityMetric::provide_xml_schema( xsd );
}

std::string
ElectrostaticComplementarityMetricCreator::keyname() const {
	return ElectrostaticComplementarityMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ElectrostaticComplementarityMetricCreator::create_simple_metric() const {
	return SimpleMetricOP( new ElectrostaticComplementarityMetric );
}


} // namespace core
} // namespace simple_metrics
} // namespace composite_metrics
