// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/simple_metrics/metrics/ShapeSimilarityMetric.cc
/// @brief    Simple metric for ShapeSimilarityCalculator
/// @details  The code modifies Luki Goldschmidt's
///           implementation of Lawrence & Coleman shape complementarity calculator
///           to allow for the comparison of similar surface shapes.
/// @author   Andreas Scheck (andreas.scheck@epfl.ch)

// Unit Headers
#include <core/simple_metrics/metrics/ShapeSimilarityMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>
#include <core/simple_metrics/util.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/ref_pose.hh>
#include <core/scoring/sc/ShapeSimilarityCalculator.hh>
#include <core/select/residue_selector/ResidueVector.hh>
#include <core/select/residue_selector/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>

// Parser headers
#include <utility/tag/Tag.hh>

// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#include <core/select/residue_selector/ResidueSelector.hh> // AUTO IWYU For ResidueSelector


//// C++ headers
static basic::Tracer tr( "core.simple_metrics.metrics.ShapeSimilarityMetric" );

namespace core {
namespace simple_metrics {
namespace metrics {

/// @brief default constructor
ShapeSimilarityMetric::ShapeSimilarityMetric():
	core::simple_metrics::RealMetric(),
	quick_( false ),
	verbose_( false ),
	reference_pose_( /* NULL */ ),
	selector1_(),
	selector2_(),
	write_int_area_( false ),
	dist_weight_( 0.5 ),
	median_( false )
{}


/// @brief Compute the overall shape similarity of two surfaces
core::Real
ShapeSimilarityMetric::calculate( pose::Pose const & pose ) const {
	try {
		scoring::sc::ShapeSimilarityCalculator ssc;

		if ( !ssc.Init() ) {
			throw CREATE_EXCEPTION(EXCN_InitFailed, "");
		}

		if ( quick_ ) {
			ssc.settings.density = 5.0;
		}

		ssc.settings.weight = dist_weight_;
		ssc.median = median_;

		ssc.Reset(); // this may not be needed anymore, but I'm leaving it here for safety
		core::Real final_ss = 0;

		if ( !reference_pose_ ) {
			tr.Error << "reference pose error" << std::endl;
			throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
		}

		if ( selector1_ && selector2_ ) {
			if ( !ssc.CalcSs(pose, *reference_pose_, selector1_, selector2_) ) {
				tr.Error << "calc error" << std::endl;
				throw CREATE_EXCEPTION(EXCN_CalcFailed, "");
			}
		}
		final_ss = ssc.CalcSs(pose, *reference_pose_, selector1_, selector2_);

		return final_ss;

	} catch( EXCN_InitFailed const & ) {
		tr.Error << "Issue initializing shape similarity calculator" << std::endl;
		return -2;
	} catch( EXCN_CalcFailed const & ) {
		tr.Error << "Issue running shape similarity calculator" << std::endl;
		return -2;
	}
}


/// @brief parse xml
void
ShapeSimilarityMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map )
{
	SimpleMetric::parse_base_tag( tag );
	verbose_ = tag->getOption<core::Size>( "verbose", false );
	quick_ = tag->getOption<core::Size>( "quick", false );
	write_int_area_ = tag->getOption<bool>( "write_int_area", false );
	dist_weight_ = tag->getOption<Real>( "dist_weight", 0.5);
	median_ = tag->getOption<bool>( "median", false );


	std::string const selector1name = tag->getOption< std::string >( "residue_selector1", "" );
	if ( !selector1name.empty() ) selector1_ = core::select::residue_selector::get_residue_selector( selector1name, data_map );

	std::string const selector2name = tag->getOption< std::string >( "residue_selector2", "" );
	if ( !selector2name.empty() ) selector2_ = core::select::residue_selector::get_residue_selector( selector2name, data_map );

	if ( quick_ ) {
		tr.Info << "Calculating shape similarity in quick mode with less accuracy." << std::endl;
	}

	if ( tag->hasOption("dist_weight") ) {
		tr.Info << "Weighting the distance of sampled dots with " << dist_weight_ << ". (default 0.5)" << std::endl;
	}

	if ( median_ ) {
		tr.Info << "Scoring method set to using the median rather than the mean." << std::endl;
	}

	if ( !selector1_ || !selector2_ ) {
		tr.Error << "Residue range selection is invalied since residues" << (selector2_ ? 1 : 2) << " is empty." << std::endl;
		tr.Error << "A valid residue selection needs to be provided - returning -2 instead." << std::endl;
	}


	if ( tag->hasOption("reference_name") ) {
		if ( tag->getOption<bool>("use_native", true) ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Cannot use use_native with reference_name in ShapeSimilarityMetric.");
		}
		reference_pose_ = core::pose::saved_reference_pose(tag, data_map, "reference_name");
		tr<<"Loaded reference pose: "<<tag->getOption< std::string >( "reference_name" )<<" with "<<reference_pose_->size()<<" residues"<<std::endl;
	} else if ( tag->getOption<bool>("use_native", true) ) {
		if ( data_map.has_resource("native_pose") ) {
			reference_pose_ = core::pose::saved_native_pose(data_map);
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Use native specified with ShapeSimilarityMetric, but no native structure has been specified.");
		}
	}
}

/// @brief Uses residue selectors to set up the ShapeSimilarityCalculator
/// @param[in]  pose Pose to be analyzed
/// @param[out] ssc Initialized, empty ShapeSimilarityCalculator, to which pose residues are added
void
ShapeSimilarityMetric::setup_from_selectors( pose::Pose const & pose, scoring::sc::ShapeSimilarityCalculator & ssc ) const
{
	using core::select::residue_selector::ResidueVector;

	ResidueVector const residues1( selector1_->apply( pose ) );
	ResidueVector const residues2( selector2_->apply( pose ) );

	// Dump information about residues
	if ( tr.Info.visible() ) {
		tr.Info << "Using residues for molecule surface (rosetta numbering):" << std::endl;
		tr.Info << "  Surface 1: ";
		for ( auto r = residues1.begin(); r != residues1.end(); ++r ) {
			tr.Info << (r == residues1.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
		tr.Info << "  Surface 2: ";
		for ( auto r = residues2.begin(); r != residues2.end(); ++r ) {
			tr.Info << (r == residues2.begin() ? "" : ", ") << *r;
		}
		tr.Info << std::endl;
	}

	for ( core::Size r : residues1 ) {
		ssc.AddResidue( 0, pose.residue( r ) );
	}

	for ( core::Size r : residues2 ) {
		ssc.AddResidue( 1, pose.residue( r ) );
	}
}


/// @brief prints results to given tracer in a human-readable format
/// @param[out] tr std::ostream object to write to
/// @param[in]  r  ShapeSimilarityCalculatorResults object containing results
void
ShapeSimilarityMetric::print_ss_result(
	std::ostream & tr,
	ShapeSimilarityCalculatorResults const & r,
	core::Real const nsubs_scalefactor ) const
{
	// Verbose view
	tr << "==================================================" << std::endl;
	tr << std::endl;
	for ( int i = 0; i <= 2; i++ ) {
		if ( i < 2 ) {
			tr << "Molecule " << (i+1) << ":" << std::endl;
		} else {
			tr << "Total/Average for both molecules:" << std::endl;
		}

		tr << "          Total Atoms: " << r.surface[i].nAtoms << std::endl;
		tr << "         Buried Atoms: " << r.surface[i].nBuriedAtoms << std::endl;
		tr << "        Blocked Atoms: " << r.surface[i].nBlockedAtoms << std::endl;
		tr << "           Total Dots: " << r.surface[i].nAllDots << std::endl;
		tr << " Trimmed Surface Dots: " << r.surface[i].nTrimmedDots << std::endl;
		tr << "         Trimmed Area: " << r.surface[i].trimmedArea << " (avg) " << std::endl;
		tr << std::endl;
	}
	tr << std::endl;

	for ( int i = 0; i <= 2; i++ ) {
		if ( i < 2 ) {
			tr << "Molecule " << (i+1) << "->" << ((i+1)%2+1) << ": " << std::endl;
		} else {
			tr << "Average for both molecules:" << std::endl;
		}
		tr << "      Mean Separation: " << r.surface[i].d_mean << std::endl;
		tr << "    Median Separation: " << r.surface[i].d_median << std::endl;
		tr << "    Mean Shape Compl.: " << r.surface[i].s_mean << std::endl;
		tr << "  Median Shape Compl.: " << r.surface[i].s_median << std::endl;
		tr << std::endl;
	}

	tr << "Shape similarity: " << r.sc << std::endl;
	tr << "Interface area: " << r.area << std::endl;
	if ( nsubs_scalefactor != 1 ) {
		tr << "Area per monomer: " << ( (core::Real) r.area / nsubs_scalefactor ) << std::endl ;
	}
	tr << "Interface seperation: " << r.distance << std::endl;
}

void ShapeSimilarityMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::attribute_w_default( "verbose" , xsct_rosetta_bool , "If true, print extra calculation details to the tracer." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "quick" , xsct_rosetta_bool , "If true, do a quicker, less accurate calculation by reducing the density." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "write_int_area" , xsct_rosetta_bool , "If true, write interface area to scorefile." , "false" )
		+ XMLSchemaAttribute( "residue_selector1" , xs_string , "Set which residues are to be considered on the target side using residue_selectors." )
		+ XMLSchemaAttribute( "residue_selector2" , xs_string , "Set which residues are to be considered on the reference side using residue_selectors." )
		+ XMLSchemaAttribute::attribute_w_default( "dist_weight" , xsct_real , "Set the weight for the distance of dots. The higher value, the stronger far dots are penalized." , "0.5" )
		+ XMLSchemaAttribute::attribute_w_default( "median" , xsct_rosetta_bool , "Set the scoring function to use the median rather than the mean." , "false" )
		+ XMLSchemaAttribute::attribute_w_default( "use_native", xsct_rosetta_bool, "Use the native if present on the cmd-line.", "true" ) ;

	core::simple_metrics::xsd_simple_metric_type_definition_w_attributes( xsd, name_static(), "Calculates the shape similarity similar to concepts introduced by Lawrence and Coleman."
		"The two surfaces have to be specified by explicitly providing lists of the residues making up each surface or using residue selectors. Selector1 specifies residues"
		"on the target protein (usually provided with -s) and selector2 the residues on the reference protein (provided either via SavePoseMover or -in:file:native flag)."
		"NO alignment of the protein is perfomed during the comparison.", attlist );
}

std::string ShapeSimilarityMetric::name() const {
	return name_static();
}

std::string ShapeSimilarityMetric::name_static() {
	return "ShapeSimilarityMetric";
}

std::string ShapeSimilarityMetric::metric() const {
	return "ss";
}

void ShapeSimilarityMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ShapeSimilarityMetric::provide_xml_schema( xsd );
}

std::string
ShapeSimilarityMetricCreator::keyname() const {
	return ShapeSimilarityMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
ShapeSimilarityMetricCreator::create_simple_metric() const {
	return SimpleMetricOP( new ShapeSimilarityMetric );
}


} // core
} // simple_metrics
} // metrics
