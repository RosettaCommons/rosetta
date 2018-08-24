// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/DensityFitResidueSelector.hh
/// @brief  Select residues that have a good fit to density.  Can invert to select bad fit.
///
/// - Original Code and Logic - Reference (eLife 2016, Dimaio)
/// @author Ray Wang (wangyr@uw.edu)
/// @author Frank Dimaio (fdimaio@gmail.com)
///
///  - Logic into Residue Selector
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Sebastian Rämisch (raemisch@scripps.edu)

// Unit headers
#include <core/select/residue_selector/DensityFitResidueSelector.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
#include <core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.hh>
#include <core/simple_metrics/util.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/ref_pose.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>
#include <numeric/zscores.hh>
#include <utility/pointer/memory.hh>

// C++ headers
#include <utility/assert.hh>
#include <cmath>

// XSD Includes
#include <core/select/residue_selector/util.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.select.residue_selector.DensityFitResidueSelector" );

namespace core {
namespace select {
namespace residue_selector {

using namespace core::simple_metrics;
using namespace core::simple_metrics::per_residue_metrics;

/// @brief Constructor.
///
DensityFitResidueSelector::DensityFitResidueSelector():
	core::select::residue_selector::ResidueSelector()
{
	den_fit_metric_ = utility::pointer::make_shared< PerResidueDensityFitMetric >();
}

DensityFitResidueSelector::DensityFitResidueSelector( PerResidueDensityFitMetricCOP den_fit_metric ):
	core::select::residue_selector::ResidueSelector()
{
	den_fit_metric_ = utility::pointer::make_shared< PerResidueDensityFitMetric >(*den_fit_metric);
}

/// @brief Destructor.
///
DensityFitResidueSelector::~DensityFitResidueSelector() {}

/// @brief Copy Constructor.
DensityFitResidueSelector::DensityFitResidueSelector(DensityFitResidueSelector const & src):
	score_cut_(src.score_cut_),
	invert_(src.invert_)
{
	den_fit_metric_ = utility::pointer::make_shared< PerResidueDensityFitMetric >( *src.den_fit_metric_);
}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
DensityFitResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new DensityFitResidueSelector(*this) );
}

void
DensityFitResidueSelector::set_den_fit_metric(simple_metrics::per_residue_metrics::PerResidueDensityFitMetricCOP den_fit_metric){
	den_fit_metric_ = utility::pointer::make_shared< PerResidueDensityFitMetric >( *den_fit_metric );
}

PerResidueDensityFitMetricOP
DensityFitResidueSelector::get_den_fit_metric(){
	return den_fit_metric_;
}

///@brief Set up options to use any data already stored in the SM cache through its apply method.
void
DensityFitResidueSelector::set_cache_options( bool use_sm_cache, std::string const & prefix, std::string const & suffix, bool fail_on_missing_cache)
{
	use_cache_ = use_sm_cache;
	prefix_ = prefix;
	suffix_ = suffix;
	fail_on_missing_cache_ = fail_on_missing_cache;
}

/// @brief XML parse.
/// @details Parse RosettaScripts tags and set up this mover.
void
DensityFitResidueSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap )
{
	set_score_cut( tag->getOption< core::Real >("cutoff", score_cut_) );
	set_invert( tag->getOption< bool >("invert", invert_) );

	// Setup Cache Options
	bool use_cache = tag->getOption< bool >("use_cache", use_cache_);
	std::string prefix = tag->getOption< std::string >("prefix", prefix_);
	std::string suffix = tag->getOption< std::string >("suffix", suffix_);
	bool fail_on_missing_cache = tag->getOption< bool >("fail_on_missing_cache", fail_on_missing_cache_);

	set_cache_options(use_cache, prefix, suffix, fail_on_missing_cache);

	if ( tag->hasOption("den_fit_metric") ) {
		SimpleMetricCOP metric = get_metric_from_datamap_and_subtags(tag, datamap, "den_fit_metric");
		PerResidueDensityFitMetric const & den_fit = dynamic_cast<PerResidueDensityFitMetric  const & >( *metric );
		den_fit_metric_ = PerResidueDensityFitMetricOP ( new PerResidueDensityFitMetric(den_fit));

		if ( den_fit_metric_->name_static() != "PerResidueDensityFitMetric" ) {
			utility_exit_with_status("DensityFitResidueSelector requires metric type of PerResidueDensityFitMetric");
		}
	} else {
		if ( tag->hasOption("mixed_sliding_window") ) {
			den_fit_metric_->set_mixed_sliding_window( tag->getOption< bool >("mixed_sliding_window") );
		}
		if ( tag->hasOption("sliding_window_size") ) {
			den_fit_metric_->set_sliding_window_size( tag->getOption< Size >("sliding_window_size") );
		}
		if ( tag->hasOption("match_res") ) {
			den_fit_metric_->set_match_mode( tag->getOption< bool >("match_res") );
		}
		if ( tag->hasOption("use_selector_as_zscore_mask") ) {
			den_fit_metric_->set_use_selector_as_zscore_mask( tag->getOption< bool >("use_selector_as_zscore_mask"));
		}

		if ( tag->hasOption("residue_selector") ) {
			den_fit_metric_->set_residue_selector(parse_residue_selector( tag, datamap ));
		}
	}


}

std::string DensityFitResidueSelector::get_name() const
{
	return DensityFitResidueSelector::class_name();
}

std::string DensityFitResidueSelector::class_name()
{
	return "DensityFitResidueSelector";
}

void DensityFitResidueSelector::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	//Syntax Example:
	using namespace utility::tag;
	AttributeList attributes;
	attributes
		+ XMLSchemaAttribute::attribute_w_default( "invert",  xsct_rosetta_bool, "Select residues that have a bad density fit instead of those with good density fit.", "false")
		+ XMLSchemaAttribute::attribute_w_default( "cutoff",  xsct_real, "Cutoff of bad match to density", "-.5")




		+ XMLSchemaAttribute::attribute_w_default( "sliding_window_size",  xsct_positive_integer, "Sliding window size for density calculation", "3")
		+ XMLSchemaAttribute::attribute_w_default( "mixed_sliding_window",  xsct_rosetta_bool, "Use a window size of 3 for protein and 1 for glycans.  May skew results.", "false")
		+ XMLSchemaAttribute::attribute_w_default( "match_res",  xsct_rosetta_bool, "Use density correlation instead of a zscore to fit to density", "false")
		+ XMLSchemaAttribute::attribute_w_default( "use_selector_as_zscore_mask",  xsct_rosetta_bool, "Use the selector as true mask to calculate the Zscore.  Otherwise, use it just as a selection for result.  Default true.", "true")


		///Cache
		+ XMLSchemaAttribute("den_fit_metric", xs_string, "Use a previously defined Per-Residue Density Fit Metric (use for caching data")

		+ XMLSchemaAttribute::attribute_w_default("use_cache", xsct_rosetta_bool, "Use any cached SM data instead of recaclulating. Data can be cached by calling RunSimpleMetrics on the metric before using this residue selector. The metric will need to be defined before this selector, perhaps in another section of SIMPLE_METRICS.", "false")

		+ XMLSchemaAttribute("prefix", xs_string, "Prefix for the cached data.")

		+ XMLSchemaAttribute("suffix", xs_string, "Suffix for the cached data.")
		+ XMLSchemaAttribute::attribute_w_default("fail_on_missing_cache", xsct_rosetta_bool, "If we are using cached data, should we fail on no cache?", "false");

	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name(attributes, "A Residue selector mask.  Used to only compute Zscore among a set of residues.  Useful for protein vs glycan density.  Since match_res is NOT a zscore, the selector acts as an AND selector, so we only compute the correlations on this set. " );
	std::string documentation = "Select residues that have a good electron density fit. (Or bad fit using the invert option). Uses internal density tools to do so.  Numbers and cutoffs match well with Coot's density fit analysis tool. Zscore uses weighted sum of density, density-compared-to-neighbors, rama (where applicable) and cart_bonded to compute)  Correlation is same values used to calculate density scores.  Zscore reference is here: eLife 2016, Dimaio";

	xsd_type_definition_w_attributes( xsd, class_name(), documentation, attributes );
}

core::select::residue_selector::ResidueSelectorOP
DensityFitResidueSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new DensityFitResidueSelector );
}

std::string
DensityFitResidueSelectorCreator::keyname() const {
	return DensityFitResidueSelector::class_name();
}

//@brief Provide XSD information, allowing automatic evaluation of bad XML.
void
DensityFitResidueSelectorCreator::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd
) const {
	DensityFitResidueSelector::provide_xml_schema( xsd );
}

void
DensityFitResidueSelector::set_invert( bool invert ){
	invert_ = invert;
}

void
DensityFitResidueSelector::set_score_cut(Real score_cut){
	score_cut_ = score_cut;
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
DensityFitResidueSelector::ResidueSubset
DensityFitResidueSelector::apply( core::pose::Pose const & pose ) const
{
	using namespace core::scoring;
	using namespace numeric;

	//Get Zscore.  Make Cut.  Flip if need be.
	core::conformation::symmetry::SymmetryInfoCOP symminfo = nullptr;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo = SymmConf.Symmetry_Info();
	}

	//This catches both code-level and RS interface.
	if ( den_fit_metric_->match_res_mode() && (( score_cut_ < 0 ) || ( score_cut_ > 1 ) ) ) {
		utility_exit_with_message("Using match_res with a the default cutoff.  This mode gives a correlation score to density and should be between 0 and 1");
	}

	//Calculate the values using the PerResidue Metric (optionally using the set cache)
	std::map < core::Size, core::Real > fit_values = den_fit_metric_->cached_calculate( pose, use_cache_, prefix_, suffix_, fail_on_missing_cache_ );

	//std::map < core::Size, core::Real > fit_values = den_fit_metric_->calculate( pose );

	//Match the values to the cutoffs, create the subset
	utility::vector1< Size > subset( pose.size(), false);
	for ( Size i = 1; i <= pose.size(); ++i ) {
		if ( ! fit_values.count(i) ) continue;
		Real fit_value = fit_values[i];

		TR << pose.pdb_info()->pose2pdb(i) << " " << fit_value << std::endl;
		if ( (! invert_) && fit_value >= score_cut_ ) {
			TR << "Good Fit" << std::endl;
			subset[ i ] = true;
		} else if ( invert_ && fit_value < score_cut_ ) {
			TR <<"Bad Fit" << std::endl;
			subset[i] = true;
		}
	}

	//Correct for symmetry:
	if ( core::pose::symmetry::is_symmetric( pose )  ) {
		for ( Size i = 1; i <= pose.size(); ++i ) {

			if ( symminfo->bb_is_independent(i) ) continue;
			else {
				subset[i] = subset[ symminfo->bb_follows(i)];
			}
		}
	}
	return subset;
}


} //core
} //select
} //residue_selector

#ifdef    SERIALIZATION



template< class Archive >
void
core::select::residue_selector::DensityFitResidueSelector::save( Archive & arc ) const {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( score_cut_ ) );
	arc( CEREAL_NVP( invert_));
	arc( CEREAL_NVP( den_fit_metric_));
	arc( CEREAL_NVP( use_cache_));
	arc( CEREAL_NVP( prefix_));
	arc( CEREAL_NVP( suffix_));
	arc( CEREAL_NVP( fail_on_missing_cache_));

}

template< class Archive >
void
core::select::residue_selector::DensityFitResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( score_cut_ );
	arc( invert_ );
	arc( den_fit_metric_);
	arc( use_cache_);
	arc( prefix_ );
	arc( suffix_ );
	arc( fail_on_missing_cache_);


}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::DensityFitResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::DensityFitResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_DensityFitResidueSelector )
#endif // SERIALIZATION


