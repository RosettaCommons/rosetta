// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.cc
/// @brief A per-residue metric that will calculate the density fit for each residue using either a correlation or a zscore.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
///
/// - Original Code zscore, matchres code, and Logic - Reference (eLife 2016, Dimaio)
/// @author Ray Wang (wangyr@uw.edu)
/// @author Frank Dimaio (fdimaio@gmail.com)

// Unit headers
#include <core/simple_metrics/per_residue_metrics/PerResidueDensityFitMetric.hh>
#include <core/simple_metrics/simple_metric_creators.hh>

// Core headers
#include <core/simple_metrics/PerResidueRealMetric.hh>
#include <core/simple_metrics/util.hh>

#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/select/residue_selector/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/ref_pose.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>
#include <numeric/zscores.hh>

// XSD Includes
#include <utility/tag/XMLSchemaGeneration.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


static basic::Tracer TR( "core.simple_metrics.per_residue_metrics.PerResidueDensityFitMetric" );


namespace core {
namespace simple_metrics {
namespace per_residue_metrics {

using namespace core::scoring::electron_density;

/////////////////////
/// Constructors  ///
/////////////////////

/// @brief Default constructor
PerResidueDensityFitMetric::PerResidueDensityFitMetric():
	core::simple_metrics::PerResidueRealMetric()
{}

////////////////////////////////////////////////////////////////////////////////
/// @brief Destructor (important for properly forward-declaring smart-pointer members)
PerResidueDensityFitMetric::~PerResidueDensityFitMetric(){}

////////////////////////////////////////////////////////////////////////////////
/// @brief Copy constructor
PerResidueDensityFitMetric::PerResidueDensityFitMetric( PerResidueDensityFitMetric const &  ) = default;

core::simple_metrics::SimpleMetricOP
PerResidueDensityFitMetric::clone() const {
	return core::simple_metrics::SimpleMetricOP(new PerResidueDensityFitMetric( *this ) );

}

std::string
PerResidueDensityFitMetric::name() const {
	return name_static();
}

std::string
PerResidueDensityFitMetric::name_static() {
	return "PerResidueDensityFitMetric";

}
std::string
PerResidueDensityFitMetric::metric() const {
	return "res_density_fit";
}


void
PerResidueDensityFitMetric::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)
{

	SimpleMetric::parse_base_tag( tag );
	PerResidueRealMetric::parse_per_residue_tag( tag, datamap );

	set_match_mode( tag->getOption< bool >("match_res", match_res_) );
	set_use_selector_as_zscore_mask( tag->getOption< bool >("use_selector_as_zscore_mask", use_selector_as_zscore_mask_));
	set_mixed_sliding_window( tag->getOption< bool >("mixed_sliding_window", mixed_sliding_window_) );
	set_sliding_window_size( tag->getOption< Size >("sliding_window_size", sliding_window_size_) );
}

void
PerResidueDensityFitMetric::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	using namespace core::select::residue_selector;

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default( "sliding_window_size",  xsct_positive_integer, "Sliding window size for density calculation", "3")
		+ XMLSchemaAttribute::attribute_w_default( "match_res",  xsct_rosetta_bool, "Use density correlation instead of a zscore to fit to density", "false")
		+ XMLSchemaAttribute::attribute_w_default( "mixed_sliding_window",  xsct_rosetta_bool, "Use a window size of 3 for protein and 1 for glycans.  May skew results.", "false")
		+ XMLSchemaAttribute::attribute_w_default( "use_selector_as_zscore_mask",  xsct_rosetta_bool, "Use the selector as true mask to calculate the Zscore.  Otherwise, use it just as a selection for computation.  Default true.", "true");


	//core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name(attlist, "A Residue selector mask.  Used to only compute Zscore among a set of residues.  Useful for protein vs glycan density.  Since match_res is NOT a zscore, the selector acts as an AND selector, so we only compute the correlations on this set. " );

	std::string documentation = "Calculate either a correlation (match_res true) or zscore related to the fit to the density.  Uses internal density tools to do so.  Numbers and cutoffs match well with Coot's density fit analysis tool. Zscore uses weighted sum of density, density-compared-to-neighbors, rama (where applicable) and cart_bonded to compute)  Correlation is same values used to calculate density scores.  Zscore reference is here: eLife 2016, Dimaio";

	core::simple_metrics::xsd_per_residue_real_metric_type_definition_w_attributes(xsd, name_static(),
		documentation, attlist);
}


void
PerResidueDensityFitMetric::set_sliding_window_size(core::Size window_size){
	sliding_window_size_ = window_size;
}

void
PerResidueDensityFitMetric::set_mixed_sliding_window(bool mixed_sliding_window){
	mixed_sliding_window_ = mixed_sliding_window;
}

void
PerResidueDensityFitMetric::set_match_mode( bool match_mode ){
	match_res_ = match_mode;
}

void
PerResidueDensityFitMetric::set_use_selector_as_zscore_mask(bool selector_as_zscore_mask){
	use_selector_as_zscore_mask_ = selector_as_zscore_mask;
}

void
PerResidueDensityFitMetric::compute_scores(
	pose::Pose & pose,

	std::map< Size, Real > & per_rsd_dens,
	std::map< Size, Real > & per_rsd_nbrdens,
	std::map< Size, Real > & per_rsd_rama,
	std::map< Size, Real > & per_rsd_geometry
) const {

	//Initialize
	per_rsd_dens.clear(); per_rsd_rama.clear();
	per_rsd_nbrdens.clear(); per_rsd_geometry.clear();

	core::conformation::symmetry::SymmetryInfoCOP symminfo = nullptr;
	core::Size n_symm_subunit = 1;

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo = SymmConf.Symmetry_Info();
		n_symm_subunit  = symminfo->score_multiply_factor();
	}

	utility::vector1< bool > rsd_mask( pose.size(), true);

	if ( get_selector() && use_selector_as_zscore_mask_ ) {
		rsd_mask = get_selector()->apply( pose );
	}

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue(i).aa() == core::chemical::aa_vrt || pose.residue_type(i).is_virtual_residue() ) continue;
		if ( symminfo && !symminfo->bb_is_independent( i ) ) continue; // only the main chain gets selected

		if ( get_selector() && (! rsd_mask[i]) ) continue;

		per_rsd_dens[i] = 0.0;
		per_rsd_nbrdens[i] = 0.0;
		per_rsd_rama[i] = 0.0;
		per_rsd_geometry[i] = 0.0;
	}


	// to catch some outliers
	calculate_density_nbr( pose, per_rsd_dens, per_rsd_nbrdens, symminfo, sliding_window_size_, mixed_sliding_window_ ); // per_rsd_dens_ and per_rsd_densnbr_

	//TR << "DENSITY" << per_rsd_dens.size() << " " <<per_rsd_nbrdens.size() << std::endl;
	for ( auto r : per_rsd_dens ) {
		std::cout << pose.pdb_info()->pose2pdb( r.first ) << " " << per_rsd_dens[r.first] <<" " << per_rsd_nbrdens[r.first] << std::endl;
	}


	calculate_rama( pose, per_rsd_rama, n_symm_subunit ); // per_rsd_rama_
	//std::cout << "RAMA" << per_rsd_rama.size() << std::endl;
	for ( auto r :per_rsd_rama ) {
		std::cout << pose.pdb_info()->pose2pdb( r.first ) << " " <<per_rsd_rama[ r.first ] << std::endl;
	}

	calculate_geometry( pose, per_rsd_geometry , n_symm_subunit); // per_rsd_geometry_
	//std::cout << "GEOMETRY" << per_rsd_geometry.size() << std::endl;
	for ( auto r: per_rsd_geometry ) {
		std::cout << pose.pdb_info()->pose2pdb( r.first ) << " " << per_rsd_geometry[ r.first ] << std::endl;
	}
	//utility_exit_with_message("Stopping");
}

std::map< core::Size, core::Real >
PerResidueDensityFitMetric::calculate(const pose::Pose & input_pose) const {

	using namespace core::scoring;
	using namespace numeric;

	std::map< core::Size, core::Real > result;

	pose::Pose pose = input_pose; //Needs to be scored


	std::map< Size, Real > per_rsd_dens, per_rsd_nbrdens, per_rsd_rama, per_rsd_geometry;
	std::map< Size, Real > zscore_dens, zscore_nbrdens, zscore_rama, zscore_geometry;
	std::map< Size, Real > scores;


	//Get Zscore.  Make Cut.  Flip if need be.
	core::conformation::symmetry::SymmetryInfoCOP symminfo = nullptr;
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		auto const & SymmConf (
			dynamic_cast<core::conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
		symminfo = SymmConf.Symmetry_Info();
	}

	utility::vector1< Size > rsd_mask = get_selector()->apply(pose);

	if ( match_res_ ) {
		TR << "Using basic correlation to density instead of Zscore." << std::endl;

		core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
		edm.setScoreWindowContext( true );
		edm.setWindow(1);  // smoother to use 3-res window

		core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
		myscore->set_weight( core::scoring::elec_dens_window, 1.0 );
		myscore->score(pose);

		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( ! rsd_mask[i] ) continue;
			Real dens_rscc = core::scoring::electron_density::getDensityMap().matchRes( i , pose.residue(i), pose, symminfo , false);

			TR.Debug << pose.pdb_info()->pose2pdb(i) << " " << dens_rscc << std::endl;
			result[i] = dens_rscc;
		}
		return result;
	}

	//ZScore Computation
	compute_scores(pose, per_rsd_dens, per_rsd_nbrdens, per_rsd_rama, per_rsd_geometry);

	calc_zscore( per_rsd_dens,     zscore_dens           );

	TR << "Calculating NBR dens " << std::endl;
	calc_zscore( per_rsd_nbrdens,  zscore_nbrdens        );
	calc_zscore( per_rsd_rama,     zscore_rama,     true );
	calc_zscore( per_rsd_geometry, zscore_geometry, true );

	for ( auto r : per_rsd_dens ) {
		if ( ! rsd_mask[ r.first ] ) continue;
		TR.Debug << per_rsd_dens[r.first] <<" " << per_rsd_nbrdens[r.first] << " " << per_rsd_rama[r.first] << " " << per_rsd_geometry[ r.first] << std::endl;
		TR.Debug << zscore_dens[r.first] << " " << zscore_nbrdens[r.first] << " "<< zscore_rama[r.first] << " " <<zscore_geometry[r.first] << std::endl;
		Real score =  0.45*zscore_dens[r.first]
			+ 0.05*zscore_nbrdens[r.first]
			+ 0.15*zscore_rama[r.first]
			+ 0.35*zscore_geometry[r.first];

		TR.Debug <<" " << pose.pdb_info()->pose2pdb( r.first) << " " << score << std::endl;
		result[r.first] = score;

		//fragbias_tr << "rsn: " << r << " fragProb: " << fragmentProbs_[r] << " score: " << score << std::endl;
	}

	//Correct for symmetry:
	if ( core::pose::symmetry::is_symmetric( pose )  ) {
		for ( Size i = 1; i <= pose.size(); ++i ) {
			if ( symminfo->bb_is_independent(i) || (!result.count( symminfo->bb_follows(i) ) ) ) {
				continue;
			} else if ( rsd_mask[i] ) {
				result[i] = result[ symminfo->bb_follows(i)];
			}
		}
	}
	return result;
}

void
PerResidueDensityFitMetricCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	PerResidueDensityFitMetric::provide_xml_schema( xsd );
}

std::string
PerResidueDensityFitMetricCreator::keyname() const {
	return PerResidueDensityFitMetric::name_static();
}

core::simple_metrics::SimpleMetricOP
PerResidueDensityFitMetricCreator::create_simple_metric() const {
	return core::simple_metrics::SimpleMetricOP( new PerResidueDensityFitMetric );

}



} //core
} //simple_metrics
} //per_residue_metrics


#ifdef    SERIALIZATION



template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueDensityFitMetric::save( Archive & arc ) const {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric>( this ) );
	arc( CEREAL_NVP( sliding_window_size_ ) );
	arc( CEREAL_NVP( mixed_sliding_window_));
	arc( CEREAL_NVP( match_res_));
	arc( CEREAL_NVP( use_selector_as_zscore_mask_));

}

template< class Archive >
void
core::simple_metrics::per_residue_metrics::PerResidueDensityFitMetric::load( Archive & arc ) {
	arc( cereal::base_class< core::simple_metrics::PerResidueRealMetric >( this ) );
	arc( sliding_window_size_ );
	arc( mixed_sliding_window_ );
	arc( match_res_);
	arc( use_selector_as_zscore_mask_);


}

SAVE_AND_LOAD_SERIALIZABLE( core::simple_metrics::per_residue_metrics::PerResidueDensityFitMetric );
CEREAL_REGISTER_TYPE( core::simple_metrics::per_residue_metrics::PerResidueDensityFitMetric )

CEREAL_REGISTER_DYNAMIC_INIT( core_simple_metrics_per_residue_metrics_PerResidueDensityFitMetric )
#endif // SERIALIZATION



