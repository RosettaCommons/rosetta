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

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/select/residue_selector/ResidueSelectorFactory.hh>
#include <core/select/residue_selector/util.hh> // for xml schema utility functions
#include <core/scoring/Energies.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/select/residue_selector/util.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>

// symmetry
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/electron_density/util.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>
#include <numeric/zscores.hh>

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
	using namespace core::scoring::electron_density;

/// @brief Constructor.
///
DensityFitResidueSelector::DensityFitResidueSelector():
	core::select::residue_selector::ResidueSelector()
{}

/// @brief Destructor.
///
DensityFitResidueSelector::~DensityFitResidueSelector() {}

/// @brief Copy Constructor.
//DensityFitResidueSelector::DensityFitResidueSelector(DensityFitResidueSelector const & src):
// core::select::residue_selector::ResidueSelector( src )
//{
//}

/// @brief Clone function.
/// @details Copy this object and return owning pointer to the copy (created on the heap).
core::select::residue_selector::ResidueSelectorOP
DensityFitResidueSelector::clone() const {
	return core::select::residue_selector::ResidueSelectorOP( new DensityFitResidueSelector(*this) );
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
	set_mixed_sliding_window( tag->getOption< bool >("mixed_sliding_window", mixed_sliding_window_) );
	set_sliding_window_size( tag->getOption< Size >("sliding_window_size", sliding_window_size_) );
	set_match_mode( tag->getOption< bool >("match_res", match_res_) );
	set_use_selector_as_mask( tag->getOption< bool >("use_selector_as_mask", use_selector_as_mask_));
	
	if ( tag->hasOption("residue_selector") ) {
		mask_ = parse_residue_selector( tag, datamap );
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
		+ XMLSchemaAttribute::attribute_w_default( "use_selector_as_mask",  xsct_rosetta_bool, "Use the selector as true mask to calculate the Zscore.  Otherwise, use it just as a selection for computation.  Default true.", "true");
	
	core::select::residue_selector::attributes_for_parse_residue_selector_default_option_name(attributes, "A Residue selector mask.  Used to only compute Zscore among a set of residues.  Useful for protein vs glycan density.  Since match_res is NOT a zscore, the selector acts as an AND selector, so we only compute the correlations on this set. " );
	std::string documentation = "Select residues that have a good electron density fit. (Or bad fit using the invert option)";
	
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
DensityFitResidueSelector::set_residue_mask( ResidueSelectorCOP selector ){
	mask_ = selector;
}

void
DensityFitResidueSelector::set_sliding_window_size(core::Size window_size){
	sliding_window_size_ = window_size;
}

void
DensityFitResidueSelector::set_mixed_sliding_window(bool mixed_sliding_window){
	mixed_sliding_window_ = mixed_sliding_window;
}

void
DensityFitResidueSelector::set_invert( bool invert ){
	invert_ = invert;
}

void
DensityFitResidueSelector::set_score_cut(Real score_cut){
	score_cut_ = score_cut;
}

void
DensityFitResidueSelector::set_match_mode( bool match_mode ){
	match_res_ = match_mode;
}

void
DensityFitResidueSelector::set_use_selector_as_mask(bool selector_as_mask){
	use_selector_as_mask_ = selector_as_mask;
}

void
DensityFitResidueSelector::compute_scores(
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
	
	if (mask_ && use_selector_as_mask_ ){
		rsd_mask = mask_->apply( pose );
	}
	
	for (core::Size i = 1; i <= pose.size(); ++i){
		if ( pose.residue(i).aa() == core::chemical::aa_vrt || pose.residue_type(i).is_virtual_residue() ) continue;
		if ( symminfo && !symminfo->bb_is_independent( i ) ) continue; // only the main chain gets selected
		
		if (mask_ && (! rsd_mask[i])) continue;
		
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
	for (auto r :per_rsd_rama){
		std::cout << pose.pdb_info()->pose2pdb( r.first ) << " " <<per_rsd_rama[ r.first ] << std::endl;
	}
	
	calculate_geometry( pose, per_rsd_geometry , n_symm_subunit); // per_rsd_geometry_
	//std::cout << "GEOMETRY" << per_rsd_geometry.size() << std::endl;
	for (auto r: per_rsd_geometry){
		std::cout << pose.pdb_info()->pose2pdb( r.first ) << " " << per_rsd_geometry[ r.first ] << std::endl;
	}
	//utility_exit_with_message("Stopping");
}

/// @brief "Apply" function.
/// @details Given the pose, generate a vector of bools with entries for every residue in the pose
/// indicating whether each residue is selected ("true") or not ("false").
DensityFitResidueSelector::ResidueSubset
DensityFitResidueSelector::apply( core::pose::Pose const & input_pose ) const
{
	using namespace core::scoring;
	using namespace numeric;
	
	pose::Pose pose = input_pose;
	
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
	
	utility::vector1< Size > subset( pose.size(), false);
	utility::vector1< Size > rsd_mask( pose.size(), false);
	
	if (mask_ && use_selector_as_mask_){
		rsd_mask = mask_->apply(pose);
	}
	if (match_res_){
		TR << "Using basic correlation to density instead of Zscore." << std::endl;
		
		core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();
		edm.setScoreWindowContext( true );
		edm.setWindow(1);  // smoother to use 3-res window
		
		core::scoring::ScoreFunctionOP myscore( new core::scoring::ScoreFunction() );
		myscore->set_weight( core::scoring::elec_dens_window, 1.0 );
		myscore->score(pose);
	
		//This catches both code-level and RS interface.
		if (score_cut_ == -.5){
			utility_exit_with_message("Using match_res with a the default cutoff.  This mode gives a correlation score to density and should be between 0 and 1");
		}
		
		for (Size i = 1; i <= pose.size(); ++i){
			if (! rsd_mask[i] ) continue;
			Real dens_rscc = core::scoring::electron_density::getDensityMap().matchRes( i , pose.residue(i), pose, symminfo , false);
			
			TR << pose.pdb_info()->pose2pdb(i) << " " << dens_rscc << std::endl;
			if ( (! invert_) && dens_rscc >= score_cut_ ){
				TR << "Good Fit" << std::endl;
				subset[ i ] = true;
			}
			else if ( invert_ && dens_rscc < score_cut_){
				TR <<"Bad Fit" << std::endl;
				subset[i] = true;
			}
		}
		return subset;
	}
	
	
	//ZScore Computation
	compute_scores(pose, per_rsd_dens, per_rsd_nbrdens, per_rsd_rama, per_rsd_geometry);
	
	calc_zscore( per_rsd_dens,     zscore_dens           );
	
	TR << "Calculating NBR dens " << std::endl;
	calc_zscore( per_rsd_nbrdens,  zscore_nbrdens        );
	calc_zscore( per_rsd_rama,     zscore_rama,     true );
	calc_zscore( per_rsd_geometry, zscore_geometry, true );
	
	
	if (mask_){
		rsd_mask = mask_->apply( pose );
	}
	for ( auto r : per_rsd_dens ) {
		if ( ! rsd_mask[ r.first ] ) continue;
		TR << per_rsd_dens[r.first] <<" " << per_rsd_nbrdens[r.first] << " " << per_rsd_rama[r.first] << " " << per_rsd_geometry[ r.first] << std::endl;
		TR << zscore_dens[r.first] << " " << zscore_nbrdens[r.first] << " "<< zscore_rama[r.first] << " " <<zscore_geometry[r.first] << std::endl;
		Real score =  0.45*zscore_dens[r.first]
			+ 0.05*zscore_nbrdens[r.first]
			+ 0.15*zscore_rama[r.first]
			+ 0.35*zscore_geometry[r.first];
	
		TR <<" " << pose.pdb_info()->pose2pdb( r.first) << " " << score << std::endl;
		if (invert_ && (score <= score_cut_)){
			TR << "Bad fit " << std::endl;
			subset[ r.first ] = true;
		}
		else if ( !invert_ && (score > score_cut_) ){
			TR << "Good fit " << std::endl;
			subset[ r.first]  = true;
		}

		//fragbias_tr << "rsn: " << r << " fragProb: " << fragmentProbs_[r] << " score: " << score << std::endl;
	}
	
	//Correct for symmetry:
	if ( core::pose::symmetry::is_symmetric( pose )  ){
		for (Size i = 1; i <= pose.size(); ++i){
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
	arc( CEREAL_NVP( sliding_window_size_ ) );
	arc( CEREAL_NVP( mixed_sliding_window_) );
	arc( CEREAL_NVP( mask_ ) );
	arc( CEREAL_NVP( match_res_ ) );
	arc( CEREAL_NVP( use_selector_as_mask_ ) );
}

template< class Archive >
void
core::select::residue_selector::DensityFitResidueSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( score_cut_ );
	arc( invert_ );
	arc( sliding_window_size_ );
	arc( mixed_sliding_window_ );
	arc( mask_ );
	arc( match_res_ );
	arc( use_selector_as_mask_ );
	
}

SAVE_AND_LOAD_SERIALIZABLE( core::select::residue_selector::DensityFitResidueSelector );
CEREAL_REGISTER_TYPE( core::select::residue_selector::DensityFitResidueSelector )

CEREAL_REGISTER_DYNAMIC_INIT( core_select_residue_selector_DensityFitResidueSelector )
#endif // SERIALIZATION
