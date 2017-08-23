// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/carbohydrates/LinkageConformerMover.cc
/// @brief This code changes all of the dihedrals of a particular glycosidic linkage based on database info,
///   esentially sampling carbohydrate dihedral conformers of two residues.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)
/// @author Labonte <JWLabonte@jhu.edu>

// Unit headers
#include <protocols/carbohydrates/LinkageConformerMover.hh>
#include <protocols/carbohydrates/LinkageConformerMoverCreator.hh>

// Package headers
#include <core/chemical/carbohydrates/CarbohydrateInfoManager.hh>
#include <core/chemical/carbohydrates/CarbohydrateInfo.hh>

// Project headers
#include <core/id/types.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/util.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/GlycanResidueSelector.hh>
#include <core/select/residue_selector/ReturnResidueSubsetSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>

#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>
#include <protocols/rosetta_scripts/util.hh>


// Numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/WeightedSampler.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <basic/basic.hh>

// Utility header
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.carbohydrates.LinkageConformerMover" );


namespace protocols {
namespace carbohydrates {
using namespace core::kinematics;
using namespace protocols::simple_moves;
using namespace protocols::simple_moves::bb_sampler;

LinkageConformerMover::LinkageConformerMover():
	protocols::moves::Mover( "LinkageConformerMover" ),
	phi_sampler_mover_(/* Null */),
	psi_sampler_mover_(/* Null */)
{
	set_defaults();
}

LinkageConformerMover::LinkageConformerMover( core::select::residue_selector::ResidueSelectorCOP selector ):
	protocols::moves::Mover( "LinkageConformerMover" ),
	phi_sampler_mover_(/* Null */),
	psi_sampler_mover_(/* Null */)
{
	set_defaults();
	set_residue_selector( selector );
}

void
LinkageConformerMover::set_defaults(){

	sample_sd_ = 1.0;
	use_sugar_bb_data_if_needed_ = true;
	idealize_torsions_ = false;
	conformer_found_ = false;
	use_sd_as_prob_ = false;
	sample_protein_linkage_ = true;
	use_conformer_population_stats_ = true;

	SugarBBSamplerOP phi_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::phi_dihedral ) );
	phi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( phi_sampler ) );

	SugarBBSamplerOP psi_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::psi_dihedral ) );
	psi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( psi_sampler ) );

	SugarBBSamplerOP omega_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::omega_dihedral ) );
	omega_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( omega_sampler ) );

}

LinkageConformerMover::~LinkageConformerMover()= default;

LinkageConformerMover::LinkageConformerMover( LinkageConformerMover const & src ):
	protocols::moves::Mover(src),
	sample_sd_(src.sample_sd_),
	use_sugar_bb_data_if_needed_(src.use_sugar_bb_data_if_needed_),
	idealize_torsions_(src.idealize_torsions_),
	conformer_found_(src.conformer_found_),
	use_sd_as_prob_(src.use_sd_as_prob_),
	sample_protein_linkage_(src.sample_protein_linkage_),
	use_conformer_population_stats_(src.use_conformer_population_stats_),
	random_sampler_(src.random_sampler_)
{
	if ( src.selector_ ) selector_ = src.selector_->clone();
	phi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( *src.phi_sampler_mover_ ) );
	psi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( *src.psi_sampler_mover_ ) );
	omega_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( *src.omega_sampler_mover_ ) );
}

void
LinkageConformerMover::set_residue_selector( core::select::residue_selector::ResidueSelectorCOP selector ){

	selector_ = selector->clone();

}

void
LinkageConformerMover::set_single_resnum( core::pose::Pose const & pose, core::Size resnum ){
	using namespace core::select::residue_selector;

	utility::vector1< bool > subset(pose.total_residue(), false);

	subset[ resnum ] = true;
	selector_ = ReturnResidueSubsetSelectorOP( new ReturnResidueSubsetSelector( subset));
}

void
LinkageConformerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& datamap,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	if ( tag->hasOption("residue_selector") ) {
		selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, datamap );
	}

	sample_sd_ = tag->getOption< core::Real >("x_sds", sample_sd_);
	use_sugar_bb_data_if_needed_ = tag->getOption< bool >( "use_sugar_bb_if_needed", use_sugar_bb_data_if_needed_);
	idealize_torsions_ = tag->getOption< bool >("idealize_torsions", idealize_torsions_);
	use_sd_as_prob_ = tag->getOption< bool >("prob_sd_sampling", use_sd_as_prob_);
	sample_protein_linkage_ = tag->getOption< bool >("sample_protein_linkage", sample_protein_linkage_);
	use_conformer_population_stats_ =
		tag->getOption< bool >( "use_conformer_population_stats", use_conformer_population_stats_ );
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

protocols::moves::MoverOP
LinkageConformerMover::clone() const{
	return protocols::moves::MoverOP( new LinkageConformerMover( *this ) );
}

/*
LinkageConformerMover & LinkageConformerMoveroperator=( LinkageConformerMover const & src){
return LinkageConformerMover( src );
}
*/


moves::MoverOP
LinkageConformerMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new LinkageConformerMover );
}

// XRW TEMP std::string
// XRW TEMP LinkageConformerMover::get_name() const {
// XRW TEMP  return "LinkageConformerMover";
// XRW TEMP }

void
LinkageConformerMover::show(std::ostream & output) const
{
	protocols::moves::Mover::show(output);
}

std::ostream &operator<< (std::ostream &os, LinkageConformerMover const &mover)
{
	mover.show(os);
	return os;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
LinkageConformerMover::set_idealize_torsions(bool idealize_torsions) {
	idealize_torsions_ = idealize_torsions;
}

void
LinkageConformerMover::set_x_standard_deviations(core::Real standard_deviation){
	sample_sd_ = standard_deviation;
}

void
LinkageConformerMover::set_use_sugar_bb_data_if_needed(bool use_sugar_bb){
	use_sugar_bb_data_if_needed_ = use_sugar_bb;
}

bool
LinkageConformerMover::conformer_found() const {
	return conformer_found_;
}

void
LinkageConformerMover::set_prob_sd_sampling( bool prob_sampling ){
	use_sd_as_prob_ = prob_sampling;
}

void
LinkageConformerMover::apply( core::pose::Pose & pose )
{
	using namespace core::pose::carbohydrates;
	using namespace core::chemical::carbohydrates;
	using namespace core::chemical;
	using namespace core::select::residue_selector;

	reset_status();

	if ( ! selector_ ) {
		TR << "No Residue Selector set.  Attempting to use all carbohydrate residues." << std::endl;
		selector_ = GlycanResidueSelectorOP( new GlycanResidueSelector() );
	}

	phi_sampler_mover_->set_residue_selector(selector_);
	psi_sampler_mover_->set_residue_selector(selector_);
	omega_sampler_mover_->set_residue_selector( selector_ );

	utility::vector1< bool > subset = selector_->apply(pose);
	utility::vector1< core::Size > movemap_residues = selection_positions( subset );
	conformer_found_ = false;


	core::Size index = numeric::random::rg().random_range2( 1, movemap_residues.size() );
	core::Size upper_resnum = movemap_residues[ index ];



	core::chemical::ResidueType const & res = pose.residue_type( upper_resnum );
	if ( ! res.is_carbohydrate() ) {
		TR << "Selected residue does not have a glycosidic linkage to its parent!  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}

	core::Size lower_resnum = pose.glycan_tree_set()->get_parent( upper_resnum );
	if ( lower_resnum == 0 ) {
		TR << "Selected residue has no parent.  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}
	core::chemical::ResidueType const & parent_res = pose.residue_type(  lower_resnum );



	std::string res_name = res.carbohydrate_info()->short_name();
	std::string parent_res_name;
	if ( parent_res.is_carbohydrate() ) {
		parent_res_name = parent_res.carbohydrate_info()->short_name();  // 3-letter code not enough

		//This is due to multiple connecting points possible. This gets the position of linking to the previous residue!
		core::Size const link_pos = pose.glycan_tree_set()->get_linkage_position( upper_resnum );

		parent_res_name[ 2 ] = '0' + link_pos;  // Set the correct connectivity.

	} else if ( sample_protein_linkage_ ) {
		parent_res_name = parent_res.name3();
	} else {
		TR << "Sampling of linkages to protein has been disabled.  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}

	set_last_move_status(protocols::moves::MS_FAIL); //Fail unless we say we are ok.



	//TR << "Upper resnum: " << upper_resnum << "  Lower resnum: " << lower_resnum << std::endl;
	TR << "Sampling " << res_name << "(?" << parent_res_name << " linkage " << std::endl;
	if ( CarbohydrateInfoManager::pair_has_linkage_statistics( parent_res_name, res_name ) ) {
		utility::vector1< LinkageConformerData > conformers =
			CarbohydrateInfoManager::linkages_from_pair( parent_res_name, res_name );
		core::Size const n_conformers = conformers.size() ;
		core::Size conformer_num;
		if ( n_conformers == 1 ) {
			conformer_num = 1;
		} else if ( use_conformer_population_stats_ ) {
			//TR.Debug << "Sampling on conformer stats" << std::endl;
			utility::vector1< core::Real > populations;
			for ( core::uint i =1 ; i <= n_conformers; ++i ) {
				populations.push_back( conformers[ i ].population );
			}
			numeric::random::WeightedSampler sampler( populations );
			conformer_num = sampler.random_sample( numeric::random::rg() );
		} else {
			//TR.Debug << "Sampling all conformers: Total="<<conformers.size() << std::endl;
			conformer_num = numeric::random::rg().random_range2( 1, conformers.size() );
		}
		TR << "Sampling conformer " << conformer_num << " which has a population of " << conformers[conformer_num].population << std::endl;

		set_dihedrals_from_linkage_conformer_data(
			pose, upper_resnum, conformers[ conformer_num ], idealize_torsions_, use_sd_as_prob_ );

		TR.Info << "Complete" << std::endl;
		conformer_found_ = true;
		set_last_move_status(protocols::moves::MS_SUCCESS);

	} else if ( use_sugar_bb_data_if_needed_ ) {

		TR << "no conformer found.  Using sugar BB sampler if possible." << std::endl;
		if ( ! parent_res.is_carbohydrate() ) {
			TR <<" Cannot run sugar BB Sampler on linkage as linkage as a protein-glycan linkage. Skipping" << std::endl;
			set_last_move_status(protocols::moves::MS_FAIL);
		}
		phi_sampler_mover_->set_single_resnum( upper_resnum );
		psi_sampler_mover_->set_single_resnum( upper_resnum );
		omega_sampler_mover_->set_single_resnum( upper_resnum );

		phi_sampler_mover_->apply(pose);
		psi_sampler_mover_->apply(pose);
		omega_sampler_mover_->apply(pose);



	} else {
		TR << "No conformer data found and use_sugar_bb_data FALSE.  Doing nothing. " << std::endl;
		set_last_move_status(protocols::moves::MS_FAIL);
	}
}


/////////////// Creator ///////////////

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP LinkageConformerMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new LinkageConformerMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LinkageConformerMoverCreator::keyname() const {
// XRW TEMP  return LinkageConformerMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP LinkageConformerMover::mover_name(){
// XRW TEMP  return "LinkageConformerMover";
// XRW TEMP }

std::string LinkageConformerMover::get_name() const {
	return mover_name();
}

std::string LinkageConformerMover::mover_name() {
	return "LinkageConformerMover";
}

void LinkageConformerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute::required_attribute("upper_resnum", xsct_non_negative_integer, "XRW TO DO")
		+ XMLSchemaAttribute("x_sds", xsct_real, "Standard deviation for sampling")
		+ XMLSchemaAttribute("use_sugar_bb_if_needed", xsct_rosetta_bool, "Use sugar backbone data if needed?")
		+ XMLSchemaAttribute("idealize_torsions", xsct_rosetta_bool, "Idealize torsion angles before run?")
		+ XMLSchemaAttribute("prob_sd_sampling", xsct_rosetta_bool, "Use standard deviation as probability")
		+ XMLSchemaAttribute("sample_protein_linkage", xsct_rosetta_bool, "Also sample linkage between glycan and protein")
		+ XMLSchemaAttribute("use_conformer_population_stats", xsct_rosetta_bool, "Use statistics about conformer populations for sampling");

	rosetta_scripts::attributes_for_parse_residue_selector( attlist );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Mover to sample glycan linkages", attlist );
}

std::string LinkageConformerMoverCreator::keyname() const {
	return LinkageConformerMover::mover_name();
}

protocols::moves::MoverOP
LinkageConformerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LinkageConformerMover );
}

void LinkageConformerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LinkageConformerMover::provide_xml_schema( xsd );
}


} //carbohydrates
} //protocols
