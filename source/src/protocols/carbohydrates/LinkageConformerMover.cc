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

#include <protocols/moves/MoverStatus.hh>
#include <protocols/simple_moves/BBDihedralSamplerMover.hh>
#include <protocols/simple_moves/bb_sampler/SugarBBSampler.hh>

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
	psi_sampler_mover_(/* Null */),
	movemap_(/* Null */)
{
	set_defaults();
}

LinkageConformerMover::LinkageConformerMover( core::kinematics::MoveMapCOP movemap ):
	protocols::moves::Mover( "LinkageConformerMover" ),
	phi_sampler_mover_(/* Null */),
	psi_sampler_mover_(/* Null */),
	movemap_(/* Null */)
{
	set_defaults();
	set_movemap( movemap );
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
	movemap_residues_.clear();

	SugarBBSamplerOP phi_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::phi_dihedral ) );
	phi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( phi_sampler ) );

	SugarBBSamplerOP psi_sampler = SugarBBSamplerOP( new SugarBBSampler( core::id::psi_dihedral ) );
	psi_sampler_mover_ = BBDihedralSamplerMoverOP( new BBDihedralSamplerMover( psi_sampler ) );


}

LinkageConformerMover::~LinkageConformerMover()= default;

LinkageConformerMover::LinkageConformerMover( LinkageConformerMover const & )= default;

void
LinkageConformerMover::set_movemap( core::kinematics::MoveMapCOP movemap ){
	using namespace core::kinematics;
	movemap_ = movemap->clone();

}

void
LinkageConformerMover::set_single_resnum( core::Size resnum ){
	movemap_residues_.clear();
	movemap_residues_.push_back( resnum );
}

void
LinkageConformerMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap& ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & pose)
{
	if ( tag->hasOption("upper_resnum") ) {
		movemap_residues_.push_back( core::pose::parse_resnum( tag->getOption< std::string >("upper_resnum"), pose ) );
	} else {
		utility_exit_with_message("Must pass upper_resnum option for LinkageConformerMover for now");
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

	reset_status();
	movemap_residues_.clear();

	if ( ! movemap_ ) {
		TR << "No Movemap Set.  Attempting to use all carbohydrate residues." << std::endl;
		movemap_ = core::kinematics::MoveMapOP( new MoveMap);
		for ( core::Size i = 1; i <= pose.size(); ++i ) {
			if ( pose.residue_type(i).is_carbohydrate() ) {
				movemap_->set_bb( i, true);
			} else {
				movemap_->set_bb( i, false);
			}
		}
	}

	phi_sampler_mover_->set_movemap(movemap_);
	psi_sampler_mover_->set_movemap(movemap_);

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( pose.residue_type( i ).is_carbohydrate() && movemap_->get_bb(i) ) {
			movemap_residues_.push_back( i );
		}
	}

	conformer_found_ = false;


	core::Size index = numeric::random::rg().random_range( 1, movemap_residues_.size() );
	core::Size upper_resnum = movemap_residues_[ index ];



	core::chemical::ResidueType const & res2 = pose.residue_type( upper_resnum );
	if ( ! res2.is_carbohydrate() ) {
		TR << "Selected residue does not have a glycosidic linkage to its parent!  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}

	core::Size lower_resnum = find_seqpos_of_saccharides_parent_residue( pose.residue( upper_resnum ) );
	if ( lower_resnum == 0 ) {
		TR << "Selected residue has no parent.  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}
	core::chemical::ResidueType const & res1 = pose.residue_type(  lower_resnum );



	std::string res2_name = res2.carbohydrate_info()->short_name();
	std::string res1_name;
	if ( res1.is_carbohydrate() ) {
		res1_name = res1.carbohydrate_info()->short_name();  // 3-letter code not enough

		//This is due to multiple connecting points possible. This gets the position of linking to the previous residue!
		core::Size const link_pos = get_linkage_position_of_saccharide_residue( pose, upper_resnum );

		res1_name[ 2 ] = '0' + link_pos;  // Set the correct connectivity.

	} else if ( sample_protein_linkage_ ) {
		res1_name = res1.name3();
	} else {
		TR << "Sampling of linkages to protein has been disabled.  Skipping..." << std::endl;
		set_last_move_status( protocols::moves::MS_FAIL );
		return;
	}

	set_last_move_status(protocols::moves::MS_FAIL); //Fail unless we say we are ok.



	//TR << "Upper resnum: " << upper_resnum << "  Lower resnum: " << lower_resnum << std::endl;
	TR << "Sampling " << res2_name << "(?" << res1_name << " linkage " << std::endl;
	if ( CarbohydrateInfoManager::pair_has_linkage_statistics( res1_name, res2_name ) ) {
		utility::vector1< LinkageConformerData > conformers =
			CarbohydrateInfoManager::linkages_from_pair( res1_name, res2_name );
		core::Size const n_conformers = conformers.size() ;
		core::Size conformer_num;
		if ( n_conformers == 1 ) {
			conformer_num = 1;
		} else if ( use_conformer_population_stats_ ) {
			utility::vector1< core::Real > populations;
			for ( core::uint i =1 ; i <= n_conformers; ++i ) {
				populations.push_back( conformers[ i ].population );
			}
			numeric::random::WeightedSampler sampler( populations );
			conformer_num = sampler.random_sample( numeric::random::rg() );
		} else {
			conformer_num = numeric::random::rg().random_range( 1, conformers.size() );
		}
		set_dihedrals_from_linkage_conformer_data(
			pose, upper_resnum, conformers[ conformer_num ], idealize_torsions_, use_sd_as_prob_ );

		conformer_found_ = true;
		set_last_move_status(protocols::moves::MS_SUCCESS);

	} else if ( use_sugar_bb_data_if_needed_ ) {

		TR << "no conformer found.  Using sugar BB on phi and [ psi ]" << std::endl;
		phi_sampler_mover_->set_single_resnum( upper_resnum );
		psi_sampler_mover_->set_single_resnum( upper_resnum );

		phi_sampler_mover_->apply(pose);

		//Remove this when needed!
		if ( ! has_exocyclic_glycosidic_linkage( pose, upper_resnum  ) ) {
			psi_sampler_mover_->apply(pose);
		} else {
			TR << upper_resnum << " has glycosidic linkage.  Skipping psi sampling." << std::endl;
		}

		set_last_move_status(protocols::moves::MS_SUCCESS);
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
