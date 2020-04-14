// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
#include <core/pose/chains_util.hh>
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/ResfileReader.cc
/// @brief  implementation of resfile reader and its command classes
/// @author Gordon Lemmon (glemmon@gmail.com)

// Unit Headers
#include <protocols/ligand_docking/Translate.hh>
#include <protocols/ligand_docking/TranslateCreator.hh>

#include <protocols/ligand_docking/RandomConformerMover.hh>
#include <protocols/ligand_docking/grid_functions.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/scoring_grid/GridSet.hh>
#include <protocols/qsar/scoring_grid/schema_util.hh>

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.io.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

namespace protocols {
namespace ligand_docking {

static basic::Tracer translate_tracer( "protocols.ligand_docking.ligand_options.translate", basic::t_debug );

void
Translate_info::set_chain_id( core::Size id ) {
	by_string_ = false;
	chain_number_ = id;
}

void
Translate_info::set_chain_letter( std::string const & str) {
	by_string_ = true;
	chain_string_ = str;
}

core::Size
Translate_info::chain_id( core::pose::Pose const & pose ) const {
	if ( ! by_string_ ) {
		return chain_number_;
	}
	utility::vector1< core::Size > chain_ids( core::pose::get_chain_ids_from_chain(chain_string_, pose) );
	if ( chain_ids.size() == 0 ) {
		utility_exit_with_message( "'Translate' mover: chain '"+chain_string_+"' does not exist.");
	} else if ( chain_ids.size() > 1 ) {
		utility_exit_with_message( "'Translate' mover: chain letter '"+chain_string_+"' represents more than one chain. Consider the 'Translates' mover (with an 's') instead.");
	}
	return chain_ids[1];
}

char
Translate_info::chain_letter( core::pose::Pose const & pose ) const {
	if ( by_string_ ) {
		if ( chain_string_.size() != 1 ) {
			utility_exit_with_message("'Translate' mover: chain designation '"+chain_string_+"' is not one character.");
		}
		return chain_string_[1];
	}

	return core::pose::get_chain_from_chain_id( chain_number_, pose );
}

core::Size
Translate_info::jump_id( core::pose::Pose const & pose ) const {
	// This indirection existed in the previous version of the mover - we probably want to short-circuit it.
	return core::pose::get_jump_id_from_chain_id(chain_id(pose), pose);
}


/// @brief
Translate::Translate():
	Mover("Translate"),
	translate_info_()
{}

Translate::Translate(Translate_info translate_info):
	Mover("Translate"),
	translate_info_(translate_info)
{}

protocols::moves::MoverOP Translate::clone() const {
	return utility::pointer::make_shared< Translate >( *this );
}

protocols::moves::MoverOP Translate::fresh_instance() const {
	return utility::pointer::make_shared< Translate >();
}


/// @brief parse XML (specifically in the context of the parser/scripting scheme)
void
Translate::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map
)
{
	if ( ! tag->hasOption("chain") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover requires chain tag");
	if ( ! tag->hasOption("distribution") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover requires distribution tag");
	if ( ! tag->hasOption("angstroms") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover requires angstroms tag");
	if ( ! tag->hasOption("cycles") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover requires cycles tag");
	//if ( ! tag->hasOption("force") ) throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover requires force tag"); optional. default is don't force, meaning ligand stays put if it can't find somewhere to go.

	// Will return a nullptr if this XML didn't define any ScoringGrids
	grid_set_prototype_ = protocols::qsar::scoring_grid::parse_optional_grid_set_from_tag( tag, data_map );

	std::string chain = tag->getOption<std::string>("chain");
	if ( chain.size() != 1 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'Translate' mover: chain must be a single letter.");
	}
	translate_info_.set_chain_letter( chain );

	std::string distribution_str= tag->getOption<std::string>("distribution");
	translate_info_.distribution= get_distribution(distribution_str);
	translate_info_.angstroms = tag->getOption<core::Real>("angstroms");
	translate_info_.cycles = tag->getOption<core::Size>("cycles");

	if ( tag->hasOption("force") ) {
		if ( tag->getOption<std::string>("force") == "true" ) {
			translate_info_.force= true;
		} else if ( tag->getOption<std::string>("force") != "false" ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "'force' option is true or false");
		}
	}

	if ( tag->hasOption("tag_along_chains") ) {
		std::string const tag_along_chains_str = tag->getOption<std::string>("tag_along_chains");
		tag_along_chains_ = utility::string_split(tag_along_chains_str, ',');
	}

}

void Translate::apply(core::pose::Pose & pose) {
	// Put chain/jump determination at the top, to quickly error out if the pose isn't suitable.
	core::Size chain_id = translate_info_.chain_id(pose);
	core::Size jump_id = translate_info_.jump_id(pose);
	char chain_letter = translate_info_.chain_letter(pose);

	core::Size const begin(pose.conformation().chain_begin(chain_id));

	utility::vector1<core::Size> chain_ids_to_exclude( chain_ids_to_exclude_ ); // these are invisible the translation grid, so ligand can land on top.
	// add this Translate's chain conditionally (for use with CompoundTranslate)
	if ( ! chain_ids_to_exclude.contains( chain_id ) ) {
		chain_ids_to_exclude.push_back(chain_id);
	}
	for ( std::string const & tag_along_chain_str: tag_along_chains_ ) {
		for ( core::Size const ta_chain_id : get_chain_ids_from_chain(tag_along_chain_str, pose) ) {
			chain_ids_to_exclude.push_back(ta_chain_id);
		}
	}

	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);

	if ( grid_set_prototype_ == nullptr ) {
		utility::pointer::shared_ptr<core::grid::CartGrid<int> > const grid = make_atr_rep_grid_without_ligands(pose, center, chain_ids_to_exclude);
		translate_ligand(grid, jump_id, pose);// move ligand to a random point in binding pocket
	} else {
		// How we treat chains here is rather poor, but it should work in the historical one-residue-ligand case.
		qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
		qsar::scoring_grid::GridSetCOP grid_set( grid_manager->get_grids( *grid_set_prototype_, pose, center, chain_letter ) );
		//TODO refactor qsar map so it works properly
		/*
		if(!grid_manager->is_qsar_map_attached())
		{
		utility::vector1<std::string> grid_names( grid_manager->get_grid_names() );
		core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
		qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
		qsar_map->fill_with_value(1,grid_names);
		grid_set->set_qsar_map(qsar_map);
		}
		*/
		translate_ligand(grid_set, jump_id, pose, begin);
	}
}

void Translate::uniform_translate_ligand(
	const utility::pointer::shared_ptr<core::grid::CartGrid<int> >  & grid,
	const core::Size jump_id,
	core::pose::Pose & pose
){
	translate_tracer.Debug<< "making a uniform translator of up to " << translate_info_.angstroms<< " angstroms" << std::endl;
	protocols::rigid::UniformSphereTransMover mover( jump_id, translate_info_.angstroms);

	core::pose::Pose const orig_pose(pose);
	core::pose::Pose best_pose;
	int best_score=100000;

	std::vector<core::Size> tag_along_jumps;
	for ( std::string const & chain: tag_along_chains_ ) {
		// A bit indirect, but matches what was used historically.
		for ( core::Size chain_id: get_chain_ids_from_chain(chain, pose) ) {
			tag_along_jumps.push_back( core::pose::get_jump_id_from_chain_id(chain_id, pose) );
		}
	}

	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	for ( core::Size cycle = 0; cycle < translate_info_.cycles; ++cycle ) {
		mover.apply(pose);
		core::Vector c = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
		// did our nbr_atom land in an empty space on the grid?
		// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
		int grid_value=best_score;
		if ( grid->is_in_grid(c.x(), c.y(), c.z()) ) {
			grid_value= grid->getValue(c.x(), c.y(), c.z());
		}
		if ( grid_value <= 0 ) {
			translate_tracer.Trace << "Accepting ligand position with nbr_atom at " << c << std::endl;
			mover.freeze();
			for ( core::Size const tag_along_jump : tag_along_jumps ) {
				translate_tracer.Trace << "moving jump " << tag_along_jump<< " the same amount"<< std::endl;
				mover.rb_jump(tag_along_jump);
				mover.apply(pose);
			}
			return;
		}
		translate_tracer.Trace << "Rejecting ligand position with nbr_atom at " << c << std::endl;

		if ( translate_info_.force && grid_value < best_score ) {
			best_score = grid_value;
			best_pose = pose;
		}
		pose = orig_pose; // reset and try again
	}
	// The code below only executes if we still haven't found a nonclashing position for the neighbor atom
	if ( translate_info_.force ) {
		translate_tracer.Trace << "Forcing neighbor atom to move (leading to a clash)" << std::endl;
		pose= best_pose;
		mover.freeze();
		for ( core::Size const tag_along_jump : tag_along_jumps ) {
			mover.rb_jump(tag_along_jump);
			mover.apply(pose);
		}
		return;
	} else {
		translate_tracer.Warning << "cannot find placement for this ligand.  Keeping original position. Use the force option to force translation"<< std::endl;
	}

}

void Translate::gaussian_translate_ligand(
	const utility::pointer::shared_ptr<core::grid::CartGrid<int> >  & grid,
	const core::Size jump_id,
	core::pose::Pose & pose
){
	translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
	protocols::rigid::RigidBodyPerturbMover mover( jump_id, 0 /*rotation*/, translate_info_.angstroms);

	core::pose::Pose const orig_pose(pose);
	core::pose::Pose best_pose;
	int best_score=100000;

	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	for ( core::Size cycle = 0; cycle < translate_info_.cycles; ++cycle ) {
		mover.apply(pose);
		core::Vector c = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
		// did our nbr_atom land in an empty space on the grid?
		// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
		int grid_value=0;
		if ( grid->is_in_grid(c.x(), c.y(), c.z()) ) {
			grid_value= grid->getValue(c.x(), c.y(), c.z()) <= 0;
		}
		if ( grid_value < 0 ) {
			translate_tracer.Trace << "Accepting ligand position with nbr_atom at " << c << std::endl;
			return;
		}
		translate_tracer.Trace << "Rejecting ligand position with nbr_atom at " << c << std::endl;

		if ( translate_info_.force && grid_value < best_score ) {
			best_score = grid_value;
			best_pose = pose;
		}
		pose = orig_pose; // reset and try again
	}
	if ( translate_info_.force ) {
		pose= best_pose;
		return;
	} else {
		translate_tracer.Warning << "cannot find placement for this ligand.  Keeping original position. Use the force option to force translation"<< std::endl;
	}
}

void Translate::translate_ligand(
	const utility::pointer::shared_ptr<core::grid::CartGrid<int> >  & grid,
	const core::Size jump_id,
	core::pose::Pose & pose
) {
	if ( translate_info_.angstroms < 0 ) utility_exit_with_message("cannot have a negative translation value");
	if ( translate_info_.angstroms == 0 ) return;

	if ( translate_info_.distribution == Uniform ) uniform_translate_ligand(grid, jump_id, pose);
	else if ( translate_info_.distribution == Gaussian ) gaussian_translate_ligand(grid, jump_id, pose);
}

void Translate::translate_ligand(qsar::scoring_grid::GridSetCOP grid_set, core::Size const jump_id,core::pose::Pose & pose, core::Size const & residue_id)
{
	if ( translate_info_.angstroms < 0 ) utility_exit_with_message("cannot have a negative translation value");
	if ( translate_info_.angstroms == 0 ) return;

	core::conformation::Residue const & residue(pose.residue(residue_id));

	protocols::rigid::RigidBodyMoverOP translate_mover;
	if ( translate_info_.distribution == Uniform ) {
		translate_tracer.Debug<< "making a uniform translator of up to "<< translate_info_.angstroms<<" angstroms"<< std::endl;
		translate_mover = utility::pointer::make_shared< protocols::rigid::UniformSphereTransMover >( jump_id, translate_info_.angstroms);
	} else if ( translate_info_.distribution == Gaussian ) {
		translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
		translate_mover = utility::pointer::make_shared< protocols::rigid::RigidBodyPerturbMover > ( jump_id, 0 /*rotation*/, translate_info_.angstroms);
	}

	RandomConformerMoverOP conformer_mover( new RandomConformerMover(residue_id) );

	core::pose::Pose  orig_pose(pose);
	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	core::Real best_score = 1000000;

	for ( core::Size cycle = 0; cycle < translate_info_.cycles; ++cycle ) {

		translate_tracer.Trace <<"Doing a translate move" <<std::endl;
		translate_mover->apply(pose);
		// }else
		// {
		//
		//  translate_tracer.Trace <<"Doing a conformer selection move" <<std::endl;
		//  conformer_mover->apply(pose);
		// }
		core::Real score(grid_set->total_score(residue));
		if ( score <= best_score ) {
			best_score = score;
			orig_pose = pose;
			translate_tracer.Trace << "accepted ligand position with score of " <<score <<std::endl;
		} else {
			pose = orig_pose;
			translate_tracer.Trace <<"rejected ligand position with score of "<<score <<std::endl;
		}
	}

}

core::Size
Translate::get_chain_id(core::pose::Pose const & pose ){
	return translate_info_.chain_id(pose);
}

void
Translate::add_excluded_chains(
	std::set<core::Size>::const_iterator begin,
	std::set<core::Size>::const_iterator end
){
	chain_ids_to_exclude_.insert( chain_ids_to_exclude_.end(), begin, end);
}

std::string Translate::get_name() const {
	return mover_name();
}

std::string Translate::mover_name() {
	return "Translate";
}

void Translate::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	// Restriction for sampling distribution tag
	XMLSchemaRestriction restriction_type;
	restriction_type.name( "distribution_string" );
	restriction_type.base_type( xs_string );
	restriction_type.add_restriction( xsr_pattern, "uniform|gaussian" );
	xsd.add_top_level_element( restriction_type );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute("chain", xs_string, "Chain ID of chain to be translated.")
		+ XMLSchemaAttribute::required_attribute("distribution", "distribution_string",
		"The random move can be chosen from a \"uniform\" or \"gaussian\" distribution.")
		+ XMLSchemaAttribute("force", xs_string, "XRW TO DO")
		+ XMLSchemaAttribute("tag_along_chains", xs_string, "XRW TO DO . Comma separated list of chain IDs to be moved together with \"chain\". ")
		+ XMLSchemaAttribute::required_attribute("angstroms", xsct_real, "Movement can be anywhere within a sphere of radius specified by \"angstroms\".")
		+ XMLSchemaAttribute::required_attribute("cycles", xsct_non_negative_integer,
		"Number of attempts to make such a movement without landing on top of another molecule.");

	protocols::qsar::scoring_grid::attributes_for_parse_optional_grid_set_from_tag( attlist, "The ScoringGrid set to use for scoring the translation. "
		"If no scoring grids (at all) are present in the XML, use a default classic grid." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Performs a coarse random movement of a small molecule in xyz-space.", attlist );
}

std::string TranslateCreator::keyname() const {
	return Translate::mover_name();
}

protocols::moves::MoverOP
TranslateCreator::create_mover() const {
	return utility::pointer::make_shared< Translate >();
}

void TranslateCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	Translate::provide_xml_schema( xsd );
}



} //namespace ligand_docking
} //namespace protocols
