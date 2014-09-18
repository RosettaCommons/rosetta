// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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

#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>

#include <core/pose/util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.io.hh>

#include <boost/foreach.hpp>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pose/Pose.hh>
using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static thread_local basic::Tracer translate_tracer( "protocols.ligand_docking.ligand_options.translate", basic::t_debug );

std::string
TranslateCreator::keyname() const
{
	return TranslateCreator::mover_name();
}

protocols::moves::MoverOP
TranslateCreator::create_mover() const {
	return new Translate;
}

std::string
TranslateCreator::mover_name()
{
	return "Translate";
}

///@brief
Translate::Translate():
		//utility::pointer::ReferenceCount(),
		Mover("Translate"),
		translate_info_(),
		chain_ids_to_exclude_(),
		tag_along_jumps_()
{}

Translate::Translate(Translate_info translate_info):
		//utility::pointer::ReferenceCount(),
		Mover("Translate"),
		translate_info_(translate_info),
		chain_ids_to_exclude_(),
		tag_along_jumps_()
{}

Translate::Translate(Translate const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		translate_info_(that.translate_info_),
		chain_ids_to_exclude_(that.chain_ids_to_exclude_),
		tag_along_jumps_(that.tag_along_jumps_)
{}

Translate::~Translate() {}

protocols::moves::MoverOP Translate::clone() const {
	return new Translate( *this );
}

protocols::moves::MoverOP Translate::fresh_instance() const {
	return new Translate;
}

std::string Translate::get_name() const{
	return "Translate";
}

///@brief parse XML (specifically in the context of the parser/scripting scheme)
void
Translate::parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & /*data_map*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose
)
{
	if ( tag->getName() != "Translate" ){
		throw utility::excn::EXCN_RosettaScriptsOption("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translate' mover requires chain tag");
	if ( ! tag->hasOption("distribution") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translate' mover requires distribution tag");
	if ( ! tag->hasOption("angstroms") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translate' mover requires angstroms tag");
	if ( ! tag->hasOption("cycles") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translate' mover requires cycles tag");
	//if ( ! tag->hasOption("force") ) throw utility::excn::EXCN_RosettaScriptsOption("'Translate' mover requires force tag"); optional. default is don't force, meaning ligand stays put if it can't find somewhere to go.

	std::string chain = tag->getOption<std::string>("chain");
	translate_info_.chain_id = core::pose::get_chain_id_from_chain(chain, pose);
	translate_info_.jump_id = core::pose::get_jump_id_from_chain_id(translate_info_.chain_id, pose);
	std::string distribution_str= tag->getOption<std::string>("distribution");
	translate_info_.distribution= get_distribution(distribution_str);
	translate_info_.angstroms = tag->getOption<core::Real>("angstroms");
	translate_info_.cycles = tag->getOption<core::Size>("cycles");

	if(tag->hasOption("force")){
		if(tag->getOption<std::string>("force") == "true")
			translate_info_.force= true;
		else if(tag->getOption<std::string>("force") != "false")
			throw utility::excn::EXCN_RosettaScriptsOption("'force' option is true or false");
	}

	if ( tag->hasOption("tag_along_chains") ){
		std::string const tag_along_chains_str = tag->getOption<std::string>("tag_along_chains");
		utility::vector1<std::string> tag_along_chain_strs= utility::string_split(tag_along_chains_str, ',');
		BOOST_FOREACH(std::string tag_along_chain_str, tag_along_chain_strs){
			utility::vector1<core::Size> chain_ids= get_chain_ids_from_chain(tag_along_chain_str, pose);
			BOOST_FOREACH( core::Size chain_id, chain_ids){
				core::Size jump_id= core::pose::get_jump_id_from_chain_id(chain_id, pose);
				chain_ids_to_exclude_.push_back(chain_id);
				tag_along_jumps_.push_back(jump_id);
			}
		}
	}

}

void Translate::apply(core::pose::Pose & pose) {
	core::Size const begin(pose.conformation().chain_begin(translate_info_.chain_id));

	{// add this Translate's chain conditionally (for use with CompoundTranslate)
		utility::vector1<core::Size>::iterator found=
				find(chain_ids_to_exclude_.begin(), chain_ids_to_exclude_.end(), translate_info_.chain_id);
		if( found ==  chain_ids_to_exclude_.end()	)
				chain_ids_to_exclude_.push_back(translate_info_.chain_id);
	}
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, translate_info_.jump_id);

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	if(grid_manager->size() == 0)
	{
		utility::pointer::owning_ptr<core::grid::CartGrid<int> > const grid = make_atr_rep_grid_without_ligands(pose, center, chain_ids_to_exclude_);
		translate_ligand(grid, translate_info_.jump_id, pose);// move ligand to a random point in binding pocket
	}else
	{
		//TODO refactor qsar map so it works properly
		/*
		if(!grid_manager->is_qsar_map_attached())
		{
			utility::vector1<std::string> grid_names( grid_manager->get_grid_names() );
			core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
			qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
			qsar_map->fill_with_value(1,grid_names);
			grid_manager->set_qsar_map(qsar_map);
		}
		*/
		grid_manager->initialize_all_grids(center);
		grid_manager->update_grids(pose,center);
		translate_ligand(translate_info_.jump_id,pose,begin);
	}
}

void Translate::uniform_translate_ligand(
		const utility::pointer::owning_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
){
	translate_tracer.Debug<< "making a uniform translator of up to " << translate_info_.angstroms<< " angstroms" << std::endl;
	protocols::rigid::UniformSphereTransMover mover( jump_id, translate_info_.angstroms);

	core::pose::Pose const orig_pose(pose);
	core::pose::Pose best_pose;
	int best_score=100000;

	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	for (Size cycle = 0; cycle < translate_info_.cycles; ++cycle) {
		mover.apply(pose);
		core::Vector c = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
		// did our nbr_atom land in an empty space on the grid?
		// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
		int grid_value=best_score;
		if ( grid->is_in_grid(c.x(), c.y(), c.z()) )
		{
			grid_value= grid->getValue(c.x(), c.y(), c.z());
		}
		if ( grid_value <= 0 ){
			translate_tracer.Trace << "Accepting ligand position with nbr_atom at " << c << std::endl;
			mover.freeze();
			BOOST_FOREACH(core::Size tag_along_jump, tag_along_jumps_){
				translate_tracer.Trace << "moving jump " << tag_along_jump<< " the same amount"<< std::endl;
				mover.rb_jump(tag_along_jump);
				mover.apply(pose);
			}
			return;
		}
		translate_tracer.Trace << "Rejecting ligand position with nbr_atom at " << c << std::endl;

		if( translate_info_.force && grid_value < best_score){
			best_score = grid_value;
			best_pose = pose;
		}
		pose = orig_pose; // reset and try again
	}
	// The code below only executes if we still haven't found a nonclashing position for the neighbor atom
	if( translate_info_.force ){
		translate_tracer.Trace << "Forcing neighbor atom to move (leading to a clash)" << std::endl;
		pose= best_pose;
		mover.freeze();
		BOOST_FOREACH(core::Size tag_along_jump, tag_along_jumps_){
			mover.rb_jump(tag_along_jump);
			mover.apply(pose);
		}
		return;
	}else{
		translate_tracer << "WARNING: cannot find placement for this ligand.  Keeping original position. Use the force option to force translation"<< std::endl;
	}

}

void Translate::gaussian_translate_ligand(
		const utility::pointer::owning_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
){
	translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
	protocols::rigid::RigidBodyPerturbMover mover( jump_id, 0 /*rotation*/, translate_info_.angstroms);

	core::pose::Pose const orig_pose(pose);
	core::pose::Pose best_pose;
	int best_score=100000;

	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	for (Size cycle = 0; cycle < translate_info_.cycles; ++cycle) {
		mover.apply(pose);
		core::Vector c = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
		// did our nbr_atom land in an empty space on the grid?
		// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
		int grid_value=0;
		if ( grid->is_in_grid(c.x(), c.y(), c.z()) )
		{
			grid_value= grid->getValue(c.x(), c.y(), c.z()) <= 0;
		}
		if ( grid_value < 0 ){
			translate_tracer.Trace << "Accepting ligand position with nbr_atom at " << c << std::endl;
			return;
		}
		translate_tracer.Trace << "Rejecting ligand position with nbr_atom at " << c << std::endl;

		if( translate_info_.force && grid_value < best_score){
			best_score = grid_value;
			best_pose = pose;
		}
		pose = orig_pose; // reset and try again
	}
	if( translate_info_.force ){
		pose= best_pose;
		return;
	}else{
		translate_tracer << "WARNING: cannot find placement for this ligand.  Keeping original position. Use the force option to force translation"<< std::endl;
	}
}

void Translate::translate_ligand(
		const utility::pointer::owning_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
) {
	if(translate_info_.angstroms < 0) utility_exit_with_message("cannot have a negative translation value");
	if(translate_info_.angstroms == 0) return;

	if(translate_info_.distribution == Uniform) uniform_translate_ligand(grid, jump_id, pose);
	else if(translate_info_.distribution == Gaussian) gaussian_translate_ligand(grid, jump_id, pose);
}

void Translate::translate_ligand(core::Size const jump_id,core::pose::Pose & pose, core::Size const & residue_id)
{
	if(translate_info_.angstroms < 0) utility_exit_with_message("cannot have a negative translation value");
	if(translate_info_.angstroms == 0) return;

	core::conformation::Residue const & residue(pose.residue(residue_id));

	protocols::rigid::RigidBodyMoverOP translate_mover;
	if(translate_info_.distribution == Uniform){
		translate_tracer.Debug<< "making a uniform translator of up to "<< translate_info_.angstroms<<" angstroms"<< std::endl;
		translate_mover= new protocols::rigid::UniformSphereTransMover( jump_id, translate_info_.angstroms);
	}
	else if(translate_info_.distribution == Gaussian){
		translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
		translate_mover= new protocols::rigid::RigidBodyPerturbMover ( jump_id, 0 /*rotation*/, translate_info_.angstroms);
	}

	RandomConformerMoverOP conformer_mover = new RandomConformerMover(residue_id);

	core::pose::Pose  orig_pose(pose);
	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	core::Real best_score = 1000000;

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	for (Size cycle = 0; cycle < translate_info_.cycles; ++cycle) {

			translate_tracer.Trace <<"Doing a translate move" <<std::endl;
			translate_mover->apply(pose);
	//	}else
	//	{
	//
	//		translate_tracer.Trace <<"Doing a conformer selection move" <<std::endl;
	//		conformer_mover->apply(pose);
	//	}
		core::Real score(grid_manager->total_score(residue));
		if(score <= best_score)
		{
			best_score = score;
			orig_pose = pose;
			translate_tracer.Trace << "accepted ligand position with score of " <<score <<std::endl;
		}else
		{
			pose = orig_pose;
			translate_tracer.Trace <<"rejected ligand position with score of "<<score <<std::endl;
		}
	}

}

core::Size
Translate::get_chain_id(core::pose::Pose const & ){
	return translate_info_.chain_id;
}

void
Translate::add_excluded_chains(
		std::set<core::Size>::const_iterator begin,
		std::set<core::Size>::const_iterator end
){
		chain_ids_to_exclude_.insert( chain_ids_to_exclude_.end(), begin, end);
}


} //namespace ligand_docking
} //namespace protocols
