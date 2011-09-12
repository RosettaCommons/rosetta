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
#include <protocols/geometry/RB_geometry.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/moves/DataMap.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <protocols/qsar/scoring_grid/GridManager.hh>

#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/tag/Tag.hh>
#include <numeric/random/random.hh>

//Auto Headers
#include <core/grid/CartGrid.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <numeric/xyzVector.io.hh>

using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace ligand_docking {

static basic::Tracer translate_tracer("protocols.ligand_docking.ligand_options.translate", basic::t_debug);

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
		chain_ids_to_exclude_()
{}

Translate::Translate(Translate const & that):
		//utility::pointer::ReferenceCount(),
		protocols::moves::Mover( that ),
		translate_info_(that.translate_info_),
		chain_ids_to_exclude_(that.chain_ids_to_exclude_)
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
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data_map,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "Translate" ){
		utility_exit_with_message("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'Translate' mover requires chain tag");
	if ( ! tag->hasOption("distribution") ) utility_exit_with_message("'Translate' mover requires distribution tag");
	if ( ! tag->hasOption("angstroms") ) utility_exit_with_message("'Translate' mover requires angstroms tag");
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'Translate' mover requires cycles tag");

	translate_info_.chain = tag->getOption<std::string>("chain");
	{
		std::string distribution_str= tag->getOption<std::string>("distribution");
		translate_info_.distribution= get_distribution(distribution_str);
	}
	translate_info_.angstroms = tag->getOption<core::Real>("angstroms");
	translate_info_.cycles = tag->getOption<core::Size>("cycles");
}

void Translate::apply(core::pose::Pose & pose) {
	assert(translate_info_.chain.size() == 1);// do the user check in parse_flags. Later remove this requirement

	core::Size const chain_id = core::pose::get_chain_id_from_chain(translate_info_.chain, pose);
	core::Size const jump_id = core::pose::get_jump_id_from_chain_id(chain_id, pose);
	core::Size const begin(pose.conformation().chain_begin(chain_id));

	{// add this Translate's chain conditionally (for use with CompoundTranslate)
		utility::vector1<core::Size>::iterator found=
				find(chain_ids_to_exclude_.begin(), chain_ids_to_exclude_.end(), chain_id);
		if( found ==  chain_ids_to_exclude_.end()	)
				chain_ids_to_exclude_.push_back(chain_id);
	}
	core::Vector const center = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	if(grid_manager->size() == 0)
	{
		utility::pointer::owning_ptr<core::grid::CartGrid<int> > const grid = make_atr_rep_grid_without_ligands(pose, center, chain_ids_to_exclude_);
		translate_ligand(grid, jump_id, pose);// move ligand to a random point in binding pocket
	}else
	{
		if(!grid_manager->is_qsar_map_attached())
		{
			utility::vector1<std::string> grid_names( grid_manager->get_grid_names() );
			core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
			qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
			qsar_map->fill_with_value(1,grid_names);
			grid_manager->set_qsar_map(qsar_map);
		}
		grid_manager->initialize_all_grids(center);
		grid_manager->update_grids(pose,center);
		translate_ligand(jump_id,pose,begin);
	}
}

void Translate::translate_ligand(
		const utility::pointer::owning_ptr<core::grid::CartGrid<int> >  & grid,
		const core::Size jump_id,
		core::pose::Pose & pose
) {
	if(translate_info_.angstroms < 0) utility_exit_with_message("cannot have a negative translation value");
	if(translate_info_.angstroms == 0) return;

	protocols::moves::RigidBodyMoverOP mover;
	if(translate_info_.distribution == Uniform){
		translate_tracer.Debug<< "making a uniform translator of up to "<< translate_info_.angstroms<<" angstroms"<< std::endl;
		mover= new protocols::moves::UniformSphereTransMover( jump_id, translate_info_.angstroms);
	}
	else if(translate_info_.distribution == Gaussian){
		translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
		mover= new protocols::moves::RigidBodyPerturbMover ( jump_id, 0 /*rotation*/, translate_info_.angstroms);
	}

	core::pose::Pose const orig_pose(pose);

	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	for (Size cycle = 0; cycle < translate_info_.cycles; ++cycle) {
		mover->apply(pose);
		core::Vector c = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
		// did our nbr_atom land in an empty space on the grid?
		// Don't want to insist the nbr_atom is in the attractive region, because it may not be!
		if (
				grid->is_in_grid(c.x(), c.y(), c.z())
				&& grid->getValue(c.x(), c.y(), c.z()) <= 0
		) {
			translate_tracer.Trace << "Accepting ligand position with nbr_atom at " << c << std::endl;
			return;
		}
		translate_tracer.Trace << "Rejecting ligand position with nbr_atom at " << c << std::endl;
		pose = orig_pose; // reset and try again
	}///TODO alert the user that we cannot find a placement for this ligand.  Do not proceed with docking this ligand.

	translate_tracer << "WARNING: cannot find placement for this ligand"<< std::endl;
}

void Translate::translate_ligand(core::Size const jump_id,core::pose::Pose & pose, core::Size const & residue_id)
{
	if(translate_info_.angstroms < 0) utility_exit_with_message("cannot have a negative translation value");
	if(translate_info_.angstroms == 0) return;

	core::conformation::Residue const & residue(pose.residue(residue_id));

	protocols::moves::RigidBodyMoverOP translate_mover;
	if(translate_info_.distribution == Uniform){
		translate_tracer.Debug<< "making a uniform translator of up to "<< translate_info_.angstroms<<" angstroms"<< std::endl;
		translate_mover= new protocols::moves::UniformSphereTransMover( jump_id, translate_info_.angstroms);
	}
	else if(translate_info_.distribution == Gaussian){
		translate_tracer.Debug<< "making a Gaussian translator of up to "<< translate_info_.angstroms<<" angstroms";
		translate_mover= new protocols::moves::RigidBodyPerturbMover ( jump_id, 0 /*rotation*/, translate_info_.angstroms);
	}

	RandomConformerMoverOP conformer_mover = new RandomConformerMover(residue_id);

	core::pose::Pose  orig_pose(pose);
	translate_tracer.Debug << "time to cycle: " << translate_info_.cycles << std::endl;
	core::Real best_score = 1000000;

	qsar::scoring_grid::GridManager* grid_manager = qsar::scoring_grid::GridManager::get_instance();
	for (Size cycle = 0; cycle < translate_info_.cycles; ++cycle) {
		core::Real pick_move(numeric::random::gaussian());
		if(pick_move >= 0.50)
		{
			translate_tracer.Trace <<"Doing a translate move" <<std::endl;
			translate_mover->apply(pose);
		}else
		{

			translate_tracer.Trace <<"Doing a conformer selection move" <<std::endl;
			conformer_mover->apply(pose);
		}
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
Translate::get_chain_id(core::pose::Pose const & pose){
	core::Size chain_id = core::pose::get_chain_id_from_chain(translate_info_.chain, pose);
	return chain_id;
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
