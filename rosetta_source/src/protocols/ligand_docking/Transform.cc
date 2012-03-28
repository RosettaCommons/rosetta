// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   /grid/src/protocols/ligand_docking/Transform.cc
/// @author Sam DeLuca


#include <protocols/ligand_docking/TransformCreator.hh>
#include <protocols/ligand_docking/Transform.hh>
#include <protocols/ligand_docking/grid_functions.hh>

#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>
#include <basic/Tracer.hh>


#include <core/pose/util.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
// AUTO-REMOVED #include <protocols/qsar/qsarMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>
//#include <utility/string_util.hh>

namespace protocols {
namespace ligand_docking {

static numeric::random::RandomGenerator RG(23459);

static basic::Tracer transform_tracer("protocols.ligand_docking.ligand_options.transform", basic::t_debug);

std::string TransformCreator::keyname() const
{
	return TransformCreator::mover_name();
}

protocols::moves::MoverOP TransformCreator::create_mover() const
{
	return new Transform;
}

std::string TransformCreator::mover_name()
{
	return "Transform";
}

Transform::Transform(): Mover("Transform"), transform_info_(),optimize_until_score_is_negative_(0.0)
{

}

Transform::Transform(
	std::string const & chain,
	core::Real const & box_size,
	core::Real const & move_distance,
	core::Real const & angle,
	core::Size const & cycles,
	core::Real const & temp
) : Mover("Transform"), transform_info_(),optimize_until_score_is_negative_(0.0)
{
	transform_info_.chain = chain;
	transform_info_.box_size = box_size;
	transform_info_.move_distance = move_distance;
	transform_info_.angle = angle;
	transform_info_.cycles = cycles;
	transform_info_.temperature = temp;
}

Transform::~Transform()
{
	//
}

protocols::moves::MoverOP Transform::clone() const
{
	return new Transform (*this);
}

protocols::moves::MoverOP Transform::fresh_instance() const
{
	return new Transform;
}

std::string Transform::get_name() const
{
	return "Transform";
}

void Transform::parse_my_tag
(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & data_map,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/
)
{
	if ( tag->getName() != "Transform" )
	{
		utility_exit_with_message("This should be impossible");
	}
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'Transform' mover requires chain tag");
	if ( ! tag->hasOption("move_distance") ) utility_exit_with_message("'Transform' mover requires move_distance tag");
	if (! tag->hasOption("box_size") ) utility_exit_with_message("'Transform' mover requires box_size tag");
	if ( ! tag->hasOption("angle") ) utility_exit_with_message("'Transform' mover requires angle tag");
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'Transform' mover requires cycles tag");
	if (!tag->hasOption("temperature")) utility_exit_with_message("'Transform' mover requires temperature tag");


	transform_info_.chain = tag->getOption<std::string>("chain");
	transform_info_.move_distance = tag->getOption<core::Real>("move_distance");
	transform_info_.box_size = tag->getOption<core::Real>("box_size");
	transform_info_.angle = tag->getOption<core::Real>("angle");
	transform_info_.cycles = tag->getOption<core::Size>("cycles");
	transform_info_.temperature = tag->getOption<core::Real>("temperature");

	optimize_until_score_is_negative_ = tag->getOption<bool>("optimize_until_score_is_negative",false);

}

void Transform::apply(core::pose::Pose & pose)
{
	qsar::scoring_grid::GridManager* grid_manager(qsar::scoring_grid::GridManager::get_instance());

	assert(transform_info_.chain.size() == 1);
	transform_info_.chain_id = core::pose::get_chain_id_from_chain(transform_info_.chain, pose);
	transform_info_.jump_id = core::pose::get_jump_id_from_chain_id(transform_info_.chain_id, pose);
	core::Size const begin(pose.conformation().chain_begin(transform_info_.chain_id));
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose, transform_info_.jump_id));
	assert(grid_manager != 0); //something has gone hopelessly wrong if this triggers


	//TODO qsar map needs to be refactored so it actually works right
	/*
	if(!grid_manager->is_qsar_map_attached())
	{
		core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
		qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
		if(!qsar_map->fill_from_mol_data(residue->type().get_mol_data()))
		{
			utility::vector1<std::string> grid_names( grid_manager->get_grid_names() );
			qsar_map->fill_with_value(1,grid_names);
		}
		//		qsar_map->fill_with_value(1);
		grid_manager->set_qsar_map(qsar_map);
	}
	*/
	core::conformation::Residue original_residue = pose.residue(begin);
	core::chemical::ResidueType residue_type = pose.residue_type(begin);
	//jd2::JobOP current_job = protocols::jd2::JobDistributor::get_instance()->current_job();
	//std::string tag = current_job->input_tag();

	grid_manager->initialize_all_grids(center);
	grid_manager->update_grids(pose,center);

	core::pose::Pose best_pose(pose);
	core::Real best_score(grid_manager->total_score(original_residue));
	core::Real temperature = transform_info_.temperature;
	core::Vector original_center(original_residue.xyz(original_residue.nbr_atom()));

	rotamers_for_trials(pose,begin,ligand_conformers_);
	transform_tracer << "Considering " << ligand_conformers_.size() << " conformers during sampling" << std::endl;
	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;

	bool not_converged = true;

	core::Size cycle = 1;
	//for(core::Size cycle = 1; cycle <= transform_info_.cycles; ++cycle)
	while(not_converged)
	{
		if(optimize_until_score_is_negative_)
		{
			if(cycle >= transform_info_.cycles && best_score <= 0.0)
			{
				not_converged= false;
			}
		}else
		{
			if(cycle >= transform_info_.cycles)
			{
				not_converged= false;
			}
		}

		cycle++;

		//during each move either move the ligand or try a new conformer (if there is more than one conformer)
		if(ligand_conformers_.size() > 1)
		{
			if(RG.uniform() >= 0.5)
			{
				transform_ligand(pose);
			}else
			{
				change_conformer(pose,begin);
			}
		}else
		{
			transform_ligand(pose);
		}
		//delete residue;
		core::conformation::Residue residue(pose.residue(begin));

		core::Real current_score = grid_manager->total_score(residue);
		core::Real const boltz_factor((best_score-current_score)/temperature);
		core::Real const probability = std::exp( boltz_factor ) ;
		core::Vector new_center(residue.xyz(residue.nbr_atom()));


		if(new_center.distance(original_center) > transform_info_.box_size) //Reject the new pose
		{
			pose = best_pose;
			rejected_moves++;
			//transform_tracer << "probability: " << probability << " rejected (out of box)"<<std::endl;
		}else if(probability < 1 && RG.uniform() >= probability)  //reject the new pose
		{
			pose = best_pose;
			rejected_moves++;
			//transform_tracer << "probability: " << probability << " rejected"<<std::endl;
		}else if(probability < 1)  // Accept the new pose
		{
			best_score = current_score;
			best_pose = pose;
			accepted_moves++;
			//transform_tracer << "probability: " << probability << " accepted"<<std::endl;
		}else  //Accept the new pose
		{
			best_score = current_score;
			best_pose = pose;
			accepted_moves++;
			//transform_tracer << "probability: " << probability << " accepted"<<std::endl;
		}
		transform_tracer << best_score << " " <<current_score <<std::endl;
	}
	pose=best_pose;

	transform_tracer <<"percent acceptance: "<< accepted_moves << " " << (core::Real)accepted_moves/(core::Real)rejected_moves <<" " << rejected_moves <<std::endl;

}

void Transform::transform_ligand(core::pose::Pose & pose)
{
	if(transform_info_.angle ==0 && transform_info_.move_distance == 0)
	{
		transform_tracer <<"WARNING: angle and distance are both 0.  Transform will do nothing" <<std::endl;
		return;
	}

	protocols::rigid::RigidBodyMoverOP mover;
	mover = new protocols::rigid::RigidBodyPerturbMover(transform_info_.jump_id,transform_info_.angle,transform_info_.move_distance);
	//core::Size chain_begin = pose.conformation().chain_begin(transform_info_.chain_id);

	mover->apply(pose);
	pose.update_actcoords();
}

void Transform::change_conformer(core::pose::Pose & pose, core::Size const & seqpos)
{
	assert(ligand_conformers_.size());
	core::Size index_to_select = RG.random_range(1,ligand_conformers_.size());
	core::conformation::ResidueOP new_residue = ligand_conformers_[index_to_select];
	pose.conformation().replace_residue(seqpos,*new_residue,false );
	pose.update_actcoords();
}


}
}
