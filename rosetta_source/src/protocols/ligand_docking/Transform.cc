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


#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>


#include <core/pose/util.hh>
#include <protocols/qsar/scoring_grid/GridManager.hh>
#include <protocols/qsar/qsarMap.hh>
#include <protocols/moves/RigidBodyMover.hh>
#include <protocols/geometry/RB_geometry.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <protocols/moves/DataMap.hh>
#include <utility/tag/Tag.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
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

Transform::Transform(): Mover("Transform"), grid_manager_(0), transform_info_()
{

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
	//tag->write(transform_tracer,1);
	if ( ! tag->hasOption("chain") ) utility_exit_with_message("'Transform' mover requires chain tag");
	//if ( ! tag->hasOption("distribution") ) utility_exit_with_message("'Translate' mover requires distribution tag");
	if ( ! tag->hasOption("move_distance") ) utility_exit_with_message("'Transform' mover requires move_distance tag");
	if (! tag->hasOption("box_size") ) utility_exit_with_message("'Transform' mover requires box_size tag");
	if ( ! tag->hasOption("angle") ) utility_exit_with_message("'Transform' mover requires angle tag");
	if ( ! tag->hasOption("cycles") ) utility_exit_with_message("'Transform' mover requires cycles tag");
	if (!tag->hasOption("temperature")) utility_exit_with_message("'Transform' mover requires temperature tag");
	if(data_map.has("scoringgrid","default"))
	{
		grid_manager_ = data_map.get<qsar::scoring_grid::GridManager *>("scoringgrid","default");
	}
	else
	{
		utility_exit_with_message("ERROR: you are using Transform without specifying a GridManger");

		//grid_manager_ = new qsar::scoring_grid::GridManager(40.0, 0.25);
	}

	transform_info_.chain = tag->getOption<std::string>("chain");
	//std::string distribution_str = tag->getOption<std::string>("distribution");
	//transform_info_.distribution = get_distribution(distribution_str);
	transform_info_.move_distance = tag->getOption<core::Real>("move_distance");
	transform_info_.box_size = tag->getOption<core::Real>("box_size");
	transform_info_.angle = tag->getOption<core::Real>("angle");
	transform_info_.cycles = tag->getOption<core::Size>("cycles");
	transform_info_.temperature = tag->getOption<core::Real>("temperature");
}

void Transform::apply(core::pose::Pose & pose)
{
	assert(transform_info_.chain.size() == 1);
	transform_info_.chain_id = core::pose::get_chain_id_from_chain(transform_info_.chain, pose);
	transform_info_.jump_id = core::pose::get_jump_id_from_chain_id(transform_info_.chain_id, pose);
	core::Size const begin(pose.conformation().chain_begin(transform_info_.chain_id));
	core::Vector const center(protocols::geometry::downstream_centroid_by_jump(pose, transform_info_.jump_id));
	assert(grid_manager_ != 0); //something has gone hopelessly wrong if this triggers


	if(!grid_manager_->is_qsar_map_attached())
	{
		core::conformation::ResidueOP residue = new core::conformation::Residue(pose.residue(begin));
		qsar::qsarMapOP qsar_map(new qsar::qsarMap("default",residue));
		if(!qsar_map->fill_from_mol_data(residue->type().get_mol_data()))
		{
			utility::vector1<std::string> grid_names( grid_manager_->get_grid_names() );
			qsar_map->fill_with_value(1,grid_names);
		}
		//		qsar_map->fill_with_value(1);
		grid_manager_->set_qsar_map(qsar_map);
	}

	core::conformation::Residue original_residue = pose.residue(begin);

	jd2::JobOP current_job = protocols::jd2::JobDistributor::get_instance()->current_job();
	std::string tag = current_job->input_tag();

	grid_manager_->initialize_all_grids(center);
	grid_manager_->update_grids(pose,center,tag);

	core::pose::Pose best_pose(pose);
	core::Real best_score(grid_manager_->total_score(original_residue));
	core::Real temperature = transform_info_.temperature;
	core::Vector original_center(original_residue.xyz(original_residue.nbr_atom()));

	core::Size accepted_moves = 0;
	core::Size rejected_moves = 0;

	for(core::Size cycle = 1; cycle <= transform_info_.cycles; ++cycle)
	{
		transform_ligand(pose);

		//delete residue;
		core::conformation::Residue residue(pose.residue(begin));

		core::Real current_score = grid_manager_->total_score(residue);
		core::Real const boltz_factor((best_score-current_score)/temperature);
		core::Real const probability = std::exp( boltz_factor ) ;
		core::Vector new_center(residue.xyz(residue.nbr_atom()));
		//std::ostringstream oss;

		//oss << std::setfill('0') << std::setw(4) << cycle;
		//std::string cycle_string(oss.str());

		//const std::string accept_path(tag+"_accept_"+cycle_string+".pdb");
		//const std::string reject_path(tag+"_reject_"+cycle_string+".pdb");
		if(new_center.distance(original_center) > transform_info_.box_size) //Reject the new pose
		{
			//pose.dump_pdb(reject_path);
			pose = best_pose;
			rejected_moves++;
			transform_tracer << "probability: " << probability << " rejected (out of box)"<<std::endl;
		}else if(probability < 1 && RG.uniform() >= probability)  //reject the new pose
		{
			//pose.dump_pdb(reject_path);
			pose = best_pose;
			rejected_moves++;
			transform_tracer << "probability: " << probability << " rejected"<<std::endl;
		}else if(probability < 1)  // Accept the new pose
		{
			//pose.dump_pdb(accept_path);
			best_score = current_score;
			best_pose = pose;
			accepted_moves++;
			transform_tracer << "probability: " << probability << " accepted"<<std::endl;
		}else  //Accept the new pose
		{
			//pose.dump_pdb(accept_path);
			best_score = current_score;
			best_pose = pose;
			accepted_moves++;
			transform_tracer << "probability: " << probability << " accepted"<<std::endl;
			//pose = best_pose;
			//rejected_moves++;
		}
	}
	pose=best_pose;

	transform_tracer <<"percent acceptance: "<< accepted_moves << " " << accepted_moves/rejected_moves <<" " << rejected_moves <<std::endl;

}

void Transform::transform_ligand(core::pose::Pose & pose)
{
	if(transform_info_.angle ==0 && transform_info_.move_distance == 0)
	{
		transform_tracer <<"WARNING: angle and distance are both 0.  Transform will do nothing" <<std::endl;
		return;
	}

	protocols::moves::RigidBodyMoverOP mover;
	mover = new protocols::moves::RigidBodyPerturbMover(transform_info_.jump_id,transform_info_.angle,transform_info_.move_distance);
	//core::Size chain_begin = pose.conformation().chain_begin(transform_info_.chain_id);

	mover->apply(pose);
	pose.update_actcoords();
}

void Transform::change_conformer
(
	core::Size const ,
	core::pose::Pose & ,
	core::Size const &
)
{
	//lets implement this later
}

}
}
