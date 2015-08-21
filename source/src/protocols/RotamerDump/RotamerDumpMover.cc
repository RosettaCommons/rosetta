// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/RotamerDump/RotamerDumpMover.cc
/// @author Sam DeLuca

#include <typeinfo>

#include <protocols/RotamerDump/RotamerDumpMover.hh>

#include <core/pack/interaction_graph/PrecomputedPairEnergiesInteractionGraph.hh>
#include <core/pack/interaction_graph/InteractionGraphFactory.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/graph/Graph.hh>
#include <basic/Tracer.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/conformation/Residue.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/rotamerdump.OptionKeys.gen.hh>

#include <utility/string_util.hh>

#include <protocols/jd2/JobDistributor.hh>

#include <core/pose/Pose.hh>
#include <protocols/jd2/Job.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace RotamerDump {

static thread_local basic::Tracer RotamerTracer( "protocols.RotamerDump.RotamerDumpMover" );

RotamerDumpMover::RotamerDumpMover(core::pack::task::TaskFactoryOP task_factory,
	core::scoring::ScoreFunctionOP score_function) :
	task_factory_(task_factory), score_function_(score_function)
{

}

void RotamerDumpMover::apply(core::pose::Pose & pose)
{

	core::pack::pack_scorefxn_pose_handshake(pose,*score_function_);

	pose.update_residue_neighbors();


	core::pack::task::PackerTaskCOP packer_task( task_factory_->create_task_and_apply_taskoperations(pose));
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets( new core::pack::rotamer_set::RotamerSets() );

	score_function_->setup_for_packing(pose, packer_task->repacking_residues(), packer_task->designing_residues());

	core::graph::GraphOP packer_neighbor_graph(core::pack::create_packer_graph(pose,*score_function_,packer_task));

	rotamer_sets->set_task(packer_task);
	rotamer_sets->build_rotamers(pose,*score_function_,packer_neighbor_graph);

	rotamer_sets->prepare_sets_for_packing(pose,*score_function_);

	RotamerTracer << "built " << rotamer_sets->nrotamers() << " rotamers at " <<rotamer_sets->nmoltenres() <<" positions" <<std::endl;

	core::pack::interaction_graph::InteractionGraphBaseOP ig(core::pack::interaction_graph::InteractionGraphFactory::create_interaction_graph(*packer_task,*rotamer_sets,pose,*score_function_));
	rotamer_sets->compute_energies(pose,*score_function_,packer_neighbor_graph,ig);

	RotamerTracer << "IG: " <<ig->getTotalMemoryUsage() << " bytes"<<std::endl;

	jd2::JobOP job(jd2::JobDistributor::get_instance()->current_job());

	std::list<std::string> job_data;

	if ( basic::options::option[basic::options::OptionKeys::rotamerdump::one_body].user() ) {
		RotamerTracer << "Generating one body energy table"<<std::endl;
		job_data.push_back(get_onebody_energy_table(ig,rotamer_sets));
	}
	if ( basic::options::option[basic::options::OptionKeys::rotamerdump::xyz].user() ) {
		RotamerTracer << "Generating XYZ rotamer table"<<std::endl;
		job_data.push_back(get_xyz_coord_table(rotamer_sets));
	}
	if ( basic::options::option[basic::options::OptionKeys::rotamerdump::two_body].user() ) {
		RotamerTracer << "Generating two body energy table"<<std::endl;
		job_data.push_back(get_twobody_energy_table(ig,rotamer_sets));
	}
	if ( basic::options::option[basic::options::OptionKeys::rotamerdump::annealer].user() ) {
		RotamerTracer <<"Running Annealer" <<std::endl;
		job_data.push_back(get_annealer_pick_table(ig,rotamer_sets,pose,packer_task));
	}
	job->add_strings(job_data);
	//RotamerTracer <<two_body_energy <<std::endl;
	//RotamerTracer << one_body_energy <<std::endl;

}

/// @details appends a line to the job in the form one_body num_items (resno, resn, rotno, energy) for each 1 body energy in the IG
std::string RotamerDumpMover::get_onebody_energy_table(core::pack::interaction_graph::InteractionGraphBaseOP ig ,
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets)
{
	std::string table_type("one_body");
	core::Size ig_size = ig->get_num_nodes();
	std::string data ="";
	core::Size elements = 0;
	//table.append(" "+utility::to_string<core::Size>(ig_size));
	for ( core::Size node_id = 1; node_id <= ig_size; ++node_id ) {
		core::Size residue_id = rotamer_sets->moltenres_2_resid(node_id);
		std::string residue_id_string = utility::to_string<core::Size>(residue_id);


		core::Size node_states = ig->get_num_states_for_node(node_id);
		for ( core::Size state_id = 1; state_id <= node_states; ++state_id ) {
			float one_body_energy = ig->get_one_body_energy_for_node_state(node_id,state_id);
			std::string one_body_energy_string = utility::to_string<float>(one_body_energy);

			core::pack::rotamer_set::RotamerSetOP rotamer_set(rotamer_sets->rotamer_set_for_moltenresidue(node_id));
			core::conformation::ResidueCOP current_residue(rotamer_set->rotamer(state_id));


			std::string residue_name(current_residue->name3());
			std::string state_id_string = utility::to_string<core::Size>(state_id);

			//RotamerTracer << residue_id_string << " " << residue_name << " " << state_id_string <<  " " << one_body_energy_string <<std::endl;
			data.append(" ("+residue_id_string+","+residue_name+","+state_id_string+","+one_body_energy_string+")");
			elements++;
		}
	}
	return table_type+" "+utility::to_string<core::Size>(elements)+" "+data;
}

/// @details appends a line to the job in the form two_body num_items (nres1,nrot1,resn1,nres2,nrot2,resn2,energy) for each 2 body energy in the IG
std::string RotamerDumpMover::get_twobody_energy_table(core::pack::interaction_graph::InteractionGraphBaseOP ig,
	core::pack::rotamer_set::RotamerSetsOP rotamer_sets)
{
	std::string table_type("two_body");
	std::string data = "";
	core::Size elements = 0;

	core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraphOP p_ig =
		utility::pointer::dynamic_pointer_cast< core::pack::interaction_graph::PrecomputedPairEnergiesInteractionGraph > ( ig );

	if ( !p_ig ) {
		utility_exit_with_message("Interaction graph is not pre-computed");
	}
	core::Size ig_size = ig->get_num_nodes();
	for ( core::Size node_1_id = 1; node_1_id <= ig_size; ++node_1_id ) {
		core::Size residue_1_id = rotamer_sets->moltenres_2_resid(node_1_id);
		std::string residue_1_id_string = utility::to_string<core::Size>(residue_1_id);


		for ( core::Size node_2_id = node_1_id+2; node_2_id <= ig_size; ++node_2_id ) {
			core::Size residue_2_id = rotamer_sets->moltenres_2_resid(node_2_id);
			std::string residue_2_id_string = utility::to_string<core::Size>(residue_2_id);


			core::Size const num_states_1(ig->get_num_states_for_node(node_1_id));
			core::Size const num_states_2(ig->get_num_states_for_node(node_2_id));
			if ( !ig->get_edge_exists(node_1_id,node_2_id) ) {
				continue;
			}

			for ( core::Size node_1_state = 1; node_1_state <= num_states_1; ++node_1_state ) {
				std::string node_1_state_string = utility::to_string<core::Size>(node_1_state);

				core::pack::rotamer_set::RotamerSetOP rotamer_set_1(rotamer_sets->rotamer_set_for_moltenresidue(node_1_id));
				core::conformation::ResidueCOP residue_1(rotamer_set_1->rotamer(node_1_state));
				std::string residue_1_name(residue_1->name3());

				for ( core::Size node_2_state = 1; node_2_state <= num_states_2; ++node_2_state ) {
					core::pack::rotamer_set::RotamerSetOP rotamer_set_2(rotamer_sets->rotamer_set_for_moltenresidue(node_2_id));
					core::conformation::ResidueCOP residue_2(rotamer_set_2->rotamer(node_2_state));
					std::string residue_2_name(residue_2->name3());

					std::string node_2_state_string = utility::to_string<core::Size>(node_2_state);
					float energy = p_ig->get_two_body_energy_for_edge(node_1_id,node_2_id,node_1_state,node_2_state);
					if ( energy != 0.0 ) {
						std::string energy_string = utility::to_string<float>(energy);
						data.append(" ("+residue_1_id_string+","+node_1_state_string+","+residue_1_name+","+residue_2_id_string+","+node_2_state_string+","+residue_2_name+","+energy_string+")" );
						elements++;
					}
				}
			}
		}
	}
	return table_type+" "+utility::to_string<core::Size>(elements)+" "+data;
}

/// @details appends a line to the job in the form xyz_coord num_items xyz_coord (resno,rotno,resn,atomName,x,y,z) for each atom in the IG
std::string RotamerDumpMover::get_xyz_coord_table(core::pack::rotamer_set::RotamerSetsOP rotamer_sets)
{
	std::string table_type("xyz_coord");
	std::string data = "";
	core::Size elements = 0;

	core::Size n_molten_res =rotamer_sets->nmoltenres();
	for ( core::Size molten_res =1; molten_res <= n_molten_res; ++molten_res ) {
		core::pack::rotamer_set::RotamerSetOP rotamer_set(rotamer_sets->rotamer_set_for_moltenresidue(molten_res));
		core::Size n_rotamers = rotamer_set->num_rotamers();
		std::string res_id_string = utility::to_string<core::Size>(rotamer_sets->moltenres_2_resid(molten_res));

		for ( core::Size rotamer_id =1; rotamer_id <= n_rotamers; ++rotamer_id ) {
			std::string rotamer_id_string = utility::to_string<core::Size>(rotamer_id);

			core::conformation::ResidueCOP current_residue(rotamer_set->rotamer(rotamer_id));
			std::string residue_name(current_residue->name3());

			core::Size n_atoms = current_residue->natoms();
			for ( core::Size atom_id = 1; atom_id <= n_atoms; ++atom_id ) {
				std::string atom_name = utility::trim(current_residue->atom_name(atom_id));

				core::Vector xyz_coords = current_residue->xyz(atom_id);
				std::string x_string = utility::to_string<core::Real>(xyz_coords.x());
				std::string y_string = utility::to_string<core::Real>(xyz_coords.y());
				std::string z_string = utility::to_string<core::Real>(xyz_coords.z());

				data.append(" ("+res_id_string+","+rotamer_id_string+","+residue_name+","+atom_name+","+x_string+","+y_string+","+z_string+")");
				elements++;
			}
		}
	}
	return table_type+" "+utility::to_string<core::Size>(elements)+" "+data;
}

/// @details appends a line to the job in the form annealer_results num_items (nres,selected_nres, selected_resn) for each atom in the IG
std::string RotamerDumpMover::get_annealer_pick_table(core::pack::interaction_graph::InteractionGraphBaseOP ig, core::pack::rotamer_set::RotamerSetsOP rotamer_sets,core::pose::Pose & pose , core::pack::task::PackerTaskCOP task)
{
	std::string table_type("annealer_results");
	std::string data = "";
	core::Size elements = 0;

	//bool start_with_current(false);
	ObjexxFCL::FArray1D_int current_rot_index(pose.total_residue(),0);
	//bool calc_rot_freq(false);
	ObjexxFCL::FArray1D<float> rot_freq(ig->get_num_total_states(),0.0);
	ObjexxFCL::FArray1D_int bestrotamer_at_seqpos(pose.total_residue());
	utility::vector0<int> rot_to_pack;
	float bestenergy = 0;

	core::pack::pack_rotamers_run(pose,task,rotamer_sets,ig,rot_to_pack,bestrotamer_at_seqpos,bestenergy);

	core::Size n_molten_res =rotamer_sets->nmoltenres();
	for ( core::Size molten_res =1; molten_res <= n_molten_res; ++molten_res ) {
		core::Size residue_id = rotamer_sets->moltenres_2_resid(molten_res);
		std::string residue_id_string = utility::to_string<core::Size>(residue_id);

		core::Size selected_rotamer = rotamer_sets->rotid_on_moltenresidue(bestrotamer_at_seqpos(residue_id));
		std::string selected_rotamer_string = utility::to_string<core::Size>(selected_rotamer);

		core::conformation::ResidueCOP selected_residue(rotamer_sets->rotamer_set_for_moltenresidue(molten_res)->rotamer(selected_rotamer));
		std::string selected_residue_name(selected_residue->name3());

		data.append(" ("+residue_id_string+","+selected_rotamer_string+","+selected_residue_name+")");
		elements++;
	}
	return table_type+" "+utility::to_string<core::Size>(elements)+" "+data;
}


std::string RotamerDumpMover::get_name() const
{
	return "RotamerDumpMover";
}

}
}
