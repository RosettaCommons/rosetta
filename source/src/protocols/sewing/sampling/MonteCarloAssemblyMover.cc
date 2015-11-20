// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MonteCarloAssemblyMover.cc
///
/// @brief Assembly substructures by MonteCarlo way
/// @author Tim Jacobs, Frank Teets
/// @modified Doonam Kim

// Unit Headers
#include <protocols/sewing/sampling/MonteCarloAssemblyMover.hh>
#include <protocols/sewing/sampling/MonteCarloAssemblyMoverCreator.hh>

// Package Headers
#include <protocols/sewing/conformation/AssemblyFactory.hh>
#include <protocols/sewing/sampling/requirements/RequirementSet.hh>
#include <protocols/sewing/scoring/AssemblyScorer.hh>

//Utility headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/sewing.OptionKeys.gen.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <utility/tag/Tag.hh>

namespace protocols {
namespace sewing  {

static basic::Tracer TR( "protocols.sewing.sampling.MonteCarloAssemblyMover" );

////////////////////////////////////////////////////////////////////////////////////
///////////////////////////  Boiler Plate Code   ///////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
protocols::moves::MoverOP
MonteCarloAssemblyMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new MonteCarloAssemblyMover );
}

std::string
MonteCarloAssemblyMoverCreator::keyname() const
{
	return MonteCarloAssemblyMoverCreator::mover_name();
}

std::string
MonteCarloAssemblyMoverCreator::mover_name()
{
	return "MonteCarloAssemblyMover";
}

protocols::moves::MoverOP
MonteCarloAssemblyMover::clone() const {
	return( protocols::moves::MoverOP( new MonteCarloAssemblyMover( *this ) ) );
}

protocols::moves::MoverOP
MonteCarloAssemblyMover::fresh_instance() const {
	return protocols::moves::MoverOP( new MonteCarloAssemblyMover );
}

std::string
MonteCarloAssemblyMover::get_name() const {
	return "MonteCarloAssemblyMover";
}

/////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////  MonteCarloAssemblyMover function   //////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////

MonteCarloAssemblyMover::MonteCarloAssemblyMover():
	cycles_(10000),
	max_temperature_(0.10),
	min_temperature_(0.05),
	min_assembly_score_(-2.0),
	use_best_assembly_score_(true),
	add_probability_(0.05),
	delete_probability_(0.005),
	switch_probability_(0.8)
{}

AssemblyOP
MonteCarloAssemblyMover::generate_assembly(){
	// TR << "MonteCarloAssemblyMover::generate_assembly()" << std::endl;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//Initialize an empty starting Assembly and add it to the list
	AssemblyOP working_assembly;
	if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "discontinuous" ) {
		working_assembly = AssemblyFactory::create_assembly("discontinuous");
	} else if ( option[ basic::options::OptionKeys::sewing::assembly_type ].value() == "continuous" ) {
		working_assembly = AssemblyFactory::create_assembly("continuous");
	}

	//Vector of all currently added Assemblies
	utility::vector1<AssemblyOP> assembly_list;

	//Best complete assembly seen throughout the simulation
	AssemblyOP best_complete_assembly;

	//set to arbitrarily high starting value
	core::Real best_score = 10000;

	//Main loop
	core::Size starttime = time(NULL);
	core::Real random_action = numeric::random::rg().uniform();
	core::Real current_add_probability = 0.0;
	core::Real current_delete_probability = 0.0;
	//core::Real current_switch_probability = 0.0;

	core::Size cur_cycle = 0;
	core::Size total_attempts_counter = 0;
	while ( cur_cycle < cycles_ ) {
		++total_attempts_counter;

		AssemblyOP pre_op_assembly = working_assembly->clone();
		utility::vector1<AssemblyOP> pre_op_assembly_list(assembly_list.size());
		for ( core::Size i=1; i <= assembly_list.size(); ++i ) {
			pre_op_assembly_list[i] = assembly_list[i]->clone();
		}

		if ( TR.Debug.visible() ) {
			TR.Debug << "Cycle " << cur_cycle << " " << working_assembly->segments().size() << " segments " << best_score << " Best_Score" << std::endl;
		}
		if ( TR.Debug.visible() && working_assembly->segments().size() == (requirement_set_->get_max_segments()-1) ) {
			TR.Debug << "Maximum length cannot be reached by SMotif nodes alone. Check starting segment number." << std::endl;
		}

		//If we have an empty assembly, we have to add an edge.
		//If we already have a completed assembly, we have to delete. Otherwise
		//we can do either
		if ( assembly_list.size() == 0 ) {
			if ( TR.Debug.visible() ) { TR.Debug << "No Edges Short-Circuit to Add" << std::endl; }
			current_add_probability = 1.0;
		} else if ( !requirement_set_->can_be_added_to(working_assembly) ) {
			if ( TR.Debug.visible() ) { TR.Debug << "Full Edges Short-Circuit to Delete/Switch" << std::endl; }
			current_add_probability = 0.0;
		} else {
			//   current_delete_probability =((core::Real)working_assembly->segments().size() / (core::Real)(requirement_set_->get_max_segments()) ) * (1.0-switch_probability_);
			//   current_add_probability = 1.0 - (current_delete_probability + switch_probability_);
			current_add_probability = add_probability_;
			current_delete_probability = delete_probability_;
		}

		random_action = numeric::random::rg().uniform();
		if ( TR.Debug.visible() ) {
			TR.Debug << "random: " << random_action
				<< " add: " << current_add_probability
				<< " delete: "<< current_delete_probability
				<< " switch: "<< switch_probability_ << std::endl;
		}

		//Add operation
		if ( random_action < current_add_probability ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << std::endl;
				TR.Debug << "ADD EDGE" << std::endl;
			}
			if ( !add_edge(pre_op_assembly, pre_op_assembly_list, working_assembly, assembly_list) ) {
				if ( TR.Debug.visible() ) { TR.Debug << "add_edge fails" << std::endl; }
				continue;
			}

			append_movie_frame(working_assembly, cur_cycle);
		} else if ( random_action < (current_add_probability + current_delete_probability) ) {
			//Delete operation
			if ( TR.Debug.visible() ) { TR.Debug << "DELETE EDGE" << std::endl; }
			delete_edge(working_assembly, assembly_list);
			append_movie_frame(working_assembly, cur_cycle);
		} else {
			//Switch operation
			if ( TR.Debug.visible() ) { TR.Debug << "SWITCH EDGE" << std::endl; }
			if ( !switch_edge(pre_op_assembly, pre_op_assembly_list, working_assembly, assembly_list) ) {
				if ( TR.Debug.visible() ) { TR.Debug << "switch_edge fails" << std::endl; }
				continue;
			}
			boltzman(
				pre_op_assembly_list,
				assembly_list,
				pre_op_assembly,
				working_assembly,
				cur_cycle,
				best_complete_assembly,
				best_score
			);
		}
		//If we haven't failed an edge addition this cycle, increment the cycles
		++cur_cycle;
	}

	core::Size endtime = time(NULL);
	TR << "Completed " << cur_cycle << " cycles (" << total_attempts_counter << " total attempts) in " << endtime - starttime << " seconds" << std::endl;
	return best_complete_assembly;
} //generate_assembly


///@details If we haven't added an edges yet, pick a random node to start from. Otherwise
///Pick a random node from the current assembly that has the ability to have edges built
///from it. This method returns true if an edge was added, false if we couldn't add an edge
///either due to no edges being in the graph, or no edges satisfying requirements.
bool
MonteCarloAssemblyMover::add_edge(
	AssemblyOP const & pre_op_assembly,
	utility::vector1<AssemblyOP> const & pre_op_assembly_list,
	AssemblyOP & assembly, //working_assembly
	utility::vector1<AssemblyOP> & assembly_list
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( TR.Debug.visible() ) { TR << "MonteCarloAssemblyMover::add_edge" << std::endl;}

	assembly_list.push_back(assembly->clone()); // so working_assembly will be pushed back to this assembly_list
	if ( TR.Debug.visible() ) { TR.Debug << "assembly_list.size(): " << assembly_list.size() << std::endl; }
	if ( assembly_list.size() == 1 ) {
		add_starting_model(assembly);
		return true;
	} else {
		if ( TR.Debug.visible() ) {
			TR.Debug << "assembly->get_next_reference_node(graph_): " << assembly->get_next_reference_node(graph_) << std::endl;
			// so it turns out that assembly->get_next_reference_node(graph_) is not model_id!
		}

		ModelNode const * reference_node = graph_->get_model_node(assembly->get_next_reference_node(graph_));



		if ( TR.Debug.visible() ) {
			TR.Debug << "edge_file_: " << edge_file_ << std::endl;
			TR.Debug << "node_id in other word, reference_node->get_node_index(): " << reference_node->get_node_index() << std::endl;
		}

		if ( TR.Debug.visible() ) { TR.Debug << "before add_edges_from_binary " << std::endl; }

		graph_->add_edges_from_binary(edge_file_, reference_node->get_node_index());
		// the 2nd parameter reference_node->get_node_index()


		if ( TR.Debug.visible() ) { TR.Debug << "after add_edges_from_binary " << std::endl; }



		//if there are no edges from this node then this ADD_EDGE operation is a no-op
		core::Size num_edges = reference_node->num_edges();

		if ( TR.Debug.visible() ) {
			TR.Debug << "num_edges: " << num_edges << std::endl;
		}

		/// num_edges: the number of edges incident on this node, which may include a loop edge
		if ( num_edges == 0 ) {
			assembly = pre_op_assembly;
			assembly_list = pre_op_assembly_list;
			TR.Debug << "rejecting add, no edges" << std::endl;
			return false;
		}

		//Pick a random edge and follow it
		core::Size edge_num = numeric::random::random_range(0, (int)num_edges-1);
		core::graph::EdgeListConstIterator edge_it = reference_node->const_edge_list_begin();
		for ( core::Size j=1; j<=edge_num; ++j ) {
			++edge_it;
		}

		//Cast the edge to a proper HashEdge and follow it
		HashEdge const * const cur_edge = static_cast< HashEdge const * >(*edge_it);

		//If we've already added this model, don't add it again. This should theoretically
		//not be a problem, but due to the way model regeneration is handled in the Assembly
		//class, you get an error. This should be fixed in Assembly soon.
		core::Size reference_node_index = reference_node->get_node_index();
		core::Size mobile_node_index = cur_edge->get_other_ind(reference_node_index);
		core::Size mobile_model_id( graph_->get_model_node(mobile_node_index)->model().model_id_ );
		std::set<core::Size> model_ids = assembly->model_ids();
		if ( model_ids.find(mobile_model_id) != model_ids.end() ) {
			assembly = pre_op_assembly;
			assembly_list = pre_op_assembly_list;
			TR.Debug << "rejecting add, model previously added" << std::endl;
			return false;
		}

		//If we are violating requirements, then return false
		if ( requirement_set_->violates(assembly) ) {
			TR.Debug << "rejecting add, violation" << std::endl;
			assembly = pre_op_assembly;
			assembly_list = pre_op_assembly_list;
			return false;
		}

		assembly->follow_edge(graph_, cur_edge, reference_node->get_node_index());

		if ( TR.Debug.visible() ) { TR.Debug << "return true in add_edge" << std::endl; }

		return true;
	}
} //add_edge



///@details Simply go back to the Assembly before
///the most recent edge addition
void
MonteCarloAssemblyMover::delete_edge(
	AssemblyOP & assembly,
	utility::vector1<AssemblyOP> & assembly_list
) const {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( TR.Debug.visible() ) { TR << "MonteCarloAssemblyMover::delete_edge" << std::endl;}

	//If we've added edges, simply go back to the
	if ( assembly_list.size() > 0 ) {
		assembly = assembly_list.back();
		assembly_list.pop_back();
	} else {
		utility_exit_with_message("Attempt to delete a non-existant edge!");
	}
}


///@details A switch is just an auto-accepted delete followed by an add, the
///reason it needs to be its own operation is due to the fact
///that our scoring doesn't properly weight the value of an empty
///Assembly, and thus removing the first node will never be accepted.
///So, if we have only one node, we only allow adding and switching.
bool
MonteCarloAssemblyMover::switch_edge(
	AssemblyOP const & pre_op_assembly,
	utility::vector1<AssemblyOP> const & pre_op_assembly_list,
	AssemblyOP & assembly,
	utility::vector1<AssemblyOP> & assembly_list
) const {
	delete_edge(assembly, assembly_list);
	//If we fail to delete then make switch a no-op
	if ( !add_edge(pre_op_assembly, pre_op_assembly_list, assembly, assembly_list) ) {
		return false;
	}
	return true;
}


void
MonteCarloAssemblyMover::boltzman(
	utility::vector1<AssemblyOP> const & pre_op_assembly_list,
	utility::vector1<AssemblyOP> & assembly_list,
	AssemblyOP const & pre_op_assembly,
	AssemblyOP & working_assembly,
	core::Size cur_cycle,
	AssemblyOP & best_complete_assembly,
	core::Real & best_score
) const{

	core::Real const old_score = assembly_scorefxn_->score(pre_op_assembly);
	core::Real const score = assembly_scorefxn_->score(working_assembly);

	core::Real const score_delta( score - old_score );
	core::Real const temp_range = max_temperature_ - min_temperature_;
	core::Real const base_range = ( 1 - 4 / (core::Real)cycles_ );
	core::Real temperature_ = temp_range * pow((core::Real)base_range,cur_cycle) + min_temperature_;
	if ( ((core::Real)cur_cycle/cycles_) > 0.9 || temperature_ < 0.0000001 )  {
		temperature_ = 0.0000001;
	}
	core::Real const boltz_factor =  -score_delta / temperature_;
	core::Real const probability = std::exp( std::min (40.0, std::max(-40.0,boltz_factor)) );

	if ( TR.Debug.visible() ) {
		TR << "min_temperature_: " << min_temperature_ << std::endl;
		TR << "max_temperature_: " << max_temperature_ << std::endl;
		TR << "temperature_: " << temperature_ << std::endl;
	}


	// if(TR.Debug.visible()) {
	//  if (working_assembly->segments().size() > 5){
	//
	//   std::string seg_tag = utility::to_string(working_assembly->segments().size());
	//
	//   //Dump the full pose
	//   core::pose::Pose to_be_dumped_full_pose_indirct = get_fullatom_pose(working_assembly);
	//   to_be_dumped_full_pose_indirct.dump_pdb( seg_tag + "_segments_long_assembly_full_from_to_pose_indirectly_called.pdb" );
	//
	//   core::pose::Pose to_be_dumped_full_pose = working_assembly->to_pose(core::chemical::FA_STANDARD);
	//   to_be_dumped_full_pose.dump_pdb( seg_tag + "_segments_long_assembly_full_from_to_pose_directly_called.pdb" );
	//
	//   //Dump the multichain pose
	//   core::pose::Pose multi_chain_pose = working_assembly->to_multichain_pose(core::chemical::FA_STANDARD);
	//   multi_chain_pose.dump_pdb( seg_tag + "_segments_long_assembly_multichain.pdb");
	//
	//   utility_exit_with_message("just dumped working_assembly as a pdb");
	//  }
	// }

	//Accept
	if ( numeric::random::rg().uniform() < probability ) {
		if ( TR.Debug.visible() ) {
			TR.Debug << "ACCEPTING " << " " << pre_op_assembly->segments().size() << " segments to " << working_assembly->segments().size() << " segments " <<
				" ( " << old_score << " -> " << score << " ) at temp " << temperature_ << std::endl;
			append_movie_frame(working_assembly, cur_cycle);
		}

		bool check_completeness_of_assembly = false;

		if ( use_best_assembly_score_ ) { // default behavior in Tim's era
			//If this is our best completed assembly, then save it
			if ( requirement_set_->satisfies(working_assembly) && score < best_score ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "SAVING BEST" << std::endl;
				}
				check_completeness_of_assembly = true;

			}
		} else {
			if ( requirement_set_->satisfies(working_assembly) && score < min_assembly_score_ ) {
				check_completeness_of_assembly = true;
			}
		}
		if ( check_completeness_of_assembly ) {
			bool this_pose_is_complete = true;
			if ( remove_cut_off_assembly_ ) { // remove incomplete (cut-off) assembly (mostly comes from inherent error in pdb file), in the future, removing incomplete model at model extration step will be pursued
				core::pose::Pose to_be_checked_pose = get_fullatom_pose(working_assembly);

				for ( Size ii=1; ii < to_be_checked_pose.total_residue(); ii++ ) {
					core::Real distance = to_be_checked_pose.residue(ii).atom("CA").xyz().distance(to_be_checked_pose.residue(ii+1).atom("CA").xyz());
					if ( distance > 5.0 ) {
						this_pose_is_complete = false;
						if ( TR.Debug.visible() ) {
							TR.Debug << "Don't use this assembly, you will throw it away anyway unless you rebuild missing region!" << std::endl;
						}
						break;
					}
				}
			}
			if ( ((remove_cut_off_assembly_) && (this_pose_is_complete)) || (!remove_cut_off_assembly_) ) {
				if ( TR.Debug.visible() ) {
					TR.Debug << "SAVING current backbone since (assembly_score: " << score << ") < (min_assembly_score: " << min_assembly_score_ << ")" << std::endl;
				}
				if ( use_best_assembly_score_ ) { // default behavior in Tim's era
					best_complete_assembly = working_assembly->clone();
					best_score = score;
				} else {
					best_complete_assembly = working_assembly->clone();
				}
			}
		}

	} else { //Reject
		if ( TR.Debug.visible() ) {
			TR.Debug << "REJECTING " << " " <<pre_op_assembly->segments().size() << " segments to " << working_assembly->segments().size() <<
				" ( " << old_score << " -> " << score << " ) at temp " << temperature_ << std::endl;
		}

		working_assembly = pre_op_assembly;
		assembly_list = pre_op_assembly_list;
	}
} //boltzman


void
MonteCarloAssemblyMover::parse_my_tag(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
){

	parent::parse_my_tag(tag, data, filters, movers, pose);

	/////////Monte-carlo cycles////////////
	if ( tag->hasOption("cycles") ) {
		cycles_ = tag->getOption<core::Size>("cycles");
	}
	if ( tag->hasOption("add_probability") ) {
		add_probability_ = tag->getOption<core::Real>("add_probability");
	}
	if ( tag->hasOption("delete_probability") ) {
		delete_probability_ = tag->getOption<core::Real>("delete_probability");
	}
	if ( tag->hasOption("switch_probability") ) {
		switch_probability_ = tag->getOption<core::Real>("switch_probability");
	}
	if ( tag->hasOption("min_temperature") ) {
		min_temperature_ = tag->getOption<core::Real>("min_temperature");
	}
	if ( tag->hasOption("max_temperature") ) {
		max_temperature_ = tag->getOption<core::Real>("max_temperature");
	}
	if ( tag->hasOption("use_best_assembly_score") ) {
		use_best_assembly_score_ = tag->getOption<bool>("use_best_assembly_score");
	}
	if ( tag->hasOption("min_assembly_score") ) {
		min_assembly_score_ = tag->getOption<core::Real>("min_assembly_score"); // for Doonam's a/b design, -2.0 is recommended
	}

	if ( tag->hasOption("remove_cut_off_assembly") ) {
		remove_cut_off_assembly_ = tag->getOption<bool>("remove_cut_off_assembly");
	}
}

} //sewing
} //protocols
