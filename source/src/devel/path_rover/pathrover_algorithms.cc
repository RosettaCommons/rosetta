
////////////////////////////////////////////////////////////////////////////////
//     Pathwaysalgorithms.cc: module for generating pathways between conformations
//
//     see pathways.h
/////////////////////////////////////////////////////////////////////////////////


// Rosetta Headers
#include "cst_set.h" // pose constraints
#include "minimize.h"
#include "pack.h"
#include "pathways.h"
#include "pathways_planners.h"
#include "pose.h"
#include "pose_io.h"
#include "pose_rotamer_trials.h"
#include "prof.h"
#include "random_numbers.h"
#include "rotamer_trials.h"
#include "score.h"
#include "score_ns.h"

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// Utility Headers
//#include "src/utility/basic_sys_util.h"
//#include "src/utility/io/ozstream.h"

// C++ Headers
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <limits>
#include <cmath>
#include <set>
#include <sstream>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include "pathways_algorithms.h"

static int allocated_nodes = 0;

namespace pathways{

//	#define MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA 2

// **********************************************
// ********* general comparator methods *********
// **********************************************

bool partial_data_distance_less(Partial_data_node_dist n, Partial_data_node_dist m){
	return (n._distance < m._distance);
}

bool energy_path_less(Path_energy p, Path_energy q)
{
	return (p._max_energy < q._max_energy);
}


//TODO: energy
bool energy_less(RRT_node* p,RRT_node* q)
{
	return (p->get_score() < q->get_score());
}


// *******************************************
// ********* Single_tree_RRT methods *********
// *******************************************


// Takes a node 'n' that contains DOFs vector, and produces,
// optimizes and scores a pose object from the node using the
// DOFs vector. The score is stored in 'n'
//
// Optimizations used: rotamer trial + BB/SC torisonal minimization
//
// Constraints: a heavy constraint does not allow big changes
// in BB torsions
//
// params -
// @param n [I/O]: node that contains the DOFs info. The new score is stored in 'n'
//
// returns: a newly allocated pose object based on info from 'n', new score is also stored in 'n'
pose_ns::Pose* Single_tree_RRT::build_and_score_pose_from_node(RRT_node* n)
{
  using namespace std;

  bool local_debug(true);
  DOFs_manager& dofs_manager = _owner->get_dofs_manager();
  pose_ns::Pose* poseP = n->produce_pose();
  std::vector<double> old_DOFs = dofs_manager.get_dofs_values_vector(*poseP); // for DEBUG purposes

  // score with optimizations (rotamer trials + min)
  if(_owner->is_full_atom())
    {
      if(local_debug)
	{
	  poseP->score(_score_weight_map);
	  double old_score = poseP->get_0D_score(pose_ns::SCORE);
	  cout << "Score before RT+minimization: " << old_score << endl;
	}

      // rot. trial
      poseP->set_allow_chi_move(true);
      poseP->rottrial(_score_weight_map, 1 /*cycles*/); // TODO: what do the cycles mean?

      if(local_debug)
	{
	  cout << "Score after RT: " << poseP->get_0D_score(pose_ns::SCORE) << endl;
	}

      // minimize (locally for backbone, unconstrained for sidechains) // TODO: flag? for tolerance? for any minimization? make same as in local-planner?
      poseP->set_allow_bb_move(true);
      double min_tolerance = 10;
      char* min_type = "lbfgs_armijo_atol"; // DFP gradient descent + absolute tolerance
      minimize_set_tolerance( min_tolerance );
      cst_set_ns::Cst_set bb_tor_cst_set;     //Backbone torsional constraints:
      _score_weight_map.set_weight( pose_ns::PHIPSI_CST, 250.0 );
      // _scoreweight_map.set_weight( pose_ns::OMEGA_CST, 1.0 );
      for( int res = 1, nres = poseP->total_residue()  ;
	   res <= nres ; ++res )
	{
	  bb_tor_cst_set.add_rosetta_torsion_constraint( res, 1, poseP->phi( res ) );
	  bb_tor_cst_set.add_rosetta_torsion_constraint( res, 2, poseP->psi( res ) );
	  // bb_tor_cst_set.add_rosetta_torsion_constraint( res, 3, pose.omega( res ) );
	}
      poseP->set_constraints( bb_tor_cst_set );
      poseP->dump_pdb("./output/before_minimize.pdb");
      poseP->main_minimize( _score_weight_map, min_type );
      poseP->score(_score_weight_map);
      if(local_debug)
	{
	  cout << "Score " << poseP->get_0D_score(pose_ns::SCORE) << " (with PHI_PSI_CST)" << endl;
	}

      _score_weight_map.set_weight( pose_ns::PHIPSI_CST, 0.0 );
      poseP->score(_score_weight_map);
      if(local_debug)
	{
	  std::vector<double> new_DOFs = dofs_manager.get_dofs_values_vector(*poseP);
	  cout << "Delta of DOFs due to minimization:" << endl << " [";
	  for(unsigned int i = 0 ; i < old_DOFs.size() ; i++)
	    cout << (new_DOFs[i] - old_DOFs[i]) << ";";
	  cout << "]" << endl;
	  cout << endl;
	}
    } else { // centroid mode
      poseP->score(_score_weight_map);
    }


  double new_score = poseP->get_0D_score(pose_ns::SCORE); // ignore the constraint weight for final weighting
  // double new_score = poseP->score(_score_weight_map);
  if(local_debug)
    {
      if(_owner->is_full_atom())
	cout << "Score " << new_score << " (without PHI_PSI_CST)" << endl;
      cout << poseP->show_scores() << endl;
    }

  n->set_score(new_score);

  return poseP;
}


// print path [0 .. path.size()]
void Single_tree_RRT::print_path(std::vector<RRT_node* >& path,int path_num,std::string dir)
{
	using namespace std;

	ostringstream ostr_energy_fname, ostr_energy_verbose_fname;
	ostr_energy_fname << dir << path_num << "/PATH_ENERGY.txt" << ends;
	ostr_energy_verbose_fname << dir << path_num << "/PATH_ENERGY_VERBOSE.txt" << ends;
	//    std::string energy_fname = ostr_energy_fname.str();
	//cout << "Print Path " << ostr_energy_fname.str() << endl;
	ofstream fout_energy(ostr_energy_fname.str().c_str(), ios::out);
	ofstream fout_energy_verbose(ostr_energy_verbose_fname.str().c_str(), ios::out);
	// also open a verbose file
	fout_energy << "% PATH ENERGY: old score / new score (old score = score when added node to tree)" << endl;
	fout_energy_verbose << "% VERBOSE PATH ENERGY" << endl;
	//cout << "Path" << path_num << std::endl;

	for(unsigned int j = 0; j < path.size(); j++){

	  double old_score = path[j]->get_score();
	  pose_ns::Pose* p_cur_pose = build_and_score_pose_from_node(path[j]);
	  double new_score = p_cur_pose->get_0D_score(pose_ns::SCORE);
	  std::cout << "Control: p_cur_pose->score = " << new_score << std::endl;
	  fout_energy<<old_score<< " " << new_score << std::endl; // use old score?
	  path[j]->set_score(old_score); // save old score, so it would be valid for subsequent paths

	  ostringstream ostr_conf_fname;
	  ostr_conf_fname << dir << path_num << "/conf_" << j << ".pdb";
	  //std::cout<<"Print Path conf "<<ostr_conf_fname.str()<<std::endl;
	  p_cur_pose->dump_pdb(ostr_conf_fname.str());

	  // print verbose score info:
	  fout_energy_verbose << "Frame #" << j << ":" << std::endl;
	  p_cur_pose->show_scores(fout_energy_verbose);
	  fout_energy_verbose << std::endl;

	  delete p_cur_pose;
	}
	fout_energy.close();
	fout_energy_verbose.close();
}


void Single_tree_RRT::reverse_path(std::vector<RRT_node* >& path)
{
	RRT_node* tmp;
	for(unsigned int j = 0,i = path.size() -1; j < i; j++,i--){
		tmp = path[j];
		path[j] = path[i];
		path[i] = tmp;
	}
}

// stop local planner:
// 1. reached max steps of teh local planner
// 2. reached to the vicinity of the random conformation
// 3. reached to the dof values torsions of the random conformation - due to different
//    bond angles - having teh same torsions values does Not garentee the same rmsd values
// 4. energy of the new node along teh interpolation is above some value - to
//    insure valid compact conformations with no clashes
bool Single_tree_RRT::fulfilled_stop_condition_local_planner
( Linear_planner_iterator& linplanner_iter,
		RRT_node* new_n,
		RRT_node* target_n)
{
	if(linplanner_iter.step_num() > _max_steps_local_planner) return true; // #1

	// BR 5/2/08 - removed #2 since even if they're close, we want to allow at least one step
	//             so it won't be counted as a consecutive fail (and #1 takes care of stopping anyway)
	//if (new_n->compute_dist_CA_rmsd(target_n) < _close_nodes_reached_random) // #2
	//	return true;

	// #3 collision
	float score = new_n->get_score();

	if(this->is_forbidden_energy(score))
		return true;

	return false;

}


Single_tree_RRT::Single_tree_RRT(pathways::Pathways* owner)
{
	using namespace std;
	bool local_debug = true;
	if(local_debug){
		cout << "*** Single_tree_RRT::Single_tree_RRT() - start ***" << endl;
	}
	_owner = owner;
	// initialize params:
	Params_handler& ph = _owner->get_params_h();
	_max_generated_nodes = ph.max_generated_nodes;
	_max_consecutive_fails = ph.max_consecutive_fails;
	_num_steps_to_generate_nodes = ph.num_steps_to_generate_nodes;
	_torsion_angle_max_step = ph.torsion_angle_max_step;
	_max_steps_local_planner = ph.max_steps_local_planner;
	_distance_offset_nearest_neighbours_CA_RMSD =
		ph.distance_offset_nearest_neighbours_CA_RMSD;
	_distance_offset_nearest_neighbours_DOFs_vector =
		ph.distance_offset_nearest_neighbours_DOFs_vector;
	_min_level_motion = ph.min_level_motion;
	_max_output_pathways = ph.max_output_pathways;
	_max_generated_pathways = ph.max_generated_pathways;
	_close_nodes_reached_random = ph.close_nodes_reached_random;
	_max_energy = ph.max_energy;
	_is_soft_max_energy = ph.is_soft_max_energy;
	_stddev_soft_max_energy = ph.stddev_soft_max_energy;
	_score_weight_map = _owner->get_score_weight_map();

	_num_tested_nodes = 0;

	// minimize parameters // TODO: flag for minimization?
	minimize_exclude_sstype( false, false );
	minimize_set_vary_phipsi( true );
	minimize_set_vary_omega( false );
	minimize_set_vary_chi( true );
	minimize_set_vary_rb_trans( false );
	minimize_set_vary_rb_angle( false );
	minimize_set_local_min( false, 0 );
	// pose_set_use_nblist( true ); // TODO: nb_list means restriction to neighbouring residue pairs in calculations. However, if this is true, it requires frequent udpates to the list - TODO: should we use it (and update the list frequently) or not?

}

void  Single_tree_RRT::run()
{
	using namespace std;
	bool local_debug = false;//true;
	_score_weight_map = _owner->get_score_weight_map();
	//char resultFileName[100];
	double score_src = _owner->get_src_pose().score(_score_weight_map );
	if(local_debug)
		cout << "Score src = " << score_src << endl;

	RRT_node* root = new RRT_node(
			_owner->get_src_pose(),
			_owner->get_src_pose() /*template*/ ,
			_owner->get_dofs_manager(),
			true /* is_score_valid*/,
			new Node_rmsd_info(_owner->get_src_pose()));
	allocated_nodes++; // DEBUG
	// add the root to the DAG (meaningless for single algo - for compatibility)
	_conf_DAG.add_node(root, NULL);
	_conf_DAG.set_active(root);

	grow_DAG();

	// find pathways to any reasonable target
	std::vector<RRT_node*> optional_targets;
	RRT_node* best_target = find_best_conformation(optional_targets);
	std::sort(optional_targets.begin(), optional_targets.end(), energy_less);
	if(local_debug) cout << "Found " << optional_targets.size() << " optional targets" << endl;

	// dump pathways
	if (optional_targets.size() == 0)
		cout << "No plausible path found" << endl;
	for(unsigned int i = 0;
	i < optional_targets.size() && i < _max_output_pathways;
	i++)
	{
		RRT_node* target = optional_targets[i];
		std::vector<RRT_node* > path;
		if(i == 0 && target != best_target && local_debug)
			cout << "DEBUG: target != best_target" << endl;
		findPath(root, target,path); // from root to target
		if(local_debug)
			cout << "Path with " << path.size() << " nodes" << endl;
		print_path(path, i ,"./output/SingleRRT/");
	}

}


int Single_tree_RRT::fulfilled_stop_condition(int num_generated_nodes,int num_conseq_fails)
{
	static int i=0;
	if (num_generated_nodes > i){
		std::cout << "(num_conseq_fails,num_generated_nodes, _MAX_GENERATED_NODES) = "
			  << num_conseq_fails << ", "<<num_generated_nodes << ", " << _max_generated_nodes << std::endl;
		i+=1000;
	}

	if (num_generated_nodes > _max_generated_nodes)
		return 1;
	if (num_conseq_fails > _max_consecutive_fails)
		return 1;
	return 0;
}


RRT_node* Single_tree_RRT::generate_random_conformation()
{
	using namespace std;
	using namespace pose_ns;
	bool local_debug = false;//true;
	if(local_debug)
		cout << "*** Single_tree_RRT::generate_random_conformation() - start ***" << endl;
	// building new pose based on root template: // TODO: can do it much more efficiently - work directly with vector and sample on it
	Pose* p_new_pose = _conf_DAG.get_active_root()->produce_pose();
	assert(p_new_pose);
	DOFs_manager& dofs_manager = _owner->get_dofs_manager();
	dofs_manager.apply_uniform_sample_all(*p_new_pose);
	float score = p_new_pose->score(_score_weight_map);
	RRT_node* n = new RRT_node(
			*p_new_pose,
			_conf_DAG.get_active_template(),
			dofs_manager,
			true /*is_score_valid*/,
			new Node_rmsd_info(*p_new_pose));
	allocated_nodes++; // DEBUG
	if(local_debug)
	{
		cout << "Score ret_value = " << score
		<< " ; pose.get_0D_score(SCORE) = " << p_new_pose->get_0D_score(SCORE)
		<< endl;
	}
	if(p_new_pose)
	{	delete p_new_pose; p_new_pose = NULL; }
	if(local_debug)
		cout << "*** Single_tree_RRT::generate_random_conformation() - finished ***" << endl;

	return n;
}


// returns vector of nearest nodes within some offeset_distance to n in the active_tree
// stores the distance in "dist"
std::vector<RRT_node *> Single_tree_RRT::find_nearest_nodes_CA_RMSD(RRT_node* n,double& dist)
{
	using namespace std;
	int local_debug = 0;//2;
	if(local_debug >= 1) cout << "::find_nearest_nodes" << endl;
	vector<RRT_node *> results;
	//RRT_node* result = NULL;
	dist = 10000;
	//double curr_dist=0;
	vector<double> distances;
	vector<RRT_node *>::iterator it;
	vector<RRT_node *>& node_list = _conf_DAG.get_active_nodes();
	for(it = node_list.begin(), end = node_list.end();++it){
		double curr_dist = n->compute_dist_CA_rmsd(*it);
		if(local_debug >= 3) cout << "Curr CA_RMSD = " << curr_dist << endl;
		distances.push_back(curr_dist);
		if(curr_dist < dist){
			//result = *it;
			dist = curr_dist;
		}
	}
	// save close nodes
	it = node_list.begin();
	for(unsigned int i = 0; i < distances.size(); i++, ++it )
	{
		if(distances[i] <= dist + _distance_offset_nearest_neighbours_CA_RMSD)
			results.push_back(*it);
	}
	// now reduce to a maximum of 2 neighbours// TODO: parametrize this number
	while( results.size() > 2 )
	{
		unsigned del_index = (unsigned)(ran3() * results.size());
		it = results.begin();
		for(unsigned i = 0; i < del_index; i++, ++it ); // get iterator to results[del_index]
		// TODO: primitive, switch to set...?
		results.erase(it);
	}

	if(local_debug >= 1)
		cout << "find_nearest_nodes_CA_RMSD() - FINISHED with "
		<< results.size() << " neighbours" << endl;

	return results;
}

// returns vector of nearest nodes within some offeset_distance to n in the active_tree
// stores the distance in "dist"
std::vector<RRT_node *>
Single_tree_RRT::find_nearest_nodes_by_DOFs_vector(
		RRT_node* n,
		double& dist )
{
	using namespace std;
	int local_debug = 0;//2;
	if(local_debug >= 1)
		cout << "Single_tree_RRT::find_nearest_nodes_by_DOFs_vector" << endl;
	std::vector<RRT_node *> results;
	//RRT_node* result = NULL;
	dist = 10000;
	//double curr_dist=0;
	std::vector<double> distances;
	std::vector<RRT_node *>::iterator it, end;
	std::vector<RRT_node *>& node_list = _conf_DAG.get_active_nodes();
	for(it = node_list.begin(), end = node_list.end(); it!= end; ++it ) {
		double curr_dist = n->compute_dist_dofs_vector(*it); // functor according to L2 norm on internal coordinates vector
		if(local_debug >= 3)
			cout << "Curr L2_DOFs_vector = " << curr_dist << endl;
		distances.push_back(curr_dist);
		if(curr_dist < dist){
			RRT_node* result = *it;
			dist = curr_dist;
		}
	}
	it = node_list.begin();
	for(unsigned int i = 0 ;
	i < distances.size() && results.size() < 10 ; i++ ){ // TODO: switch to random version below for limiting #
		if(distances[i] <= dist + _distance_offset_nearest_neighbours_DOFs_vector)
		{
			results.push_back(*it);
			if(local_debug >= 2){
				cout << "Found near neighbour (dist = " << dist << ")" << endl;
			}
		}
		++it;
	}
	// now stay with maximum 2 // TODO: parametrize this number
	while( results.size() > 2 )
	{
		unsigned del_index = (unsigned)(ran3() * results.size());
		it = results.begin();
		for(unsigned i = 0; i < del_index; i++, ++it); // get iterator to results[del_index]
		// TODO: primitive, switch to set...?
		results.erase(it);
	}
	if(local_debug >= 1)
		cout << "find_nearest_nodes_by_DOFs_vector() - FINISHED with "
		<< results.size() << " neighbours" << endl;

	return results;
}


// new_generated_nodes gets nodes by order in which they're added
int Single_tree_RRT::add_nodes_with_local_planner(RRT_node* from,RRT_node* to,
		std::vector<RRT_node*>& new_generated_nodes)
{
	using namespace std;
	using namespace pose_ns;

	int local_debug = 2; //3;//2;
	int num_added_nodes = 0;
	DOFs_manager& dofs_manager = _owner->get_dofs_manager();
	RRT_node* new_n = NULL,
	* curr_parent = from,
	* prev_n = NULL;
	pose_ns::Pose debug_pose;
	//new_generated_nodes.clear();
	Linear_planner_iterator linplanner_dofs_iter (
			from->get_dofs_vector(),
			to->get_dofs_vector(),
			_torsion_angle_max_step );

	if(local_debug >= 2){
		std::cout << endl;
		cout << "Local planner FROM: [";
		for(unsigned int i = 0 ; i < from->get_dofs_vector().size() ; i++)
			cout << from->get_dofs_vector()[i] << ";";
		cout << "]" << endl;
		std::cout << endl;
		std::cout << endl;
		cout << "Local planner TO: [";
		for(unsigned int i = 0 ; i < to->get_dofs_vector().size() ; i++)
			cout << to->get_dofs_vector()[i] << ";";
		cout << "]" << endl;
		cout << "dist_CA_rmsd(FROM,TO) = " <<
		  from->compute_dist_CA_rmsd( to ) << " A" << endl;
		std::cout << endl;
	}


	//char resultFileName[100];
	if (fulfilled_stop_condition_local_planner(linplanner_dofs_iter, from, to) )
	{	  return 0;} // stop local planner}

	linplanner_dofs_iter++; // skip "from"

	//ANGELA DEBUG
	//Pose* pp_new_pose = to->produce_pose();
	//std::sprintf(resultFileName, "./output/local_planner_rand_conf.pdb");
	//pp_new_pose->dump_pdb(resultFileName);
	//delete pp_new_pose;

	// REPACK ARRAYS CODE:
	//  // TODO:: Barak - is it ok? the side chains are also symetric
	//  // we want to move all of them right?
	//  FArray1D_bool allow_repack(from->get_pose().total_residue(),false);
	//  for( int ii = 1; ii <= from->get_pose().total_residue(); ++ii ) {
	//    allow_repack(ii) = true;
	//  }

	bool added_in_prev_round = false;
	while( !linplanner_dofs_iter.is_end() ) // TODO: energy, close to the random node
	{
		// use DOFs from iterator to build new node pose from template
		int step = linplanner_dofs_iter.step_num();
		new_n = new RRT_node(
				*linplanner_dofs_iter,
				_conf_DAG.get_active_template(),
				dofs_manager,
				false /* is_score_valid */);
		allocated_nodes++;
		_num_tested_nodes++;
		if(local_debug >= 2){
			cout << endl;
			cout << "Local planner _v_dofs #" << step << " [";
			for(unsigned int i = 0 ; i < (*linplanner_dofs_iter).size() ; i++)
				cout << (*linplanner_dofs_iter)[i] << ";";
			cout << "]" << endl;
			if(_num_tested_nodes % 100 == 0)
				cout << "Num tested nodes: " << _num_tested_nodes << endl;
			cout << endl;
		}

		// build (minimized) new pose and update node;
		// TODO: maybe replace rotamer trials with full side-chain repack? every few cycles?

		Pose* p_new_pose = this->build_and_score_pose_from_node( new_n );
		new_n->set_custom_info(
		  new Node_rmsd_info(*p_new_pose));

		if(local_debug)
		  {
		    if(_owner->is_full_atom())
		      cout << "DEBUG: Score " << new_n->get_score() << " (without PHI_PSI_CST) + distances" << endl;
		    cout << p_new_pose->show_scores() << endl;
		    cout << "dist_CA_rmsd(new_n,FROM) = " <<
		      new_n->compute_dist_CA_rmsd( from ) << " A ; ";
		    cout << "dist_CA_rmsd(new_n,TO) = " <<
		      new_n->compute_dist_CA_rmsd( to ) << " A"
			 << endl;
		    cout << "dist_CA_rmsd(new_n,SOURCE_NODE) = " <<
		      new_n->compute_dist_CA_rmsd( _conf_DAG.get_active_root() ) << endl;

		}

		if(p_new_pose){
		  delete p_new_pose; p_new_pose = NULL;
		}

		// breakout condition
		if (fulfilled_stop_condition_local_planner (
		      linplanner_dofs_iter, new_n, to) )
		  {

		    if(new_n) {delete new_n; new_n = NULL;}
		    allocated_nodes--;
		    // add last good node to tree
		    if ( !added_in_prev_round && prev_n)
		      {
			_conf_DAG.add_node(prev_n, curr_parent);
			new_generated_nodes.push_back(prev_n);
			num_added_nodes++;
		      }
		    if(local_debug >= 2)
		      cout << "add_nodes_with_local_planner(): *** Stopping local_planner ***" << endl;
		    break; // stop local planner
		  }


		if (!added_in_prev_round && prev_n){
		  delete prev_n;
		  //allocated_nodes--;
		}
		prev_n = new_n;

		// all is good -> so add node every NUM_STEPS_TO_GENERATE_NODE steps
		if (step % _num_steps_to_generate_nodes == 0)
		  {
		    added_in_prev_round = true;
		    _conf_DAG.add_node(new_n, curr_parent);
		    // book keeping:
		    new_generated_nodes.push_back(new_n);
		    curr_parent = new_n;
		    num_added_nodes++;
		    if(local_debug >= 2)
		      cout << "add_nodes_with_local_planner(): adding node, score = "/*"after rot trial = "*/
			   << new_n->get_score() << endl;
		    if(local_debug >= 3)
		      {
			// dump PDB of every k'th = 100 nodes for each tree
			const int k = 1000;
			int active_tree = _conf_DAG.get_active_id();
			int tree_size = _conf_DAG.get_tree_size(active_tree);
			std::cout << "add_nodes_with_local_planner() - Active tree: "
				  << active_tree << " Tree size: "
				  << tree_size << std::endl;
			if(tree_size % k == 0)
			  { // only if node is going to be accepted
			    std::stringstream sstr_fname;
			    sstr_fname << "./output/node_" << active_tree << "_" << tree_size << ".pdb" << std::ends;
			    debug_pose.dump_pdb(sstr_fname.str());
			  }
		      }
		  }
		else
		  added_in_prev_round = false;

		linplanner_dofs_iter++;
	}

	return num_added_nodes;
}


void Single_tree_RRT::grow_DAG()
{
	using namespace std;
	bool local_debug = true;
	if(local_debug) cout << "*** Single_tree_RRT::grow_DAG() - start ***" << endl;
	//assert(_owner);
	//RRT_node* n=NULL;//,*neighbour;
	int num_added_nodes;
	int num_conseq_fails =0,num_generated_nodes=1;
	std::vector<RRT_node *> generated_nodes;
	std::vector<RRT_node *> neighboring_nodes;
	double dist;

	while(!fulfilled_stop_condition(num_generated_nodes, num_conseq_fails)){

		RRT_node* n = generate_random_conformation();
		//n->get_pose().dump_pdb("./output/last_rand_conf.pdb");
		neighboring_nodes.clear();
		neighboring_nodes = find_nearest_nodes_by_DOFs_vector(n, dist); // return to general nearest nodes stuff
		// apply local-planner to all nearest neighbours
		for(unsigned int i = 0; i < neighboring_nodes.size(); i++){
			generated_nodes.clear();
			num_added_nodes =
				add_nodes_with_local_planner(neighboring_nodes[i],
						n, generated_nodes);
			if(local_debug)
			{
				cout << "grow_DAG(): Finished neighbor #" << i << std::endl;
				cout << "Added " << num_added_nodes << " nodes" << endl;
				cout << "So far, allocated " << allocated_nodes << " nodes" << endl;
				cout << "grow_DAG(): generated " << num_generated_nodes
				<< " nodes, num_conseq_fails = " << num_conseq_fails << endl;
			}
			// book-keeping
			num_generated_nodes += num_added_nodes;
			if(num_added_nodes==0)
				num_conseq_fails++;
			else
				num_conseq_fails = 0;
			if (fulfilled_stop_condition(num_generated_nodes, num_conseq_fails))
			{
				cout << "fulfilled_stop_condition(num_generated_nodes"<<endl;
				break;
			}
		}
		if(n){
			delete n;
			allocated_nodes--;
		}
	}
	if(local_debug)
		cout << "Single_tree::grow_DAG() - finished" << endl;
}


// meet_goal_criterions()
//
// does node meet most basic criteria for putative targets
int Single_tree_RRT::meet_goal_criterions(RRT_node* n){
	double energy = n->get_score();
	bool good_energy = (energy <= _max_energy);
	bool good_num_frames = (n->getLevel() >= _min_level_motion); // distance from root = # of frames
	if (!good_energy || !good_num_frames)
		return 0;

	// check if better to go thru children
	// TODO: isn't it more efficient to scan from the children in first place?
	std::vector< RRT_node* >& children = n->get_children();
	for(unsigned int i = 0; i < children.size();i++){
		if (children[i]->get_score() < energy)
			return 0;
	}
	return 1;
}


// findPath()
//
// outputs a path from "ancestor" to "descendant"
// NOTE: assumes that all nodes along the path are with the same tree_id,
//       except perhaps the descendant itself // TODO: generalize this using dijkstra?
//
// Params:
//  ancestor, descendant - first and last nodes in query for directed path
//  path - [output] - the path from ancestor to descendant
//
// Returns:
//  true if succesful
int Single_tree_RRT::findPath(RRT_node* ancestor, RRT_node* descendant, std::vector<RRT_node* >& path){
	using namespace std;
	int local_debug = 1;
	assert(ancestor);
	assert(descendant);
	path.clear();
	int ancestor_tree = ancestor->get_orig_src_root_id();

	// climb the parents that match ancestor's tree_id only // TODO: change to full dijkstra?
	// (build path along the way)
	RRT_node* cur = descendant;
	while(cur != ancestor){
		path.push_back(cur);
		std::vector<RRT_node*>& cur_parents = cur->get_parents();
		bool is_good_parent = false;
		cout << "Descendant had " << cur_parents.size() << " parents" << endl;
		for(unsigned i = 0; i < cur_parents.size(); i++){
			if(cur_parents[i]->get_orig_src_root_id() == ancestor_tree)
			{
				is_good_parent = true;
				cur = cur_parents[i];
				break;
			}
		}
		if(!is_good_parent) // dead-end
		{
			if(local_debug >= 1)
				cout << "findPath() - Didn't find path!" << endl;
			path.clear();
			return false;
		}
	}
	path.push_back(ancestor);
	reverse_path(path); // from "ancestor" to "descendant"
	if(local_debug >= 1)
		cout << "findPath() - found path length = " << path.size() << endl;
	return true;
}


//   // outputs a reverse directed path from "to" to "from" // TODO: counterintuitive... reverse output...
//   //
//   // from, to - "from" and "to" nodes
//   // path - [output] the path from "to" to "from" // TODO... reverse...
//   //
//   // returns true if successful
//   int Single_tree_RRT::findPath(RRT_node* from, RRT_node* to,std::vector<RRT_node* >& path){
//     using namespace std;
//     static int debug_depth = 0;
//     static bool did_it = false;
//     if(debug_depth % 5000 == 0)
//       cout << "findPath() RECURSION - depth = " << debug_depth;
//     debug_depth++;
//     assert(from);
//     assert(to);
//     path.clear();
//     if(from == to){
//       cout << "findPath() - FOUND******* depth=" << debug_depth << endl;
//       path.push_back(to);
//       debug_depth--;
//       return 1;
//     }
//     std::vector< RRT_node* >& children = from->get_children();
//     for(unsigned i = 0; i < children.size(); i++){
//       if(!children[i]){
// 	cout << "findPath() NULL child!!! i = " << i << ", num_children = " << children.size() << endl;
// 	assert(children[i]);
//       }
//       if (findPath(children[i], to, path)){
// 	cout << "findPath() - going back, depth = " << debug_depth << endl;
// 	path.push_back(from);
// 	//	std::cout << "findPath() - END path.size()=" << path.size() << std::endl;
// 	debug_depth--;
// 	return 1;
//       }
//     }
//     //std::cout << "findPath() - BAD END path.size()=" << path.size() << std::endl;
//     debug_depth--;
//     return 0;
//   }


RRT_node * Single_tree_RRT::find_best_conformation(std::vector<RRT_node *>& v_res_nodes){

	v_res_nodes.clear();
	std::vector<RRT_node*>::iterator it, end;
	std::vector<RRT_node*>& node_list = _conf_DAG.get_active_nodes();
	RRT_node* min_energy_node= NULL;

	//int level =0;
	double energy, min_energy = 1000000;
	for(it = node_list.begin()+1, end = node_list.end(); it != end; ++it ) { // TODO: indeed begin + 1?
		if (!(meet_goal_criterions(*it)))
			continue;
		v_res_nodes.push_back(*it);
		energy = (*it)->get_score();

		if (energy < min_energy){
			min_energy = energy;
			min_energy_node = *it;
		}
		/*if (level < (*it)->getLevel()){
	level = (*it)->getLevel();
	min_energy_node = *it;
	}*/


	}
	return min_energy_node;
}


double Single_tree_RRT::energy_function_between_flex_units(pose_ns::Pose* pose){
	double score_all_units = pose->score(_score_weight_map),score_seperated_units=0,score_interface=0;

	int num_jumps=7;
	bool local_debug= false;
	std::vector<pose_ns::Jump> vj;
	for(int i=1; i <= num_jumps; i++){
		vj.push_back(pose->get_jump(i));
	}
	std::cout<<"vj="<<vj.size()<<std::endl;

	int step_size = 200;
	int dock_jump = 1;

	int const pos1( pose->fold_tree().get_jump_point()(1, dock_jump) );
	int const pos2( pose->fold_tree().get_jump_point()(2, dock_jump) );

	const FArray3D_float & Epos( pose->Eposition() );
	numeric::xyzVector_double trans_axis (
			numeric::xyzVector_double( &Epos(1,2,pos2) ) -
			numeric::xyzVector_double( &Epos(1,2,pos1) ) );

	for(unsigned int i=0; i < vj.size(); i++ , step_size+=200){
		vj[i].translation_along_axis(Epos(1,1,pos1), trans_axis, step_size);
		pose->set_jump( i+1/*dock_jump*/, vj[i] );
	}

	if (local_debug){
		char resultFileName[100];
		std::sprintf(resultFileName, "./output/after_translation_conf.pdb");
		pose->dump_pdb(resultFileName);

		exit(0);
	}


	score_seperated_units = pose->score(_score_weight_map);
	score_interface  = score_all_units - score_seperated_units;

	std::cout<<"\nEnergy all = "<<score_all_units<<" seperated = "<<score_seperated_units<<" Interface = "<<score_interface<<std::endl;


	return score_interface;

}


// compute_energy_of_path()
//
// return the energy barrier in the path
double
Single_tree_RRT::compute_energy_of_path(
		std::vector<RRT_node* >& path)
{
	double max_energy =  path[0]->get_score();//,energy=0;

	for(unsigned int i = 1; i < path.size()-1;i++){
		double energy = path[i]->get_score();
		if(energy > max_energy)
			max_energy = energy;
	}
	return max_energy;
}


// **************************************************************************
// Bi-tree
// **************************************************************************

// ctr for growing 2 RRT trees, one from source and one from target,
// connect them and find a valid pathway from src to target
biRRT::biRRT(Pathways* owner)
: Single_tree_RRT(owner)
{
	// algo params
	_num_steps_to_connect_trees = _owner->get_params_h().num_steps_to_connect_trees;
	_max_rmsd_to_connect_trees = _owner->get_params_h().max_rmsd_to_connect_trees;
	_distance_offset_nearest_neighbours_CA_RMSD_connect_trees =
		_owner->get_params_h().distance_offset_nearest_neighbours_CA_RMSD_connect_trees;
}


// run()
//
// run the bi-RRT protocol for simultaneously growing a SOURCE
// and a TARGET tree till they connect
//
void biRRT::run()
{
	using namespace std;
	int local_debug = 1; // NOTE: if >=2, pre-exit on path finding stage below
	if(local_debug)
		cout << "*** biRRT::run() - start ***" << endl;


	pose_ns::Pose& template_pose_S = _owner->get_src_pose();
	pose_ns::Pose& template_pose_T = _owner->get_trg_pose();

	//char resultFileName[100];
	// initialize source & target
	float score_S =
		template_pose_S.score(_score_weight_map);
	if(local_debug)
		cout << "Score source: " << score_S << endl;

	template_pose_S.show_scores(std::cout);

	float score_T =
		template_pose_T.score(_score_weight_map);
	if(local_debug)
		cout << "Score target: " << score_T << endl;

	template_pose_T.show_scores(std::cout);


	RRT_node* root_S =
		new RRT_node(
				_owner->get_src_pose(),
				template_pose_S,
				_owner->get_dofs_manager() ,
				true /*is_score_valid */,
				new Node_rmsd_info(_owner->get_src_pose()) );
	RRT_node* root_T =
		new RRT_node(
				_owner->get_trg_pose(),
				template_pose_T,
				_owner->get_dofs_manager(),
				true /*is_score_valid*/,
				new Node_rmsd_info(_owner->get_trg_pose()) );

	//	pose_ns::Pose* tmp_pose = root_S->produce_pose();
	//  root_S->set_score(energy_function_between_flex_units(tmp_pose));
	//  delete tmp_pose;
	//
	//  tmp_pose = root_T->produce_pose();
	//  root_T->set_score(energy_function_between_flex_units(tmp_pose));
	//  delete tmp_pose;
	//
	//  if(local_debug){
	//    cout << "Score interface source: " << root_S->get_score() << endl;
	//    cout << "Score interface target: " << root_T->get_score() << endl;
	//  }
	//  //exit(0);

	allocated_nodes += 2;
	// add roots to DAG
	int id_root_S = _conf_DAG.add_node(root_S, NULL /*parent*/);
	int id_root_T = _conf_DAG.add_node(root_T, NULL /*parent*/);
	{
		//DEBUG
		pose_ns::Pose* ppose = NULL;
		ppose = root_S->produce_pose();
		ppose->dump_pdb("./output/root_S.pdb");
		delete ppose;
		ppose = root_T->produce_pose();
		ppose->dump_pdb("./output/root_T.pdb");
		delete ppose;
		ppose = NULL;
	}
	_conf_DAG.set_active(id_root_S);
	if(local_debug)
		cout << "biRRT::run(): Src_root_id = " << id_root_S
		<< " ; Trg_root_id = " << id_root_T << endl;

	// grow the trees
	vector< RRT_node* > connecting_nodes;
	grow_DAG(connecting_nodes);
	if(connecting_nodes.size() == 0)
	{
		std::cout<<"biRRT::run(): No path was found to connect the 2 confs :("<<std::endl;
		return;
	}

	dump_pathways(connecting_nodes, root_S, root_T);

	if(local_debug)
		cout << "*** biRRT::run(): END ***" << endl;
}


// dump_pathways()
//
// Dumps the pathways between root_S and root_T in the DAG
// assuming nodes in connecting_nodes appears on the
// interface between the two trees (src tree & trg tree)
void
biRRT::dump_pathways(
		std::vector< RRT_node* > const& connecting_nodes,
		RRT_node* root_S,
		RRT_node* root_T)
{
	using namespace std;
	int local_debug = 2;
	std::vector<Path_energy> paths; // a vector of paths with computed energy for path

	if(local_debug)
		// find a path via each connecting_node
		for(unsigned i=0; i < connecting_nodes.size();i++)
		{
			if(local_debug >= 2){
				cout << "Connecting node #" << i << endl;
			}
			std::vector<RRT_node* > path;
			std::vector<RRT_node* > rev_path;
			cout << "invoke findPath to connecting_node[i] = " << connecting_nodes[i] << endl;
			bool is_ok;
			is_ok = findPath(root_S, connecting_nodes[i], path);
			cout << "is_ok = " << is_ok << endl;
			if(!is_ok) exit(1);
			if(path[0] != root_S) exit(1);
			is_ok = findPath(root_T, connecting_nodes[i], rev_path);
			cout << "is_ok = " << is_ok << endl; if(!is_ok) exit(1);
			if(rev_path[0] != root_T) exit(1);
			reverse_path(rev_path); // go from connection to T
			if(*(--rev_path.end()) != root_T)  exit(1);
			for(unsigned int j = 1; j < rev_path.size(); j++) // skip connecting node (already in path)
			{
				path.push_back(rev_path[j]);
			}
			assert(*(--path.end()) == root_T);
			if (local_debug >= 2)
			{
				std::cout<<"Reverse-path size = " << rev_path.size()<<std::endl;
				std::cout<<"Compute energy of path and push it..."<<std::endl;
			}
			paths.push_back(
					Path_energy(path, compute_energy_of_path(path)));

		} // for i

	std::cout<<":)"<<std::endl;
	//sort according to energy buriers
	std::sort(paths.begin(), paths.end(), energy_path_less);
	std::cout<<":)"<<std::endl;


	for(unsigned int i = 0; i < paths.size() && i < _max_output_pathways; i++){
		print_path(paths[i]._path,i,"./output/BiRRT/");
	}
}


// grow_DAG()
//
// TODO: barak - maybe each node should contain a parent list too i.e. each edge = bidirected edge
//
// General:
// grow_DAG() grows trees by stochastically adding nodes to existing nodes, until a given stop
// condition is fulfilled.
//
// Specifics for biRRT:
// grow_DAG() grows a src and trg node, tries to connect them and stops when ...(?) // TODO: complete this
//
// params:
// connecting_nodes - // TODO: write this
//
void biRRT::grow_DAG(std::vector<RRT_node*>& connecting_nodes)
{
	using namespace std;
	bool local_debug = true;
	if(local_debug)
		cout << "biRRT::grow_DAG() - start" << endl;
	assert(_owner);
	RRT_node* n = NULL;//,*neighbour;
	int num_added_nodes;
	int num_conseq_fails = 0,
	num_generated_nodes = 2;
	//std::vector<RRT_node *> generated_nodes;
	std::vector<std::vector<RRT_node *> > generated_nodes_from_all_trees;
	std::vector<RRT_node *> neighboring_nodes;
	double dist=0;
	unsigned int i;
	int num_iters=1;

	generated_nodes_from_all_trees.resize(_conf_DAG.get_num_roots());

	while(!fulfilled_stop_condition(num_generated_nodes,num_conseq_fails) &&
			connecting_nodes.size() < _max_generated_pathways)
	{

		// one cycle -> go over all trees
		_conf_DAG.initialize_active();
		while(_conf_DAG.set_next_active()){
			// starting from first tree till last
			n = generate_random_conformation();
			neighboring_nodes = find_nearest_nodes_by_DOFs_vector(n, dist); // return to general nearest nodes stuff
			//generated_nodes.clear(); // per one tree

			for(i = 0; i < neighboring_nodes.size(); i++){
				num_added_nodes = add_nodes_with_local_planner(neighboring_nodes[i],
						n,
						generated_nodes_from_all_trees[_conf_DAG.get_active_id()]);
				//generated_nodes);
				num_generated_nodes += num_added_nodes;
				if(num_added_nodes==0)
					num_conseq_fails++;
				else
					num_conseq_fails = 0;
				if(local_debug)
					cout << "grow_DAG(): generated " << num_generated_nodes
					<< " nodes, num_conseq_fails = " << num_conseq_fails << endl;
			}
			//generated_nodes_from_all_trees.push_back(generated_nodes);
		}
		// try to connect trees every few cycles:
		if (num_iters % _num_steps_to_connect_trees == 0){
			try_to_connect_trees(generated_nodes_from_all_trees /*input/output*/,
					connecting_nodes/*output*/, num_generated_nodes/*output*/);//TODO
			if(local_debug)
				cout << "After trying to connect the trees - apply again the main procedure"<<endl;
		}
		num_iters++;
	} // while (fulfilling stop condition)

	if(local_debug){
		cout << "biTree::grow_DAG(): *** finished ***" << endl;
		cout << "                    -> tested " << _num_tested_nodes << " nodes" << endl;
	}

}


// find_nearest_nodes_by_DOFs_vector()
//
// Synopsis:
// go over all newly generated nodes,
// and find which one gets closest to the other tree
//
// Params:
// <generated_nodes_from_all_trees> - [input] - a list of new nodes from each tree
// <output_pairs> - put all close pairs in <output_pairs>
void biRRT::find_nearest_nodes_between_trees(
		std::vector<std::vector<RRT_node *> >& generated_nodes_from_all_trees,
		std::vector< Node_pair >& output_pairs)
{
	using namespace std;
	unsigned num_trees =
		generated_nodes_from_all_trees.size();
	t_map_node_pairs_by_dist map_dist2pairs_CA;
	t_map_node_pairs_by_dist map_dist2pairs_DOFs_vec;
	vector<RRT_node *> neighboring_nodes_CA;
	//vector<RRT_node *> neighboring_nodes_DOFs_vector;
	RRT_node *n=NULL;
	double dist=0;
	int local_debug=1;

	if (local_debug > 0)
		std::cout<<"Start::find_nearest_nodes_between trees ..."<<std::endl;

	// go over all trees
	for( unsigned curr_tree=0;
	curr_tree < num_trees;
	curr_tree++){
		// go over all the new generated nodes for each tree
		unsigned num_nodes = generated_nodes_from_all_trees[curr_tree].size();
		for( unsigned i_node = 0 ;
		i_node < num_nodes ;
		i_node++)
		{
			n = generated_nodes_from_all_trees[curr_tree][i_node];
			assert(n);
			// compare n against all trees (except for its own tree)
			for( unsigned i_tree=0;
			i_tree < num_trees;
			i_tree++)
			{
				if (i_tree == curr_tree)
					continue; // connect i_node to other trees
				if (local_debug >= 2)
					cout<<"curr_tree " << curr_tree << " vs. i_tree " <<i_tree << endl;
				if(local_debug == 0)
					_conf_DAG.print();
				// calc by the two criteria: CA & DOFs_vector
				_conf_DAG.set_active(i_tree); // (test tree #i_tree against n)
				vector< RRT_node* >::iterator it_neighb, end;
				neighboring_nodes_CA = find_nearest_nodes_CA_RMSD(n, dist); // TODO: currenly works on active tree, perhaps this should be a function parameter
				if(local_debug >= 2)
					cout << "Nearest dist RMSD CA = " << dist << endl;
				for(it_neighb = neighboring_nodes_CA.begin(), end = neighboring_nodes_CA.end();
				it_neighb != end; ++it_neighb )
				{
					double dist_CA = n->compute_dist_CA_rmsd(*it_neighb);
					map_dist2pairs_CA.insert(
							std::make_pair( dist_CA, Node_pair(*it_neighb, n) ) );
				}
				// 	      neighboring_nodes_DOFs_vector = find_nearest_nodes_by_DOFs_vector(n, dist);
				// 	      if(local_debug >= 2)
				// 		cout << "Nearest dist by DOFs_vector = " << dist << endl;
				// 	      for(it_neighb = neighboring_nodes_DOFs_vector.begin();
				// 		  it_neighb != neighboring_nodes_DOFs_vector.end();
				// 		  it_neighb++)
				// 		{
				// 		  double dist_DOFs = n->compute_dist_dofs_vector(*it_neighb);
				// 		  map_dist2pairs_DOFs_vec.insert(
				// 		    std::make_pair( dist_DOFs, Node_pair(*it_neighb, n) ) );
				// 		}
			} // for i_tree
		} // for i_node
	} // for curr_tree

	// now save only pairs s.t. dist <= [minimal distance + _distance_offset_to_refer_as_neighbours]
	double min_dist_CA = map_dist2pairs_CA.begin()->first;
	//    double min_dist_DOFs = map_dist2pairs_DOFs_vec.begin()->first;
	if(local_debug >= 1)
	{
		std::cout << "min_dist_CA = " << min_dist_CA << std::endl;
		//std::cout << "min_dist_DOFs = " << min_dist_DOFs << std::endl;
	}
	double thresh_CA =
		min_dist_CA +
		_distance_offset_nearest_neighbours_CA_RMSD_connect_trees;
	//    double thresh_DOFs =
	//min_dist_DOFs + _distance_offset_nearest_neighbours_DOFs_vector;
	cout << "dist thresholds = (" << thresh_CA << /*"," << thresh_DOFs <<*/ ")" << std::endl;
	output_pairs.clear();
	t_map_node_pairs_by_dist::const_iterator
	it_pairs, it_last;
	it_last = map_dist2pairs_CA.upper_bound(thresh_CA);
	for(it_pairs = map_dist2pairs_CA.begin();
	it_pairs != it_last ; // TODO: more efficient
	++it_pairs )
	{
		output_pairs.push_back(it_pairs->second);
	}
	//     it_last = map_dist2pairs_DOFs_vec.upper_bound(thresh_DOFs);
	//     for(it_pairs = map_dist2pairs_DOFs_vec.begin();
	// 	it_pairs != it_last ;
	// 	it_pairs++)
	//       {
	// 	output_pairs.push_back(it_pairs->second);
	//       }

	// now reduce to a maximum of 5 neighbours (but always include the best...) // TODO: parametrize this number
	while( output_pairs.size() > 5 )
	{
		unsigned del_index = (unsigned)(ran3() * output_pairs.size());
		if(del_index == 0) continue; // always skip best
		std::vector<Node_pair>::iterator it;
		it = output_pairs.begin();
		for(unsigned i = 0; i < del_index; i++, ++it); // get iterator to results[del_index]
		// TODO: primitive, switch to set...?
		output_pairs.erase(it);
	}


	if (local_debug >= 1)
		std::cout<<"END::find_nearest_nodes_between_trees " << output_pairs.size() << " final pairs" << std::endl;

}


// try_to_connect_trees()
//
// params:
// <generated_nodes_from_all_trees> - [input] - a list of new nodes from each tree (e.g. from previous runs)
//                                    [output] - only the nodes that were created in this run	//
// <connecting nodes> - [output] insertion of new nodes that connect trees if any
void biRRT::try_to_connect_trees(
		std::vector<std::vector<RRT_node *> >& generated_nodes_from_all_trees,
		std::vector<RRT_node*>& connecting_nodes,
		int& num_generated_nodes)
{
	using namespace std;
	vector< std::vector<RRT_node *> >
	new_nodes( _conf_DAG.get_num_roots() );
	vector< Node_pair > v_close_pairs; // closest inter-tree pairs
	int local_debug=1;

	if(local_debug >= 1)
		cout << endl << "*** try_to_connect_trees(): START ***" << endl;

	find_nearest_nodes_between_trees(
			generated_nodes_from_all_trees /*input*/,
			v_close_pairs /*output*/ );

	// now go over all these pairs and apply local planner
	for(unsigned i = 0; i < v_close_pairs.size(); i++)
	{
		if(local_debug >= 1){
			cout << "try_to_connect_trees(): trying pair " << i << "/" << v_close_pairs.size() << std::endl;
		}
		Node_pair& cur_pair = v_close_pairs[i];
		int n1_tree_id = cur_pair.n1->get_orig_src_root_id() ;
		vector<RRT_node *> cycle_new_nodes; // new nodes in one cycle of local planner
		int num_added_nodes;

		_conf_DAG.set_active( n1_tree_id ); // local planners works on active tree! // TODO: maybe should be in cmd-line of local planner?
		num_added_nodes = add_nodes_with_local_planner(	// grow "n1" towards "n2"
				cur_pair.n1,
				cur_pair.n2,
				cycle_new_nodes);
		// keep track of all new nodes of this run:
		num_generated_nodes += num_added_nodes;
		cout << "try_to_connect_trees(): Added " << num_added_nodes << " nodes with local planner" << endl;
		cout << "new_nodes[" << n1_tree_id << "].size() BEFORE = " << new_nodes[n1_tree_id].size() << endl;
		new_nodes[ n1_tree_id ].insert(
				new_nodes[ n1_tree_id ].end(),
				cycle_new_nodes.begin(),
				cycle_new_nodes.end() );
		cout << "new_nodes[" << n1_tree_id << "].size() AFTER = " << new_nodes[n1_tree_id].size() << endl;

		// check if n2 is close enough to any of the new nodes coming from n1 towards it, including n1
		// TODO: how do we settle the difference between DOFs and RMSD measures?
		RRT_node* closest_node = cur_pair.n1; // a-priori
		double min_dist_CA = cur_pair.n2->compute_dist_CA_rmsd( cur_pair.n1 );
		cout << "start with <cur_pair.n1> - dist_CA = " << min_dist_CA << endl;
		vector<RRT_node *>::const_iterator node_it, end;
		static bool did_it = false; static int fname_i = 0; // DEBUG
		for(node_it = cycle_new_nodes.begin(); end = cycle_new_nodes.end(); node_it != end; ++node_it)
		{
			double dist_CA =
				cur_pair.n2->compute_dist_CA_rmsd( *node_it );
			if(!did_it) // DEBUG
			{
				pose_ns::Pose* p = (*node_it)->produce_pose();
				std::ostringstream fname;
				fname << "./output/test_" << fname_i++ << "_" << dist_CA << ".pdb" << std::ends;
				p->dump_pdb(fname.str());
				fname_i++;
				delete p;
			}
			cout << "next dist_CA = " << dist_CA << endl;
			if(dist_CA < min_dist_CA)
			{
				min_dist_CA = dist_CA;
				closest_node = *node_it;
			}
		}
		did_it = true; // DEBUG
		if(local_debug >= 1)
			std::cout << "[try_connect] 'min_dist_CA'(RMSD) after local-planner from cur_pair.n2 = "
			<< min_dist_CA << std::endl;
		if(min_dist_CA < _max_rmsd_to_connect_trees)
		{   // connect two trees - line n2 to last generated node
			std::cout<< "Connecting :) *****************************************" <<std::endl;
			_conf_DAG.add_node(cur_pair.n2, closest_node /*parent*/);
			std::cout << "Connecting details: (n2,closest_node) = " << cur_pair.n2
			<< "," << closest_node << endl;
			//TODO: in case of many trees - I am not sure
			// union the two trees to one in _conf_DAG
			// need to split each edge to two
			connecting_nodes.push_back(cur_pair.n2);
		}

		if (fulfilled_stop_condition(num_generated_nodes, 0) ||
				connecting_nodes.size() > _max_output_pathways)
		{
			generated_nodes_from_all_trees = new_nodes; // output, for next round
			if (local_debug > 0)
				std::cout<<"END::Try to connect trees - stopped..."<<std::endl;
			return;
		}
	} // for i

	if (local_debug > 0)
	{
		generated_nodes_from_all_trees = new_nodes; // output, for next round
		std::cout<<"*** try_to_connect_trees(): END ***"<<std::endl;
	}
}


// void biRRT::try_to_connect_trees(std::vector<std::vector<RRT_node *> >& generated_nodes_from_all_trees,
//   std::vector<RRT_node*>& connecting_nodes,int& num_generated_nodes){
//   unsigned int curr_tree=0, i_trees;//,try_to_expand_node,i;
//   std::vector<RRT_node *> neighboring_nodes;
//   std::vector<RRT_node *> neighboring_nodes_dofs_vector;
//   std::vector<RRT_node *> generated_nodes;
//   unsigned int i_node;
//   RRT_node *n=NULL;
//   double dist=0;
//   int local_debug=1;
//   int num_added_nodes;

//   if (local_debug > 0)
//     std::cout<<"Start::Try to connect trees ..."<<std::endl;
//   for(curr_tree=0; curr_tree < generated_nodes_from_all_trees.size(); curr_tree++){
//     // gp over all trees
//     for(i_node=0; i_node < generated_nodes_from_all_trees[curr_tree].size();
// 	 i_node++)
//       {
// 	// go over all the new generated nodes for each tree
// 	n = generated_nodes_from_all_trees[curr_tree][i_node];
// 	assert(n);
// 	// this node is now regarded as the random node that
// 	// we want to move towards it
// 	for(i_trees=0; i_trees < generated_nodes_from_all_trees.size();
// 	    i_trees++)
// 	  {
// 	    if (local_debug > 0)
// 	      std::cout<<"Try to connect tree " << curr_tree << " to tree "<<i_trees<<std::endl;
// 	    if(local_debug == 0)
// 	      _conf_DAG.print();
// 	    // go over all trees - except the tree that n belongs to
// 	    if (i_trees == curr_tree)
// 	      continue; // connect i_node to other trees
// 	    _conf_DAG.set_active(i_trees); // grow tree i_tree towards n
// 	    // by the two criteria...
// 	    neighboring_nodes = find_nearest_nodes_CA_RMSD(n, dist);
// 	    neighboring_nodes_dofs_vector = find_nearest_nodes_by_dofs_vector(n, dist);
// 	    neighboring_nodes.insert(
// 	      neighboring_nodes.end(),
// 	      neighboring_nodes_dofs_vector.begin(),
// 	      neighboring_nodes_dofs_vector.end());
// 	    if (local_debug > 0){
// 	      std::cout<<"Num neighbors =  "<<neighboring_nodes.size()<<std::endl;
// 	      std::cout<<"i_trees =  ("<<i_trees<<","<<generated_nodes_from_all_trees.size()<<")"<<std::endl;
// 	      std::cout<<"i_node =  ("<<i_node<<","<<generated_nodes_from_all_trees[curr_tree].size()<<") "<<curr_tree<<std::endl;
// 	    }
// 	    std::vector< RRT_node* >::iterator it_neighb;
// 	    for(it_neighb = neighboring_nodes.begin(); it_neighb != neighboring_nodes.end(); it_neighb++)
// 	      {
// 		generated_nodes.clear();
// 		num_added_nodes = add_nodes_with_local_planner(*it_neighb, n,
// 		  generated_nodes);
// 		num_generated_nodes += num_added_nodes;

// 		if (generated_nodes.size() > 0)
// 		  {
// 		    double dist = n->compute_dist_CA_rmsd(
// 		      generated_nodes[generated_nodes.size()-1]);
// 		    if(local_debug > 0)
// 		      std::cout << "[try_connect] 'dist' (after local-planner) from cur_neighbour = "
// 				<< dist << std::endl;
// 		    if(dist < _max_rmsd_to_connect_trees)
// 		      {   // connect two trees
// 			std::cout<< "Connecting :) *****************************************" <<std::endl;
// 			RRT_node* parent = generated_nodes[generated_nodes.size()-1];
// 			_conf_DAG.add_node(n, parent);
// 			//TODO: in case of many trees - I am not sure
// 			// union the two trees to one in _conf_DAG
// 			// need to split each edge to two
// 			connecting_nodes.push_back(n);
// 		      }
// 		  }
// 		if (fulfilled_stop_condition(num_generated_nodes,0) ||
// 		  connecting_nodes.size() > _max_output_pathways)
// 		  {
// 		    if (local_debug > 0)
// 		      std::cout<<"END::Try to connect trees ..."<<std::endl;
// 		    return;
// 		  }
// 	      }
// 	  }
//       }
//   }
//   if (local_debug > 0)
//     std::cout<<"END::Try to connect trees ..."<<std::endl;
// }

// ------------------------------------------------------------------------
// ------------------ Towards_partial_data_RRT ----------------------------
// ------------------------------------------------------------------------

// in addition to the usual checks (see single-tree) check also if you
// move towards the partial data, if it is not the case for
// X consecutive steps - stop the local planner
// X = MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA

bool Towards_partial_data_RRT::fulfilled_stop_condition_local_planner
( Linear_planner_iterator& linplanner_iter,
		RRT_node* new_n,
		RRT_node* target_n)
{
  bool local_debug = false;
  static int num_local_planners = 0;
  bool flag_res = false;
  static bool flag_init = true;
  static int num_steps_move_away_from_partial_data=0;
  static double pre_dist=0; // previouse distance from partial data info
  double curr_dist = measure_distance_from_partial_data(new_n); // curr dist
  //std::cout<<"------------------------- pre_distance = "<<pre_dist<<" curr_dist = "<<curr_dist<<std::endl;
  //std::cout<<"------------------------- "<<num_steps_move_away_from_partial_data<<std::endl;
  if(flag_init){ // first iteration
    pre_dist = curr_dist;
    flag_init = false;
    num_local_planners++;
    if (local_debug){
      std::cout<<"********************* INIT"<<std::endl;
    }
  }
  else{
    if (std::fabs(pre_dist - curr_dist) < 0.001){
      num_steps_move_away_from_partial_data = 0;
    }
    else if (pre_dist < curr_dist){ // i.e., moving away from target
      if (local_debug){
	std::cout<<"Pre = "<<pre_dist<<" Curr = "<<curr_dist<<std::endl;
      }
      num_steps_move_away_from_partial_data++; // count num of 'bad' consecutive steps
    }
    else{
      num_steps_move_away_from_partial_data = 0; // productive move - intialize bad moves
      //std::cout << "Partial data fulfilled_stop_condiotion(): improved from " << pre_dist << " to " << curr_dist << std::endl;
    }
    pre_dist = curr_dist;
  }
  // every second local-planner iteration, invoke partial data condition
  if ((num_local_planners % 2 == 0) &&
    num_steps_move_away_from_partial_data >=  _max_steps_moving_away_from_partial_data){
    if (local_debug){
      std::cout<<"Stop local planner due to max steps moving away from partial data**"<<std::endl;
    }
    flag_res = true;
  }

  if(linplanner_iter.step_num() > _max_steps_local_planner){
    flag_res = true;  // #1
    if (local_debug)
      std::cout<<"Stop local planner due to linplanner_iter.step_num()"<<std::endl;
  }

  if (new_n->compute_dist_CA_rmsd(target_n) < _close_nodes_reached_random){ // #2
    flag_res = true;
    if (local_debug)
      std::cout<<"Stop local planner due to _close_nodes_reached_random = "<<_close_nodes_reached_random<<" "<<new_n->compute_dist_CA_rmsd(target_n)<<std::endl;
  }

  float score = new_n->get_score();
  if(is_forbidden_energy(score))
    flag_res = true;

  // reset flags for next time
  if  (flag_res){
    num_steps_move_away_from_partial_data=0;
    pre_dist=0;
    flag_init = true;
  }
  if (local_debug)
    std::cout<<"Flg res = "<<flag_res<<" Stop local planner - due to fulfilled stopping criterion"<<std::endl;

  return flag_res;

}


Towards_partial_data_RRT::Towards_partial_data_RRT(Pathways* owner)
: Single_tree_RRT(owner)
{_partial_data = owner->get_partial_data();
  _min_partial_target =100;
  Params_handler& ph = _owner->get_params_h();
  _max_steps_moving_away_from_partial_data = ph.max_steps_moving_away_from_partial_data;
}

void Towards_partial_data_RRT::run()
{
  using namespace std;
  bool local_debug = true;
  //char resultFileName[100];
  //double score_src = _owner->get_src_pose().score(_score_weight_map );
  if(local_debug)
    cout << "Towards_partial_data_RRT::Score src = " << score_src << endl;

  RRT_node* root = new RRT_node(
    _owner->get_src_pose(),
    _owner->get_src_pose() /*template*/ ,
    _owner->get_dofs_manager(),
    true /* is_score_valid*/,
    new Node_rmsd_info(_owner->get_src_pose()) );
  allocated_nodes++; // DEBUG
  // add the root to the DAG (meaningless for single algo - for compatibility)
  _conf_DAG.add_node(root, NULL);
  _conf_DAG.set_active(root);

  grow_DAG();

  std::cout<<":)"<<std::endl;

  // find pathways to any reasonable target
  std::vector<RRT_node*> optional_targets;
  find_close_to_partial_data_conformations(optional_targets);
  //TODO: we want to sort the paths according to the energy barriers or
  // according to the energy of the last conformation???????????????
  // current - according to energy of the target
  //std::sort(optional_targets.begin(), optional_targets.end(), energy_less); // already sorted

  std::vector<RRT_node* > path;
  std::vector<Path_energy> paths;

  std::cout<<":)"<<std::endl;

  if(local_debug)
    cout << "Found " << optional_targets.size() << " optional targets" << endl;

  for(unsigned int i=0; i < optional_targets.size();i++)
    {
      path.clear();
      findPath(root, optional_targets[i], path); // already finds a path from root to target
      //reverse_path(path);
      std::cout<<"Compute energy of path..."<<std::endl;
      paths.push_back(Path_energy(path, compute_energy_of_path(path)));
    } // for i


  for(unsigned int i = 0; i < paths.size() && i < _max_output_pathways; i++){
    int path_size = paths[i]._path.size();
    std::cout<<"Path "<<i<<" - distance from target "<<measure_distance_from_partial_data(paths[i]._path[path_size-1]) << std::endl;
    print_path(paths[i]._path,i,"./output/Partial_Data_RRT/");
  }


}


// find the optional targets sorted by their closenest to the partial data
void Towards_partial_data_RRT::find_close_to_partial_data_conformations(std::vector<RRT_node*>& optional_targets)
{
	//double curr_dist;
	//std::vector<RRT_node *>::iterator it;
	//std::cout<<"find_close_to_partial_data_conformations..."<<std::endl;
	std::vector<RRT_node *>& node_list = _conf_DAG.get_active_nodes();
	std::vector<Partial_data_node_dist> node_dist_list;
	unsigned int i;
	for(i = 0; i <  node_list.size(); i++){
		//TODO - eliminate nodes that are too far from the partial data
		// or maybe it is better to have some results instead of none????
		node_dist_list.push_back(Partial_data_node_dist(node_list[i],
				measure_distance_from_partial_data(node_list[i])));
	}
	//std::cout<<"find_close_to_partial_data_conformations...2"<<std::endl;
	std::sort(
			node_dist_list.begin(), node_dist_list.end(),
			partial_data_distance_less);
	optional_targets.clear();
	for(i = 0; i <  node_dist_list.size() && i < _max_output_pathways; i++){
		optional_targets.push_back(node_dist_list[i]._n);
		std::cout<<"Path "<<i<<" - distance from target "<<node_dist_list[i]._distance<<std::endl;
	}
	//std::cout<<"find_close_to_partial_data_conformations...3"<<std::endl;
}


double Towards_partial_data_RRT::measure_distance_from_partial_data(RRT_node* new_n)
{
  pose_ns::Pose* p_cur_pose = new_n->produce_pose();
  p_cur_pose->score(_score_weight_map);
  double res = _partial_data.compute_match(*p_cur_pose);
  delete p_cur_pose;
  if (res < _min_partial_target){
    std::cout<<"------------------------------------------ curr measure from target = "<<res<<std::endl;
    _min_partial_target = res;

    }
  return res;

}


// returns vector of nearest nodes within some offeset_distance to n in the active_tree
// stores the distance in "dist"
std::vector<RRT_node *>
Towards_partial_data_RRT::find_nearest_nodes_by_DOFs_vector(
		RRT_node* n,
		double& dist )
{
	using namespace std;
	int local_debug = 0;
	if(local_debug >= 1)
		cout << "Single_tree_RRT::find_nearest_nodes_by_DOFs_vector" << endl;
	std::vector<RRT_node *> results;
	//RRT_node* result = NULL;
	dist = 10000;
	std::vector<double> distances;
	std::vector<RRT_node *>::iterator it;
	std::vector<RRT_node *>& node_list = _conf_DAG.get_active_nodes();
	static int from_index_to_check =0;
	int curr_index=0;

	static double partial_target_dist =100000;
	double partial_target_dist_tmp;
	static RRT_node* result_partial_data = NULL;

	for(it = node_list.begin(), end = node_list.end(); it != end; ++it, ++curr_index ){
	  double curr_dist = n->compute_dist_dofs_vector(*it); // functor according to L2 norm on internal coordinates vector
	  if (curr_index >= from_index_to_check){
	    partial_target_dist_tmp = measure_distance_from_partial_data(*it);
	    if(partial_target_dist_tmp < partial_target_dist){
	      result_partial_data = *it;
	      partial_target_dist = partial_target_dist_tmp;
	    }
	  }

	  if(local_debug >= 3)
	    cout << "Curr L2_DOFs_vector = " << curr_dist << endl;
	  distances.push_back(curr_dist);
	  if(curr_dist < dist){
	    //result = *it;
	    dist = curr_dist;
	  }
	}

	from_index_to_check = curr_index;

	if (result_partial_data)
	  results.push_back(result_partial_data);

	it = node_list.begin();
	for(unsigned int i = 0 ;
	    i < distances.size() && results.size() < 10 ; i++ ){ // TODO: switch to random version below for limiting #
	  if(distances[i] <= dist + _distance_offset_nearest_neighbours_DOFs_vector)
	    {
	      results.push_back(*it);
	      if(local_debug >= 2){
		cout << "Found near neighbour (dist = " << dist << ")" << endl;
	      }
		}
	  ++it;
	}
	// now stay with maximum 2 // TODO: parametrize this number
	while( results.size() > 2 )
	  {
	    unsigned del_index = (unsigned)(ran3() * results.size());
	    it = results.begin();
	    for(unsigned i = 0; i < del_index; i++, ++it ); // get iterator to results[del_index]
	    // TODO: primitive, switch to set...?
	    results.erase(it);
	  }

	if(local_debug >= 1)
	  cout << "find_nearest_nodes_by_DOFs_vector() - FINISHED with "
		<< results.size() << " neighbours" << endl;

	return results;
}

/*****************************************************************************/
/******************************* private methods *****************************/
/*****************************************************************************/

// determines if an energy score should be considered in C_Feasible
// (this depends on parameters - '_max_energy', '_is_soft_max_energy' and '_stddev_soft_max_energy'.
//  The later leads to a deterministic choice)
bool
Single_tree_RRT::is_forbidden_energy(float score)
{
	double cur_max_energy;
	if(!_is_soft_max_energy)
		cur_max_energy = _max_energy;
	else
		cur_max_energy = _max_energy + gaussian() * _stddev_soft_max_energy;
	return score > cur_max_energy;
}


} // namespace pathways


