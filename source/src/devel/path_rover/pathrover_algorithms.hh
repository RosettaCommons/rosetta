////////////////////////////////////////////////////////////////////////////////
//     pathways_algorithms.h:
//     This class is for algorithms for generating conformation trees, especially
//     using RRT based algorithms, but with flexibility to allow other algorithms
//     and modifications.
//     The idea is to use DOFs_manager objects as a wrapper for maniplulating pose
//     DOFs, so the same algorithm can manipulate different type of DOFs
//
//	   Author: Barak Raveh & Angela Enosh
//     Created: 03/12/2007
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_algorithms_hh
#define INCLUDED_devel_path_rover_pathrover_algorithms_hh


// ObjexxFCL Headers
#include <ObjexxFCL/ObjexxFCL.hh>
#include "pathways.h"
// #include "pathways_DOFs_manager.h"
#include "pathways_planners.h"
#include "pathways_RRT_conformation_DAG.h"
#include "pathways_partial_data_manager.h"
#include "pose.h"
#include "score_data.h"
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <map>
#include <ostream>
#include <cstdlib>
#include <cstdio>

#define MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA 2


namespace pathways {

class Pathways;

class Path_energy{
public:
  double _max_energy;
  std::vector<RRT_node* > _path;
  Path_energy(std::vector<RRT_node* >& path,double max_energy):
  _max_energy( max_energy ),
  _path( path )
  {
  }
};

class Partial_data_node_dist
{
public:
  RRT_node* _n;
  double _distance;
  Partial_data_node_dist(RRT_node* n,double distance){
    _distance = distance;
    _n = n;
  }
};


// ordered pair
 struct Node_pair{
   RRT_node* n1;
   RRT_node* n2;

   // ctr:  ordered construction
   Node_pair(RRT_node* n1_, RRT_node* n2_)
     { n1 = n1_; n2 = n2_; }
 };

 typedef std::multimap< double, Node_pair > t_map_node_pairs_by_dist;


class Single_tree_RRT
{
protected:

	Pathways* _owner;
        RRT_conformation_DAG _conf_DAG;
	//int _active_root_id;
	pose_ns::Score_weight_map _score_weight_map; // current weight map

	// algorithm parameters
	int _max_generated_nodes; // maximum total number of nodes in all trees
	int _max_consecutive_fails; // how many consecutive local-planner failues are allowed
	int _num_steps_to_generate_nodes; // for local-planner
	double _torsion_angle_max_step; // maximal local-planner perturbation at any given angle
	int _max_steps_local_planner; // for one run of local planner
	double _distance_offset_nearest_neighbours_CA_RMSD; // for find_nearest_node_CA_RMSD()
	double _distance_offset_nearest_neighbours_DOFs_vector; // for find_nearest_node_DOFs_vector()
	int _min_level_motion;
	unsigned int _max_output_pathways;
	unsigned int _max_generated_pathways;
	double _close_nodes_reached_random;
	double _max_energy; // max energy for inserting nodes to a tree (i.e. C-feasible), see soft_max_energy in pathways_parameters.h
	bool _is_soft_max_energy; // see pathways_parameters.h
	double _stddev_soft_max_energy; // see pathways_parameters.h
	unsigned _num_tested_nodes; // number of nodes that were tested for
			       // addition thruought the run
	                       // TODO: make sure I count in all relevant places

public:
	// TODO: documentation...
        Single_tree_RRT(Pathways* owner);

	// Takes a node 'n' that contains DOFs vector, and produces,
	// optimizes and scores a pose object from the node using the
	// DOFs vector. The score is stored in 'n'
	//
	// Optimizations used: rotamer trial + BB/SC torisonal minimization
	//
	// Constraints: a heavy constraint does not allow big changes
	// in BB torsions
	//
	// returns: a newly allocated pose object based on info from 'n'
	pose_ns::Pose* build_and_score_pose_from_node(RRT_node* n);

	// print PDBs of motion path 'path', which is a sequence of
	// nodes
	// 'path_num' is used as a serial number for PDB file name
	// 'dir' - target directory for output PDBs
	void print_path(std::vector<RRT_node* >& path,int path_num,std::string dir);

	// reverses the list of nodes 'path'
	void reverse_path(std::vector<RRT_node* >& path);

	virtual void run();

	virtual ~Single_tree_RRT() = default;

	double energy_function_between_flex_units(pose_ns::Pose* pose);

	// when to stop the groth of the tree/s
        int fulfilled_stop_condition(int num_generated_nodes,int num_conseq_fails);


	// generates a random conformation based on the dof set by the owner
	// returns a dummy node for the random conformation
        RRT_node* generate_random_conformation();


	// returns vector of nearest nodes within some offeset_distance to n in the active_tree
	// stores the distance in "dist"
        std::vector<RRT_node *> find_nearest_nodes_CA_RMSD(RRT_node* n,double& dist);
	std::vector<RRT_node *> find_nearest_nodes_by_DOFs_vector(RRT_node* n,double& dist);

	// stop local planner:
	// 1. reached max steps of teh local planner
	// 2. reached to the vicinity of the random conformation
	// 3. reached to the dof values torsions of the random conformation - due to different
	//    bond angles - having teh same torsions values does Not garentee the same rmsd values
	// 4. energy of the new node along teh interpolation is above some value - to
	//    insure valid compact conformations with no clashes
	virtual bool fulfilled_stop_condition_local_planner
	  ( Linear_planner_iterator& linplanner_iter,
	    RRT_node* new_n,
	    RRT_node* target_n);


	int add_nodes_with_local_planner(RRT_node* from,RRT_node* to,
					 std::vector<RRT_node*>& new_generated_nodes);


	virtual void grow_DAG();

	int meet_goal_criterions(RRT_node* n);

	int findPath(RRT_node* from, RRT_node* to,std::vector<RRT_node* >& path);

	RRT_node * find_best_conformation(std::vector<RRT_node *>& v_res_nodes);

       double compute_energy_of_path(std::vector<RRT_node* >& path);
protected:
	// determines if an energy score should be considered in C_Feasible
	// (this depends on parameters - '_max_energy', '_is_soft_max_energy' and '_stddev_soft_max_energy'.
	//  The later leads to a deterministic choice)
	bool is_forbidden_energy(float score);

}; //end class single RRT


class biRRT : Single_tree_RRT
{
protected:
	int _num_steps_to_connect_trees; // for more than one tree

	double _max_rmsd_to_connect_trees; // for more than one tree - connect only nodes with less rmsd

	// when connecting trees, consider only pairs of nodes
	// with such offset RMSD from the closest pair
	double _distance_offset_nearest_neighbours_CA_RMSD_connect_trees;

public:
	biRRT(Pathways* owner);

    void run() override;

  //double compute_energy_of_path(std::vector<RRT_node* >& path);


 protected:
    // grows a DAG of the two trees, trying to connect them,
    //
    // Params:
    // connecting_nodes - [output] nodes that connect the two trees
    //                             after the run
    virtual void grow_DAG( std::vector<RRT_node*>& connecting_nodes );

    // go over all newly generated nodes,
    // and find which one gets closest to the other tree
    // put all close pairs in <output_pairs>
    void find_nearest_nodes_between_trees
      ( std::vector<std::vector<RRT_node *> >& generated_nodes_from_all_trees,
	std::vector< Node_pair >& output_pairs);


    void try_to_connect_trees
      ( std::vector<std::vector<RRT_node *> >& generated_nodes_from_all_trees,
	std::vector<RRT_node*>& connecting_nodes,
	int& num_generated_nodes );

    // dumps the pathways following the growing of the tree
    // following a run of grow_DAG, between <root_S>y & <root_T>
    // and going thru <connecting_nodes>
    //
    // IMPLEMENTATION NOTE: assumses a DAG is ready and that
    // _connecting_nodes contains nodes in the interface between
    // the two trees
    virtual void dump_pathways(
       std::vector<RRT_node*> const& connecting_nodes,
       RRT_node* root_S,
       RRT_node* root_T);

}; //end class biRRT

class Towards_partial_data_RRT : Single_tree_RRT
{

	Partial_Data _partial_data;

        double _min_partial_target;
public:
	Towards_partial_data_RRT(Pathways* owner);

	virtual bool fulfilled_stop_condition_local_planner
	  ( Linear_planner_iterator& linplanner_iter,
	    RRT_node* new_n,
	    RRT_node* target_n);
    void run() override;
    void find_close_to_partial_data_conformations(std::vector<RRT_node*>& optional_targets);
	double measure_distance_from_partial_data(RRT_node* new_n);
    std::vector<RRT_node *> find_nearest_nodes_by_DOFs_vector(RRT_node* n, double& dist );


    int  _max_steps_moving_away_from_partial_data;


}; //end class Towards_partial_data_RRT


}

#endif
