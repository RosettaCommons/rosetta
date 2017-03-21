////////////////////////////////////////////////////////////////////////////////
//     pathways_parameters.h: module for handling run-time parameters of
//                            Rosetta Pathways
//
//	   Author: Barak Raveh & Angela Enosh
//     Created: 12/12/2007
//
// methods list:
//     TODO: fill this in
//
/////////////////////////////////////////////////////////////////////////////////

#ifndef INCLUDED_devel_path_rover_pathrover_parameters_hh
#define INCLUDED_devel_path_rover_pathrover_parameters_hh

#define DEFAULT_DOF_STD_DEV 100 // std. dev for gaussian sampling of DOF values
#define DEFAULT_UNI_MAX_DEV 180 // maximal range for uniform sampling of DOF values

// ObjexxFCL Headers
#include <ObjexxFCL/ObjexxFCL.hh>
#include <fstream>
#include <iostream>
#include <string>


#include "pathways_DOFs_manager.h"

namespace pathways{

// TODO: add as user parameters
#define DEFAULT_ENERGY_FUNCTION_NAME "score12"
#define DEFAULT_FULL_ATOM_MODE true
#define DEFAULT_MAX_GENERATED_NODES 3000//1000
#define DEFAULT_MAX_CONSECUTIVE_FAILS 100
#define DEFAULT_NUM_STEPS_TO_GENERATE_NODE 1
#define DEFAULT_NUM_STEPS_TO_CONNECT_TREES 5
#define DEFAULT_TORSION_ANGLE_MAX_STEP 2.0
#define DEFAULT_MAX_STEPS_LOCAL_PLANNER 50
#define DEFAULT_DISTANCE_OFFSET_NEAREST_NEIGHBOURS_CA_RMSD 0.01
#define DEFAULT_DISTANCE_OFFSET_NEAREST_NEIGHBOURS_CA_RMSD_CONNECT_TREES 1 // BiTree - more liberal since we want to connect a lot of nodes
#define DEFAULT_DISTANCE_OFFSET_NEAREST_NEIGHBOURS_DOFS_VECTOR 0.01
#define DEFAULT_MIN_LEVEL_MOTION 10
#define DEFAULT_MAX_OUTPUT_PATHWAYS 10
#define DEFAULT_MAX_GENERATED_PATHWAYS 30
#define DEFAULT_MAX_RMSD_TO_CONNECT_TREES 2.5 //2.5 // rmsd <= between two nodes of different tree -> connect trees
#define DEFAULT_CLOSE_NODES_RECHEAD_RANDOM 0.8
#define DEFAULT_MAX_ENERGY 300 // maximum energy score for C-feasible space (also see is stochastic max energy)
#define DEFAULT_IS_SOFT_MAX_ENERGY false // softens the max-energy term (and C-feasible space) into a probabilistic space that depends on energy (specifically a conformation is in C-feasible if energy < gaussian around MAX-ENERGY, satdev = STDDEV_SOFT_MAX_ENERGY)
#define DEFAULT_STDDEV_SOFT_MAX_ENERGY 25 // see IS_SOFT_MAX_ENERGY
#define DEFAULT_MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA 3

/***************************************/
/***** specific params structures: *****/
/***************************************/
// partial_data_file
    // MATCH target A 15-20 B 14 16 source A 20-30 2 4
    // LINES_ANGLE A 15 20 A 20 30 = 45 (degrees) usually will follow the below line
    // CENTROIDS_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
    // LINES_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
    // FORM_ALPHA B 15 20

// LINES_ANGLE <chain1> <from1> <to1> <chain2> <from2> <to2> = <angle(degree)>
struct PD_line_angle{ //PD - partial data
	char chain1;
	char chain2;
	int from1,from2,to1,to2;
	double angle;
};

std::istream& operator >>(std::istream &is,PD_line_angle &p);


// LINES_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
struct PD_line_distance{
	char chain1;
	char chain2;
	int from1,from2,to1,to2;
	double distance;
};

std::istream& operator >>(std::istream &is,PD_line_distance &p);

// CENTROIDS_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
struct PD_centroid_distance{
	char chain1;
	char chain2;
	int from1,from2,to1,to2;
	double distance;
};

std::istream& operator >>(std::istream &is,PD_centroid_distance &p);


// FORM_ALPHA <chain> <from> <to>
struct PD_form_alpha{
	char chain;
	int from,to;
};

std::istream& operator >>(std::istream &is,PD_form_alpha &p);

// MATCH target A 15 20  source A 20 30
// CENTROIDS_DISTANCE <type> <from1> <to1> <type> <chain2> <from2> <to2>
struct PD_match{
	char chain_s;
	char chain_t;
	int from_s,from_t,to_s,to_t;
	int id;
};

std::istream& operator >>(std::istream &is,PD_match &p);


// MATCH_RMSD target A 15 20  source A 20 30
// CENTROIDS_DISTANCE <type> <from1> <to1> <type> <chain2> <from2> <to2>
struct PD_match_rmsd{
	char chain_s;
	char chain_t;
	int from_s,from_t,to_s,to_t;
  //int id;
};

std::istream& operator >>(std::istream &is,PD_match_rmsd &p);


struct DOF_param{
	std::string s_dof_type;
	Pdbres_id pdb_res;
	double std_dev, uni_dev;

	// by chain,res and then dof_type
	bool operator<(DOF_param const&other) const;
};

// DOF: "<type> <chain> <residue>"
// type: "phi" or "psi" or "residue" (i.e. phi+psi)
// chain: PDB chain, "_" for null chain
// residue: PDB residue
std::istream& operator >>(std::istream &is,DOF_param &p);

struct DOF_range_param{
	std::string s_dof_type;
	char pdb_chain;
	int from_pdb_res, to_pdb_res;
	double std_dev, uni_dev;

	bool operator<(DOF_range_param const&other) const;
};


// DOF_RANGE: "<type> <chain> <from_res> <to_res>"
// see "DOF" comments
std::istream& operator >>(std::istream &is,DOF_range_param &p);

// Parameter that describes a single resiude
struct Pdbres_param{
	Pdbres_id pdbres;

	bool operator<(Pdbres_param const&other) const;
};

// "<chain> <pdb_res>"
std::istream& operator >>(std::istream &is,Pdbres_param &p);


// a fold-tree jump from one residue to another
struct Jump_param{
	Pdbres_id from_pdbres;
	Pdbres_id to_pdbres;

	bool operator<(Jump_param const&other) const;
};

// "JUMP <from_chain> <from_res> <to_chain> <to_res>
std::istream& operator >>(std::istream &is,Jump_param &p);


// defines one chain as symmetric to another chain
// (TODO: all internal jumps should be duplicated as well)
struct Symmetric_param{
	std::string symm_type;
	char from_pdb_chain;
	std::set<char> to_pdb_chains; // set of target symmetries
	bool active;

	bool operator<(Symmetric_param const&other) const;
};

// SYMMETRIC "<from_chain> <to_chain_list>", e.g. "A BCD"
// i.e. <to_chain_list> are symmetric to a reference chain <from_chain>
std::istream& operator >>(std::istream &is,Symmetric_param &p);


/***************************************/
// main class for reading and handling pathways parameters file
// information about the running -
// * the list of DOFs to change while running
// * predicates
// * energy parameters to use
// * characteristics of algorithms
// * pathway alignment algorithms
// * rigid bodies and fold tree info ?
/***************************************/
class Params_handler{
public:
	typedef std::map<std::string, std::string> t_map_str2str;
public:
	t_map_str2str pdbs; // a list of pdbs (key = type, i.e. "src", value = path to pdb-file)
	bool full_atom;

	std::set<DOF_param> dofs; // DOFs defined over a single residue
	std::set<DOF_range_param> dof_ranges; // DOFs defined over a range of residues
	Pdbres_param fold_root_res; // root for fold-tree
	std::set<Pdbres_param> cuts; // cuts of fold-tree
	std::set<Jump_param> jumps; // jumps of fold-tree
	Symmetric_param symmetry;
	bool is_symmetry; // is symmetry active for molecule
	std::string algo_name; // algorithm for generating pathways // TODO: implement
	std::string energy_func_name; // from our set of custom energy functions // TODO: implement
	std::map<std::string, double> reweight_energy_map; // reweight mapping on energy terms // TODO: implement
	double vdw_scale_factor; // some scaling // TODO: implement

	int max_generated_nodes; // maximum total number of nodes in all trees
	int max_consecutive_fails; // how many consecutive local-planner failues are allowed
	int num_steps_to_generate_nodes; // for local-planner
	int num_steps_to_connect_trees; // for more than one tree
	double torsion_angle_max_step; // maximal local-planner perturbation at any given angle
	int max_steps_local_planner; // for one run of local planner
	double distance_offset_nearest_neighbours_CA_RMSD; // for find_nearest_node_CA_RMSD()
	double distance_offset_nearest_neighbours_CA_RMSD_connect_trees; // for BiRRT::find_nearest_pairs()
	double distance_offset_nearest_neighbours_DOFs_vector; // for find_nearest_node_DOFs_vector()
	int min_level_motion;
	unsigned int max_output_pathways;
	unsigned int max_generated_pathways;
	double max_rmsd_to_connect_trees; // for more than one tree - connect only nodes with less rmsd
	double close_nodes_reached_random;
	double max_energy; // max energy for inserting nodes to a tree (C-feasible)
	bool is_soft_max_energy; // softens the max_energy term (and C-feasible space) into a probabilistic space that depends on energy (specifically a conformation is in C-feasible iff a energy < gaussian around MAX-ENERGY, with std_dev of STDDEV_SOFT_MAX_ENERGY)
	double stddev_soft_max_energy; // see is_soft_max_energy

	int max_steps_moving_away_from_partial_data; // how many steps are allowed to be made without improving the partial data continuous predicate
	// Different types of partial data:
	// partial_data_file
    // MATCH target A 15-20 B 14 16 source A 20-30 2 4
    // LINES_ANGLE A 15 20 A 20 30 = 45 (degrees) usually will follow the below line
    // CENTROIDS_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
    // LINES_DISTANCE A 15 20  B 20 30 = 5 (angstrom)
    // FORM_ALPHA B 15 20
	std::vector<PD_line_angle> line_angle_recs;
	std::vector<PD_line_distance> line_distance_recs;
	std::vector<PD_centroid_distance> centroid_distance_recs;
	std::vector<PD_form_alpha> form_alpha_recs;
	std::vector<PD_match> match_recs;
	std::vector<PD_match_rmsd> match_rmsd_recs;
    bool form_beta, form_alpha;
    bool partial_target_available;
    bool full_target_available;


	// LINES_ANGLE <chain1> <from1> <to1> <chain2> <from2> <to2> = <angle(degree)>
// LINES_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
// CENTROIDS_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
// FORM_ALPHA <chain> <from> <to>
// MATCH target A 15 20  source A 20 30


public:
	Params_handler(std::string fname);

	std::string const& get_pdb_path(std::string type) const
	{
		t_map_str2str::const_iterator iter = pdbs.find(type);
		debug_assert(iter != pdbs.end());
		return iter->second;
	}
};

/* function wrapping the macro */
// NOTE: this is just a small utility stuff to wrap the macro "toupper",
// so that it can be passed as function ptr to std::transform().
// TODO: move somewhere else, to a utility header of some form
inline char toupper_wrapper (const char c) { return toupper(c); }

/* function wrapping the macro */
// NOTE: this is just a small utility stuff to wrap the macro "toupper",
// so that it can be passed as function ptr to std::transform().
// TODO: move somewhere else, to a utility header of some form
inline char tolower_wrapper (const char c) { return tolower(c); }


} // namespace pathways

#endif
