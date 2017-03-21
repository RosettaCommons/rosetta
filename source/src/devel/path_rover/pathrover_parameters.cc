////////////////////////////////////////////////////////////////////////////////
//     pathways_parameters.cc: handling parameters file for Rosetta Pathways
//
//     see pathways_parameters.h
/////////////////////////////////////////////////////////////////////////////////


#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include "pathways_parameters.h"

namespace pathways {

std::string str2upper(std::string& str)
{
	std::string ret_value;
	ret_value.resize(str.size());
	std::transform(str.begin(), str.end(), ret_value.begin(), toupper_wrapper); // convert to uppercase
	return ret_value;
}

std::string str2lower(std::string& str)
{
	std::string ret_value;
	ret_value.resize(str.size());
	std::transform(str.begin(), str.end(), ret_value.begin(), tolower_wrapper); // convert to lowercase
	return ret_value;
}


// ctr
// read parameters file (case insensitive) according to format:
Params_handler::Params_handler(std::string fname)
: is_symmetry(false)
{
	using namespace std;
	using namespace pose_ns;
	//  bool local_debug = false;

	// defaults
	this->symmetry.active = false; // default
	this->full_atom = DEFAULT_FULL_ATOM_MODE; // default
	this->energy_func_name = DEFAULT_ENERGY_FUNCTION_NAME;
	this->max_generated_nodes = DEFAULT_MAX_GENERATED_NODES;
	this->max_consecutive_fails = DEFAULT_MAX_CONSECUTIVE_FAILS;
	this->num_steps_to_generate_nodes = DEFAULT_NUM_STEPS_TO_GENERATE_NODE;
	this->num_steps_to_connect_trees = DEFAULT_NUM_STEPS_TO_CONNECT_TREES;
	this->torsion_angle_max_step = DEFAULT_TORSION_ANGLE_MAX_STEP;
	this->max_steps_local_planner = DEFAULT_MAX_STEPS_LOCAL_PLANNER;
	this->distance_offset_nearest_neighbours_CA_RMSD =
		DEFAULT_DISTANCE_OFFSET_NEAREST_NEIGHBOURS_CA_RMSD;
	this->distance_offset_nearest_neighbours_DOFs_vector =
		DEFAULT_DISTANCE_OFFSET_NEAREST_NEIGHBOURS_DOFS_VECTOR;
	this->min_level_motion = DEFAULT_MIN_LEVEL_MOTION;
	this->max_output_pathways = DEFAULT_MAX_OUTPUT_PATHWAYS;
	this->max_generated_pathways = DEFAULT_MAX_GENERATED_PATHWAYS;
	this->max_rmsd_to_connect_trees = DEFAULT_MAX_RMSD_TO_CONNECT_TREES;
	this->close_nodes_reached_random = DEFAULT_CLOSE_NODES_RECHEAD_RANDOM;
	this->max_energy = DEFAULT_MAX_ENERGY;
	this->is_soft_max_energy = DEFAULT_IS_SOFT_MAX_ENERGY;
	this->stddev_soft_max_energy = DEFAULT_STDDEV_SOFT_MAX_ENERGY;
	this->max_steps_moving_away_from_partial_data = DEFAULT_MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA;
	this->form_beta=false;
	this->form_alpha = false;
	std::ifstream data ( fname.c_str() );
	std::string line;
    this->partial_target_available = false;
    this->full_target_available = false;
	while( getline(data, line) ) {
		if(line.size() == 0 || line.at(0) == '#') // comments or blacks TODO: is this ok to check at(0) at the same close with line length check? depends on compiler? make sure no out_of_range exception
			continue;
		//  	char a = tolower('A');
		std::cout << "PARAMS line: '" << line << "'" << std::endl;
		std::istringstream line_stream(str2upper(line));
		std::string tag;
		line_stream >> tag;
		// TODO: switch to XML and SCHEME?
		if( tag == "DOF"){
			DOF_param dof_p;
			line_stream >> dof_p; // TODO: format checking?
			this->dofs.insert(dof_p);
		}
		else if( tag == "DOF_RANGE"){
			DOF_range_param dof_range_p;
			line_stream >> dof_range_p; // TODO: format checking?
			this->dof_ranges.insert(dof_range_p);
		}
		else if( tag == "FOLD_ROOT_RESIDUE"){
			line_stream >> this->fold_root_res; // TODO: format checking?
		}
		else if (tag == "CUT_CHAIN"){
			Pdbres_param cut_chain_p;
			line_stream >> cut_chain_p;
			this->cuts.insert(cut_chain_p);
		}
		else if (tag == "JUMP"){
			Jump_param jump_p;
			line_stream >> jump_p;
			this->jumps.insert(jump_p);
		}
		else if(tag == "PDB"){
			// NOTE: read pdb path without application of uppercase! this is a bit cumbersome TODO: something smarter for non-uppercase params?
			std::string type, pdb_path, tmp;
			line_stream >> type;
			std::istringstream line_stream_noupper(line); // reprocess line "PDB TYPE <pdb_path>"
			line_stream_noupper >> tmp >> tmp >> pdb_path;
			debug_assert((type == "SRC" || type == "TRG" || type == "TRG_PARTIAL" ));
			this->pdbs[type] = pdb_path;
            if (type == "TRG_PARTIAL") partial_target_available = true;
            if (type == "TRG") full_target_available = true;
		}
		else if (tag == "SYMMETRIC"){
			line_stream >> this->symmetry;
			this->symmetry.active = true;
			std::cout << "SYMMETRIC: " << symmetry.symm_type
			<< " from " << symmetry.from_pdb_chain
			<< " to " << symmetry.to_pdb_chains.size() << " clones" << std::endl;
		}
		else if (tag == "ALGORITHM"){
			// "ALGORITHM <algo_name>"
			line_stream >> this->algo_name;
			debug_assert((algo_name == "SINGLE_TREE") ||
					(algo_name == "BI_TREE") ||
					(algo_name == "RRT_PARTIAL_TARGET_INFO"));
		}
		else if (tag == "ENERGY_FUNC" ){
			// "ENERGY_FUNC <func_name>" // TODO: separate FA from Centroid function? several choices?
			line_stream >> this->energy_func_name;
			debug_assert(energy_func_name == "SCORE12" ||
					energy_func_name == "SCORE12_PATHWAYS" ||
					energy_func_name == "SCORE4" ||
					energy_func_name == "SCORE_PATHWAY_CENTROID" ||
			                energy_func_name == "SCORE6");
			energy_func_name = str2lower(energy_func_name);
		}
		else if (tag == "ENERGY_REWEIGHT" ){
			// "ENERGY_REWEIGHT <term-name> <reweight_factor>"
			// Reweights terms in the energy function
			string energy_term_name;
			double reweight_factor;
			line_stream >> energy_term_name >> reweight_factor;
			this->reweight_energy_map[ energy_term_name ] = reweight_factor;
			std::cout << "*** ENERGY_REWEIGHT param not implemented yet ***" << std::endl;
			exit(1);
		}
		else if (tag == "VDW_SCALING" ){
			// "VDW_SCALING <scale_factor>"
			// Rescales the radii of atoms by a certain factor (1 = 100%)
			line_stream >> this->vdw_scale_factor;
			std::cout << "*** VDW_SCALING param not implemented yet ***" << std::endl;
			exit(1);
		}
		else if(tag == "MAX_ENERGY"){
			line_stream >> max_energy;
		}
		else if(tag == "SOFT_MAX_ENERGY"){
			line_stream >> is_soft_max_energy;
		}
		else if(tag == "STDDEV_SOFT_MAX_ENERGY"){
			line_stream >> stddev_soft_max_energy;
		}
		else if (tag == "MAX_GENERATED_NODES"){
			line_stream >> max_generated_nodes;
		}
		else if (tag == "MAX_STEPS_MOVING_AWAY_FROM_PARTIAL_DATA"){
			line_stream >> max_steps_moving_away_from_partial_data;
		}
		else if (tag == "MAX_CONSECUTIVE_FAILS"){
			line_stream >> max_consecutive_fails;
		}
		else if (tag == "NUM_STEPS_TO_GENERATE_NODE"){
			line_stream >> num_steps_to_generate_nodes;
		}
		else if(tag == "NUM_STEPS_TO_CONNECT_TREES"){
			line_stream >> num_steps_to_connect_trees;
		}
		else if(tag == "TORSION_ANGLE_MAX_STEP"){
			line_stream >> torsion_angle_max_step;
		}
		else if(tag == "MAX_STEPS_LOCAL_PLANNER"){
			line_stream >> max_steps_local_planner;
		}
		else if(tag == "DISTANCE_OFFSET_NEAREST_NEIGHBOURS_CA_RMSD"){
			line_stream >> distance_offset_nearest_neighbours_CA_RMSD;
		}
		else if(tag == "DISTANCE_OFFSET_NEAREST_NEIGHBOURS_CA_RMSD_CONNECT_TREES"){
			line_stream >> distance_offset_nearest_neighbours_CA_RMSD_connect_trees;
		}
		else if(tag == "DISTANCE_OFFSET_NEAREST_NEIGHBOURS_DOFS_VECTOR"){
			line_stream >> distance_offset_nearest_neighbours_DOFs_vector;
		}
		else if(tag == "MIN_LEVEL_MOTION"){
			line_stream >> min_level_motion;
		}
		else if(tag == "MAX_OUTPUT_PATHWAYS"){
			line_stream >> max_output_pathways;
		}
		else if(tag == "MAX_GENERATED_PATHWAYS"){
			line_stream >> max_generated_pathways;
		}
		else if(tag == "MAX_RMSD_TO_CONNECT_TREES"){
			line_stream >> max_rmsd_to_connect_trees;
		}
		else if(tag == "CLOSE_NODES_RECHEAD_RANDOM"){
			line_stream >> close_nodes_reached_random;
		}
		else if(tag == "FULL_ATOM"){
			std::string fa_string;
			line_stream >> fa_string;
			debug_assert(fa_string == "TRUE" || fa_string == "T" || fa_string == "1" ||
					fa_string == "FALSE" || fa_string == "F" || fa_string == "0" );
			this->full_atom = (fa_string == "TRUE" || fa_string == "T" || fa_string == "1");
			std::cout << "*** Full-Atom mode: " << full_atom << std::endl;
		}
		// Partial Data -
		// ------ begin --------
		else if( tag == "LINES_ANGLE"){
			PD_line_angle line_angle;
			line_stream >> line_angle; // TODO: format checking?
			this->line_angle_recs.push_back(line_angle);
		}
		else if( tag == "LINES_DISTANCE"){
			PD_line_distance line_distance;
			line_stream >> line_distance; // TODO: format checking?
			this->line_distance_recs.push_back(line_distance);
		}
		else if( tag == "CENTROIDS_DISTANCE"){
			PD_centroid_distance center_distance;
			line_stream >> center_distance; // TODO: format checking?
			this->centroid_distance_recs.push_back(center_distance);
		}
		else if( tag == "FORM_ALPHA"){
			//PD_form_alpha form_alpha;
			//line_stream >> form_alpha; // TODO: format checking?
			//this->form_alpha_recs.push_back(form_alpha);
			this->form_alpha = true;
		}
		else if( tag == "FORM_BETA"){
			this->form_beta=true;
		}
		else if( tag == "MATCH"){
			PD_match match;
			line_stream >> match; // TODO: format checking?
			this->match_recs.push_back(match);
		}
		else if( tag == "MATCH_RMSD"){
			PD_match_rmsd match;
			line_stream >> match; // TODO: format checking?
			this->match_rmsd_recs.push_back(match);
		}
		// Partial data -
		// ------ end --------
		else
		{
			std::cout << "Unknown parameters line tag: '" << tag << "'" << std::endl;
			exit(1);
		}
	}
}

/***********************************************************/
/**************** Params Structures Methods ****************/
/***********************************************************/

bool DOF_param::operator<(DOF_param const&other) const
{
	// first by res, then by type
	if(pdb_res != other.pdb_res)
		return pdb_res < other.pdb_res;
	return s_dof_type < other.s_dof_type;
}


// DOF: "<type> <chain> <residue>"
// type: "phi" or "psi" or "residue" (i.e. phi+psi)
// chain: PDB chain, "_" for null chain
// residue: PDB residue
std::istream& operator >>(std::istream &is,DOF_param &p)
{
	is >> p.s_dof_type >> p.pdb_res.chain >> p.pdb_res.res_id;
	p.uni_dev = DEFAULT_UNI_MAX_DEV;
	p.std_dev = DEFAULT_DOF_STD_DEV;
	if(!is.eof()){
	  is >> p.uni_dev;
	}
	if(!is.eof()){
		is >> p.std_dev;
	}
	std::cout << "DOF: " << p.s_dof_type << ", "
		  << p.pdb_res.chain << ", " << p.pdb_res.res_id << ", unifrom: "
		  << p.uni_dev << ", gauss-dev: " << p.std_dev
		  << std::endl;
	return is;
}

bool DOF_range_param::operator<(DOF_range_param const&other) const
{
	// first by from res, then to res, then by type
	if(pdb_chain != other.pdb_chain)
		return pdb_chain < other.pdb_chain;
	if(from_pdb_res != other.from_pdb_res)
		return from_pdb_res < other.from_pdb_res;
	if(to_pdb_res != other.to_pdb_res)
		return to_pdb_res < other.to_pdb_res;
	return s_dof_type < other.s_dof_type;
}

// DOF_RANGE: "<type> <chain> <from_res> <to_res>"
// see "DOF" comments
std::istream& operator >>(std::istream &is,DOF_range_param &p)
{
	is >> p.s_dof_type >> p.pdb_chain >> p.from_pdb_res >> p.to_pdb_res;
	p.uni_dev = DEFAULT_UNI_MAX_DEV;
	p.std_dev = DEFAULT_DOF_STD_DEV;
	if(!is.eof())
		is >> p.uni_dev;
	if(!is.eof())
		is >> p.std_dev;
	std::cout << "DOF_RANGE>>(): " << p.s_dof_type << " " << p.pdb_chain << " " << p.from_pdb_res << "->" << p.to_pdb_res << " uniform-dev: " << p.uni_dev << " gauss_dev: " << p.std_dev << std::endl;
	return is;
}

// cutting a chain after a certain res
bool Pdbres_param::operator<(Pdbres_param const& other) const
{
	return pdbres < other.pdbres;
}

// "<chain> <pdb_res>"
std::istream& operator >>(std::istream &is, Pdbres_param &p)
{
	is >> p.pdbres.chain >> p.pdbres.res_id;
	return is;
}

bool Jump_param::operator<(Jump_param const&other) const
{
	// by from res, then to res
	if(from_pdbres != other.from_pdbres)
		return from_pdbres < other.from_pdbres;
	return to_pdbres < other.to_pdbres;
}

// "JUMP <from_chain> <from_res> <to_chain> <to_res>
std::istream& operator >>(std::istream &is,Jump_param &p)
{
	is >> p.from_pdbres.chain >> p.from_pdbres.res_id >> p.to_pdbres.chain >> p.to_pdbres.res_id;
	return is;
}

bool Symmetric_param::operator<(Symmetric_param const&other) const
{
	// from_pdb_chain should uniquely define a symmetry
	return from_pdb_chain < other.from_pdb_chain;
}

// SYMMETRIC "<from_chain> <to_chain_list>", e.g. "A BCD"
// i.e. <to_chain_list> are symmetric to a reference chain <from_chain>
std::istream& operator >>(std::istream &is,Symmetric_param &p)
{
	is >> p.symm_type >> p.from_pdb_chain;
	while(!is.eof()){
		char ch;
		is >> ch;
		p.to_pdb_chains.insert(ch);
	}
	return is;
}


// LINES_ANGLE <chain1> <from1> <to1> <chain2> <from2> <to2> = <angle(degree)>
std::istream& operator >>(std::istream &is,PD_line_angle &p){
	std::string eq;
	is >> p.chain1 >> p.from1 >> p.to1 >> p.chain2 >>  p.from2 >> p.to2 >> eq >> p.angle;
	if (p.from1 > p.to1){
		int tmp = p.from1;
		p.from1 = p.to1;
		p.to1 = tmp;//p.from1;
	}
	if (p.from2 > p.to2){
		int tmp = p.from2;
		p.from2 = p.to2;
		p.to2 = tmp;//p.from2;
	}

	if (eq != "="){
		std::cout<<"bad command line: LINES_ANGLE <chain1> <from1> <to1> <chain2> <from2> <to2> = <angle(degree)>"<<std::endl;
		exit(0);
	}
	return is;
}


// LINES_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
std::istream& operator >>(std::istream &is,PD_line_distance &p){
	std::string eq;
	is >> p.chain1 >> p.from1 >> p.to1 >> p.chain2 >>  p.from2 >> p.to2 >> eq >> p.distance;
	//int tmp;
	if (p.from1 > p.to1){
		int tmp = p.from1;
		p.from1 = p.to1;
		p.to1 = tmp;//p.from1;
	}
	if (p.from2 > p.to2){
		int tmp = p.from2;
		p.from2 = p.to2;
		p.to2 = tmp;//p.from2;
	}

	if (eq != "="){
		std::cout<<"bad command line: LINES_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>"<<std::endl;
		exit(0);
	}
	return is;
}

// CENTROIDS_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>
std::istream& operator >>(std::istream &is,PD_centroid_distance &p){
	std::string eq;
	is >> p.chain1 >> p.from1 >> p.to1 >> p.chain2 >>  p.from2 >> p.to2 >> eq >> p.distance;
	if (p.from1 > p.to1){
		int tmp = p.from1;
		p.from1 = p.to1;
		p.to1 = tmp;//p.from1;
	}
	if (p.from2 > p.to2){
		int tmp = p.from2;
		p.from2 = p.to2;
		p.to2 = tmp;//p.from2;
	}

	if (eq != "="){
		std::cout<<"bad command line: CENTROIDS_DISTANCE <chain1> <from1> <to1> <chain2> <from2> <to2> = <distance>"<<std::endl;
		exit(0);
	}
	return is;
}


// FORM_ALPHA <chain> <from> <to>
std::istream& operator >>(std::istream &is,PD_form_alpha &p){

	is >> p.chain >> p.from >> p.to;
	if (p.from > p.to){
		int tmp = p.from;
		p.from = p.to;
		p.to = tmp;//p.from;
	}

	return is;
}

// MATCH target A 15 20  source A 20 30 id
std::istream& operator >>(std::istream &is,PD_match &p){
	std::string s1,s2;
	char chain1,chain2;
	int from1,from2,to1,to2;
	is >> p.id >> s1 >> chain1 >> from1 >> to1 >> s2 >> chain2 >>  from2 >> to2;
	if (from1 > to1){
		int tmp = from1;
		from1 = to1;
		to1 = tmp;//from1;
	}
	if (from2 > to2){
		int tmp = from2;
		from2 = to2;
		to2 = tmp;//from2;
	}
	s1 = str2upper(s1);
	s2 = str2upper(s2);
	if ((s1 == "TARGET") && (s2 == "SOURCE")){
		p.chain_s = chain2;
		p.chain_t = chain1;
		p.from_s = from2;
		p.to_s = to2;
		p.from_t = from1;
		p.to_t = to1;
	}
	else if ((s2 == "TARGET") && (s1 == "SOURCE")){
		p.chain_s = chain1;
		p.chain_t = chain2;
		p.from_s = from1;
		p.to_s = to1;
		p.from_t = from2;
		p.to_t = to2;
	}

	else{
		std::cout<<"bad command line: MATCH target <chain> <from> <to> <source> <chain> <from> <to> "<<std::endl;
		exit(0);
	}
	return is;


}


// MATCH_rmsd target A 15 20  source A 20 30 id
std::istream& operator >>(std::istream &is,PD_match_rmsd &p){
	std::string s1,s2;
	char chain1,chain2;
	int from1,from2,to1,to2;
	is >> s1 >> chain1 >> from1 >> to1 >> s2 >> chain2 >>  from2 >> to2;
	if (from1 > to1){
		int tmp = from1;
		from1 = to1;
		to1 = tmp;//from1;
	}
	if (from2 > to2){
		int tmp = from2;
		from2 = to2;
		to2 = tmp;//from2;
	}
	s1 = str2upper(s1);
	s2 = str2upper(s2);
	if ((s1 == "TARGET") && (s2 == "SOURCE")){
		p.chain_s = chain2;
		p.chain_t = chain1;
		p.from_s = from2;
		p.to_s = to2;
		p.from_t = from1;
		p.to_t = to1;
	}
	else if ((s2 == "TARGET") && (s1 == "SOURCE")){
		p.chain_s = chain1;
		p.chain_t = chain2;
		p.from_s = from1;
		p.to_s = to1;
		p.from_t = from2;
		p.to_t = to2;
	}

	else{
		std::cout<<"bad command line: MATCH_RMSD target <chain> <from> <to> <source> <chain> <from> <to> "<<std::endl;
		exit(0);
	}
	return is;


}


} // namespace pathways

