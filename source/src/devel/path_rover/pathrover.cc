////////////////////////////////////////////////////////////////////////////////
//     pathways.cc: module for generating pathways between conformations
//
//     Originally Created in R++: 29/11/2007 (Barak Raveh & Angela Enosh)
//     Transformed to mini: 27/10/2009 (Barak Raveh)
//
//     see pathways.hh
/////////////////////////////////////////////////////////////////////////////////


// Rosetta Headers
# #include "misc.h"
# #include "pack.h"
# #include "pdb.h"
#include <core/pose/Pose.hh>
# #include "pose_io.h"
# #include "pose_rotamer_trials.h"
# #include "random_numbers.h"
# #include "rotamer_trials.h"
# #include "runlevel.h"
# #include "score.h"
# #include "score_ns.h"
# #include "symmetry_info.h"
# #include "util_vector.h"

#include "pathrover.hh"
#include "pathrover_parameters.hh"

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
#include <cctype>
#include <fstream>
#include <iostream>
#include <map>
#include <cmath>
#include <set>
#include <string>
#include <cstdio>
#include <vector>

namespace protocols{
namespace pathrover{


//////////////////////////////////////////////////////////
//
// A global function for invoking the flexible peptide docking
// protocol based on the "misc" global data structure
// (although the actual implementation is via Pose)
// TODO: stop using misc at all
//
// output parameter:
//   fail - signifies any major failure in the run of the docking
///////////////////////////////////////////////////////////
// utility
// invoke pathways protocol for generating conformation pathways
void pathrover_generator_main(bool &failed)
{
	// initialize
	failed = true; // a-priori
	Pathways pw;
	if(pw.fail()){
		// TODO: print something?
		failed = true;
		return;
	}
	// run
	pw.run();
	if(pw.fail()){
		// TODO: print something?
		failed=true;
		return;
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief ctr of protocol for generating pathways of conformations with RRT class algorithms,
//         etc.
///
/// @details
/// ctr of // TODO: fill this up
///
/// @global_read
/// standard,runlevel (in runlevel.h) - documentation is runlevel dependent
///
/// @global_write
/// score (in misc.h) - // TODO: fill this up
///
/// @remarks
/// many indirect changes to global variables (score and structure related).
///
/// @references
///
/// @author Barak & Angela 29/Nov/2007
///
/////////////////////////////////////////////////////////////////////////////////
path_rover::PathRover()
	: 	 _params_h("pathways.params"), // TODO: run-time definition of file name
		_is_fail(false) // a-priori

{
	using namespace misc;
	using namespace runlevel_ns;
	using namespace pose_ns;

	bool local_debug = true;
	if ( local_debug ) std::cout << "Pathways::Pathways()" << std::endl;
	set_pose_flag( true ); // we use pose TODO: but this has a global effect - what to do? // TODO: perhaps this is done already in options.cc

	// process parameters
	initialize_from_params(_params_h);

//	// ***** Initialize local variables and data structures: *****
//
//	// set up structure pose for docking, and tree of folding hierarchy
//	//copy misc to pose
//	// (misc was already loaded, so we copy PDB structures and so on into the docking_pose data structure)
//	ideal_pos = false;
//  	coords_init = true;
//  	this->full_atom = docking::docking_fullatom_flag; // TODO: change the name of the flag
//  	if(this->full_atom){
//  		std::cout << "Pathway: <Full atom mode>" << std::endl;
//  	}
//  pose_from_misc( this->pose, fullatom, ideal_pos, coords_init );


}


/////////////////////////////////////////////////
// create a pose, a fold tree and a DOFs manager
// according to parameters
/////////////////////////////////////////////////
void
path_rover::initialize_from_params(Params_handler& param_h)
{
	typedef std::map<char, Chain_boundaries> Chain_boundaries_map;
	bool local_debug = true;
	if(local_debug){
		std::cout << "Pathways::initialize_from_params() - start" << std::endl;
		if(param_h.full_atom)
			std::cout << "FULL_ATOM mode" << std::endl;
		else
			std::cout << "CENTROID mode" << std::endl;
	}

		std::cout << "CENTROID mode" << std::endl;
	_full_atom = param_h.full_atom;

std::cout << "CENTROID mode" << std::endl;
	std::cout << "Src pdb: " << _params_h.get_pdb_path("SRC");
	bool is_ok;
	is_ok = pose_from_file(
		_src_pose,
		_params_h.get_pdb_path("SRC"), /* filename */
		_full_atom,  /* fullatom*/
		false, /* ideal_pose */
		true,  /* read_all_chains, // = false */
		'-',   /* chain */
		false, /* skip_missing, // = false */
		false  /* allow_missing // = false */
	);
	assert(is_ok);
	// read target
	if(_params_h.full_target_available)
	{
		is_ok = pose_from_file(
			_trg_pose,
			_params_h.get_pdb_path("TRG"), /* filename */
			_full_atom,  /* fullatom*/
			false, /* ideal_pose */
			true,  /* read_all_chains, // = false */
			'-',   /* chain */
			false, /* skip_missing, // = false */
			false  /* allow_missing // = false */
		);
		assert(is_ok);
		_trg_pose.dump_pdb("./output/test1.pdb");
	}
	else
	  _trg_pose = 	_src_pose; //TODO:ignore target when it doesnt exsists
	// read partial target
    if (_params_h.partial_target_available)
	{
		bool is_ok = pose_from_file(
			_trg_pose_partial,
			_params_h.get_pdb_path("TRG_PARTIAL"), /* filename */
			_full_atom,  /* fullatom*/
			false, /* ideal_pose */
			true,  /* read_all_chains, // = false */
			'-',   /* chain */
			false, /* skip_missing, // = false */
			false  /* allow_missing // = false */
		);
		assert(is_ok);
	}

	Pdb_info& pdb_info = _src_pose.pdb_info();

	// initialize DOFs managing based on loaded pose
	// TODO: initialize should clear maps and so on
	_dofs_manager.initialize(_src_pose, false /* read_pose_allow_xxx */);

	// calculate chain boundaries
	Chain_boundaries_map chain_boundaries_map;
	{
		char cur_ch = pdb_info.res_chain(1);
		std::cout << "current chain: " <<  cur_ch << std::endl;
		int cur_ch_start = 1;
		for(int i = 2; i <= _src_pose.size(); i++)
		{
			char prev_ch = cur_ch;
			cur_ch = pdb_info.res_chain(i);
			if(prev_ch != cur_ch) // changed chain
			{
				std::cout << "from chain " << prev_ch << " to " << cur_ch << std::endl;
				chain_boundaries_map.insert(
						std::make_pair(prev_ch, Chain_boundaries(cur_ch_start, i-1))
				);
				cur_ch_start = i; // update for next time
			}
		}
		// last chain
		chain_boundaries_map.insert(
				std::make_pair(cur_ch, Chain_boundaries(cur_ch_start, _src_pose.size()))
		);
	}

	// Symmetery auxilary variables:
	Symmetric_param& symm_p = param_h.symmetry;
	int src_unit_first_res = -1; // first residue of symmetry unit
	int src_unit_last_res = -1; // last...
	std::vector<int> src_unit_jumps; // list of jumps to be duplicated
	int n_clones = 0; // # of times symmetry unit is duplicated
	std::set<char>& clones_pdb_chains = symm_p.to_pdb_chains;
	if(symm_p.active)
		n_clones = clones_pdb_chains.size();
	std::vector<int> clones_first_res_list;
	std::vector< std::vector<int> > clones_jumps;

	// set first residue of symmetry units (src & clones)
	if(symm_p.active && symm_p.symm_type == "INTERNAL_ONLY")
	{
		assert(chain_boundaries_map.find(symm_p.from_pdb_chain) != chain_boundaries_map.end()); // valid chain
		assert(chain_boundaries_map[ symm_p.from_pdb_chain ].valid);
		src_unit_first_res = chain_boundaries_map[ symm_p.from_pdb_chain ].first;
		src_unit_last_res = chain_boundaries_map[ symm_p.from_pdb_chain ].last;
		for(std::set<char>::const_iterator it = clones_pdb_chains.begin(), it = clones_pdb_chains.end(); it != end; ++it )
		{
			assert(chain_boundaries_map.find(*it) != chain_boundaries_map.end()); // valid chain
			assert(chain_boundaries_map[ *it ].valid);
			clones_first_res_list.push_back(chain_boundaries_map[*it].first);
		}
	}

	// add DOFs to params handler
	{ for(std::set<DOF_param>::const_iterator it = param_h.dofs.begin(), end = param_h.dofs.end(); it != end; ++it )
	{
		if(it->s_dof_type == "PHI")
			_dofs_manager.set_phi_dof(it->pdb_res, it->std_dev, it->uni_dev);
		else if(it->s_dof_type == "PSI")
			_dofs_manager.set_psi_dof(it->pdb_res, it->std_dev, it->uni_dev);
		else if(it->s_dof_type == "RESIDUE")
			_dofs_manager.set_residue_dofs(it->pdb_res, it->std_dev, it->uni_dev);
		else{
			std::cout << "Pathways::initialize_from_params(): Bad DOF_type '" << it->s_dof_type << "'" << std::endl;
			assert(false);
		}
	} // for
	}

	// add DOF residues range to params handler
	{
	  for(std::set<DOF_range_param>::const_iterator it = param_h.dof_ranges.begin(), end = param_h.dof_ranges.end(); it != end; ++it )
	    {
	      if(it->s_dof_type == "PHI")
		for(int k = it->from_pdb_res; k <= it->to_pdb_res; k++)
		  {
		    _dofs_manager.set_phi_dof(it->pdb_chain, k, it->std_dev, it->uni_dev);
		  }
	      else if(it->s_dof_type == "PSI")
		for(int k = it->from_pdb_res; k <= it->to_pdb_res; k++)
		  {
		    _dofs_manager.set_psi_dof(it->pdb_chain, k, it->std_dev, it->uni_dev);
		  }
	      else if(it->s_dof_type == "RESIDUE")
		for(int k = it->from_pdb_res; k <= it->to_pdb_res; k++)
		  {
		    _dofs_manager.set_residue_dofs(it->pdb_chain, k, it->std_dev, it->uni_dev);
		  }
	      else{
		std::cout << "Pathways::initialize_from_params(): Bad DOF_type '" << it->s_dof_type << "'" << std::endl;
		assert(false);
	      }
	    }
	}

	// build Rosetta fold-tree according to params file:
  // (formally, a directed tree that describes the hierarchy in which
  //  torsion angles and rigid body transformations are propogated into
  //  cartesian coordinates, starting from a start resiude which is fixed)
  pose_ns::Fold_tree ft;
  typedef std::set< int > Int_set;
  Int_set vertices, cuts;
  // insert root of fold-tree:
  int root_pose_res = _dofs_manager.pdbres_to_poseres(param_h.fold_root_res.pdbres);
  assert(root_pose_res != -1);
  vertices.insert(root_pose_res);
  // add cuts for chain boundaries
  {
    vertices.insert(1); // first residue in pose
    for(Chain_boundaries_map::const_iterator it = chain_boundaries_map.begin(), end = chain_boundaries_map.end(); it != end; ++it )
      {
	cuts.insert(it->second.last);
      }
  }
  // add cuts from params file to list
  {
    for(std::set<Pdbres_param>::const_iterator it = param_h.cuts.begin(), end = param_h.cuts.end(); it != end; ++it )
      {
	int cut_pose_res = _dofs_manager.pdbres_to_poseres(it->pdbres);
	assert(cut_pose_res != -1);
	cuts.insert(cut_pose_res);
	if(local_debug)
	  std::cout << "adding params file cut: pose_res #" << cut_pose_res << " ; pdbres: " << it->pdbres << std::endl;
      }
  }
  // add all cuts as vertices (i, i+1)
  {
    for(Int_set::const_iterator it = cuts.begin(), end = cuts.end(); it != end; ++it )
      {
	vertices.insert(*it);
	vertices.insert(*it + 1);
      }
  }
  // add jumps vertices
  {
    for(std::set<Jump_param>::const_iterator it = param_h.jumps.begin(), end = param_h.jumps.end(); it != end; ++it )
      {
	int from_pose_res = _dofs_manager.pdbres_to_poseres(it->from_pdbres);
	int to_pose_res = _dofs_manager.pdbres_to_poseres(it->to_pdbres);
	assert(from_pose_res != -1 && to_pose_res != -1);
	vertices.insert(from_pose_res);
	vertices.insert(to_pose_res);
	if(local_debug)
	  std::cout << "adding params file jump: pose_res [from,to] = [" << from_pose_res << ", " << to_pose_res << "] "
		    << " ; pdbres [from,to] = [" << it->from_pdbres << ", " << it->to_pdbres << "]" << std::endl;
      }
  }
  // duplicate symmetry vertices & cuts (if needed) for all clones
  if(symm_p.active && symm_p.symm_type == "INTERNAL_ONLY")
    {
      Int_set::const_iterator it_v;
      for(it_v = ++vertices.lower_bound(src_unit_first_res); // this is first iter s.t. *it >= src_unit_first_res
	  *it_v <= src_unit_last_res ; ++it_v ) 
	{
	  if(local_debug)
	    std::cout << "duplicating " << *it_v << std::endl;
	  for(int i = 0 ; i < n_clones; i++)
	    {
	      int offset = clones_first_res_list[i]
		- src_unit_first_res;
	      vertices.insert(*it_v + offset);
	      if(local_debug)
		std::cout << "clone #" << i << "; offset " << offset << std::endl;
	      if(cuts.find(*it_v) != cuts.end()) // NOTE: a cut must be a vertex so this is complete
		{
		  cuts.insert(*it_v + offset);
		  if(local_debug)
		    std::cout << "[+ cut]" << std::endl;
		}
	    }
	}
    }

  if(local_debug)
    {
      std::cout << "Vertices for fold-tree: " << std::endl;
      for(Int_set::const_iterator it = vertices.begin(), end = vertices.end(); it != end; ++it)
	std::cout << *it << "; ";
      std::cout << std::endl;
    }

  // build all peptide edges (direction doesn't matter, will be set by "REORDER"
  {
    // go over (prev_it, cur_it) pairs
    Int_set::const_iterator	cur_it = vertices.begin();
    for(Int_set::const_iterator	prev_it = cur_it++;
	cur_it != vertices.end();
	++prev_it, ++cur_it )

      {
	if(cuts.find(*prev_it) == cuts.end()) // add only if not a cut after prev_it
	  {
	    ft.add_edge(*prev_it, *cur_it, pose_ns::Fold_tree::PEPTIDE);
	    if(local_debug)
	      std::cout << "Adding peptide edge: " << *prev_it << " , " << *cur_it << std::endl;
	  }
      }
  }
  // add jumps edges, including symmetry
  {
    int num_jumps = 0;
    clones_jumps.resize(n_clones);
    for(std::set<Jump_param>::const_iterator it = param_h.jumps.begin(), end = param_h.jumps.end(); it != end; ++it)
      {
	int from_pose_res = _dofs_manager.pdbres_to_poseres(it->from_pdbres);
	int to_pose_res = _dofs_manager.pdbres_to_poseres(it->to_pdbres);
	ft.add_edge(from_pose_res, to_pose_res, ++num_jumps);
	if(local_debug)
	  std::cout << "Adding jump edge: " << from_pose_res << " , " << to_pose_res << std::endl;
	// if needed, add symmetry jumps, and save mapping from src unit
	// (internal = jumps within the source unit)
	if(symm_p.active && symm_p.symm_type == "INTERNAL_ONLY"
	  && it->from_pdbres.chain == symm_p.from_pdb_chain
	  && it->to_pdbres.chain == symm_p.from_pdb_chain )
	  {
	    src_unit_jumps.push_back(num_jumps);
	    for(int i = 0 ; i < n_clones; i++) // for each clone
	      {
		num_jumps++;
		int offset = clones_first_res_list[i]
		  - src_unit_first_res;
		ft.add_edge(from_pose_res + offset,
		  to_pose_res + offset,
		  num_jumps);
		clones_jumps[i].push_back(num_jumps);
		if(local_debug)
		  std::cout << "Adding jump edge for symmetry clone: "
			    << from_pose_res + offset << " , "
			    << to_pose_res + offset << std::endl;
	      }
	  }
      }
  }
  // now set edge directions starting from root
  // (which results in direction of "folding", i.e. directions of DOF propogation)
  if(local_debug)
    std::cout << "Foldtree b4 reorder: " << ft << std::endl;
  ft.reorder(root_pose_res);
  if(local_debug)
    std::cout << "Foldtree after reorder(" << root_pose_res << "): " << ft << std::endl;
  assert(ft.check_fold_tree());
  // finally set fold tree
  _src_pose.set_fold_tree(ft);
  if(_params_h.algo_name == "BI_TREE")
	  _trg_pose.set_fold_tree(ft);

  // create and set final symmetry info object
  if(symm_p.active) // TODO: does it need to depend on symmetry type???
    {
      int n_res_src_unit = src_unit_last_res - src_unit_first_res + 1;
      pose_ns::Symmetry_info
	s(src_unit_first_res, n_res_src_unit, src_unit_jumps,
	  clones_first_res_list, clones_jumps);
      _src_pose.setup_symm_info(s);
      if(_params_h.algo_name == "BI_TREE")
    	  _trg_pose.setup_symm_info(s);
    }

  if(local_debug)
    std::cout << "Pathways::initialize_from_params()" << std::endl;

}


void
path_rover::run(){
	using namespace pose_ns;
	int local_debug = 2;

	std::cout << "[PathRover::run()] *** starting ***" << std::endl;

	for(int i=261; i <= 280; i++)
	  {
	    if(i==271) i++; // skip C31
	    double phi_offset = ran3()*30 - 15;
	    double psi_offset = ran3()*30 - 15;
	    double new_phi = periodic_range( _src_pose.phi(i) + phi_offset, 360.0f );
	    double new_psi = periodic_range( _src_pose.psi(i) + psi_offset, 360.0f );
	    _src_pose.set_phi(i,new_phi);
	    _src_pose.set_psi(i,new_psi);
	  }
	//	_src_pose.set_segment_extended(261,263);
	//	_src_pose.set_segment_extended(272,280);

	_src_pose.dump_pdb("./output/run_start_SRC.pdb"); // DEBUG
	_trg_pose.dump_pdb("./output/run_start_TRG.pdb"); // DEBUG

	if(local_debug >= 2){
		std::cout << "PHI/PSI sample - SOURCE structure : ";
		int i;
		for ( i=1; i <= _src_pose.size(); i += 1 /*10 sample only...*/) {
			std::cout << "[res #" << i
				  << " phi " << _src_pose.phi(i)
				  << " psi = " << _src_pose.psi(i) << "] ";
		}
		std::cout << std::endl;
		std::cout << "PHI/PSI sample - TARGET structure : ";
		for ( i=1; i <= _trg_pose.size(); i += 1 /*10 sample only...*/) {
			std::cout << "[res #" << i
				  << " phi " << _trg_pose.phi(i)
				  << " psi = " << _trg_pose.psi(i) << "] ";
		}
		std::cout << std::endl;
		std::cout << "PHI/PSI sample - Delta(SRC, TRG) : ";
                for ( i=1; i <= _trg_pose.size(); i += 1 /*10 sample only...*/) {
		  std::cout << "[res #" << i
			    << " phi " << _trg_pose.phi(i) - _src_pose.phi(i)
			    << " psi = " << _trg_pose.psi(i) - _src_pose.psi(i) << "] "
		            << std::endl;
                }

		std::cout << std::endl;

		// apply dofs of each pose to the other pose and see what we get...
		pose_ns::Pose
		  pose_src_with_trg_dofs,
		  pose_trg_with_src_dofs;
		pose_src_with_trg_dofs = _src_pose;
		pose_trg_with_src_dofs = _trg_pose;
		std::vector<double> src_dofs, trg_dofs;
		src_dofs = _dofs_manager.get_dofs_values_vector(_src_pose);
		trg_dofs = _dofs_manager.get_dofs_values_vector(_trg_pose);
		_dofs_manager.apply_dofs_values_vector(
		  pose_src_with_trg_dofs, trg_dofs);
		_dofs_manager.apply_dofs_values_vector(
		  pose_trg_with_src_dofs, src_dofs);
		pose_src_with_trg_dofs.dump_pdb("./output/SRC_with_TRG_DOFS.pdb");
		pose_trg_with_src_dofs.dump_pdb("./output/TRG_with_SRC_DOFS.pdb");
	}

// 	// idealize the bond angles & length of the peptide, to allow fragment insertion
// 	_src_pose.insert_ideal_bonds( 1, _src_pose.size() );
// 	_trg_pose.insert_ideal_bonds( 1, _trg_pose.size() );
// 	_src_pose.dump_pdb("./output/run_start_SRC_ideal_bonds.pdb"); // DEBUG
// 	_trg_pose.dump_pdb("./output/run_start_TRG_ideal_bonds.pdb"); // DEBUG

// 	if(local_debug){
// 		std::cout << "PHI/PSI sample - ideal SOURCE structure : ";
// 		int i;
// 		for ( i=1; i <= _src_pose.size(); i += 10 /*sample only...*/) {
// 			std::cout << "[res #" << i << " phi " << _src_pose.phi(i) << " psi = " << _src_pose.psi(i) << "] ";
// 		}
// 		std::cout << std::endl;
// 		std::cout << "PHI/PSI sample - ideal TARGET structure : ";
// 		for ( i=1; i <= _trg_pose.size(); i += 10 /*sample only...*/) {
// 			std::cout << "[res #" << i << " phi " << _trg_pose.phi(i) << " psi = " << _trg_pose.psi(i) << "] ";
// 		}
// 		std::cout << std::endl;
// 	}

	// FIRST OPTIMIZE START STRUCTURE
	setup_score_weight_map(
		_score_weight_map,
		score_function_from_string( _params_h.full_atom , _params_h.energy_func_name )
	);
	std::cout << "Energy function: " << _params_h.energy_func_name << std::endl;

	if(local_debug){
	  // DEBUG: print scoring
	  _src_pose.score( _score_weight_map );
	  std::cout << "pose <" << _params_h.energy_func_name << "> with source start structure: " << std::endl;
	  _src_pose.show_scores(std::cout); std::cout << std::endl;

	  if(_params_h.algo_name == "BI_TREE")
	  {
		  _trg_pose.score( _score_weight_map );
		  std::cout << "pose <" << _params_h.energy_func_name << "> with target start structure:" << std::endl;
		  _trg_pose.show_scores(std::cout); std::cout << std::endl;
	  }
  	}

	if(_params_h.full_atom)
	{
		// full repack of src & trg pose
		FArray1D_bool allow_repack(_src_pose.size(), false);
		for( int i = 1; i <= _src_pose.size(); ++i ) {
			allow_repack(i) = true;
		}
		_src_pose.repack( allow_repack, true/*include_current*/ );
		_trg_pose.repack( allow_repack, true/*include_current*/ );

		if(_params_h.algo_name == "BI_TREE")
		{
			// repack trg pose
			FArray1D_bool allow_repack_t(_trg_pose.size(), false);
			for(int i = 1; i <= _trg_pose.size(); ++i ) {
				allow_repack(i) = true;
			}
			_trg_pose.repack( allow_repack, true/*include_current*/ );
		}

		if(local_debug){
			// DEBUG: print scoring
			_src_pose.score( _score_weight_map );
			std::cout << "pose <" << _params_h.energy_func_name << "> with source start structure +REPACK: " << std::endl;
			_src_pose.show_scores(std::cout); std::cout << std::endl;

			if(_params_h.algo_name == "BI_TREE")
			{
			 _trg_pose.score( _score_weight_map );
			 std::cout << "pose <" << _params_h.energy_func_name << "> with target start structure + REPACK:" << std::endl;
			 _trg_pose.show_scores(std::cout); std::cout << std::endl;
			}
		}
	}

	// ********** Now really run.... **************

	std::cout<<"CHOOSE ALGORITHM -------------------------"<<std::endl;
	if(_params_h.algo_name == "SINGLE_TREE")
	{
		Single_tree_RRT rrt(this);
		rrt.run();
	}
	if(_params_h.algo_name == "BI_TREE")
	{
		std::cout<<"Start bi-RRT"<<std::endl;
		biRRT rrt(this);
		rrt.run();
	}
	if(_params_h.algo_name == "RRT_PARTIAL_TARGET_INFO")

	{
		    int count_flags=0;
		    std::cout<<"Start RRT_PARTIAL_TARGET_INFO"<<std::endl;
		    _partial_data = Partial_Data(this);
                        if (_params_h.partial_target_available){
                            _partial_data.set_partial_data();
                        }
			if (_params_h.line_angle_recs.size() > 0){
				count_flags++;
			    _partial_data.set_line_angle_recs(_params_h.line_angle_recs);
			}
			if (_params_h.line_distance_recs.size() > 0){
				count_flags++;
			    _partial_data.set_line_distance_recs(_params_h.line_distance_recs);
			}
			if (_params_h.centroid_distance_recs.size() > 0){
				count_flags++;
				_partial_data.set_centroid_distance_recs(_params_h.centroid_distance_recs);

			}
			if (_params_h.form_alpha_recs.size() > 0){
				count_flags++;
			   _partial_data.set_form_alpha_recs(_params_h.form_alpha_recs);
			}
			if (_params_h.match_recs.size() > 0){
				count_flags++;
				_partial_data.set_match_recs(_params_h.match_recs);
			}
                        if (_params_h.match_rmsd_recs.size() > 0){
				count_flags++;
				_partial_data.set_match_rmsd_recs(_params_h.match_rmsd_recs);
			}
                        if (_params_h.form_alpha){
                          count_flags++;
			  if (!_params_h.full_atom){
			    std::cout<<"H-bonds are considered in Full Atom mode"<<std::endl;
                            exit(0);
			  }
			    _partial_data.set_alpha_flag();
			}
                        if (_params_h.form_beta){
                            count_flags++;
			    _partial_data.set_beta_flag();
			}
			if(count_flags==0){
			  std::cout<<"Missing Partial Data"<<std::endl;
				exit(0);
			}

		Towards_partial_data_RRT rrt(this);
		rrt.run();
	}
}


} // namespace pathways
} // namespace protocols
