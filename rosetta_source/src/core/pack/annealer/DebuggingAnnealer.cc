// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/DebuggingAnnealer.cc
/// @brief  Debugging annealer class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/DebuggingAnnealer.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>

// Rosetta Headers
//#include "after_opts.h"
//#include "DebuggingAnnealer.h"
//#include "InteractionGraphBase.h"
//#include "param.h"
//#include "random_numbers.h"
//#include "RotamerSet.h"

#include <ObjexxFCL/Fmath.hh>

#include <fstream>
#include <istream>
#include <iostream>

namespace core{
namespace pack{
namespace annealer{

////////////////////////////////////////////////////////////////////////////////
/// @begin DebuggingAnnealer::DebuggingAnnealer()
///
/// @brief
/// constructor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors
///
/// @last_modified

////////////////////////////////////////////////////////////////////////////////
DebuggingAnnealer::DebuggingAnnealer(
	utility::vector0< int > & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	pack::InteractionGraphBase * ig,
	const RotamerSet * p_rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	rot_to_pack,
	(int) rot_to_pack.size(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	p_rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
}

DebuggingAnnealer::DebuggingAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	pack::InteractionGraphBase * ig,
	const RotamerSet * p_rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq
):
	RotamerAssigningAnnealer(
	(ig->get_num_total_states()),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	p_rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
}

////////////////////////////////////////////////////////////////////////////////
/// @begin DebuggingAnnealer::~DebuggingAnnealer()
///
/// @brief
/// virtual destructor
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
DebuggingAnnealer::~DebuggingAnnealer()
{}

////////////////////////////////////////////////////////////////////////////////
/// @begin FixbbSimAnnealer::run
///
/// @brief
/// sim_annealing for fixed backbone design mode
///
/// @detailed
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors
///
/// @last_modified
////////////////////////////////////////////////////////////////////////////////
void DebuggingAnnealer::run()
{
	using namespace param;

	//--------------------------------------------------------------------
	//internal variables

	int const nmoltenres = ig_->get_num_nodes();
	float currentenergy;
	//FArray1D_int list( p_rotamer_set_->nrotamers() );
	//FArray1D_int state_on_node( nmoltenres,0 );
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_int current_state( nmoltenres, 0 );
	//FArray1D_float loopenergy(maxouteriterations,0.0);

	//kwk internal rot_index expanded to handled non amino acid molten residues
	//kwk rot_to_moltenres(1,X)=moltenres_id
	//kwk rot_to_moltenres(2,X)=moltenres_state
	FArray2D_int rot_2_moltenres( 2,ig_->get_num_total_states(),0);

	//bk variables for calculating rotamer frequencies during simulation
	//int nsteps = 0;

	//--------------------------------------------------------------------
	//initialize variables

	currentenergy = 0.0;

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	int rot_tmp=1;
	for( int ii =1; ii<=nmoltenres; ++ii){
		for(int jj =1; jj<=ig_->get_num_states_for_node( ii); ++jj){
			rot_2_moltenres(1,rot_tmp)=ii;
			rot_2_moltenres(2,rot_tmp)=jj;
			rot_tmp++;
		}
	}

	//--------------------------------------------------------------------
	if ( num_of_rot_to_pack_ == 0 ) return;

	setup_iterations();

	FArray1D_float previous_nsteps_for_rot(p_rotamer_set_->nrotamers(), 0.0);

	//int outeriterations = get_outeriterations();

	std::string annealer_file( stringafteroption( "db_annealer_file" ));
	std::ifstream instructions( annealer_file.c_str() );

	while ( instructions )
	{
		int node_to_change;
		instructions >> node_to_change;

		if ( node_to_change <= 0 || node_to_change > nmoltenres )
		{
			std::cerr << "Error in instructions for Debugging Annealer: node_to_change out of range." << std::endl;
			std::cerr << "Requires 0 < node_in_range <= nmoltenres.  nmoltenres: " << nmoltenres << " node_to_change: " << node_to_change << std::endl;
			return;
		}

		int new_state_for_node;
		instructions >> new_state_for_node;

		if ( new_state_for_node <= 0 || new_state_for_node > ig_->get_num_states_for_node( node_to_change ) )
		{
			std::cerr << "Error in instructions for Debugging Annealer: new_state_for_node out of range." << std::endl;
			std::cerr << "Requires 0 < new_state_for_node <= num_states_for_node( " <<  node_to_change << ");" << std::endl;
			std::cerr << "new_state_for_node: " << new_state_for_node << " num_states_for_node( " <<  node_to_change << "): ";
			std::cerr << ig_->get_num_states_for_node( node_to_change ) << std::endl;
			return;
		}

		char decision;
		instructions >> decision;

		if ( decision != 'A' && decision != 'R' && decision != 'M' )
		{
			std::cerr << "Error in instructions for Debugging Annealer: invalid decision" << std::endl;
			std::cerr << "Requires decision to be either 'A', or 'R' (for accept or reject)" << std::endl;
			std::cerr << "decision: " << decision << std::endl;
			return;
		}

		float deltaE( 0.0f ), previous_energy_for_node( 0.0f );
		ig_->consider_substitution( node_to_change, new_state_for_node, deltaE, previous_energy_for_node );

		std::cout << "mres: " << node_to_change << ", state: ";
		std::cout << new_state_for_node << ", deltaE: " << deltaE << std::endl;


		if ( decision == 'A' )
		{
			ig_->commit_considered_substitution();
			current_state( node_to_change ) = new_state_for_node;
		}

	}

	best_state_on_node = current_state;

	//convert best_state_on_node into best_rotamer_at_seqpos
	for (int ii = 1;ii <= nmoltenres; ++ii){
		int iiresid = p_rotamer_set_->moltenres_2_resid(ii);
		bestrotamer_at_seqpos_(iiresid) = best_state_on_node(ii) + p_rotamer_set_->rotindex_offsets(iiresid);
	}

}

}//end namespace annealer
}//end namespace pack
}//end namespace core
