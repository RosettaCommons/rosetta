// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbLinkingRotamerSimAnnealer.cc
/// @brief  Packer's standard simulated annealing class implementation
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/FixbbLinkingRotamerSimAnnealer.hh>

// Package Headers
#include <core/pack/rotamer_set/RotamerSets.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerLinks.hh>
#include <core/pack/interaction_graph/InteractionGraphBase.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <basic/Tracer.hh>

//#include "after_opts.h"
//#include "FixbbLinkingRotamerSimAnnealer.h"
//#include "RotamerAssigningAnnealer.h"
//#include "random_numbers.h"
//#include "param.h"
//#include "RotamerSet.h"

#include <utility/exit.hh>
#include <numeric/random/random.hh>

// AUTO-REMOVED #include <ObjexxFCL/Fmath.hh>

#include <iostream>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

static thread_local basic::Tracer TR( "core.pack.annealer.FixbbLinkingRotamerSimAnnealer", basic::t_info );

////////////////////////////////////////////////////////////////////////////////
/// @begin FixbbLinkingRotamerSimAnnealer::FixbbLinkingRotamerSimAnnealer()
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
FixbbLinkingRotamerSimAnnealer::FixbbLinkingRotamerSimAnnealer(
	utility::vector0<int> & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::InteractionGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq,
	RotamerLinksCOP rotamer_links
):
	RotamerAssigningAnnealer(
	rot_to_pack,
	(int) rot_to_pack.size(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	rotamer_sets,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
	setup_rotamer_links( rotamer_links );
}

FixbbLinkingRotamerSimAnnealer::FixbbLinkingRotamerSimAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::InteractionGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D< core::PackerEnergy > & rot_freq,
	RotamerLinksCOP rotamer_links
):
	RotamerAssigningAnnealer(
	(ig->get_num_total_states()),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	rotamer_set,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	), ig_(ig)
{
	setup_rotamer_links( rotamer_links );
}

// rotamer links provide a map for each position and all the positions linking
// together, including itself.

void
FixbbLinkingRotamerSimAnnealer::setup_rotamer_links(
	RotamerLinksCOP rotamer_links
)
{
	// setup the rotamer links
	// for each moltenres
	rotamer_links_ = new rotamer_set::RotamerLinks();
	rotamer_links_->resize( rotamer_sets()->nmoltenres() );

	//std::cout << "nmoltenres " << rotamer_sets()->nmoltenres()<< std::endl;

	for ( Size moltenres_id=1; moltenres_id<= rotamer_sets()->nmoltenres(); ++moltenres_id ) {
		uint const resid( rotamer_sets()->moltenres_2_resid( moltenres_id ) );
		//std::cout << "RESID: " << resid << std::endl;
		//init anything linking to it locally
		//iterate over the associated set to check the positions, if molten
		//std::cout << "linking positions: ";
		if (rotamer_links->has(resid)){ // not a null
		//std::cout << " linked IN " << std::endl;
			utility::vector1<int> copies = rotamer_links->get_equiv(resid);
			for (Size i = 1; i <= copies.size(); ++i ){
				if ( rotamer_sets()->resid_2_moltenres( copies[i] )) {
	//				rotamer_links_->set_equiv(resid, copies[i]);
					rotamer_links_->set_equiv(moltenres_id, rotamer_sets()->resid_2_moltenres(copies[i]));
					//std::cout << copies[i];
					//std::cout << ":" << rotamer_sets()->resid_2_moltenres(copies[i]) << std::endl;
				}
				else {
					std::cout << "a position in the link isn't to be changed" << std::endl;
				}
			}
		}
		else {
			std::cout << "singular unlinked position" << std::endl;
		}
	}
}


/// @brief virtual destructor
FixbbLinkingRotamerSimAnnealer::~FixbbLinkingRotamerSimAnnealer()
{}


void FixbbLinkingRotamerSimAnnealer::run()
{
	using namespace core::conformation;

	core::Size const nmoltenres = ig_->get_num_nodes();

	FArray1D_int state_on_node( nmoltenres, 0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres, 0 );
	FArray1D< core::PackerEnergy > loopenergy(maxouteriterations,0.0);

	//bk variables for calculating rotamer frequencies during simulation
	core::Size nsteps = 0;
	FArray1D_int nsteps_for_rot( ig_->get_num_total_states(), 0 );

	//--------------------------------------------------------------------
	//initialize variables

	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	// Detect the quasisymmetrical case by checking to see that there are
	// residues with only 1 link to themselves as well as residues with multiple
	// links.
	bool flag1 = false; bool flag2 = false;
	for ( core::Size i = 1; i <= nmoltenres; ++i ) {
		utility::vector1<Size> these_links = rotamer_links_->get_equiv(i);
		if ( flag1 && flag2 ) break;
		if ( these_links.size() == 1 && these_links[1] == i )
			flag1 = true;
		if ( these_links.size() > 1 )
			flag2 = true;
	}

	core::Size totalrot = 0;
	// totalrot needs to be calculated differently for the quasisymmetrical case
	if ( flag1 && flag2 ) {
		for ( core::Size res=1; res<=nmoltenres; ++res ) {
			totalrot += rotamer_sets()->nrotamers_for_moltenres( res );
		}
	} else {

		//experimental
		utility::vector1<Size> segmentTest = rotamer_links_->get_equiv(nmoltenres);
		// get the first element of the last repeat.  it should be segment length
		//std::cout<< "SEGMENTLENGTH from ROTAMER LINK" << segmentTest[1] << std::endl;
		//Size repeat_number = segmentTest.back()/segmentTest[1];
		//std::cout<< "number of repeats" << repeat_number << std::endl;

		for (core::Size res = segmentTest[1]; res <= segmentTest[1]*2 ; res++){
			totalrot += rotamer_sets()->nrotamers_for_moltenres(res);
		}
			//std::cout << "TOTAL ROTAMER " << totalrot << std::endl;
	} // end quasisymmetric if-else


	//setup_iterations();
	//setup_iterations(num_rots_to_pack()/repeat_number);
	setup_iterations(totalrot*2);

	FArray1D< core::PackerEnergy > previous_nsteps_for_rot( rotamer_sets()->nrotamers(), 0.0);

	int outeriterations = get_outeriterations();


	//std::ofstream annealer_trajectory;
	//static bool const record_annealer_trajectory( truefalseoption("record_annealer_trajectory") ); //look up once
	//if ( record_annealer_trajectory )
	//{
	//	std::string trajectory_file_name( stringafteroption("record_annealer_trajectory" ) );
	//	annealer_trajectory.open(trajectory_file_name.c_str() );
	//}


			// some rotamer may not exist on other repeats, and use a new vector to
			// iterate the "good" rotamers


		int allrot = rotamer_sets()->nrotamers();
		utility::vector1<bool> rot_valid(allrot, true);
		//std::cout << "outer iteration: " << outeriterations << std::endl;


	//outer loop
	for (int nn = 1; nn <= outeriterations; ++nn ){
		setup_temperature(loopenergy,nn);
		if ( quench() ){
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}
		//rh std::cout << "Sim Annealer Temperature: " << get_temperature() << std::endl;

		int inneriterations = get_inneriterations();
		//std::cout << "inner iteration: " << inneriterations << std::endl;

		core::PackerEnergy treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//inner loop
		for (int n = 1; n <= inneriterations; ++n ){

			int ranrotamer = -1;
			bool invalid_rotamer = false;
			while (!invalid_rotamer){
				ranrotamer = static_cast<int>( numeric::random::rg().random_range(1, allrot ));
				if (rot_valid[ ranrotamer ]){
					invalid_rotamer = true;
				}
			}

			//int const ranrotamer = pick_a_rotamer( n );
			if (ranrotamer == -1) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int const rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if (rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			core::PackerEnergy previous_energy_for_node, delta_energy;

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
																	delta_energy, previous_energy_for_node);

			// specialize to the case of coupled pairs in this first pass implementation

			// assume all couplings are between moltenres -- couplings between fixed and moltenres could
			// have been used to trim the residueset
			utility::vector1<int> linked_residues = rotamer_links_->get_equiv(moltenres_id);

			RotamerSetCOP rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( moltenres_id ) );
			ResidueCOP new_rotamer( rotamer_set->rotamer( rotamer_state_on_moltenres ) );

			//core::PackerEnergy tmp_currentenergy;

			//setup energy to keep track of changes
			core::PackerEnergy delta_energy_accumulated=0, previous_energy_for_node_accumulated=0;

			std::map<Size, Size> resid_states;

			//record the seeding position
			resid_states[moltenres_id]  = rotamer_state_on_moltenres;

			int other_prevrotamer_state(0);

			bool found_rotamer = false;
			Size num_linked_res =0;
			for (utility::vector1<int>::iterator itr = linked_residues.begin(), ite = linked_residues.end(); itr != ite; itr++ ){
				num_linked_res++;
				if ( (*itr != 0) && (*itr != moltenres_id )){
					TR.Trace << "moltenres_id " << moltenres_id << " coupled to moltenres_id " << *itr << std::endl;

					//try multiple substitutions

					TR.Trace << "Picked rotamer incompatible, trying multiple substitution" << std::endl;

					other_prevrotamer_state = state_on_node(*itr);

					//pick a rotamer at a linked position
					RotamerSetCOP other_rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( *itr ) );
					ResidueCOP other_rotamer( other_prevrotamer_state == 0 ? ResidueCOP(0) : other_rotamer_set->rotamer( other_prevrotamer_state ) );

					int other_rotamer_state(0);
					int const other_nrotamers( other_rotamer_set->num_rotamers() );
					int tries = other_nrotamers;
					 found_rotamer = false;
					while ( tries ) {
						// pick a rotamer at the other position
						other_rotamer_state = tries;

						other_rotamer = other_rotamer_set->rotamer(other_rotamer_state);
						--tries;

						// For quasisymmetric case, check for the same AA at linked positions,
						// but not the same rotamer.
						if ( flag1 && flag2 ) {
							if ( new_rotamer->is_similar_aa( *other_rotamer ) ) {
								found_rotamer = true;
								resid_states[*itr] = other_rotamer_state;
								break;
							}
						} else { // not psuedosymmetric
							if ( new_rotamer->is_similar_rotamer( *other_rotamer ) ) { //found the same rotamer, move on
							//if ( new_rotamer->is_similar_aa( *other_rotamer ) ) { //found the same rotamer, move on
							//std::cout << "found the same rotamer for " << moltenres_id << " and " <<  *itr << "of types " << new_rotamer->aa() << " and " << other_rotamer->aa() << std::endl;
								found_rotamer = true;

								// record the state
								resid_states[*itr] = other_rotamer_state;

								break;
							}
						}
					}
					if (!found_rotamer){ // any of the linked position without the same rotamer should be passed
						//std::cout << "same rotamer not found for " << moltenres_id << " and " <<  *itr << std::endl;
						break;
					}

				}
        // For quasisymmetrical design, one might have designable positions that are not
        // quasisymmetrical, and therefore do not need RotamerLinks to other residues. But in
        // order to get around a segfault in src/core/pack/rotamer_sets/RotamerSets.cc,
        // these residues must have RotamerLinks to only themselves; they therefore fail the
        // if above, leading to a continue directly below at if (!foundrotamer). So, here we
        // will detect these positions and set found_rotamer to true.
        else if ( (*itr == moltenres_id) && (linked_residues.size() == 1) ) {
          found_rotamer = true;
					other_prevrotamer_state = state_on_node(*itr);
        }
			} // for linked residues

			if (!found_rotamer){ // any of the linked position without the same rotamer should be passed
				//invalidate all the linked positions
				//std::cout << "invalidate " ;
				for (std::map<Size, Size>::iterator it = resid_states.begin(), ite = resid_states.end(); it != ite; it++){
					rot_valid[ rotamer_sets()->moltenres_rotid_2_rotid( (*it).first, (*it).second ) ] = false;
					//	std::cout << (*it).first << "(" << (*it).second << ")"  ;
				}
				//std::cout << std::endl;
				continue;

			}

			//score the good rotamers and pass through metropolis
			//std::cout << "summing energies: " ;

			for (std::map<Size, Size>::iterator it = resid_states.begin(), ite = resid_states.end(); it != ite; it++){

			//std::cout << (*it).first << "(" << (*it).second << ")"  ;

				core::PackerEnergy delta_energy_temp, previous_energy_for_node_temp;
				ig_->consider_substitution( (*it).first,  (*it).second,
																		delta_energy_temp, previous_energy_for_node_temp );

				//tmp_currentenergy = currentenergy;  // set but never used ~Labonte
				currentenergy = ig_->commit_considered_substitution();

				{ // debugging
						//Real const dev( std::abs( currentenergy - tmp_currentenergy - delta_energy_temp ) );
			//std::cout << (*it).first << "(" << (*it).second << ")"  ;
				//if ( dev > 0.01 ) {
				//	std::cout  << "equal2? " << dev << ' ' << currentenergy << " " <<  tmp_currentenergy << ' ' << delta_energy_temp <<
				//	'\n';
				//}
				} // scope


				delta_energy_accumulated += delta_energy_temp;
				previous_energy_for_node_accumulated += previous_energy_for_node_temp;

			//	std::cout << "delta_energy_temp " << delta_energy_temp << " previous_energy_for_node_temp " << previous_energy_for_node_temp ;

			}
			//std::cout << " accumulated " << previous_energy_for_node_accumulated << "  delta_accumulated" << delta_energy_accumulated;
			//std::cout << std::endl;


			//core::PackerEnergy previous_energy_average = ( previous_energy_for_node + previous_energy_for_node_accumulated )/(num_linked_res+1);
			//core::PackerEnergy delta_energy_average = ( delta_energy + delta_energy_accumulated )/(num_linked_res+1);
			core::PackerEnergy previous_energy_average = ( previous_energy_for_node + previous_energy_for_node_accumulated );
			core::PackerEnergy delta_energy_average = ( delta_energy + delta_energy_accumulated );

			//std::cout << prevrotamer_state << " " << other_prevrotamer_state << " " <<  previous_energy_average << " " << delta_energy_average << std::endl;

			if ( prevrotamer_state == 0 || other_prevrotamer_state == 0 ||
						pass_metropolis( previous_energy_average, delta_energy_average ) ) {
				// accept !!!!!!!
				TR.Trace << "accepting multiple rotamer substitution" << std::endl;


			//std::cout << "ACCEPT: substitution on ";
			//set state
				for (std::map<Size, Size>::iterator it = resid_states.begin(), ite = resid_states.end(); it != ite; it++){
	/*				core::PackerEnergy dE, oldE;
					ig_->consider_substitution( (*it).first,  (*it).second,
																			dE, oldE );
																			tmp_currentenergy = currentenergy;
					currentenergy = ig_->commit_considered_substitution();

					{ // debugging
						Real const dev( std::abs( currentenergy - tmp_currentenergy - dE ) );
				std::cout << (*it).first << "(" << (*it).second << ")"  ;
				//		if ( dev > 0.01 ) {
							std::cout  << "equal2? " << dev << ' ' << currentenergy << " " <<  tmp_currentenergy << ' ' << dE <<
								'\n';
				//		}
					} // scope

					std::cout << "RUNNING CURRENT " << currentenergy << std::endl;
				//	std::cout << (*it).first << "(" << (*it).second << ")"  ;
*/
					state_on_node( (*it).first ) = (*it).second;
				}
			//std::cout << "ACCEPT: currentenergy " << currentenergy << std::endl;
			//std::cout << std::endl;

			//std::cout << "CURRENT: " << currentenergy << "BEST: " << bestenergy() << std::endl;
				if ( ( prevrotamer_state == 0 ) || ( other_prevrotamer_state == 0 ) || ( currentenergy <= bestenergy() )) {
					bestenergy() = currentenergy;
					best_state_on_node = state_on_node;
					if ( false ) { // hacking ///////////////////////////////////////////
						std::cout << "best-accept: ";
						for ( Size i=1; i<= Size(nmoltenres); ++i ) {
							if ( state_on_node( i ) == 0 ) {
								std::cout << '.';
							} else {
								RotamerSetCOP rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( i ) );
								conformation::ResidueCOP rotamer( rotamer_set->rotamer( state_on_node( i ) ) );
								if ( rotamer->is_DNA() ) TR << rotamer->name1();
							}
						}
						std::cout << ' ' << nn << ' ' << n << ' ' << currentenergy << '\n';
					} // end hacking ///////////////////////////////////////
				}

			} else {
				// reject
				TR.Trace << "rejecting multiple rotamer substitution" << std::endl;
				//revert changes:

				for (std::map<Size, Size>::iterator it = resid_states.begin(), ite = resid_states.end(); it != ite; it++){

					core::PackerEnergy dE, oldE;
					ig_->consider_substitution( (*it).first,  state_on_node((*it).first),
																			dE, oldE );
					//tmp_currentenergy = currentenergy;  // set but never used ~Labonte
					currentenergy = ig_->commit_considered_substitution();
				}
			//std::cout << "REJECT: currentenergy " << currentenergy << std::endl;
			//	{ // debugging
			//		Real const dev( std::abs( tmp_currentenergy - currentenergy ) );
			//		if ( dev > 0.01 ) {
			//			TR << "equal3? " << dev << ' ' << tmp_currentenergy << ' ' << currentenergy << '\n';
			//		}
			//	} // scope
			//	currentenergy = previous_energy_average + delta_energy_average;
			} // accept or reject?

			loopenergy(nn) = currentenergy;

			assert( !calc_rot_freq() );
			continue; // skip the logic below for single-rotamer substitution ////////////////////////////////////




			//std::cerr << "mres: " << moltenres_id << ", state: ";
			//std::cerr << rotamer_state_on_moltenres << ", deltaE: " << delta_energy << std::endl;

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) )
			{

				//std::cerr << "Accepted" << std::endl;
//				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ((prevrotamer_state == 0)||(currentenergy < bestenergy() ))
				{
					bestenergy() = currentenergy;
					best_state_on_node = state_on_node;
				}
			} // accept

			loopenergy(nn) = currentenergy;
			core::PackerEnergy const temperature = get_temperature();

			if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) )
			{
				++nsteps;
				for (Size ii = 1; ii <= nmoltenres; ++ii )
				{
					int iistate = state_on_node(ii);
					if (iistate != 0)
					{
						++nsteps_for_rot( rotamer_sets()->moltenres_rotid_2_rotid(ii, iistate) );
					}
				}
			}


		} // end of inneriteration loop
	} //end of outeriteration loop

	if ( ig_->any_vertex_state_unassigned() )
	{
		std::cerr << "Critical error -- In FixbbLinkingRotamerSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
		std::cerr << "Critical error -- assignment and energy of assignment meaningless" << std::endl;

		FArray1D_int nstates_for_moltenres( rotamer_sets()->nmoltenres(), 0 );
		for ( uint ii = 0; ii < num_rots_to_pack(); ++ii)
		{
			++nstates_for_moltenres( rotamer_sets()->res_for_rotamer( rot_to_pack()[ ii ] ) );
		}

		for ( uint ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii)
		{
			if ( best_state_on_node( ii ) == 0 )
			{
				std::cout << "Molten res " << ii << " (residue " << rotamer_sets()->moltenres_2_resid( ii );
				std::cout << " ) assigned state 0 despite having " << nstates_for_moltenres( ii ) << " states to choose from" << std::endl;
			}
		}
		assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();


	}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for (Size ii = 1; ii <= nmoltenres; ++ii){
		Size const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}


}

}//end namespace annealer
}//end namespace pack
}//end namespace core
