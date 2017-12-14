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
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <basic/Tracer.hh>

#include <utility>
#include <utility/exit.hh>
#include <numeric/random/random.hh>

#include <iostream>

#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

using namespace ObjexxFCL;

namespace core {
namespace pack {
namespace annealer {

static basic::Tracer TR( "core.pack.annealer.FixbbLinkingRotamerSimAnnealer", basic::t_info );

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
/// constructor
///
/// @details
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author


////////////////////////////////////////////////////////////////////////////////
FixbbLinkingRotamerSimAnnealer::FixbbLinkingRotamerSimAnnealer(
	utility::vector0<int> & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
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
	), ig_(std::move(ig))
{
	setup_rotamer_links( rotamer_links );
}

FixbbLinkingRotamerSimAnnealer::FixbbLinkingRotamerSimAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	core::PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
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
	rotamer_links_ = RotamerLinksOP( new rotamer_set::RotamerLinks() );
	rotamer_links_->resize( rotamer_sets()->nmoltenres() );

	TR.Debug << "nmoltenres " << rotamer_sets()->nmoltenres()<< std::endl;

	for ( Size moltenres_id=1; moltenres_id<= rotamer_sets()->nmoltenres(); ++moltenres_id ) {
		uint const resid( rotamer_sets()->moltenres_2_resid( moltenres_id ) );
		//init anything linking to it locally
		//iterate over the associated set to check the positions, if molten

		if ( ! rotamer_links->has(resid) ) {
			TR.Debug << "singular unlinked position" << std::endl;
		} else {
			utility::vector1<int> copies = rotamer_links->get_equiv(resid);
			TR.Debug << "copies" << copies.size() << std::endl;
			for ( Size i = 1; i <= copies.size(); ++i ) {
				if ( rotamer_sets()->resid_2_moltenres( copies[i] ) ) {
					TR.Debug << "setting equivalent" << moltenres_id << "==" << rotamer_sets()->resid_2_moltenres(copies[i]) << std::endl;
					rotamer_links_->set_equiv(moltenres_id, rotamer_sets()->resid_2_moltenres(copies[i]));
				} else {
					TR.Debug << "a position in the link isn't to be changed" << std::endl;
				}
			}
		}
	}
}


/// @brief virtual destructor
FixbbLinkingRotamerSimAnnealer::~FixbbLinkingRotamerSimAnnealer() = default;


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

	/*// Detect the quasisymmetrical case by checking to see that there are
	// residues with only 1 link to themselves as well as residues with multiple
	// links.
	//
	// Repeat proteins are also qua
	bool flag1 = false; bool flag2 = false;
	for ( core::Size i = 1; i <= nmoltenres; ++i ) {

	utility::vector1<Size> these_links = rotamer_links_->get_equiv(i);
	if ( flag1 && flag2 ) break;
	if ( these_links.size() == 1 && these_links[1] == i ) {
	flag1 = true;
	TR.Debug << "quasisymmetric annealer flag 1: ON" << std::endl;
	}
	if ( these_links.size() > 1 ) {
	flag2 = true;
	TR.Debug << "quasisymmetric annealer flag 2: ON" << std::endl;
	}
	}*/

	bool quasiflag = false;
	//this quasisymmetry flag turns on a lot of bypasses below
	if ( basic::options::option[ basic::options::OptionKeys::packing::quasisymmetry]() == true ) {
		quasiflag = true;
		TR << "NOTICE: QUASISYMMETRIC PACKING IS TURNED ON in FixbbLinkingRotamerSimAnnealer. (quasiflag = " << quasiflag  << ")" << std::endl;
	}

	core::Size totalrot = 0;
	// totalrot needs to be calculated differently for the quasisymmetrical case
	if ( quasiflag ) {
		for ( core::Size res=1; res<=nmoltenres; ++res ) {
			totalrot += rotamer_sets()->nrotamers_for_moltenres( res );

		}
		TR << "QUASIBYPASS: quasisymmetric totalrot = " << totalrot << std::endl;
	} else {

		//experimental
		utility::vector1<Size> segmentTest = rotamer_links_->get_equiv(nmoltenres);

		for ( Size ii=1; ii<=segmentTest.size(); ++ii ) {
			// get the first element of the last repeat.  it should be segment length why? :: Bad logic.
			//std::cout<< "SEGMENTLENGTH from ROTAMER LINK" << segmentTest[1] << std::endl;
			//Size repeat_number = segmentTest.back()/segmentTest[1];
			//std::cout<< "number of repeats" << repeat_number << std::endl;
			for ( core::Size res = segmentTest[1]; res <= segmentTest[1]*2 ; res++ ) {
				totalrot += rotamer_sets()->nrotamers_for_moltenres(res);
			}
		}
		//std::cout << "TOTAL ROTAMER " << totalrot << std::endl;
	} // end calculate totalrot


	//setup_iterations();
	//setup_iterations(num_rots_to_pack()/repeat_number);
	setup_iterations(totalrot*2);

	FArray1D< core::PackerEnergy > previous_nsteps_for_rot( rotamer_sets()->nrotamers(), 0.0);

	int outeriterations = get_outeriterations();

	// some rotamer may not exist on other repeats, and use a new vector to
	// iterate the "good" rotamers
	int allrot = rotamer_sets()->nrotamers();
	utility::vector1<bool> rot_valid(allrot, true);

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		setup_temperature(loopenergy,nn);
		if ( quench() ) {
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
		for ( int n = 1; n <= inneriterations; ++n ) {

			int ranrotamer = -1;
			bool invalid_rotamer = false;
			while ( !invalid_rotamer ) {
				ranrotamer = static_cast<int>( numeric::random::rg().random_range(1, allrot ));
				if ( rot_valid[ ranrotamer ] ) {
					invalid_rotamer = true;
				}
			}

			//int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int const rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

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

			for ( auto itr = linked_residues.begin(), ite = linked_residues.end(); itr != ite; ++itr ) { // go through each linked residue
				num_linked_res++;
				TR.Debug << "analyzing num_linked_res: " << num_linked_res << std::endl;
				if ( ( *itr != 0 ) && ( *itr != moltenres_id ) ) { //skip this step if residue is self
					//try multiple substitutions
					if ( TR.Trace.visible() ) {
						TR.Trace << "moltenres_id " << moltenres_id << " coupled to moltenres_id " << *itr << std::endl;
					}

					other_prevrotamer_state = state_on_node(*itr);

					//pick a rotamer at a linked position
					RotamerSetCOP other_rotamer_set( rotamer_sets()->rotamer_set_for_moltenresidue( *itr ) );
					//ResidueCOP other_rotamer( other_prevrotamer_state == 0 ? ResidueCOP(0) : other_rotamer_set->rotamer( other_prevrotamer_state ) );
					ResidueCOP other_rotamer( ResidueCOP(nullptr) );

					utility::vector1<int> passed_states;

					int const other_nrotamers( other_rotamer_set->num_rotamers() );
					int tries = other_nrotamers;
					found_rotamer = false;
					while ( tries ) {

						// pick a rotamer at the other position
						int other_rotamer_state = tries;

						other_rotamer = other_rotamer_set->rotamer(other_rotamer_state);
						--tries;

						// For quasisymmetric case, check for the same AA at linked positions,
						// but not the same rotamer.
						if ( quasiflag ) {
							if ( new_rotamer->is_similar_aa( *other_rotamer ) ) {
								TR.Debug << "QUASIBYPASS: similar AAs found: " << moltenres_id << " (" << new_rotamer->aa() << ") and "<< *itr << " (" << other_rotamer->aa() << ")";
								if ( new_rotamer->is_similar_rotamer( *other_rotamer ) ) { //found the same rotamer
									TR.Debug << " [IDENTICAL ROTAMERS]" << std::endl;
								} else {
									TR.Debug << std::endl;
								}
								//generate new list of all AAs that pass
								passed_states.push_back( other_rotamer_state );
								continue;
							} else {
								TR.Debug << "QUASIBYPASS: AAs not similar, tries remaining: " << tries << std::endl;
							}

						} else { // not quasisymmetric
							if ( new_rotamer->is_similar_rotamer( *other_rotamer ) ) { //found the same rotamer, move on
								//if ( new_rotamer->is_similar_aa( *other_rotamer ) ) { //found the same rotamer, move on
								//std::cout << "found the same rotamer for " << moltenres_id << " and " <<  *itr << "of types " << new_rotamer->aa() << " and " << other_rotamer->aa() << std::endl;
								found_rotamer = true;
								// record the state
								resid_states[*itr] = other_rotamer_state;

								break;
							}
						}
					} //tries

					if ( quasiflag ) { //if quasisymmetric case, RNG pick a rotamer from the newly compiled list of similar AAs
						TR.Debug << "QUASIBYPASS: number of states with similar AA: " << passed_states.size() << std::endl;
						auto ranrotamer2 = static_cast<int>( numeric::random::rg().random_range(1, passed_states.size() ));
						found_rotamer = true; //flags this rotamer to be "same" for the sake of code downstream
						resid_states[*itr] = passed_states[ ranrotamer2 ]; //record state of "other rotamer"
					}

					if ( !found_rotamer ) { // any of the linked position without the same rotamer should be passed
						TR.Debug << "same rotamer not found for " << moltenres_id << " and " <<  *itr << std::endl;
						break;
					}

				} else if ( (*itr == moltenres_id) && (linked_residues.size() == 1) ) {
					// For quasisymmetrical design, one might have designable positions that are not
					// quasisymmetrical, and therefore do not need RotamerLinks to other residues. But in
					// order to get around a segfault in src/core/pack/rotamer_sets/RotamerSets.cc,
					// these residues must have RotamerLinks to only themselves; they therefore fail the
					// if above, leading to a continue directly below at if (!foundrotamer). So, here we
					// will detect these positions and set found_rotamer to true.
					TR.Debug << "QUASIBYPASS: SELF-LINKED AA found: " << moltenres_id << " (" << new_rotamer->aa() << ")" << std::endl;
					found_rotamer = true;
					other_prevrotamer_state = state_on_node(*itr);
				}

			} //for linked residues

			if ( ( !found_rotamer ) && ( ! ( quasiflag ) ) ) { // any of the linked position without the same rotamer should be passed (don't do this for quasisymmetry)
				//invalidate all the linked positions
				TR.Debug << "invalidate " ;
				for ( auto & resid_state : resid_states ) {
					rot_valid[ rotamer_sets()->moltenres_rotid_2_rotid( resid_state.first, resid_state.second ) ] = false;
					TR.Debug << resid_state.first << "(" << resid_state.second << ")"  ;
				}
				TR.Debug << std::endl;
				continue;
			}

			//score the good rotamers and pass through metropolis

			core::PackerEnergy totalenergy = 0.0; //reset totalenergy
			for ( auto & resid_state : resid_states ) {
				core::PackerEnergy delta_energy_temp( 0.0 ), previous_energy_for_node_temp( 0.0 ); //initialize to zero?
				ig_->consider_substitution( resid_state.first, resid_state.second, delta_energy_temp, previous_energy_for_node_temp );

				currentenergy = ig_->commit_considered_substitution();
				TR.Debug << "current energy(" << resid_state << "): " << currentenergy << std::endl;
				totalenergy += currentenergy;
				TR.Debug << "total energy: " << totalenergy << std::endl;
				delta_energy_accumulated += delta_energy_temp;
				previous_energy_for_node_accumulated += previous_energy_for_node_temp;
			}
			core::PackerEnergy avgenergy = ( totalenergy / resid_states.size() ); //calculate avgenergy
			TR.Debug << "average energy: " << totalenergy << "/" << resid_states.size() << " = " << avgenergy << std::endl;

			core::PackerEnergy previous_energy_average = ( previous_energy_for_node + previous_energy_for_node_accumulated );
			core::PackerEnergy delta_energy_average = ( delta_energy + delta_energy_accumulated );

			if ( prevrotamer_state == 0 || other_prevrotamer_state == 0 || pass_metropolis( previous_energy_average, delta_energy_average ) ) {

				// accept !!!!!!!
				if ( TR.Trace.visible() ) {
					TR.Trace << "accepting multiple rotamer substitution (pass_metropolis)" << std::endl;
				}

				//set state
				for ( auto & resid_state : resid_states ) {
					TR.Debug << "resid_states " << resid_state << " first(residue): " << resid_state.first << " second(rotamer): " << resid_state.second << std::endl;
					state_on_node( resid_state.first ) = resid_state.second;
				}

				TR.Debug << "current (avg) energy: " << avgenergy << " best energy: " << bestenergy() << std::endl;
				//if ( ( prevrotamer_state == 0 ) || ( other_prevrotamer_state == 0 ) || ( currentenergy <= bestenergy() ) ) { //prevrotamerstate == 0 means this position has not taken a new rotamer before
				if ( ( prevrotamer_state == 0 ) || ( other_prevrotamer_state == 0 ) || ( avgenergy <= bestenergy() ) ) { //prevrotamerstate == 0 means this position has not taken a new rotamer before
					TR.Debug << "accepted rotamer stored in best (linked)" << std::endl;
					//bestenergy() = currentenergy;
					bestenergy() = avgenergy;
					best_state_on_node = state_on_node;
				}

			} else {
				// reject
				if ( TR.Trace.visible() ) {
					TR.Trace << "rejecting multiple rotamer substitution (fail_metropolis)" << std::endl;
				}
				//revert changes:

				for ( auto & resid_state : resid_states ) {

					core::PackerEnergy dE, oldE;
					ig_->consider_substitution( resid_state.first,  state_on_node(resid_state.first),
						dE, oldE );
					currentenergy = ig_->commit_considered_substitution();
				}
			} // accept or reject?

			loopenergy(nn) = currentenergy;

			debug_assert( !calc_rot_freq() );
			// skip the logic below for single-rotamer substitution ////////////////////////////////////
			continue;

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				TR.Debug << "entering bk outer loop" << std::endl;
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					TR.Debug << "entering bk inner loop" << std::endl;
					bestenergy() = currentenergy;
					best_state_on_node = state_on_node;
				}
			} // accept

			loopenergy(nn) = currentenergy;
			core::PackerEnergy const temperature = get_temperature();

			if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
				++nsteps;
				for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
					int iistate = state_on_node(ii);
					if ( iistate != 0 ) {
						++nsteps_for_rot( rotamer_sets()->moltenres_rotid_2_rotid(ii, iistate) );
					}
				}
			}
		} // end of inneriteration loop
	} //end of outeriteration loop

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In FixbbLinkingRotamerSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
		std::cerr << "Critical error -- assignment and energy of assignment meaningless" << std::endl;

		FArray1D_int nstates_for_moltenres( rotamer_sets()->nmoltenres(), 0 );
		for ( uint ii = 0; ii < num_rots_to_pack(); ++ii ) {
			++nstates_for_moltenres( rotamer_sets()->res_for_rotamer( rot_to_pack()[ ii ] ) );
		}

		for ( uint ii = 1; ii <= rotamer_sets()->nmoltenres(); ++ii ) {
			if ( best_state_on_node( ii ) == 0 ) {
				std::cout << "Molten res " << ii << " (residue " << rotamer_sets()->moltenres_2_resid( ii );
				std::cout << " ) assigned state 0 despite having " << nstates_for_moltenres( ii ) << " states to choose from" << std::endl;
			}
		}
		debug_assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();
	}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
		Size const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}
}

}//end namespace annealer
}//end namespace pack
}//end namespace core
