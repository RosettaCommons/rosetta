// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/flexpack/annealer/FlexbbSimAnnealer.cc
/// @brief  Annealer for optimizing both backbone and sidechain conformations implementation.
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

/// Unit Headers
#include <protocols/flexpack/annealer/FlexbbSimAnnealer.hh>

/// Package headers
#include <protocols/flexpack/rotamer_set/FlexbbRotamerSets.hh>
#include <protocols/flexpack/interaction_graph/FlexbbInteractionGraph.hh>

/// Project headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/flexpack.OptionKeys.gen.hh>

/// ObjexxFCL headers

/// Utility headers
#include <utility/exit.hh>

/// Numeric headers
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

/// C++ headers
#include <iostream>

#include <utility/vector0.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {
namespace annealer {


static THREAD_LOCAL basic::Tracer TR( "protocols.flexpack.annealer.FlexbbSimAnnealer" );

FlexbbSimAnnealer::FlexbbSimAnnealer(
	ObjexxFCL::FArray1D_int & bestrotamer_at_seqpos,
	PackerEnergy & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::FlexbbInteractionGraphOP ig,
	rotamer_set::FlexbbRotamerSetsCOP rotsets,
	ObjexxFCL::FArray1D_int & current_rot_index,
	bool calc_rot_freq,
	ObjexxFCL::FArray1D< PackerEnergy > & rot_freq
):
	parent(
	ig->get_num_total_states(),
	bestrotamer_at_seqpos,
	bestenergy,
	start_with_current, // start simulation with current rotamers
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	ig_(ig),
	rotsets_( rotsets )
{}

FlexbbSimAnnealer:: ~FlexbbSimAnnealer()
{}

/// @details The FlexbbSimAnnealer operates with three submodes:
/// mvsc_only :    same as the fixbb mode, only consider moving the side chains
/// mvbb_only :    no movement with the side chain, but only move the backbone
/// both_sc_and bb:default mode would be a combination of move side chain and
///                backbone fragments
void FlexbbSimAnnealer::run()
{
	using namespace interaction_graph;


	TR << "Beginning flexible backbone simulated annealing" << std::endl;

	//--------------------------------------------------------------------

	Size ranrotamer,ranbbfrag;//,rannum;
	Size nmoltenres = ig_->get_num_nodes();
	//Size num_accessible, num_accessible_bb, num_backbones_available;
	Size moltenres_id, rotamer_state_on_moltenres, prevrotamer_state;
	Size const nrotamers( rotsets_->nrotamers() );
	Size cycle1, cycle2, cycle3;
	bool valid_bb_move; //apl if a node on a flexible fragment is in state 0 no meaningful "backbone move" exists
	PackerEnergy currentenergy, previous_energy_for_node, delta_energy;
	PackerEnergy previous_energy_for_bbfrag;

	utility::vector1< Size > accessible_state_list( nrotamers, 0 );
	utility::vector1< Size > accessible_state_list_bbfrag( rotsets_->nbackbone_conformations(), 0 );//max number of backbone fragment

	ObjexxFCL::FArray1D< PackerEnergy > loopenergy( maxouteriterations, 0.0);

	//bk variables for calculating rotamer frequencies during simulation
	Size nsteps = 0;
	ObjexxFCL::FArray1D_int nsteps_for_rot( nrotamers, 0 );

	ObjexxFCL::FArray1D_int state_on_node( nmoltenres, 0);
	ObjexxFCL::FArray1D_int best_state_on_node( nmoltenres, 0);
	utility::vector1< int > moltenres_rotoffsets( nmoltenres, 0 );
	for ( Size ii = 1; ii <= nmoltenres; ++ii ) { moltenres_rotoffsets[ ii ] = rotsets_->nrotamer_offset_for_moltenres( ii ); }

	//--------------------------------------------------------------------
	//initialize variables

	cycle1 = 130;
	cycle2 = 10;
	cycle3 = 10;

	currentenergy = 0.0;

	//variables to keep track of the frequency of rotamers
	for ( Size i = 1; i <= rotsets_->nrotamers(); ++i ) {
		accessible_state_list[ i ] = i;
		nsteps_for_rot(i) = 0;
	}

	ig_->prepare_for_simulated_annealing();

	set_lowtemp( 0.2 ); // go colder than regular SA!

	if ( start_with_current() ) {
		for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
			state_on_node(ii) = current_rot_index()( rotsets_->moltenres_2_resid( ii ));
		}
		ig_->set_network_state( state_on_node );
		currentenergy = ig_->get_energy_current_state_assignment();
		bestenergy() = currentenergy;
		best_state_on_node = state_on_node;
	} else {
		ig_->blanket_assign_state_0();
		state_on_node = 0;
	}

	//outer iterations and inner iterations
	setup_iterations();

	Size outeriterations = get_outeriterations();
	Size inneriterations_usual = get_inneriterations();

	/// Reduce the number of inner iterations to reflect the three
	/// inner-inner loops that extend the number of rotamer substitutions.
	/// Note, the base class  assigns 5x as many iterations
	/// if start_with_current is false.
	inneriterations_usual /= (start_with_current() ? 75 : 75 * 5 );

	using namespace basic::options;
	using namespace basic::options::OptionKeys::flexpack::annealer;

	if ( option[ inner_iteration_scale ].user() ) {
		inneriterations_usual *=  static_cast< Size > (inneriterations_usual * option[ inner_iteration_scale ]);
	}
	if ( option[ outer_iteration_scale ].user() ) {
		outeriterations = static_cast< Size > (outeriterations * option[ outer_iteration_scale ]);
	}
	if ( option[ fixbb_substitutions_scale ].user() ) {
		cycle1 = static_cast< Size > ( cycle1 * option[ fixbb_substitutions_scale ]);
	}
	if ( option[ pure_movebb_substitutions_scale ].user() ) {
		cycle2 = static_cast< Size > ( cycle2 * option[ pure_movebb_substitutions_scale ]);
	}
	if ( option[ rotsub_movebb_substitutions_scale ].user() ) {
		cycle3 = static_cast< Size > ( cycle3 * option[ rotsub_movebb_substitutions_scale ]);
	}

	//outeriterations *= 2;
	//inneriterations_usual /= 25;

	//collect the list of all accessible states once
	utility::vector1< Size > list_of_all( nrotamers, 0 );

	ig_->get_accessible_states(
		FlexbbInteractionGraph::BOTH_SC_AND_BB,
		list_of_all
	);


	bool sconly_move_accessible_state_list_current( false );

	//outer loop
	PackerEnergy last_temperature = 100000;
	for ( Size nn = 1; nn <= outeriterations; ++nn ) {

		Size num_fixbb_move_accepts = 0;
		Size num_fixbb_move_valid = 0;
		Size num_bb_move_accepts = 0;
		Size num_bb_move_valid = 0;
		Size num_bb_sub_accepts = 0;
		Size num_bb_sub_valid = 0;

		//set up the temperature
		setup_temperature(loopenergy,nn);
		Size inneriterations = inneriterations_usual;

		if ( nn > 1 && get_temperature() > last_temperature ) {
			/// short circuit -- do not return to high temperature.
			//assert( outeriterations > 0 ); // don't wrap to 4 billion
			//nn = outeriterations - 1; // force quence next round
			//continue;
			set_temperature( 0.5 ); // don't return to high temperatures /// ? maybe?
		}
		last_temperature = get_temperature();

		if ( get_temperature() > 5.0f ) {
			inneriterations /= 2;
		} else if ( quench() ) {
			inneriterations *= 6;
		}

		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}

		//default mode, which is a combination of move side chain, move backbone
		//need to test how many cycles for each mode
		PackerEnergy temperature = get_temperature();
		for ( Size j = 1; j <= inneriterations; ++j ) {

			if ( ! sconly_move_accessible_state_list_current ) {
				ig_->get_accessible_states( FlexbbInteractionGraph::SC_ONLY, accessible_state_list );
				sconly_move_accessible_state_list_current = true;
			}

			//     std::cout<<"num_accessible = "<<num_accessible<<"totalinnner = "<<inneriterations<<'\n'
			for ( Size n = 1; n <= cycle1; ++n ) {//move side chain only
				//pick random rotamers from the rotamer list,this subroutine should go into the base class
				//  ranrotamer = pick_a_rotamer(list,n,nn);


				ranrotamer = pick_a_rotamer( j, n, cycle1, accessible_state_list);

				if ( quench() ) {
					//std::cout<<"inneriteration = "<<j<<"cycle1 = "<<n<<"rannum="<<rannum<<"ranrot="<<ranrotamer<<'\n';
				}

				//ranrotamer = accessible_state_list( static_cast< int > (ran3() * num_accessible) + 1 );
				moltenres_id               = rotsets_->moltenres_for_rotamer( ranrotamer );
				//seqpos                     = rotsets_->moltenres_2_resid( moltenres_id );
				rotamer_state_on_moltenres = rotsets_->local_rotid_for_rotamer_on_moltenres( moltenres_id, ranrotamer );
				prevrotamer_state = state_on_node( moltenres_id );

				if ( rotamer_state_on_moltenres == prevrotamer_state ) continue;// rotamer is already assigned; skip iteration

				++num_fixbb_move_valid;
				//std::cout << "State on node: ";
				//for ( Size ii = 1; ii <= state_on_node.size(); ++ii ) { std::cout << state_on_node( ii ) << " ";}
				//std::cout << std::endl;
				//std::cout << "Consider fixed backbone rotamer substitution" << std::endl;
				ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
					delta_energy, previous_energy_for_node);

				if ( (prevrotamer_state == 0)||pass_metropolis(previous_energy_for_node,delta_energy) ) {
					//std::cerr << "pass_metropolis sc only: " << delta_energy << ", " << alternate_energy_total << ", " << bestenergy_ << std::endl;
					++num_fixbb_move_accepts;
					currentenergy = ig_->commit_considered_substitution();
					state_on_node( moltenres_id ) = rotamer_state_on_moltenres;
					if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
						bestenergy() = currentenergy;
						best_state_on_node = state_on_node;
					}

					if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
						++nsteps;
						for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
							Size iistate = state_on_node(ii);
							if ( iistate != 0 ) {
								++nsteps_for_rot( iistate + moltenres_rotoffsets[ii] );
							}
						}
					}
				}
			}

			//switch to moving the backbone only
			ig_->get_backbone_list( accessible_state_list_bbfrag );
			for ( Size n = 1; n <= cycle2; ++n ) {
				//pick a random backbone fragment from the list
				//   ranbbfrag = accessible_state_list_bbfrag( static_cast< int > (ran3() * num_accessible_bb) + 1 );

				ranbbfrag = pick_a_rotamer(j,n,cycle2,accessible_state_list_bbfrag);

				//figure out the old energy and newenergy, at this circustance, the
				//rotamers are unchanged or at least very close to the old rotamer
				//  frag_index = bbfrag_index(ranbbfrag);
				//   std::cout<<"ranbbfrag="<<ranbbfrag<<'\n';
				if ( ig_->get_backbone_currently_assigned( ranbbfrag ) ) continue;
				int num_changing_nodes = 0;

				//std::cout << "State on node: ";
				//for ( Size ii = 1; ii <= state_on_node.size(); ++ii ) { std::cout << state_on_node( ii ) << " ";}
				//std::cout << std::endl;


				//std::cout << "Consider backbone substitution " << ranbbfrag << " ";
				ig_->consider_backbone_move( ranbbfrag, delta_energy,
					previous_energy_for_bbfrag, valid_bb_move, num_changing_nodes);
				//std::cout << valid_bb_move << std::endl;
				if ( ! valid_bb_move ) continue;

				++num_bb_move_valid;
				if ( pass_metropolis_multiple_nodes_changing(
						previous_energy_for_bbfrag,delta_energy, num_changing_nodes) ) {
					++num_bb_move_accepts;
					sconly_move_accessible_state_list_current = false;
					//std::cerr << "pass_metropolis: bb(2) " << delta_energy << ", " << alternate_energy_newbbfrag << ", " << bestenergy_ << std::endl;

					//switch from prevbbfrag to ranbbfrag, make sure update all the coupled residues
					currentenergy = ig_->commit_considered_backbone_move(state_on_node);
					//since multiple sequence changed the backbone conformation
					//bbfrag_for_fragindex( frag_index ) = ranbbfrag;
					if ( currentenergy < bestenergy() ) {
						bestenergy() = currentenergy;
						best_state_on_node = state_on_node;
					}
				}

				if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
					++nsteps;
					for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
						Size iistate = state_on_node(ii);
						if ( iistate != 0 ) {
							++nsteps_for_rot( iistate + moltenres_rotoffsets[ii] );
						}
					}
				}
			}

			//at some point, switch to move side/backbone simutaneously
			//get this list once at the beginning of simulation and reuse it.
			//ig_->get_accessible_states( FlexbbInteractionGraph::BOTH_SC_AND_BB, list, num_accessible);
			//all states are accessible;
			for ( Size n = 1; n <= cycle3; ++n ) {
				//return a list of all rotamers,not sure this is fair enough, because
				//those residues with alternate backbone conformations defitenly have higher
				//probability to get picked, but at least in the quench step, the annealer needs to go through
				//the entire rotamer list exhaustively

				ranrotamer = pick_a_rotamer(j,n,cycle3,list_of_all);

				moltenres_id               = rotsets_->moltenres_for_rotamer( ranrotamer );
				//seqpos                     = rotsets_->moltenres_2_resid( moltenres_id );
				rotamer_state_on_moltenres = rotsets_->local_rotid_for_rotamer_on_moltenres( moltenres_id, ranrotamer );
				prevrotamer_state = state_on_node( moltenres_id );

				if ( prevrotamer_state == rotamer_state_on_moltenres ) continue;

				int num_changing_nodes = 0;

				//std::cout << "State on node: ";
				//for ( Size ii = 1; ii <= state_on_node.size(); ++ii ) { std::cout << state_on_node( ii ) << " ";}
				//std::cout << std::endl;


				//std::cout << "Consider backbone substitution with rotamer substitution: " << moltenres_id << " " << rotamer_state_on_moltenres << " ";
				ig_->consider_bbmove_w_state_substitution(
					moltenres_id, rotamer_state_on_moltenres,
					delta_energy, previous_energy_for_bbfrag,
					valid_bb_move, num_changing_nodes );
				//std::cout << valid_bb_move << std::endl;
				if ( ! valid_bb_move ) continue;

				++num_bb_sub_valid;
				if ( (prevrotamer_state == 0) ||
						pass_metropolis_multiple_nodes_changing(
						previous_energy_for_bbfrag, delta_energy, num_changing_nodes) ) {
					++num_bb_sub_accepts;

					/// If we accept a change in backbone conformation, then mark the state list for SCONLY moves as out-of-date.
					if ( ! rotsets_->rotamers_on_same_bbconf( moltenres_id, rotamer_state_on_moltenres, prevrotamer_state ) ) {
						sconly_move_accessible_state_list_current = false;
					}

					//std::cerr << "pass_metropolis: bb & sc " << delta_energy << ", " << alternate_energy_total << ", " << bestenergy_ << std::endl;
					currentenergy = ig_->commit_considered_backbone_move(state_on_node);
					if ( (prevrotamer_state == 0)||(currentenergy < bestenergy()) ) {
						bestenergy() = currentenergy;
						best_state_on_node = state_on_node;
					}
				}              // end Metropolis criteria

				if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
					++nsteps;
					for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
						Size iistate = state_on_node(ii);
						if ( iistate != 0 ) {
							++nsteps_for_rot( iistate + moltenres_rotoffsets[ii] );
						}
					}
				}
			}//end of both_sc_and_bb
			loopenergy(nn) = currentenergy;
		}

		/*std::cerr << "Finished Inner Iterations, temperature = " << get_temperature() << std::endl;
		std::cerr << "fixbb moves attempted: " << num_fixbb_move_valid << " accepted: " << num_fixbb_move_accepts;
		if ( num_fixbb_move_valid > 0 ) std::cerr << " (" << (double) num_fixbb_move_accepts / num_fixbb_move_valid << ")";
		std::cerr << std::endl;
		std::cerr << "bb moves attempted: " << num_bb_move_valid << " accepted: " << num_bb_move_accepts;
		if ( num_bb_move_valid > 0 ) std::cerr << " (" << (double) num_bb_move_accepts / num_bb_move_valid << ")";
		std::cerr << std::endl;
		std::cerr << "bb substitutions attempted: " << num_bb_sub_valid << " accepted: " << num_bb_sub_accepts;
		if ( num_bb_sub_accepts > 0 ) std::cerr << " (" << (double) num_bb_sub_accepts / num_bb_sub_valid << ")";
		std::cerr << std::endl;
		std::cerr << "bestenergy: " << bestenergy() << " currentenergy: " << currentenergy;
		std::cerr << std::endl << "-------------" << std::endl;*/
	}//end of outeriterations

	//std::cerr << "bestenergy after quench: " << bestenergy() << std::endl;
	//apl now convert best_state_on_node into bestrotamer_at_seqpos

	for ( Size ii = 1; ii <= nmoltenres; ++ii ) {
		Size iiresid = rotsets_->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = best_state_on_node( ii ) + rotsets_->nrotamer_offset_for_moltenres( ii );
	}

	//std::cerr << "bestenergy_ " << bestenergy() << std::endl;
	//ig_->print_current_state_assignment();

	return;
}

core::Size
FlexbbSimAnnealer::pick_a_rotamer(
	Size outercycle,
	Size innercycle,
	Size inner_loop_iteration_limit, // ?
	utility::vector1< Size > & accessible_state_list
) const
{

	Size num = 0;

	if ( quench() ) {
		num = numeric::mod<Size>( (outercycle - 1) * inner_loop_iteration_limit + innercycle - 1, accessible_state_list.size() ) + 1;
		if ( num == 1 ) {
			numeric::random::random_permutation( accessible_state_list, numeric::random::rg() );
		}
	} else {
		num = static_cast< Size > ( numeric::random::rg().uniform() * accessible_state_list.size() ) + 1;
	}
	return accessible_state_list[ num ];
}

bool
FlexbbSimAnnealer::pass_metropolis_multiple_nodes_changing(
	PackerEnergy previous_fragmentE,
	PackerEnergy deltaE,
	Size num_changing_nodes
) const
{
	assert( num_changing_nodes > 0 );
	if ( get_temperature() > 0.5 ) {
		return pass_metropolis( previous_fragmentE / num_changing_nodes , deltaE / ( num_changing_nodes > 4 ? 4 : num_changing_nodes)  );
	} else {
		return pass_metropolis( previous_fragmentE, deltaE );
	}
}


}
}
}
