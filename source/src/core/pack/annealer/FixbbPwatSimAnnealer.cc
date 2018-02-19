// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/annealer/FixbbPwatSimAnnealer.cc
/// @brief  fixed temperature annealer for packing with statistical PWAT water model
/// @author modified from FixbbSimAnnealer.cc by Ryan Pavlovicz (rpavlov@uw.edu)
/// @author FixbbSimAnnealer.cc originally by Andrew Leaver-Fay (aleaverfay@gmail.com)

// Unit Headers
#include <core/pack/annealer/FixbbPwatSimAnnealer.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetsBase.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>

#include <iostream>
#include <fstream>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>

#include <basic/options/option.hh>

// option key includes
#include <basic/options/keys/corrections.OptionKeys.gen.hh>

using namespace ObjexxFCL;

//#ifndef NDEBUG
static basic::Tracer TR("core.pack.annealer.FixbbPwatSimAnnealer");
//#endif

namespace core {
namespace pack {
namespace annealer {

using namespace basic::options;

////////////////////////////////////////////////////////////////////////////////
/// @begin FixbbPwatSimAnnealer::FixbbPwatSimAnnealer()
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
FixbbPwatSimAnnealer::FixbbPwatSimAnnealer(
	utility::vector0< int > & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq,
	utility::vector1< PointDwell > & all_rot
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
	rot_freq//,
	// all_rot
	),
	ig_(ig),
	all_rot_(all_rot)
{
}

/// @brief virtual destructor
FixbbPwatSimAnnealer::~FixbbPwatSimAnnealer()
{}

bool sortDwell(const PointDwell &i, const PointDwell &j) { return i.dwell > j.dwell; }

/// @brief sim_annealing for fixed backbone design mode
void FixbbPwatSimAnnealer::run()
{
	bool record_dwell = option[ OptionKeys::corrections::water::record_dwell ].value();
	bool print_dwell = option[ OptionKeys::corrections::water::print_dwell ].value();
	bool spike_anneal = option[ OptionKeys::corrections::water::spike_anneal ].value();
	Real spike_steps_scale  = option[ OptionKeys::corrections::water::spike_steps_scale ].value();
	int spike_cycles  = option[ OptionKeys::corrections::water::spike_cycles ].value();
	Real dwell_cutoff = option[ OptionKeys::corrections::water::dwell_cutoff ].value();

	Real record_count = 0;

	int const nmoltenres = ig_->get_num_nodes();
	int totalrot = 0;
	for ( int resi = 1; resi <= nmoltenres; resi++ ) {
		totalrot += rotamer_sets()->nrotamers_for_moltenres(resi);
	}

	// problem if annealer needs to be called before water packing
	// for example in case of missing side chains in input pdb
	// iterate through moltenres and make sure PWAT is found
	// if not, turn record_dwell off to allow normal packing
	bool packwat = false;
	for ( int resi = 1; resi <= nmoltenres; resi++ ) {
		std::string resname = rotamer_sets()->rotamer_set_for_moltenresidue(resi)->rotamer(1)->name();
		if ( resname == "PWAT" ) {
			packwat = true;
			break;
		}
	}
	if ( packwat == false ) {
		TR << "turning record_dwell off!" << std::endl;
		record_dwell = false;
	}

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//--------------------------------------------------------------------
	//initialize variables

	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	setup_iterations();

	int outeriterations = get_outeriterations();
	if ( record_dwell ) {
		if ( spike_anneal ) {
			outeriterations = spike_cycles; // spiking cycles only works with record_dwell
		} else {
			outeriterations = 1;  // run single temp step if recording annealer trajectory
		}
	}

	std::ofstream dwell_file;
	uint nmolres = rotamer_sets()->nmoltenres();
	utility::vector1< utility::vector1<int> > rot_states(nmolres);
	if ( record_dwell ) {
		if ( print_dwell ) {
			std::string const dwell_file_name = option[ OptionKeys::corrections::water::dwell_name ].value();
			dwell_file.open(dwell_file_name.c_str() );
		}
		for ( uint roti=1; roti<=nmolres; roti++ ) {
			rot_states[roti].resize(rotamer_sets()->nrotamers_for_moltenres(roti),0);
		}
	}

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		if ( record_dwell ) {
			set_temperature( option[ OptionKeys::corrections::water::pack_temp ].value() );
		} else {
			setup_temperature(loopenergy,nn);
		}

		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}

		int inneriterations = get_inneriterations();

		float treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//high temperature scrambling of states
		if ( spike_anneal ) {    // start with high temp if spike_anneal
			int spikeaccepts = 0;
			set_temperature( 100.0 );
			int spike_steps = round(totalrot*spike_steps_scale);
			for ( int n = 1; n <= spike_steps; ++n ) {
				int const ranrotamer = pick_a_rotamer( n );
				if ( ranrotamer == -1 ) continue;

				int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
				int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
				int const prevrotamer_state = state_on_node(moltenres_id);

				if ( rotamer_state_on_moltenres == prevrotamer_state ) continue;

				core::conformation::Residue curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
				if ( curres.aa() == core::chemical::aa_h2o ) {
					double rand_num = numeric::random::rg().uniform();
					double nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
					if ( rand_num < (nrot/2.0-1)/nrot ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				}
				core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );
				ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
					delta_energy, previous_energy_for_node);
				if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
					spikeaccepts++;
					currentenergy = ig_->commit_considered_substitution();
					state_on_node(moltenres_id) = rotamer_state_on_moltenres;
					if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
						best_state_on_node = state_on_node;
						bestenergy() = currentenergy;
					}
				}
			}
		}
		if ( spike_anneal ) {
			// set temperature back to record temp
			set_temperature( option[ OptionKeys::corrections::water::pack_temp ].value() );
			inneriterations = 6*totalrot;   // increased from 5 to 6 to allow for quiet period of totalrot
		}
		int accepts = 0;
		for ( int n = 1; n <= inneriterations; ++n ) {  // normal inneriterations
			int const ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int const moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// for waters, set to virtual 50% of the time
			core::conformation::Residue curres( *rotamer_sets()->rotamer_for_moltenres(moltenres_id, rotamer_state_on_moltenres) );
			if ( curres.aa() == core::chemical::aa_h2o ) {
				double rand_num = numeric::random::uniform();
				// last rotamer is the virtual state
				double nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				if ( rand_num < (nrot/2.0-1)/nrot ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			// initializing to zero but should be updated below.
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);

			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.

			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {
				accepts++;
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;
					bestenergy() = currentenergy;
				}

			}

			// iterate through nmoltenres and count the current state to accumulate dwell times
			if ( n > totalrot ) {  // allow for totalrot steps of equilibration before recording
				record_count++;
				if ( record_dwell ) {
					for ( uint roti=1; roti <= nmolres; ++roti ) {
						if ( state_on_node(roti) == 0 ) {
							continue;
						} else {
							rot_states[roti][state_on_node(roti)]++;
						}
					}
				}
			}

			loopenergy(nn) = currentenergy;

		} // end of inneriteration loop
	} //end of outeriteration loop

	if ( record_dwell ) {
		int watcount = 0;
		utility::vector1< PointDwell > rot_cutoff;
		for ( Size x = 1; x <= rot_states.size(); ++x ) {
			std::string resname = rotamer_sets()->rotamer_set_for_moltenresidue(x)->rotamer(1)->name();
			if ( resname == "PWAT" ) {
				watcount++;
				if ( print_dwell ) dwell_file << rot_states[x] << std::endl;
				utility::vector1<Real> normvect(rot_states[x].size()-1); // -1 to skip virtual state
				for ( Size y = 1; y <= normvect.size(); ++y ) {
					normvect[y] = rot_states[x][y]/record_count;
					PointDwell current_rotamer;
					core::conformation::ResidueCOP pwat_rot = rotamer_sets()->rotamer_set_for_moltenresidue(x)->rotamer(y);
					current_rotamer.xyz = pwat_rot->xyz("O");
					current_rotamer.dwell = normvect[y];
					all_rot().push_back(current_rotamer);
					if ( normvect[y] >= dwell_cutoff ) {
						PointDwell new_cutoff;
						core::conformation::ResidueCOP keeprot = rotamer_sets()->rotamer_set_for_moltenresidue(x)->rotamer(y);
						new_cutoff.xyz = keeprot->xyz("O");
						new_cutoff.dwell = normvect[y];
						rot_cutoff.push_back(new_cutoff);
					}
				}
			}
		}
	}


	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In FixbbPwatSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
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
		assert( ! ig_->any_vertex_state_unassigned() );
		utility_exit();
	}

	//convert best_state_on_node into best_rotamer_at_seqpos
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}

	//std::cout << "Annealing ends" << std::endl;
}

utility::vector1< PointDwell > & FixbbPwatSimAnnealer::all_rot() { return all_rot_; }

}//end namespace annealer
}//end namespace pack
}//end namespace core
