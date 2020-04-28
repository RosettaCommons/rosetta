// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/annealer/SmartFixbbSimAnnealer.cc
/// @brief  Derivation of the FixbbSimAnnealer with some fancy Tensorflow logic to adaptively decrease sample space during the run
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit Headers
#include <core/pack/annealer/SmartFixbbSimAnnealer.hh>

// Package Headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack_basic/RotamerSetsBase.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/interaction_graph/AnnealableGraphBase.hh>
#include <basic/Tracer.hh>

#include <numeric/random/random.hh>
#include <utility>
#include <utility/exit.hh>

#include <iostream>
#include <fstream>
#include <math.h>

#include <basic/citation_manager/CitationManager.hh>
#include <basic/citation_manager/UnpublishedModuleInfo.hh>
#include <basic/citation_manager/CitationCollection.fwd.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/corrections.OptionKeys.gen.hh>
#include <basic/database/open.hh>

#include <core/pack/rotamer_set/FixbbRotamerSets.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#ifdef USE_TENSORFLOW
#include <tensorflow/c/c_api.h>
#include <basic/tensorflow_manager/RosettaTensorflowManager.hh>
#include <basic/tensorflow_manager/RosettaTensorflowSessionContainer.tmpl.hh>
#include <basic/tensorflow_manager/RosettaTensorflowTensorContainer.tmpl.hh>
#endif

using namespace ObjexxFCL;

static basic::Tracer TR( "core.pack.annealer.SmartFixbbSimAnnealer" );

namespace core {
namespace pack {
namespace annealer {

#ifdef USE_TENSORFLOW
using namespace basic::tensorflow_manager;

std::string
make_key( std::string const & model_name ){
	return "smart_annealer_" + model_name;
}

RosettaTensorflowSessionContainerCOP
setup_TF_session(
	std::string const & model_name
){
	std::string const key = make_key( model_name );

	RosettaTensorflowManager * manager = RosettaTensorflowManager::get_instance();
	if( manager->session_exists( key ) ){
		return manager->get_session( key );
	}

	//read from database
	std::string path_to_saved_model = basic::database::find_database_path( "protocol_data/tensorflow_graphs/smart_annealer/" + model_name + "/", "saved_model.pb" );
	//remove "saved_model.pb"
	// "saved_model.pb".size() == 14
	path_to_saved_model = path_to_saved_model.substr( 0, path_to_saved_model.size() - 14 );
	TR << "Loading network from: " << path_to_saved_model << std::endl;

	return manager->get_session( path_to_saved_model, "serve", key );
}

float
TF_run_single(
	std::array< float, 35 > const & input_data,
	RosettaTensorflowSessionContainer const & session,
	std::string const in_operation = "in",
	std::string const output_operation = "output/Sigmoid"
){
	constexpr core::Size in_size = 35;
	RosettaTensorflowTensorContainer< float > input_tensor( {1,in_size} );
	input_tensor.copy_data_from_container( input_data );

	RosettaTensorflowTensorContainer< float > output_tensor( {1,1} );

	std::chrono::duration< double, std::micro > runtime;

	session.run_session( in_operation, output_operation, input_tensor, output_tensor, runtime );

	return output_tensor( 1, 1 );
}

constexpr core::Size window_size = 5;
constexpr core::Size array_len = (window_size * 3) + 20;

#endif //USE_TENSORFLOW

void
register_unpublished_author_info() {
	using namespace basic::citation_manager;

	UnpublishedModuleInfoOP info( utility::pointer::make_shared< UnpublishedModuleInfo >(
		"SmartFixbbSimAnnealer",
		CitedModuleType::Mover,
		"Jack Magure",
		"University of North Carloina",
		"jackmaguire1444@gmail.com",
		"Created the SmartFixbbSimAnnealer."
		) );

	utility::vector1< UnpublishedModuleInfoCOP > infovect{ info };
	CitationManager::get_instance()->add_unpublished_modules( infovect );
}


////////////////////////////////////////////////////////////////////////////////
SmartFixbbSimAnnealer::SmartFixbbSimAnnealer(
	utility::vector0< int > & rot_to_pack,
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_sets,
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
	rotamer_sets,
	current_rot_index,
	calc_rot_freq,
	rot_freq
	),
	ig_(std::move(ig)),
	record_annealer_trajectory_( false )
{
	register_unpublished_author_info();
}

SmartFixbbSimAnnealer::SmartFixbbSimAnnealer(
	FArray1D_int & bestrotamer_at_seqpos,
	float & bestenergy,
	bool start_with_current, // start simulation with current rotamers
	interaction_graph::AnnealableGraphBaseOP ig,
	FixbbRotamerSetsCOP rotamer_set,
	FArray1_int & current_rot_index,
	bool calc_rot_freq,
	FArray1D_float & rot_freq
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
	),
	ig_(ig),
	record_annealer_trajectory_( false )
{
	register_unpublished_author_info();
}

/// @brief virtual destructor
SmartFixbbSimAnnealer::~SmartFixbbSimAnnealer() = default;

struct Key {
	core::Size mres;
	char AA;

	//Key() = default;
	//Key( Key const & ) = default;
	//Key( Key && ) = default;

	bool operator == ( Key const & o ) const {
		return o.mres == mres && o.AA == AA;
	}
};

//https://stackoverflow.com/questions/17016175/c-unordered-map-using-a-custom-class-type-as-the-key
struct KeyHasher {
	std::size_t operator() ( Key const & k ) const {
		return ((std::hash< char >()( k.AA )
			^ (std::hash< Size >()( k.mres ) << 1)) >> 1);
	}
};

struct Value {
	core::Size num_accept_to_AA = 0;
	core::Size num_total_to_AA = 0;
	core::Real frac_to_AA() const {
		if ( num_total_to_AA == 0 ) return -1;
		return core::Real( num_accept_to_AA ) / core::Real ( num_total_to_AA );
	}

	core::Size num_accept_from_AA = 0;
	core::Size num_total_from_AA = 0;
	core::Real frac_from_AA() const {
		if ( num_total_from_AA == 0 ) return -1;
		return core::Real( num_accept_from_AA ) / core::Real ( num_total_from_AA );
	}

	core::Size num_accept_intra_AA = 0;
	core::Size num_total_intra_AA = 0;
	core::Real frac_intra_AA() const {
		if ( num_total_intra_AA == 0 ) return -1;
		return core::Real( num_accept_intra_AA ) / core::Real ( num_total_intra_AA );
	}

	std::string to_string() const {
		return std::to_string( frac_to_AA() ) + "_" +
			std::to_string( frac_from_AA() ) + "_" +
			std::to_string( frac_intra_AA() );
	}

};

namespace {
#ifdef USE_TENSORFLOW

int
name1_to_index( char const name1 ){//using 0-indexing
	switch( name1 ){
	case('A'): return 0;
	case('C'): return 1;
	case('D'): return 2;
	case('E'): return 3;
	case('F'): return 4;
	case('G'): return 5;
	case('H'): return 6;
	case('I'): return 7;
	case('K'): return 8;
	case('L'): return 9;
	case('M'): return 10;
	case('N'): return 11;
	case('P'): return 12;
	case('Q'): return 13;
	case('R'): return 14;
	case('S'): return 15;
	case('T'): return 16;
	case('V'): return 17;
	case('W'): return 18;
	case('Y'): return 19;
	default: return -1;
	}
}

std::unordered_map< Key, bool, KeyHasher >
determine_likely_mutations(
	std::unordered_map< Key, utility::vector1< Value >, KeyHasher > const & my_data_map,
	core::Real const cutoff,
	std::string const & model_name
){
	std::unordered_map< Key, bool, KeyHasher > likely; //this can just be a std::set
	likely.max_load_factor( 0.1 );

	RosettaTensorflowSessionContainerCOP const session =
		setup_TF_session( model_name );

	for( auto const & pair : my_data_map ){
		core::Real my_cutoff = cutoff;
		Key const key = pair.first;
		char const AA = key.AA;

		switch( AA ){//ex1ex2 rots:
			case( 'A' )://1 rot
			case( 'G' )://1 rot
			case( 'P' )://5 rots
			//case( 'V' )://2 rots //Need to double-check this before including valine. 2 seems too small
				//cheap to keep, not many rotamers, let's just keep them
				likely[ key ] = true;
				continue;
			break;//redundant break, but just being safe
			case( 'W' )://59 rots
			case( 'F' )://32 rots
				//high false negative rate
				//too many rotamers to just blindly keep
				//let's just decrease the cutoff a hair
				my_cutoff *= 0.75;
			break;
		default:
			break;
		}//switch

		utility::vector1< Value > const & values = pair.second;
		std::array< float, array_len > network_input;//Zero Indexed!
		network_input.fill( 0.0 );

		for( core::Size ii = 0; ii < window_size; ++ii ){
			network_input[ ii * 3 ]       = values[ ii + 1 ].frac_to_AA();
			network_input[ (ii * 3) + 1 ] = values[ ii + 1 ].frac_from_AA();
			network_input[ (ii * 3) + 2 ] = values[ ii + 1 ].frac_intra_AA();
		}

		int const n1_index = name1_to_index( AA );
		debug_assert( n1_index != -1 );
		network_input[ 15 + n1_index ] = 1.0;

		float const prediction = TF_run_single( network_input, *session );
		likely[ key ] = prediction > my_cutoff;
	}

	return likely;
}

#endif//USE_TENSORFLOW
}

/// @brief sim_annealing for fixed backbone design mode
void
SmartFixbbSimAnnealer::run(){

	if ( ! has_been_initialized_from_task_ ) {
		utility_exit_with_message( "initialize_from_task() needs to be called before running this annealer. This is a developer error, you likely did nothing wrong if you are seeing this message in a production run." );
	}

#ifndef USE_TENSORFLOW
	utility_exit_with_message( "Thanks for using the SmartFixbbSimAnnealer, I really hope you like it. Unfortunately, it requires that you compile with extras=tensorflow. Once you do that, please make sure you run a version of rosetta that has 'tensorflow' in the name (rosetta_scripts.mpitensroflow.linuxgccrelease, for example)." );
#else

	TR << "Running with settings: " << std::endl;
	TR << "model: " << model_ << std::endl;
	TR << "cutoff: " << cutoff_ << std::endl;
	TR << "pick_again: " << pick_again_ << std::endl;

	perform_validation_test();
	//get_correct_io_layers_for_model();

	int const nmoltenres = ig_->get_num_nodes();

	FArray1D_int state_on_node( nmoltenres,0 ); // parallel representation of interaction graph's state
	FArray1D_int best_state_on_node( nmoltenres,0 );
	FArray1D_float loopenergy(maxouteriterations,0.0);

	//bk variables for calculating rotamer frequencies during simulation
	int nsteps = 0;
	FArray1D_int nsteps_for_rot( ig_->get_num_total_states(), 0 );

	//--------------------------------------------------------------------
	//initialize variables

	core::PackerEnergy currentenergy = 0.0;

	ig_->prepare_graph_for_simulated_annealing();
	ig_->blanket_assign_state_0();

	//--------------------------------------------------------------------
	if ( num_rots_to_pack() == 0 ) return;

	setup_iterations();

	FArray1D_float previous_nsteps_for_rot( rotamer_sets()->nrotamers(), 0.0);

	int outeriterations = get_outeriterations();

	/// if pose has water molecules then check if virtualization is allowed
	/// (default=true) so that water molecules may be virtualized 50% of the time
	bool include_vrt = basic::options::option[ basic::options::OptionKeys::corrections::water::include_vrt ].value();

	std::ofstream annealer_trajectory;
	if ( record_annealer_trajectory_ ) {
		annealer_trajectory.open(trajectory_file_name_.c_str() );
	}

	std::unordered_map< Key, utility::vector1< Value >, KeyHasher > my_data_map;
	my_data_map.max_load_factor( 0.1 );

	Key previous_key;
	previous_key.mres = 9999;
	previous_key.AA = '=';//using starting values that are unlikely to match
	for ( core::Size rotid = 1; rotid <= rotamer_sets()->nrotamers(); ++rotid ) {
		Key key;
		key.mres = rotamer_sets()->moltenres_for_rotamer( rotid );
		key.AA = rotamer_sets()->rotamer( rotid )->name1();

		if ( key == previous_key ) continue;
		else previous_key = key;

		utility::vector1< Value > & vec = my_data_map[ key ];
		if ( vec.empty() ) {
			vec.resize( window_size );
		}
	}

	std::unordered_map< Key, bool, KeyHasher > likely_mutations;//this will become populated after the window_size-th outer iteration

	//outer loop
	for ( int nn = 1; nn <= outeriterations; ++nn ) {
		setup_temperature(loopenergy,nn);

		if ( quench() ) {
			currentenergy = bestenergy();
			state_on_node = best_state_on_node;
			ig_->set_network_state( state_on_node );
		}

		int inneriterations = get_inneriterations();

		float treshold_for_deltaE_inaccuracy = std::sqrt( get_temperature() );
		ig_->set_errorfull_deltaE_threshold( treshold_for_deltaE_inaccuracy );

		//inner loop
		for ( int n = 1; n <= inneriterations; ++n ) {
			int ranrotamer = pick_a_rotamer( n );
			if ( ranrotamer == -1 ) continue;

			int moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );

			if ( core::Size( nn ) > window_size && ! quench() ) {
				Key key;
				key.mres = moltenres_id;
				key.AA = rotamer_sets()->rotamer( ranrotamer )->name1();
				debug_assert( likely_mutations.find( key ) != likely_mutations.end() );

				bool const bypass_because_of_quenching = disable_during_quench_ && quench();
				if( ! bypass_because_of_quenching ){
					if ( ! pick_again_ ) {
						if ( ! likely_mutations[ key ] ) {
							continue;// for int n
						}
					} else { //keep sampling
						core::Size attempts = 0;
						constexpr core::Size max_attempts = 10;

						while ( ! likely_mutations[ key ] ) {
							++attempts;
							if( attempts > max_attempts ) break;//just to prevent an infinite loop

							ranrotamer = pick_a_rotamer( n );
							moltenres_id = rotamer_sets()->moltenres_for_rotamer( ranrotamer );
							key.mres = moltenres_id;
							key.AA = rotamer_sets()->rotamer( ranrotamer )->name1();
						}
						if( attempts == max_attempts ) continue;
					}
				}
			}


			/// removed const from rotamer_state_on_moltenres for code that virtualizes waters half the time
			int rotamer_state_on_moltenres = rotamer_sets()->rotid_on_moltenresidue( ranrotamer );
			int const prevrotamer_state = state_on_node(moltenres_id);

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// for waters, set to virtual 50% of the time
			if ( ( include_vrt ) && ( rotamer_sets()->rotamer_for_moltenres( moltenres_id, rotamer_state_on_moltenres )->type().is_virtualizable_by_packer() ) ) {  // don't want to do this for waters without a virtual state, like TP3
				core::Size const nrot = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
				/// last state is the virtual state
				if ( numeric::random::rg().uniform() < (static_cast<core::Real>(nrot)/2.0-1)/static_cast<core::Real>(nrot) ) rotamer_state_on_moltenres = rotamer_sets()->nrotamers_for_moltenres(moltenres_id);
			}

			if ( rotamer_state_on_moltenres == prevrotamer_state ) continue; //skip iteration

			// initializing to zero but should be updated below.
			core::PackerEnergy previous_energy_for_node( 0.0 ), delta_energy( 0.0 );

			ig_->consider_substitution( moltenres_id, rotamer_state_on_moltenres,
				delta_energy, previous_energy_for_node);

			Key k_to;
			Key k_from;
			if ( prevrotamer_state != 0 && core::Size( nn ) <= window_size ) {
				k_to.mres = moltenres_id;
				k_to.AA = rotamer_sets()->rotamer_for_moltenres( moltenres_id, rotamer_state_on_moltenres )->name1();
				utility::vector1< Value > & inner_map_to = my_data_map[ k_to ];
				Value & val_to = inner_map_to[ nn ];

				k_from.mres = moltenres_id;
				k_from.AA = rotamer_sets()->rotamer_for_moltenres( moltenres_id, prevrotamer_state )->name1();
				utility::vector1< Value > & inner_map_from = my_data_map[ k_from ];
				Value & val_from = inner_map_from[ nn ];

				if ( k_from.AA == k_to.AA ) {
					++val_to.num_total_intra_AA;
				} else {
					++val_to.num_total_to_AA;
					++val_from.num_total_from_AA;
				}
			}


			//bk keep new rotamer if it is lower in energy or accept it at some
			//bk probability if it is higher in energy, if it is the first
			//bk rotamer to be tried at this position automatically accept it.
			if ( (prevrotamer_state == 0) || pass_metropolis(previous_energy_for_node,delta_energy) ) {

				if ( prevrotamer_state != 0 && core::Size( nn ) <= window_size ) {
					utility::vector1< Value > & inner_map_to = my_data_map[ k_to ];
					Value & val_to = inner_map_to[ nn ];

					utility::vector1< Value > & inner_map_from = my_data_map[ k_from ];
					Value & val_from = inner_map_from[ nn ];

					if ( k_from.AA == k_to.AA ) {
						++val_to.num_accept_intra_AA;
					} else {
						++val_to.num_accept_to_AA;
						++val_from.num_accept_from_AA;
					}

				}


				//std::cout << " accepted\n";
				currentenergy = ig_->commit_considered_substitution();
				state_on_node(moltenres_id) = rotamer_state_on_moltenres;
				if ( (prevrotamer_state == 0)||(currentenergy < bestenergy() ) ) {
					best_state_on_node = state_on_node;
					bestenergy() = currentenergy;
				}

				if ( record_annealer_trajectory_ ) {
					annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " A\n";
				}

			} else if ( record_annealer_trajectory_ ) {
				annealer_trajectory << moltenres_id << " " << rotamer_state_on_moltenres << " R\n";
			}

			loopenergy(nn) = currentenergy;
			float const temperature = get_temperature();

			if ( calc_rot_freq() && ( temperature <= calc_freq_temp ) ) {
				++nsteps;
				for ( int ii = 1; ii <= nmoltenres; ++ii ) {
					int iistate = state_on_node(ii);
					if ( iistate != 0 ) {
						++nsteps_for_rot( rotamer_sets()->moltenres_rotid_2_rotid(ii, iistate) );
					}
				}
			}

		} // end of inneriteration loop

		if ( core::Size( nn ) == window_size ) {
			likely_mutations = determine_likely_mutations( my_data_map, cutoff_, model_ );
#ifndef NDEBUG
			core::Size nrots_kept = 0;
			for ( core::Size rotid = 1; rotid <= rotamer_sets()->nrotamers(); ++rotid ) {
				Key key;
				key.mres = rotamer_sets()->moltenres_for_rotamer( rotid );
				key.AA = rotamer_sets()->rotamer( rotid )->name1();
				if( likely_mutations[ key ] ){
					++nrots_kept;
				}
			}
			TR << "Kept " << nrots_kept << " of " << rotamer_sets()->nrotamers() << " rotamers" << std::endl;
#endif
		}

	} //end of outeriteration loop: nn

#ifndef NDEBUG
	TR << "pack_rotamers run final: ";
	for ( Size ii=0; ii < best_state_on_node.size(); ii++ ) {
		if ( best_state_on_node[ii] != 0 ) {
			TR << rotamer_sets()->rotamer_set_for_moltenresidue( ii+1 )->rotamer( best_state_on_node[ii] )->name1();
		} else {
			TR << '-';
		}
	}
	TR << ", best_energy: " << bestenergy() << std::endl;
#endif

	if ( ig_->any_vertex_state_unassigned() ) {
		std::cerr << "Critical error -- In SmartFixbbSimAnnealer, one or more vertex states unassigned at annealing's completion." << std::endl;
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
	for ( int ii = 1; ii <= nmoltenres; ++ii ) {
		int const iiresid = rotamer_sets()->moltenres_2_resid( ii );
		bestrotamer_at_seqpos()( iiresid ) = rotamer_sets()->moltenres_rotid_2_rotid( ii, best_state_on_node(ii));
	}

#endif
}

void SmartFixbbSimAnnealer::record_annealer_trajectory( bool setting ) {
	record_annealer_trajectory_ = setting;
}

void SmartFixbbSimAnnealer::trajectory_file_name( std::string const & setting ) {
	trajectory_file_name_ = setting;
}

core::Real
expected_output_for_model_name( std::string const & model_name ){
	//Input: L,0.932039,0.923077,-1,0.828571,0.90625,1,0.586207,0.809524,-1,0.461538,0.5,0,0.210526,0.246914,1,0.0
	if ( model_name == "generation1" ) {
		return 0.92780346;
	} else if ( model_name == "generation2" ) {
		return 0.961535;
	} else {
		return -1;
	}
}

void
SmartFixbbSimAnnealer::perform_validation_test() const {
#ifdef USE_TENSORFLOW
	core::Real const expected_output = expected_output_for_model_name( model_ );
	if( expected_output < 0 ){
		TR << "expected output is not registered for model name " << model_ << ". Skipping validation test." << std::endl;
		return;
	}

	auto const tf_info = setup_TF_session( model_ );
	std::array< float, array_len > network_input;
	{
		//Input: L,0.932039,0.923077,-1,0.828571,0.90625,1,0.586207,0.809524,-1,0.461538,0.5,0,0.210526,0.246914,1,0.0
		network_input.fill( 0.0 );

		network_input[ 0 ] = 0.932039;
		network_input[ 1 ] = 0.923077;
		network_input[ 2 ] = -1;
		network_input[ 3 ] = 0.828571;
		network_input[ 4 ] = 0.90625;
		network_input[ 5 ] = 1;
		network_input[ 6 ] = 0.586207;
		network_input[ 7 ] = 0.809524;
		network_input[ 8 ] = -1;
		network_input[ 9 ] = 0.461538;
		network_input[ 10 ] = 0.5;
		network_input[ 11 ] = 0;
		network_input[ 12 ] = 0.210526;
		network_input[ 13 ] = 0.246914;
		network_input[ 14 ] = 1;

		int const n1_index = name1_to_index( 'L' );
		debug_assert( n1_index != -1 );
		network_input[ 15 + n1_index ] = 1.0;
		float const prediction = TF_run_single( network_input, *tf_info );
		TR << "Validation Test - prediction: " << prediction << ", expected: " << expected_output << std::endl;
	}
#endif
}

void
SmartFixbbSimAnnealer::initialize_from_task(
	core::pack::task::PackerTask const & task
) {
	model_ = task.smart_annealer_model();
	cutoff_ = task.smart_annealer_cutoff();
	pick_again_ = task.smart_annealer_pick_again();
	disable_during_quench_ = task.smart_annealer_disable_during_quench();
	has_been_initialized_from_task_ = true;
}


/*
I'm leaving this here because some poor soul might find it helpful someday when they can't figure out which layers to use as input and output for their model.

void
SmartFixbbSimAnnealer::get_correct_io_layers_for_model() const {
#ifdef USE_TENSORFLOW
auto const tf_info( setup_TF_session( model_ ) );
utility::vector1< std::string > operations;
{
size_t pos = 0;
TF_Operation* oper;
TR << "All Operations: " << std::endl;
while ((oper = TF_GraphNextOperation( tf_info->graph(), &pos )) != nullptr) {
operations.emplace_back( TF_OperationName( oper ) );
TR << TF_OperationName( oper ) << std::endl;
}
}

std::array< float, array_len > network_input;
{//debug.h5
//Input: L,0.932039,0.923077,-1,0.828571,0.90625,1,0.586207,0.809524,-1,0.461538,0.5,0,0.210526,0.246914,1,0.0
//Output: 0.97097194
network_input.fill( 0.0 );

network_input[ 0 ] = 0.932039;
network_input[ 1 ] = 0.923077;
network_input[ 2 ] = -1;
network_input[ 3 ] = 0.828571;
network_input[ 4 ] = 0.90625;
network_input[ 5 ] = 1;
network_input[ 6 ] = 0.586207;
network_input[ 7 ] = 0.809524;
network_input[ 8 ] = -1;
network_input[ 9 ] = 0.461538;
network_input[ 10 ] = 0.5;
network_input[ 11 ] = 0;
network_input[ 12 ] = 0.210526;
network_input[ 13 ] = 0.246914;
network_input[ 14 ] = 1;

network_input[ 15 + name1_to_index( 'L' ) ] = 1.0;
//float const prediction = TF_run_single( network_input, *tf_info );
//TR << "prediction: " << prediction << ", expected: 0.97097194" << std::endl;
}

for( core::Size ii = 1; ii < operations.size(); ++ii ){
for( core::Size jj = operations.size(); jj > ii; --jj ){
try{
float const p = TF_run_single( network_input, * tf_info, operations[ ii ], operations[ jj ] );
if( p > 0.95 && p < 1.0 ){
TR << operations[ ii ] << " " << operations[ jj ] << " " << p << " 0.97097194" << std::endl;
}
} catch( ... ){}
}//jj
}//ii

#endif
}
*/

}//end namespace annealer
}//end namespace pack
}//end namespace core
