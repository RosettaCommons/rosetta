// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseProteinBackboneSampler
/// @brief Makes a list of (phi, psi, omega) at moving_residues that
///              could be useful for full-atom packing
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/StepWiseProteinBackboneSampler.hh>
#include <protocols/stepwise/modeler/protein/util.hh>
#include <protocols/stepwise/modeler/protein/MainChainTorsionSet.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/rna/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoreType.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/fixedsizearray1.hh>

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise modeler of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.protein.StepWiseProteinBackboneSampler" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseProteinBackboneSampler::StepWiseProteinBackboneSampler(  working_parameters::StepWiseWorkingParametersCOP working_parameters ):
	working_parameters_( working_parameters ),
	moving_residues_input_( working_parameters_->working_moving_res_list() ),
	n_sample_( 18 /* Corresponds to 20 degree bins */ ),
	n_sample_beta_( 6 /* Corresponds to 60 degree bins, for betas */ ),
	rmsd_cutoff_( -1.0 ),
	silent_file_( "" ),
	filter_native_big_bins_( false ),
	centroid_screen_( false ),
	centroid_score_ref_( 999999999999999.99),
	centroid_score_diff_cut_( 20.0 ),
	apply_vdw_cut_( false ),
	centroid_vdw_ref_( 9999999999999999.999 ),
	nstruct_centroid_( 0 ),
	ramachandran_( core::scoring::ScoringManager::get_instance()->get_Ramachandran() ),
	ghost_loops_( false ),
	is_pre_proline_( working_parameters_->is_pre_proline() ),
	expand_loop_takeoff_( false ) // may switch to true soon. connects all psi,omega,phi in CA-to-CA connections.
{
	initialize_is_fixed_res();
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseProteinBackboneSampler::~StepWiseProteinBackboneSampler()
{}

/////////////////////
std::string
StepWiseProteinBackboneSampler::get_name() const {
	return "StepWiseProteinBackboneSampler";
}


//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::apply( core::pose::Pose & pose )
{
	Size which_res( 1 );
	Size count( 1 );

	clock_t const time_start( clock() );

	define_moving_res( pose );
	setup_torsion_sets( pose );
	centroid_scores_.clear();

	// convert to centroid
	pose::Pose pose_save = pose;
	if ( centroid_screen_ ) convert_to_centroid( pose );
	if ( ghost_loops_ ) prepare_ghost_pose( pose ); // pose with loops excised

	sample_residues_recursively( which_res, count, pose );

	TR.Debug << "Total time in StepWiseProteinBackboneSampler: " <<
		static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

	pose = pose_save;

	if ( nstruct_centroid_ > 0 && centroid_screen_ ) filter_main_chain_torsion_sets();
}


///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::define_moving_res( pose::Pose const & pose ) {

	is_fixed_res_ = is_fixed_res_input_;
	moving_residues_ = moving_residues_input_;
	runtime_assert( is_fixed_res_.size() == pose.size() );

	if ( !expand_loop_takeoff_ ) return;

	// expand_loop_takeoff makes protein modeler similar to RNA -- sample psi, omega, *and* phi
	// at CA-to-CA connections going into and out of moving residues, just like in RNA,
	// we sample epsilon, zeta, alpha, beta, and gamma in each sugar-to-sugar connection.
	for ( Size const moving_res : moving_residues_input_ ) {
		Size const takeoff_res = pose.fold_tree().get_parent_residue( moving_res );
		if ( takeoff_res == 0 ) continue;
		if ( pose.fold_tree().jump_nr( moving_res, takeoff_res ) > 0 ) continue; // jump
		runtime_assert( (moving_res == takeoff_res + 1) || (moving_res == takeoff_res - 1 ) );
		if ( !moving_residues_.has_value( takeoff_res ) ) {
			moving_residues_.push_back( takeoff_res );
			is_fixed_res_[ takeoff_res ] = true;
		}
	}
	std::sort( moving_residues_.begin(), moving_residues_.end() );

	runtime_assert( is_fixed_res_.size() == pose.size() );
}

///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::setup_torsion_sets( core::pose::Pose const & pose ) {

	main_chain_torsion_set_for_moving_residues_.clear();
	for ( Size const moving_res : moving_residues_ ) {
		//
		//}Size n = 1; n <= moving_residues_.size(); n++ )  {
		// Preserve alpha behavior rather than trying to unify.
		if ( pose.residue_type( moving_res ).is_alpha_aa() ) {
			main_chain_torsion_set_for_moving_residues_.emplace_back( 0.0, 0.0, 0.0 );
		} else if ( pose.residue_type( moving_res ).is_beta_aa() ) {
			// Obviously this will eventually be general!
			utility::fixedsizearray1< Real, 3 > dihs;
			main_chain_torsion_set_for_moving_residues_.emplace_back( dihs, 0.0 );
		}
	}

	main_chain_torsion_sets_for_moving_residues_.clear();
}


////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::sample_residues_recursively(
	Size const which_res,
	Size & count,
	core::pose::Pose & pose )
{

	using namespace core::chemical;
	using namespace core::scoring;

	if ( moving_residues_.size() < 1 ) return; //nothing to do.

	Size const n = moving_residues_[ which_res ];

	// todo: make a vec of vecs of Real.
	MainChainTorsionSetList main_chain_torsion_set_list;
	get_main_chain_torsion_set_list( n, pose, main_chain_torsion_set_list );

	for ( auto const & main_chain_torsion_set : main_chain_torsion_set_list ) {

		if ( pose.residue_type( n ).is_alpha_aa() ) {
			// For now don't risk changing alpha behavior!
			pose.set_phi( n, main_chain_torsion_set.phi() );
			pose.set_psi( n, main_chain_torsion_set.psi() );
			pose.set_omega( n, main_chain_torsion_set.omega() );
		} else if ( pose.residue_type( n ).is_beta_aa() ) {
			// change to loop over MC torsion index for rsd n
			// Ugh these could be unequal.
			//for ( Size ii = 1; ii <= pose.residue( n ).mainchain_torsions(); ++ii ) {
			for ( Size ii = 1; ii <= main_chain_torsion_set.mainchain_dihedral_values().size(); ++ii ) {
				pose.set_torsion( core::id::TorsionID( n, id::BB, ii ), main_chain_torsion_set.mainchain_dihedral_values()[ ii ] );
			}
		}

		main_chain_torsion_set_for_moving_residues_[ which_res ] = main_chain_torsion_set;

		if ( which_res == moving_residues_.size() ) {

			count++;
			std::string const tag = "S_"+ ObjexxFCL::lead_zero_string_of( count, 5 );
			filter_and_save( pose, tag );

		} else {
			sample_residues_recursively( which_res+1, count, pose );
		}
	}
}


//////////////////////////////////////////////////////////////////////
// This could even be its own class...
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list(
	Size const n,
	pose::Pose const & pose,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	using namespace core::chemical;
	utility::vector1< Size > const fixed_domain_map = core::pose::full_model_info::get_fixed_domain_from_full_model_info_const( pose );
	main_chain_torsion_set_list.clear();

	// following is totally hacky... would be much better to
	// figure out which bins would sum to give probability > 99%, or something like that.
	Real best_energy_cutoff = 0.8;

	// A few special cases first.
	if ( n == 1
			&& !pose.residue( 1 ).has_variant_type( core::chemical::N_ACETYLATION )
			//&& !pose.residue( 1 ).has_variant_type( "LOWER_TERMINUS" )  /*new!*/
			) {

		// If we're at the N-terminal, there's no point in modeler phi -- just a dinky hydrogen...
		// or three hydrogens.
		// basically slave the phi to the psi.
		get_main_chain_torsion_set_list_n_terminus( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

	} else if ( n == pose.size()
			&& !pose.residue( n ).has_variant_type( core::chemical::C_METHYLAMIDATION )
			&& !pose.residue( n ).has_variant_type( core::chemical::UPPER_TERMINUS_VARIANT ) /*new!*/
			) {
		// DEFUNCT: If we're at the C-terminal, there's no point in modeler psi -- just an oxygen.
		// basically slave the psi to the phi.
		get_main_chain_torsion_set_list_c_terminus( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

	} else if ( is_fixed_res_[ n ] &&
			(n == pose.size() || is_fixed_res_[ n+1 ]) &&
			(n > 1 && ( !is_fixed_res_[ n-1] ||
			( fixed_domain_map[n-1] != fixed_domain_map[n] ) ) ) ) {

		get_main_chain_torsion_set_list_sample_phi_only( n, pose,  best_energy_cutoff, main_chain_torsion_set_list );
		// TR << "Sampling phi only for: " << n << ' ' << main_chain_torsion_set_list.size() << std::endl;

	} else if ( is_fixed_res_[ n ] &&
			(n == 1 || is_fixed_res_[ n-1 ]) &&
			(n < pose.size() && ( !is_fixed_res_[n+1] ||
			( fixed_domain_map[n] != fixed_domain_map[n+1] ) ) ) ) {

		get_main_chain_torsion_set_list_sample_psi_only( n, pose, best_energy_cutoff, main_chain_torsion_set_list );
		//    TR << "Sampling psi only for: " << n << ' ' << main_chain_torsion_set_list.size() << std::endl;

	} else if ( n > moving_residues_[ 1 ] &&
			n < moving_residues_[ moving_residues_.size() ] ) {

		if ( expand_loop_takeoff_ ) {
			get_main_chain_torsion_set_list_full( n, pose, best_energy_cutoff, main_chain_torsion_set_list );
		} else {
			// Trying coarse sample for internal residues -- otherwise number of conformations really blows up.
			get_main_chain_torsion_set_list_coarse( n, pose, main_chain_torsion_set_list );
		}

	} else {

		/////////////////
		// General case.
		/////////////////
		get_main_chain_torsion_set_list_full( n, pose, best_energy_cutoff, main_chain_torsion_set_list );

	}

	//TR << "Filtering..." << std::endl;
	if ( filter_native_big_bins_ )  {
		filter_native_BIG_BINS( n, main_chain_torsion_set_list );
	} else {
		filter_based_on_desired_secstruct( pose.secstruct( n ), main_chain_torsion_set_list );
	}

	if ( is_pre_proline_[ n ] && ( n == pose.size() || !is_fixed_res_[ n+1 ] ) ) {
		TR.Debug << "-----  SAMPLING CIS OMEGA --------" << std::endl;
		sample_cis_omega( main_chain_torsion_set_list );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_coarse(
	Size const n,
	pose::Pose const & pose,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	using namespace core::chemical;

	// AMW TODO this.
	// Later, can individually dial in cluster centers based on careful
	// inspection of residue-specific rama plots. For now, this is basically a quick hack.
	if ( pose.aa( n ) == aa_pro ) {
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, -30.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet(  60.0, 140.0 ) );
	} else if ( pose.aa( n ) == aa_gly ) {
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( -80.0, -10.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, 160.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet(  90.0,  10.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( 100.0, 180.0 ) );
	} else if ( pose.residue_type( n ).is_alpha_aa() ) {
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, -30.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet( -70.0, 140.0 ) );
		main_chain_torsion_set_list.push_back( MainChainTorsionSet(  50.0,  50.0 ) );
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		// AMW: since we haven't made the MainChainTorsionSet a thin wrapper to a vector
		// of reals yet, this acts like setting OMEGA (but will be well interpreted
		// when applied to betas, once we have the chance)
		utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
		dihs[ 1 ] = -140.0;
		dihs[ 2 ] =   60.0;
		dihs[ 3 ] = -120.0;
		main_chain_torsion_set_list.emplace_back( dihs );//-140.0, 60.0, -120.0 ) );
	}
}


/////////////////////////////////////////////////////////////
Real
StepWiseProteinBackboneSampler::get_rotamer_angle( core::Size const i, core::Size const N_SAMPLE ){
	return  ( -180.0 + ( 360.0 / N_SAMPLE ) * i + 0.001 );
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_full( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
	MainChainTorsionSetList & main_chain_torsion_set_list )

{
	//TR << "get_main_chain_torsion_set_list_full rsd " << n << " (" << pose.residue_type( n ).name() << ")" << std::endl;
	//generic -- a residue in the middle of the loop
	if ( pose.residue_type( n ).is_alpha_aa() ) {
		for ( Size i = 1; i <= n_sample_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_ );
			for ( Size j = 1; j <= n_sample_; j++ ) {
				Real const psi = get_rotamer_angle( j, n_sample_ );

				if ( ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  )  > best_energy_cutoff ) {
					continue;
				}

				main_chain_torsion_set_list.emplace_back( phi, psi );
			}
		}
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		// AMW TODO: eventually, just loop over ALL MAINCHAIN INDICES, unifying these two loops. But, rama care.
		for ( Size i = 1; i <= n_sample_beta_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_beta_ );//n_sample_ );
			for ( Size j = 1; j <= n_sample_beta_; j++ ) {
				Real const tht = get_rotamer_angle( j, n_sample_beta_ );//n_sample_ );
				for ( Size k = 1; k <= n_sample_beta_; k++ ) {
					Real const psi = get_rotamer_angle( k, n_sample_beta_ );//n_sample_ );

					// Can't apply a rama filter, oy vey!
					utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
					dihs[ 1 ] = phi;
					dihs[ 2 ] = tht;
					dihs[ 3 ] = psi;
					main_chain_torsion_set_list.emplace_back( dihs );//phi, tht, psi ) );
				}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_n_terminus( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	debug_assert( n == 1);
	if ( pose.residue_type( n ).is_alpha_aa() ) {
		for ( Size j = 1; j <= n_sample_; j++ ) {
			Real const psi = get_rotamer_angle( j, n_sample_ );
			Real best_phi = 0.0;
			Real best_energy = 999999.9999;
			for ( Size i = 1; i <= n_sample_; i++ ) {
				Real const phi = get_rotamer_angle( i, n_sample_ );
				Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				if ( temp_rama  < best_energy ) {
					best_energy = temp_rama;
					best_phi = phi;
				}
			}
			if ( best_energy < best_energy_cutoff ) {
				main_chain_torsion_set_list.emplace_back( best_phi, psi );
			}
		}
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		for ( Size j = 1; j <= n_sample_beta_; j++ ) {
			Real const tht = get_rotamer_angle( j, n_sample_beta_ );
			for ( Size k = 1; k <= n_sample_beta_; k++ ) {
				Real const psi = get_rotamer_angle( k, n_sample_beta_ );

				// AMW TODO: lack of rama can't get best. Just guess.
				Real best_phi = -120.0;
				//Real best_energy = 999999.9999;
				//for ( Size i = 1; i <= n_sample_; i++ ) {
				// Real const phi = get_rotamer_angle( i, n_sample_ );
				// Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				// if ( temp_rama  < best_energy ) {
				//  best_energy = temp_rama;
				//  best_phi = phi;
				// }
				//}
				//if ( best_energy < best_energy_cutoff ) {
				utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
				dihs[ 1 ] = best_phi;
				dihs[ 2 ] = tht;
				dihs[ 3 ] = psi;
				main_chain_torsion_set_list.emplace_back( dihs );//best_phi, tht, psi ) );
				//}
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_c_terminus( core::Size const n, core::pose::Pose const & pose, core::Real const best_energy_cutoff,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	if ( pose.residue_type( n ).is_alpha_aa() ) {
		for ( Size i = 1; i <= n_sample_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_ );

			Real best_psi = 0.0;
			Real best_energy = 999999.9999;
			for ( Size j = 1; j <= n_sample_; j++ ) {
				Real const psi = get_rotamer_angle( j, n_sample_ );
				Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				if ( temp_rama  < best_energy ) {
					best_energy = temp_rama;
					best_psi = psi;
				}
			}
			if ( best_energy < best_energy_cutoff ) {
				main_chain_torsion_set_list.emplace_back( phi, best_psi );
			}
		}
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		for ( Size i = 1; i <= n_sample_beta_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_beta_ );
			for ( Size j = 1; j <= n_sample_beta_; j++ ) {
				Real const tht = get_rotamer_angle( j, n_sample_beta_ );

				// Without rama can't do this. Guess!
				Real best_psi = -120.0;
				//Real best_energy = 999999.9999;
				//for ( Size j = 1; j <= n_sample_; j++ ) {
				// Real const psi = get_rotamer_angle( j, n_sample_ );
				// Real temp_rama = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
				// if ( temp_rama  < best_energy ) {
				//  best_energy = temp_rama;
				//  best_psi = psi;
				// }
				//}
				//if ( best_energy < best_energy_cutoff ) {
				utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
				dihs[ 1 ] = phi;
				dihs[ 2 ] = tht;
				dihs[ 3 ] = best_psi;
				main_chain_torsion_set_list.emplace_back( dihs );
				//}
			}
		}
	}
}

// AMW TODO: have yet to understand the purpose of these functions truly, perhaps there
// will never be a need for sample_theta only. But really we should say
// 'sample lowermost only'

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_sample_phi_only( core::Size const n,
	core::pose::Pose const & pose,
	core::Real const energy_cutoff,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	//TR << "get_main_chain_torsion_set_list_sample_phi_only" << std::endl;
	if ( pose.residue_type( n ).is_alpha_aa() ) {
		// we are prepending and this is the junction residue. sample phi only!
		Real best_rama_energy( 99999.999 ), best_phi( 0.0 );
		Real const psi = pose.psi( n );
		Real const omega = pose.omega( n );
		for ( Size i = 1; i <= n_sample_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_ );
			Real const rama_energy = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
			if ( rama_energy < best_rama_energy || i == 1 ) {
				best_rama_energy = rama_energy;
				best_phi = phi;
			}
			if ( rama_energy  > energy_cutoff )  continue;
			main_chain_torsion_set_list.emplace_back( phi, psi, omega );
		}

		// Make sure to return something at least...
		if ( main_chain_torsion_set_list.size() == 0 )  main_chain_torsion_set_list.push_back( MainChainTorsionSet( best_phi, psi, omega ) );
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		// we are prepending and this is the junction residue. sample phi only!
		// AMW can't judge rama.
		//Real best_rama_energy( 99999.999 ), best_phi( 0.0 );
		Real const psi = pose.psi( n );
		Real const tht = pose.theta( n );
		Real const omega = pose.omega( n );
		for ( Size i = 1; i <= n_sample_beta_; i++ ) {
			Real const phi = get_rotamer_angle( i, n_sample_beta_ );
			//Real const rama_energy = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
			//if ( rama_energy < best_rama_energy || i == 1 ) {
			// best_rama_energy = rama_energy;
			// best_phi = phi;
			//}
			//if ( rama_energy  > energy_cutoff )  continue;
			utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
			dihs[ 1 ] = phi;
			dihs[ 2 ] = tht;
			dihs[ 3 ] = psi;
			main_chain_torsion_set_list.emplace_back( dihs, omega );
		}

		// Make sure to return something at least...
		if ( main_chain_torsion_set_list.size() == 0 )  {
			utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
			dihs[ 1 ] = -140.0;
			dihs[ 2 ] = tht;
			dihs[ 3 ] = psi;
			main_chain_torsion_set_list.emplace_back( dihs, omega );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::get_main_chain_torsion_set_list_sample_psi_only( core::Size const n,
	core::pose::Pose const & pose,
	core::Real const energy_cutoff,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	//TR << "get_main_chain_torsion_set_list_sample_psi_only" << std::endl;
	if ( pose.residue_type( n ).is_alpha_aa() ) {
		//  TR << "JUNCTION RESIDUE --> APPEND " << n << std::endl;

		// we are appending and this is the junction residue. sample psi only!
		Real best_rama_energy( 99999.999 ), best_psi( 0.0 );
		Real const phi = pose.phi( n );
		for ( Size j = 1; j <= n_sample_; j++ ) {
			Real const psi = get_rotamer_angle( j, n_sample_ );
			Real const rama_energy = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
			if ( rama_energy < best_rama_energy || j == 1 ) {
				best_rama_energy = rama_energy;
				best_psi = psi;
			}
			if ( rama_energy  > energy_cutoff )  continue;
			main_chain_torsion_set_list.emplace_back( phi, psi );
			//main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, psi ) );
		}

		// Make sure to return something at least...
		if ( main_chain_torsion_set_list.size() == 0 )  main_chain_torsion_set_list.push_back( MainChainTorsionSet( phi, best_psi ) );
	} else if ( pose.residue_type( n ).is_beta_aa() ) {
		//  TR << "JUNCTION RESIDUE --> APPEND " << n << std::endl;

		// we are appending and this is the junction residue. sample psi only!
		// AMW can't judge rama.
		//Real best_rama_energy( 99999.999 ), best_psi( 0.0 );
		Real const phi = pose.phi( n );
		Real const tht = pose.theta( n );
		for ( Size j = 1; j <= n_sample_beta_; j++ ) {
			Real const psi = get_rotamer_angle( j, n_sample_beta_ );
			//Real const rama_energy = ramachandran_.eval_rama_score_residue( pose.aa( n ), phi, psi  );
			//if ( rama_energy < best_rama_energy || j == 1 ) {
			// best_rama_energy = rama_energy;
			// best_psi = psi;
			//}
			//if ( rama_energy  > energy_cutoff )  continue;
			utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
			dihs[ 1 ] = phi;
			dihs[ 2 ] = tht;
			dihs[ 3 ] = psi;
			main_chain_torsion_set_list.emplace_back( dihs );
			//main_chain_torsion_set_list.push_back( MainChainTorsionSet( dihs ) );
		}

		// Make sure to return something at least...
		if ( main_chain_torsion_set_list.size() == 0 ) {
			utility::fixedsizearray1< Real, 3 > dihs;//{ phi, tht, psi };
			dihs[ 1 ] = phi;
			dihs[ 2 ] = tht;
			dihs[ 3 ] = -120.0;
			main_chain_torsion_set_list.emplace_back( dihs );
			//main_chain_torsion_set_list.push_back( MainChainTorsionSet( dihs ) );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// Coarse -- see set_ss_from_phipsi in core/pose/variant_util.cc.
//////////////////////////////////////////////////////////////////////////////
Size
StepWiseProteinBackboneSampler::get_big_bin( Real const phi, Real const psi ) const
{
	if ( phi < -20.0 && psi > -90.0 && psi < -10.0 ) {
		return 1;
	} else if ( phi < -20.0 && (psi > 20.0 || psi < -170.0) ) {
		return 2;
	}
	return 3;
}


///////////////////////////////////////////////////////////////
// Check set_ss_from_phipsi -- pretty coarse. Could
// later use some guess for strand/helix propensity as an
// energy term or something.
///////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::filter_native_BIG_BINS(
	Size const n,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{
	if ( get_native_pose() == 0 ) {
		utility_exit_with_message(  "Filter based on native big bins, but not native specified?" );
	}

	//big bin not really specified for end residues...
	pose::Pose const & native_pose( *get_native_pose() );

	if ( n == 1 ) return;
	if ( n == native_pose.size() ) return;

	Size const big_bin =  get_big_bin( native_pose.phi(n), native_pose.psi(n) );
	if ( big_bin == 3 ) return; // If loop, no constraints.

	filter_big_bin( big_bin, main_chain_torsion_set_list );
}

//////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::filter_based_on_desired_secstruct(
	char const secstruct,
	MainChainTorsionSetList & main_chain_torsion_set_list )
{

	Size big_bin( 3 );
	if ( secstruct == 'H' ) big_bin = 1;
	if ( secstruct == 'E' ) big_bin = 2;

	// default secstruct is 'L' --> no filtering.
	if ( big_bin == 3 ) return; // If loop, no constraints.

	filter_big_bin( big_bin, main_chain_torsion_set_list );
}

//////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::filter_big_bin( Size const big_bin,
	MainChainTorsionSetList & main_chain_torsion_set_list ) {

	MainChainTorsionSetList main_chain_torsion_set_list_new;

	for ( auto const & main_chain_torsion_set : main_chain_torsion_set_list ) {
		Real const phi = main_chain_torsion_set.phi();
		Real const psi = main_chain_torsion_set.psi();
		if ( ( big_bin == 1 && get_big_bin( phi, psi ) == 1 ) ||
				( big_bin == 2 && get_big_bin( phi, psi ) == 2 ) ) {
			main_chain_torsion_set_list_new.push_back( main_chain_torsion_set );
		}
	}

	// If the list is empty, don't chuck it in!
	if ( !main_chain_torsion_set_list_new.empty() ) {
		main_chain_torsion_set_list = main_chain_torsion_set_list_new;
	}
}

///////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::filter_main_chain_torsion_sets(){

	std::list< std::pair< Real, Size > > energy_index_list;
	for ( Size n = 1; n <= centroid_scores_.size(); n++ ) {
		energy_index_list.push_back( std::make_pair( centroid_scores_[n], n ) ) ;
	}
	energy_index_list.sort();

	utility::vector1< MainChainTorsionSetList > main_chain_torsion_sets_for_moving_residues_new;

	for ( Size n = 1; n <= centroid_scores_.size(); n++ ) {
		if ( n > nstruct_centroid_ ) break;
		main_chain_torsion_sets_for_moving_residues_new.push_back(
			main_chain_torsion_sets_for_moving_residues_[ n ] );
	}

	main_chain_torsion_sets_for_moving_residues_ =  main_chain_torsion_sets_for_moving_residues_new;
}

/////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::sample_cis_omega( MainChainTorsionSetList & main_chain_torsion_set_list ){
	MainChainTorsionSetList main_chain_torsion_set_list_new;

	for ( auto const & main_chain_torsion_set : main_chain_torsion_set_list ) {
		// AMW: this is a bad criterion -- if you were starting from an alpha
		// sample that was EXACTLY cis, this would get fooled. But that's unlikely,
		// right?!
		if ( main_chain_torsion_set.mainchain_dihedral_values()[ 4 ] == 0.0 ) {
			Real const phi = main_chain_torsion_set.phi();
			Real const psi = main_chain_torsion_set.psi();
			main_chain_torsion_set_list_new.push_back( MainChainTorsionSet( phi, psi, 180.0 ) );
			main_chain_torsion_set_list_new.push_back( MainChainTorsionSet( phi, psi,   0.0 ) );
		} else {
			utility::fixedsizearray1< Real, 3 > new_dihs;
			new_dihs[ 1 ] = main_chain_torsion_set.mainchain_dihedral_values()[ 1 ];
			new_dihs[ 2 ] = main_chain_torsion_set.mainchain_dihedral_values()[ 2 ];
			new_dihs[ 3 ] = main_chain_torsion_set.mainchain_dihedral_values()[ 3 ];
			main_chain_torsion_set_list_new.push_back( MainChainTorsionSet( new_dihs, 180.0 ) );
			main_chain_torsion_set_list_new.push_back( MainChainTorsionSet( new_dihs,   0.0 ) );
		}
	}

	main_chain_torsion_set_list = main_chain_torsion_set_list_new;

}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_silent_file( std::string const & silent_file ){
	silent_file_ = silent_file;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_rmsd_cutoff( core::Real const setting ){
	rmsd_cutoff_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_n_sample( core::Size const setting ){
	n_sample_ = setting;
	// This is the only way to manipulate n_sample_beta
	// because it must be always one-third the above.
	n_sample_beta_ = setting/3;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_nstruct_centroid( core::Size const setting ){
	nstruct_centroid_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_centroid_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
	centroid_scorefxn_ = scorefxn;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_filter_native_big_bins( bool const setting ){
	filter_native_big_bins_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_centroid_screen( bool const setting ){
	centroid_screen_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_ghost_loops( bool const setting ){
	ghost_loops_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_apply_vdw_cut( bool const setting ){
	apply_vdw_cut_ = setting;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_centroid_score_diff_cut( core::Real const setting ){
	centroid_score_diff_cut_ = setting;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< MainChainTorsionSetList > const &
StepWiseProteinBackboneSampler::main_chain_torsion_set_lists() const
{
	return main_chain_torsion_sets_for_moving_residues_;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< utility::vector1< core::Real > >
StepWiseProteinBackboneSampler::main_chain_torsion_set_lists_real( Pose const & sampler_pose ) const
{
	utility::vector1< utility::vector1< core::Real > > output_list;
	for ( auto const & main_chain_torsion_set_list : main_chain_torsion_sets_for_moving_residues_ ) {

		utility::vector1< core::Real > output;
		Size moving_res_index = 1;
		for ( auto const & main_chain_torsion_set : main_chain_torsion_set_list ) {

			if ( sampler_pose.residue_type( moving_residues_[ moving_res_index++ ] ).is_beta_aa() ) {
				// AMW: this breaks with alphas, because it tries to interpret it wrong -- because these are
				// always 4 long...
				for ( Real const dih : main_chain_torsion_set.mainchain_dihedral_values() )  {
					output.push_back( dih );
				}
			} else {
				output.push_back( main_chain_torsion_set.phi() );
				output.push_back( main_chain_torsion_set.psi() );
				output.push_back( main_chain_torsion_set.omega() );
			}
		}
		output_list.push_back( output );
	}
	return output_list;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< id::TorsionID >
StepWiseProteinBackboneSampler::which_torsions( Pose const & sampler_pose )
{
	using namespace core::id;
	utility::vector1< id::TorsionID > which_torsions;
	for ( Size const moving_res : moving_residues_ ) {
		// Push back as many bb torsions as there are in this residue type
		// Note: this will get screwed up for terminal types with extra variants, MAYBE.
		// AMW: we can use mainchain_torsions; I'm dumb.
		for ( Size n = 1; n <= sampler_pose.residue( moving_res ).mainchain_torsions().size(); ++n ) {
			which_torsions.push_back( TorsionID( moving_res, BB, n ) );
		}
	}
	return which_torsions;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::initialize_is_fixed_res(){
	set_fixed_residues( working_parameters_->working_fixed_res() );
}

/////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_moving_residues( utility::vector1< Size > const & moving_res ){
	moving_residues_ = moving_res;
}

/////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::set_fixed_residues( utility::vector1< Size > const & fixed_res ){
	is_fixed_res_input_.clear();
	Size const nres = core::pose::rna::remove_bracketed( working_parameters_->working_sequence() ).size();
	for ( Size n = 1; n <= nres; ++n ) {
		is_fixed_res_input_.push_back( false );
	}
	for ( Size i = 1; i <= fixed_res.size(); i++ ) {
		is_fixed_res_input_[ fixed_res[i] ] = true;
	}
}

/////////////////////////////////////////////////////////////
// This is basically deprecated... may not work anymore.
/////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::copy_coords( pose::Pose & pose, pose::Pose const & template_pose, ResMap const & ghost_map ) const
{
	using namespace core::id;
	template_pose.residue( 1 ); // force a refold.

	for ( auto const & elem : ghost_map ) {
		Size const res( elem.first );
		Size const template_res( elem.second );

		for ( Size j = 1; j <= pose.residue_type( res ).natoms(); j++ ) {
			if ( pose.residue_type( res ).atom_name( j ) !=
					template_pose.residue_type( template_res ).atom_name( j ) ) {
				TR << "PROBLEM! " << res << " " << pose.residue( res ).atom_name( j ) <<  "  !=  " <<
					template_res << " " << template_pose.residue( template_res ).atom_name( j ) <<  std::endl;
				utility_exit_with_message( "mismatch in ghost pose" );
			}

			pose.set_xyz( AtomID( j, res ), template_pose.xyz( AtomID( j, template_res ) ) );
		}
	}

	pose.residue( 1 ); // force a refold.
}

///////////////////////////////////////////////////////////////////////////
core::kinematics::FoldTree
StepWiseProteinBackboneSampler::figure_out_fold_tree( ResMap const & ghost_map ) const {

	Size prev_res( 0 ), total_res( 0 );
	utility::vector1< Size > cutpoints;

	for ( auto const & elem : ghost_map ) {
		Size const res( elem.first );
		Size const template_res( elem.second );

		TR << "MAPPING " << res << " --> " << template_res << std::endl;

		if ( (template_res-1) != prev_res  ) cutpoints.push_back( res-1 );
		prev_res = template_res;

		total_res = res;
	}

	core::kinematics::FoldTree f( total_res );

	for ( Size const cut : cutpoints ) {
		TR << "Adding jump across: " << cut << std::endl;
		f.new_jump( cut, cut + 1, cut );
	}

	return f;
}

void
StepWiseProteinBackboneSampler::setup_centroid_screen(
	Real const centroid_score_diff_cut,
	std::string const & centroid_weights,
	Size const nstruct_centroid,
	bool const ghost_loops) {

	using namespace core::scoring;

	ScoreFunctionOP centroid_scorefxn = ScoreFunctionFactory::create_score_function( centroid_weights );
	set_centroid_screen( true );
	set_centroid_score_diff_cut( centroid_score_diff_cut );

	if ( ghost_loops ) {
		// Trying a mode where loops are not included in scoring.
		// the idea was to cut out poses where secondary structure elements
		// are in contact -- energy less than a reference pose in which
		// the secondary structure elements are really far apart.
		set_ghost_loops( true );
		set_centroid_score_diff_cut( 0.0 );
		// if reference is "expanded", rg doesn't make sense.
		centroid_scorefxn->set_weight( rg, 0.0 );
		// disallow steric clashes beyond what it in expanded pose.
		set_apply_vdw_cut( true );
	}

	set_centroid_scorefxn( centroid_scorefxn );
	set_nstruct_centroid( nstruct_centroid );
}


///////////////////////////////////////////////////////////////////////////
// pose with loops excised. centroid level.
void
StepWiseProteinBackboneSampler::prepare_ghost_pose( core::pose::Pose const & pose ){

	using namespace core::kinematics;
	using namespace core::scoring;

	// Figure out sequence of a pose without the loops (i.e., the "moving residues" );
	std::string const full_sequence = pose.sequence();
	ObjexxFCL::FArray1D<bool> moving_array( full_sequence.size(), false );
	for ( Size n = 1; n <= moving_residues_.size(); n++ ) moving_array( moving_residues_[n] ) = true;

	ghost_map_.clear();
	Size count( 0 );
	std::string desired_sequence  = "";
	for ( Size n = 1; n <= full_sequence.size(); n++ ) {
		if ( moving_array(n) ) continue;
		count++;
		desired_sequence += full_sequence[n-1];
		ghost_map_[ count ] = n;
	}

	FoldTree f = figure_out_fold_tree( ghost_map_ ); //chainbreaks, etc.

	TR << "FULL SEQUENCE " << full_sequence << std::endl;
	TR << "DESIRED SEQUENCE " << desired_sequence << std::endl;

	initialize_ghost_pose( ghost_pose_, desired_sequence, pose, ghost_map_, f);
	ghost_pose_->dump_pdb( "GHOST_START.pdb" );

	if ( get_native_pose() ) {
		Pose native_centroid_pose =  *get_native_pose();
		convert_to_centroid( native_centroid_pose );
		initialize_ghost_pose( ghost_native_pose_, desired_sequence, native_centroid_pose, ghost_map_, f);
		ghost_native_pose_->dump_pdb( "GHOST_NATIVE_START.pdb" );
	}

	// Also need to figure out reference energy.
	Pose ghost_pose_blowup = *ghost_pose_;
	for ( Size n = 1; n <= ghost_pose_blowup.num_jump(); n++ ) {
		Jump jump( ghost_pose_blowup.jump( n ) );
		jump.set_translation( Vector( 100.0, 0.0, 0.0 ) ); //This is a little dangerous.
		ghost_pose_blowup.set_jump( n, jump );
	}
	ghost_pose_blowup.dump_pdb( "GHOST_BLOWUP.pdb" );
	runtime_assert( centroid_scorefxn_ != 0 );
	centroid_score_ref_ = (*centroid_scorefxn_)( ghost_pose_blowup );
	centroid_vdw_ref_ = ghost_pose_blowup.energies().total_energies()[ vdw ];

	centroid_scorefxn_->show( TR, ghost_pose_blowup );
	TR << " REFERENCE SCORES " << centroid_score_ref_ << " " << centroid_vdw_ref_ << std::endl;
}

///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::initialize_ghost_pose(
	core::pose::PoseOP & ghost_pose,
	std::string const & desired_sequence,
	core::pose::Pose const & template_pose,
	ResMap const & ghost_map,
	core::kinematics::FoldTree f )
{
	using namespace core::chemical;
	using namespace core::pose;
	static const ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( CENTROID ) );

	ghost_pose = core::pose::PoseOP( new Pose );
	make_pose_from_sequence( *ghost_pose, desired_sequence, *rsd_set );
	copy_coords( *ghost_pose, template_pose, ghost_map );
	ghost_pose->fold_tree( f );
}


///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::convert_to_centroid( core::pose::Pose & pose ) {

	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::N_ACETYLATION, 1 );
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );

	core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::C_METHYLAMIDATION, pose.size() );
	core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, pose.size() );

	core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID_t );

	// get DSSP, assign secondary structure
	core::scoring::dssp::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	//pose.dump_pdb( "CENTROID.pdb" );
}


///////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinBackboneSampler::filter_and_save( core::pose::Pose & pose,
	std::string const & tag  )
{
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::pose;

	if ( get_native_pose() && rmsd_cutoff_ > 0.0 ) {
		Real const rmsd = rmsd_with_super( pose, *get_native_pose(), is_protein_backbone_including_O );
		TR << "CHECKING: decoy " << tag << ": " <<  rmsd << " vs. filter " << rmsd_cutoff_ << std::endl;
		if ( rmsd > rmsd_cutoff_ ) return;
	}

	Real centroid_score( 0.0 );
	if ( centroid_screen_ ) {
		if ( ghost_loops_ ) {
			//this may be a little inefficient, since all internal dofs needs to be recalculated. Its safe, though, I think.
			copy_coords( *ghost_pose_, pose, ghost_map_ );

			centroid_score = (*centroid_scorefxn_)( *ghost_pose_ ); // pose with loops excised

			Real const centroid_vdw = ghost_pose_->energies().total_energies()[ vdw ];
			if ( apply_vdw_cut_ && centroid_vdw > centroid_vdw_ref_ ) return;

		} else {
			centroid_score = (*centroid_scorefxn_)( pose );
			// Keep running tabs on best score seen so far.
			if ( centroid_score < centroid_score_ref_ ) centroid_score_ref_ = centroid_score;
		}

		Real const centroid_score_diff = centroid_score - centroid_score_ref_;
		if (  centroid_score_diff >  centroid_score_diff_cut_ )  return;

		//TR << "COMPARING " << centroid_score << " " << centroid_score_ref_ << std::endl;
	}

	main_chain_torsion_sets_for_moving_residues_.push_back( main_chain_torsion_set_for_moving_residues_ );
	centroid_scores_.push_back( centroid_score );
}


} //protein
} //modeler
} //stepwise
} //protocols
