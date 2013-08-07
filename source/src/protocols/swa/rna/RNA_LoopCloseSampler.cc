// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA_LoopCloseSampler
/// @brief Loop Close Sampler...
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/RNA_LoopCloseSampler.hh>
#include <protocols/swa/rna/RNA_AnalyticLoopCloser.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Conformation.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh> // !! for monte carlo, new 2013.
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>

#include <time.h>

#include <string>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>


using namespace core;
using core::Real;

static numeric::random::RandomGenerator RG( 199123 );  // <- Magic number, do not change it!
static basic::Tracer TR( "protocols.swa.rna.RNA_LoopCloseSampler" );

namespace protocols {
namespace swa {
namespace rna {

//////////////////////////////////////////////////////////////////////////
//constructor!
RNA_LoopCloseSampler::RNA_LoopCloseSampler ( Size const moving_suite, Size const chainbreak_suite ) :
	moving_suite_ ( moving_suite ),
	chainbreak_suite_ ( chainbreak_suite ),
	scorefxn_ ( core::scoring::getScoreFunction() ),
	bin_size_ ( 20 ),
	n_construct_ ( 0 ),
	rep_cutoff_ ( 0.1 ),
	torsion_range_ ( 20.0 ),
	torsion_increment_ ( 5.0 ),
	epsilon_range_ ( 40.0 ),
	just_output_score_ ( false ),
	sample_only_ ( false ),
	include_current_( false ),
	sample_native_torsion_ ( false ),
	choose_random_( false ){
	RNA_AnalyticLoopCloser rna_analytic_loop_closer ( moving_suite, chainbreak_suite );
	initialize_rep_scorefxn();
}

//////////////////////////////////////////////////////////////////////////
//destructor
RNA_LoopCloseSampler::~RNA_LoopCloseSampler()
{}

////////////////////////////////////////////////////////////////////////////////////////////
// iterate over 4 dofs:
//     epsilon1, zeta1, alpha1;     alpha2
//
// 	solve for the other 6 by chain closure:
//     beta1, gamma1;               epsilon2, zeta2, beta2, gamma2
//
// This could probably be encapsulated into a PoseSampleGenerator or something similar.
////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::apply ( core::pose::Pose & pose ) {
	using namespace pose;
	using namespace chemical;
	using namespace id;
	using namespace chemical::rna;
	using namespace numeric::conversions;

	Real const fa_rep_score_baseline = initialize_fa_rep ( pose, utility::tools::make_vector1 ( moving_suite_ ), rep_scorefxn_ );

	PuckerState pucker_state1 = Get_residue_pucker_state ( pose, moving_suite_ );

	// following is only used if we are estimating jacobians (numerically).
	utility::vector1< utility::vector1< utility::vector1< Real > > > perturbed_solution_torsions;
	utility::vector1< utility::vector1< Real > > J; // 6 x 6 Jacobian.
	utility::vector1< Real > six_zeros;
	for ( int i = 1; i <= 6; i++ ) six_zeros.push_back ( 0.0 );
	for ( int i = 1; i <= 6; i++ ) J.push_back ( six_zeros );

	// The closer.
	RNA_AnalyticLoopCloser rna_analytic_loop_closer ( moving_suite_, chainbreak_suite_ );

	// set up backbones to close.
	// epsilon1, zeta1, alpha1; and alpha2 are 'driver torsions' -- deltas will be fixed,
	//  and the other six torsions will be solved by analytical loop closure.
	utility::vector1< TorsionID > driver_torsion_IDs;
	driver_torsion_IDs.push_back( TorsionID ( moving_suite_ , id::BB, EPSILON ) );
	driver_torsion_IDs.push_back( TorsionID ( moving_suite_ , id::BB, ZETA ) );
	driver_torsion_IDs.push_back( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ) );
	driver_torsion_IDs.push_back( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ) );

	utility::vector1< utility::vector1< Real > > driver_torsion_sets;

	Real epsilon1_center = ( pucker_state1 == NORTH ) ? -150.17 : -98.45;
	Real epsilon1_min = epsilon1_center - epsilon_range_;
	Real epsilon1_max = epsilon1_center + epsilon_range_;
	if ( pucker_state1 == SOUTH && epsilon1_min > -178 ) epsilon1_min = -178.45;
	Real epsilon1_increment = bin_size_;
	Real zeta1_min = 20.0;
	Real zeta1_max = 340.0 - bin_size_;
	Real zeta1_increment = bin_size_;
	Real alpha1_min = 20.0;
	Real alpha1_max = 340.0 - bin_size_;
	Real alpha1_increment = bin_size_;
	Real alpha2_min = 20.0;
	Real alpha2_max = 340.0 - bin_size_;
	Real alpha2_increment = bin_size_;

	for ( Real epsilon1 = epsilon1_min; epsilon1 <= epsilon1_max; epsilon1 += epsilon1_increment ) {
		for ( Real zeta1 = zeta1_min; zeta1 <= zeta1_max; zeta1 += zeta1_increment ) {
			for ( Real alpha1 = alpha1_min; alpha1 <= alpha1_max; alpha1 += alpha1_increment ) {
				for ( Real alpha2 = alpha2_min; alpha2 <= alpha2_max; alpha2 += alpha2_increment ) {
					utility::vector1< Real > driver_torsion_set = utility::tools::make_vector1( epsilon1, zeta1, alpha1, alpha2 );
					driver_torsion_sets.push_back( driver_torsion_set );
				}
			}
		}
	}

	if ( include_current_ ) add_driver_torsion_set( driver_torsion_sets, driver_torsion_IDs, pose );
	if ( sample_native_torsion_ && get_native_pose() ) add_driver_torsion_set( driver_torsion_sets, driver_torsion_IDs, *get_native_pose() );

	PoseCOP native_pose = get_native_pose();
	if ( sample_native_torsion_ ) utility_exit_with_message( "sample_native_torsion_ not set up properly in RNA_LoopCloseSampler.cc!" );

	for ( Size count = 1; count <= driver_torsion_sets.size(); count++ ){

		utility::vector1< Real > driver_torsion_set = driver_torsion_sets[ count ];
		if ( choose_random_ ) driver_torsion_set = RG.random_element( driver_torsion_sets );

		runtime_assert( driver_torsion_set.size() == 4 );
		for ( Size  k = 1; k <= driver_torsion_set.size(); k++ ) pose.set_torsion( driver_torsion_IDs[k], driver_torsion_set[k] );

		//close loop.
		rna_analytic_loop_closer.apply ( pose );

		// iterate over solutions -- anything with OK repulsive?
		for ( Size n = 1; n <= rna_analytic_loop_closer.nsol(); n++ ) {
			rna_analytic_loop_closer.fill_solution ( pose, n );

			if ( !torsion_angles_within_cutoffs ( pose, moving_suite_, chainbreak_suite_ ) ) continue;

			if ( !sample_only_ ) {
				if ( !check_clash ( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn_ ) ) continue;
			}

			// save data
			n_construct_++;
			torsion_info_.clear();
			torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ , id::BB, EPSILON ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ , id::BB, ZETA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, BETA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( moving_suite_ + 1, id::BB, GAMMA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ , id::BB, EPSILON ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ , id::BB, ZETA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, BETA ) ) );
			torsion_info_.push_back ( pose.torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, GAMMA ) ) );
			all_torsion_info_.push_back ( torsion_info_ );

			if ( choose_random_ && all_torsion_info_.size() > 1 ) break; // our work is done!

		}

		if ( choose_random_ && all_torsion_info_.size() > 1 ) break; // our work is done!
	}

}

///////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::add_driver_torsion_set(
														utility::vector1< utility::vector1< Real > > & driver_torsion_sets,
														utility::vector1< core::id::TorsionID > const & driver_torsion_IDs,
														pose::Pose const & pose ) {
	utility::vector1< Real > current_driver_torsion_set;
	for ( Size k = 1; k <= driver_torsion_IDs.size(); k++ )	 current_driver_torsion_set.push_back( pose.torsion( driver_torsion_IDs[k] ) );
	driver_torsion_sets.push_back( current_driver_torsion_set );
}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::fill_pose ( pose::Pose & pose, Size construct_number ) {
	using namespace pose;
	using namespace id;
	using namespace chemical::rna;
	utility::vector1< Real > & torsion_info = all_torsion_info_[construct_number];
	pose.set_torsion ( TorsionID ( moving_suite_ , id::BB, EPSILON ), torsion_info[1] );
	pose.set_torsion ( TorsionID ( moving_suite_ , id::BB, ZETA ), torsion_info[2] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, ALPHA ), torsion_info[3] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, BETA ), torsion_info[4] );
	pose.set_torsion ( TorsionID ( moving_suite_ + 1, id::BB, GAMMA ), torsion_info[5] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ , id::BB, EPSILON ), torsion_info[6] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ , id::BB, ZETA ), torsion_info[7] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, ALPHA ), torsion_info[8] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, BETA ), torsion_info[9] );
	pose.set_torsion ( TorsionID ( chainbreak_suite_ + 1, id::BB, GAMMA ), torsion_info[10] );
	return;
}
///////////////////////////////
void
RNA_LoopCloseSampler::clear_all() {
	all_torsion_info_.clear();
	n_construct_ = 0;
	return;
}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::initialize_rep_scorefxn() {
	using namespace core::scoring;
	rep_scorefxn_ = new ScoreFunction;
	rep_scorefxn_->set_weight ( fa_rep, 0.12 );
}

///////////////////////////////////////////////////////////////
Real
RNA_LoopCloseSampler::initialize_fa_rep ( pose::Pose const & pose,
    utility::vector1< Size > const & moving_suites,
    scoring::ScoreFunctionOP rep_scorefxn ) {
	using namespace pose;
	using namespace kinematics;
	using namespace scoring;
	using namespace protocols::swa;

	if ( sample_only_ ) {
		return 0;
	}

	Pose pose_expand = pose;

	for ( Size n = 1; n <= moving_suites.size(); n++ ) {
		Size const jump_at_moving_suite = make_cut_at_moving_suite ( pose_expand, moving_suites[n] );
		Jump j = pose_expand.jump ( jump_at_moving_suite );
		j.set_translation ( Vector ( 1.0e4 * n, 0.0, 0.0 ) );
		pose_expand.set_jump ( jump_at_moving_suite, j );
	}

	( *rep_scorefxn ) ( pose_expand );
	EnergyMap const & energy_map = pose_expand.energies().total_energies();
	return energy_map[ fa_rep ] * rep_scorefxn->get_weight ( fa_rep );
}


///////////////////////////////////////////////////////////////
bool
RNA_LoopCloseSampler::check_clash ( pose::Pose & pose,
                                    Real const & fa_rep_score_baseline,
                                    Real const & rep_cutoff_,
                                    scoring::ScoreFunctionOP rep_scorefxn ) {
	using namespace scoring;
	( *rep_scorefxn ) ( pose );
	EnergyMap const & energy_map = pose.energies().total_energies();
	Real const fa_rep_score = energy_map[ fa_rep ] * rep_scorefxn->get_weight ( fa_rep );
	//	TR << fa_rep_score << " " << fa_rep_score_baseline << std::endl;

	if ( ( fa_rep_score - fa_rep_score_baseline ) > rep_cutoff_ ) return false;

	static Real const tolerance ( 1.0e-3 );

	if ( ( fa_rep_score - fa_rep_score_baseline ) < -1.0 * tolerance ) {
		TR << fa_rep_score << " " << fa_rep_score_baseline << std::endl;
		//		utility_exit_with_message( "Weird fa_rep?" );
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////
void
RNA_LoopCloseSampler::set_scorefxn ( core::scoring::ScoreFunctionOP const & scorefxn ) {
	scorefxn_ = scorefxn;
}

//////////////////////////////////////////////////////////////////////////
bool
RNA_LoopCloseSampler::torsion_angles_within_cutoffs ( pose::Pose const & pose,
    Size const moving_suite,
    Size const chainbreak_suite ) {
	using namespace id;
	using namespace core::chemical::rna;
	using namespace protocols::swa::rna;

	static Size const epsilonmin_n = 155, epsilonmax_n = 310;
	static Size const epsilonmin_s = 175, epsilonmax_s = 310;
	static Size const gammapmin  =  20, gammapmax  =  95;
	static Size const gammatmin  = 140, gammatmax  = 215;
	static Size const gammammin  = 260, gammammax  = 335;
	static Size const betamin    =  50, betamax    = 290;
	static Size const alphamin   =  25, alphamax   = 335;
	static Size const zetamin    =  25, zetamax    = 335;

	//Beta
	Real beta1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite + 1, id::BB, BETA ) ) );
	if ( beta1 < 0 ) beta1 += 360;
	if ( beta1 < betamin || beta1 > betamax ) return false;

	Real beta2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite + 1, id::BB, BETA ) ) );
	if ( beta2 < 0 ) beta2 += 360;
	if ( beta2 < betamin || beta2 > betamax ) return false;

	//Gamma
	Real gamma1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite + 1, id::BB, GAMMA ) ) );
	if ( gamma1 < 0 ) gamma1 += 360;
	if ( ( gamma1 < gammapmin || gamma1 > gammapmax ) &&
			 ( gamma1 < gammatmin || gamma1 > gammatmax ) &&
			 ( gamma1 < gammammin || gamma1 > gammammax ) ) {
		return false;
	}

	Real gamma2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite + 1, id::BB, GAMMA ) ) );
	if ( gamma2 < 0 ) gamma2 += 360;
	if ( ( gamma2 < gammapmin || gamma2 > gammapmax ) &&
			 ( gamma2 < gammatmin || gamma2 > gammatmax ) &&
			 ( gamma2 < gammammin || gamma2 > gammammax ) ) {
		return false;
	}

	//Epsilon
	PuckerState pucker_state1 = Get_residue_pucker_state ( pose, chainbreak_suite );
	Real epsilon1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite, id::BB, EPSILON ) ) );
	if ( epsilon1 < 0 ) epsilon1 += 360;
	if ( pucker_state1 == NORTH ) {
		if ( epsilon1 < epsilonmin_n || epsilon1 > epsilonmax_n ) return false;
	} else {
		if ( epsilon1 < epsilonmin_s || epsilon1 > epsilonmax_s ) return false;
	}

	PuckerState pucker_state2 = Get_residue_pucker_state ( pose, chainbreak_suite );
	Real epsilon2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite, id::BB, EPSILON ) ) );
	if ( epsilon2 < 0 ) epsilon2 += 360;
	if ( pucker_state2 == NORTH ) {
		if ( epsilon2 < epsilonmin_n || epsilon2 > epsilonmax_n ) return false;
	} else {
		if ( epsilon2 < epsilonmin_s || epsilon2 > epsilonmax_s ) return false;
	}

	//Zeta
	Real zeta1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite, id::BB, ZETA ) ) );
	if ( zeta1 < 0 ) zeta1 += 360;
	if ( zeta1 < zetamin || zeta1 > zetamax ) return false;

	Real zeta2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite, id::BB, ZETA ) ) );
	if ( zeta2 < 0 ) zeta2 += 360;
	if ( zeta2 < zetamin || zeta1 > zetamax ) return false;

	//Alpha
	Real alpha1 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( moving_suite + 1, id::BB, ALPHA ) ) );
	if ( alpha1 < 0 ) alpha1 += 360;
	if ( alpha1 < alphamin || alpha1 > alphamax ) return false;

	Real alpha2 = numeric::principal_angle_degrees ( pose.torsion ( TorsionID ( chainbreak_suite + 1, id::BB, ALPHA ) ) );
	if ( alpha2 < 0 ) alpha2 += 360;
	if ( alpha2 < alphamin || alpha2 > alphamax ) return false;

	return true;
}


}
}
}
