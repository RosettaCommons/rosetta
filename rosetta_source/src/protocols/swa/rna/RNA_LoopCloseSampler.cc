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
#include <protocols/swa/rna/determinants.hh>
#include <protocols/swa/rna/rigid_body_settings.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_TorsionPotential.hh>
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
// AUTO-REMOVED #include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <numeric/conversions.hh>
// AUTO-REMOVED #include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/exit.hh>
// AUTO-REMOVED #include <time.h>

#include <string>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <sstream>
#include <map>

//Auto Headers
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <utility/vector1.hh>



using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.rna.rna_loop_close_sampler" ) ;

namespace protocols {
namespace swa {
namespace rna {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  RNA_LoopCloseSampler::RNA_LoopCloseSampler( Size const moving_suite, Size const chainbreak_suite ):
		moving_suite_( moving_suite ),
		chainbreak_suite_( chainbreak_suite ),
		scorefxn_( core::scoring::getScoreFunction() ),
		bin_size_( 20 ),
		epsilon_range_( 40.0 ),
		rep_cutoff_( 0.1 ),
		torsion_range_( 20.0 ),
		torsion_increment_( 5.0 ),
		center_around_native_( false ),
		calculate_jacobian_( false ),
		save_torsion_info_( false ),
		just_output_score_( false ),
		silent_file_( "" )
  {
		RNA_AnalyticLoopCloser rna_analytic_loop_closer( moving_suite, chainbreak_suite );
		initialize_rep_scorefxn();
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  RNA_LoopCloseSampler::~RNA_LoopCloseSampler()
  {}

	/////////////////////
	std::string
	RNA_LoopCloseSampler::get_name() const {
		return "RNA_LoopCloseSampler";
	}


	//////////////////////////////////////////////////////////////////////////
	void
	RNA_LoopCloseSampler::apply( core::pose::Pose & pose ) {

		using namespace pose;
		using namespace scoring;
		using namespace io::silent;
		using namespace chemical;
		using namespace id;
		using namespace scoring::rna;
		using namespace protocols::swa;
		using namespace protocols::swa::rna;
		using namespace numeric::conversions;

		Real const fa_rep_score_baseline = initialize_fa_rep( pose, utility::tools::make_vector1( moving_suite_ ), rep_scorefxn_ );

		sfd_ = new core::io::silent::SilentFileData;

		// iterate over 4 dofs:
		//     epsilon1, zeta1, alpha1;    alpha2
		// 	solve for the other 6 by chain closure:
		//     beta1, gamma1;               epsilon2, rho2, beta2, gamma2
		// This should probably be encapsulated into a PoseSampleGenerator or something similar
		int const bins1_ = 360/bin_size_ ; //This is total bins, default is 18
		int const bins2_ = bins1_/2; //This is total bins divided by 2; default is 9
		int const bins3_ = bins1_/3; //This is total bins divided by 3; default is 6
		int const bins4_ = 1 + 40/bin_size_; //This is the bin for chi and episilon, these two torsion angles vary from -20+mean to 20+mean
		PuckerState pucker_state1 = Get_residue_pucker_state( pose, moving_suite_ );
		Real epsilon1, zeta1, alpha1( 0.0 ), alpha2( 0.0 ), perturb_epsilon1, perturb_zeta1, perturb_alpha1, perturb_alpha2;
		Real beta1, beta2, epsilon2;
		Size count( 0 );

		// following is only used if we are estimating jacobians (numerically).
		utility::vector1< utility::vector1< utility::vector1< Real > > > perturbed_solution_torsions;
		utility::vector1< utility::vector1< Real > > J; // 6 x 6 Jacobian.
		utility::vector1< Real > six_zeros;
		for ( int i = 1; i <= 6; i++ ) six_zeros.push_back( 0.0 );
		for ( int i = 1; i <= 6; i++ ) J.push_back( six_zeros );
		Real const perturbation_size( 1.0e-5 );
		Real const detJ_cutoff_( 1.0 );

		Real epsilon1_center = (pucker_state1 == NORTH) ? -150.17 : -98.45;
		//		std::cout << "PUCKER_STATE1 " << pucker_state1 << ' ' << pose.torsion( id::TorsionID( moving_suite_, id::BB, scoring::rna::DELTA ) ) <<  "  IDEAL EPSILON: " << epsilon1_center << " +/- " << epsilon_range_ << std::endl;
		Real epsilon1_min = epsilon1_center - epsilon_range_;
		Real epsilon1_max = epsilon1_center + epsilon_range_;
		Real epsilon1_increment = bin_size_;

		Real alpha1_min = 0.0;
		Real alpha1_max = 360.0-bin_size_;
		Real alpha1_increment = bin_size_;

		Real alpha2_min = 0.0;
		Real alpha2_max = 360.0-bin_size_;
		Real alpha2_increment = bin_size_;

		Real zeta1_min = 0.0;
		Real zeta1_max = 360.0-bin_size_;
		Real zeta1_increment = bin_size_;

		if ( calculate_jacobian_ && rbs_new_pair_.size() != 6 ) {
			utility_exit_with_message( "To calculate jacobian, for now need to specify 'rigid body setting' -- to be fixed in the future." );
		}

		// move this to its own function?
		PoseCOP native_pose = get_native_pose();
		if ( center_around_native_ ){
			if (!native_pose) utility_exit_with_message( "must supply -native" );

			epsilon1_center            = pose.torsion( TorsionID( moving_suite_  , id::BB, EPSILON ) );
			Real const zeta1_center    = pose.torsion( TorsionID( moving_suite_  , id::BB, ZETA ) );
			Real const alpha1_center   = pose.torsion( TorsionID( moving_suite_+1, id::BB, ALPHA ) );
			Real const alpha2_center   = pose.torsion( TorsionID( chainbreak_suite_+1, id::BB, ALPHA ) );

			epsilon1_min = epsilon1_center - torsion_range_;
			epsilon1_max = epsilon1_center + torsion_range_;
			epsilon1_increment = torsion_increment_;

			alpha1_min = alpha1_center - torsion_range_;
			alpha1_max = alpha1_center + torsion_range_;
			alpha1_increment = torsion_increment_;

			zeta1_min = zeta1_center - torsion_range_;
			zeta1_max = zeta1_center + torsion_range_;
			zeta1_increment = torsion_increment_;

			alpha2_min = alpha2_center - torsion_range_;
			alpha2_max = alpha2_center + torsion_range_;
			alpha2_increment = torsion_increment_;
		}

		RNA_AnalyticLoopCloser rna_analytic_loop_closer( moving_suite_, chainbreak_suite_ );

		for (Real epsilon1 = epsilon1_min; epsilon1 <= epsilon1_max; epsilon1 += epsilon1_increment ){
			for (Real alpha1 = alpha1_min; alpha1 <= alpha1_max; alpha1 += alpha1_increment ){
				for (Real alpha2 = alpha2_min; alpha2 <= alpha2_max; alpha2 += alpha2_increment ){
					for (Real zeta1 = zeta1_min; zeta1 <= zeta1_max; zeta1 += zeta1_increment ){


						pose.set_torsion( TorsionID( moving_suite_,       id::BB, EPSILON ), epsilon1 );
						pose.set_torsion( TorsionID( moving_suite_,       id::BB, ZETA ),    zeta1 );
						pose.set_torsion( TorsionID( moving_suite_+1,     id::BB, ALPHA ),   alpha1  );
						pose.set_torsion( TorsionID( chainbreak_suite_+1, id::BB, ALPHA ),   alpha2  );

						//close loop.
						rna_analytic_loop_closer.apply( pose );

						bool perturb_solutions_calculated = false;
						// iterate over solutions -- anything with OK repulsive?
						for ( Size n = 1; n <= rna_analytic_loop_closer.nsol(); n++ ){

							rna_analytic_loop_closer.fill_solution( pose, n );
							utility::vector1< Real > const & solution_torsions = rna_analytic_loop_closer.get_torsions( n );

							if ( !torsion_angles_within_cutoffs( pose, moving_suite_, chainbreak_suite_, bin_size_, bins2_ ) ) continue;

							if ( !check_clash( pose, fa_rep_score_baseline, rep_cutoff_, rep_scorefxn_ ) ) continue;

							if ( center_around_native_  && !check_close_to_native( pose, *native_pose, moving_suite_, chainbreak_suite_, torsion_range_) ) continue;

							// if OK, score and save it.
							Real const score = ( *scorefxn_ )( pose );

							// Jacobian --> phase space transformation factor.
							// Calculate torsions upon six perturbations of rigid body settings
							Real detJ( 0.0 ), volume( 0.0 );

							// move Jacobian stuff to its own function?
							if ( calculate_jacobian_ ) {
								if ( !perturb_solutions_calculated ){
									calculate_perturbed_solutions( pose,
																								 moving_suite_, chainbreak_suite_,
																								 rbs_new_pair_,
																								 perturbation_size,
																								 perturbed_solution_torsions );
									bool ok( true );
									for ( Size m = 1; m <= 6; m++ ){
										if ( perturbed_solution_torsions[m].size() != rna_analytic_loop_closer.nsol() ) {
											std::cout << "WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!" << std::endl;
											std::cout << "Number of solutions ORIGINAL: " <<  rna_analytic_loop_closer.nsol() << std::endl;
											std::cout << "After perturbation in direction " << m << ": " << perturbed_solution_torsions[m].size() << std::endl;
										std::cout << " --> skipping this one." << std::endl;
										std::cout << "WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!" << std::endl;
										ok = false;
										//													utility_exit_with_message( "Numerical problem in calculating jacobian" );
										}
									}
									if ( !ok ) break;

									perturb_solutions_calculated = true;
								}

								// for this particular solution, find which of the perturbed solutions are little nudges away.
								// calculate determinant of jacobian, and then fold this into the volume element.
								get_jacobian( solution_torsions, perturbed_solution_torsions, perturbation_size, J );

								detJ =  std::abs( get_determinant( J ) );
								if( detJ > detJ_cutoff_ ){
									std::cout << "WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!" << std::endl;
									std::cout << " det J " << detJ << " is higher than cutoff " << detJ_cutoff_ << " --> skipping it. " << std::endl;
									std::cout << "WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!WARNING!" << std::endl;
									continue;
								}
								volume = detJ * radians(epsilon1_increment) * radians(zeta1_increment) * radians(alpha1_increment) * radians(alpha2_increment);
								//						std::cout << "detJ " << detJ << "   bin_size " << bin_size_ << "    volume " << volume << std::endl;
							}

							// save data to outfile
							count++;

							if ( save_torsion_info_ ){
								torsion_info_.clear();
								torsion_info_.push_back( pose.torsion( TorsionID( moving_suite_  , id::BB, EPSILON ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( moving_suite_  , id::BB, ZETA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( moving_suite_+1, id::BB, ALPHA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( moving_suite_+1, id::BB, BETA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( moving_suite_+1, id::BB, GAMMA ) ) );

								torsion_info_.push_back( pose.torsion( TorsionID( chainbreak_suite_  , id::BB, EPSILON ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( chainbreak_suite_  , id::BB, ZETA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( chainbreak_suite_+1, id::BB, ALPHA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( chainbreak_suite_+1, id::BB, BETA ) ) );
								torsion_info_.push_back( pose.torsion( TorsionID( chainbreak_suite_+1, id::BB, GAMMA ) ) );

								torsion_info_.push_back( score );
								torsion_info_.push_back( volume );

								all_torsion_info_.push_back( torsion_info_ );
							}


							std::string const tag = "S_" + ObjexxFCL::lead_zero_string_of( count, 6 );
							BinaryRNASilentStruct s( pose, tag );

							if ( save_torsion_info_ ) {
								s.add_energy( "epsilon1", torsion_info_[1] );
								s.add_energy( "zeta1"   , torsion_info_[2] );
								s.add_energy( "alpha1"  , torsion_info_[3] );
								s.add_energy( "beta1"   , torsion_info_[4] );
								s.add_energy( "gamma1"  , torsion_info_[5] );
								s.add_energy( "epsilon2", torsion_info_[6] );
								s.add_energy( "zeta2"   , torsion_info_[7] );
								s.add_energy( "alpha2"  , torsion_info_[8] );
								s.add_energy( "beta2"   , torsion_info_[9] );
								s.add_energy( "gamma2"  , torsion_info_[10] );
								s.add_energy( "log_vol" , log( volume ) );
							}

							if ( native_pose ) s.add_energy( "all_rms", all_atom_rmsd( pose,*native_pose) );

							if ( silent_file_.size() > 0 ) sfd_->write_silent_struct( s, silent_file_, just_output_score_ );

							sfd_->add_structure( s );

						}
					}
				}
			}
		}


	}


	//////////////////////////////////////////////////////////////////////////
	void
	RNA_LoopCloseSampler::initialize_rep_scorefxn(){
		using namespace core::scoring;
		rep_scorefxn_ = new ScoreFunction;
		rep_scorefxn_->set_weight( fa_rep, 0.12 );
	}

	///////////////////////////////////////////////////////////////
	Real
	RNA_LoopCloseSampler::initialize_fa_rep( pose::Pose const & pose,
																					 utility::vector1< Size > const & moving_suites,
																					 scoring::ScoreFunctionOP rep_scorefxn ) {

		using namespace pose;
		using namespace kinematics;
		using namespace scoring;
		using namespace protocols::swa;

		Pose pose_expand = pose;

		for ( Size n = 1; n <= moving_suites.size(); n++ ){
			Size const jump_at_moving_suite = make_cut_at_moving_suite( pose_expand, moving_suites[n] );
			Jump j = pose_expand.jump( jump_at_moving_suite );
			j.set_translation( Vector( 1.0e4 * n, 0.0, 0.0 ) );
			pose_expand.set_jump( jump_at_moving_suite, j );
		}

		(*rep_scorefxn)( pose_expand );
		EnergyMap const & energy_map=pose_expand.energies().total_energies();
		return energy_map[ fa_rep ] * rep_scorefxn->get_weight( fa_rep );

	}


	///////////////////////////////////////////////////////////////
	bool
	RNA_LoopCloseSampler::check_clash( pose::Pose & pose,
																		 Real const & fa_rep_score_baseline,
																		 Real const & rep_cutoff_,
																		 scoring::ScoreFunctionOP rep_scorefxn ){

		using namespace scoring;

		(*rep_scorefxn)( pose );
		EnergyMap const & energy_map=pose.energies().total_energies();
		Real const fa_rep_score = energy_map[ fa_rep ] * rep_scorefxn->get_weight( fa_rep );

		//	std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;

		if ( (fa_rep_score - fa_rep_score_baseline) > rep_cutoff_ ) return false;

		static Real const tolerance( 1.0e-3 );
		if ( (fa_rep_score - fa_rep_score_baseline) < -1.0 * tolerance ) {
		std::cout << fa_rep_score << " " << fa_rep_score_baseline << std::endl;
		//		utility_exit_with_message( "Weird fa_rep?" );
		}

		return true;
	}

  //////////////////////////////////////////////////////////////////////////
  void
  RNA_LoopCloseSampler::set_silent_file( std::string const & silent_file){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	RNA_LoopCloseSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP
	RNA_LoopCloseSampler::silent_file_data(){
		return sfd_;
	}

	///////////////////////////////////////////////////////////////
	bool
	RNA_LoopCloseSampler::torsion_angles_within_cutoffs( pose::Pose const & pose,
																											 Size const moving_suite,
																											 Size const chainbreak_suite,
																											 int const bin_size_,
																											 int const bins2_ ) {

		using namespace id;
		using namespace core::scoring::rna;
		using namespace protocols::swa::rna;

		// Quick check on range of epsilon, beta. Other torsion angles span full 360.0 range.
		Real const beta1 = numeric::principal_angle_degrees( pose.torsion( TorsionID( moving_suite+1, id::BB, BETA ) ) );
		if ( ( beta1 < (  180.0 - bin_size_ * ( bins2_/2.0 - 0.5 ) ) ) &&
				 ( beta1 > ( -180.0 + bin_size_ * ( bins2_/2.0 - 0.5 ) ) ) ) {
			return false;
		}

		Real const beta2 = numeric::principal_angle_degrees( pose.torsion( TorsionID( chainbreak_suite+1, id::BB, BETA ) ) );
		if ( ( beta2 < (  180.0 - bin_size_ * ( bins2_/2.0 - 0.5 ) ) ) &&
				 ( beta2 > ( -180.0 + bin_size_ * ( bins2_/2.0 - 0.5 ) ) ) ) {
			return false;
		}

		PuckerState pucker_state1 = Get_residue_pucker_state( pose, moving_suite );
		Real const epsilon1 = numeric::principal_angle_degrees( pose.torsion( TorsionID( moving_suite, id::BB, EPSILON ) ) );
		Real const epsilon1_ideal= (pucker_state1 == NORTH) ? -150.17 : -98.45;
		if ( ( epsilon1 < (epsilon1_ideal - epsilon_range_ ) )||
				 ( epsilon1 > (epsilon1_ideal + epsilon_range_ ) ) ) {
			return false;
		}

		PuckerState pucker_state2 = Get_residue_pucker_state( pose, chainbreak_suite );
		Real const epsilon2 = numeric::principal_angle_degrees( pose.torsion( TorsionID( chainbreak_suite, id::BB, EPSILON ) ) );
		Real const epsilon2_ideal= (pucker_state2 == NORTH) ? -150.17 : -98.45;
		if ( ( epsilon2 < (epsilon2_ideal - epsilon_range_ ) )||
				 ( epsilon2 > (epsilon2_ideal + epsilon_range_ ) ) ) {
			return false;
		}

		return true;

	}


///////////////////////////////////////////////////////////////
	bool
	RNA_LoopCloseSampler::check_close_to_native( pose::Pose const & pose1, pose::Pose const & pose2, Size const suite1, Size const suite2, Real const angle_range ){

		using namespace id;
		using namespace scoring::rna;
		using namespace numeric::conversions;

		utility::vector1< TorsionID> torsion_ids;
		torsion_ids.push_back( TorsionID( suite1,   id::BB, EPSILON ) ) ;
		torsion_ids.push_back( TorsionID( suite1,   id::BB, ZETA ) ) ;
		torsion_ids.push_back( TorsionID( suite1+1, id::BB, ALPHA ) ) ;
		torsion_ids.push_back( TorsionID( suite1+1, id::BB, BETA ) ) ;
		torsion_ids.push_back( TorsionID( suite1+1, id::BB, GAMMA ) ) ;

		torsion_ids.push_back( TorsionID( suite2,   id::BB, EPSILON ) ) ;
		torsion_ids.push_back( TorsionID( suite2,   id::BB, ZETA ) ) ;
		torsion_ids.push_back( TorsionID( suite2+1, id::BB, ALPHA ) ) ;
		torsion_ids.push_back( TorsionID( suite2+1, id::BB, BETA ) ) ;
		torsion_ids.push_back( TorsionID( suite2+1, id::BB, GAMMA ) ) ;

		for ( Size n = 1; n <= torsion_ids.size(); n++ ) {
			if ( std::abs( numeric::principal_angle_degrees( pose1.torsion( torsion_ids[n] ) - pose2.torsion( torsion_ids[n]) ) ) > angle_range ) return false;
		}
		return true;

	}

	///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////
	// All this Jacobian stuff might be better in a util file.
	///////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////
	void
	RNA_LoopCloseSampler::calculate_perturbed_solutions( pose::Pose const & pose,
																											 Size const moving_suite,
																											 Size const chainbreak_suite,
																											 utility::vector1< Real > const & rbs,
																											 Real const perturbation_size,
																											 utility::vector1< utility::vector1< utility::vector1< Real > > >  & perturbed_solution_torsions ){

		using namespace pose;
		using namespace protocols::swa::rna;
		using namespace numeric::conversions;

		perturbed_solution_torsions.clear();
		Pose pose_scratch = pose; // This costs a lot, but should not happen that often.

		RNA_AnalyticLoopCloser rna_analytic_loop_closer( moving_suite, chainbreak_suite );

		utility::vector1< Real > how_much_to_perturb;
		how_much_to_perturb.push_back( degrees( perturbation_size ) ); //alpha ( euler )
		Real const cos_beta= cos( radians( rbs[2] ) );
		//	std::cout << "DEL(COSBETA) in DEGREES: " << degrees(acos(cos_beta+perturbation_size) -acos( cos_beta)) << std::endl;
		how_much_to_perturb.push_back( degrees( acos( cos_beta + perturbation_size )  - acos( cos_beta ) )  ); //beta ( euler )
		how_much_to_perturb.push_back( degrees( perturbation_size ) ); //gamma ( euler )
		how_much_to_perturb.push_back( perturbation_size  ); //x
		how_much_to_perturb.push_back( perturbation_size  ); //y
		how_much_to_perturb.push_back( perturbation_size  ); //z

		for( Size i = 1; i <= 6; i++ ){
			utility::vector1< Real > rbs_perturb;
			rbs_perturb = rbs;
			rbs_perturb[ i ] += how_much_to_perturb[ i ];
			apply_rigid_body_settings( pose_scratch, rbs_perturb, moving_suite+1, chainbreak_suite );
			rna_analytic_loop_closer.apply( pose_scratch );
			utility::vector1< utility::vector1< Real > > const torsions_for_all_solutions =
				rna_analytic_loop_closer.get_torsions_for_all_solutions();
			perturbed_solution_torsions.push_back( torsions_for_all_solutions );
		}

	}

	///////////////////////////////////////////////////////////////
	utility::vector1< Real >
	RNA_LoopCloseSampler::find_closest_solution(
																							utility::vector1< utility::vector1< Real > >  const & perturbed_torsions,
																							utility::vector1< Real > const & solution_torsions ){

		using namespace numeric::conversions;

		Size closest_set( 0 );
		Real closest_deviation( 0.0 );
		for ( Size i = 1; i <= perturbed_torsions.size(); i++ ){
			Real deviation( 0.0 );
			for ( Size j = 1; j <= solution_torsions.size(); j++ ){
				deviation += std::abs( numeric::principal_angle( radians( solution_torsions[j] - perturbed_torsions[i][j] ) ) );
			}
			if ( deviation < closest_deviation || i == 1 ){
				closest_deviation = deviation;
				closest_set = i;
			}
		}

		if( closest_set == 0 ) utility_exit_with_message( "Not enough solutions in find_closest_solution?" );

		return perturbed_torsions[ closest_set ];
	}


	///////////////////////////////////////////////////////////////
	void
	RNA_LoopCloseSampler::get_jacobian( utility::vector1< Real > const & solution_torsions,
																			utility::vector1< utility::vector1< utility::vector1< Real > > >  const & perturbed_solution_torsions,
																			Real const perturbation_size,
																			utility::vector1< utility::vector1< Real > > & J ){

		using namespace numeric::conversions;

		for ( Size i = 1; i <= 6; i++ ){
			utility::vector1< Real > closest_perturbed_torsions = find_closest_solution( perturbed_solution_torsions[ i ], solution_torsions );
			for ( Size j = 1; j <= 6; j++ ){
				J[i][j] = numeric::principal_angle( radians( closest_perturbed_torsions[j] - solution_torsions[ j ] ) )/ perturbation_size;
				//			std::cout << "Filling jacobian " << i << ' ' << j << " --> " << closest_perturbed_torsions[ j ] << ' ' << solution_torsions[ j ] << " " << perturbation_size << "  --> " << J[i][j] << std::endl;
			}
		}
	}



}
}
}
