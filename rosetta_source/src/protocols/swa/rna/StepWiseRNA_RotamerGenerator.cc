// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_RotamerGenerator
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh>

#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using core::Real;
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_rotamer_generator" );

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_RotamerGenerator::StepWiseRNA_RotamerGenerator( Size const moving_suite, PuckerState const & pucker1, PuckerState const & pucker2 ):
		moving_suite_( moving_suite ),
		pucker1_specified_( pucker1 ),
		pucker2_specified_( pucker2 ),
		sample_extra_rotamers_( true ),
		fast_( false ),
		bin_size_( 20 ), // must be 20, 10, or 5
		bins1_( 360/bin_size_ ), //This is total bins, default is 18
		bins2_( bins1_/2 ), //This is total bins divided by 2; default is 9
		bins3_( bins1_/3  ), //This is total bins divided by 3; default is 6
		bins4_( 1 + 40/bin_size_ ) //This is the bin for chi and episilon, these two torsion angles vary from -20+mean to 20+mean
	{
		initialize_rotamers();
		initialize_extra_rotamer_perturbations();
		reset();
	}

	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGenerator::~StepWiseRNA_RotamerGenerator(){}


	////////////////////////////////////////////////////////////////////////
	utility::vector1< core::id::TorsionID > const &
	StepWiseRNA_RotamerGenerator::torsion_ids() const {
		return torsion_ids_;
	}

	////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_RotamerGenerator::has_another_rotamer() const{
		return ( group_rotamer_ <= rotamer_centers_.size() &&
						 (!fast_ || group_rotamer_ <= 2 ) /*premature stop*/ );
	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::reset() {
		group_rotamer_ = 1;
		subgroup_rotamer_ = 1;
	}

	////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Real > const &
	StepWiseRNA_RotamerGenerator::get_next_rotamer(){

		using namespace core::id;

		assert( has_another_rotamer() );

		rotamer_values_ =  rotamer_centers_[ group_rotamer_ ];

		if ( sample_extra_rotamers_ ) {

			for ( Size i = 1; i <= perturb_torsion_ids_.size(); i++ ){
				TorsionID const & torsion_id = perturb_torsion_ids_[ i ] ;

				// Need to figure out which torsion calue to add on to.
				for ( Size n = 1; n <= torsion_ids_.size(); n++ ){
					if ( torsion_ids_[n] == torsion_id ){
						rotamer_values_[ n ] += extra_rotamer_perturbations_[ subgroup_rotamer_ ][ i ];
					}
				}

			}

			subgroup_rotamer_++;
			if ( subgroup_rotamer_ > extra_rotamer_perturbations_.size() ) {
				std::cout << "Done with rotamer ... " << group_rotamer_ << " out of " << rotamer_centers_.size() << std::endl;
				subgroup_rotamer_ = 1;
				group_rotamer_++;
			}


		} else {
			group_rotamer_++;
		}

		return rotamer_values_;

	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::initialize_rotamers(){

		using namespace core::id;

		///////////////////////////////////////////////////////////////////////////////////////
		torsion_ids_.clear();
		if ( !pucker1_specified_ /* sample it */ ) {
			torsion_ids_.push_back( TorsionID( moving_suite_    , BB, 4 ) ); //delta1
			torsion_ids_.push_back( TorsionID( moving_suite_    , CHI, 1 ) ); // chi1
			torsion_ids_.push_back( TorsionID( moving_suite_    , CHI, 2 ) ); // nu2 [1]
			torsion_ids_.push_back( TorsionID( moving_suite_    , CHI, 3 ) ); // nu1 [1]
		}

		torsion_ids_.push_back( TorsionID( moving_suite_    , BB, 5 ) ); // epsilon1
		torsion_ids_.push_back( TorsionID( moving_suite_    , BB, 6 ) ); // zeta1
		torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 1 ) ); // alpha2
		torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 2 ) ); // beta2
		torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 3 ) ); // gamma2

		if ( !pucker2_specified_ /* sample it */) {
			torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 4 ) ); //delta2
			torsion_ids_.push_back( TorsionID( moving_suite_ + 1, CHI, 1 ) ); // chi2
			torsion_ids_.push_back( TorsionID( moving_suite_ + 1, CHI, 2 ) ); // nu2 [2]
			torsion_ids_.push_back( TorsionID( moving_suite_ + 1, CHI, 3 ) ); // nu1 [2]
		}

		rotamer_centers_.clear();

		core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;

		Real delta1( 0.0 ),	chi_1( 0.0 ),	nu2_1( 0.0 ),	nu1_1( 0.0 ),	epsilon1( 0.0 ),	zeta1( 0.0 ),	alpha2( 0.0 ),	beta2( 0.0 ),	gamma2( 0.0 ),	delta2( 0.0 ),	chi_2( 0.0 ),	nu2_2( 0.0 ),	nu1_2( 0.0 );

		for (int d1 = 1; d1 <= 2; d1++ ) {
			for (int c1 = 1; c1 <= bins4_; c1++ ) {

				if ( pucker1_specified_ > 0 /*pucker is specified*/	 &&
						 ( (pucker1_specified_ != d1) /* sample the right pucker */ ||
							 (c1 > 1) /*sample one dummy chi angle*/  ) )  continue;

				if (d1 == 1) {
					delta1 = rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center;
					chi_1 = rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center + bin_size_*(c1-1) - 20;
					nu2_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center;
					nu1_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center;
				}	else {
					delta1 = rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center;
					chi_1 = rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center + bin_size_*(c1-1) - 20;
					nu2_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center;
					nu1_1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center;
				}

				for (int d2 = 1; d2 <= 2; d2++ ) {
					for (int c2 = 1; c2 <= bins4_; c2++ ) {

						//			std::cout << " d2= " << d2;
						//			std::cout << " c2= " << c2;
						//			std::cout << " (c2-1)-20)= " << bin_size_*(c2-1)-20;
						//			std::cout << " Size((c2-1)-20)))= " << ((int)(bin_size_*(c2-1)-20)) << std::endl;

						if ( pucker2_specified_ > 0 /*pucker is specified*/	 &&
								 ( (pucker2_specified_ != d2) /* sample the right pucker */ ||
									 (c2 > 1) /*sample one dummy chi angle*/ )   ) continue;

						if (d2 == 1) {
							delta2 = rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center;
							chi_2 = rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center + bin_size_*(c2-1) - 20;
							nu2_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center;
							nu1_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center;
						}	else {
							delta2 = rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center;
							chi_2 = rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center + bin_size_*(c2-1) - 20;
							nu2_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center;
							nu1_2 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center;
						}

						std::cout << " 	delta1= " <<  F(8, 3, delta1);
						std::cout << " 	chi_1= " <<  F(8, 3, chi_1);
						std::cout << " 	nu2_1= " <<  F(8, 3, nu2_1);
						std::cout << " 	nu1_1= " <<  F(8, 3, nu1_1);
						std::cout << " 	delta2= " <<  F(8, 3, delta2);
						std::cout << " 	chi_2= " <<  F(8, 3, chi_2);
						std::cout << " 	nu2_2= " <<  F(8, 3, nu2_2);
						std::cout << " 	nu1_2= " <<  F(8, 3, nu1_2);
						std::cout << std::endl;

						//epsilon depends on delta of same residue...
						for (int e1 = 1; e1 <= 1; e1++ ) {
								for (int z1 = 1; z1 <= 2; z1++ ) {
									// zeta depends on alpha of next residue... note that in original nested loop
									//   we looped on alpha2 first, then zeta1, since they are now sampled completely,
									//   in 20 degree bins.
									for (int a2 = 1; a2 <= 3; a2++ ) {
										for (int b2 = 1; b2 <= 1; b2++ ) {
											for (int g2 = 1; g2 <= 3; g2++ ) {

											// Why do these deviate from RNA_FittedTorsionInfo?
											if (d1 == 1) epsilon1 = -150.17;
											else         epsilon1 = -98.45;

											if(a2==1) alpha2= 240.0;
											if(a2==2) alpha2=  0.0;
											if(a2==3) alpha2= 120.0;

											if(z1==1) zeta1= 0.0;
											if(z1==2) zeta1= 180.0;

											beta2 = 0.0;

											if(g2==1) gamma2= 240.0;
											if(g2==2) gamma2= 0.0;
											if(g2==3) gamma2= 120.0;

											utility::vector1 < Real >  rotamer_center;

											if ( !pucker1_specified_ ) {
												rotamer_center.push_back( delta1 );    //1
												rotamer_center.push_back( chi_1 );     //2
												rotamer_center.push_back( nu2_1 );     //3
												rotamer_center.push_back( nu1_1 );     //4
											}

											rotamer_center.push_back( epsilon1 );  //5
											rotamer_center.push_back( zeta1 );     //6
											rotamer_center.push_back( alpha2 );    //7
											rotamer_center.push_back( beta2 );     //8
											rotamer_center.push_back( gamma2 );    //9

											if ( !pucker2_specified_ ) {
												rotamer_center.push_back( delta2 );    //10
												rotamer_center.push_back( chi_2 );     //11
												rotamer_center.push_back( nu2_2 );     //12
												rotamer_center.push_back( nu1_2 );     //13
											}

											rotamer_centers_.push_back( rotamer_center );
										}
									}
								}
							}
						}
					}
				}
			}
		}


	}


	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::initialize_extra_rotamer_perturbations(){

		using namespace core::id;

		///////////////////////////////////////////////////////////////////////////////////////
		perturb_torsion_ids_.clear();

		perturb_torsion_ids_.push_back( TorsionID( moving_suite_    , BB, 5 ) ); // epsilon1
		perturb_torsion_ids_.push_back( TorsionID( moving_suite_    , BB, 6 ) ); // zeta1
		perturb_torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 1 ) ); // alpha2
		perturb_torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 2 ) ); // beta2
		perturb_torsion_ids_.push_back( TorsionID( moving_suite_ + 1, BB, 3 ) ); // gamma2

		Real perturb_epsilon1( 0.0 ),	perturb_zeta1( 0.0 ),	perturb_alpha2( 0.0 ),	perturb_beta2( 0.0 ),	perturb_gamma2( 0.0 );

		for ( int e1_std = 1; e1_std <= static_cast<int>(bins4_); e1_std++ ) {

			perturb_epsilon1 = bin_size_*(e1_std-1) - 20.0;

			for ( int z1_std = 1; z1_std <= bins2_; z1_std++ ) {

					perturb_zeta1 = bin_size_ * z1_std;

					for ( int a2_std = 1; a2_std <= bins3_; a2_std++ ) {

						perturb_alpha2 = bin_size_ * a2_std;

						for ( int b2_std = 1; b2_std <= bins1_; b2_std++ ) {

							perturb_beta2 = bin_size_ * b2_std;

							for ( int g2_std = 1; g2_std <= bins3_; g2_std++ ) {

								perturb_gamma2 = bin_size_ * g2_std;

								utility::vector1 < Real >  rotamer_perturbation;

								rotamer_perturbation.push_back( perturb_epsilon1 );
								rotamer_perturbation.push_back( perturb_zeta1 );
								rotamer_perturbation.push_back( perturb_alpha2 );
								rotamer_perturbation.push_back( perturb_beta2 );
								rotamer_perturbation.push_back( perturb_gamma2 );

								extra_rotamer_perturbations_.push_back( rotamer_perturbation );
						}
					}
				}
			}
		}
	}

	core::Size const &
	StepWiseRNA_RotamerGenerator::group_rotamer(){ return group_rotamer_; }
	core::Size const &
	StepWiseRNA_RotamerGenerator::subgroup_rotamer(){ return subgroup_rotamer_; }

}
}
}
