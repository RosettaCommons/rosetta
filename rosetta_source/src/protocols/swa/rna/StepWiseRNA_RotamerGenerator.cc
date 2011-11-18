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
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.hh>
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>
//////////////////////////////////

#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
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


	////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_RotamerGenerator::StepWiseRNA_RotamerGenerator(
							 core::pose::Pose const & pose,
               Size const moving_suite,
							 PuckerState const & pucker1,
							 PuckerState const & pucker2,
							 bool const & Is_bulge,
							 Real const bin_size ):
		moving_suite_( moving_suite ),
		Is_bulge_(Is_bulge),
		pucker1_specified_( pucker1 ),
		pucker2_specified_( pucker2 ),
		sample_extra_rotamers_( true ),
		fast_( false ),
		bin_size_( bin_size ), // must be 20, 10, or 5
		bins1_( 360/bin_size_ ), //This is total bins, default is 18
		bins2_( bins1_/2 ), //This is total bins divided by 2; default is 9
		bins3_( bins1_/3  ), //This is total bins divided by 3; default is 6
		bins4_( 1 + 40/bin_size_ ), //This is the bin for chi and episilon, these two torsion angles vary from -20+mean to 20+mean
		verbose_(true)
	{
		Output_title_text("Enter StepWiseRNA_RotamerGenerator Constructor");
		initialize_syn_chi( pose );
		initialize_rotamers();
		initialize_extra_rotamer_perturbations();
		reset();
		Output_title_text("Exit StepWiseRNA_RotamerGenerator Constructor");
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
	PuckerState
	StepWiseRNA_RotamerGenerator::pucker_state(std::string const which_sugar){

		using namespace core::id;

		//A variable declared static in a function retains its state between calls to that function. Parin Feb 6, 2010
		static scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;
		Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );

		TorsionID const & torsion_id= (which_sugar=="lower") ? TorsionID( moving_suite_ , BB, 4 ) : TorsionID( moving_suite_+1 , BB, 4 );

		Real delta_value;

		for(Size n=1; n<=torsion_ids_.size(); n++){
			if(torsion_ids_[n]==torsion_id){
				delta_value=rotamer_centers_[group_rotamer_][n];

				std::cout << "  DELTA_CUTOFF angle=" << DELTA_CUTOFF << "  delta angle=" << delta_value << std::endl;
				if (delta_value <= DELTA_CUTOFF) {
					return NORTH;
				} else {
					return SOUTH;
				}
			}
		}

		utility_exit_with_message( "torsion_ids does not contain the desire delta torsion id" );
		exit (1); //Prevent compiler warning!

	}
	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::update_to_next_rotamer(){

		using namespace core::id;

		assert( has_another_rotamer() );
		rotamer_list_.clear();

		for(Size n=1; n<=torsion_ids_.size(); n++){
			Torsion_Info torsion_info;
			torsion_info.id=torsion_ids_[n];
			torsion_info.value=rotamer_centers_[group_rotamer_][n];
			rotamer_list_.push_back(torsion_info);
		}

//		rotamer_values_ =  rotamer_centers_[ group_rotamer_ ];

		if ( sample_extra_rotamers_ ) {

			for ( Size i = 1; i <= perturb_torsion_ids_.size(); i++ ){
				TorsionID const & torsion_id = perturb_torsion_ids_[ i ] ;

				// Need to figure out which torsion calue to add on to.
				for ( Size n = 1; n <= rotamer_list_.size(); n++ ){
					if ( rotamer_list_[n].id == torsion_id ){
						rotamer_list_[ n ].value += extra_rotamer_perturbations_[ subgroup_rotamer_ ][ i ];
					}
				}
			}

			subgroup_rotamer_++;
			if ( subgroup_rotamer_ > extra_rotamer_perturbations_.size() ) {
				if(verbose_) std::cout << "Done with rotamer ... " << group_rotamer_ << " out of " << rotamer_centers_.size() << std::endl;
				subgroup_rotamer_ = 1;
				group_rotamer_++;
			}


		} else {
			group_rotamer_++;
		}

	}

	////////////////////////////////////////////////////////////////////////
/*
	void
	StepWiseRNA_RotamerGenerator::initialize_puckers(
								 core::pose::Pose const & pose,
								 bool const & sample_sugar_and_base1,
								 bool const & sample_sugar_and_base2 )
	{
		pucker1_specified_ = ALL;
		pucker2_specified_ = ALL;

		if ( !sample_sugar_and_base1 ) pucker1_specified_ = Get_residue_pucker_state( pose, moving_suite_ );
		if ( !sample_sugar_and_base2 ) pucker2_specified_ = Get_residue_pucker_state( pose, moving_suite_ + 1 );
	}
*/
	////////////////////////////////////////////////////////////////////////

	//This should be call in initialize_rotamers() function
	void
	StepWiseRNA_RotamerGenerator::initialize_syn_chi(
								 core::pose::Pose const & pose )
	{
		using namespace core::scoring::rna;
		sample_syn_chi1_ = is_purine( pose.residue( moving_suite_ ) );
		sample_syn_chi2_ = is_purine( pose.residue( moving_suite_+1 ) );
	}

	////////////////////////////////////////////////////////////////////////
	//This check for repeats. Parin Feb 4, 2010
	void
	StepWiseRNA_RotamerGenerator::add_torsion_id(core::id::TorsionID const torsion_id){

		// Check that the torsion_id is not already in the TorsionID list
		for ( Size n = 1; n <= torsion_ids_.size(); n++ ){
			if ( torsion_ids_[n] == torsion_id ){
				utility_exit_with_message( "Error: torsion_id is already in the torsion_id_list!" );
			}
		}

		torsion_ids_.push_back(torsion_id);

	}


	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::print_pucker_state(PuckerState const & pucker_state, std::string const tag) const {
		std::cout << tag;
		if(pucker_state==ALL) std::cout << "ALL ";
		else if(pucker_state==NORTH) std::cout << "NORTH ";
		else std::cout << "SOUTH ";
	}

	void
	StepWiseRNA_RotamerGenerator::print_base_state(BaseState const & base_state, std::string const tag) const {
		std::cout << tag;
		if(base_state==BOTH) std::cout << "BOTH ";
		else if(base_state==ANTI) std::cout << "ANTI ";
		else std::cout << "NONE ";
	}

	void
	StepWiseRNA_RotamerGenerator::initialize_rotamers(){

		using namespace core::id;

		BaseState lower_base_state=NONE;
		BaseState upper_base_state=NONE;

		core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;

		torsion_ids_.clear();
		rotamer_centers_.clear();

		///////////////////////////////////////////////////////////////////////////////////////
		if ( !pucker1_specified_ /* sample it */ ) {
			add_torsion_id(TorsionID( moving_suite_    , BB, 4 ));  //delta1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 1 )); //chi_1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 2 )); //nu2_1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 3 )); //nu2_2

			lower_base_state = (sample_syn_chi1_) ? BOTH: ANTI;
			if(Is_bulge_) lower_base_state= NONE ;
		}

		add_torsion_id( TorsionID( moving_suite_    , BB, 5 ) ); // epsilon1
		add_torsion_id( TorsionID( moving_suite_    , BB, 6 ) ); // zeta1
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 1 ) ); // alpha2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 2 ) ); // beta2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 3 ) ); // gamma2

		if ( !pucker2_specified_ /* sample it */) {
			add_torsion_id( TorsionID( moving_suite_ + 1, BB, 4 ) );  //delta1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 1 ) ); //chi_1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 2 ) ); //nu2_1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 3 ) ); //nu2_2

			upper_base_state = (sample_syn_chi2_) ? BOTH: ANTI;
			if(Is_bulge_) upper_base_state = NONE ;
		}

		if(verbose_){
			std::cout << "moving_suite_= " << moving_suite_ << " Is_bulge_= "; Output_boolean(Is_bulge_);
			std::cout << " bin_size= " << bin_size_ << " bins1_= " << bins1_ << " bins2_= " << bins2_ << " bins3_= " << bins3_ << " bins4_= " << bins4_ << std::endl;
			std::cout << "sample_syn_chi1_= "; Output_boolean(sample_syn_chi1_); std::cout << " sample_syn_chi2_= "; Output_boolean(sample_syn_chi2_);
			std::cout << " sample_extra_rotamers_= "; Output_boolean(sample_extra_rotamers_);
			print_pucker_state(pucker1_specified_, " pucker1_specified= "); print_pucker_state(pucker2_specified_ , " pucker2_specified_= ");
			print_base_state(lower_base_state, " lower_base_state= "); print_base_state(upper_base_state, " upper_base_state= ");
			std::cout << std::endl;
		}

		Real epsilon1( 0.0 ),	zeta1( 0.0 ),	alpha2( 0.0 ),	beta2( 0.0 ),	gamma2( 0.0 );

		StepWiseRNA_Base_Sugar_RotamerOP lower_base_sugar_rotamer = new StepWiseRNA_Base_Sugar_Rotamer( lower_base_state, pucker1_specified_, rna_fitted_torsion_info, bin_size_, bins4_ );
		StepWiseRNA_Base_Sugar_RotamerOP upper_base_sugar_rotamer = new StepWiseRNA_Base_Sugar_Rotamer( upper_base_state, pucker2_specified_, rna_fitted_torsion_info, bin_size_, bins4_ );

		while(lower_base_sugar_rotamer->get_next_rotamer()){

			upper_base_sugar_rotamer->reset();

			while(upper_base_sugar_rotamer->get_next_rotamer()){
//			std::cout << "Exit upper_base_sugar_rotamer->get_next_rotamer()" << std::endl;

				std::cout << " 	delta1= " <<  F(8, 3, lower_base_sugar_rotamer->delta()) << " 	chi_1= " <<  F(8, 3, lower_base_sugar_rotamer->chi());
				std::cout << " 	nu2_1= " <<  F(8, 3, lower_base_sugar_rotamer->nu2())    << " 	nu1_1= " <<  F(8, 3, lower_base_sugar_rotamer->nu1());
				std::cout << " 	delta2= " <<  F(8, 3, upper_base_sugar_rotamer->delta()) << " 	chi_2= " <<  F(8, 3, upper_base_sugar_rotamer->chi());
				std::cout << " 	nu2_2= " <<  F(8, 3, upper_base_sugar_rotamer->nu2())    << " 	nu1_2= " <<  F(8, 3, upper_base_sugar_rotamer->nu1());
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
									epsilon1= (lower_base_sugar_rotamer->current_pucker_state() == NORTH) ? -150.17 : -98.45;

									if(z1==1) zeta1= 0.0;
									if(z1==2) zeta1= 180.0;

									if(a2==1) alpha2= 240.0;
									if(a2==2) alpha2=  0.0;
									if(a2==3) alpha2= 120.0;

									beta2 = 180.0;

									if(g2==1) gamma2= 240.0;
									if(g2==2) gamma2= 0.0;
									if(g2==3) gamma2= 120.0;

									utility::vector1 < Real >  rotamer_center;

									if ( !pucker1_specified_ ) {
										rotamer_center.push_back( lower_base_sugar_rotamer->delta());    //1
										rotamer_center.push_back( lower_base_sugar_rotamer->chi());     //2
										rotamer_center.push_back( lower_base_sugar_rotamer->nu2() );     //3
										rotamer_center.push_back( lower_base_sugar_rotamer->nu1() );     //4
									}

									rotamer_center.push_back( epsilon1 );  //5
									rotamer_center.push_back( zeta1 );     //6
									rotamer_center.push_back( alpha2 );    //7
									rotamer_center.push_back( beta2 );     //8
									rotamer_center.push_back( gamma2 );    //9

									if ( !pucker2_specified_ ) {
										rotamer_center.push_back( upper_base_sugar_rotamer->delta() );    //10
										rotamer_center.push_back( upper_base_sugar_rotamer->chi() );     //11
										rotamer_center.push_back( upper_base_sugar_rotamer->nu2() );     //12
										rotamer_center.push_back( upper_base_sugar_rotamer->nu1());     //13
									}

									rotamer_centers_.push_back( rotamer_center );
								}
							}
						}
					}
				}
			}
		}

		std::cout << "rotamer_centers_[1].size()= " << rotamer_centers_[1].size() << "  torsion_ids_.size()= " << torsion_ids_.size() << std::endl;
		if(rotamer_centers_[1].size()!=torsion_ids_.size()){
			utility_exit_with_message( "rotamer_centers_[1].size()!=torsion_ids_.size()" );
		}

	}


	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::initialize_extra_rotamer_perturbations(){

		using namespace core::id;

		///////////////////////////////////////////////////////////////////////////////////////
		perturb_torsion_ids_.clear();
		extra_rotamer_perturbations_.clear();

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

						for ( int b2_std = 1; b2_std <= bins2_; b2_std++ ) {

							perturb_beta2 = bin_size_ * static_cast<Real>( b2_std - (bins2_+1)/2 );

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

		std::cout << "extra_rotamer_perturbations_[1].size()= " << extra_rotamer_perturbations_[1].size() << "  perturb_torsion_ids_.size()= " << perturb_torsion_ids_.size() << std::endl;
		if(extra_rotamer_perturbations_[1].size()!=perturb_torsion_ids_.size()){
			utility_exit_with_message( "extra_rotamer_perturbations_[1].size()!=perturb_torsion_ids_.size()" );
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size
	StepWiseRNA_RotamerGenerator::num_rotamer_centers(){ return rotamer_centers_.size() ; }

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size const &
	StepWiseRNA_RotamerGenerator::group_rotamer(){ return group_rotamer_; }

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size const &
	StepWiseRNA_RotamerGenerator::subgroup_rotamer(){ return subgroup_rotamer_; }

	//////////////////////////////////////////////////////////////////////////////////////
	core::Size const &
	StepWiseRNA_RotamerGenerator::moving_suite(){ return moving_suite_; }

}
}
}
