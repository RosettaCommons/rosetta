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
#include <protocols/swa/rna/StepWiseRNA_BaseSugarRotamer.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseSugarRotamer.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
//////////////////////////////////

#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

#include <string>

using namespace core;
using core::Real;
using ObjexxFCL::fmt::F;

static numeric::random::RandomGenerator RG ( 26699940 ); // <- Magic number, do not change it!

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_rotamer_generator" );

namespace protocols {
namespace swa {
namespace rna {


	////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_RotamerGenerator::StepWiseRNA_RotamerGenerator(
               Size const moving_suite,
							bool const sample_lower_sugar_and_base,
							bool const sample_upper_sugar_and_base,
							PuckerState const pucker1,
							PuckerState const pucker2):
		moving_suite_( moving_suite ),
		sample_lower_sugar_and_base_(sample_lower_sugar_and_base),
		sample_upper_sugar_and_base_(sample_upper_sugar_and_base),
		pucker1_specified_( pucker1 ),
		pucker2_specified_( pucker2 ),
		sample_extra_rotamers_( true ),
		fast_( false ),
		bin_size_( 20 ), // must be 20, 10, or 5
		verbose_(true),
		sample_chi_torsion_(true), //Oct 2, 2010
		include_syn_chi_(true),
		extra_epsilon_(false), //Aug 30, 2010
		extra_beta_(false), //Aug 30, 2010
		extra_anti_chi_(false), //Aug 30, 2010
		extra_syn_chi_(false), //Aug 30, 2010
		exclude_alpha_beta_gamma_sampling_(false),
		allow_syn_pyrimidine_(false), //Nov 14, 2010
		lower_base_state_( NONE ), //Will be initialized by initialize_sample_base_states, April 29, 2011
		upper_base_state_( NONE ),  //Will be initialized by initialize_sample_base_states, April 29, 2011
		choose_random_( false )
	{

		//These vectors should be empty to begin with, but not harm to ensure this.
		force_syn_chi_res_list_.clear();
	}

	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGenerator::~StepWiseRNA_RotamerGenerator(){}
	////////////////////////////////////////////////////////////////////////

	//This function is called only after all the parameter setting functions (set_bin_size, set_include_syn_chi, set_sample_extra_rotamers, set_fast and etc..) had been called!...
	void
	StepWiseRNA_RotamerGenerator::initialize_rotamer_generator(core::pose::Pose const & pose){
		Output_title_text("Enter StepWiseRNA_RotamerGenerator::initialize_rotamer_generator", TR );

		bins1_= 360/bin_size_;     //This is total bins, default is 18
		bins2_= bins1_/2;          //This is total bins divided by 2; default is 9
		bins3_= bins1_/3; 				 //This is total bins divided by 3; default is 6

		//beta_bins_ = (extra_beta_) ?  (bins1_) : ((bins1_/2)+2) ; //Before Mar 14, 2011
		beta_bins_ = (extra_beta_) ?  (bins1_) : ( (200/bin_size_)+1 ) ; //Mar 14, 2011
		eps_bins_= (extra_epsilon_) ? (1+120/bin_size_) : (1 +40/bin_size_) ;

		initialize_sample_base_states( pose );
		initialize_rotamers();
		initialize_extra_rotamer_perturbations();
		reset();
		Output_title_text("Exit StepWiseRNA_RotamerGenerator::initialize_rotamer_generator", TR );
	}

	////////////////////////////////////////////////////////////////////////
	utility::vector1< core::id::TorsionID > const &
	StepWiseRNA_RotamerGenerator::torsion_ids() const {
		return torsion_ids_;
	}

	////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_RotamerGenerator::has_another_rotamer() const{

		if ( choose_random_ ) return true; // can always get another random rotamer!

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
			if (torsion_ids_[n]==torsion_id){
				delta_value=rotamer_centers_[group_rotamer_][n];

				TR << "  DELTA_CUTOFF angle=" << DELTA_CUTOFF << "  delta angle=" << delta_value << std::endl;

				if ((delta_value>1.0 && delta_value<179.00)==false){
					utility_exit_with_message( "delta angle out of range!" );
				}

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

		if ( choose_random_ ) update_to_random_rotamer();
		else runtime_assert( has_another_rotamer() ); // update to next rotamer occurs below...

		rotamer_list_.clear();

		for (Size n = 1; n <= torsion_ids_.size(); n++){
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

		}

		//Added in Aug 30 to ensure that dihedral torsion angle is in [-180:180] range.
		for(Size n=1; n<=rotamer_list_.size(); n++){
			rotamer_list_[ n ].value=numeric::principal_angle_degrees(rotamer_list_[ n ].value);
		}

		///////////////////////////////////////////////////////////////////////////////
		// time to actually update the rotamer to get ready for the next request...
		if ( !choose_random_ ){
			if ( sample_extra_rotamers_ ) {
				subgroup_rotamer_++;
				if ( subgroup_rotamer_ > extra_rotamer_perturbations_.size() ) {
					if (verbose_) TR << "Done with rotamer ... " << group_rotamer_ << " out of " << rotamer_centers_.size() << std::endl;
					subgroup_rotamer_ = 1;
					group_rotamer_++;
				}
			} else {
				group_rotamer_++;
			}
		}

	}

	////////////////////////////////////////////////////////////////////////
	// slightly "off-label" use of RotamerGenerator.
	void
	StepWiseRNA_RotamerGenerator::update_to_random_rotamer(){
		runtime_assert( choose_random_ );
		group_rotamer_    = RG.random_range( 1, rotamer_centers_.size() );
		subgroup_rotamer_ = RG.random_range( 1, extra_rotamer_perturbations_.size() );
	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator::initialize_sample_base_states(core::pose::Pose const & pose ){
		using namespace core::scoring::rna;

		/* BEFORE APRIL 29, 2011
		if (include_syn_chi_){
			if (allow_syn_pyrimidine_){
				sample_syn_chi1_= true;
				sample_syn_chi2_= true;
			}else{
				sample_syn_chi1_ = is_purine( pose.residue( moving_suite_ ) );
				sample_syn_chi2_ = is_purine( pose.residue( moving_suite_+1 ) );
			}
		}else{
			sample_syn_chi1_= false;
			sample_syn_chi2_= false;
		}
		*/


		lower_base_state_=ANTI;
		upper_base_state_=ANTI;

		if (include_syn_chi_==true){
			if ( is_purine( pose.residue( moving_suite_   ) ) || allow_syn_pyrimidine_)	lower_base_state_=BOTH;
			if ( is_purine( pose.residue( moving_suite_+1 ) ) || allow_syn_pyrimidine_)	upper_base_state_=BOTH;
		}

		if ( force_syn_chi_res_list_.has_value(moving_suite_  ) ){ //lower base
			if ( (is_purine( pose.residue( moving_suite_   ) )==false) && (allow_syn_pyrimidine_==false) ){
				utility_exit_with_message("forcing working_res= "+ ObjexxFCL::string_of(moving_suite_  ) +" to be syn chi but residue is pyrimidine and allow_syn_pyrimidine_==false!");
			}
			lower_base_state_=SYN;
		}

		if ( force_syn_chi_res_list_.has_value(moving_suite_+1) ){ //upper base
			if ( (is_purine( pose.residue( moving_suite_+1 ) )==false) && (allow_syn_pyrimidine_==false) ){
				utility_exit_with_message("forcing working_res= "+ ObjexxFCL::string_of(moving_suite_+1) +" to be syn chi but residue is pyrimidine and allow_syn_pyrimidine_==false!");
			}
			upper_base_state_=SYN;
		}

		if ( sample_chi_torsion_==false){
			lower_base_state_=NONE;
			upper_base_state_=NONE;
		}


		if ( sample_lower_sugar_and_base_==false  /* DON'T sample sugar+base */) lower_base_state_=NONE;
		if ( sample_upper_sugar_and_base_==false  /* DON'T sample sugar+base */) upper_base_state_=NONE;

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
	StepWiseRNA_RotamerGenerator::initialize_rotamers(){

		using namespace core::id;



		core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;

		torsion_ids_.clear();
		rotamer_centers_.clear();

		///////////////////////////////////////////////////////////////////////////////////////
		if ( sample_lower_sugar_and_base_ /* sample it */ ){
			add_torsion_id(TorsionID( moving_suite_    , BB, 4 ));  //delta1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 2 )); //nu2_1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 3 )); //nu1_1
		}

		if (lower_base_state_!=NONE){
			add_torsion_id(TorsionID( moving_suite_    , CHI, 1 )); //chi_1
		}

		add_torsion_id( TorsionID( moving_suite_    , BB, 5 ) ); // epsilon1
		add_torsion_id( TorsionID( moving_suite_    , BB, 6 ) ); // zeta1
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 1 ) ); // alpha2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 2 ) ); // beta2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 3 ) ); // gamma2

		if ( sample_upper_sugar_and_base_ /* sample it */){
			add_torsion_id( TorsionID( moving_suite_ + 1, BB, 4 ) );  //delta2
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 2 ) ); //nu2_2 (this is the one that include the delta bond atom)
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 3 ) ); //nu1_2 (this is the one that include the base nitrogen atom), actually nearer to chi
		}

		if (upper_base_state_!=NONE){
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 1 ) ); //chi_2
		}

		if (verbose_){
			TR << "moving_suite_= " << moving_suite_ << " bin_size= " << bin_size_ << " bins1_= " << bins1_ << " bins2_= " << bins2_ << " bins3_= " << bins3_ << std::endl;
			Output_boolean("extra_beta_: ", extra_beta_, TR );
			Output_boolean(" extra_epsilon_: ", extra_epsilon_, TR);
			Output_boolean(" extra_anti_chi_: ", extra_anti_chi_, TR );
			Output_boolean(" extra_syn_chi_: ", extra_syn_chi_, TR);
			TR << " beta_bins_= " << beta_bins_ << " eps_bins_ = " << eps_bins_ << std::endl;
			Output_boolean("sample_chi_torsion_= ", sample_chi_torsion_, TR );
			Output_boolean(" sample_extra_rotamers_= ", sample_extra_rotamers_, TR );
			Output_boolean(" exclude_alpha_beta_gamma_sampling_: ", exclude_alpha_beta_gamma_sampling_, TR );
			TR << std::endl;
			Output_boolean("sample_lower_sugar_and_base_= " , sample_lower_sugar_and_base_, TR );
			print_ribose_pucker_state(" pucker1_specified= ", pucker1_specified_);
			print_base_state(" lower_base_state= ", lower_base_state_);
			Output_boolean(" sample_upper_sugar_and_base_= " , sample_upper_sugar_and_base_, TR );
			print_ribose_pucker_state(" pucker2_specified_= ", pucker2_specified_);
			print_base_state(" upper_base_state= ", upper_base_state_ );
			TR << std::endl;
		}

		if (pucker1_specified_==ALL && sample_lower_sugar_and_base_==false){
			utility_exit_with_message("pucker1_specified==ALL but sample_lower_sugar_and_base_==false");
		}

		if (pucker2_specified_==ALL && sample_upper_sugar_and_base_==false){
			utility_exit_with_message("pucker2_specified==ALL but sample_upper_sugar_and_base_==false");
		}

		Real epsilon1( 0.0 ),	zeta1( 0.0 ),	alpha2( 0.0 ),	beta2( 0.0 ),	gamma2( 0.0 );

		StepWiseRNA_BaseSugarRotamerOP lower_base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( lower_base_state_, pucker1_specified_, rna_fitted_torsion_info, bin_size_);
		StepWiseRNA_BaseSugarRotamerOP upper_base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( upper_base_state_, pucker2_specified_, rna_fitted_torsion_info, bin_size_);

		lower_base_sugar_rotamer->set_extra_anti_chi(extra_anti_chi_);
		upper_base_sugar_rotamer->set_extra_anti_chi(extra_anti_chi_);

		lower_base_sugar_rotamer->set_extra_syn_chi(extra_syn_chi_);
		upper_base_sugar_rotamer->set_extra_syn_chi(extra_syn_chi_);


		while (lower_base_sugar_rotamer->get_next_rotamer()){

			upper_base_sugar_rotamer->reset();

			while (upper_base_sugar_rotamer->get_next_rotamer()){
//			TR << "Exit upper_base_sugar_rotamer->get_next_rotamer()" << std::endl;

				TR << " 	delta1= " <<  F(8, 3, lower_base_sugar_rotamer->delta()) << " 	chi_1= " <<  F(8, 3, lower_base_sugar_rotamer->chi());
				TR << " 	nu2_1= " <<  F(8, 3, lower_base_sugar_rotamer->nu2())    << " 	nu1_1= " <<  F(8, 3, lower_base_sugar_rotamer->nu1());
				TR << " 	delta2= " <<  F(8, 3, upper_base_sugar_rotamer->delta()) << " 	chi_2= " <<  F(8, 3, upper_base_sugar_rotamer->chi());
				TR << " 	nu2_2= " <<  F(8, 3, upper_base_sugar_rotamer->nu2())    << " 	nu1_2= " <<  F(8, 3, upper_base_sugar_rotamer->nu1());
				TR << std::endl;

				//epsilon depends on delta of same residue...
				for (int e1 = 1; e1 <= 1; e1++ ) {
					for (int z1 = 1; z1 <= 2; z1++ ) {
						// zeta depends on alpha of next residue... note that in original nested loop
						//   we looped on alpha2 first, then zeta1, since they are now sampled completely,
						//   in 20 degree bins.
						for (int a2 = 1; a2 <= 3; a2++ ) {
							for (int b2 = 1; b2 <= 1; b2++ ) {
								for (int g2 = 1; g2 <= 3; g2++ ) {

									if (exclude_alpha_beta_gamma_sampling_==true){
										if (a2!=1 || b2!=1 || g2!=1) continue;
									}

									// -150.17 and -98.45 MIGHT but match the value in RNA_FittedTorsionInfo? If that is the case, then value came from Parin calculatation which exclude A-form
									//The mean is shifted in the extra_epsilon mode since the distribution is for both North and South are uneven.
									//For North the tail is longer in the right side. For South the tail is longer in the left side
									if (extra_epsilon_){ //+/- 60 degrees, added on Aug 30 2010 to fix BB torsion between Bulged G and the nucleotide 3' of it. Epsilon=-167.7 Delta=SOUTH
										epsilon1= (lower_base_sugar_rotamer->current_pucker_state() == NORTH) ? -130.17 : -118.45;
									}else{ // +/- 20 degrees standard mode
										epsilon1= (lower_base_sugar_rotamer->current_pucker_state() == NORTH) ? -150.17 : -98.45;
									}

									if (z1==1) zeta1= 0.0;
									if (z1==2) zeta1= 180.0;

									if (a2==1) alpha2= 240.0;
									if (a2==2) alpha2=  0.0;
									if (a2==3) alpha2= 120.0;

									beta2 = 180.0;

									if (g2==1) gamma2= 240.0;
									if (g2==2) gamma2= 0.0;
									if (g2==3) gamma2= 120.0;

									utility::vector1 < Real >  rotamer_center;

									if ( sample_lower_sugar_and_base_ /* sample it */) {
										rotamer_center.push_back( lower_base_sugar_rotamer->delta());    //1
										rotamer_center.push_back( lower_base_sugar_rotamer->nu2() );     //3
										rotamer_center.push_back( lower_base_sugar_rotamer->nu1() );     //4
									}

									if (lower_base_state_!=NONE){
										rotamer_center.push_back( lower_base_sugar_rotamer->chi());     //2
									}

									rotamer_center.push_back( epsilon1 );  //5
									rotamer_center.push_back( zeta1 );     //6
									rotamer_center.push_back( alpha2 );    //7
									rotamer_center.push_back( beta2 );     //8
									rotamer_center.push_back( gamma2 );    //9

									if ( sample_upper_sugar_and_base_ /* sample it */) {
										rotamer_center.push_back( upper_base_sugar_rotamer->delta() );    //10
										rotamer_center.push_back( upper_base_sugar_rotamer->nu2() );     //12
										rotamer_center.push_back( upper_base_sugar_rotamer->nu1());     //13
									}

									if (upper_base_state_!=NONE){
										rotamer_center.push_back( upper_base_sugar_rotamer->chi() );     //11
									}

									rotamer_centers_.push_back( rotamer_center );
								}
							}
						}
					}
				}
			}
		}

		TR << "rotamer_centers_[1].size()= " << rotamer_centers_[1].size() << "  torsion_ids_.size()= " << torsion_ids_.size() << std::endl;
		if (rotamer_centers_[1].size()!=torsion_ids_.size()){
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

		for ( int e1_std = 1; e1_std <= static_cast<int>(eps_bins_); e1_std++ ) {

			if (extra_epsilon_){
				perturb_epsilon1 = bin_size_*(e1_std-1) - 60.0;
			}else{
				perturb_epsilon1 = bin_size_*(e1_std-1) - 20.0;
			}


			for ( int z1_std = 1; z1_std <= bins2_; z1_std++ ) {

					perturb_zeta1 = bin_size_ * z1_std;

					for ( int a2_std = 1; a2_std <= bins3_; a2_std++ ) {

						perturb_alpha2 = bin_size_ * a2_std;

						for ( int b2_std = 1; b2_std <= beta_bins_; b2_std++ ) {

							if (extra_beta_){
								perturb_beta2 = bin_size_ * static_cast<Real>( b2_std - 1 ); //perturb: (0, +340) actual: (180,520) --> (-180, -160,-140,-120,-100,-80,-60.....140, 160 ) :18 bins
							}else{
								perturb_beta2 = bin_size_ * static_cast<Real>( b2_std - (beta_bins_+1)/2 );	// beta_bins=11 perturb: (-100,100) actual: (80, -80) --> (80, 100, 120, 140, .....-140, -120, -100, -80) : 9 bins

//Before Oct 12, 2010 beta_bins=9 perturb: (-80,80) actual: (100, -100) --> (100, 120, 140, .....-140, -120, -100) : 9 bins
							}

							for ( int g2_std = 1; g2_std <= bins3_; g2_std++ ) {

								perturb_gamma2 = bin_size_ * g2_std;

								if (exclude_alpha_beta_gamma_sampling_==true){
									if (a2_std!=1 || b2_std!=1 || g2_std!=1) continue;
								}


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

		TR << "extra_rotamer_perturbations_[1].size()= " << extra_rotamer_perturbations_[1].size() << "  perturb_torsion_ids_.size()= " << perturb_torsion_ids_.size() << std::endl;
		if (extra_rotamer_perturbations_[1].size()!=perturb_torsion_ids_.size()){
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
	/////////////////////////////////////////////////////////////////////////

}
}
}
