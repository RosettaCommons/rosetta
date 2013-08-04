// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts = 2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_BaseSugarRotamer
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_BaseSugarRotamer.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/random/random.hh>
#include <utility/tools/make_vector1.hh>

#include <string>

using namespace core;
using core::Real;
static numeric::random::RandomGenerator RG( 2952388 );  // <- Magic number, do not change it!
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_BaseSugarRotamer" );

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_BaseSugarRotamer::StepWiseRNA_BaseSugarRotamer(
											BaseState const & base_state,
											PuckerState const & pucker_state,
											core::scoring::rna::RNA_FittedTorsionInfo const & rna_fitted_torsion_info,
											core::Size const bin_size ): //This is to determine the bin value for chi
		base_state_( base_state ),
		pucker_state_( pucker_state ),
		rna_fitted_torsion_info_( rna_fitted_torsion_info ),
		inputted_bin_size_( bin_size ), // must be 20, 10, or 5
		extra_anti_chi_( false ),
		extra_syn_chi_( false ),
		choose_random_( false ),
		rotamer_count_( 1 )
	{
		reset();

//		num_base_std_ID_act_ = (1 + 40/bin_size_)


		base_state_list_.clear(); //April 30, 2011
		if ( base_state_ == BOTH ){
			base_state_list_.push_back( ANTI );
			base_state_list_.push_back( SYN );
		} else if ( base_state_ == ANTI ){
			base_state_list_.push_back( ANTI );
		} else if ( base_state == SYN ){
			base_state_list_.push_back( SYN );
		} else if ( base_state == NONE ){
			//TR << "BLAH: base_state == NONE" << std::endl;
		} else {
			utility_exit_with_message( "Invalid base_state_ = " + ObjexxFCL::string_of( base_state_ ) );
		}

		//WARNING..DO NOT TRY TO INTRODUCE "NONE" value. pucker_state needs alway be defined even if not sampled, since EPSILON (in Rotamer Generator) depends on value of pucker_state! April 30, 2011
		pucker_state_list_.clear();
		if ( pucker_state_ == ALL ){
			pucker_state_list_.push_back( NORTH );
			pucker_state_list_.push_back( SOUTH );
		} else if ( pucker_state_ == NORTH ){
			pucker_state_list_.push_back( NORTH );
		} else if ( pucker_state_ == SOUTH ){
			pucker_state_list_.push_back( SOUTH );
		} else {
			utility_exit_with_message( "Invalid pucker_state_ = " + ObjexxFCL::string_of( pucker_state_ ) );
		}

	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_BaseSugarRotamer::~StepWiseRNA_BaseSugarRotamer()
  {}

  //////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseSugarRotamer::get_next_rotamer(){

		if ( master_rotamer_list_.size() == 0 ) initialize_master_rotamer_list();

		if ( rotamer_count_ > master_rotamer_list_.size() ) return false;

		utility::vector1< Size > & rotamer = master_rotamer_list_[ rotamer_count_ ];
		if ( choose_random_ ) rotamer = RG.random_element( master_rotamer_list_ );

		pucker_ID_   = rotamer[1];
		base_ID_     = rotamer[2];
		base_std_ID_ = rotamer[3];

		//		std::cerr << "STATE: " << pucker_ID_ << " " << base_ID_ << " " << base_std_ID_ << std::endl;
		get_next_rotamer_original(); // this is a silly hack to actually fill chi, etc.

		rotamer_count_++;
		return true;

	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseSugarRotamer::initialize_master_rotamer_list(){
		runtime_assert( rotamer_count_ == 1 );
		runtime_assert( pucker_ID_ == 1 );
		// new -- keep track of all possibilities in a master list.
		// this permits selection of random state.
		reset();
		master_rotamer_list_.clear();
		while ( get_next_rotamer_original() ){
			// The "old" values are the ones consistent with the current pose, before incrementing IDs to next rotamer.
			// this is confusing, but should be fixable if we refactor a bit.
			master_rotamer_list_.push_back( utility::tools::make_vector1( pucker_ID_old_, base_ID_old_, base_std_ID_old_ ) );
		}
		// TR <<  "Number of states in StepWiseRNA_BaseSugarRotamer = " + ObjexxFCL::string_of( master_rotamer_list_.size() );
		rotamer_count_ = 1;
		reset();
	}

  //////////////////////////////////////////////////////////////////////////
	// original behavior -- does not make out a full list of possibilities. works
	// great (and not memory intensive), but makes it harder to get random states.
  //////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseSugarRotamer::get_next_rotamer_original(){

		if ( pucker_ID_ > pucker_state_list_.size() ) return false;

		//old_ID is always consistent the stored chi_/delta_ value
		pucker_ID_old_ = pucker_ID_;
		base_ID_old_ = base_ID_;
		base_std_ID_old_ = base_std_ID_;

		PuckerState const & curr_pucker_state = current_pucker_state();

		Size base_center_ID = 99;
		if ( base_state_list_.size() != 0 ){

			if ( ( extra_anti_chi_ && ( base_state_list_[base_ID_] == ANTI ) ) || ( extra_syn_chi_ && ( base_state_list_[base_ID_] == SYN ) ) ){
	//			total_variation_ = 60; // + -30 Aug_29 to Sept 15 2010 //Mod out on May 06, 2011
	//			total_variation_ = 100; // + -50 //testing on Sept 15 2010
	//			bin_size_ = std::min(int(inputted_bin_size_), 10); // Aug_29 to Sept 15 2010 //Mod out on May 06, 2011

				total_variation_ = 120;
				bin_size_ = inputted_bin_size_;
				num_base_std_ID_ = Size( 1 + ( total_variation_/bin_size_ ) ); //-40,-20,0,20,40
			} else {
				total_variation_ = 40;	// + -20
				bin_size_ = inputted_bin_size_;
				num_base_std_ID_ = Size( 1 + ( total_variation_/bin_size_ ) ); //-20,0,20
			}

			if ( base_state_list_[base_ID_] == ANTI ){
				base_center_ID = 1;
			} else if ( base_state_list_[base_ID_] == SYN ){
				base_center_ID = 2;
			} else {
				utility_exit_with_message( "Invalid current_base_state" + ObjexxFCL::string_of( base_state_list_[base_ID_] ) );
			}
		}

		if ( curr_pucker_state == NORTH ) {

			if ( base_state_list_.size() != 0 ){
				chi_ = rna_fitted_torsion_info_.gaussian_parameter_set_chi_north()[ base_center_ID ].center + bin_size_*( base_std_ID_ - 1 ) - ( total_variation_/2 );
			}

			delta_ = rna_fitted_torsion_info_.ideal_delta_north();
			nu2_ = rna_fitted_torsion_info_.ideal_nu2_north();
			nu1_ = rna_fitted_torsion_info_.ideal_nu1_north();
		}	else if ( curr_pucker_state == SOUTH ) {

			if ( base_state_list_.size() != 0 ){
				chi_ = rna_fitted_torsion_info_.gaussian_parameter_set_chi_south()[ base_center_ID ].center + bin_size_*( base_std_ID_ - 1 ) - ( total_variation_/2 );
			}

			delta_ = rna_fitted_torsion_info_.ideal_delta_south(); //default
			nu2_ = rna_fitted_torsion_info_.ideal_nu2_south(); //default
			nu1_ = rna_fitted_torsion_info_.ideal_nu1_south(); //default

//			delta_ = rna_fitted_torsion_info_.gaussian_parameter_set_delta_south()[ 1 ].center;
//			nu2_ = rna_fitted_torsion_info_.gaussian_parameter_set_nu2_south()[ 1 ].center;
//			nu1_ = rna_fitted_torsion_info_.gaussian_parameter_set_nu1_south()[ 1 ].center;
//			nu2_ = -38.9; //2eew res A_19 torsion value
//			nu1_ = 161.8; //2eew res A_19 torsion value
//	    delta_ = 146.8; //1q9a res G_9 torsion value
//			nu2_ = -41.6; //1q9a res G_9 torsion value
//			nu1_ = 167.2; //1q9a res G_9 torsion value

		} else {
			utility_exit_with_message( "Invalid current_pucker_state!" + ObjexxFCL::string_of( curr_pucker_state ) );
		}

		base_std_ID_++;

		if ( base_state_list_.size() != 0 ){
			if ( base_std_ID_ > num_base_std_ID_ ){
				base_std_ID_ = 1;
				base_ID_++;
			}
		}

		if ( base_ID_ > base_state_list_.size() ){
			base_ID_ = 1;
			pucker_ID_++;
		}

		return true;
	}

  //////////////////////////////////////////////////////////////////////////
	PuckerState const &
	StepWiseRNA_BaseSugarRotamer::current_pucker_state() const {

//		TR << "pucker_ID_old_ = " << pucker_ID_old_ << std::endl;
//		TR << "pucker_state_list_.size() = " << pucker_state_list_.size() << std::endl;

		PuckerState const & pucker_state = pucker_state_list_[pucker_ID_old_];

		if ( pucker_state != NORTH && pucker_state != SOUTH ) utility_exit_with_message( "pucker_state should equal NORTH or SOUTH!" );

		return pucker_state;
	}

  //////////////////////////////////////////////////////////////////////////

	std::string const
	StepWiseRNA_BaseSugarRotamer::current_base_state() const {

		using namespace ObjexxFCL;

		Real const principal_chi = numeric::principal_angle_degrees( chi_ );
		std::string base_state_string;

		if ( base_state_list_.size() == 0 ){
			base_state_string = "ZZ";
		} else {

			if ( principal_chi > 0.0 ){
				if ( base_state_list_[base_ID_old_] != ANTI ) utility_exit_with_message( "principal_chi > 0 but base_state_list_( base_ID_old_ ) != ANTI" );
				base_state_string = "A"; //anti
			} else {
				if ( base_state_list_[base_ID_old_] != SYN ) utility_exit_with_message( "principal_chi <= 0 but base_state_list_( base_ID_old_ ) != SYN" );
				base_state_string = "S"; //syn
			}

			base_state_string += string_of( base_std_ID_old_ );
		}

		return base_state_string;
	}

	//////////////////////////////////////////////////////////////////////////
	std::string const
	StepWiseRNA_BaseSugarRotamer::current_tag() const {

		std::string tag = current_base_state();

		if ( current_pucker_state() == NORTH ){
			tag += "N";
		} else if ( current_pucker_state() == SOUTH ){
			tag += "S";
		} else {
			utility_exit_with_message( "Invalid current_pucker_state!" );
		}

		return tag;

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseSugarRotamer::reset() {

		rotamer_count_ = 1;

		pucker_ID_ = 1;
		base_ID_ = 1;
		base_std_ID_ = 1;

		pucker_ID_old_ = pucker_ID_;
		base_ID_old_ = base_ID_;
		base_std_ID_old_ = base_std_ID_;

		chi_ = 0.0;
		delta_ = 0.0;
		nu2_ = 0.0;
		nu1_ = 0.0;

	}

  //////////////////////////////////////////////////////////////////////////



}
}
}
