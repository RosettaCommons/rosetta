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
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGeneratorWrapper.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
//////////////////////////////////

#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using namespace core::chemical::rna;
using core::Real;
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.StepWiseRNA_RotamerGeneratorWrapper" );

//////////////////////////////////////////////////////////////////////////////////////////
// This chains a bunch of rotamer generators together -- particularly useful in cases
// where there are more than one nucleotide suite is to be sampled.
//
//  In the case of multiple nucleotide suites, the first suite supplied in moving_suite_list
//  is assumed to be most distal in sequence from the fixed point. [see setting of
//  Is_prepend below.]
//
//////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGeneratorWrapper::StepWiseRNA_RotamerGeneratorWrapper(
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & moving_suite_list,
								bool const & sample_sugar_and_base1,
								bool const & sample_sugar_and_base2 ):
		pose_( pose ),
		moving_suite_list_( moving_suite_list ),
		sample_sugar_and_base1_( sample_sugar_and_base1 ), // sugar/base of 5' most nucleoside
		sample_sugar_and_base2_( sample_sugar_and_base2 ), // sugar/base of 3' most nucleoside
		sample_extra_rotamers_( true ),
		fast_( false ),
		verbose_( false ),
		sample_chi_torsion_( true ), //Oct 2, 2010
		include_syn_chi_( true ),
		bin_size_( 20 ),
		extra_epsilon_( false ), //Aug 30, 2010
		extra_beta_( false ), //Aug 30, 2010
		extra_anti_chi_( false ), //Aug 30, 2010
		extra_syn_chi_( false ), //Aug 30, 2010
		exclude_alpha_beta_gamma_sampling_( false ),
		allow_syn_pyrimidine_( false ),
		choose_random_( false ),
		rotamer_generator_list_( moving_suite_list_.size(), NULL )
	{

		/////////////////////////Check that moving_suite_list_ is correctly order/////////////////////////
		bool const can_prepend = check_can_prepend( moving_suite_list_ ); //[12,13]
		bool const can_append  = check_can_append( moving_suite_list_ ); //[13,12]

		if ( !can_prepend && !can_append ){
			Output_seq_num_list( "working_moving_suite_list_:", moving_suite_list_, TR.Debug );
			utility_exit_with_message( "Cannot prepend or append residue in moving_suite_list_" );
		}

		if ( moving_suite_list_.size() > 1 ){
			if ( can_prepend && can_append ){
				Output_seq_num_list( "working_moving_suite_list_:", moving_suite_list_, TR.Debug );
				utility_exit_with_message( "moving_suite_list_.size() > 1 but WHATEVER can_prepend = true && can_append == true!" );
			}
		}

		Is_prepend_ = can_prepend; //WARNING THIS VARIABLE ONLY HAVE MEANING IF moving_suite_list_.size()>1.
		/////////////////////////Check that moving_suite_list_ is correctly order/////////////////////////

		//These vectors should be empty to begin with, but not harm to ensure this.
		force_syn_chi_res_list_.clear();
		force_north_ribose_list_.clear();
		force_south_ribose_list_.clear();

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGeneratorWrapper::~StepWiseRNA_RotamerGeneratorWrapper(){}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGeneratorWrapper::initialize_rotamer_generator_list(){

		Output_title_text( "Enter StepWiseRNA_RotamerGenerator::initialize_rotamer_generator_list", TR.Debug );
		Output_boolean( "sample_sugar_and_base1_ =  ", sample_sugar_and_base1_, TR.Debug );
		Output_boolean( "  sample_sugar_and_base2_ =  ", sample_sugar_and_base2_, TR.Debug );
		Output_boolean( "  Is_prepend_( only_meaningful_if_#res>1) =  ", Is_prepend_, TR.Debug );
		Output_boolean( "  sample_extra_rotamers_ =  ", sample_extra_rotamers_, TR.Debug );
		Output_boolean( "  include_syn_chi_ =  ", include_syn_chi_, TR.Debug );
		Output_boolean( "  allow_syn_pyrimidine_ =  ", allow_syn_pyrimidine_, TR.Debug );
		TR.Debug << std::endl;
		Output_boolean( "exclude_alpha_beta_gamma_sampling_: ", exclude_alpha_beta_gamma_sampling_, TR.Debug );
		Output_boolean( "  extra_beta_: ", extra_beta_, TR.Debug );
		Output_boolean( " extra_epsilon_: ", extra_epsilon_, TR.Debug );
		Output_boolean( " extra_anti_chi_: ", extra_anti_chi_, TR.Debug );
		Output_boolean( " extra_syn_chi_: ", extra_syn_chi_, TR.Debug );
		TR.Debug << std::endl;
		Output_seq_num_list( "working_moving_suite_list_:", moving_suite_list_, TR.Debug );
		Output_seq_num_list( "working_force_syn_chi_res_list_:", force_syn_chi_res_list_, TR.Debug );
		Output_seq_num_list( "working_force_north_ribose_list_:", force_north_ribose_list_, TR.Debug );
		Output_seq_num_list( "working_force_south_ribose_list_:", force_south_ribose_list_, TR.Debug );

		for ( Size n = 1; n <= force_north_ribose_list_.size(); n++ ){
			if ( force_south_ribose_list_.has_value( force_north_ribose_list_[n] ) ){
				utility_exit_with_message( "seq_num =  " + ObjexxFCL::string_of( force_north_ribose_list_[n] ) + " is in both force_north_ribose_list_ and force_south_ribose_list_! " );
			}
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////
		for ( Size list_position = rotamer_generator_list_.size(); list_position >= 1; list_position-- ){
			TR.Debug << "list_position =  " << list_position << " working_moving_suite =  " << moving_suite_list_[list_position] << std::endl;
			rotamer_generator_list_[list_position] = setup_rotamer_generator( list_position ); //This assumes that rotamer_generator[list_position+1] is setup and properly initialized
		}

		Output_title_text( "Exit StepWiseRNA_RotamerGenerator::initialize_rotamer_generator_list", TR.Debug );

	}

	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGeneratorOP const
	StepWiseRNA_RotamerGeneratorWrapper::setup_rotamer_generator( Size const list_position ){

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool sample_lower_sugar_and_base, sample_upper_sugar_and_base;

		if ( moving_suite_list_.size() == 1 ){
			sample_lower_sugar_and_base = sample_sugar_and_base1_;
			sample_upper_sugar_and_base = sample_sugar_and_base2_;
		} else{
			if ( Is_prepend_ ){
				if ( list_position == 1 ){ //5' most sample res.
					sample_lower_sugar_and_base = sample_sugar_and_base1_;
					sample_upper_sugar_and_base = false;
				} else if ( list_position == moving_suite_list_.size() ) { //3' most sample res.
					sample_lower_sugar_and_base = true;
					sample_upper_sugar_and_base = sample_sugar_and_base2_;
				} else{
					sample_lower_sugar_and_base = true;
					sample_upper_sugar_and_base = false;
				}
			} else{ //Append
				if ( list_position == 1 ){ //3' most sample res.
					sample_lower_sugar_and_base = false;
					sample_upper_sugar_and_base = sample_sugar_and_base2_;
				} else if ( list_position == moving_suite_list_.size() ) { //5' most sample res.
					sample_lower_sugar_and_base = sample_sugar_and_base1_;
					sample_upper_sugar_and_base = true;
				} else{
					sample_lower_sugar_and_base = false;
					sample_upper_sugar_and_base = true;
				}
			}
		}

		// Following function does not necessarily look in the pose for the information it needs -- it can look at previously set up
		// rotamer generators!
		core::Size const lower_res_puckerstate  =  Get_residue_pucker_state_internal( pose_, list_position, "lower", sample_lower_sugar_and_base );
		core::Size const upper_res_puckerstate  =  Get_residue_pucker_state_internal( pose_, list_position, "upper", sample_upper_sugar_and_base );

		///////////////////////////////////////////////////////////////////////////////////////////////////////////////
		bool const Is_bulge  =  ( list_position == 1 ) ? false : true;
		TR.Debug << "list_position = " << list_position; Output_boolean( " Is_bulge =  ", Is_bulge, TR.Debug ); TR.Debug << std::endl;

		bool const sample_extra_rotamers =  ( Is_bulge ) ? false : sample_extra_rotamers_;
		bool const fast  =  ( list_position == 1 ) ? fast_ : false; //Only make the sampling suite fast..

		//StepWiseRNA_RotamerGeneratorOP rotamer_generator  =  new StepWiseRNA_RotamerGenerator( moving_suite_list_[list_position], lower_res_puckerstate, upper_res_puckerstate, Is_bulge);

		StepWiseRNA_RotamerGeneratorOP rotamer_generator  =  new StepWiseRNA_RotamerGenerator( moving_suite_list_[list_position],
																																										sample_lower_sugar_and_base,
																																										sample_upper_sugar_and_base,
																																										lower_res_puckerstate,
																																										upper_res_puckerstate );

		if ( Is_bulge ){
			rotamer_generator->set_sample_chi_torsion( false );
		} else{
			rotamer_generator->set_sample_chi_torsion( sample_chi_torsion_ );
		}

		rotamer_generator->set_fast( fast );
		rotamer_generator->set_sample_extra_rotamers( sample_extra_rotamers );

		rotamer_generator->set_include_syn_chi( include_syn_chi_ );
		rotamer_generator->set_allow_syn_pyrimidine( allow_syn_pyrimidine_ );
		rotamer_generator->set_force_syn_chi_res_list( force_syn_chi_res_list_ );

		rotamer_generator->set_bin_size( bin_size_ );
		rotamer_generator->set_extra_epsilon(   (  ( Is_bulge ) ? false : extra_epsilon_  )  );
		rotamer_generator->set_extra_beta(      (  ( Is_bulge ) ? false : extra_beta_     )  );
		rotamer_generator->set_extra_anti_chi(  (  ( Is_bulge ) ? false : extra_anti_chi_ )  );
		rotamer_generator->set_extra_syn_chi(   (  ( Is_bulge ) ? false : extra_syn_chi_  )  );
		rotamer_generator->set_exclude_alpha_beta_gamma_sampling( exclude_alpha_beta_gamma_sampling_ );

		rotamer_generator->set_choose_random( choose_random_ );

		rotamer_generator->initialize_rotamer_generator( pose_ );

		return rotamer_generator;
	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGeneratorWrapper::set_fast( bool const & fast ){
		fast_ = fast;

	}

	////////////////////////////////////////////////////////////////////////
	//Rewrote this function on April 30th, 2011
	//The returned pucker_state value serves two roles
	//First 		if sample_sugar_pucker == True,  then pucker_state then specifies the pucker to be sampled by rotamer_generator_list_[list_position]
	//Second		if sample_sugar_pucker == False, then pucker_state then specifies the pucker determine from input pose or rotamer_generator_list_[list_position+1]

	core::Size
	StepWiseRNA_RotamerGeneratorWrapper::Get_residue_pucker_state_internal( core::pose::Pose const & pose, Size const list_position, std::string const which_sugar, bool sample_sugar_pucker ) const{



		Size const working_moving_suite  =  moving_suite_list_[list_position];
		Size const working_moving_pucker =  ( which_sugar == "lower" ) ? working_moving_suite : working_moving_suite + 1; //April 29, 2011


		if ( sample_sugar_pucker ){

			if ( force_north_ribose_list_.has_value( working_moving_pucker ) ) return NORTH;

			if ( force_south_ribose_list_.has_value( working_moving_pucker ) ) return SOUTH;

			return WHATEVER;

		} else{

			core::Size pucker_state;

			bool Is_first_of_multiple_res_plus_internal_case = false;

			if ( ( list_position == 1 && moving_suite_list_.size() > 1 ) ){
				if ( which_sugar == "lower" && Is_prepend_ )  Is_first_of_multiple_res_plus_internal_case = true; //Must be internal since prepend, since lower sugar is not sampled.
				if ( which_sugar == "upper" && !Is_prepend_ ) Is_first_of_multiple_res_plus_internal_case = true; //Must be internal since append, since lower sugar is not sampled.
			}

			if ( list_position == rotamer_generator_list_.size() || Is_first_of_multiple_res_plus_internal_case ){ //the ribose belongs to a prexisting input pose.

				pucker_state = Get_residue_pucker_state( pose, working_moving_pucker, true );

			} else{ //these correspond to sugars that were sampled by the rotamer_generator_list[list_position+1]

				if ( moving_suite_list_.size() < 2 ) utility_exit_with_message( "moving_suite_list_.size() < 2" );


				if ( Is_prepend_ ){

					if ( which_sugar != "upper"	)	utility_exit_with_message( "which_sugar != \"upper\"" );		//upper sugar pucker is the one not sampled if prepend

					if ( ( working_moving_suite + 1 )  !=  rotamer_generator_list_[list_position + 1]->moving_suite() ){
						utility_exit_with_message( "( moving_suite + 1 )  !=  rotamer_generator_list_[list_position + 1]->moving_suite()" );
					}

					pucker_state = rotamer_generator_list_[list_position + 1]->pucker_state( "lower" ); //if prepend, then lower of rotamer_generator_list_[list_position+1] is upper of rotamer_generator_list_[list_position]


				} else{ //Append

					if ( which_sugar != "lower" ) utility_exit_with_message( "which_sugar != \"lower\"" );		//lower sugar pucker is the one not sampled if append

					if ( ( working_moving_suite - 1 )  !=  rotamer_generator_list_[list_position + 1]->moving_suite() ){
						utility_exit_with_message( "( moving_suite - 1 )  !=  rotamer_generator_list_[list_position + 1]->moving_suite()" );
					}

					pucker_state = rotamer_generator_list_[list_position + 1]->pucker_state( "upper" ); //if append upper of rotamer_generator_list_[list_position+1] is lower of rotamer_generator_list_[list_position]
				}
			}

			if ( force_north_ribose_list_.has_value( working_moving_pucker ) && pucker_state != NORTH ){
				utility_exit_with_message( "force_north_ribose_list_.has_value( working_moving_pucker ) && pucker_state != NORTH, working_moving_pucker = " + ObjexxFCL::string_of( working_moving_pucker ) );
			}

			if ( force_south_ribose_list_.has_value( working_moving_pucker ) && pucker_state != SOUTH ){
				utility_exit_with_message( "force_south_ribose_list_.has_value( working_moving_pucker ) && pucker_state != SOUTH, working_moving_pucker = " + ObjexxFCL::string_of( working_moving_pucker ) );
			}

			return pucker_state;

		}
	}


	////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_RotamerGeneratorWrapper::has_another_rotamer() const{

		if ( choose_random_ ) return true; // can always get another random rotamer!

		for ( Size list_position = rotamer_generator_list_.size(); list_position >= 2; list_position-- ){
			if ( rotamer_generator_list_[list_position]->num_rotamer_centers()  !=  rotamer_generator_list_[list_position]->group_rotamer() ) return true;
		}

		return ( rotamer_generator_list_[1]->has_another_rotamer() );
	}


	////////////////////////////////////////////////////////////////////////
	utility::vector1< Torsion_Info >
	StepWiseRNA_RotamerGeneratorWrapper::get_next_rotamer(){

		using namespace core::id;

		assert( has_another_rotamer() );

		utility::vector1< bool > need_initialization_list( moving_suite_list_.size(), false );

		if ( choose_random_ ){

			for ( Size n = rotamer_generator_list_.size(); n >= 1; n-- ){ //Important to initialize in this order.

				// following gives correct behavior, I think, but is insanely expensive -- will need to optimize if/when we really get to dinucleotide case!
				if ( moving_suite_list_.size() > 1 ) rotamer_generator_list_[n] = setup_rotamer_generator( n ); //This assumes that rotamer_generator[n+1] is setup and initialized

				rotamer_generator_list_[n]->update_to_next_rotamer(); // will actually update_to_random_rotamer.
			}

		} else {

			//Update rotamer_generator_list_.
			for ( Size n = 1; n <= rotamer_generator_list_.size(); n++ ){
				if ( rotamer_generator_list_[n]->has_another_rotamer() ){
					rotamer_generator_list_[n]->update_to_next_rotamer();
					break;
				} else {
					need_initialization_list[n] = true;
					TR.Debug << "Bulge_rotamer_ID =  " << rotamer_generator_list_[2]->group_rotamer() << std::endl;
				}
			}

			//Need to reinitialize rotamer_generator since the pucker of other suites might change.
			for ( Size n = rotamer_generator_list_.size(); n >= 1; n-- ){ //Important to initialize in this order.
				if ( need_initialization_list[n] ){
					rotamer_generator_list_[n] = setup_rotamer_generator( n ); //This assumes that rotamer_generator[n+1] is setup andinitialized
				}
			}

		}

		utility::vector1< Torsion_Info > all_rotamer_list;
		all_rotamer_list.clear();

		//Update all_rotamer_list
		for ( Size n = 1; n <= moving_suite_list_.size(); n++ ){
			utility::vector1< Torsion_Info > suite_rotamer_list = rotamer_generator_list_[n]->get_current_rotamer();
			for ( Size i = 1; i <= suite_rotamer_list.size(); i++ ){
				all_rotamer_list.push_back( suite_rotamer_list[i] );
			}
		}

//		TR.Debug << "all_rotamer_list.size() =  " << all_rotamer_list.size() << std::endl;

		return all_rotamer_list;
	}

	/////////////////////////////////////////////////////////////////////////
	core::Size
	StepWiseRNA_RotamerGeneratorWrapper::group_rotamer( core::Size const list_position ) {
		return rotamer_generator_list_[list_position]->group_rotamer();
	}

	/////////////////////////////////////////////////////////////////////////
	core::Size
	StepWiseRNA_RotamerGeneratorWrapper::subgroup_rotamer( core::Size const list_position ) {
		return rotamer_generator_list_[list_position]->subgroup_rotamer();
	}

	/////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGeneratorWrapper::set_choose_random( bool const  setting ){
		choose_random_ = setting;
	}

}
}
}
