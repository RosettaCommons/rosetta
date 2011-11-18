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
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator_Wrapper.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/StepWiseUtil.hh>
//////////////////////////////////

#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using core::Real;
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_rotamer_generator_wrapper" );

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////


	StepWiseRNA_RotamerGenerator_Wrapper::StepWiseRNA_RotamerGenerator_Wrapper(
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & moving_suite_list,
								bool const & sample_sugar_and_base1,
								bool const & sample_sugar_and_base2,
								Real const bin_size ):
		pose_(pose),
		moving_suite_list_( moving_suite_list ),
		sample_sugar_and_base1_(sample_sugar_and_base1),
		sample_sugar_and_base2_(sample_sugar_and_base2),
		sample_extra_rotamers_( true ),
		fast_( false ),
		verbose_(false),
		rotamer_generator_list_(moving_suite_list_.size(), NULL),
		bin_size_( bin_size )
	{
		Output_title_text("Enter StepWiseRNA_RotamerGenerator_Wrapper Constructor");
		Output_boolean("sample_sugar_and_base1_= " , sample_sugar_and_base1_);
		Output_boolean("  sample_sugar_and_base2_= " , sample_sugar_and_base2_);
		Output_boolean("  sample_extra_rotamers_= " , sample_extra_rotamers_); std::cout << std::endl;
		Output_seq_num_list("working_moving_suite_list_:" , moving_suite_list_);

		initialize_rotamer_generator_list();

		//		if ( bin_size_ != 20.0 && bin_size_ != 10.0 && bin_size_ != 5.0 ) utility_exit_with_message( "Disallowed bin_size for rotamer generator" ) ;

		Output_title_text("Exit StepWiseRNA_RotamerGenerator_Wrappper Constructor");
	}
	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGenerator_Wrapper::~StepWiseRNA_RotamerGenerator_Wrapper(){}

	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_RotamerGeneratorOP const
	StepWiseRNA_RotamerGenerator_Wrapper::setup_rotamer_generator(Size const list_position){

		PuckerState const lower_res_puckerstate = Get_residue_pucker_state_internal(pose_, list_position, "lower");
		PuckerState const upper_res_puckerstate = Get_residue_pucker_state_internal(pose_, list_position, "upper");

		bool const Is_bulge = (list_position==1) ? false : true;
//		Output_boolean("Is_bulge= " , Is_bulge);
		bool const sample_extra_rotamers= (Is_bulge) ? false : sample_extra_rotamers_;
		bool const fast = (list_position==1) ? fast_ : false; //Only make the sampling suite fast..

		StepWiseRNA_RotamerGeneratorOP rotamer_generator = new StepWiseRNA_RotamerGenerator( pose_, moving_suite_list_[list_position], lower_res_puckerstate, upper_res_puckerstate, Is_bulge, bin_size_ );
		rotamer_generator->set_fast( fast );
		rotamer_generator->set_sample_extra_rotamers(sample_extra_rotamers);

		return rotamer_generator;
	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator_Wrapper::set_fast( bool const & fast ){
		fast_=fast;

		//Need this since the first time the rotamer_generator_list is initialize happens below set_fast() is called
		rotamer_generator_list_[1]->set_fast( fast_ ); //Only make the sampling suite fast..
	}

	////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_RotamerGenerator_Wrapper::initialize_rotamer_generator_list(){

		for(Size list_position=rotamer_generator_list_.size(); list_position>=1; list_position--){
			std::cout << "list_position= " << list_position << " working_moving_suite= " << moving_suite_list_[list_position] << std::endl;
			rotamer_generator_list_[list_position]=setup_rotamer_generator(list_position); //This assumes that rotamer_generator[list_position+1] is setup and probably initialized
		}

	}

	////////////////////////////////////////////////////////////////////////
  //Might need to modify this slightly for it to work in the INTERNAL + more than one moving res case.
	PuckerState
	StepWiseRNA_RotamerGenerator_Wrapper::Get_residue_pucker_state_internal( core::pose::Pose const & pose, Size const list_position, std::string const which_sugar) const{

		if(which_sugar=="lower"){ //lower
			if(sample_sugar_and_base1_) return ALL;
		}else{
			if(sample_sugar_and_base2_) return ALL;
		}

		Size const moving_suite = moving_suite_list_[list_position];

		if(list_position==rotamer_generator_list_.size()){ //Last moving_suite_list_position, the suite that is connected to the prexisting structure
			if(which_sugar=="lower"){
				return Get_residue_pucker_state( pose, moving_suite );
			}else{
				return Get_residue_pucker_state( pose, moving_suite + 1);
			}
		}

		//n!=1 case
		if(which_sugar=="lower"){ //Append case
			if( (moving_suite-1) != rotamer_generator_list_[list_position+1]->moving_suite()){
				utility_exit_with_message( "(moving_suite-1) != rotamer_generator_list_[list_position+1]->moving_suite()" );
			}
			return rotamer_generator_list_[list_position+1]->pucker_state("upper");
		}else{ //Prepend case
			if( (moving_suite+1) != rotamer_generator_list_[list_position+1]->moving_suite()){
				utility_exit_with_message( "(moving_suite+1) != rotamer_generator_list_[list_position+1]->moving_suite()" );
			}
			return rotamer_generator_list_[list_position+1]->pucker_state("lower");
		}
	}


	////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_RotamerGenerator_Wrapper::has_another_rotamer() const{

		for(Size list_position=rotamer_generator_list_.size(); list_position>=2; list_position--){
			if(rotamer_generator_list_[list_position]->num_rotamer_centers() != rotamer_generator_list_[list_position]->group_rotamer()) return true;
		}

		return ( rotamer_generator_list_[1]->has_another_rotamer() );
	}


	////////////////////////////////////////////////////////////////////////
	utility::vector1< Torsion_Info >
	StepWiseRNA_RotamerGenerator_Wrapper::get_next_rotamer(){

		using namespace core::id;

		assert( has_another_rotamer() );

		utility::vector1< bool > need_initialization_list(moving_suite_list_.size(), false);

		//Update rotamer_generator_list_.
		for(Size n=1; n<=rotamer_generator_list_.size(); n++){
			if(rotamer_generator_list_[n]->has_another_rotamer()){
				rotamer_generator_list_[n]->update_to_next_rotamer();
				break;
			}else{
				need_initialization_list[n]=true;
				std::cout << "Bulge_rotamer_ID= " << rotamer_generator_list_[2]->group_rotamer() << std::endl;
			}
		}

		//Need to reinitialize rotamer_generator since the pucker of other suites might change.
		for(Size n=rotamer_generator_list_.size(); n>=1; n--){ //Important to initialize in this order.
			if(need_initialization_list[n]==true){
				rotamer_generator_list_[n]=setup_rotamer_generator(n); //This assumes that rotamer_generator[n+1] is setup and probably initialized
			}
		}

		utility::vector1< Torsion_Info > all_rotamer_list;
		all_rotamer_list.clear();

		//Update all_rotamer_list
		for(Size n=1; n<=moving_suite_list_.size(); n++){
			utility::vector1< Torsion_Info > suite_rotamer_list=rotamer_generator_list_[n]->get_current_rotamer();
			for(Size i=1; i<=suite_rotamer_list.size(); i++){
				all_rotamer_list.push_back(suite_rotamer_list[i]);
			}
		}

//		std::cout << "all_rotamer_list.size()= " << all_rotamer_list.size() << std::endl;

		return all_rotamer_list;
	}

	core::Size
	StepWiseRNA_RotamerGenerator_Wrapper::group_rotamer(core::Size const list_position) {
		return rotamer_generator_list_[list_position]->group_rotamer();
	}

	core::Size
	StepWiseRNA_RotamerGenerator_Wrapper::subgroup_rotamer(core::Size const list_position) {
		return rotamer_generator_list_[list_position]->subgroup_rotamer();
	}

}
}
}
