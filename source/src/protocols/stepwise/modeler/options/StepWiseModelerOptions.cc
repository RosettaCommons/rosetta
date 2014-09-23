// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/options/StepWiseModelerOptions.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.hh>

#include <basic/Tracer.hh>

//Req'd on WIN32
#include <utility/tag/Tag.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.options.StepWiseModelerOptions" );

using namespace basic::options;
using namespace basic::options::OptionKeys;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace options {

	/////////////////////////////////////////////////////////////////////////////////////
	//Constructor
	StepWiseModelerOptions::StepWiseModelerOptions()
	{
		initialize_variables();
	}

	/////////////////////////////////////////////////////////////////////////////////////
	//Destructor
	StepWiseModelerOptions::~StepWiseModelerOptions()
	{
	}

	/// @brief copy constructor
	StepWiseModelerOptions::StepWiseModelerOptions( StepWiseModelerOptions const & src ) :
		ResourceOptions( src ),
		StepWiseBasicModelerOptions( src ),
		StepWiseProteinModelerOptions( src ),
		StepWiseRNA_ModelerOptions( src )
	{
		*this = src;
	}


	/////////////////////////////////////////////////////////////////////////////////////
	// If you add a variable, initialize it here, and include in operator = definition below!
	void
	StepWiseModelerOptions::initialize_variables(){
		StepWiseBasicModelerOptions::initialize_variables();
		StepWiseProteinModelerOptions::initialize_variables();
		StepWiseRNA_ModelerOptions::initialize_variables();
	}

	/// @brief clone the options
	StepWiseModelerOptionsOP
	StepWiseModelerOptions::clone() const
	{
		return StepWiseModelerOptionsOP( new StepWiseModelerOptions( *this ) );
	}

	/////////////////////////////////////////////////////////////////////////////////////
	StepWiseModelerOptions &
	StepWiseModelerOptions::operator = ( StepWiseModelerOptions const &  )
	{
		return *this;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseModelerOptions::initialize_from_command_line(){

		StepWiseBasicModelerOptions::initialize_from_command_line();
		StepWiseProteinModelerOptions::initialize_from_command_line();
		StepWiseRNA_ModelerOptions::initialize_from_command_line();

		if ( integration_test_mode() ) set_num_pose_minimize( 1 );
		if ( option[ out::nstruct].user() )	utility_exit_with_message( "Can not use -out::nstruct anymore, use -sampler_num_pose_kept and/or -num_pose_minimize for fine control." );

	}

	////////////////////////////////////////////////////////////////
	void
	StepWiseModelerOptions::setup_options_for_VDW_bin_checker( rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker ) const {
		user_input_VDW_bin_checker->set_VDW_rep_alignment_RMSD_CUTOFF ( VDW_rep_alignment_RMSD_CUTOFF() );
		user_input_VDW_bin_checker->set_VDW_rep_delete_matching_res( VDW_rep_delete_matching_res() );
		user_input_VDW_bin_checker->set_physical_pose_clash_dist_cutoff( VDW_rep_screen_physical_pose_clash_dist_cutoff() );
		user_input_VDW_bin_checker->set_output_pdb( dump_ );
	}

	///////////////////////////////////////////////////////////////////////////////////
	StepWiseModelerOptionsOP
	StepWiseModelerOptions::get_sampler_options() const {
		StepWiseModelerOptionsOP sampler_options = clone();
		if ( sampler_options->choose_random() )	sampler_options->set_cluster_rmsd( 0.0 ); // don't cluster.
		if ( sampler_options->integration_test_mode() ){
			sampler_options->set_sampler_num_pose_kept( 2 );
		 sampler_options->set_rmsd_screen( 1.0 ); // RNA_ConnectionSampler will initially have this off, but toggle true later in integration test.
		}
		return sampler_options;
	}

} //options
} //modeler
} //stepwise
} //protocols
