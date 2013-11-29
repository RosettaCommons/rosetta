// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/StepWiseResidueSampler.hh>
#include <protocols/swa/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/Ramachandran.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>

//GreenPacker
#include <protocols/simple_moves/GreenPacker.hh>
#include <protocols/simple_moves/GreenPacker.fwd.hh>

// AUTO-REMOVED #include <numeric/angle.functions.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/exit.hh>
// AUTO-REMOVED #include <time.h>

#include <string>

#include <utility/vector1.hh>

//Auto Headers



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

static basic::Tracer TR( "protocols.swa.StepWiseResidueSampler" ) ;

namespace protocols {
namespace swa {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseResidueSampler::StepWiseResidueSampler(
																								 utility::vector1< Size > const & moving_residues,
																								 utility::vector1< MainChainTorsionSetList > const & main_chain_torsion_set_lists ):
		moving_residues_( moving_residues ),
		main_chain_torsion_set_lists_( main_chain_torsion_set_lists ),
		scorefxn_( core::scoring::getScoreFunction() ),
		green_packer_( new protocols::simple_moves::GreenPacker ),
		silent_file_( "" ),
		sfd_( new core::io::silent::SilentFileData)
  {
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseResidueSampler::~StepWiseResidueSampler()
  {}

  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseResidueSampler::apply( core::pose::Pose & pose )
	{

		//Size which_res( 1 );
		//		Size count( 1 );

		clock_t const time_start( clock() );

		initialize_green_packer( pose.total_residue() );

		sample_residues( pose );

		std::cout << "Total time in StepWiseResidueSampler: " << clock() - time_start / CLOCKS_PER_SEC
							<< std::endl;

	}


	std::string
	StepWiseResidueSampler::get_name() const {
		return "StepWiseResidueSampler";
	}

	////////////////////////////////////////////////////////////////////////////
	void
	StepWiseResidueSampler::sample_residues( core::pose::Pose & pose )
	{

		 using namespace core::chemical;
		 using namespace core::scoring;
		 using namespace core::pose;

		 for ( Size k = 1; k <= main_chain_torsion_set_lists_.size(); k++ ) {

			 MainChainTorsionSetList const & main_chain_torsion_set_list( main_chain_torsion_set_lists_[ k ] );

			 for ( Size i = 1; i <= moving_residues_.size(); i++ ) {

				 MainChainTorsionSet const & main_chain_torsion_set( main_chain_torsion_set_list[ i ] );

				 Size const n = moving_residues_[ i ];
				 pose.set_phi( n,  main_chain_torsion_set.phi() );
				 pose.set_psi( n,  main_chain_torsion_set.psi() );
				 //Probably need to sample omega=0.0 for proline... easy fix, do it later.
				 pose.set_omega( n,  main_chain_torsion_set.omega() );

			 }

			 std::string const tag = "S_"+ lead_zero_string_of( k-1, 5 );
			 TR << " Decoy " << tag << " : " << k << " out of " << main_chain_torsion_set_lists_.size() << std::endl;

			 green_packer_->apply( pose );
			 (*scorefxn_)( pose );

			 output_silent_struct( pose, get_native_pose(), silent_file_, tag, sfd_ );

		 }

	}


	////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseResidueSampler::initialize_green_packer( Size const & nres )
	{
		using namespace protocols::moves;
		using namespace core::pack;
		using namespace core::pack::task;
		using namespace core::pack::task::operation;

		protocols::simple_moves::UserDefinedGroupDiscriminatorOP user_defined_group_discriminator( new UserDefinedGroupDiscriminator);
		utility::vector1< Size > group_ids;

		Size current_group = 0;
		Size spectator_group = 0;
		for (Size i = 1; i <= nres; i++ ) {
			bool found_it( false );
			for (Size k = 1; k <= moving_residues_.size(); k++ ) {
				if ( i == moving_residues_[k] ) {
					found_it = true;
					break;
				}
			}
			if (found_it ) {
				current_group = 0;
				TR << "GREENPACKER SAMPLER " << i << std::endl;
			} else {
				if ( current_group == 0 ) spectator_group++;
				current_group = spectator_group;
				TR << "GREENPACKER SPECTATOR   " << i <<  " --> group " << spectator_group << std::endl;
			}
			group_ids.push_back( current_group );
		}

		user_defined_group_discriminator->set_group_ids( group_ids );
		green_packer_->set_scorefunction( *scorefxn_ );
		green_packer_->set_group_discriminator( user_defined_group_discriminator );

		TaskFactoryOP initial_task_factory( new TaskFactory );
		initial_task_factory->push_back( new InitializeFromCommandline );
		initial_task_factory->push_back( new RestrictToRepacking );
		green_packer_->set_reference_round_task_factory( initial_task_factory );

		TaskFactoryOP general_task_factory( new TaskFactory );
		general_task_factory->push_back( new InitializeFromCommandline );
		general_task_factory->push_back( new RestrictToRepacking );
		green_packer_->set_task_factory( general_task_factory );

		//green_packer_->reset();
	}


  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseResidueSampler::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseResidueSampler::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		scorefxn_ = scorefxn;
	}

  //////////////////////////////////////////////////////////////////////////
	core::io::silent::SilentFileDataOP &
	StepWiseResidueSampler::silent_file_data(){
		return sfd_;
	}


}
}
