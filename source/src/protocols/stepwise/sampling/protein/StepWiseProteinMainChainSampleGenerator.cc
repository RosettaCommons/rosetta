// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinMainChainSampleGenerator
/// @brief Subclass of StepWisePoseSampleGenerator
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/sampling/protein/StepWiseProteinMainChainSampleGenerator.hh>
#include <core/id/TorsionID.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <utility/exit.hh>

//Auto Headers
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseProteinMainChainSampleGenerator::StepWiseProteinMainChainSampleGenerator(
							utility::vector1< core::id::TorsionID > const & which_torsions,
							utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists ):
		which_torsions_( which_torsions ),
		main_chain_torsion_set_lists_( main_chain_torsion_set_lists ),
		count_( 0 ),
		num_samples_( main_chain_torsion_set_lists_.size() )
  {
  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinMainChainSampleGenerator::reset(){
		count_ = 0;
	}

  //////////////////////////////////////////////////////////////////////////
	bool
  StepWiseProteinMainChainSampleGenerator::has_another_sample(){
		return (count_ < num_samples_);
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinMainChainSampleGenerator::get_next_sample( core::pose::Pose & pose )
	{

		count_++;

		if ( count_ > num_samples_ ) utility_exit_with_message( "Asked StepWiseProteinMainChainSampleGenerator for another sample but it does not have one!" );

		utility::vector1< Real > const & main_chain_torsion_set_list( main_chain_torsion_set_lists_[ count_ ] );

		for ( Size i = 1; i <= which_torsions_.size(); i++ ) {
			//std::cout << "SETTING TORSION " << which_torsions_[ i ] << "  to " << main_chain_torsion_set_list[ i ] << std::endl;
			pose.set_torsion( which_torsions_[ i ], main_chain_torsion_set_list[ i ] );
		}

	}

	///////////////////////////////////////////////////////
	Size
	StepWiseProteinMainChainSampleGenerator::size() const{
		return num_samples_;
	}

} //protein
} //sampling
} //stepwise
} //protocols

