// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWisePoseMinimizer
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/StepWisePoseMinimizer.hh>
#include <protocols/swa/StepWiseUtil.hh>

//////////////////////////////////
#include <core/types.hh>
//#include <core/io/silent/BinaryProteinSilentStruct.hh>
//#include <core/io/silent/SilentFileData.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/kinematics/MoveMap.hh>

//#include <utility/basic_sys_util.hh>
#include <time.h>


#include <string>

using namespace core;
using core::Real;

namespace protocols {
namespace swa {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWisePoseMinimizer::StepWisePoseMinimizer( PoseList & pose_list, utility::vector1< Size > const & moving_residues ):
    pose_list_( pose_list ),
    moving_residues_( moving_residues ),
    fa_scorefxn_( core::scoring::getScoreFunction() )
  {
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWisePoseMinimizer::~StepWisePoseMinimizer()
  {}

  //////////////////////////////////////////////////////////////////////////
  void
  StepWisePoseMinimizer::apply( core::pose::Pose & pose )
  {
    using namespace core::optimization;
    using namespace core::scoring;
    using namespace core::scoring::constraints;
    using namespace core::pose;

		clock_t const time_start( clock() );

		ConstraintSetOP cst_set = pose.constraint_set()->clone();

    AtomTreeMinimizer minimizer;
    float const dummy_tol( 0.00000025);
    bool const use_nblist( true );
    MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
    options.nblist_auto_update( true );

    kinematics::MoveMap mm;
    mm.set_jump( true );
    mm.set_bb( true );
    mm.set_chi( true );

    // Limited minimize... not as effective.
    if ( false ) limit_minimize( mm, pose.total_residue() );

    Size count( 1 );

    for ( PoseList::iterator iter = pose_list_.begin(); iter != pose_list_.end(); iter++ ) {

      PoseOP & pose_op( iter->second );
      pose = *pose_op; // This copy is to allow for easy graphic visualization.

			pose.constraint_set( cst_set );

			Real const score_original = (*fa_scorefxn_)( pose );

			//Should we also save the original constraints and put them back in later?
			core::scoring::constraints::add_coordinate_constraints( pose );

      std::cout << "Minimizing decoy " << count++ << " out of " << pose_list_.size() << std::endl;
      fa_scorefxn_->set_weight( coordinate_constraint, 0.1 );
      minimizer.run( pose, mm, *fa_scorefxn_, options );
      fa_scorefxn_->set_weight( coordinate_constraint, 0.0 );
      minimizer.run( pose, mm, *fa_scorefxn_, options );

			setPoseExtraScores( pose, "score_orig", score_original );

      std::string const & tag( iter->first );

      protocols::swa::output_silent_struct( pose, get_native_pose(), silent_file_, tag );

			// Running into file locking issues
			//			utility::sys_sleep( 0.5 );
			//exit( 0 );

      // Might was well replace pose in original list.
      *pose_op = pose;
    }

		std::cout << "Total time in StepWisePoseMinimizer: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

  }

std::string
StepWisePoseMinimizer::get_name() const {
	return "StepWisePoseMinimizer";
}

  ///////////////////////////////////////////////////////////////////////////
  void
  StepWisePoseMinimizer::limit_minimize( core::kinematics::MoveMap & mm, Size const & nres ) {

    mm.set_jump( false );
    mm.set_bb( false );
    mm.set_chi( false );

    for ( Size n = 1; n <= moving_residues_.size(); n++ ) {
      mm.set_bb( n, true );
      mm.set_chi( n, true );
      if ( n > 1 ){
				mm.set_bb ( n-1, true );
				mm.set_chi( n-1, true );
      }
      if ( n < nres ){
				mm.set_bb ( n+1, true );
				mm.set_chi( n+1, true );
      }
    }

  }

  //////////////////////////////////////////////////////////////////////////
  void
  StepWisePoseMinimizer::set_silent_file( std::string const & silent_file ){
    silent_file_ = silent_file;
  }

  //////////////////////////////////////////////////////////////////////////
	void
	StepWisePoseMinimizer::set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn ){
		fa_scorefxn_ = scorefxn;
	}

	//	void
	//	StepWisePoseMinimizer::set_constraint_set( core::scoring::constraints::ConstraintSetOP const & cst_set ){
	//		cst_set_ = cst_set;
	//	}


}
}
