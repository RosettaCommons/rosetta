// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseFilterer
/// @brief Not particularly fancy, just filters a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/enumerate/protein/StepWiseFilterer.hh>

//////////////////////////////////
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <list>

#include <utility/vector1.hh>

//Auto Headers



using namespace core;
using core::Real;

namespace protocols {
namespace stepwise {
namespace enumerate {
namespace protein {


  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseFilterer::StepWiseFilterer():
		final_number_( 400 )
  {
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseFilterer::~StepWiseFilterer()
  {}

  //////////////////////////////////////////////////////////////////////////
	void
  StepWiseFilterer::filter( PoseList & pose_list, PoseList & filter_pose_list ) const
	{

		using namespace core::scoring;
		using namespace core::pose;

		static ScoreFunctionOP fa_scorefxn = core::scoring::getScoreFunction();

		typedef std::list < std::pair< Real, std::string > >  ScoreList;
		ScoreList score_list;

		for ( PoseList::iterator iter = pose_list.begin(); iter != pose_list.end(); iter++ ) {
			PoseOP & pose_op( iter->second );
			Pose & pose( *pose_op );
			Real score = (*fa_scorefxn)( pose );

			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////
			score -=  pose.energies().total_energies()[ fa_rep ]; /*downweight fa_rep */
			/////////////////////////////////////////////////////////////////////////////
			/////////////////////////////////////////////////////////////////////////////

			std::string const & tag( iter->first );

			score_list.push_back( std::make_pair( score, tag ) );
		}

		//	std::cout << "LIST OF POSES " << pose_list.size() << " " << score_list.size() << std::endl;

		score_list.sort();

		Size count( 0 );
		for ( ScoreList::const_iterator iter=score_list.begin(); iter != score_list.end(); iter++ ) {
			count++;
			if ( count > final_number_ ) break;
			std::string const  & tag = iter->second;
			filter_pose_list[ tag ] = pose_list[ tag ];
		}

	}

} //protein
} //enumerate
} //stepwise
} //protocols
