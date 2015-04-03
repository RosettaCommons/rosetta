// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file looplooprelax_protocols.hh
/// @brief
/// @details
///
/// @author Srivatsan Raman


#ifndef INCLUDED_devel_ssrbrelax_SSrbrelax_hh
#define INCLUDED_devel_ssrbrelax_SSrbrelax_hh
#include <devel/ssrbrelax/SSRbClass.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
//#include <core/chemical/ResidueTypeSet.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.hh>
#include <protocols/frags/TorsionFragment.hh>
#include <protocols/frags/VallData.hh>

#include <numeric/xyzVector.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

namespace devel {
namespace ssrbrelax {

class RbRelax: public protocols::moves::Mover {
public:
	//	LoopRelax(
	//		core::scoring::ScoreFunctionOP scorefxn_in,
	//		core::pose::PoseOP pose_in
	//		) : Mover(), pose_(pose_in), scorefxn_(scorefxn_in)
	RbRelax():Mover()

	{}
//		set_default();
//		}

	/*
	void set_default() {
		set_default_mc();
		set_default_move_map();
	}
	*/
	virtual RbRelax* clone() const	{
		return new RbRelax(*this);
	}

	protocols::moves::MonteCarloOP get_mc();
	//	void set_default_move_map();
	//	void set_default_mc();

	void apply( core::pose::Pose & pose );
	void apply();


private:
	// protocol stuff
	core::pose::PoseOP pose_;
	core::scoring::ScoreFunctionOP scorefxn_;
	core::kinematics::MoveMapOP movemap_;
	protocols::moves::MonteCarloOP mc_;

	// default temperature for monte carlo
	core::Real m_Temperature_;

	void rbrelax_main(
										core::pose::Pose & pose
										);

	void segment_rb_move(
											 core::pose::Pose & pose,
											 devel::ssrbrelax::RbSegments & rbsegments,
											 std::map< Size, devel::frags::TorsionFragmentLibraryOP > & frag_libs
											 );

	void perturb_segment_and_close_loops(
																			 core::pose::Pose & pose,
																			 devel::ssrbrelax::RbSegments & rbsegments,
																			 devel::ssrbrelax::RbSegments & this_segment,
																			 std::map< Size, devel::frags::TorsionFragmentLibraryOP > & frag_libs
																			 );

	void perturb_segment(
											 core::pose::Pose & pose,
											 devel::ssrbrelax::RbSegments & this_segment,
											 int const & dof,
											 float const & stddev
											 );

	void set_rbrelax_allow_move_map(
																	core::kinematics::MoveMap & rb_move_map,
																	int const & flexible_jump
																	);

	void initialize_fragments(
														std::map< Size, devel::frags::TorsionFragmentLibraryOP > & frag_libs
														);

	void close_both_loops(
												core::pose::Pose & pose,
												devel::ssrbrelax::RbSegments & this_segment,
												std::map< Size, devel::frags::TorsionFragmentLibraryOP > & frag_libs
												);

	void refine_segment(
											core::pose::Pose & pose,
											devel::ssrbrelax::RbSegments & this_segment
											);


}; // class RbRelax

} // SSrbrelax
} // devel

#endif
