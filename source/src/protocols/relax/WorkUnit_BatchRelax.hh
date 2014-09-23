// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/realx/WorkUnit_BatchRelax.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_relax_WorkUnit_BatchRelax_hh
#define INCLUDED_protocols_relax_WorkUnit_BatchRelax_hh

#include <protocols/wum/WorkUnitBase.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace relax {



class WorkUnit_BatchRelax;
typedef utility::pointer::shared_ptr< WorkUnit_BatchRelax > WorkUnit_BatchRelaxOP;
typedef utility::pointer::shared_ptr< WorkUnit_BatchRelax const > WorkUnit_BatchRelaxCOP;

class WorkUnit_BatchRelax : public protocols::wum::WorkUnit_SilentStructStore {
public:
	WorkUnit_BatchRelax();
	~WorkUnit_BatchRelax();

	// @brief Run the workunit - overloaded by children of this class
	virtual void run();

	// @brief Hook for post processing such as rescoring etc.
	virtual void pre_process();

	// @brief Hook for post processing such as rescoring etc.
	virtual void post_process();

	virtual protocols::wum::WorkUnitBaseOP clone() const;

	void
	set_native_pose( core::pose::PoseCOP native_pose);

	protected:

		core::pose::PoseCOP native_pose_;
		core::scoring::ScoreFunctionOP scorefxn_;
	private:
};


class WorkUnit_BatchRelax_and_PostRescore : public WorkUnit_BatchRelax{
public:
	WorkUnit_BatchRelax_and_PostRescore();
	~WorkUnit_BatchRelax_and_PostRescore();

	virtual protocols::wum::WorkUnitBaseOP clone() const;

	void set_defaults();

	// @brief Hook for post processing such as rescoring etc.
	virtual void pre_process();

	// @brief Hook for post processing such as rescoring etc.
	virtual void post_process();

	// @brief Hook for post processing such as rescoring etc.
	virtual void rescore_all_decoys();

	protected:

	void trim();

	private:

	core::Real trim_proportion_;
	core::pose::Pose ref_pose_ ;
};




}
}

#endif
