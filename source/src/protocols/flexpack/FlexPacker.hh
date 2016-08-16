// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/flexpack/
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_flexpack_FlexPacker_hh
#define INCLUDED_protocols_flexpack_FlexPacker_hh

// Unit Headers
#include <protocols/flexpack/FlexPacker.fwd.hh>

#include <protocols/moves/Mover.hh>

// Project headers
#include <core/fragment/Frame.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

/// Utility headers

#include <utility/vector1.hh>


namespace protocols {
namespace flexpack {

class FlexPacker : public moves::Mover
{
public:
	//FlexPacker();


	FlexPacker(
		core::pack::task::PackerTaskCOP task,
		utility::vector1< core::fragment::FrameOP > const & frames,
		core::scoring::ScoreFunctionCOP scorefxn
	);


	virtual ~FlexPacker();

	virtual
	void
	apply( core::pose::Pose & pose );

	virtual std::string get_name() const;

	void
	set_sfxn( core::scoring::ScoreFunctionCOP sfxn );

	void
	set_task( core::pack::task::PackerTaskCOP task );

	void
	set_frames( utility::vector1< core::fragment::FrameOP > const & frames );

	void
	set_task_factory( core::pack::task::TaskFactoryOP factory );

private:

	utility::vector1< core::fragment::FrameCOP > frames_;

	core::pack::task::TaskFactoryOP factory_;

	core::pack::task::PackerTaskCOP task_;
	core::scoring::ScoreFunctionCOP scorefxn_;

};

}
}

#endif
