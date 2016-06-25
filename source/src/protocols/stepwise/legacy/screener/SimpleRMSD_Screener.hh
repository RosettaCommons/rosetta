// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/screener/SimpleRMSD_Screener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_SimpleRMSD_Screener_HH
#define INCLUDED_protocols_stepwise_screener_SimpleRMSD_Screener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/legacy/screener/SimpleRMSD_Screener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/types.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

class SimpleRMSD_Screener: public stepwise::screener::StepWiseScreener {

public:

	//constructor
	SimpleRMSD_Screener( pose::Pose const & pose,
		utility::vector1< Size > const & calc_rms_res,
		core::pose::PoseCOP native_pose,
		core::Real const rmsd_cutoff,
		bool const force_align = false );

	//destructor
	~SimpleRMSD_Screener();

public:

	std::string
	name() const { return "SimpleRMSD_Screener"; }

	stepwise::screener::StepWiseScreenerType
	type() const { return stepwise::screener::SIMPLE_RMSD; }

	bool
	check_screen();

private:

	void
	initialize_corresponding_atom_id_map( core::pose::Pose const & pose );

private:

	pose::Pose const & pose_;

	utility::vector1< Size > calc_rms_res_;
	core::pose::PoseCOP native_pose_;
	core::Real const rmsd_cutoff_;
	bool const force_align_;
	bool const cluster_by_all_atom_rmsd_;

	std::map< core::id::AtomID, core::id::AtomID > corresponding_atom_id_map_;

};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
