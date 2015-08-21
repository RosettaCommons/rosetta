// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/BOINCJobDistributor.hh
/// @brief  implementation of BOINCJobDistributor
/// @author Mike Tyka

#ifndef INCLUDED_protocols_jd2_BOINCJobDistributor_hh
#define INCLUDED_protocols_jd2_BOINCJobDistributor_hh

// Unit headers
#include <protocols/jd2/ShuffleJobDistributor.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.fwd.hh>

#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

class BOINCJobDistributor : public ShuffleFileSystemJobDistributor
{
protected:
	/// @brief ctor is protected; singleton pattern
	BOINCJobDistributor();

	virtual
	void
	job_failed( core::pose::Pose & pose, bool will_retry );

	virtual
	void
	job_succeeded( core::pose::Pose & pose, core::Real run_time, std::string const & tag );

	virtual
	void
	mark_current_job_id_for_repetition();

	virtual
	void
	begin_critical_section();

	virtual
	void
	end_critical_section();

	core::Size total_completed_nstruct_; // for graphics display

public:
	virtual ~BOINCJobDistributor();
	void checkpoint_read();
	void checkpoint_write();
	void checkpoint_clear();

	friend class JobDistributorFactory; //ctor access

	virtual
	void
	go( protocols::moves::MoverOP mover );

	virtual
	core::Size
	get_new_job_id();
};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_BOINCJobDistributor_HH
