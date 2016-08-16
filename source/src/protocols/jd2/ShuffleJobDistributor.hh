// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/ShuffleJobDistributor.hh
/// @brief  implementation of ShuffleFileSystemJobDistributor
/// @author Mike Tyka

#ifndef INCLUDED_protocols_jd2_ShuffleJobDistributor_hh
#define INCLUDED_protocols_jd2_ShuffleJobDistributor_hh

// Unit headers
#include <protocols/jd2/FileSystemJobDistributor.hh>

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

class ShuffleFileSystemJobDistributor : public FileSystemJobDistributor
{
protected:
	/// @brief ctor is protected; singleton pattern
	ShuffleFileSystemJobDistributor();

	virtual void handle_interrupt() {}

public:
	virtual ~ShuffleFileSystemJobDistributor();


	friend class JobDistributorFactory; //ctor access

	virtual
	core::Size
	get_new_job_id();

	virtual
	void
	mark_current_job_id_for_repetition();

	core::Size                     next_random_job() { return next_random_job_; }

private:
	utility::vector1< core::Size > scrambled_job_order_;
	core::Size                     next_random_job_;
};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_ShuffleFileSystemJobDistributor_HH
