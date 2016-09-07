// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/MPI_LoopHashRefine_Emperor.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_MPI_LoopHashRefine_Emperor_hh
#define INCLUDED_protocols_loophash_MPI_LoopHashRefine_Emperor_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/wum/MPI_WorkUnitManager.hh>
#include <protocols/loophash/MPI_LoopHashRefine.hh>

#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <string>
#include <vector>

namespace protocols {
namespace loophash {


class MPI_LoopHashRefine_Emperor: public MPI_LoopHashRefine {

public:

	MPI_LoopHashRefine_Emperor():
		MPI_LoopHashRefine( 'E' )
	{
		set_defaults();
	}

	~MPI_LoopHashRefine_Emperor() override= default;

	void set_defaults();

public:

	void go() override;

protected:

	void init() override;

	void process_inbound_wus() override;

	void process_outbound_wus() override;

	// adding arriving structures to library
	bool add_structures_to_library( protocols::wum::SilentStructStore &new_structs, std::string add_algorithm = "" ) override;
private:

	core::Size max_emperor_lib_round_;
};


} // namespace loops
} // namespace protocols


#endif


