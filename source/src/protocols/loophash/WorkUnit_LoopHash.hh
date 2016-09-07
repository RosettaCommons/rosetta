// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/WorkUnit_LoopHash.hh
/// @brief
/// @author Mike Tyka

#ifndef INCLUDED_protocols_loophash_WorkUnit_LoopHash_hh
#define INCLUDED_protocols_loophash_WorkUnit_LoopHash_hh

#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/loophash/LoopHashLibrary.hh>

#include <utility/vector1.hh>

namespace protocols {
namespace loophash {


class WorkUnit_LoopHash;
typedef utility::pointer::shared_ptr< WorkUnit_LoopHash > WorkUnit_LoopHashOP;
typedef utility::pointer::shared_ptr< WorkUnit_LoopHash const > WorkUnit_LoopHashCOP;

class WorkUnit_LoopHash: public protocols::wum::WorkUnit_SilentStructStore {
public:
	// initialize only via this
	WorkUnit_LoopHash( core::Size start_ir=0, core::Size end_ir=0, core::Size ssid=0 );

	// @brief Run the workunit - overloaded by children of this class
	void run() override;

	protocols::wum::WorkUnitBaseOP clone() const override {
		runtime_assert( library_ != nullptr );
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_LoopHash( *this ) );
	}

	void init_from_cmd( const core::Size mpi_rank );
	void set_start( core::Size start_ir ){ header.extra_data_1_ = start_ir; }
	void set_end( core::Size end_ir ){ header.extra_data_2_ = end_ir; }
	void set_ssid( core::Size ssid ){ header.extra_data_3_ = ssid; }

protected:

	core::Size get_start(){ return header.extra_data_1_; }
	core::Size get_end(){ return header.extra_data_2_; }
	core::Size get_ssid(){ return header.extra_data_3_; }

	void set_defaults();
private:
	LoopHashLibraryOP library_;
};


}
}

#endif

