

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/WorkUnit_Sampler.hh
/// @brief
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

#ifndef INCLUDED_protocols_mpi_refinement_WorkUnit_Loop_hh
#define INCLUDED_protocols_mpi_refinement_WorkUnit_Loop_hh

#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <utility/vector1.hh>
#include <set>  // required for std::set

namespace protocols {
namespace mpi_refinement {

////////////////////////////////////////////
//////// WorkUnit loophash
class WorkUnit_LoopHash : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_LoopHash( core::Size start_ir = 0, core::Size end_ir = 0, core::Size ssid = 0,
		core::Size is_global = 0);

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_LoopHash( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	virtual void init_from_cmd( const core::Size );

	void set_start( core::Size start_ir ){ header.extra_data_1_ = start_ir; }
	void set_end( core::Size end_ir ){ header.extra_data_2_ = end_ir; }
	void set_ssid( core::Size ssid ){ header.extra_data_3_ = ssid; }
	void set_global( core::Size is_global ){ header.extra_data_4_ = is_global; }

protected:

	core::Size get_start()  const { return header.extra_data_1_; }
	core::Size get_end()    const { return header.extra_data_2_; }
	core::Size get_ssid()   const { return header.extra_data_3_; }
	core::Size get_global() const { return header.extra_data_4_; }

	void set_defaults();

private:
	protocols::loophash::LoopHashLibraryOP library_;

};

////////////////////////////////////////////
//////// WorkUnit FragInsert
class WorkUnit_FragInsert : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_FragInsert( core::Size const nsteps = 0,
		core::Size const scoretype = 0,
		core::Size const res1 = 0,
		core::Size const res2 = 0,
		bool const fullatom = false
	);
	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_FragInsert( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	//void set_nstruct( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_nsteps( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_scoretype( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_res1( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_res2( core::Size const setting ){ header.extra_data_4_ = setting; }

protected:

	//core::Size get_nstruct() const { return header.extra_data_1_; }
	core::Size get_nsteps() const { return header.extra_data_1_; }
	core::Size get_scoretype() const { return header.extra_data_2_; }
	core::Size get_res1() const { return header.extra_data_3_; }
	core::Size get_res2() const { return header.extra_data_4_; }

	void set_defaults();

private:
	//std::set< core::Size > pos_from_cmd_;
	bool fullatom_;

};

////////////////////////////////////////////
//////// WorkUnit KicCloser
class WorkUnit_KicCloser : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_KicCloser( core::Size const nsteps = 0,
		core::Size const scoretype = 0,
		core::Size const res1 = 0,
		core::Size const res2 = 0,
		bool const kicclose = true
	);
	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_KicCloser( *this ));
	}

	// Pure virtual functions
	virtual void run();

	void set_nsteps( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_scoretype( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_res1( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_res2( core::Size const setting ){ header.extra_data_4_ = setting; }

protected:
	core::Size get_nsteps() const { return header.extra_data_1_; }
	core::Size get_scoretype() const { return header.extra_data_2_; }
	core::Size get_res1() const { return header.extra_data_3_; }
	core::Size get_res2() const { return header.extra_data_4_; }

	void set_defaults();

private:
	bool kicclose_;

};

////////////////////////////////////////////
//////// WorkUnit PartialAbinitio
class WorkUnit_PartialAbinitio : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_PartialAbinitio( core::Size const nsteps = 0,
		bool const reconstruct = false
	);
	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_PartialAbinitio( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	void set_nsteps( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_res1( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_res2( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_reconstruct( core::Size const setting ){ header.extra_data_4_ = setting; }

protected:
	core::Size get_nsteps() const { return header.extra_data_1_; }
	core::Size get_res1() const { return header.extra_data_2_; }
	core::Size get_res2() const { return header.extra_data_3_; }
	core::Size get_reconstruct() const { return header.extra_data_4_; }

	void set_defaults();

private:
	/*
	void
	setup_score( core::scoring::ScoreFunctionOP score0,
	core::scoring::ScoreFunctionOP score1,
	core::scoring::ScoreFunctionOP score2,
	core::scoring::ScoreFunctionOP score3,
	core::scoring::ScoreFunctionOP score5
	) const;
	*/

	core::scoring::ScoreFunctionOP
	setup_score( std::string sfxn_name ) const;

};

}
}

#endif

