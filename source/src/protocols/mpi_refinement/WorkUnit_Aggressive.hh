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
/// @author Mike Tyka
/// @author Hahnbeom Park: Generalized as a "Sampler" from "Loop Hasher"

#ifndef INCLUDED_protocols_mpi_refinement_WorkUnit_Aggressive_hh
#define INCLUDED_protocols_mpi_refinement_WorkUnit_Aggressive_hh

#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <protocols/wum/WorkUnitBase.hh>
#include <protocols/normalmode/NormalModeRelaxMover.hh>
#include <protocols/normalmode/NormalModeRelaxMover.fwd.hh>

#include <utility/vector1.hh>
#include <set>  // required for std::set

namespace protocols {
namespace mpi_refinement {

////////////////////////////////////////////
//////// WorkUnit Combine
class WorkUnit_CombinePose : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	WorkUnit_CombinePose( core::Size nstruct = 10,
		bool const cartesian = false );

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_CombinePose( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	virtual void init_from_cmd( const core::Size );

	void set_maxfrac( core::Real const setting ){ maxfrac_ = setting; }
	void set_minfrac( core::Real const setting ){ minfrac_ = setting; }

protected:

	void set_nstruct( core::Size const nstruct ){ header.extra_data_1_ = nstruct; }
	void set_cartesian( core::Size const cartesian ){ header.extra_data_2_ = cartesian; }

	core::Size get_nstruct() const { return header.extra_data_1_; }
	core::Size get_cartesian() const { return header.extra_data_2_; }
	void set_defaults();

private:

	core::Real maxfrac_;
	core::Real minfrac_;

};

////////////////////////////////////////////
//////// WorkUnit NormalMode
class WorkUnit_NormalMode : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	WorkUnit_NormalMode( core::Size const nmodes = 0,
		core::Size const nmtype = 0,
		core::Size const relaxtype = 1,
		core::Real const maxscale = 0.0
	);

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_NormalMode( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	virtual void init_from_cmd( const core::Size );

protected:

	void set_nmodes(  core::Size nmodes ){ header.extra_data_1_ = nmodes; }
	void set_nmtype( core::Size relaxtype ){ header.extra_data_2_ = relaxtype; }
	void set_relaxtype( core::Size relaxtype ){ header.extra_data_3_ = relaxtype; }
	void set_maxscale( core::Real setting ){ header.extra_data_4_ = (core::Size)(setting*100.0); }

	core::Size get_nmodes()      const { return header.extra_data_1_; }
	core::Size get_nmtype()    const { return header.extra_data_2_; }
	core::Size get_relaxtype() const { return header.extra_data_3_; }
	core::Real get_maxscale() const { return 0.01*header.extra_data_4_; }
	void set_defaults();

private:

	inline
	protocols::normalmode::NormalModeRelaxMoverOP
	get_NMmover( core::pose::Pose pose,
		core::scoring::ScoreFunctionCOP sfxn_loc,
		core::kinematics::MoveMapOP mm,
		core::Real distcut,
		std::string relaxmode,
		bool const cart )
	{
		if ( cart ) {
			return protocols::normalmode::NormalModeRelaxMoverOP(
				new protocols::normalmode::CartesianNormalModeMover( pose, sfxn_loc, mm, "CA", distcut, relaxmode ) );
		} else {
			return protocols::normalmode::NormalModeRelaxMoverOP(
				new protocols::normalmode::TorsionNormalModeMover( pose, sfxn_loc,  mm, "CA", distcut, relaxmode ) );
		}
	}

};

////////////////////////////////////////////
//////// WorkUnit RamaPerturber
class WorkUnit_RamaPerturber : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_RamaPerturber( core::Size const nsteps = 0,
		core::Size const res1 = 0,
		core::Size const res2 = 0,
		core::Real const kT = 0.0
	);
	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_RamaPerturber( *this ) );
	}

	// Pure virtual functions
	virtual void run();

	virtual void init_from_cmd( const core::Size );

	void set_nsteps( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_res1( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_res2( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_kT( core::Real const setting ){ header.extra_data_4_ = (core::Size)(setting*1000); }

protected:
	core::Size get_nsteps() const { return header.extra_data_1_; }
	core::Size get_res1() const { return header.extra_data_2_; }
	core::Size get_res2() const { return header.extra_data_3_; }
	core::Real get_kT() const { return header.extra_data_4_*0.001; }

	void set_defaults();

private:

};

}
}

#endif

