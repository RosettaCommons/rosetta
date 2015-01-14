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

#ifndef INCLUDED_protocols_mpi_refinement_WorkUnit_Relax_hh
#define INCLUDED_protocols_mpi_refinement_WorkUnit_Relax_hh

#include <protocols/mpi_refinement/WorkUnit_Sampler.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace mpi_refinement {

////////////////////////////////////////////
//////// WorkUnit bbGauss
class WorkUnit_bbGauss : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_bbGauss( core::Size const nstruct = 0,
										core::Real const kT = 0.5,
										bool const centroid = false,
										bool const on_defined_segment = false );

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_bbGauss( *this ) );
	}

	// Pure virtual functions
	virtual	void run();

	void set_nstruct( core::Size setting ){ header.extra_data_1_ = setting; }
	void set_kT( core::Real const setting ){ header.extra_data_2_ = (core::Size)(setting*1000); }
	void set_centroid( core::Size setting ){ header.extra_data_3_ = setting; }
	void set_segdef( core::Size setting ){ header.extra_data_4_ = setting; }

protected:

	core::Size get_nstruct(){ return header.extra_data_1_; }
	core::Real get_kT(){ return 0.001*header.extra_data_2_; }
	core::Size get_centroid(){ return header.extra_data_3_; }
	core::Size get_segdef(){ return header.extra_data_4_; }

	void set_defaults();

private:

};

////////////////////////////////////////////
//////// WorkUnit MD
class WorkUnit_MD : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_MD( core::Size const relaxtype = 0,
							 core::Size const scoretype = 0,
							 core::Size const nstruct = 10,
							 core::Real const cstweight = 0.0,
							 bool const looponly = false );

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_MD( *this ) );
	}

	// Pure virtual functions
	virtual	void run();

	void set_relaxtype( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_scoretype( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_nstruct  ( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_cstweight( core::Real const setting ){ header.extra_data_4_ = (core::Size)(setting*100); }

protected:

	core::Size get_relaxtype() const { return header.extra_data_1_; }
	core::Size get_scoretype() const { return header.extra_data_2_; }
	core::Size get_nstruct  () const { return header.extra_data_3_; }
	core::Real get_cstweight() const { return 0.01*header.extra_data_4_; }

	void set_defaults();

private:

};

////////////////////////////////////////////
//////// WorkUnit Relax
class WorkUnit_Relax : public protocols::mpi_refinement::WorkUnit_Sampler
{

public:
	// initialize only via this
	WorkUnit_Relax( 
									core::Size const relaxtype = 0,
									core::Size const scoretype = 0,
									core::Size const nrepeat = 0,
									core::Real const cstweight = 0.0 );

	virtual protocols::wum::WorkUnitBaseOP clone() const {
		return protocols::wum::WorkUnitBaseOP( new WorkUnit_Relax( *this ) );
	}

	// Pure virtual functions
	virtual	void run();

	void set_relaxtype( core::Size const setting ){ header.extra_data_1_ = setting; }
	void set_scoretype( core::Size const setting ){ header.extra_data_2_ = setting; }
	void set_nrepeat  ( core::Size const setting ){ header.extra_data_3_ = setting; }
	void set_cstweight( core::Real const setting ){ header.extra_data_4_ = (core::Size)(setting*100); }

protected:

	core::Size get_relaxtype() const { return header.extra_data_1_; }
	core::Size get_scoretype() const { return header.extra_data_2_; }
	core::Size get_nrepeat()   const { return header.extra_data_3_; }
	core::Real get_cstweight() const { return 0.01*header.extra_data_4_; }

	void set_defaults();

private:
	std::vector< std::string > set_relax_schedule() const;
};

}
}

#endif

