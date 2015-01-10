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

#ifndef INCLUDED_protocols_mpi_refinement_WorkUnit_Sampler_hh
#define INCLUDED_protocols_mpi_refinement_WorkUnit_Sampler_hh

#include <protocols/wum/WorkUnitBase.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace mpi_refinement {

class WorkUnit_Sampler;
typedef utility::pointer::shared_ptr< WorkUnit_Sampler > WorkUnit_SamplerOP;
	typedef utility::pointer::shared_ptr< WorkUnit_Sampler const > WorkUnit_SamplerCOP;

class WorkUnit_Sampler: public protocols::wum::WorkUnit_SilentStructStore {

public:
	// initialize only via this
	WorkUnit_Sampler():	WorkUnit_SilentStructStore(){}

	// @brief Run the workunit - overloaded by children of this class
	virtual void run() = 0;

	virtual void init_from_cmd( const core::Size mpi_rank );

protected:

	core::kinematics::MoveMapOP
	get_movemap( core::pose::Pose const &pose,
							 std::string const mode,
							 bool const nonideal ) const;

	void
	store_to_decoys( core::io::silent::SilentStructCOP start_struct,
									 core::pose::Pose const pose,
									 std::string const additional_tag = "" );

	void
	store_to_decoys( core::io::silent::SilentStructCOP start_struct,
									 core::io::silent::SilentStructOP ss,
									 std::string const additional_tag = "" );

	void
	repack( core::pose::Pose &pose,
					core::scoring::ScoreFunctionOP sfxn );

	void 
	ramp_minpack_loop2( core::pose::Pose &pose,
											utility::vector1< core::Size > const loopres, 
											core::scoring::ScoreFunctionCOP sfxn,
											bool const nonideal = true,
											bool const ramp = true,
											bool const efficient = false,
											core::Real dist_cut = 0.0
											);

	void 
	superimpose_to_ref( core::pose::Pose const &pose_ref,
											core::pose::Pose &pose_work,
											utility::vector1< core::Size > exclude_res 
											= utility::vector1< core::Size >(0) ) const;


	core::scoring::ScoreFunctionOP 
	get_energy( std::string const sfxn_name,
							bool const softpack = false,
							core::Real const weight_coord_cst = 0.0 ) const;

	void
	revert_facts_params() const;

};


}
}

#endif

