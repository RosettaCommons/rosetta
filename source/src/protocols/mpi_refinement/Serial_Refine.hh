// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/mpi_refinement/MPI_Refine.hh
/// @brief
/// @author Mike Tyka
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_Serial_Refine_hh
#define INCLUDED_protocols_mpi_refinement_Serial_Refine_hh

#include <protocols/wum/SilentStructStore.hh>
#include <protocols/mpi_refinement/MultiObjective.hh>
#include <protocols/mpi_refinement/Scheduler.hh>

#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <string>
#include <vector>

namespace protocols {
namespace mpi_refinement {

class Serial_Refine {
public:
	Serial_Refine();
	 ~Serial_Refine();

	void set_defaults();

public:
	core::Real
	apply( core::pose::Pose &pose, 
				 utility::vector1< core::Size > fixres
				 );

	void init();

private:

	void
	load_structures_from_cmdline_into_library(
															core::pose::Pose const & pose,
															protocols::wum::SilentStructStore &library );

	protocols::wum::SilentStructStore 
	perturb( MethodParams const &params,
					 const core::io::silent::SilentStructOP &start_struct );

	void
	dump_structures( protocols::wum::SilentStructStore const &new_structs, 
									 bool score_only,
									 std::string prefix ) const;

	core::pose::Pose
	get_average_structure( protocols::wum::SilentStructStore &decoys,
												 utility::vector1< core::Size > const touse,
												 std::string const columnname, 
												 bool const minimize ) const;

private:

	// library
	protocols::wum::SilentStructStore library_central_;
	protocols::wum::SilentStructStore library_ref_;

	// native
	bool native_given_;
	core::pose::Pose native_pose_;

	// Scheduler
	Scheduler scheduler_;

	// Objfunction
	MultiObjectiveOP fobj_;

	// structure id based on run name
	std::map< std::string, core::Size > ssids_by_name_;

};

} // namespace loops
} // namespace protocols



#endif
