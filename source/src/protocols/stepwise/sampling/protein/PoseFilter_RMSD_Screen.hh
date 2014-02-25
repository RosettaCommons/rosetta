// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_stepwise_PoseFilter_RMSD_Screen_HH
#define INCLUDED_protocols_stepwise_PoseFilter_RMSD_Screen_HH

#include <protocols/stepwise/sampling/protein/PoseFilter.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>
#include <map>

//Auto Headers
#include <core/id/AtomID.fwd.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace protein {

	class PoseFilter_RMSD_Screen: public PoseFilter {
	public:

    PoseFilter_RMSD_Screen( utility::vector1< Size > const calc_rms_res,
														core::pose::PoseCOP native_pose,
														core::Real const rmsd_cutoff,
														bool const force_align = false );

		bool passes_filter( core::pose::Pose & pose );

	private:

		void
		initialize_corresponding_atom_id_map( core::pose::Pose const & pose );

		utility::vector1< Size > calc_rms_res_;
		core::pose::PoseCOP native_pose_;
		core::Real const rmsd_cutoff_;
		bool const force_align_;
		bool const cluster_by_all_atom_rmsd_;

		std::map< core::id::AtomID, core::id::AtomID > corresponding_atom_id_map_;

  };

} //protein
} //sampling
} //stepwise
} //protocols

#endif

