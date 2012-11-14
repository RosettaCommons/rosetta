// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author 

#ifndef INCLUDED_protocols_hotspot_hashing_filters_PlaceMinimizeSearch_hh
#define INCLUDED_protocols_hotspot_hashing_filters_PlaceMinimizeSearch_hh


#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

#include <protocols/moves/Mover.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>

#include <core/kinematics/Stub.hh>
#include <protocols/hotspot_hashing/SearchPattern.hh>

namespace protocols {
namespace hotspot_hashing {

class PlaceMinimizeSearch : public utility::pointer::ReferenceCount
{
	public:
    PlaceMinimizeSearch(
			core::pose::Pose const & target_pose,
			core::conformation::ResidueCOP target_residue,
			SearchPatternOP search_pattern,
			protocols::moves::MoverOP relax_mover,
			protocols::filters::FilterOP triage_filter,
			std::string output_tag);

		PlaceMinimizeSearch(
			core::pose::Pose target_pose,
			core::conformation::ResidueCOP target_residue,
			SearchPatternOP search_pattern,
			protocols::moves::MoverOP relax_mover,
			protocols::filters::FilterOP triage_filter,
			std::string output_tag,
			protocols::jd2::JobOP current_job,
			protocols::jd2::JobOutputterOP current_job_outputter);
		
		void execute();

		friend class HotspotHashingTests;

		static void placeResidueAtTransform( core::pose::Pose & pose, core::conformation::Residue const & residue, core::kinematics::Stub transform, core::Size & residuejumpindex, core::Size & residueindex );

		static core::kinematics::Stub residueStubCentroidTransform(core::conformation::Residue const & residue);
		static Vector residueStubCentroid(core::conformation::Residue const & residue);
		static void placeResidueOnPose(core::pose::Pose & pose, core::conformation::Residue const & residue);

	protected:
		void logPreMinPose(core::pose::Pose & pose, core::Size transformindex, core::Size residuejumpindex, core::Size residueindex);
		void logPostMinPose(core::pose::Pose & pose, core::Size transformindex, core::Size residuejumpindex, core::Size residueindex);

  private:
		core::pose::Pose target_pose_;
		core::conformation::ResidueCOP target_residue_;

		protocols::hotspot_hashing::SearchPatternOP search_pattern_;

    protocols::moves::MoverOP relax_mover_;
    protocols::filters::FilterOP triage_filter_;

		std::string output_tag_;

		protocols::jd2::JobOP current_job_;
		protocols::jd2::JobOutputterOP current_job_outputter_;
};

}
}

#endif 
