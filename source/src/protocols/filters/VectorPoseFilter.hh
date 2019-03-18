// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/filters/VectorPoseFilter.hh
/// @brief Designates a filter that can be passed multiple poses by the VectorPoseJobDistributor
/// Any filters deriving from this subclass can then act on all of the input poses simultaneously
/// Only accessible through recon application.
/// @author Alex Sevy (alex.sevy@gmail.com)

#ifndef INCLUDED_PROTOCOLS_FILTERS_VectorPoseFilter_HH
#define INCLUDED_PROTOCOLS_FILTERS_VectorPoseFilter_HH

#include <protocols/filters/VectorPoseFilter.fwd.hh>
#include <protocols/filters/Filter.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace filters {

/// @brief Designates a filter that can be passed multiple poses by the VectorPoseJobDistributor
/// Any filters deriving from this subclass can then act on all of the input poses simultaneously
/// Only accessible through recon application.
class VectorPoseFilter : public Filter {

public:
	/// @brief Constructor
	VectorPoseFilter();

	virtual ~VectorPoseFilter();

	VectorPoseFilter( std::string const& name );

	VectorPoseFilter( VectorPoseFilter const & other );

	virtual bool apply ( core::pose::Pose const & pose ) const = 0;

	/// @brief Apply function to run under MPI
	virtual bool apply_mpi ( core::pose::Pose const & pose ) const = 0;

	/// @brief Set the vector of poses for the mover to act upon
	void set_poses( utility::vector1< core::pose::PoseOP > const & poses );

protected:
	utility::vector1< core::pose::PoseOP > poses_;
};

}  // namespace filters
}  // namespace protocols


#endif  // PROTOCOLS_FILTERS_VectorPoseFilter_HH_
