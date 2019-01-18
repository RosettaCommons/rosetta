// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd3/job_results/PoseJobResult.hh
/// @brief  A standard JobResult that takes and stores a pose.
///   Derived classes may derive additional objects to store from the set pose.
///
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)

#ifndef INCLUDED_protocols_jd3_job_results_PoseJobResult_HH
#define INCLUDED_protocols_jd3_job_results_PoseJobResult_HH

// Unit headers
#include <protocols/jd3/job_results/PoseJobResult.fwd.hh>

// Package headers
#include <protocols/jd3/JobResult.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION


namespace protocols {
namespace jd3 {
namespace job_results {

///@brief
/// A typical PoseJobResult. Stores a pose.
///
///@details
///
/// Can be derived in order give additional stored variables
/// from the pose.  IE - override the pose( pose ) method and grab whatever you need from the pose (sequence, subpose, fragment, etc.)
///
class PoseJobResult : public JobResult
{
public:
	PoseJobResult();
	PoseJobResult( core::pose::PoseOP pose );
	~PoseJobResult() override;

	JobStatus status() const override;

	core::pose::PoseOP pose();
	core::pose::PoseCOP pose() const;

	///@brief Set and Store the pose.  Derive any additional results FROM this pose.
	virtual
	void pose( core::pose::PoseOP setting );

private:

	core::pose::PoseOP pose_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} //job_results
} //jd3
} //protocols


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_jd3_job_results_PoseJobResult )
#endif // SERIALIZATION

#endif //INCLUDED_protocols_jd3_job_results_PoseJobResult
