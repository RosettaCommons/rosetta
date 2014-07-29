// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/CentroidDistanceScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_CentroidDistanceScreener_HH
#define INCLUDED_protocols_stepwise_screener_CentroidDistanceScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/CentroidDistanceScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class CentroidDistanceScreener: public StepWiseScreener {

	public:

		//constructor
		CentroidDistanceScreener(  pose::Pose & screening_pose,
															 Size const moving_res,
															 Vector const & reference_centroid,
															 core::Real const max_distance_squared );
		//destructor
		~CentroidDistanceScreener();

	public:

		virtual
		bool check_screen();

		virtual
		std::string
		name() const { return "CentroidDistanceScreener"; }

		virtual
		StepWiseScreenerType
		type() const { return CENTROID_DISTANCE; }

		virtual
		void
		fast_forward( sampler::StepWiseSamplerBaseOP sampler );

	private:

		pose::Pose & screening_pose_;
		Size const moving_res_;
		Vector const & reference_centroid_;
		core::Real const max_distance_squared_;

	};

} //screener
} //stepwise
} //protocols

#endif
