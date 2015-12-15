// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/PoseMetricCalculators/SurfaceCalculator.hh
/// @brief A Pose metric which will keep track of the surface energy of a Pose object
/// @author Ron Jacak


#ifndef INCLUDED_protocols_toolbox_pose_metric_calculators_SurfaceCalculator_hh
#define INCLUDED_protocols_toolbox_pose_metric_calculators_SurfaceCalculator_hh

#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <basic/MetricValue.fwd.hh>

#include <utility/vector1.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace toolbox {
namespace pose_metric_calculators {


class SurfaceCalculator : public core::pose::metrics::StructureDependentCalculator {

public:

	SurfaceCalculator( bool remove_nonprotein_res = false );
	~SurfaceCalculator();

public:

	core::pose::metrics::PoseMetricCalculatorOP clone() const { return core::pose::metrics::PoseMetricCalculatorOP( new SurfaceCalculator( remove_nonprotein_res_ ) ); };

protected:

	virtual void lookup( std::string const & key, basic::MetricValueBase* valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::Real total_surface_energy_;
	bool remove_nonprotein_res_;
	utility::vector1< core::Real > residue_surface_energy_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


} // namespace pose_metric_calculators
} // namespace toolbox
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_toolbox_pose_metric_calculators_SurfaceCalculator )
#endif // SERIALIZATION


#endif
