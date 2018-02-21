// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetric_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetric_HH

#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

#include <platform/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

class ClusterMetric : public utility::pointer::ReferenceCount
{
public:
	ClusterMetric();
	~ClusterMetric();

	///@brief Measures the distance between this metric and another of the same type.
	/// You can assume that other has the same type as this
	virtual
	platform::Real
	distance( ClusterMetric const & other ) const = 0;

	virtual
	void
	parse_my_tag (
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) = 0;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // cluster
} // multistage_rosetta_scripts
} // protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_ClusterMetric )
#endif // SERIALIZATION


#endif
