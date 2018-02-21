// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetricRegistrator.hh
/// @brief  Declaration of the template class for registrating ClusterMetricCreators with
///         the ClusterMetricFactory
/// @author Jack Maguire, jackmaguire1444@gmail.com


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricRegistrator_hh
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricRegistrator_hh

// Package headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricFactory.hh>
#include <utility/factory/WidgetRegistrator.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

/// @brief This templated class will register an instance of an
/// ClusterMetricCreator (class T) with the ClusterMetricFactory.  It will ensure
/// that no ClusterMetricCreator is registered twice, and centralizes
/// this registration logic so that thread safety issues can be handled in
/// one place
template < class T >
class ClusterMetricRegistrator : public utility::factory::WidgetRegistrator< ClusterMetricFactory, T >
{
public:
	typedef utility::factory::WidgetRegistrator< ClusterMetricFactory, T > parent;
public:
	ClusterMetricRegistrator() : parent() {}
};

}
}
}

#endif
