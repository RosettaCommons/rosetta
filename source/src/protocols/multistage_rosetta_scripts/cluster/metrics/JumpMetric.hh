// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.hh
/// @brief  Stores the 6-dimensions of a given jump
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetric_hh
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetric_hh

// unit headers
#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// package headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>
#include <string>
#include <utility/vector1.hh>
#include <core/types.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#include <utility/vector1.srlz.hh>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

class JumpMetric : public ClusterMetric {

public:

	///@brief This enum is not explicitly used but these match the index-numbering I'm using
	enum DOF {
		TRANS_X = 1,
		TRANS_Y,
		TRANS_Z,
		ROT_X,
		ROT_Y,
		ROT_Z
	};

public:
	JumpMetric();

	JumpMetric( utility::vector1< core::Real > const & dofs );

	~JumpMetric() override;

	///@brief just counts number of mutations.
	/// Might someday use BLOSUM62 to create different weights
	platform::Real distance( JumpMetric const & other ) const;

	///@brief just counts number of mutations.
	/// Might someday use BLOSUM62 to create different weights
	platform::Real distance( ClusterMetric const & other ) const override {
		return distance( static_cast< JumpMetric const & > ( other ) );
	}

	/// @brief measure and store the jump dimensions
	void
	analyze( core::pose::Pose const &, unsigned int jump );

public://XML stuff
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag (
		core::pose::Pose const &,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datacache
	) override;

public://getters/setters
	utility::vector1< core::Real > const & dofs() const {
		return dofs_;
	}

private:
	utility::vector1< core::Real > dofs_;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace metrics
} // namespace cluster
} // namespace multistage_rosetta_scripts
} // namespace protocols

#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_metrics_JumpMetric )
#endif // SERIALIZATION


#endif
