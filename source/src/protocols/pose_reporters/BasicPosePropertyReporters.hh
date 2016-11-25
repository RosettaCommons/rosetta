// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/pose_reporters/BasicPosePropertyReporters.hh
/// @brief  Collection of simple pose property reporters
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_pose_reporters_BasicPosePropertyReporters_hh
#define INCLUDED_protocols_pose_reporters_BasicPosePropertyReporters_hh

// Unit Headers
#include <protocols/rosetta_scripts/PosePropertyReporter.hh>

// Project headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <list>

namespace protocols {
namespace pose_reporters {

////////////////////////////////////////////////////////////////////////
// EnergyReporter

/// @brief Report the pose score or a specific energy term

class EnergyReporter : public protocols::rosetta_scripts::PosePropertyReporter {

public:
	EnergyReporter();

	static std::string name() { return "EnergyReporter"; }
	std::string get_name() const override { return name(); }

	core::Real report_property( core::pose::Pose & p ) const override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

protected:

private:
	core::scoring::ScoreFunctionOP scorefxn_;
	core::scoring::ScoreType scoretype_;

};

////////////////////////////////////////////////////////////////////////
// FilterReporter

/// @brief Use a RosettaScripts filter as a reporter

class FilterReporter : public protocols::rosetta_scripts::PosePropertyReporter {

public:
	FilterReporter();

	static std::string name() { return "FilterReporter"; }
	std::string get_name() const override { return name(); }

	core::Real report_property( core::pose::Pose & p ) const override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

protected:

private:
	protocols::filters::FilterOP filter_;

};

////////////////////////////////////////////////////////////////////////
// RMSD reporter

/// @brief Report the RMSD between two poses

class RMSDReporter : public protocols::rosetta_scripts::PosePropertyReporter {

	enum {
		MODE_NONE,
		MODE_CA_rmsd,
		MODE_all_atom_rmsd
	};

public:
	RMSDReporter();

	static std::string name() { return "RMSDReporter"; }
	std::string get_name() const override { return name(); }
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	core::Real report_property( core::pose::Pose & p1, core::pose::Pose & p2 ) const override;


	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	) override;

protected:

private:
	core::Size mode_;
	std::list<Size> residues_;

};

////////////////////////////////////////////////////////////////////////

} // pose_reporters
} // protocols

#endif //INCLUDED_protocols_pose_reporters_BasicPosePropertyReporters_HH
