// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/pose_selectors/BasicPoseSelectors.hh
/// @brief  Collection of simple pose selectors
/// @author Luki Goldschmidt <lugo@uw.edu>

#ifndef INCLUDED_protocols_pose_selectors_BasicPoseSelectors_hh
#define INCLUDED_protocols_pose_selectors_BasicPoseSelectors_hh

// Unit Headers
#include <protocols/rosetta_scripts/PoseSelector.fwd.hh>
#include <protocols/rosetta_scripts/PoseSelector.hh>
#include <protocols/rosetta_scripts/PosePropertyReporter.fwd.hh>

// Project headers
#include <basic/datacache/DataMap.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/filters/Filter.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <vector>

namespace protocols {
namespace pose_selectors {

/// @brief Logical boolean selector base class

class LogicalSelector : public protocols::rosetta_scripts::PoseSelector {

protected:
	LogicalSelector();
	virtual ~LogicalSelector();

public:
	static std::string name() {	return "LogicalSelector"; }
	virtual std::string get_name() const { return name(); }

	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	utility::vector1<bool> select_poses( utility::vector1< core::pose::PoseOP > poses ) const;

protected:
	virtual inline bool selection_operation( bool a, bool b ) const { return false; }
	virtual inline bool get_default() const { return false; };
	
private:
	std::vector < protocols::rosetta_scripts::PoseSelectorOP > selectors_;

}; // LogicalSelector


/// @brief AND Selector: select poses that were selected by all child selectors

class AndSelector : public LogicalSelector {

public:
	AndSelector() {}
	~AndSelector() {}

	static std::string name() {	return "AndSelector"; }
	std::string get_name() const { return name(); }

protected:
	inline bool selection_operation( bool a, bool b ) const { return a && b; }
	inline bool get_default() const { return true; }

}; // AndSelector


/// @brief OR Selector: select poses that were selected by any child selectors

class OrSelector : public LogicalSelector {

public:
	OrSelector() {}
	~OrSelector() {}

	static std::string name() {	return "OrSelector"; }
	std::string get_name() const { return name(); }

protected:
	inline bool selection_operation( bool a, bool b ) const { return a || b; }
	inline bool get_default() const { return false; }

}; // OrSelector


/// @brief Select top N poses by a specific pose property

class TopNByProperty : public protocols::rosetta_scripts::PoseSelector {

public:
	TopNByProperty();
	~TopNByProperty() {}

	static std::string name() {	return "TopNByProperty"; }
	std::string get_name() const { return name(); }

	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	utility::vector1<bool> select_poses( utility::vector1< core::pose::PoseOP > poses ) const;

protected:
	
private:
	core::Size n_;
	protocols::rosetta_scripts::PosePropertyReporterOP reporter_;
	int order_;

}; // TopNByProperty


/// @brief Use existing RosettaScripts filter as a post selector

class Filter : public protocols::rosetta_scripts::PoseSelector {

public:
	Filter();
	~Filter() {}

	static std::string name() { return "Filter"; }
	std::string get_name() const { return name(); }

	virtual
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose
	);

	utility::vector1<bool> select_poses( utility::vector1< core::pose::PoseOP > poses ) const;

protected:
	
private:
	protocols::filters::FilterOP filter_;

}; // Filter

} // pose_selectors
} // protocols

#endif //INCLUDED_protocols_pose_selectors_BasicPoseSelectors_hh
