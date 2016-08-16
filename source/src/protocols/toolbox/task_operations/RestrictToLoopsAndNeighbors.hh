// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.hh
/// @brief  This class allows the selection for packing (and possibly design) of residues contained in a Loops object as
///         well as the neighbors within a specified cutoff distance, with the default and maximum cutoff at 10.0 A.
/// @author Brian D. Weitzner (brian.weitzner@gmail.com)

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToLoopsAndNeighbors_HH
#define INCLUDED_protocols_toolbox_task_operations_RestrictToLoopsAndNeighbors_HH

// Unit headers
#include <protocols/toolbox/task_operations/RestrictToLoops.hh>
#include <protocols/toolbox/task_operations/RestrictToLoopsAndNeighbors.fwd.hh>

// Package headers
#include <protocols/loops/Loops.fwd.hh>

// Project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @brief  This class allows the selection for packing (and possibly design) of residues contained in a Loops object as
///         well as the neighbors within a specified cutoff distance, with the default and maximum cutoff at 10.0 A.
class RestrictToLoopsAndNeighbors : public RestrictToLoops
{
public:
	typedef RestrictToLoops parent;

public:
	/// @brief Default constructor.
	RestrictToLoopsAndNeighbors();

	/// @brief Copy constructor.
	RestrictToLoopsAndNeighbors(RestrictToLoopsAndNeighbors const & src);

	/// @brief Default destructor.
	virtual ~RestrictToLoopsAndNeighbors();

	/// @brief Assignment operator.
	RestrictToLoopsAndNeighbors & operator=( RestrictToLoopsAndNeighbors const & rhs );

	/// @copydoc RestrictToLoops::clone
	virtual core::pack::task::operation::TaskOperationOP clone() const;

	/// @copydoc RestrictToLoops::parse_tag
	void parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToLoopsAndNeighbors"; }

	bool include_neighbors() const;
	void set_include_neighbors( bool include_neighbors );

	bool design_neighbors() const;
	void set_design_neighbors( bool design_neighbors );

	core::Real cutoff_distance() const;
	void set_cutoff_distance( core::Real cutoff_distance );

private:
	void init();
	void init_for_copy( RestrictToLoopsAndNeighbors & lhs, RestrictToLoopsAndNeighbors const & rhs);

private:
	bool include_neighbors_;
	bool design_neighbors_;
	core::Real cutoff_distance_;

};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif // INCLUDED_protocols_toolbox_task_operations_RestrictToLoopsAndNeighbors_HH
