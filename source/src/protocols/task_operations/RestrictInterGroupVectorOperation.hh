// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/RestrictInterGroupVectorOperation.hh
/// @brief restricts design to only those residues between two groups of structures
/// @author Ben Stranges stranges@unc.edu

#ifndef INCLUDED_protocols_task_operations_RestrictInterGroupVectorOperation_hh
#define INCLUDED_protocols_task_operations_RestrictInterGroupVectorOperation_hh

// unit headers
#include <protocols/task_operations/RestrictInterGroupVectorOperation.fwd.hh>

//package headers
#include <core/pack/task/operation/TaskOperation.hh>

//project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

//C++ headers
#include <set>
#include <utility> //pair


namespace protocols {
namespace task_operations {

class RestrictInterGroupVectorOperation : public core::pack::task::operation::TaskOperation {
public:
	typedef core::pack::task::operation::TaskOperation parent;
	typedef std::set< core::Size > one_group;
	typedef std::pair< one_group, one_group > group_pair;
	typedef utility::vector1< group_pair > group_vector;
public:

	/// @brief default constructor
	RestrictInterGroupVectorOperation();

	/// @brief full constructor
	RestrictInterGroupVectorOperation(
		group_vector const & group,
		core::Real CB_dist_cutoff,
		core::Real nearby_atom_cutoff,
		core::Real vector_angle_cutoff,
		core::Real vector_dist_cutoff);

	//@brief convienience contstuctor for one pair
	RestrictInterGroupVectorOperation(
		group_pair const & one_group,
		core::Real CB_dist_cutoff,
		core::Real nearby_atom_cutoff,
		core::Real vector_angle_cutoff,
		core::Real vector_dist_cutoff);

	/// @brief destructor
	~RestrictInterGroupVectorOperation();

	/// @brief make clone
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictInterGroupVectorOperation"; }

public:
	/// @brief parse_tag
	//parse_tag(
	//  utility::tag::TagCOP tag,
	// core::pose::Pose const & pose);

	/// @brief apply
	virtual void apply( Pose const & pose, core::pack::task::PackerTask & task ) const;

	// /// @brief parse_tag
	// void parse_tag(utility::tag::TagCOP tag);

	/// @brief setters for member data
	void insert_pair( group_pair  pair);
	void CB_dist_cutoff( core::Real  CB_dist_cutoff);
	void nearby_atom_cutoff(core::Real  nearby_atom_cutoff);
	void vector_angle_cutoff(core::Real  vector_angle_cutoff);
	void vector_dist_cutoff(core::Real  vector_dist_cutoff);

private:
	group_vector pair_vector_;
	core::Real CB_dist_cutoff_;
	core::Real nearby_atom_cutoff_;
	core::Real vector_angle_cutoff_;
	core::Real vector_dist_cutoff_;
};


} // TaskOperations
} // protocols
#endif
