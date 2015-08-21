// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Kale Kundert (kale.kundert@ucsf.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictToLoops_HH
#define INCLUDED_protocols_toolbox_task_operations_RestrictToLoops_HH

// Unit headers
#include <protocols/toolbox/task_operations/RestrictToLoops.fwd.hh>

// Core headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Protocol headers
#include <protocols/loops/Loops.fwd.hh>

// Utility headers
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

// C++ headers
#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

class RestrictToLoops : public core::pack::task::operation::TaskOperation {

public:
	typedef core::pack::task::operation::TaskOperation parent;

public:

	/// @brief Default constructor.
	RestrictToLoops();

	/// @brief Copy constructor.
	RestrictToLoops( RestrictToLoops const & src );

	/// @brief Default destructor.
	virtual ~RestrictToLoops();

	/// @brief Assignment operator.
	RestrictToLoops & operator= ( RestrictToLoops const & rhs );

	/// @brief Return a deep-copied OP.
	core::pack::task::operation::TaskOperationOP clone() const;

	/// @brief Configure from a RosettaScripts XML tag.
	void parse_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & );

protected:

	/// @brief Help construct instances of this class.
	virtual void init();

	/// @brief Help copy instances of this class.
	virtual void copy( RestrictToLoops & lhs, RestrictToLoops const & rhs );

public:

	/// @brief Apply this operation to the packer task.
	void apply(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask & task ) const;

	/// @brief Return true if design is allowed.
	bool design_loop() const;

	/// @brief Specify whether or not design is allowed.
	void set_design_loop( bool design_loop );

	/// @brief Return the loops allowed to pack.
	loops::LoopsCOP loops() const;

	/// @brief Specify the loops that will be allowed to pack.
	void set_loops( loops::LoopsCOP loops );

	/// @brief Specify the loops that will be allowed to pack.
	void set_loops_from_file( std::string loops_file );

	/// @brief Return true if we are restricting to only design.
	///  AKA RestrictDesignToLoops.
	bool restrict_only_design_to_loops() const;

	/// @brief Specify whether to restrict only design to loops/neighbors
	///  AKA RestrictDesignToLoops.  Does not disable packing for any residue.
	///  Implies and sets design_loop to true.
	void set_restrict_only_design_to_loops( bool restrict_only_design );

protected:

	/// @brief Helper function to prevent code duplication in subclasses.
	void apply_helper(
		core::pose::Pose const & pose,
		core::pack::task::PackerTask & task,
		bool include_neighbors,
		core::Real cutoff_distance,
		bool design_neighbors ) const;

private:
	bool design_loops_;
	bool restrict_only_design_;
	loops::LoopsCOP loops_;

};

}
}
}

#endif
