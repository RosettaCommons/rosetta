// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/InteractingRotamerExplosion.hh
/// @brief  TaskOperation class that oversamples rotamers interacting with a certain target
/// @author Florian Richter, florian.richter.1@hu-berlin.de, june '14


#ifndef INCLUDED_protocols_task_operations_InteractingRotamerExplosion_hh
#define INCLUDED_protocols_task_operations_InteractingRotamerExplosion_hh

// Unit Headers
#include <protocols/task_operations/InteractingRotamerExplosion.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh> // abstract base classes

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace task_operations {

/// @brief: taskop to oversample rotamers that make good interactions
///         with a (set of) specified target residue(s)
///         similar to Rosetta 2.x rotamer explosion
class InteractingRotamerExplosion : public core::pack::task::operation::TaskOperation
{

public:
	typedef core::pack::task::operation::TaskOperation parent;

	InteractingRotamerExplosion();
	virtual ~InteractingRotamerExplosion();
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	virtual void apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const;
	virtual void parse_tag( TagCOP tag, DataMap & datamap );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "InteractingRotamerExplosion"; }

private:
	std::string string_target_seqpos_; // this can only be parsed at apply time, when the pose is available;
	core::Real score_cutoff_;
	core::Size ex_level_;
	bool debug_;
	core::Real exclude_radius_;

};

} //namespace protocols
} //namespace task_operations

#endif
