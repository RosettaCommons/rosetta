// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/NoRepackDisulfides.hh
/// @brief  prevent disulfides from being repacked; assumes disulfide info in
///         Pose is up-to-date
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_pack_task_operation_NoRepackDisulfides_hh
#define INCLUDED_core_pack_task_operation_NoRepackDisulfides_hh

// unit headers
#include <core/pack/task/operation/NoRepackDisulfides.fwd.hh>

// project headers
#include <core/pack/task/operation/TaskOperation.hh>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {


/// @brief prevent disulfides from being repacked; assume disulfides info in
///  Pose is up-to-date
class NoRepackDisulfides : public core::pack::task::operation::TaskOperation {

private: // typedefs

	typedef core::pack::task::operation::TaskOperation Super;

public: // typedefs

	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef core::pose::Pose Pose;

public: // construct/destruct

	/// @brief default constructor
	NoRepackDisulfides();

	/// @brief copy constructor
	NoRepackDisulfides( NoRepackDisulfides const & rval );

	/// @brief default destructor
	virtual ~NoRepackDisulfides();

public: // virtual constructors

	/// @brief clone this object
	virtual TaskOperationOP clone() const;

public: // methods

	/// @brief apply operations to PackerTask
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname();

};


} // namespace operation
} // namespace task
} // namespace pack
} // namespace core


#endif /* INCLUDED_core_pack_task_operation_NoRepackDisulfides_HH */
