// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/matdes/RestrictIdentitiesOperation.hh
/// @brief  TaskOperation class that restricts a vector of Size defined residues to repacking
///   when parsed, it takes in a string and splits by ","
/// @author Neil King (neilking@uw.edu)

#ifndef INCLUDED_protocols_toolbox_task_operations_RestrictIdentitiesOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_RestrictIdentitiesOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/RestrictIdentitiesOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictIdentitiesOperationCreator.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace protocols {
namespace toolbox {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class RestrictIdentitiesOperation : public core::pack::task::operation::TaskOperation {
public:

	RestrictIdentitiesOperation();
	RestrictIdentitiesOperation( utility::vector1 < std::string > identities, bool prevent_repacking );
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	virtual ~RestrictIdentitiesOperation();

	// @brief getters
	utility::vector1< std::string > identities() const;
	bool prevent_repacking() const;

	// @brief setters
	void identities( utility::vector1 < std::string > residues_vec );
	void prevent_repacking( bool const prevent_repacking);

	virtual void apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;
	virtual void parse_tag( TagCOP, DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictIdentities"; }

private:
	std::string unparsed_identities_;
	utility::vector1 < std::string > identities_;
	bool prevent_repacking_;

};

} //namespace task_operations
} //namespace toolbox
} //namespace protocols

#endif // INCLUDED_devel_matdes_RestrictIdentitiesOperation_HH
