// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictToTerminiOperation.fwd.hh
/// @brief  Restrict to packing only the residues at either or both termini
/// @author Arpit Tandon
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#ifndef INCLUDED_protocols_task_operations_RestrictToTerminiOperation_hh
#define INCLUDED_protocols_task_operations_RestrictToTerminiOperation_hh

// Unit Headers
#include <protocols/task_operations/RestrictToTerminiOperation.fwd.hh>
#include <core/pack/task/operation/TaskOperation.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>


namespace protocols {
namespace task_operations {

class RestrictToTerminiOperation : public core::pack::task::operation::TaskOperation {
public:

	RestrictToTerminiOperation();
	RestrictToTerminiOperation(
		core::Size chain,
		bool restrict_n_terminus,
		bool restrict_c_terminus);

	RestrictToTerminiOperation(
		RestrictToTerminiOperation const & src);

	~RestrictToTerminiOperation();

	core::pack::task::operation::TaskOperationOP
	clone() const;

	void
	apply(
		core::pose::Pose const &,
		core::pack::task::PackerTask & ) const;

	void
	parse_tag( utility::tag::TagCOP, DataMap & );
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "RestrictToTermini"; }

private:
	core::Size chain_;

	//TODO generalize this to alternate definitions of termini
	//(e.g. last k-residues etc.)
	bool repack_n_terminus_;
	bool repack_c_terminus_;


};

} //namespace
} //namespace

#endif // include guard
