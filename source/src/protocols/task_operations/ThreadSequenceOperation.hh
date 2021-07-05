// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/ThreadSequenceOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_task_operations_ThreadSequenceOperation_hh
#define INCLUDED_protocols_task_operations_ThreadSequenceOperation_hh

// Unit Headers
#include <protocols/task_operations/ThreadSequenceOperation.fwd.hh>
#include <protocols/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace task_operations {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class ThreadSequenceOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	ThreadSequenceOperation();
	ThreadSequenceOperation( std::string const & seq );

	~ThreadSequenceOperation() override;

	TaskOperationOP clone() const override;

	void target_sequence( std::string const & seq );
	std::string target_sequence() const;

	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const override;
	static void restrict_aas(core::Size const resi, char const aa, core::pose::Pose const & pose, core::pack::task::PackerTask & task);

	void parse_tag( TagCOP, DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ThreadSequence"; }

	core::Size start_res() const;
	void start_res( core::Size const s );
	bool allow_design_around() const{ return allow_design_around_;}
	void allow_design_around( bool const b ){ allow_design_around_ = b ; }
	core::Size chain_num() const{ return chain_num_; }
	void chain_num( core::Size const s ) { chain_num_ = s;}
	bool filter_non_aas() const{ return filter_non_aas_;}
	void filter_non_aas( bool const b ){ filter_non_aas_ = b ; }

private:
	std::string target_sequence_;
	core::Size start_res_; // dflt 1; which residue number to start threading (useful for partial threads)
	bool allow_design_around_; //dflt true; if false restricts rest of the pose to repakcing
	core::Size chain_num_; //dflt 0; which chain to start threading on. 0 equals entire pose
	bool filter_non_aas_; //dflt false; if true restricts threading to amino acids only
};

} //namespace protocols
} //namespace task_operations

#endif // INCLUDED_protocols_TaskOperations_ThreadSequenceOperation_HH
