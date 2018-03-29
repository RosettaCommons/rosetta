// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/DesignInterfacesOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_devel_splice_DesignInterfacesOperation_hh
#define INCLUDED_devel_splice_DesignInterfacesOperation_hh

// Unit Headers
#include <devel/splice/DesignInterfacesOperation.fwd.hh>
#include <protocols/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace devel {
namespace splice {

/// @details this class is a TaskOperation to prevent repacking of residues not near an interface.
class DesignInterfacesOperation : public protocols::task_operations::RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	DesignInterfacesOperation();

	void design_shell( core::Real const radius );
	core::Real design_shell() const{ return design_shell_; }

	~DesignInterfacesOperation() override;

	TaskOperationOP clone() const override;


	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const override;

	void parse_tag( TagCOP, DataMap & ) override;
	static std::string keyname() { return "DesignInterfacesOperation"; }
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	void repack_shell( core::Real const repack_shell) {
		repack_shell_ = repack_shell;
		if ( repack_shell <= design_shell() ) {
			design_shell_ = repack_shell;
		}
	}
	core::Real repack_shell() const{ return repack_shell_; }

	bool restrict_to_repacking_chain1() const { return restrict_to_repacking_chain1_;}
	void restrict_to_repacking_chain1( bool const c ){ restrict_to_repacking_chain1_ = c; }

	bool restrict_to_repacking_chain2() const { return restrict_to_repacking_chain2_;}
	void restrict_to_repacking_chain2( bool const c ){ restrict_to_repacking_chain2_ = c; }
private:
	core::Real design_shell_, repack_shell_; //dflt 6 and 8
	bool restrict_to_repacking_chain1_, restrict_to_repacking_chain2_; //dflt false, true
};

} //namespace splice
} //namespace devel

#endif // INCLUDED_protocols_toolbox_TaskOperations_DesignInterfacesOperation_HH
