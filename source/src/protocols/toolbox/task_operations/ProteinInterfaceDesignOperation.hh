// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh
/// @brief  TaskOperation class that restricts a chain to repacking
/// @author Sarel Fleishman sarelf@uw.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_ProteinInterfaceDesignOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_ProteinInterfaceDesignOperation_hh

// Unit Headers
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictOperationsBase.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <utility/tag/Tag.fwd.hh>

// Utility Headers
#include <core/types.hh>

// C++ Headers
// AUTO-REMOVED #include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace toolbox {
namespace task_operations {

///@details this class is a TaskOperation to prevent repacking of residues not near an interface.
class ProteinInterfaceDesignOperation : public RestrictOperationsBase
{
public:
	typedef RestrictOperationsBase parent;

	ProteinInterfaceDesignOperation();
	void repack_chain1( bool const repack );
	void repack_chain2( bool const repack );
	void design_chain1( bool const design );
	void design_chain2( bool const design );
	///@brief allow all amino acids to be designed at all positions, do not exclude C, G, P
	void allow_all_aas( bool const allow  );
	///@brief allow design of all residues on input pose, do not exclude G,P
	void design_all_aas( bool const design_all  );
	///@brief distance cutoff for atom distance in an interface.
	/// All residues an atoms less than the distance cutoff from an atom in the other chain are
	/// defined as interface.
	void interface_distance_cutoff( core::Real const dist );
	void jump( core::Size const j );
	core::Size jump() const;

	virtual ~ProteinInterfaceDesignOperation();

	virtual TaskOperationOP clone() const;

	virtual
	void
	apply( core::pose::Pose const &, core::pack::task::PackerTask & ) const;

	virtual void parse_tag( TagCOP, DataMap & );

	bool modify_before_jump() const{ return modify_before_jump_; }
	void modify_before_jump( bool const m ) { modify_before_jump_ = m; }
	bool modify_after_jump() const{ return   modify_after_jump_; }
	void modify_after_jump( bool const m ) { modify_after_jump_ = m; }

private:
	bool repack_chain1_, repack_chain2_;
	bool design_chain1_, design_chain2_;
	bool allow_all_aas_, design_all_aas_;
	core::Real interface_distance_cutoff_;
	core::Size jump_;//dflt 1; the jump below which it is chain1, and above which it is chain2
	bool modify_before_jump_, modify_after_jump_; //dflt true,true; set to false to prevent this taskoperation from modifying one of the sides of the jump. This is useful, if, e.g., you want to only restrict the interface on chain2
};

} //namespace protocols
} //namespace toolbox
} //namespace task_operations

#endif // INCLUDED_protocols_toolbox_TaskOperations_ProteinInterfaceDesignOperation_HH
