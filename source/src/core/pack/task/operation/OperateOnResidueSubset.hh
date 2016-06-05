// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/operation/OperateOnResidueSubset.hh
/// @brief  Class, much like the OperateOnCertainResidues task operation, to apply a particular
///         residue-level task operation on the residues identified by a ResidueSelector (Or a pre-selected set of residues).
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com) - Generalized for any vector1< bool >, flip_subset.

#ifndef INCLUDED_core_pack_task_operation_OperateOnResidueSubset_hh
#define INCLUDED_core_pack_task_operation_OperateOnResidueSubset_hh

// Unit Headers
#include <core/pack/task/operation/OperateOnResidueSubset.fwd.hh>

// Package headers
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ Headers

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace task {
namespace operation {

/// @brief  Class, much like the OperateOnCertainResidues task operation, to apply a particular
///         residue-level task operation on the residues identified by a ResidueSelector (Or a pre-selected subset).
class OperateOnResidueSubset : public TaskOperation
{
public:
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;
	typedef pose::Pose Pose;

public:
	OperateOnResidueSubset();
	OperateOnResidueSubset( ResLvlTaskOperationCOP, core::select::residue_selector::ResidueSelectorCOP, bool flip_subset = false );
	OperateOnResidueSubset( ResLvlTaskOperationCOP, utility::vector1< bool > const & subset);

	OperateOnResidueSubset( OperateOnResidueSubset const & );
	OperateOnResidueSubset & operator = ( OperateOnResidueSubset const & );
	virtual ~OperateOnResidueSubset();

	virtual TaskOperationOP clone() const;
	virtual void apply( Pose const &, PackerTask & ) const;

public:

	/// @brief sets the ResLvlTaskOperation that will be applied to residues
	void op( ResLvlTaskOperationCOP );

	/// @brief sets the ResidueSelector that will be used to determine which residues to apply the RLTOP to
	void selector( core::select::residue_selector::ResidueSelectorCOP );

	/// @brief sets a specific subset that will be operated on instead of a selector itself.
	void subset( utility::vector1< bool > const & subset_residues);

	/// @brief Flip the subset on apply.  False becomes True, True becomes False.  Easier to interact with certain RLTOP.
	void flip_subset( bool flip);

public:

	/// @brief Used to parse an xml-like tag to construct the ResLvlTaskOperation and the ResFilter
	virtual void parse_tag( TagCOP, DataMap & );

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	ResLvlTaskOperationCOP op_;
	core::select::residue_selector::ResidueSelectorCOP residue_selector_;
	utility::vector1< bool > user_provided_subset_;
	bool flip_subset_;

};


} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
