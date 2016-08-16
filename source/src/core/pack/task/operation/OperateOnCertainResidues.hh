// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/operation/OperateOnCertainResidues.hh
/// @brief  class to support the general case of configuring PackerTask at the ResidueLevelTask level in some way, with an optional filter to limit the effects to certain residues
/// @author ashworth

// do not add any derived classes to this file, unless they are generalized abstract base classes and do not actually 'do any work'

#ifndef INCLUDED_core_pack_task_operation_OperateOnCertainResidues_hh
#define INCLUDED_core_pack_task_operation_OperateOnCertainResidues_hh

// Unit Headers
#include <core/pack/task/operation/OperateOnCertainResidues.fwd.hh>

#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/ResLvlTaskOperation.fwd.hh>
#include <core/pack/task/operation/ResFilter.fwd.hh>

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

class OperateOnCertainResidues : public TaskOperation
{
public:
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;
	typedef pose::Pose Pose;
	typedef utility::vector1< Size > ResidueIndices;

public:
	OperateOnCertainResidues();
	OperateOnCertainResidues( ResLvlTaskOperationOP, ResFilterOP );
	OperateOnCertainResidues( OperateOnCertainResidues const & );
	OperateOnCertainResidues & operator = ( OperateOnCertainResidues const & );
	virtual ~OperateOnCertainResidues();

	virtual TaskOperationOP clone() const;
	virtual void apply( Pose const &, PackerTask & ) const;

	/// @brief supports direct limitation of residues to be affected, without the need for a filter
	void residue_indices( ResidueIndices const & );
	ResidueIndices & residue_indices() { return residue_indices_; }
	ResidueIndices const & residue_indices() const { return residue_indices_; }

	/// @brief sets the ResLvlTaskOperation that will be applied to residues
	void op( ResLvlTaskOperationCOP );

	/// @brief sets an optional filter that is applied to each individual residue
	void filter( ResFilterCOP );

	/// @brief Used to parse an xml-like tag to construct the ResLvlTaskOperation and the ResFilter
	virtual void parse_tag( TagCOP, DataMap & );

	static std::string keyname();
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	ResidueIndices residue_indices_;
	ResLvlTaskOperationOP op_;
	ResFilterOP filter_;
};

// do not add any derived classes to this file, unless they are generalized abstract base classes and do not actually 'do any work'

} //namespace operation
} //namespace task
} //namespace pack
} //namespace core

#endif
