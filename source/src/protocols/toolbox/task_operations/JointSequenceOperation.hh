// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/JointSequenceOperation.hh
/// @brief set every position to be designable to residues observed in a set of structures
/// @author Rocco Moretti, rmoretti@u.washington.edu

#ifndef INCLUDED_protocols_toolbox_task_operations_JointSequenceOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_JointSequenceOperation_hh

// unit headers
#include <protocols/toolbox/task_operations/JointSequenceOperation.fwd.hh>

//package headers
#include <core/pack/task/operation/TaskOperation.hh>

//project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/sequence/Sequence.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>


#include <core/pack/rotamer_set/UnboundRotamersOperation.fwd.hh>


using namespace core::pack::task;

namespace protocols {
namespace toolbox {
namespace task_operations {

class JointSequenceOperation : public core::pack::task::operation::TaskOperation {
public:

	typedef std::string String;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagCOP TagCOP;

public:

	/// @brief default constructor
	JointSequenceOperation();

	/// @brief destructor
	~JointSequenceOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	void parse_tag( TagCOP tag , DataMap & );

	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	/// @brief Add the sequence from the given filename to the set of allowed aas.
	void add_pdb( std::string filename );

	/// @brief Add the sequence from the given pose to the set of allowed aas.
	void add_pose( Pose const & pose );

	/// @brief Add the sequence from the given filename to the set of allowed aas
	/// and add the rotamers to the set of possible rotamers
	void add_native_pdb( std::string filename );

	/// @brief Add the sequence from the given pose to the set of allowed aas
	/// and add the rotamers to the set of possible rotamers
	void add_native_pose( core::pose::PoseCOP posecop );

	/// @brief Should the current pose (pose supplied to apply) be used in addition to the other ones?
	void use_current_pose( bool ucp );

	/// @brief Should the rotamers for the native poses be used?
	void use_natro( bool unr );

	/// @brief Return which chain to operate on with 0 meaning all
	void set_chain(core::Size chain);

	/// @brief Add the sequence from the given fasta filename to the set of allowed aas
	void add_native_fasta( std::string fasta_file );

	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "JointSequence"; }

private:

	bool use_current_pose_;
	bool use_natro_; // set only with use_natro(), bookkeeping of ubr_
	bool use_fasta_;
	core::pack::rotamer_set::UnboundRotamersOperationOP ubr_;
	std::vector<core::sequence::SequenceOP> sequences_;
	/// @brief Which chain to operate on
	core::Size chain_;
	/// @brief from which chain to extract the sequence information
	core::Size use_chain_;
};


} // TaskOperations
} // toolbox
} // protocols
#endif
