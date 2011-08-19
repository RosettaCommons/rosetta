// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/task_operations/SeqprofConsensus.hh
/// @brief set every position to be designable to residues observed in sequence profile
/// @author Florian Richter, floric@u.washington.edu, april 2011

#ifndef INCLUDED_protocols_toolbox_task_operations_SeqprofConsensusOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_SeqprofConsensusOperation_hh

// unit headers
#include <protocols/toolbox/task_operations/SeqprofConsensusOperation.fwd.hh>

//package headers
#include <core/pack/task/operation/TaskOperation.hh>

//project headers
#include <core/chemical/AA.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/sequence/SequenceProfile.fwd.hh>
#include <core/types.hh>
#include <protocols/ddg/ddGData.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

#include <iostream>

using namespace core::pack::task;

namespace protocols{
namespace toolbox{
namespace task_operations{

class SeqprofConsensusOperation : public core::pack::task::operation::TaskOperation {
public:

	typedef std::string String;
	typedef core::Real Real;
	typedef core::pose::Pose Pose;
	typedef core::pack::task::PackerTask PackerTask;
	typedef core::pack::task::operation::TaskOperation TaskOperation;
	typedef core::pack::task::operation::TaskOperationOP TaskOperationOP;
	typedef TaskOperation parent;
	typedef utility::tag::TagPtr TagPtr;


public:

	/// @brief default constructor
	SeqprofConsensusOperation();

	/// @brief destructor
	 ~SeqprofConsensusOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	void parse_tag( TagPtr tag );

	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	core::sequence::SequenceProfileCOP
	seqprof() const;

	void
	set_seqprof(
		core::sequence::SequenceProfileCOP seqprof );

private:

	std::string seqprof_filename_;
	core::sequence::SequenceProfileCOP seqprof_;
	core::Real min_aa_probability_; // mininum probability that an aa must have in the sequence profile to be considered
	bool prob_larger_current_; //whether probability of a given aa to be included needs to be higher than the probability of the aa in the input pose

};

/// @brief a Task operation that will check whether the amino acid at a
/// position is conserved in the sequence profile and has an unfavorable
/// ddG when mutated to ala. all positions that match this criterion will
/// get set to repacking.
/// @details wt ala positions are set to repacking based on seqprof criterion
/// only.
/// If the input pose contains a forbidden (i.e. non wildtype ) residue
/// at an untouchable position, the residue currently in the pose is
/// also allowed.
class RestrictConservedLowDdgOperation : public SeqprofConsensusOperation {

public:
	typedef SeqprofConsensusOperation Parent;

	RestrictConservedLowDdgOperation();

	~RestrictConservedLowDdgOperation();

	virtual TaskOperationOP clone() const;

	void parse_tag( TagPtr tag );

	virtual void apply( Pose const & pose, PackerTask & task ) const;

	/// @brief returns true if seqpos has a sequence profile
	/// frequency > conservation_cutoff_ and an X->A ddG of >
	/// ddG_cutoff_
	bool
	position_untouchable(
		core::Size seqpos,
		core::chemical::AA seqprof_wt
	) const;

	//convenience function that returns the wild type residue
	//in the pssm file at seqpos
	core::chemical::AA
	seqprof_wt_aa( core::Size seqpos ) const;

	/// @brief convenience function to query
	/// what the ddG is for a to ala mutation
	/// at a certain position
	core::Real
	position_ala_ddG( core::Size seqpos ) const;

	bool
	verbose() const {
		return verbose_;}

private:
	std::string ddG_predictions_filename_;
	std::map< core::Size, ddG::PositionDdGInfoOP > position_ddGs_;
	core::Real conservation_cutoff_; //how freqeunt a residue must be in the sequence profile to count as conserved
	core::Real ddG_cutoff_; //how favorable the ddG at a position has to be for the residue to be considered untouchable
	bool verbose_; //spit out information about untouchable residues

};


} // TaskOperations
} // toolbox
} // protocols
#endif
