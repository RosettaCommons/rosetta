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
#include <core/io/ddg/PositionDdGInfo.fwd.hh>
#include <protocols/toolbox/task_operations/RestrictToAlignedSegments.fwd.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// AUTO-REMOVED #include <iostream>

#include <map>


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
	typedef utility::tag::TagCOP TagCOP;


public:

	/// @brief default constructor
	SeqprofConsensusOperation();

	/// @brief destructor
	 ~SeqprofConsensusOperation();

	/// @brief make clone
	virtual TaskOperationOP clone() const;

public:

	void parse_tag( TagCOP tag , DataMap & );

	/// @brief apply
	virtual void apply( Pose const & pose, PackerTask & task ) const;

	core::sequence::SequenceProfileCOP
	seqprof() const;

	void
	set_ignore_pose_profile_length_mismatch( bool const setting ) {
		ignore_pose_profile_length_mismatch_ = setting; }

	/// @brief Set the sequence profile. If reweight is true, convert the profile into per-residue probabilities first
	void
	set_seqprof( core::sequence::SequenceProfileOP seqprof, bool reweight = false );
	void convert_scores_to_probabilities( bool const c ){ convert_scores_to_probabilities_ = c; }
	bool convert_scores_to_probabilities() const{ return convert_scores_to_probabilities_;}
	RestrictToAlignedSegmentsOperationOP restrict_to_aligned_segments() const;
	void restrict_to_aligned_segments( RestrictToAlignedSegmentsOperationOP rtas );
	core::Real conservation_cutoff_aligned_segments() const { return conservation_cutoff_aligned_segments_; }
	void conservation_cutoff_aligned_segments( core::Real const c ) { conservation_cutoff_aligned_segments_ = c; }
	ProteinInterfaceDesignOperationOP protein_interface_design() const;
	void protein_interface_design( ProteinInterfaceDesignOperationOP pido );
	core::Real conservation_cutoff_protein_interface_design() const{ return conservation_cutoff_protein_interface_design_; }
	void conservation_cutoff_protein_interface_design( core::Real const c ){ conservation_cutoff_protein_interface_design_ = c; }
	void debug( bool const b ){ debug_ = b; }
	bool debug() const{ return debug_; }
private:

	std::string seqprof_filename_;
	/// @details Stored as a per-position probability weighted value
	core::sequence::SequenceProfileOP seqprof_;
	/// @brief mininum probability that an aa must have in the sequence profile to be considered
	core::Real min_aa_probability_;
	bool keep_native_;
	/// @brief whether probability of a given aa to be included needs to be higher than the probability of the aa in the input pose
	bool prob_larger_current_;

	/// @default false. if true, every pose seqpos that is bigger
	/// than the size of the sequence_profile will be set to repacking
	/// allows using this taskop in situations where one wants to
	/// do consensus design on one chain of a protein/protein interface
	bool ignore_pose_profile_length_mismatch_;
	bool convert_scores_to_probabilities_; //dflt 1; after reading the pssm, should we convert the scores to probabilities?

	RestrictToAlignedSegmentsOperationOP restrict_to_aligned_segments_; // dflt NULL; this is used to define which residues are considered to be aligned for conservation_cutoff_aligned_segments_
	core::Real conservation_cutoff_aligned_segments_; // dflt -100000;

	ProteinInterfaceDesignOperationOP protein_interface_design_; //dflt NULL; this is used to define which residues are considered to be interface, for conservation_cutoff_interface_design_
	core::Real conservation_cutoff_protein_interface_design_; // dflt -100000;
	core::Size chain_num_; //dflt set to 1

	bool debug_; // dflt false; if true be more chatty
	bool restrict_to_repacking_; // dflt true; if set to false then all the protein can go through design and not just residues with seqeunce constratins. Gideon Oct14
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

	void parse_tag( TagCOP tag , DataMap & );

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

	bool verbose() const { return verbose_;}

private:
	std::string ddG_predictions_filename_;
	std::map< core::Size, core::io::PositionDdGInfo::PositionDdGInfoOP > position_ddGs_;
	core::Real conservation_cutoff_; //how freqeunt a residue must be in the sequence profile to count as conserved
	core::Real ddG_cutoff_; //how favorable the ddG at a position has to be for the residue to be considered untouchable
	bool verbose_; //spit out information about untouchable residues
};


} // TaskOperations
} // toolbox
} // protocols
#endif
