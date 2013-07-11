// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/design/ConservativeDesignOperation.hh
/// @brief TaskOperation to allow only conservative mutations at designable residues
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_design_ConservativeDesignOperation_hh
#define INCLUDED_protocols_antibody_design_ConservativeDesignOperation_hh

#include <protocols/antibody/design/ConservativeDesignOperationCreator.hh>
#include <protocols/antibody/design/ConservativeDesignOperation.fwd.hh>

#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>


namespace protocols {
namespace antibody {
namespace design {
	
	using utility::vector1;
	using core::chemical::AA;

	
///@brief A TaskOperation that sets the allowed amino acids of designable residues to the native amino acid's conservative mutations.
///
///@details Default is to act on all designable residues.  Use limit_to_positions to limit this.
/// Default is to replace the allowed_aas with these conservative mutations. 
///
class ConservativeDesignOperation : public core::pack::task::operation::TaskOperation {

public:
	
	///@brief Default constructor.  Will use native aa from apply (Changes each pack if tf is passed).
	ConservativeDesignOperation();
	
	virtual ~ConservativeDesignOperation();
	
	ConservativeDesignOperation(ConservativeDesignOperation const & src);
	
	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;
	
	
public:
	
	////////////////////////////////////////////////////////////////////////////
	// Class Options
	//
	//
	

	
	///@brief Limit to a subset of residue positions, already set to designable.
	void
	limit_to_positions(vector1< Size > const positions);
	
	///@brief Clear any set positions.
	void
	clear_positions();
	
	///@brief Add to the allowed amino acids list instead of replacing it.  Default false.
	void
	add_to_allowed_aas(bool const & setting);
	
	///@brief Include native amino acid in the allowed_aas list.  Default true.
	void
	include_native_aa(bool const & setting);

	///@brief Will use native residues from the this pose to determine conserved aa.
	void
	use_pose_sequence_as_native(core::pose::Pose const & pose);
	
public:
	ConservativeDesignOperation & operator =( ConservativeDesignOperation const & rhs);
	
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	
private:

	void load_data_from_db();
	
	void load_data_from_file();
	
	void init_for_equal_operator_and_copy_constructor( ConservativeDesignOperation & lhs, ConservativeDesignOperation const & rhs);
	
private:
	
	vector1< vector1<bool> > conserved_mutations_; //AA index in alphabetical order -> conservative mutations.
	vector1< Size > positions_;
	
	bool add_to_allowed_aas_;
	bool include_native_aa_;
	
	std::string pose_sequence_;
	
};
	
	
	
	
	
	
} //task_operations
} //toolbox
} //protocols



#endif	//INCLUDED_protocols_antibody_design_ConservativeDesignOperation_hh

