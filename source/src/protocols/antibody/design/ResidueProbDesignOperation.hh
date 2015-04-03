// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/design/ResidueProbDesignOperation.hh
/// @brief Uses a probability distribution to sample allowed residues at a position, adding those types to the packer -  only if they are designable. 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_antibody_design_ResidueProbDesignOperation_hh
#define INCLUDED_protocols_antibody_design_ResidueProbDesignOperation_hh

#include <protocols/antibody/design/ResidueProbDesignOperationCreator.hh>
#include <protocols/antibody/design/ResidueProbDesignOperation.fwd.hh>

#include <core/chemical/AA.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <map>
//#include <core/pose/PDBInfo.hh>


//Will move to protocols/toolbox after thorough testing.
namespace protocols {
namespace antibody {
namespace design {
	using core::Real;
	using core::Size;
	using utility::vector1;
	
/// @brief A TaskOperation that allows amino acids at designable positions through a set of aa probabilities.
///
/// @details Default is for the TaskOp to control the allowed amino acids only if the class has probabilities to use for that resnum. 
/// If no probabilities were set, it will not do anything.
/// Probabilities should add to one, but if not, they will act as weights.
///
class ResidueProbDesignOperation : public core::pack::task::operation::TaskOperation {
	typedef std::map< core::chemical::AA, Real > AAProbabilities; //Map of an amino acid and it's probability.
	typedef std::map< Size, AAProbabilities > PerResidueAAProbSet; //Amino acid probabilities for a particular residue number.
	
public:
	
	ResidueProbDesignOperation();
	
	//ResidueProbDesignOperation(PerResidueAAProbSet aa_probs);
	
	ResidueProbDesignOperation(ResidueProbDesignOperation const & src);
	
	virtual ~ResidueProbDesignOperation();
	
	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;
	

public:
	
	////////////////////////////////////////////////////////////////////////////
	// Probability Settings
	//
	//
	
	/// @brief Sets the whole aa probability set.
	void
	set_aa_probability_set( PerResidueAAProbSet per_res_aa_probs);
	
	/// @brief Sets aa probabilities for a particular resnum
	void
	set_aa_probabilities( Size const resnum, AAProbabilities aa_probs );

	/// @brief Set the probability for a particular aa @ a particular resnum.  
	//set_aa_probability( Size resnum, core::chemical::AA const amino_acid, Real probability );

	/// @brief Sets aa probabilities for all positions not given specifically. All designable residues will use this if not specified through other functions.
	void
	set_overall_aa_probabilities( AAProbabilities aa_probs );
	
	
public:

	////////////////////////////////////////////////////////////////////////////
	// Reset Probabilities
	//
	//
	
	/// @brief Clears the probability matrix
	void
	clear_prob_set();
	
	/// @brief Clears the overall probabilities if any are set.
	void
	clear_overall_prob_set();
	
	
public:
	
	////////////////////////////////////////////////////////////////////////////
	// General Options
	//
	//
	
	/// @brief Include current residuetype in packer in addition to selected mutation. Default is true.
	void
	set_include_native_restype(bool include);
	
	bool
	include_native_restype() const;
	
	
	/// @brief Sample any zero probabilities at this probability instead of not using them.
	void
	set_sample_zero_probs_at(core::Real const probability);
	
	Real
	sample_zero_probs_at() const;
	
	
	/// @brief Max number of times to sample from the probability distribution.  Default is once. 
	/// @details This is to increase variability.
	/// Too many rounds will even out the distribution.
	/// Can also run the PackRotamers mover or equivalent multiple times for sampling.  
	void
	set_picking_rounds(core::Size picking_rounds);
	
	Size
	picking_rounds() const;
	
	
	/// @brief Add the aa types chosen from the probabilities to the  allowed aa at the position instead of overwriting it.
	/// @details  Will result in very different results.  Used one wants to use this Op to sample other residues in addition within a specific probability.  
	void
	set_keep_task_allowed_aas(bool setting);
	
	bool
	keep_task_allowed_aas() const;
	
	
	bool
	resnum_exists_in_set(core::Size const resnum) const;
	
public:

	////////////////////////////////////////////////////////////////////////////
	// Secondary Structure Probabilities + Options
	//
	//
	
	
public:
	//ResidueProbDesignOperation & operator =( ResidueProbDesignOperation const & rhs);
	
	virtual core::pack::task::operation::TaskOperationOP clone() const;
	
	//virtual void parse_tag( TagCOP, DataMap & );  //Not implemented - Not currently a general way to load data from a file or db.  
	
private:
	void init();
	void init_for_equal_operator_and_copy_constructor( ResidueProbDesignOperation & lhs, ResidueProbDesignOperation const & rhs);
	
	
	/// @brief  Create full weight set for all amino acids from a probability set.
	vector1<Real>
	get_weights(AAProbabilities & aa_prob_set) const;
	
	bool aa_exists(core::chemical::AA amino_acid, AAProbabilities const & aa_prob_set) const;
	
	PerResidueAAProbSet prob_set_;
	AAProbabilities overall_prob_set_;
	
	bool include_native_restype_;
	bool keep_task_allowed_aa_;
	Real zero_probs_overwrite_;
	Size picking_rounds_;
	
};

}
}
}
#endif
