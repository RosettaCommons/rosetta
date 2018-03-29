// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/task_operations/ResidueProbDesignOperation.hh
/// @brief Uses a probability distribution to sample allowed residues at a position, adding those types to the packer -  only if they are designable.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_task_operations_ResidueProbDesignOperation_hh
#define INCLUDED_protocols_task_operations_ResidueProbDesignOperation_hh

#include <protocols/task_operations/ResidueProbDesignOperationCreator.hh>
#include <protocols/task_operations/ResidueProbDesignOperation.fwd.hh>

#include <core/chemical/AA.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <map>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
//#include <core/pose/PDBInfo.hh>


namespace protocols {
namespace task_operations {


/// @brief A TaskOperation that allows amino acids at designable positions through a set of aa probabilities (aa profile).
///
///    Each time a task is generated, it will choose one amino acid from the set for that position and add it to (or replace)
///     the set of packable residues chosen for that position.
///     Decreases number of rotamers for packing and the space for design.  Instead of using energy for profiles, we use selection of residues through
///     probabilities at each task generation.  If picking rounds is higher can result in more than one additional residue in the task from the native.
///
/// @details Default is for the TaskOp to control the allowed amino acids only if the class has probabilities to use for that resnum.
/// If no probabilities were set in general, it will not do anything.
/// Probabilities should add to one, but if not, they will act as weights.
/// Does NOT control which positions are set to design or repack, only what can be designed into.
///
/// If probabilities/weights for a particular residue are not present in the map, will not modify that residues task.
///
class ResidueProbDesignOperation : public core::pack::task::operation::TaskOperation {
	typedef std::map< core::chemical::AA, core::Real > AAProbabilities; //Map of an amino acid and it's probability.
	typedef std::map< core::Size, AAProbabilities > PerResidueAAProbSet; //Amino acid probabilities for a particular residue number.

public:

	ResidueProbDesignOperation();

	//ResidueProbDesignOperation(PerResidueAAProbSet aa_probs);

	ResidueProbDesignOperation(ResidueProbDesignOperation const & src);

	virtual ~ResidueProbDesignOperation();

	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const override;
	void parse_tag( TagCOP tag , DataMap & ) override;
	static void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
	static std::string keyname() { return "ResidueProbDesignOperation"; }


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
	set_aa_probabilities( core::Size const resnum, AAProbabilities aa_probs );

	/// @brief Reads weights for each posenum and residue type. Format for each line in the file is: POSENUM RESIDUETYPE WEIGHT. With weight as a real between 0 and 1
	void
	set_aa_probabilities_from_file( const std::string& weights_file);

	/// @brief Set the probability for a particular aa @ a particular resnum.
	//set_aa_probability( Size resnum, core::chemical::AA const amino_acid, Real probability );

	/// @brief Sets aa probabilities for all positions not given specifically. All designable residues will use these probabilities this if not specified through other functions.
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

	core::Real
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
	/// @details Will result in very different results.
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

	core::pack::task::operation::TaskOperationOP clone() const override;

	void set_defaults();

private:

	void init_for_equal_operator_and_copy_constructor( ResidueProbDesignOperation & lhs, ResidueProbDesignOperation const & rhs);


	/// @brief  Create full weight set for all amino acids from a probability set.
	utility::vector1<core::Real>
	get_weights(AAProbabilities & aa_prob_set) const;

private:

	bool aa_exists(core::chemical::AA amino_acid, AAProbabilities const & aa_prob_set) const;

	PerResidueAAProbSet prob_set_;
	AAProbabilities overall_prob_set_;

	bool include_native_restype_;
	bool keep_task_allowed_aa_;
	core::Real zero_probs_overwrite_;
	core::Size picking_rounds_;

};

}
}
#endif
