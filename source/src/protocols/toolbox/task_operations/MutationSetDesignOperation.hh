// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file /MutationSetDesignOperation.hh
/// @brief Sample a set of mutations each time packer is generated.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_toolbox_task_operations_MutationSetDesignOperation_hh
#define INCLUDED_protocols_toolbox_task_operations_MutationSetDesignOperation_hh

#include <protocols/toolbox/task_operations/MutationSetDesignOperation.fwd.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <map>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace toolbox {
namespace task_operations {


typedef std::map<core::Size, core::chemical::AA> MutationSet;

/// @brief  Sample a set of mutations each time packer is generated.
/// A MutationSet is a simple map of resnum:aa.
/// Each apply will sample a set either at random, or with a set of weights.
///
/// @details Does not enable/disable packing or design by default.  Simply controls what the design set will be.
/// Typically each set would be of the same resnums, but this is not necessary.
/// Iterative Sampling can be achieved somewhat by setting the sample index with each pack.
///  Does not touch any resnums not in a sampled mutation set.
///
class MutationSetDesignOperation : public core::pack::task::operation::TaskOperation {
public:

	MutationSetDesignOperation();
	MutationSetDesignOperation(
		utility::vector1< MutationSet > mutation_sets);

	MutationSetDesignOperation(
		utility::vector1< MutationSet > mutation_sets,
		utility::vector1< core::Real > weights);

	MutationSetDesignOperation(MutationSetDesignOperation const & src);

	virtual
	void
	apply(core::pose::Pose const & pose, core::pack::task::PackerTask & task) const;


	virtual ~MutationSetDesignOperation();

	/// @brief Set the mutation sets.
	/// Each MutationSet is a std::map of resnum:aa
	/// Will use a weight of 1 for each.
	void
	set_mutation_sets( utility::vector1< MutationSet > mutation_sets);

	/// @brief Set the mutation sets and corresponding weights
	/// Each MutationSet is a std::map of resnum:aa
	void
	set_mutation_sets(
		utility::vector1< MutationSet > mutation_sets,
		utility::vector1< core::Real > mutation_set_weights);

	/// @brief add a mutation set with weight 1.
	void
	add_mutation_set( MutationSet mutation_set);

	/// @brief add a mutation set with a weight.
	void
	add_mutation_set( MutationSet mutation_set, core::Real weight);

	/// @brief Clear any stored mutation sets and weights.
	void
	clear_mutation_sets();


	///////// Sampling Control /////////////

	/// @brief Add to the allowed amino acids list instead of replacing them.  Default false.
	void
	add_to_allowed_aas(bool const & setting);

	/// @brief Include native amino acid in the allowed_aas list.  Default False.
	void
	include_native_aa(bool const & setting);

	/// @brief Number of times we sample from our sets. Default 1/apply.
	/// The more rounds, the closer the final amino acid set for each position will be to the full profile for that position.
	/// If using weights for the residue sets, this would increase variability.
	void
	set_picking_rounds(core::Size picking_rounds);


	///////// Iterative Sampling Mode ///////////  This can't work as apply is const. Fuck.

	/// @brief Set the class to iterate through the sets instead of sampling from them.
	//void
	//use_iterative_mode(bool setting);

	/// @brief Get the sampling number if using iterative mode.  Starts at zero, resets with size of mutation set.
	//core::Size
	//get_sample_number() const;

	/// @brief Used to sample particular indexes, since we can't do this iteratively due to const apply. 1 through n
	/// Set to 0 in order to sample from the weights.
	void
	set_sample_index( core::Size sample_index);

	core::Size
	get_sample_index() const;

	/// @brief Get the total number of mutation sets housed.
	core::Size
	get_total_mutation_sets() const;

	/// @brief Set the sample index back to zero - which means we will sample from all of them according to weights.
	void
	reset_sample_index();

public:

	virtual core::pack::task::operation::TaskOperationOP clone() const;

private:
	//void init_for_equal_operator_and_copy_constructor( MutationSetDesignOperation & lhs, MutationSetDesignOperation const & rhs);

	void set_defaults();

private:
	//bool iterative_mode_;
	bool add_to_allowed_aas_;
	bool include_native_aa_;

	utility::vector1< MutationSet > mutation_sets_;
	utility::vector1< core::Real > weights_;
	core::Size picking_rounds_;
	//core::Size sample_number_;

	core::Size sample_index_;


};

} //design
} //antibody
} //protocols


#endif //INCLUDED_protocols_antibody_design_MutationSetDesignOperation_hh



