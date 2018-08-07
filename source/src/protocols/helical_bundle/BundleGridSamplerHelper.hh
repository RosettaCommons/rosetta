// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/helical_bundle/BundleGridSamplerHelper.hh
/// @brief  Prototypes and method declarations for the BundleGridSamplerHelper class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_hh
#define INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_hh


// Unit headers
#include <protocols/helical_bundle/BundleGridSamplerHelper.fwd.hh>
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>

// Package headers
#include <core/conformation/Residue.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers


namespace protocols {
namespace helical_bundle {

/// @brief  BundleGridSamplerHelper class, which stores options for the PerturbBundle mover.
///
class BundleGridSamplerHelper : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< BundleGridSamplerHelper >
{
public:

	/// @brief constructors
	///
	BundleGridSamplerHelper();

	BundleGridSamplerHelper( BundleGridSamplerHelper const & src );

	~BundleGridSamplerHelper() override;

	/// @brief Copy this BundleGridSamplerHelper object( allocate actual memory for it )
	///
	virtual
	BundleGridSamplerHelperOP clone() const;

	/// @brief Reset this BundleGridSamplerHelper object.
	/// @details Clears all internal data.
	void reset();

	/// self pointers
	inline BundleGridSamplerHelperCOP get_self_ptr() const { return shared_from_this(); }
	inline BundleGridSamplerHelperOP get_self_ptr() { return shared_from_this(); }
	inline BundleGridSamplerHelperCAP get_self_weak_ptr() const { return BundleGridSamplerHelperCAP( shared_from_this() ); }
	inline BundleGridSamplerHelperAP get_self_weak_ptr() { return BundleGridSamplerHelperAP( shared_from_this() ); }

public: //Calculators

	/// @brief Perform the pre-calculation that sets up the lists of parameter values to be sampled.
	///
	void initialize_samples();

public: //Getters

	/// @brief Returns the number of different DoFs that the BundleGridSampler will sample over.
	///
	core::Size nDoFs() const { return nDoFs_; }

	/// @brief Returns the helix index of the specified DoF.
	///
	core::Size DoF_helix_index( core::Size const &index ) const {
		runtime_assert_string_msg( index <= allowed_dof_helix_indices_.size(),
			"In BundleGridSamplerHelper::DoF_helix_index(): index out of range!" );
		return allowed_dof_helix_indices_[index];
	}

	/// @brief Returns the DoF type of the specified DoF.
	///
	BPC_Parameters DoF_type( core::Size const &index ) const {
		runtime_assert_string_msg( index <= allowed_dof_types_.size(),
			"In BundleGridSamplerHelper::DoF_type(): index out of range!" );
		return allowed_dof_types_[index];
	}

	/// @brief Returns the current sample value for the specified DoF.
	///
	core::Real DoF_sample_value( core::Size const &index ) const {
		runtime_assert_string_msg( index <= dof_sample_vals_.size(),
			"In BundleGridSamplerHelper::DoF_sample_value(): index out of range!" );
		runtime_assert_string_msg( cur_indices_[index] <= dof_sample_vals_[index].size(),
			"In BundleGridSamplerHelper::DoF_sample_value(): sample index out of range!  Have the samples been initialized properly?" );
		return dof_sample_vals_[index][cur_indices_[index]];
	}


public: //Setters

	/// @brief Add an allowed DoF that the BundleGridSampler will sample over.
	/// @details Must specify the DoF type, the helix index, the number of samples,
	/// and the lower and upper bounds of the range to be sampled.
	void add_DoF(
		BPC_Parameters const doftype,
		core::Size const helix_index,
		core::Size const n_samples,
		core::Real const &lower_val,
		core::Real const &upper_val
	) {
		++nDoFs_; //Increment the number of DoFs
		allowed_dof_types_.push_back(doftype);
		allowed_dof_helix_indices_.push_back(helix_index);
		dof_samples_.push_back(n_samples);
		dof_lower_vals_.push_back(lower_val);
		dof_upper_vals_.push_back(upper_val);
		return;
	}

	/// @brief RECURSIVE function that increments the current sample index.
	/// @details This function adds 1 to the last DoF index.  If the last DoF
	/// index value exceeds the number of samples for that DoF, it resets that
	/// DoF index and increments the second-last DoF index by 1 by recursively
	/// calling itself.  This function is overloaded; the public version by
	/// default tries to increment the last index, while the private version
	/// can be called recursively.
	void increment_cur_indices();

	/// @brief Return the name of a DoF type given its enum.
	///
	inline
	std::string const &
	DoF_name(
		BPC_Parameters const &type
	) const {
		return BundleParametrizationCalculator::parameter_name_from_enum( type );
	}

private:

	/// @brief RECURSIVE function that increments the current sample index.
	/// @details This function adds 1 to the last DoF index.  If the last DoF
	/// index value exceeds the number of samples for that DoF, it resets that
	/// DoF index and increments the second-last DoF index by 1 by recursively
	/// calling itself.
	void increment_cur_indices( core::Size const index_to_increment );

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief Number of DoFs that will be sampled by the BundleGridSampler mover.
	///
	core::Size nDoFs_;

	/// @brief Allowed DoF types.
	/// @details This class stores a list of allowed DoFs that will be sampled by the BundleGridSampler
	/// mover.  It does this in several vectors.  This one lists the DoF type (r0, omega0, etc.).
	utility::vector1 < BPC_Parameters  > allowed_dof_types_;

	/// @brief Allowed DoF helix indices.
	/// @details This class stores a list of allowed DoFs that will be sampled by the BundleGridSampler
	/// mover.  It does this in several vectors.  This one lists the relevant helix index.
	utility::vector1 < core::Size  > allowed_dof_helix_indices_;

	/// @brief Number of samples for each DoF.
	/// @details This class stores a list of allowed DoFs that will be sampled by the BundleGridSampler
	/// mover.  It does this in several vectors.  This one lists the number of samples for each DoF.
	utility::vector1 < core::Size  > dof_samples_;

	/// @brief Lower bound of the range sampled for each DoF.
	/// @details This class stores a list of allowed DoFs that will be sampled by the BundleGridSampler
	/// mover.  It does this in several vectors.  This one lists the lower bound of the range sampled for
	/// each DoF.
	utility::vector1 < core::Real  > dof_lower_vals_;

	/// @brief Upper bound of the range sampled for each DoF.
	/// @details This class stores a list of allowed DoFs that will be sampled by the BundleGridSampler
	/// mover.  It does this in several vectors.  This one lists the upper bound of the range sampled for
	/// each DoF.
	utility::vector1 < core::Real  > dof_upper_vals_;

	/// @brief The vector of vectors of sample values for each DoF.
	/// @details This is pre-computed by the initialize_samples() function.
	utility::vector1 < utility::vector1 < core::Real > > dof_sample_vals_;

	/// @brief The vector indicating the current set of values that we're sampling.
	/// @details Entries in this vector represent indices in the dof_sample_vals_ vector.
	utility::vector1 < core::Size > cur_indices_;

}; //class BundleGridSamplerHelper

} // namespace helical_bundle
} // namespace protocols

#endif //INCLUDED_protocols_helical_bundle_BundleGridSamplerHelper_hh
