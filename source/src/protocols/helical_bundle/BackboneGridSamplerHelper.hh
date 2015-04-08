// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/helical_bundle/BackboneGridSamplerHelper.hh
/// @brief  A class that stores the various mainchain torsion values that will be sampled by the
/// BackboneGridSampler mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_hh
#define INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_hh


// Unit headers
#include <protocols/helical_bundle/BackboneGridSamplerHelper.fwd.hh>

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

		///@brief  BackboneGridSamplerHelper class, which stores options for the PerturbBundle mover.
		///
		class BackboneGridSamplerHelper : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< BackboneGridSamplerHelper >
		{
			public:

				/// @brief constructors
				///
				BackboneGridSamplerHelper();

				BackboneGridSamplerHelper( BackboneGridSamplerHelper const & src );

				~BackboneGridSamplerHelper();

				/// @brief Copy this BackboneGridSamplerHelper object( allocate actual memory for it )
				///
				virtual
				BackboneGridSamplerHelperOP clone() const;

				/// @brief Reset this BackboneGridSamplerHelper object.
				/// @details Clears all internal data.
				void reset();

				/// self pointers
				inline BackboneGridSamplerHelperCOP get_self_ptr() const { return shared_from_this(); }
				inline BackboneGridSamplerHelperOP get_self_ptr() { return shared_from_this(); }
				inline BackboneGridSamplerHelperCAP get_self_weak_ptr() const { return BackboneGridSamplerHelperCAP( shared_from_this() ); }
				inline BackboneGridSamplerHelperAP get_self_weak_ptr() { return BackboneGridSamplerHelperAP( shared_from_this() ); }

			private: //Calculators

			/// @brief Perform the pre-calculation that sets up the lists of mainchain torsion values to be sampled.
			/// @details Called by the initialize_data() function.
			void initialize_samples();

			public: //Getters

				/// @brief Returns the number of different mainchain torsions that the BackboneGridSampler will sample over.
				///
				core::Size n_torsions() const { return n_torsions_; }
				
				/// @brief Returns the torsion ID for the specified torsion index.
				core::Size torsion_id( core::Size const &index ) const {
					if(index==0 || index > allowed_torsion_indices_.size()) { utility_exit_with_message(
							"In BackboneGridSamplerHelper::torsion_id(): index out of range!" );
					}
					return allowed_torsion_indices_[index];
				}

				/// @brief Returns the current sample value for the specified torsion.
				///
				core::Real torsion_sample_value( core::Size const &index ) const {
					runtime_assert_string_msg( index <= torsion_sample_vals_.size(),
						"In BackboneGridSamplerHelper::torsion_sample_value(): index out of range!" );
					runtime_assert_string_msg( cur_indices_[index] <= torsion_sample_vals_[index].size(),
						"In BackboneGridSamplerHelper::torsion_sample_value(): sample index out of range!  Have the samples been initialized properly?" );
					return torsion_sample_vals_[index][cur_indices_[index]];
				}


			public: //Setters

				/// @brief RECURSIVE function that increments the current sample index.
				/// @details This function adds 1 to the last torsion index.  If the last torsion
				/// index value exceeds the number of samples for that torsion, it resets that
				/// torsion index and increments the second-last torsion index by 1 by recursively
				/// calling itself.  This function is overloaded; the public version by
				/// default tries to increment the last index, while the private version
				/// can be called recursively.
				void increment_cur_indices();

				/// @brief Initialize this object from an object passed from the BackboneGridSampler mover.
				/// 
				void initialize_data(
					utility::vector1 <std::pair <core::Size /*mainchain torsion index*/, std::pair < std::pair < core::Real /*start of range to sample*/, core::Real /*end of range to sample*/ >, core::Size /*samples*/ > > > const &torsions_to_sample
				);

			private:

				/// @brief RECURSIVE function that increments the current sample index.
				/// @details This function adds 1 to the last torsion index.  If the last torsion
				/// index value exceeds the number of samples for that torsion, it resets that
				/// torsion index and increments the second-last torsion index by 1 by recursively
				/// calling itself.
				void increment_cur_indices( core::Size const index_to_increment );

			/********************************************************************************
						PRIVATE DATA
			*********************************************************************************/

			/// @brief Number of mainchain torsions that will be sampled by the BackboneGridSampler mover.
			///
			core::Size n_torsions_;

			/// @brief Allowed mainchain torsion indices.
			/// @details This stores a list of mainchain torsion indices that will be sampled over.
			utility::vector1 < core::Size  > allowed_torsion_indices_;

			/// @brief Number of samples for each mainchain torsion.
			/// @details This class stores a list of allowed mainchain torsions that will be sampled by the BackboneGridSampler
			/// mover.  It does this in several vectors.  This one lists the number of samples for each mainchain torsion.
			utility::vector1 < core::Size  > torsion_samples_;

			/// @brief Lower bound of the range sampled for each mainchain torsion torsion.
			/// @details This class stores a list of allowed mainchain torsion that will be sampled by the BackboneGridSampler
			/// mover.  It does this in several vectors.  This one lists the lower bound of the range sampled for
			/// each mainchain torsion.
			utility::vector1 < core::Real  > torsion_lower_vals_;

			/// @brief Upper bound of the range sampled for each mainchain torsion.
			/// @details This class stores a list of allowed mainchain torsions that will be sampled by the BackboneGridSampler
			/// mover.  It does this in several vectors.  This one lists the upper bound of the range sampled for
			/// each mainchain torsion.
			utility::vector1 < core::Real  > torsion_upper_vals_;

			/// @brief The vector of vectors of sample values for each mainchain torsion.
			/// @details This is pre-computed by the initialize_samples() function.  Order of the outer vector
			/// is the order of mainchain torsion indices in the allowed_torsion_indices_ vector.
			utility::vector1 < utility::vector1 < core::Real > > torsion_sample_vals_;

			/// @brief The vector indicating the current set of values that we're sampling.
			/// @details Entries in this vector represent indices in the torsion_sample_vals_ vector.
			utility::vector1 < core::Size > cur_indices_;

		}; //class BackboneGridSamplerHelper

	} // namespace helical_bundle
} // namespace protocols

#endif //INCLUDED_protocols_helical_bundle_BackboneGridSamplerHelper_hh
