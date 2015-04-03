// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/parameters/ParametersSet.hh
/// @brief  Prototypes and method declarations for the ParametersSet class, a class for holding sets of parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_ParametersSet_hh
#define INCLUDED_core_conformation_parametric_ParametersSet_hh


// Unit headers
#include <core/conformation/parametric/ParametersSet.fwd.hh>

// Package headers
#include <core/conformation/parametric/Parameters.fwd.hh>
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/Conformation.fwd.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// Numeric headers

// C++ headers


namespace core {
	namespace conformation {
		namespace parametric {

			/// @brief  ParametersSet class, used to store sets of parameters for parametric backbone generation.
			///
			class ParametersSet : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< ParametersSet >
			{
				public:

					/// @brief constructors
					///
					ParametersSet();

					ParametersSet( ParametersSet const & src );

					~ParametersSet();

					/// @brief Copy this residue( allocate actual memory for it )
					///
					virtual
					ParametersSetOP clone() const;

					/// self pointers
					inline ParametersSetCOP get_self_ptr() const { return shared_from_this(); }
					inline ParametersSetOP get_self_ptr() { return shared_from_this(); }
					inline ParametersSetCAP get_self_weak_ptr() const { return ParametersSetCAP( shared_from_this() ); }
					inline ParametersSetAP get_self_weak_ptr() { return ParametersSetAP( shared_from_this() ); }

					/// @brief Delete all owning pointers in the parameters_ list and reset the list.
					///
					virtual
					void clear_parameters_list() { parameters_.clear(); return; }

					/// @brief Add a Parameters object to the set included in this ParametersSet object.
					///
					virtual
					void add_parameters( ParametersOP new_parameters ) { parameters_.push_back(new_parameters); return; }

					/// @brief Only for copying Conformation objects, this ensures that the new ParametersSet object's
					/// Parameters objects have lists of ResidueOPs that point to residues in the new Conformation object,
					/// rather than to residues that only exist in the Parameters objects.
					void update_residue_links( core::conformation::Conformation &new_conf );

					/// @brief Get the number of Parameters objects associated with this ParametersSet.
					///
					core::Size n_parameters() const { return parameters_.size(); }

					/// @brief Get a Parameters object by index.
					///
					ParametersOP parameters( core::Size const index ) {
						runtime_assert_string_msg( index>0 && index<=parameters_.size(),
							"In core::conformation::parametric::ParametersSet::parameters() : Index out of range.  Expect 0 < index <= number of Parameters objects." );
						return parameters_[index];
					}

					/// @brief Get a Parameters object by index (const-access).
					///
					ParametersCOP parameters( core::Size const index ) const {
						runtime_assert_string_msg( index>0 && index<=parameters_.size(),
							"In core::conformation::parametric::ParametersSet::parameters() : Index out of range.  Expect 0 < index <= number of Parameters objects." );
						return parameters_[index];
					}


				private:

				/********************************************************************************
							PRIVATE DATA
				*********************************************************************************/

				/// @brief List of pointers to parameters objects associated with this ParametersSet.
				/// @details Each Parameters object has parameters describing, for example, one piece of
				/// secondary structure (e.g. one helix in a helical bundle).  Each ParametersSet object
				/// is a collection of Parameters objects, and can describe an assemblage of secondary
				/// structure elements (e.g. a helical bundle).  A Conformation can link to many
				/// ParametersSet objects. 
				utility::vector1 < ParametersOP > parameters_;

			}; //class ParametersSet

		} // namespace parametric
	} // namespace conformation
} // namespace core

#endif
