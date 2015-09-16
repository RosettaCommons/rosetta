// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/parameters/Parameters.hh
/// @brief  Prototypes and method declarations for the Parameters class, a class for holding parameters for parametric backbone generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)


#ifndef INCLUDED_core_conformation_parametric_Parameters_hh
#define INCLUDED_core_conformation_parametric_Parameters_hh


// Unit headers
#include <core/conformation/parametric/Parameters.fwd.hh>

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


namespace core {
namespace conformation {
namespace parametric {

/// @brief  Parameters class, used to store sets of parameters for parametric backbone generation.
///
class Parameters : public utility::pointer::ReferenceCount, public utility::pointer::enable_shared_from_this< Parameters >
{
public:

	/// @brief constructors
	///
	Parameters();

	Parameters( Parameters const & src );

	~Parameters();

	/// @brief Copy this residue( allocate actual memory for it )
	///
	virtual
	ParametersOP clone() const;

	/// self pointers
	inline ParametersCOP get_self_ptr() const { return shared_from_this(); }
	inline ParametersOP get_self_ptr() { return shared_from_this(); }
	inline ParametersCAP get_self_weak_ptr() const { return ParametersCAP( shared_from_this() ); }
	inline ParametersAP get_self_weak_ptr() { return ParametersAP( shared_from_this() ); }

	/// @brief Add a residue to the list of residues that these parameters describe.
	///
	virtual
	void add_residue( core::conformation::ResidueOP residue ) { residue_list_.push_back(residue); return; }

	/// @brief Get an owning pointer to a residue in the list of residues that these parameters describe.
	///
	virtual
	core::conformation::ResidueOP residue( core::Size const index ) {
		debug_assert(index <= residue_list_.size());
		return residue_list_[index];
	}

	/// @brief Get a const owning pointer to a residue in the list of residues that these parameters describe.
	///
	virtual
	core::conformation::ResidueCOP residue_cop( core::Size const index ) const {
		debug_assert(index <= residue_list_.size());
		return residue_list_[index];
	}

	/// @brief Get a const owning pointer to the first residue in the list of residues that these parameters describe.
	/// @details  Note that this might not be the first residue in linear sequence, if the residues were put in in non-
	/// sequential order or the residue numbering has changed.
	virtual
	inline core::conformation::ResidueOP first_residue() const {
		debug_assert( residue_list_.size() >= 1);
		return residue_list_[1];
	}

	/// @brief Get a const owning pointer to the last residue in the list of residues that these parameters describe.
	/// @details  Note that this might not be the last residue in linear sequence, if the residues were put in in non-
	/// sequential order or the residue numbering has changed.
	virtual
	inline core::conformation::ResidueOP last_residue() const {
		debug_assert( residue_list_.size() >= 1);
		return residue_list_[residue_list_.size()];
	}

	/// @brief Get the index of the first residue.
	/// @details MUST BE REWRITTEN when the OP issue is resolved.
	virtual
	inline core::Size first_residue_index() const {
		debug_assert( residue_list_.size() >=1 );
		return residue_list_[1]->seqpos();
	}

	/// @brief Get the index of the last residue.
	/// @details MUST BE REWRITTEN when the OP issue is resolved.
	virtual
	inline core::Size last_residue_index() const {
		debug_assert( residue_list_.size() >=1 );
		return residue_list_[residue_list_.size()]->seqpos();
	}

	/// @brief Returns the number (count) of residues that these parameters describe.
	///
	virtual
	inline core::Size n_residue() const {
		return residue_list_.size();
	}

	/// @brief Clears the list of residues that these parameters describe.
	///
	virtual
	inline void reset_residue_list() { residue_list_.clear(); return; }

	/// @brief Assign an element in the residue list to be an owning pointer to an existing residue.
	///
	virtual
	inline void set_residue( core::Size const index, core::conformation::ResidueOP existing_residue ) {
		debug_assert( index > 0 && index <= residue_list_.size() );
		residue_list_[index] = existing_residue;
		return;
	}

private:

	/********************************************************************************
	PRIVATE DATA
	*********************************************************************************/

	/// @brief A list of the residues in this conformation that are described by the parameters in this Parameters object.
	/// @details This is a list of owning pointers so that the residue indices don't mess things up.  (That is, as residue
	/// indices change, these should still point to the proper residues).
	utility::vector1 < core::conformation::ResidueOP > residue_list_;

}; //class Parameters

} // namespace parametric
} // namespace conformation
} // namespace core

#endif
