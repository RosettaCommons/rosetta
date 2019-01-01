// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/pcs/PCSMultiSet.hh
/// @brief   class that stores and handles data for one tagging site which can contain datasets
///          for multiple lanthanides (i.e. multiple PCSSingleSet objects)
/// @details last Modified: 06/30/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_pcs_PCSMultiSet_HH
#define INCLUDED_core_scoring_nmr_pcs_PCSMultiSet_HH

// Unit headers
#include <core/scoring/nmr/pcs/PCSMultiSet.fwd.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingleSet.fwd.hh>
#include <core/scoring/nmr/pcs/PCSTensor.fwd.hh>
#include <core/scoring/nmr/NMRGridSearch.fwd.hh>
#include <core/scoring/nmr/NMRSpinlabel.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// C++ headers
#include <iostream>
#include <algorithm>

namespace core {
namespace scoring {
namespace nmr {
namespace pcs {

class PCSMultiSet {

public: // Methods

	/// @brief construct from a vector of PCSSingleSets and the
	///        position of the tagging site. Set a default NMRSpinlabel.
	PCSMultiSet(
		utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
		Size const spinlabel_site,
		Real const weight = 1.0
	);

	/// @brief construct from a vector of PCSSingleSets, the position
	///        of the tagging site and an NMRSpinlabel object
	PCSMultiSet(
		utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
		Size const spinlabel_site,
		NMRSpinlabelOP spinlabel_ptr,
		Real const weight = 1.0
	);

	/// @brief construct from a vector of PCSSingleSets, the position
	///        of the tagging site and an NMRGridSearch object
	PCSMultiSet(
		utility::vector1<PCSSingleSetOP> const & pcs_singleset_vec,
		Size const spinlabel_site,
		NMRGridSearchOP gridsearch_ptr,
		Real const weight = 1.0
	);

	/// @brief copy constructor
	PCSMultiSet(PCSMultiSet const & other);

	/// @brief assignment operator
	PCSMultiSet&
	operator=(PCSMultiSet const & rhs);

	/// @brief destructor
	~PCSMultiSet();

	/// @brief compute the score of all PCS datasets given the input pose.
	///        the PCS tensor is solved by SVD or NLS. In the first case,
	///        the metal ion position can be determined by a grid search.
	///        Alternatively, if the spinlabel data member has been set,
	///        the metal ion position can be inferred from the conformational
	///        ensemble of the spinlabel. Here, the metal coordinates are
	///        calculated as the weighted average across the ensemble.
	Real
	compute_score(
		pose::Pose & pose,
		utility::vector1< PCSTensorCOP > & singleset_tensors
	);

	// Getters and Setters
	utility::vector1<PCSSingleSetOP> & get_pcs_singleset_vec() { return pcs_singleset_vec_; }
	utility::vector1<PCSSingleSetOP> const & get_pcs_singleset_vec() const { return pcs_singleset_vec_; }
	Size get_number_metal_ions() const { return number_metal_ions_; }
	Size get_tag_residue_number() const { return tag_residue_number_; }
	Real get_weight() const { return weight_; }
	NMRGridSearchOP get_gridsearch_iterator() { return gridsearch_iterator_; }
	NMRSpinlabelOP get_spinlabel() { return spinlabel_; }
	bool tensors_fixed() const { return fixed_tensors_; }
	Size get_total_number_pcs() const;

	void set_weight(Real weight) { weight_ = weight; }
	void set_gridsearch_iterator(NMRGridSearchOP gridsearch_ptr) { gridsearch_iterator_ = gridsearch_ptr; }
	void set_spinlabel(NMRSpinlabelOP spinlabel_ptr) { spinlabel_ = spinlabel_ptr; }
	void fix_tensors() { fixed_tensors_ = true; }
	void unfix_tensors() { fixed_tensors_ = false; }

	void show(std::ostream & TR) const;

private: // Methods

	/// @brief default constructor
	PCSMultiSet();

	void deep_copy_pcs_singleset_vec(utility::vector1<PCSSingleSetOP> const & other_vec);

	/// @brief creates MTSL spinlabel as default
	void set_default_spinlabel();

private: // Data

	utility::vector1<PCSSingleSetOP> pcs_singleset_vec_;
	Size number_metal_ions_;
	Size tag_residue_number_;
	Real weight_;
	NMRGridSearchOP gridsearch_iterator_;
	NMRSpinlabelOP spinlabel_;
	// Perform PCS calculation from fixed input tensor values
	// and keep tensor values fixed during the score calculation
	bool fixed_tensors_;

};

} // namespace pcs
} // namespace nmr
} // namespace scoring
} // namespace core

#endif // INCLUDED_core_scoring_nmr_pcs_PCSMultiSet_HH
