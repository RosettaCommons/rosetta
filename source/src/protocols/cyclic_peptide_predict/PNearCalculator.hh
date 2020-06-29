// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/cyclic_peptide_predict/PNearCalculator.hh
/// @brief A class to compute P_Near and DG_folding.
/// @details PNear = sum_i( exp(-(rmsd_i/lambda_i)^2 * exp(-E_i/kbt) ) / sum_j( -exp(-E_j/kbt) ).
/// DG_folding = -kbt*ln( PNear/( 1 - PNear ) )
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_cyclic_peptide_predict_PNearCalculator_hh
#define INCLUDED_protocols_cyclic_peptide_predict_PNearCalculator_hh

#include <protocols/cyclic_peptide_predict/PNearCalculator.fwd.hh>

// Core headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/VirtualBase.hh>

namespace protocols {
namespace cyclic_peptide_predict {

/// @brief A class to compute P_Near and DG_folding.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)
class PNearCalculator : public utility::VirtualBase {

public:

	/// @brief Default constructor (deleted).
	PNearCalculator() = delete;

	/// @brief Options constructor.
	PNearCalculator( core::Real const lambda, core::Real const kbt );

	/// @brief Destructor.
	~PNearCalculator() override;

	/// @brief Clone operation: make a copy of this object, and return an owning pointer to the copy.
	PNearCalculatorOP clone() const;

public:

	/// @brief Add a data point to the set used to compute PNear.
	void add_data_point( core::Real const energy, core::Real const rmsd );

	/// @brief Compute PNear and DeltaG_folding given the data points that have been
	/// added previously using add_data_point().
	/// @details The pnear, Keq, and dgfolding values are overwritten by this operation.  (These
	/// are the outputs.)
	void compute_pnear_and_dgfolding( core::Real &pnear, core::Real &Keq, core::Real &dgfolding ) const;

private:

	/// @brief Accumulator for the numerator.
	core::Real numerator_ = 0.0;

	/// @brief Accumulator for the denominator.
	core::Real denominator_ = 0.0;

	/// @brief Lambda parameter (see PNear equation).
	/// @details Defaults to 0.5.
	core::Real lambda_ = 0.5;

	/// @brief k_B*t parameter (see PNear equation).
	/// @details Defaults to 1.0.
	core::Real kbt_ = 1.0;

};

} //cyclic_peptide_predict
} //protocols

#endif //INCLUDED_protocols_cyclic_peptide_predict_PNearCalculator_hh
