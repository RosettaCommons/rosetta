// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/nmr/pcs/PCSTensorOptimizer.hh
/// @brief   class that performs optimization of the lanthanide coordinates and other components of the PCSTensor
///          in the laboratory frame. Used in combination with PCS grid search and SVD.
/// @details last Modified: 07/05/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_HH
#define INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_HH

// Unit headers
#include <protocols/nmr/pcs/PCSTensorOptimizer.fwd.hh>

// Package headers
#include <core/scoring/nmr/pcs/PCSSingleSet.fwd.hh>

// Project headers
#include <core/optimization/Multifunc.hh>
#include <core/optimization/types.hh>
#include <core/types.hh>

// Basic headers
#include <basic/Tracer.fwd.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers
#include <numeric/xyzVector.fwd.hh>

// Objexx headers
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1D.hh>

// C++ headers
#include <cmath>

namespace protocols {
namespace nmr {
namespace pcs {

class PCSTensorOptimizer : public core::optimization::Multifunc {

public: // Types

	typedef core::Real Real;
	typedef core::Size Size;
	typedef numeric::xyzVector<core::Real> Vector;
	typedef core::scoring::nmr::pcs::PCSSingleSetOP PCSSingleSetOP;
	typedef core::optimization::Multivec Multivec;

public: // Methods

	/// @brief constructor with a vector of PCSSingleSet pointers as argument
	PCSTensorOptimizer(utility::vector1< PCSSingleSetOP > const & singleset_vec);

	/// @brief destructor
	~PCSTensorOptimizer() override;

	/// @brief error function used in optimization of the PCS tensor parameter
	Real
	operator()(
		Multivec const & tensor_params
	) const override;

	/// @brief gradient function used in optimization of the PCS tensor parameter
	void
	dfunc(
		Multivec const & tensor_params,
		Multivec & dPCS_dparams
	) const override;

private: // Methods

	/// @brief default constructor, should not be called without a valid data object
	PCSTensorOptimizer();

private: // Data

	utility::vector1< PCSSingleSetOP > singleset_vec_;

};

} // namespace pcs
} // namespace nmr
} // namespace protocols

#endif // INCLUDED_protocols_nmr_pcs_PCSTensorOptimizer_HH
