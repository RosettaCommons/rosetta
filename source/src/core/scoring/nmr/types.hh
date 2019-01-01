// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/scoring/nmr/types.hh
/// @brief   type definitions and enums for classes in core/scoring/nmr
/// @details last Modified: 11/24/18
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_scoring_nmr_types_HH
#define INCLUDED_core_scoring_nmr_types_HH

// Utility headers
#include <core/types.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/HomogeneousTransform.fwd.hh>
#include <utility/fixedsizearray1.hh>

namespace core {
namespace scoring {
namespace nmr {

// Types and Enums
typedef numeric::xyzMatrix<core::Real> Matrix;
typedef numeric::HomogeneousTransform<core::Real> HT;
typedef utility::fixedsizearray1<core::Real,5> Vec5;

/// @brief possible sequences of euler rotation axes
enum EULER_CONVENTION {
	ZXZ_CONVENTION = 1,
	ZYZ_CONVENTION = 2,
	ZYX_CONVENTION = 3
};

/// @brief type of weighting of single PCS/RDC/PRE values
/// @details CONST = constant weighting, 1.0 for every value
///          SIGMA = weighting based on error, 1.0/(error * error)
///          OBSIG = weighting based on observed NMR value and error
///                  val_obs/val_max/ (error * error)
enum SINGLE_NMR_VALUE_WEIGHTING {
	CONST = 1,
	SIGMA = 2,
	OBSIG = 3
};

/// @brief the type how equivalent, ambiguous spins are averaged
enum NMR_VALUE_AVERAGING_TYPE {
	SUM = 1,
	MEAN = 2
};

/// @brief the type of the residual dipolar coupling measured between two nuclei
enum RDC_TYPE {
	RDC_TYPE_NH   = 1,
	RDC_TYPE_NCO  = 2,
	RDC_TYPE_CAHA = 3,
	RDC_TYPE_CACO = 4,
	RDC_TYPE_CACB = 5,
	RDC_TYPE_NCA  = 6,
	RDC_TYPE_CAHN = 7,
	RDC_TYPE_COHN = 8
};

/// @brief type of RDC normalization
///        NORM_TYPE_NH   = relative to N-H RDCs
///        NORM_TYPE_CH   = relative to CA-HA RDCs
///        NORM_TYPE_NONE = it is assumed that input RDCs have been normalized to N-H RDCs in beforehand
///                         (which is the NMR convention). Thus no additional scaling will be applied here
enum RDC_NORM_TYPE {
	NORM_TYPE_NH   = 1,
	NORM_TYPE_CH   = 2,
	NORM_TYPE_NONE = 3
};

/// @brief type of PRE nuclear spin
enum PRE_SPIN_TYPE {
	PRE_SPIN_TYPE_H = 1,
	PRE_SPIN_TYPE_N = 2,
	PRE_SPIN_TYPE_C = 3
};

/// @brief type of the paramagnetic relaxation rate
enum PRE_RATE_TYPE {
	R2_PARA = 1,
	R1_PARA = 2
};

} // nmr
} // scoring
} //core

#endif /* INCLUDED_core_scoring_nmr_types_HH */
