// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

 //////////////////////////////////////////////
 /// @begin
 ///
 /// @file protocols/scoring/methods/pcs2/TensorsOptimizerFix.hh
 ///
 /// @brief
 ///
 /// @detailed
 ///
 /// @param
 ///
 /// @return
 ///
 /// @remarks
 ///
 /// @references
 ///
 /// @authorsv Christophe Schmitz
 ///
 /// @last_modified February 2010
 ////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_TensorsOptimizerFix_hh
#define INCLUDED_protocols_scoring_methods_pcs2_TensorsOptimizerFix_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsDataCenter.fwd.hh>
// Project headers
#include <core/optimization/Multifunc.hh>

// Utility headers
#include <utility/vector1.hh>

// Numeric headers

// Objexx headers

// C++ headers


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

class TensorsOptimizerFix : public core::optimization::Multifunc {

public:
	PcsDataCenter const & pcs_d_c_;
	/*
	core::Real const xM_;
	core::Real const yM_;
	core::Real const zM_;
	*/

	/*
	utility::vector1< core::Real > Xxx_coef_vect_;
	utility::vector1< core::Real > Xxy_coef_vect_;
	utility::vector1< core::Real > Xxz_coef_vect_;
	utility::vector1< core::Real > Xyy_coef_vect_;
	utility::vector1< core::Real > Xyz_coef_vect_;
	*/

private:

	/// @brief No default constructor: You must provide a PcsDataCenter object when initializing.
	// The unimplemented private default constructor inhibits autogeneration of a default constructor.
  TensorsOptimizerFix();

public:

  TensorsOptimizerFix(PcsDataCenter const & pcs_d_c/*,
											core::Real xM,
											core::Real yM,
											core::Real zM*/);

  virtual
  ~TensorsOptimizerFix();

  // @brief OptE func
  virtual
  core::Real
  operator ()( core::optimization::Multivec const & vars ) const;

  core::Real
  func( core::optimization::Multivec const & vars ) const;


  /// @brief OptE dfunc
  virtual
  void
  dfunc(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
	) const;

	/// @brief exact derivative (fast)
  void
  dfunc_exact(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
	) const;

	/// @brief numeric derivative (slow)
  void
  dfunc_numeric(core::optimization::Multivec const & vars,
	core::optimization::Multivec & dE_dvars
	) const;

private:


};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
