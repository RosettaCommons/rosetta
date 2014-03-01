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
 /// @file protocols/scoring/TensorsOptimizer.hh
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
 /// @last_modified June 2009
 ////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcs_TensorsOptimizer_hh
#define INCLUDED_protocols_scoring_methods_pcs_TensorsOptimizer_hh

// Package headers
#include <protocols/scoring/methods/pcs/PseudocontactShiftData.fwd.hh>
// Project headers
#include <core/optimization/Multifunc.hh>

#include <utility/vector1.hh>


// Utility headers

// Numeric headers

// Objexx headers

// C++ headers


namespace protocols{
namespace scoring{
namespace methods{
namespace pcs{


class TensorsOptimizer : public core::optimization::Multifunc {

public:
	PCS_data const & pcs_d_;

private:
	/// @brief No default constructor: You must provide a PCS_data object when constructing.
	// The unimplemented private constructor turns off automatic generation
  TensorsOptimizer();

public:
  TensorsOptimizer(PCS_data const & pcs_d);

  virtual
  ~TensorsOptimizer();

  // @brief OptE func
  virtual
  core::Real
  operator ()( core::optimization::Multivec const & vars ) const;

  /// @brief OptE dfunc
  virtual
  void
  dfunc(core::optimization::Multivec const & vars,
		core::optimization::Multivec & dE_dvars
	) const;


	//	void
	//	dfunc_test(optimization::Multivec const & vars) const;

private:


}; // TensorsOptimizer



}//namespace pcs
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
