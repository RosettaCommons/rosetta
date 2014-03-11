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
 /// @authorsv Christophe Schmitz //kalabharath
 ///
 /// @last_modified Aug 2011
 ////////////////////////////////////////////////


#ifndef INCLUDED_protocols_scoring_methods_pcsTs9_TensorsOptimizer_hh
#define INCLUDED_protocols_scoring_methods_pcsTs9_TensorsOptimizer_hh

// Package headers
#include <protocols/scoring/methods/pcsTs9/PseudocontactShiftData.fwd.hh>
// Project headers
#include <core/optimization/Multifunc.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers


namespace protocols{
namespace scoring{
namespace methods{
namespace pcsTs9{


class TensorsOptimizer_Ts9 : public core::optimization::Multifunc {

public:
	PCS_data_Ts9 const & pcs_d_;

  TensorsOptimizer_Ts9();

  TensorsOptimizer_Ts9(PCS_data_Ts9 const & pcs_d);

  virtual
  ~TensorsOptimizer_Ts9();

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


}; // TensorsOptimizer_Ts9



}//namespace pcsTs9
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
