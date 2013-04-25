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
 /// @file protocols/scoring/methods/pcs2/PcsDataCenterManagerSingleton.hh
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

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenterManagerSingleton_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenterManagerSingleton_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>
#include <protocols/scoring/methods/pcs2/PcsEnergyParameterManager.fwd.hh>

// Project headers
#include <core/types.hh>

#include <utility/vector1.hh>


// Utility headers

// Numeric headers

// ObjexxFCL headers

// c++ headers

namespace protocols{
namespace scoring{
namespace methods{
namespace pcs2{

class PcsDataCenterManagerSingleton {
private:

  utility::vector1<PcsDataCenter> PCS_data_all_;
	static PcsDataCenterManagerSingleton * instance_;

public:
	static PcsDataCenterManagerSingleton *
	get_instance(PcsEnergyParameterManager & pcs_e_p_m);

  PcsDataCenterManagerSingleton(PcsEnergyParameterManager & pcs_e_p_m); //construct


	/// @brief Give me the number of lanthanide
  core::Size
  get_n_multi_data() const;

	/// @brief Give me the vector of PcsDataCenter
  utility::vector1<PcsDataCenter>&
  get_PCS_data_all() ;

	/// @brief print me
  friend std::ostream &
  operator<<(std::ostream& out, const PcsDataCenterManagerSingleton & me);

};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
