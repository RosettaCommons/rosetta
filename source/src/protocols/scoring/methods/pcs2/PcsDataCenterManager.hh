// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsDataCenterManager.hh
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenterManager_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsDataCenterManager_hh

// Package headers
#include <protocols/scoring/methods/pcs2/PcsDataCenter.hh>

// Project headers
#include <basic/datacache/CacheableData.hh>
#include <core/types.hh>

#include <utility/vector1.hh>


// Utility headers

// Numeric headers

// ObjexxFCL headers

// c++ headers

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsDataCenterManager : public basic::datacache::CacheableData {
private:

	utility::vector1<PcsDataCenter> PCS_data_all_;

public:

	PcsDataCenterManager(); //construct

	~PcsDataCenterManager(); //destruct

	PcsDataCenterManager(PcsDataCenterManager const &other); //copy

	PcsDataCenterManager & // =
	operator=( PcsDataCenterManager const & src );

	virtual basic::datacache::CacheableDataOP //clone
	clone() const;

	/// @brief Give me the number of lanthanide
	core::Size
	get_n_multi_data() const;

	/// @brief Give me the vector of PcsDataCenter
	utility::vector1<PcsDataCenter>&
	get_PCS_data_all() ;

	/// @brief print me
	friend std::ostream &
	operator<<(std::ostream& out, const PcsDataCenterManager & me);

};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
#endif
