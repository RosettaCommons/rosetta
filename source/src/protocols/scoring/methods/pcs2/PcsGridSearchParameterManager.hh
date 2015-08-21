// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

//////////////////////////////////////////////
/// @file protocols/scoring/methods/pcs2/PcsGridSearchParameterManager.hh
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////

#ifndef INCLUDED_protocols_scoring_methods_pcs2_PcsGridSearchParameterManager_hh
#define INCLUDED_protocols_scoring_methods_pcs2_PcsGridSearchParameterManager_hh

// Unit Headers

// Package Headers
#include <protocols/scoring/methods/pcs2/PcsGridSearchParameter.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

class PcsGridSearchParameterManager : public utility::SingletonBase< PcsGridSearchParameterManager >
{
public:
	friend class utility::SingletonBase< PcsGridSearchParameterManager >;

private:

	utility::vector1<PcsGridSearchParameter> grid_s_p_all_;

	PcsGridSearchParameterManager();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static PcsGridSearchParameterManager * create_singleton_instance();

public:

	/// @brief give me the number of paramagnetic center
	core::Size
	get_n_multi_data() const;

	/// @brief Incremente the number of paramagnetic center
	void
	incremente_n_multi_data();

	/// @ Give me the PcsGridSearchParamater of the center i_multi_data
	PcsGridSearchParameter &
	get_grid_search_parameters(core::Size i_multi_data);

	/// @ Re init the singleton to default value
	void
	re_init();

};

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols

#endif
