// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/packing/HolesParamsRes.cc
/// @brief  Packing Score Params
/// @author Will Sheffler

#ifndef INCLUDED_core_scoring_packing_HolesParamsRes_hh
#define INCLUDED_core_scoring_packing_HolesParamsRes_hh

#include <core/types.hh>




#include <utility/vector1.hh>
#include <map>

namespace core {
namespace scoring {
namespace packing {


class HolesParamsRes {
public:

	HolesParamsRes( std::string const & fname = std::string("") ) {
		read_data_file( fname );
	}

	void read_data_file( std::string const & fname );

	utility::vector1<core::Real> const & param(std::string res) const {
		return params_.find(res)->second;
	}

	core::Real rho() const { return rho_; }

	core::Real rho( std::string res ) const {
		if ( have_params(res) ) {
			return params_.find(res)->second[params_.find(res)->second.size()];
		} else return 0;
	}

	std::string residues( core::Size i ) const {
		return residues_[i];
	}

	bool have_params( std::string res ) const {
		return params_.find(res) != params_.end();
	}

private:
	utility::vector1<std::string> residues_;
	std::map<std::string,core::Size> res_size_;
	std::map<std::string,utility::vector1<core::Real> > params_;
	core::Real rho_;
};

} // end packing
} // end scoring
} // end core

#endif // INCLUDED_core_scoring_packing_HolesParamsRes_HH
