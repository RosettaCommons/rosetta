// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/dfire/DFIRE_Potential.hh
/// @brief  DFIRE Potential class declaration
/// @author James Thompson

#ifndef INCLUDED_core_scoring_methods_dfire_DFIRE_Potential_HH
#define INCLUDED_core_scoring_methods_dfire_DFIRE_Potential_HH

// core

// Unit Headers
#include <core/scoring/methods/dfire/DFIRE_Potential.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <boost/unordered_map.hpp>

// Utility headers

#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace methods {
namespace dfire {

class DFIRE_Potential : public utility::pointer::ReferenceCount {
public:
	DFIRE_Potential();

	virtual ~DFIRE_Potential() ; // auto-removing definition from header{}

	core::Real
	eval_dfire_pair_energy(
		core::conformation::Residue const & rsd1,
		core::conformation::Residue const & rsd2
	) const;

	bool is_loaded() const;

	void
	read_potential(std::string const & fn);

private:
	core::Size res_index ( std::string const & res_name ) const;
	core::Size atom_index( std::string const & atom_name ) const;

	bool potential_is_loaded_;

	boost::unordered_map< std::string, core::Size > atom_res_idx_;
	utility::vector1< utility::vector1< core::Real > > potential_;
};

DFIRE_Potential & get_DFIRE_potential();

} // dfire
} // methods
} // scoring
} // core

#endif
