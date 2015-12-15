// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/PocketExemplarMultifunc.hh
/// @brief  Pocket comparison multifunction class
/// @author David Johnson


#ifndef INCLUDED_protocols_pockets_PocketExemplarMultifunc_hh
#define INCLUDED_protocols_pockets_PocketExemplarMultifunc_hh

// Package headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Multifunc.hh>
#include <protocols/pockets/PocketGrid.hh>

#include <protocols/pockets/PocketExemplarMultifunc.fwd.hh>
#include <utility/vector1.hh>


using namespace core;
using namespace std;


namespace protocols {
namespace pockets {

/// @brief Pocket multifunction class, does objective function of optimization
class PocketExemplarMultifunc : public core::optimization::Multifunc {

public: // Constructor/Destructor

	// Constructor
	PocketExemplarMultifunc(std::string const input_pdb_name, std::string const resid, core::Real const c_rad, core::Real const rep_weight, utility::vector1<core::Real>& p_min, utility::vector1<core::Real>& p_max);

	/// @brief Destructor
	inline
	virtual
	~PocketExemplarMultifunc()
	{}

public: // Methods

	// objective func
	virtual
	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	// dfunc
	virtual
	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const;

	//using core::optimization::Multifunc::dump;

	/// @brief Error state reached; dump out current pdb.
	virtual
	void
	dump( core::optimization::Multivec const & vars, core::optimization::Multivec const &) const;

private:

	core::pose::Pose input_pose;
	std::vector< conformation::ResidueCOP > residues;
	protocols::pockets::PocketGrid pg;
	core::Real cRad;
	core::Real repW;
	core::Real optimal;
	core::Real vdWpen;

}; // PocketExemplarMultifunc


} // namespace Pockets
} // namespace protocols


#endif // INCLUDED_protocols_pockets_PocketExemplarMultifunc_HH
