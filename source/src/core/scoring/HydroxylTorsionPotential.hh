// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/HydroxylTorsionPotential.hh
/// @brief  Term for proton_chi on Ser/Thr/Tyr residues
/// @author Hahnbeom Park (hahnbeom@gmail.com)

#ifndef INCLUDED_core_scoring_HydroxylTorsionPotential_hh
#define INCLUDED_core_scoring_HydroxylTorsionPotential_hh

// Unit Headers
#include <core/scoring/HydroxylTorsionPotential.fwd.hh>

// Package headers
#include <core/chemical/ResidueType.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/scoring/DerivVectorPair.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

struct TorsionParams
{
	utility::vector1< Size > atm;
	Real n, k, delta;
};

class HydroxylTorsionPotential : public utility::pointer::ReferenceCount
{
public:
	typedef boost::unordered_multimap< std::string, TorsionParams >::const_iterator tors_iterator;

public:
	HydroxylTorsionPotential();
	~HydroxylTorsionPotential() override = default;

	Real
	eval_residue_energy(
		conformation::Residue const & rsd
	) const;

	void
	eval_residue_derivative(
		conformation::Residue const & rsd,
		utility::vector1< DerivVectorPair > & atom_derivs
	) const;


private:
	std::string get_restag( core::chemical::ResidueType const & restype ) const;

	void read_database( std::string filename );
	boost::unordered_multimap< std::string, TorsionParams > torsion_params_;

};

} // namespace core
} // namespace scoring

#endif // INCLUDED_core_scoring_HydroxylTorsionPotential_hh
