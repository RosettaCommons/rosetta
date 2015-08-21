// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/vall_lookback/VallLookbackEnergy.hh
/// @brief  lookback to native scores.
/// @author TJ Brunette (brunette@gmail.com)

/// This energy function requires A. the VALL lookback database gets added to the potential
/// B. For efficiency when locations are sampled they need to be marked changed in the VallLookback datacache object

#ifndef INCLUDED_core_scoring_methods_FragVallLookbackEnergy_hh
#define INCLUDED_core_scoring_methods_FragVallLookbackEnergy_hh

#include <basic/datacache/CacheableData.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/methods/vall_lookback/VallLookbackPotential.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

// Utility headers


namespace core {
namespace scoring {
namespace methods {

class VallLookbackEnergy : public WholeStructureEnergy {
public:
	typedef WholeStructureEnergy  parent;

	VallLookbackEnergy();

	/// clone
	virtual EnergyMethodOP clone() const;

	void
	finalize_total_energy(
		pose::Pose & pose,
		ScoreFunction const &,
		EnergyMap & totals
	) const;

	void
	indicate_required_context_graphs(
		utility::vector1< bool > & /*context_graphs_required*/
	) const {}

private:

	VallLookbackPotential const & potential_;
	virtual core::Size version() const;

};


} // methods
} // scoring
} // core

#endif

