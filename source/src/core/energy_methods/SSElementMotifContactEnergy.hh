// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/motif/MotifScore.cc
/// @brief  Will's motif score to determine how well packed the protein core is implemented as energy function
/// @author TJ Brunette


#ifndef INCLUDED_core_energy_methods_SSElementMotifContactEnergy_hh
#define INCLUDED_core_energy_methods_SSElementMotifContactEnergy_hh


// Package headers

#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/motif/motif_hash_stuff.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>


// Utility headers
#include <utility/vector1.hh>

#include <set> // AUTO IWYU For multiset

namespace core {
namespace energy_methods {


class SSElementMotifContactEnergy : public core::scoring::methods::WholeStructureEnergy  {
public:
	typedef core::scoring::methods::WholeStructureEnergy parent;

public:
	SSElementMotifContactEnergy();

	core::scoring::methods::EnergyMethodOP
	clone() const override {
		return utility::pointer::make_shared< SSElementMotifContactEnergy >();
	}

	utility::vector1<std::pair<Size,Size> > get_ss_elements(const pose::Pose & pose) const;

	Size which_ssElement(Size res,utility::vector1<std::pair<Size,Size> > ssElements) const;

	Size get_ssElements_in_contact_w_threshold(std::multiset<Size> ssElements_in_contact) const;
	Size get_SSelements_in_contact(Size element,utility::vector1<std::pair<Size,Size> > ssElements, const pose::Pose & pose) const;


	/// @brief Called at the end of the energy evaluation.
	void finalize_total_energy( pose::Pose & pose, core::scoring::ScoreFunction const &, core::scoring::EnergyMap & totals ) const override;


	void indicate_required_context_graphs( utility::vector1< bool > & ) const override {};


	core::Size version() const override;

private:
	core::scoring::motif::MotifHashManager *mman_;
};

} // scoring
} // core

#endif
