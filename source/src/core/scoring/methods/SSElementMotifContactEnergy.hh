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


#ifndef INCLUDED_core_scoring_methods_SSElementMotifContactEnergy_hh
#define INCLUDED_core_scoring_methods_SSElementMotifContactEnergy_hh


// Package headers
#include <basic/datacache/CacheableData.hh>

#include <core/scoring/methods/WholeStructureEnergy.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/motif/motif_hash_stuff.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>


// Utility headers
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace methods {

using utility::vector1;
using core::Size;
using core::Real;
using std::set;

class SSElementMotifContactEnergy : public WholeStructureEnergy  {
public:
    typedef WholeStructureEnergy parent;

public:
    SSElementMotifContactEnergy();

    virtual
    EnergyMethodOP
    clone() const {
        return EnergyMethodOP( new SSElementMotifContactEnergy );
    }

    vector1<std::pair<Size,Size> > get_ss_elements(const pose::Pose & pose) const;

    Size which_ssElement(Size res,vector1<std::pair<Size,Size> > ssElements) const;

    Size get_ssElements_in_contact_w_threshold(std::multiset<Size> ssElements_in_contact) const;
    Size get_SSelements_in_contact(Size element,vector1<std::pair<Size,Size> > ssElements, const pose::Pose & pose) const;


    /// @brief Called at the end of the energy evaluation.
    virtual void finalize_total_energy( pose::Pose & pose, ScoreFunction const &, EnergyMap & totals ) const;


    virtual void indicate_required_context_graphs( utility::vector1< bool > & ) const {};


    virtual core::Size version() const;

private:
    core::scoring::motif::MotifHashManager *mman_;
};

} // motif
} // scoring
} // core

#endif
