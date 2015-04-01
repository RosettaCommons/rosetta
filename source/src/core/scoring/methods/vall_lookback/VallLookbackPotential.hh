// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/VallLookbackPotenial.fwd 
/// @brief  threshold based distance metric that measures from any fragment to the structure in the vall
/// @author TJ Brunette (tjbrunette@gmail.com)

#ifndef INCLUDED_core_scoring_methods_VallLookbackPotential_hh
#define INCLUDED_core_scoring_methods_VallLookbackPotential_hh

#include <core/scoring/methods/vall_lookback/VallLookbackPotential.fwd.hh>

#include <core/pose/Pose.fwd.hh>
#include <core/indexed_structure_store/FragmentStore.hh>

#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <map>

using core::Size;
using core::Real;

namespace core {
namespace scoring {
namespace methods {
    class VallLookbackPotential: public utility::pointer::ReferenceCount {
public:
    VallLookbackPotential();
    Real lookback(pose::Pose & pose, Size resid) const; //note only reports rmsd if above threshold. 
    Real lookback(const pose::Pose & pose, Size resid) const;  //const by only using the datacache object if it exists.
    Real lookbackRegion(pose::Pose & pose, Size startRes, Size endRes) const;
    Real lookbackRegion(const pose::Pose & pose, Size startRes, Size endRes) const;
    Real lookback(pose::Pose & pose) const; //note only reports rmsd if above threshold. 
    Real lookback(const pose::Pose & pose) const; //note only reports rmsd if above threshold. 
    Size fragLengthInDB() const;
private:
    Real lookback_db(pose::Pose & pose, Size resid) const;
    Real lookback_db(const pose::Pose & pose, Size resid) const;
    std::map<Size, core::indexed_structure_store::FragmentStoreOP> abegoHashedFragmentStore_;
};
}// end methods
}// end scoring
} //end core
#endif

