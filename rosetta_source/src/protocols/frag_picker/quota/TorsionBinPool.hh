// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/quota/TorsionBinPool.hh
/// @brief a quota pool based on torsion bin  prediction
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#ifndef INCLUDED_protocols_frag_picker_quota_TorsionBinPool_hh
#define INCLUDED_protocols_frag_picker_quota_TorsionBinPool_hh

#include <protocols/frag_picker/QuotaPool.hh>

// package headers
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.fwd.hh>

#include <core/fragment/SecondaryStructure.hh>

// utility headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

namespace protocols {
namespace frag_picker {
namespace quota {

/// @brief represents a single pool used by quota selector
class TorsionBinPool : public QuotaPool {
public:
	/// @brief Creates a pool of a given size and name
	/// @param size - how many fragments may fit into this pool
	/// @param name - name assigned to this pool. This in general may be any string that
	///	later allows one control pool's behavior from a flag file
	TorsionBinPool(Size size,std::string pool_name,
		core::fragment::SecondaryStructureOP prediction) :
		QuotaPool(size,pool_name) {
		prediction_ = prediction;
	}

	virtual ~TorsionBinPool();

	void show_availability(std::ostream &) const;

	/// @brief try to insert a given fragment candidate into this pool
	bool try_fragment(ScoredCandidate & candidate);

	void restart(Size,Size);

private:
        core::fragment::SecondaryStructureOP prediction_;
	Size nh_, ne_, nl_;
	Size max_h_,max_l_,max_e_;
	Size frag_size_;
        inline Size round(Real x) { return Size(x > 0.0 ? x + 0.5 : x - 0.5); }
};

} // quota
} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_quota_TorsionBinPool_HH */
