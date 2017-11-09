// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file LegacyMotifScorer.hh
///
/// @brief
/// @author Tim Jacobs

#ifndef INCLUDED_protocols_legacy_sewing_scoring_LegacyMotifScorer_hh
#define INCLUDED_protocols_legacy_sewing_scoring_LegacyMotifScorer_hh

//Unit headers
#include <protocols/legacy_sewing/scoring/LegacyMotifScorer.fwd.hh>
#include <protocols/legacy_sewing/scoring/LegacyAssemblyScorer.hh>

//Package headers
#include <protocols/legacy_sewing/conformation/Assembly.hh>

//Core headers
#include <core/scoring/motif/motif_hash_stuff.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/types.hh>

//Utility headers
#include <utility/vector1.hh>
#include <numeric/xyzTransform.hh>

namespace protocols {
namespace legacy_sewing  {
namespace scoring {

class LegacyMotifScorer : public LegacyAssemblyScorer {

public:

	///@brief default construct
	LegacyMotifScorer();

	virtual ~LegacyMotifScorer(){}

	virtual
	core::Real
	score(
		AssemblyCOP assembly
	);

	core::Real
	norm_motif_score(
		AssemblyCOP assembly
	);

	core::Real
	full_motif_score(
		AssemblyCOP assembly
	);

	numeric::xyzTransform<core::Real>
	get_stub(
		utility::vector1<SewSegment> const & segments,
		core::Size segment_num,
		core::Size resnum
	) const;

	core::Real
	get_score(
		numeric::xyzTransform<core::Real> stub1,
		char ss1,
		char aa1,
		numeric::xyzTransform<core::Real> stub2,
		char ss2,
		char aa2
	) const;

	void
	dump_motif(
		AssemblyCOP assembly
	) const;

protected:

	core::scoring::motif::MotifHashManager & mman_;
	core::chemical::ResidueTypeSetCOP res_type_set_;

};


} //scoring namespace
} //legacy_sewing namespace
} //protocols namespace

#endif
