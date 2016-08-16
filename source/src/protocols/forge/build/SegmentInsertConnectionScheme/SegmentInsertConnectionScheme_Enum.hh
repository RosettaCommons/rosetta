// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_forge_build_SegmentInsert_SegmentInsertConnectionScheme_SegmentInsertConnectionScheme_Enum_hh
#define INCLUDED_protocols_forge_build_SegmentInsert_SegmentInsertConnectionScheme_SegmentInsertConnectionScheme_Enum_hh

namespace protocols {
namespace forge {
namespace build {
// description namespace to hold enum
namespace SegmentInsertConnectionScheme {
/// @brief connect insertion on its N-side, C-side, or decide randomly
///  between the two
enum Enum {
	RANDOM_SIDE,
	N,
	C
};

}
}
}
}

#endif
