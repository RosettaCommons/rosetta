// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/rna/denovo/fragments/TorsionSet.hh>
#include <protocols/rna/denovo/fragments/FragmentLibrary.hh>
#include <core/chemical/rna/util.hh>
#include <ObjexxFCL/StaticIndexRange.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace fragments {

using namespace core;

TorsionSet::TorsionSet( Size const size, Size const index_in_vall ){
	torsions.dimension( core::chemical::rna::NUM_RNA_TORSIONS, SRange(0, size) );
	torsion_source_name.dimension( SRange(0, size), std::string( 4, ' ' )  );
	secstruct.dimension( SRange(0, size), 'L' );
	non_main_chain_sugar_coords_defined = false;
	size_ = size;
	index_in_vall_ = index_in_vall;
}

} //fragments
} //denovo
} //rna
} //protocols
