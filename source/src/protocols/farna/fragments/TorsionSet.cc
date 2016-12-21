// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/farna/fragments/TorsionSet.hh>
#include <protocols/farna/fragments/FragmentLibrary.hh>
#include <core/chemical/rna/util.hh>
#include <ObjexxFCL/StaticIndexRange.hh>

namespace protocols {
namespace farna {
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

//////////////////////////////////////////////////////////////////////
TorsionSet &
TorsionSet::operator =(
	TorsionSet const & src
) {
	size_ = src.size_;

	for ( Size offset = 0; offset < size_; offset++ ) {
		for ( Size j = 1; j <= core::chemical::rna::NUM_RNA_TORSIONS; j++ ) {
			torsions( j, offset) = src.torsions( j, offset);
		}
		torsion_source_name( offset ) = src.torsion_source_name( offset );
		secstruct( offset ) = src.secstruct( offset );
	}

	non_main_chain_sugar_coords_defined = src.non_main_chain_sugar_coords_defined;

	if ( non_main_chain_sugar_coords_defined ) {
		non_main_chain_sugar_coords.dimension( SRange(0,size_), 3, 3 );
		for ( Size offset = 0; offset < size_; offset++ ) {
			for ( Size j = 1; j <= 3; j++ ) {
				for ( Size k = 1; k <= 3; k++ ) {
					non_main_chain_sugar_coords( offset, j, k ) =
						src.non_main_chain_sugar_coords( offset, j, k );
				}
			}
		}
	}

	index_in_vall_ = src.index_in_vall_;

	return *this;
}

}
}
}
