// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.


// Rosetta Headers
#include <protocols/rna/denovo/fragments/FragmentLibrary.hh>
#include <protocols/rna/denovo/fragments/TorsionSet.hh>

#include <protocols/rna/denovo/fragments/RNA_Fragments.hh>
#include <protocols/rna/denovo/fragments/FullAtomRNA_Fragments.hh>
#include <core/chemical/rna/util.hh>
#include <core/types.hh>

#include <ObjexxFCL/StaticIndexRange.hh>

namespace protocols {
namespace rna {
namespace denovo {
namespace fragments {

using namespace core;

FragmentLibrary::~FragmentLibrary() = default;

Real FragmentLibrary::get_fragment_torsion( Size const num_torsion,  Size const which_frag, Size const offset ){
	return align_torsions_[ which_frag - 1 ].torsions( num_torsion, offset) ;
}

/////////////////////////////////////////////////////////////////////////////////////////////////
TorsionSet const &
FragmentLibrary::get_fragment_torsion_set( Size const which_frag ) const
{
	return align_torsions_[ which_frag - 1 ];
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void  FragmentLibrary::add_torsion( TorsionSet const & torsion_set ){
	align_torsions_.push_back( torsion_set );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
void  FragmentLibrary::add_torsion(
	FullAtomRNA_Fragments const & vall,
	Size const position,
	Size const size
) {
	TorsionSet torsion_set( size, position );

	for ( Size offset = 0; offset < size; offset++ ) {
		for ( Size j = 1; j <= core::chemical::rna::NUM_RNA_TORSIONS; j++ ) {
			torsion_set.torsions( j, offset) = vall.torsions( j, position+offset);
		}
		torsion_set.torsion_source_name( offset ) = vall.name( position+offset );
		torsion_set.secstruct( offset ) = vall.secstruct( position+offset );

		//Defined non-ideal geometry of sugar ring -- to keep it closed.
		if ( vall.non_main_chain_sugar_coords_defined() ) {
			torsion_set.non_main_chain_sugar_coords_defined = true;
			torsion_set.non_main_chain_sugar_coords.dimension( SRange(0,size), 3, 3 );
			for ( Size j = 1; j <= 3; j++ ) {
				for ( Size k = 1; k <= 3; k++ ) {
					torsion_set.non_main_chain_sugar_coords( offset, j, k ) =
						vall.non_main_chain_sugar_coords( position+offset, j, k );
				}
			}
		} else {
			torsion_set.non_main_chain_sugar_coords_defined = false;
		}

	}

	align_torsions_.push_back( torsion_set );
}

/////////////////////////////////////////////////////////////////////////////////////////////////
Size FragmentLibrary::get_align_depth() const {
	return align_torsions_.size();
}

} //fragments
} //denovo
} //rna
} //protocols
