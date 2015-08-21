// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief secondary structure will hold statistics about secondary structure predictions
/// sources can be from
///      - fragments
///      - psipred files ? other stuff
///
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange


// Unit Headers
#include <protocols/jumping/SheetBuilder.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>

// Package Headers
//#include <protocols/jumping/SecondaryStructure.hh>
//#include <protocols/jumping/SameStrand.hh>
//#include <core/scoring/dssp/PairingsList.hh>

// Project Headers
#include <core/types.hh>

//#include <core/kinematics/FoldTree.hh>


// Utility headers
#include <basic/Tracer.hh>
//#include <utility/io/izstream.hh>
#include <utility/io/util.hh>

// ObjexxFCL Headers
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray1A.hh>
//#include <ObjexxFCL/FArray2A.hh>
//#include <ObjexxFCL/FArray3D.hh>

//#include <ObjexxFCL/StaticIndexRange.hh>

//#include <ObjexxFCL/format.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
//#include <numeric/numeric.functions.hh>

//// C++ headers
#include <cstdlib>

#include <core/fragment/SecondaryStructure.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//#include <string>
//#include <vector>


static thread_local basic::Tracer tr( "protocols.jumping" );

namespace protocols {
namespace jumping {

using namespace ObjexxFCL;
using namespace core;


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// NOTES ( old comment copied from jumping_pairings )
//
////////////////////////
// PAIRINGS FILE FORMAT:
//
// one line per pairing, 4 numbers per line, white-space delimited
//
// format: "pos1 pos2 orientation pleating"
//
// pos1 and pos2 are the sequence numbers of the beta-paired positions
//
// orientation is the beta-strand orientation:
// "1" for antiparallel, "2" for parallel
// or just "A" or "P"
//
// pleating determines the pleat of the beta-carbons
// "1" if the N and O of pos1 are pointed away from pos2
// "2" if the N and O of pos1 are pointed toward pos2
//
// eg, in the antiparallel case, a pleat of 2 would mean that
// there are two backbone-backbone hydrogen bonds between pos1 and pos2
// In the parallel case a pleat of 2 means that pos1 is Hbonding with
// pos2-1 and pos2+1
//
// if you check out rosetta_benchmarks, in the directory 1d3z/
// is a pairing file "pairings.dat" with two pairings in it.
// These are native pairings, ie they match the pdb structure
// 1d3z.pdb in the same directory. In Fortran you had to tell
// Rosetta how many pairings to expect; that's why the first line
// has a "2" on it. This is no longer necessary.

//////////////////////
// COMMAND LINE SYNTAX
//
// for ab initio folding with pairings: the usual "xx 1xyz _ -silent"
// plus:
//
// -pairing_file <pairings-file>
//
// with no other arguments and it will try to build decoys with
// *all* the pairings in the file. You can also specify what
// kind of sheet topology Rosetta should try to construct:
//
// -sheet1 <N1> -sheet2 <N2> ... -sheetk <Nk>
//
// Here Nj is the number of strands in sheet number j. So the
// number of forced pairings will be (Nj-1). The sheet can
// get other strands during the folding simulation -- this is
// just specifying how many Rosetta should actual build from
// the start using the broken chain stuff.
//
// So the total number of forced pairings will be:
//  N1-1 + N2-1 + N3-1 + ... + Nk-1
//
// For example, to specify two strand pairings in two different
// sheets, use the args "-sheet1 2 -sheet2 2"
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////


RandomSheetBuilder::RandomSheetBuilder( core::fragment::SecondaryStructureOP ss, core::scoring::dssp::PairingsList const& pairings, SheetTopology const& sheet_topol) :
	SheetBuilder( ss, pairings, sheet_topol ),
	input_sheet_sizes_( sheet_topol )
{
	//...
}


//default do nothing always use input_sheet_sizes_ as sheet_sizes_.
SheetBuilder::SheetTopology RandomSheetBuilder::create_new_random_topol() const
{
	Size num_sheets = std::max( 1, static_cast< int >( numeric::random::rg().uniform() * (input_sheet_sizes_.size() + 1) ) );
	tr.Debug << "random choice: num_sheets: " << num_sheets << std::endl;

	// generate random sequence from 1 .. N
	utility::vector1< Size > strand_ids;
	for ( Size i = 1; i <= input_sheet_sizes_.size(); i ++ ) {
		strand_ids.push_back( i );
	}

	numeric::random::random_permutation( strand_ids, numeric::random::rg() );
	numeric::random::random_permutation( strand_ids, numeric::random::rg() ); //want it really random

	tr.Debug << "strand_ids.size(): "
		<< strand_ids.size() << " ";
	utility::io::write_vector( tr.Debug, strand_ids );
	tr.Debug << std::endl;

	SheetTopology new_sheet_sizes;
	Size trials = 20;
	while ( new_sheet_sizes.size() < num_sheets && trials-- > 0 ) {
		for ( Size i = 1; i<=strand_ids.size() && new_sheet_sizes.size()<num_sheets; i++ ) {
			int nr( static_cast< int >( numeric::random::rg().uniform() * input_sheet_sizes_[ strand_ids[ i ] ] ) + 1 );
			if ( nr > 0 ) {
				new_sheet_sizes.push_back( nr );
			}
		}
	}

	tr.Debug << "chosen sheet parameters: ";
	utility::io::write_vector( tr.Debug, new_sheet_sizes );
	tr.Debug << std::endl;

	return new_sheet_sizes;
}

} //protocols
} //jumping
