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
/// @detailed
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange
/// @author Christopher Miles (cmiles@uw.edu)

// Unit Headers
#include <protocols/jumping/SheetBuilder.hh>

// Package Headers
#include <core/fragment/SecondaryStructure.hh>
#include <protocols/jumping/SameStrand.hh>
#include <core/scoring/dssp/PairingsList.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/FoldTree.hh>

// Utility headers
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/StaticIndexRange.hh>
#include <ObjexxFCL/format.hh>

// numeric headers
#include <numeric/random/random.hh>
#include <numeric/numeric.functions.hh>

//// C++ headers
#include <cstdlib>

#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "protocols.jumping" );

namespace protocols {
namespace jumping {

using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
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
// Here Nj is the number of strand pairs in sheet number j. So the
// number of forced pairings will be (Nj-1). The sheet can
// get other strands during the folding simulation -- this is
// just specifying how many Rosetta should actual build from
// the start using the broken chain stuff.
//
// So the total number of forced pairings will be:
//  N1 + N2 + N3 + ... + Nk
//
// For example, to specify two strand pairings in two different
// sheets, use the args "-sheet1 2 -sheet2 2"
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
SheetBuilder::SheetBuilder( core::fragment::SecondaryStructureOP ss, core::scoring::dssp::PairingsList const& pairings, SheetTopology const& sheet_topol) :
  total_residue_( ss->total_residue() ),
  pairings_( pairings ),
  same_strand_( new SameStrand( ss ) ),
  secondary_structure_( ss ),
	sheet_sizes_( sheet_topol ),
  bForceSingleSheet_( true ) // this should be set by an option or so
{}

//copy c'stor
SheetBuilder::SheetBuilder( SheetBuilder const& other )
	: BaseJumpSetup(other),
		total_residue_( other.total_residue_ ),
		pairings_( other.pairings_ ),
		same_strand_( other.same_strand_ ),
		secondary_structure_( other.secondary_structure_ ),
		sheet_sizes_( other.sheet_sizes_ ),
		bForceSingleSheet_( other.bForceSingleSheet_ ) // this should be set by an option or so
{}

//d'stor
SheetBuilder::~SheetBuilder() {}

///@brief simply random choice of pairing from pool
void
SheetBuilder::choose_next_pairing( FArray3D_int& sheet_pairing, Size pairing, Size sheet ) const {
  int const	p = static_cast< int >( numeric::random::rg().uniform() * pairings_.size() ) + 1;
	//	tr.Trace << "Picked pairing " << p << " out of " << pairings_.size() << std::endl;
	runtime_assert( p>=1 && p<= (int) pairings_.size() );
  // you should replace sheet_pairing with array of Pairings!!!!!
  // then just use operator=
	sheet_pairing( 1, pairing, sheet ) = pairings_[ p ].Pos1();
	sheet_pairing( 2, pairing, sheet ) = pairings_[ p ].Pos2();
	sheet_pairing( 3, pairing, sheet ) = pairings_[ p ].Orientation();
	sheet_pairing( 4, pairing, sheet ) = pairings_[ p ].Pleating();
}

///////////////////////////////////////////////////////////////////////////////
///@detail return true if parings have a strand in common
bool
SheetBuilder::check_pairing_intersect (
  FArray1A_int p1,
  FArray1A_int p2
) const
{
  for ( int i = 1; i <= 2; ++i ) {
    for ( int j = 1; j <= 2; ++j ) {
      if ( same_strand_->eval( p1(i), p2(j) ) ) return true;
    }
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
bool
SheetBuilder::check_sheet_pairings(
	FArray2A_int pairing_list,
	const int last_pairing,
	const bool force_single_sheet
) const
{
	if ( last_pairing < 2 ) return true; // first strand is OK by default
	//	tr.Trace << "check_sheet_pairings ... last_pairing: " << last_pairing << std::endl;
	const int max_sheet_size( 40 ); // can we get rid of this ? XXXX
	pairing_list.dimension( 4, max_sheet_size );

	int total_common_strands( 0 );
	for ( int i = 1; i <= last_pairing-1; ++i ) {
		int common_strands(0);

		//count comman_strands and check if pleating is compatible
		const bool ok( check_two_pairings( pairing_list(1,last_pairing),
																			 pairing_list(1,i), common_strands ) );

		//				tr.Trace << "( " << i << "-"<<last_pairing <<"  common_strands: " << common_strands << ") " << ( ok ? "good" : "bad" ) << " pleating" << std::endl;
		//wrong pleating or two common_strands doesn't help to build beta-sheets
		if ( !ok || common_strands > 1 ) return false;

		//???
		total_common_strands += common_strands;
	}

	return ( ( total_common_strands >= 1 || !force_single_sheet ) &&
					 ( total_common_strands < 2 ) );
}

///////////////////////////////////////////////////////////////////////////////
///@detail count how many strands are in common : 0, 1, or 2
/// if 1 or 2 strands are in common check if pleating of the two pairings is compatible
bool
SheetBuilder::check_two_pairings(
	FArray1A_int pairing1,
	FArray1A_int pairing2,
	int & common_strands
) const
{
	pairing1.dimension(4);
	pairing2.dimension(4);

	common_strands = 0;

	core::scoring::dssp::Pairing p1 ( pairing1 );
	core::scoring::dssp::Pairing p2 ( pairing2 );

	if ( same_strand_->eval( p1.Pos1(), p2.Pos1() ) && same_strand_->eval( p1.Pos2(), p2.Pos2() ) ) {
		common_strands = 2;
	} else if ( same_strand_->eval(p1.Pos1(), p2.Pos2()) && same_strand_->eval(p1.Pos2(), p2.Pos1()) ) {
		common_strands = 2;
		// could reverse either one:
		p2.reverse();
	} else if ( same_strand_->eval(p1.Pos1(), p2.Pos1()) ) {
		common_strands = 1;
	} else if ( same_strand_->eval(p1.Pos1(), p2.Pos2()) ) {
		p2.reverse();
		common_strands = 1;
	} else if ( same_strand_->eval(p1.Pos2(), p2.Pos1()) ) {
		p1.reverse();
		common_strands = 1;
	} else if ( same_strand_->eval(p1.Pos2(), p2.Pos2()) ) {
		p1.reverse();
		p2.reverse();
		common_strands = 1;
	}

	// now we have set things up so that p1.pos1 and p2.pos1
	// are in the same strand.
	//
	// and maybe also p1.pos2 and p2.pos2, if common_strands == 2
	runtime_assert ( common_strands == 0 || same_strand_->eval(p1.Pos1(), p2.Pos1()) );

	const int seqsep_mod2( numeric::mod( std::abs( (int) p1.Pos1() - (int) p2.Pos1() ), 2) );

	return ( ( common_strands == 0 ) ||
					 ( common_strands == 1 && seqsep_mod2 == 0 && p1.Pleating() != p2.Pleating() ) ||
					 ( common_strands == 1 && seqsep_mod2 == 1 && p1.Pleating() == p2.Pleating() ) ||
					 ( common_strands == 2 && seqsep_mod2 == 0 && p1.Pleating() == p2.Pleating() ) ||
					 ( common_strands == 2 && seqsep_mod2 == 1 && p1.Pleating() != p2.Pleating() ) );
}

bool
SheetBuilder::check_next_pairing( FArray3D_int& sheet_pairings, Size pairing, Size sheet ) const {
	for ( Size ii = 1; ii<=2; ii++ ) { //check if both residues of the pairing make sense at all
		Size pos = sheet_pairings( ii, pairing, sheet );
		if ( pos == 1 || pos >= total_residue_ ) return false;
	}

  // full compatibility chec
	//check1:  dont want to intersect previous sheets
	for ( int prev_sheet = 1; prev_sheet <= (int) sheet - 1; ++prev_sheet ) {
		int const num_pairings = sheet_sizes_[prev_sheet];
    for ( int prev_pairing = 1;  prev_pairing <= num_pairings; ++prev_pairing ) {
      if ( check_pairing_intersect( sheet_pairings( 1, prev_pairing, prev_sheet ),
					sheet_pairings( 1, pairing, sheet ) ) ) {
				//				tr.Debug << "paring intersect : check failed for pairing" << pairing << " sheet: " << sheet << std::endl;
				return false;
      }
    }
  }
	//check2:
	return check_sheet_pairings(
			sheet_pairings( 1, 1, sheet ),
			pairing, bForceSingleSheet_ );
}

JumpSample SheetBuilder::create_jump_sample() const{
	sheet_sizes_ = create_new_random_topol();
	FArray1A_int cuts;
	for ( int trial = 1; trial < 30; trial ++ ) {
		core::scoring::dssp::PairingsList jump_pairings;
		bool success = builder_loop( jump_pairings );
		if ( !success ) {
			tr.Warning << "redo builder_loop failed, will try " << 10-trial << " more times" << std::endl;
			sheet_sizes_ = create_new_random_topol();
			continue;
		}
		JumpSample jumps( total_residue_, jump_pairings, *secondary_structure_ );
		//		tr.Debug << " created jump sample " << jumps << std::endl;
		/*	if ( tr.Debug.visible() && jumps.is_valid() ) { //debugging output
			//			tr.Debug << "same_strand: ";
			typedef utility::vector1< int > Int_List;
			Int_List jres;
			for ( core::scoring::dssp::PairingsList::const_iterator it = jump_pairings.begin(),
							eit = jump_pairings.end(); it != eit; ++it ) {
				jres.push_back( it->Pos1() );
				jres.push_back( it->Pos2() );
			}
			for ( Int_List::iterator it=jres.begin(), eit=jres.end(); it!=eit; ++it ) {
				for ( Int_List::iterator sit=jres.begin(), seit=jres.end(); sit!=seit; ++sit ) {
					if ( *it > *sit ) {
						//						tr.Debug << *it << ":"<< *sit << " " << (same_strand_->eval(*it,*sit ) ? "same" : "diff") << " ";
					}
				}
			}
			//			tr.Debug << std::endl;
		} // is_valid()
		*/
		if ( jumps.is_valid() ) return jumps;
		//		tr.Debug << " ...which had corrupted fold-tree. Try again! " << std::endl;
	}
	//utility_exit_with_message( "impossible to find a valid fold-tree with given sheet-topology and pairings" );
	return JumpSample(); //to make compiler happy
}

bool
SheetBuilder::builder_loop( core::scoring::dssp::PairingsList& jump_pairings ) const {
  Size const max_sheet_size( 40 );
  Size const max_sheets( 40 ); // can we determine these from something ?
	int num_sheets( sheet_sizes_.size() );
  //// fill the sheet_pairing array:
  FArray3D_int sheet_pairing( 4, max_sheet_size, max_sheets );
	tr.Info << "Start Sheet Building for " << num_sheets << " sheets of size ";
	for ( int ii = 1; ii <= num_sheets; ii++ ) tr.Info << sheet_sizes_[ ii ] << " ";
	tr.Info << std::endl;
  bool success = false;
	int failed_once = 0;
  for ( int tries1 = 0; tries1 < 30 && !success; ++tries1 ) { // redo_same_strand: try different sheet boundaries
		same_strand_->redo();
    for ( int tries2 = 0; tries2 < 30 && !success; ++tries2 ) { // sheet_fail
			tr.Debug << "SheetBuilder-loop round: " << tries1 << "/" << tries2 << std::endl;
      // here we fill the array "sheet_pairing" that has all the pairings
      // organized by sheets

      int tries3 ( 0 );
      int const max_tries3 ( 1000 );

      // for each sheet
      for ( int sheet = 1; sheet <= num_sheets; ++sheet ) {
				// select sheet_size-1 pairings
				for ( int pairing = 1; pairing <= (int) sheet_sizes_[ sheet ];	++pairing ) {
					choose_next_pairing( sheet_pairing, pairing, sheet ); //draw randomly from pool
									for ( int s = 1; s<=sheet; s++ ) {
											for ( int p =1; p<=pairing; p ++ ) {
												//tr.Debug << s << ":" << p << " pairs: " << sheet_pairing(1, p, s) << "-" <<sheet_pairing(2, p, s) << std::endl;
											}
									}


					while ( !( success=check_next_pairing( sheet_pairing, pairing, sheet ) )  && tries3 < max_tries3 ) {
						choose_next_pairing( sheet_pairing, pairing, sheet );
						for ( int s = 1; s<=sheet; s++ ) {
							for ( int p =1; p<=pairing; p ++ ) {
								//							tr.Debug << s << ":" << p << " pairs: " << sheet_pairing(1, p, s) << "-" <<sheet_pairing(2, p, s) << std::endl;
							}
							}

						++tries3;
					} // loop to choose pairings
					//					tr.Trace << "Pairing " << pairing << " chosen for sheet " << sheet << std::endl;
				} // for pairings per sheet
      } // for sheets
    } // for tries2 ( sheet_fail )
    // if not successful try different sheet boundaries
    if ( !success ) {
			failed_once=tries1;
			tr.Info << "redo same_strand: too many tries with this one!" <<
      				format::SS( tries1 ) << std::endl;
      // makes stochastic decisions about strand boundaries:
      // some decisions may not be compatible with the desired pairings
      // and our logic for choosing sheets
		}
  } // for tries1 ( redo_same_strand)

  if ( !success ) {
		return false;
		//  utility_exit_with_message( "problem in SheetBuilder::builder_loop(): tries1 > 100." );
  }
	if ( failed_once ) tr.Warning << "figured out a valid strand after " << failed_once << " outer loop iterations " << std::endl;
  // now load these pairings into the jump_pairings list:
  jump_pairings.clear();
  for ( int sheet = 1; sheet <= num_sheets; ++sheet ) {
    for ( int pairing = 1, pe = sheet_sizes_[sheet]; pairing <= pe; ++pairing ) {

			core::scoring::dssp::Pairing p( sheet_pairing( 1,pairing,sheet) );
			jump_pairings.push_back( p );

			//      tr.Debug << "sheet_pairing:" << SS(sheet) << SS(pairing) << ' ' <<
			//				p.Pos1() << ' ' << p.Pos2() << ' ' <<
			//				p.Orientation() << ' ' << p.Pleating() << std::endl;
    }
  }
	return true;
}

} //protocols
} //jumping
