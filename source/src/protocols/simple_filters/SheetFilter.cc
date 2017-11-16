// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file SheetFilter.cc
/// @brief runs reject or accept filters on pose
/// @details
///   Contains currently: SheetFilter
///
///
/// @author Robert Vernon

// Unit Headers
#include <protocols/simple_filters/SheetFilter.hh>

// Package Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/xyzVector.hh>

//// C++ headers
#include <cstdlib>
#include <string>

#include <utility/vector1.hh>


static basic::Tracer tr( "protocols.simple_filters.SheetFilter" );

using core::Real;
using namespace core;
using namespace basic;
using namespace ObjexxFCL;

namespace protocols {
namespace simple_filters {


bool SheetFilter::apply( core::pose::Pose const & pose ) const {
	if ( !AbinitioBaseFilter::apply( pose ) ) return false;
	if ( sstype_ == "fail" ) return true;
	//car sheet_filter disable
	if ( beta_ < 20 || beta_ratio_ <= 0.2 ) return true;

	//  now the ss-type is read in apply():
	// and the disable decision is made here
	// std::string sstype = get_protein_sstype( pose );
	//  if ( disable_sheet_filter_ ) {
	//   return true;
	//  }

	//SHEET FILTER
	float const distance_cutoff = { 6.5 };
	//std::string protein_sstype;

	//bool sheet_filter (true);

	float min_dotprod;
	float max_distance;
	//float handedness_score_;
	int nstrands,nsheets;
	int local_pairs,nonlocal_pairs;

	//protein_sstype = get_protein_sstype();

	FArray1D_int isN;
	FArray1D_int isCA;
	FArray1D_int isC;
	FArray1D_int res_num;

	int const max_pos(5);
	int const max_res(pose.size());

	FArray2D_float position( 3, max_pos*max_res );

	for ( int i = 1; i <= max_res; i += 1 ) {
		int const itemp = max_pos*(i-1)+1;
		pose.residue(i).atom("N").xyz().x();
		position(1,itemp)   = pose.residue(i).xyz("N").x();
		position(2,itemp)   = pose.residue(i).xyz("N").y();
		position(3,itemp)   = pose.residue(i).xyz("N").z();
		position(1,itemp+1) = pose.residue(i).xyz("CA").x();
		position(2,itemp+1) = pose.residue(i).xyz("CA").y();
		position(3,itemp+1) = pose.residue(i).xyz("CA").z();
		//position(1,itemp+2) = pose.residue(i).xyz("CB").x();
		//position(2,itemp+2) = pose.residue(i).xyz("CB").y();
		//position(3,itemp+2) = pose.residue(i).xyz("CB").z();
		position(1,itemp+3) = pose.residue(i).xyz("C").x();
		position(2,itemp+3) = pose.residue(i).xyz("C").y();
		position(3,itemp+3) = pose.residue(i).xyz("C").z();
		position(1,itemp+4) = pose.residue(i).xyz("O").x();
		position(2,itemp+4) = pose.residue(i).xyz("O").y();
		position(3,itemp+4) = pose.residue(i).xyz("O").z();
	}

	res_num = FArray1D_int( max_pos * max_res );
	isN = FArray1D_int( max_pos * max_res );
	isC = FArray1D_int( max_pos * max_res );
	isCA = FArray1D_int( max_pos * max_res );

	for ( int i = 1, e = max_pos*max_res; i<= e; i+= max_pos ) {
		isN(i)   = 1;
		isN(i+1) = 0;
		isN(i+2) = 0;
		isN(i+3) = 0;
		isN(i+4) = 0;

		isCA(i)   = 0;
		isCA(i+1) = 1;
		isCA(i+2) = 0;
		isCA(i+3) = 0;
		isCA(i+4) = 0;

		isC(i)   = 0;
		isC(i+1) = 0;
		isC(i+2) = 0;
		isC(i+3) = 1;
		isC(i+4) = 0;

		for ( int j = 0; j < max_pos; ++j ) {
			res_num(i+j) = (i/max_pos) + 1;
		}
	}

	FArray1D_int ss( max_res );

	for ( int i = 1; i <= max_res; ++i ) {

		if ( pose.secstruct(i) == 'H' ) {
			ss(i) = 1;
		} else if ( pose.secstruct(i) == 'E' ) {
			ss(i) = 2;
		} else {
			ss(i) = 3;
		}
	}

	handedness_score_ = 0.0;

	int result = 0;

	ingo_sheet_stuff(max_res, ss, max_res*max_pos, position, isN, isCA, isC,
		res_num, distance_cutoff, nstrands, nsheets, max_distance,
		min_dotprod, local_pairs, nonlocal_pairs, result);

	//$$$ if ( result > 0 ) sheet_filter = false;
	//$$$ if ( result == 2 ) sheet_filter = true;  // pass toblerone (not possible)
	//$$$ if ( result == 5 ) sheet_filter = true;  // pass dotprod (not possible)
	//$$$ if ( result == 8 ) sheet_filter = true;  // pass large barrels

	if ( result == 3 ) return false;
	if ( result == 4 ) return false;
	if ( result == 6 ) return false;
	if ( result == 7 ) return false;
	if ( result == 9 ) return false;

	return true;
}

// std::string SheetFilter::get_protein_sstype() {
//  return "misc::sstype::sstype??";
// }


////////////////////////
/// Private: Methods ///
////////////////////////

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM IDENTIFIES THE STRAND DIAMERS
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  strnm - [in/out]? -
/// @param  natm - [in/out]? -
/// @param  atmps - [in/out]? -
/// @param  rnm - [in/out]? -
/// @param  indC - [in/out]? -
/// @param  indN - [in/out]? -
/// @param  strdm - [in/out]? -
/// @param  inddm - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_diamers(
	int const & nres,
	FArray1A_int strnm,
	int natm,
	FArray2D_float const & atmps,
	FArray1D_int & rnm,
	FArray1D_int & indC,
	FArray1D_int & indN,
	FArray2A_float strdm,
	FArray1A_int inddm
) const
{
	strnm.dimension( nres );
	// atmps.dimension( 3, natm );
	rnm.dimension( natm );
	indC.dimension( natm );
	indN.dimension( natm );
	strdm.dimension( nres, 3 );
	inddm.dimension( nres );


	for ( int k = 1; k <= nres; ++k ) {
		strdm(k,1) = 0.0;
		strdm(k,2) = 0.0;
		strdm(k,3) = 0.0;
		inddm(k) = 0;
	}

	for ( int k = 1; k <= natm; ++k ) {
		if ( rnm(k) < nres ) {
			if ( ( strnm(rnm(k)) > 0 ) && ( strnm(rnm(k)+1) > 0 ) ) {
				if ( indN(k) == 1 ) {
					float x1 = atmps(1,k);
					float y1 = atmps(2,k);
					float z1 = atmps(3,k);
					int next_carbon = 0;
					int i = 3;
					while ( next_carbon == 0 ) {
						if ( indC(k+i) == 1 ) {
							float x2 = atmps(1,k+i);
							float y2 = atmps(2,k+i);
							float z2 = atmps(3,k+i);
							strdm(rnm(k),1) = 0.5*(x1+x2);
							strdm(rnm(k),2) = 0.5*(y1+y2);
							strdm(rnm(k),3) = 0.5*(z1+z2);
							inddm(rnm(k)) = 1;
							next_carbon = 1;
						} else {
							++i;
						}
					}
				}
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM FINDS THE DIRECTIONS OF THE STRANDS IN THE SHEETS
///
/// @details
///
/// @param  shnm - [in/out]? -
/// @param  hm - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  slct - [in/out]? -
/// @param  order - [in/out]? -
/// @param  strlbl - [in/out]? -
/// @param  strdr - [in/out]? -
/// @param  proper - [in/out]? - Not used
/// @param  directions - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_find_dir(
	int shnm,
	int hm,
	int & nstr,
	FArray1A_int slct,
	FArray1A_int order,
	FArray1A_int strlbl,
	FArray2A_float strdr,
	// int & proper, // Not used
	FArray1A_int directions
) const
{

	slct.dimension( hm );
	order.dimension( hm );
	strlbl.dimension( nstr );
	strdr.dimension( nstr, 3 );
	directions.dimension( hm );

	FArray2D_float centers( max_nstr, 3 );
	int cnt = 0;

	//js finding strands that are in the sheet

	for ( int j = 1; j <= nstr; ++j ) {
		if ( strlbl(j) == shnm ) {
			//js count them
			++cnt;
			//js put them in an array
			slct(cnt) = j;
		}
	}
	//js centers holds the vectors from the start to stop
	//js of each strand.
	for ( int i = 1; i <= hm; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			centers(i,j) = strdr(slct(i),j);
		}
	}

	//js pretty sure this order thing is screwing things up in ingo_hand.
	//js The directions array is later accessed according the strand
	//js number, but it is being filled according to the order of the strands
	//js in the sheet.
	//js Easiest fix seems to me to be to put the order(i) in reference to
	//js directions here.  Then later calls to directions should do the
	//js right thing?
	//rhiju Double checked -- jack's fix look ok!
	directions(order(1)) = 0;

	for ( int i = 2; i <= hm; ++i ) {
		float dot_prod =
			centers(order(i),1) * centers(order(i-1),1) +
			centers(order(i),2) * centers(order(i-1),2) +
			centers(order(i),3) * centers(order(i-1),3);
		if ( dot_prod < 0.0 ) {
			int const dir1 = directions(order(i-1)) - 1;
			directions(order(i)) = dir1 * dir1;
		} else {
			directions(order(i)) = directions(order(i-1));
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM FINDS THE ORDER OF THE STRANDS IN THE SHEETS
///
/// @details
///
/// @param  shnm - [in/out]? -
/// @param  hm - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  strlbl - [in/out]? -
/// @param  dstrmin - [in/out]? -
/// @param  ngbhct - [in/out]? -
/// @param  order - [in/out]? -
/// @param  sequence - [in/out]? -
/// @param  slct - [in/out]? -
/// @param  rubbish - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_find_ord(
	int shnm,
	int hm,
	int & nstr,
	FArray1A_int strlbl,
	FArray2A_float dstrmin,
	float ngbhct,
	FArray1A_int order,
	FArray1A_int sequence,
	FArray1A_int slct,
	int & rubbish
) const
{

	strlbl.dimension( nstr );
	dstrmin.dimension( nstr, nstr );
	order.dimension( hm );
	sequence.dimension( hm );
	slct.dimension( hm );

	//---- local:
	int cnt,cnt1,cnt2,i,j,proper;//,token;
	FArray1D_int nghbrs_cnt( max_nstr );
	FArray1D_int taken( max_nstr );
	FArray2D_float nghbrs( max_nstr, max_nstr );
	FArray2D_float nghbrs_d( max_nstr, max_nstr );
	//------------------------------------------------------------------------------

	cnt = 0;
	for ( j = 1; j <= nstr; ++j ) {
		if ( strlbl(j) == shnm ) {
			++cnt;
			slct(cnt) = j;
		}
	}

	//     std::cout << "slc " << slct << std::endl;

	for ( i = 1; i <= hm; ++i ) {
		nghbrs_cnt(i) = 0;
		taken(i) = 0;
		for ( j = 1; j <= hm; ++j ) {
			nghbrs_d(i,j) = dstrmin(slct(i),slct(j));
			nghbrs(i,j) = 0;
			if ( nghbrs_d(i,j) <= ngbhct ) {
				nghbrs(i,j) = 1;
				++nghbrs_cnt(i);
			}
		}
	}

	cnt1 = 0;
	cnt2 = 0;
	proper = 0;
	for ( j = 1; j <= hm; ++j ) {
		if ( nghbrs_cnt(j) == 1 ) {
			++cnt1;

		} else if ( nghbrs_cnt(j) == 2 ) {
			++cnt2;
		}
	}

	if ( cnt1 == 2 && cnt2 == hm-2 ) {
		cnt = 0;
		int token = 0;
		proper = 1;
		while ( token == 0 ) {
			++cnt;
			if ( nghbrs_cnt(cnt) == 1 ) token = 1;
		}
		order(1) = cnt;
		taken(cnt) = 1;
		while ( token < hm ) {
			for ( j = 1; j <= hm; ++j ) {
				if ( nghbrs(order(token),j) == 1 && taken(j) == 0 ) {
					taken(j) = 1;
					++token;
					order(token) = j;
				}
			}
		}
		token = 0;
		while ( token < hm ) {
			++token;
			for ( j = 1; j <= hm; ++j ) {
				if ( order(j) == token ) sequence(token) = j;
			}
		}
	}

	//     std::cout << "ord" << std::endl;
	//     std::cout << order << std::endl;
	//     std::cout << sequence << std::endl;

	if ( proper != 1 ) rubbish = 7;
	//  addition to distinguish between barrel-types
	//car a quick fix per ingo:
	//car assigns the value 8 to the rubbish variable for all non-open sheets
	//car (ie they don't have 2 strands with only one neighbour, and the other
	//car strands all have 2).
	if ( ( proper != 1 ) && ( hm > 4 ) ) rubbish = 8;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM CHECKS THE HANDEDNESS OF A SHEET
///
/// @details
///
/// @param  k - [in/out]? -
/// @param  shnm - [in/out]? -
/// @param  hm - [in/out]? -
/// @param  stpppt - [in/out]? -
/// @param  strtpt - [in/out]? -
/// @param  str1 - [in/out]? -
/// @param  str2 - [in/out]? -
/// @param  order - [in/out]? -
/// @param  strdr - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  nres - [in/out]? -
/// @param  sequence - [in/out]? -
/// @param  directions - [in/out]? -
/// @param  scstr - [in/out]? -
/// @param  lctn - [in/out]? -
/// @param  rubbish - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_hand(
	int const k,
	int const hm,
	FArray1A_int stpppt,
	FArray1A_int strtpt,
	int const str1,
	int const str2,
	FArray1A_int order,
	FArray2A_float strdr,
	int const nstr,
	int const nres,
	FArray1A_int sequence,
	FArray1A_int directions,
	FArray1A_int scstr,
	FArray2A_float lctn,
	bool const use_whole_helix,
	int & rubbish
) const
{
	stpppt.dimension( nstr );
	strtpt.dimension( nstr );
	order.dimension( hm );
	strdr.dimension( nstr, 3 );
	sequence.dimension( hm );
	directions.dimension( hm );
	scstr.dimension( nres );
	lctn.dimension( nres, 3 );

	//---- local:
	int cnt, /*hel_ib,*/ i, j, ll, str_ib;
	float dot_prod;
	//      float dpr,ln1,ln2;
	FArray1D_float midpoint( 3 );
	FArray1D_float midpoint_loop( 3 );
	FArray1D_float vector0( 3 );
	FArray1D_float vector1( 3 );
	FArray1D_float vector2( 3 );
	FArray1D_float vector3( 3 );
	//------------------------------------------------------------------------------


	//  -  loop length between strands has to be at least 5 residues
	ll = strtpt(str2) - stpppt(str1) - 1;


	if ( ll < 5 ) return;

	//     std::cout << ll << std::endl;


	//  -  strands must be parallel
	//      std::cout << SS( shnm ) << SS( k ) << SS( directions(sequence(k-1)) ) <<
	//       SS( directions(sequence(k)) ) << SS( sequence(k-1) ) <<
	//       SS( sequence(k) ) << std::endl;
	//      if ( directions(sequence(k-1))/=directions(sequence(k)) ) return;
	if ( directions(k-1) != directions(k) ) return;

	//      if ( directions(order(k-1))/=directions(order(k)) ) return;

	//  -  straight dot-product between the two strands
	//     did not work well to define handedness
	//      goto L2222:
	//      dpr = strdr(str1,1)*strdr(str2,1)+strdr(str1,2)*strdr(str2,2)+
	//       strdr(str1,3)*strdr(str2,3);
	//      ln1 = std::sqrt(square( strdr(str1,1) )+square( strdr(str1,2) )+square( strdr(str1,3) ));
	//      ln2 = std::sqrt(square( strdr(str2,1) )+square( strdr(str2,2) )+square( strdr(str2,3) ));
	//      dpr /= (ln1*ln2);
	// L2222:;

	//  -  check for helix or strands from other sheets in between strands
	//     must be consecutive strands
	//hel_ib = 0; // js helix in between
	str_ib = 0; // js strand in between
	for ( j = stpppt(str1) + 1; j <= strtpt(str2) - 1; ++j ) {
		//if ( scstr(j) == 1 ) hel_ib = 1;  // set but never used ~Labonte
		if ( scstr(j) == 2 ) str_ib = 1;
	}
	if ( str_ib == 1 ) return;

	//js finding position of the intervening segment

	if ( !use_whole_helix ) {
		//default behavior is to just look at regions near strands ("takeoff" points)
		for ( i = 1; i <= 3; ++i ) {
			midpoint_loop(i) = 0.0;
			cnt = 0;
			for ( j = stpppt(str1) + 1; j <= stpppt(str1) + 3; ++j ) {
				//js first 3 residues after first strand
				midpoint_loop(i) += lctn(j,i);
				++cnt;
			}
			//js last 3 residues before the second strand
			for ( j = strtpt(str2) - 3; j < strtpt(str2); ++j ) {
				midpoint_loop(i) += lctn(j,i);
				++cnt;
			}
			midpoint_loop(i) /= cnt;
		}
	} else {
		for ( i = 1; i <= 3; ++i ) {
			midpoint_loop(i) = 0.0;
			cnt = 0;
			// Look at whole intervening loop (helix).
			for ( j = stpppt(str1) + 1; j <= stpppt(str2) - 1; ++j ) {
				midpoint_loop(i) += lctn(j,i);
				++cnt;
			}
			midpoint_loop(i) /= cnt;
		}
	}

	for ( i = 1; i <= 3; ++i ) {
		//js  last residue of the first strand,
		//js  first residue of the second strand
		midpoint(i) = lctn(strtpt(str1),i) + lctn(strtpt(str2),i) +
			lctn(stpppt(str1),i) + lctn(stpppt(str2),i);
		midpoint(i) /= 4.0;
		//js vector0 is the difference between the center of mass of
		//js the CA postions of 3 residues at either end of the loop (midpoint_loop)
		//js and the center of mass of the CA positions of the start and stop
		//js residues of each strand.
		vector0(i) = midpoint_loop(i)-midpoint(i);
	}

	for ( i = 1; i <= 3; ++i ) {
		vector1(i) = lctn(stpppt(str1),i) + lctn(stpppt(str2),i) -
			lctn(strtpt(str1),i) - lctn(strtpt(str2),i);
		vector1(i) /= 2.0;
		vector2(i) = lctn(strtpt(str2),i) + lctn(stpppt(str2),i) -
			lctn(strtpt(str1),i) - lctn(stpppt(str1),i);
		vector2(i) /= 2.0;
	}

	vector3(1) = vector2(2)*vector1(3)-vector2(3)*vector1(2);
	vector3(2) = vector2(3)*vector1(1)-vector2(1)*vector1(3);
	vector3(3) = vector2(1)*vector1(2)-vector2(2)*vector1(1);

	dot_prod = vector0(1)*vector3(1)+vector0(2)*vector3(2)+vector0(3)*vector3(3);

	// hand = 0;
	// if ( dot_prod > 0 ) hand = 1;
	if ( dot_prod < 0 ) rubbish = 6;

	if ( dot_prod < 0 && use_whole_helix ) {
		rubbish = 9;
		handedness_score_ += 0.1*std::fabs(dot_prod); // Arbitrary penalty for wrong handedness.
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM IDENTIFIES THE SHEETS
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  dstrmin - [in/out]? -
/// @param  ngbhct - [in/out]? -
/// @param  strprs - [in/out]? -
/// @param  nsht - [in/out]? -
/// @param  strsht - [in/out]? -
/// @param  strlbl - [in/out]? -
/// @param  rubbish - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_ident_sheets(
	int & nstr,
	FArray2A_float dstrmin,
	float ngbhct,
	FArray2A_int strprs,
	int & nsht,
	FArray1A_int strsht,
	FArray1A_int strlbl,
	int & rubbish
) const
{

	dstrmin.dimension( nstr, nstr );
	strprs.dimension( nstr, nstr );
	strsht.dimension( nstr );
	strlbl.dimension( nstr );

	//---- local:
	int cnt,lb1,lb2,lbmn,memo,strsm;//,smm;
	FArray1D_int fnlbl( max_nstr, 0 );
	FArray1D_int fncnt( max_nstr, 0 );
	FArray1D_int strlblord( nstr, 0 );
	//------------------------------------------------------------------------------

	strprs = 0;

	for ( int i = 1; i <= nstr; ++i ) {
		strlbl(i) = i;
		for ( int j = 1; j <= nstr; ++j ) {
			if ( dstrmin(i,j) < ngbhct ) strprs(i,j) = 1;
		}
	}

	rubbish = 0;
	for ( int i = 1; i <= nstr; ++i ) {
		//int smm = 0;
		//for ( int j = 1; j <= nstr; ++j ) {
		// smm += strprs(i,j);
		//}
		//js commenting out filter number 2.  It checks whether strands
		//js have more than 3 neighbors.  This should not be a filter,
		//js because native proteins can have strands with 3 neighbors.
		//js *Dimers* shouldn't be able to have more than 2 strands as
		//js neighbors.  That would require more work to fix.  For now
		//js we just comment out this filter.
		//js         if ( smm > 2 ) {
		//js            rubbish = 2;
		//js            goto L1234;
		//js         }
	}


	//      rubbish = 0
	//      for ( int i = 1; i <= nstr; ++i ) {
	//         if ( SUM(strprs(i,1:nstr)) > 2 ) {
	//            rubbish = 2;
	//            goto L1234;
	//         }
	//      }

	nsht = nstr;
	for ( int j = 1; j <= nstr; ++j ) {
		strsht(j) = 1;
	}

	for ( int i = 1; i <= (nstr-1); ++i ) {
		for ( int j = (i+1); j <= nstr; ++j ) {
			if ( strprs(i,j) == 1 ) {
				if ( strlbl(i) != strlbl(j) ) {
					--nsht;
					strsm = strsht(i)+strsht(j);
					lb1 = strlbl(i);
					lb2 = strlbl(j);
					lbmn = std::min(strlbl(i),strlbl(j));
					for ( int k = 1; k <= nstr; ++k ) {
						if ( strlbl(k) == lb1 || strlbl(k) == lb2 ) {
							strlbl(k) = lbmn;
							strsht(k) = strsm;
						}
					}
				}
			}
		}
	}


	// --- SORTING

	for ( int j = 1; j <= nstr; ++j ) {
		fnlbl(j) = 0;
		fncnt(j) = 0;
	}
	for ( int j = 1; j <= nstr; ++j ) {
		fnlbl(strlbl(j)) = 1;
		++fncnt(strlbl(j));
	}

	cnt = 0;
	for ( int j = 1; j <= nstr; ++j ) {
		strsht(j) = 0;
		if ( fnlbl(j) == 1 ) {
			++cnt;
			strsht(cnt) = fncnt(j);
		}
	}

	cnt = 0;
	for ( int j = 1; j <= nstr; ++j ) {
		strlblord(j) = 0;
	}

	for ( int j = 1; j <= (nstr-1); ++j ) {
		if ( (strlbl(j) > 0) && (strlblord(j) == 0) ) {
			++cnt;
			strlblord(j) = cnt;
			for ( int k = j+1; k <= nstr; ++k ) {
				if ( strlbl(j) == strlbl(k) ) {
					strlblord(k) = cnt;
				}
			}
		}
	}

	memo = 0;
	if ( nsht > 1 ) {
		cnt = 1;
		while ( cnt < nsht ) {
			for ( int j = 1; j <= nstr; ++j ) {
				if ( strlbl(j) > cnt ) {
					memo = strlbl(j);
					break;
				}
			}
			++cnt;
			for ( int j = 1; j <= nstr; ++j ) {
				if ( strlbl(j) == memo ) strlbl(j) = cnt;
			}
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM DETERMINES THE NUMBER OF LOCALS AND NONLOCAL
///     STRAND PAIRS
///
/// @details
///
/// @param  hm - [in/out]? -
/// @param  order - [in/out]? -
/// @param  nloc - [in/out]? -
/// @param  nnloc - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_lnl(
	int hm,
	FArray1A_int order,
	int & nloc,
	int & nnloc
) const
{
	order.dimension( hm );

	for ( int i = 2; i <= hm; ++i ) {
		if ( std::abs(order(i)-order(i-1)) > 1 ) {
			++nnloc;
		} else {
			++nloc;
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS FUNCTION DETERMINES THE CA LOCATION OF EACH RESIDUE
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  natm - [in/out]? -
/// @param  indCA - [in/out]? -
/// @param  atmps - [in/out]? -
/// @param  lctn - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_locations(
	int const & nres,
	int natm,
	FArray1D_int & indCA,
	FArray2D_float const & atmps,
	FArray2A_float lctn
) const
{
	indCA.dimension( natm );
	// atmps.dimension( 3, natm );
	lctn.dimension( nres, 3 );

	int cnt = 0;
	for ( int j = 1; j <= natm; ++j ) {
		if ( indCA(j) == 1 ) {
			++cnt;
			for ( int k = 1; k <= 3; ++k ) {
				lctn(cnt,k) = atmps(k,j);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS FUNCTION COUNTS THE NUBER OF STRANDS IN THE
///     PROTEIN/DECOY
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  scstr - [in/out]? -
/// @param  nstr - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_number_of_strands(
	int const & nres,
	FArray1A_int scstr,
	int & nstr
) const
{
	scstr.dimension( nres );

	nstr = 0;
	for ( int j = 1; j <= (nres-1); ++j ) {
		if ( scstr(j) == 2 && scstr(j+1) != 2 ) ++nstr;
	}
	if ( scstr(nres) == 2 ) ++nstr;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM CHECKS THE CONFIGURATION OF THE STRANDS FROM
///     THE SECONDARY STRUCTURE PREDICTION
///
/// @details
///
/// @param  nstr - [in/out]? -
/// @param  nsht - [in/out]? -
/// @param  ngbhct - [in/out]? -
/// @param  dotcut - [in/out]? -
/// @param  dstrmin - [in/out]? -
/// @param  strdr - [in/out]? -
/// @param  strsht - [in/out]? -
/// @param  strlbl - [in/out]? -
/// @param  rubbish - [in/out]? -
/// @param  maxdist - [in/out]? -
/// @param  mindotprodabs - [in/out]? -
/// @param  strtpt - [in/out]? -
/// @param  stpppt - [in/out]? -
/// @param  locdsm - [in/out]? -
/// @param  lctn - [in/out]? -
/// @param  nres - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_proper_sheets(
	int & nstr,
	int & nsht,
	float ngbhct,
	FArray2A_float dstrmin,
	FArray2A_float strdr,
	FArray1A_int strsht,
	FArray1A_int strlbl,
	int & rubbish,
	float & maxdist,
	float & mindotprodabs,
	FArray1A_int strtpt,
	FArray1A_int stpppt,
	FArray2A_int locdsm,
	FArray2A_float lctn,
	int const & nres
) const
{
	dstrmin.dimension( nstr, nstr );
	strdr.dimension( nstr, 3 );
	strsht.dimension( nstr );
	strlbl.dimension( nstr );
	strtpt.dimension( nstr );
	stpppt.dimension( nstr );
	locdsm.dimension( nstr, nstr );
	lctn.dimension( nres, 3 );

	//float dot_prod, ln1, ln2;
	//int pick1,pick2;
	FArray1D_float dir1( 3 );
	FArray1D_float dir2( 3 );

	//------------------------------------------------------------------------------

	rubbish = 0;

	//      std::cout << SS( nstr ) << SS( nsht ) << std::endl;
	if ( nstr < 2 ) rubbish = 3;

	//         std::cout << strsht(1:nsht) << std::endl;

	for ( int j = 1; j <= nsht; ++j ) {
		if ( strsht(j) < 2 ) rubbish = 4;
	}

	if ( rubbish == 0 ) {
		mindotprodabs = 1.0;
		maxdist = 0.0;

		for ( int i = 1; i <= nsht; ++i ) {
			int lbbl = i;
			for ( int j = 1; j <= (nstr-1); ++j ) {
				if ( strlbl(j) != lbbl ) continue;
				if ( stpppt(j)-strtpt(j) < 6 ) continue;
				for ( int k = (j+1); k <= nstr; ++k ) {
					if ( strlbl(k) != lbbl ) continue;
					if ( stpppt(k)-strtpt(k) < 6 ) continue;
					if ( dstrmin(j,k) < ngbhct ) {
						maxdist = std::max(dstrmin(j,k),maxdist);
						float dot_prod = 0.0;
						for ( int l = 1; l <= 3; ++l ) {
							//                 std::cout << SS( locdsm(j,k) ) << SS( nres ) << std::endl;
							//                 std::cout << SS( locdsm(k,j) ) << SS( nres ) << std::endl;
							int pick1 = locdsm(j,k);
							int pick2 = locdsm(k,j);
							if ( locdsm(j,k) <= 3 ) pick1 = 4;
							if ( locdsm(j,k) >= nres-3 ) pick1 = nres-4;
							if ( locdsm(k,j) <= 3 ) pick2 = 4;
							if ( locdsm(k,j) >= nres-3 ) pick2 = nres-4;
							dir1(l) = lctn(pick1+3,l)-lctn(pick1-3,l);
							dir2(l) = lctn(pick2+3,l)-lctn(pick2-3,l);
						}
						dot_prod = dir1(1)*dir2(1) + dir1(2)*dir2(2) + dir1(3)*dir2(3);
						float ln1 = std::sqrt(
							( dir1(1) * dir1(1) ) +
							( dir1(2) * dir1(2) ) +
							( dir1(3) * dir1(3) ) );
						float ln2 = std::sqrt(
							( dir2(1) * dir2(1) ) +
							( dir2(2) * dir2(2) ) +
							( dir2(3) * dir2(3) ) );
						dot_prod /= (ln1*ln2);
						if ( std::abs(dot_prod) < mindotprodabs ) {
							mindotprodabs = std::abs(dot_prod);
						}
						//                 std::cout << SS( dot_prod ) << std::endl;
						//                 if ( dot_prod > -dotcut && dot_prod < dotcut ) rubbish = 5;
					}
				}
			}
		}

	}

	//      std::cout << "rub " << rubbish << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS FUNCTION SPITS OUT THE BETA SHEET INFORMATION
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  scstr - [in/out]? -
/// @param  natm - [in/out]? -
/// @param  atmps - [in/out]? -
/// @param  indN - [in/out]? -
/// @param  indCA - [in/out]? -
/// @param  indC - [in/out]? -
/// @param  rnm - [in/out]? -
/// @param  ngbhct - [in/out]? -
/// @param  dotcut - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  nsht - [in/out]? -
/// @param  maxdt - [in/out]? -
/// @param  mindp - [in/out]? -
/// @param  nloc - [in/out]? -
/// @param  nnloc - [in/out]? -
/// @param  rubbish - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_sheet_stuff(
	int const & nres,
	FArray1A_int scstr,
	int natm,
	FArray2D_float const & atmps,
	FArray1D_int & indN,
	FArray1D_int & indCA,
	FArray1D_int & indC,
	FArray1D_int & rnm,
	float ngbhct,
	int & nstr,
	int & nsht,
	float & maxdt,
	float & mindp,
	int & nloc,
	int & nnloc,
	int & rubbish
) const
{
	//using namespace param;

	scstr.dimension( nres );
	// atmps.dimension( 3, natm );
	indN.dimension( natm );
	indCA.dimension( natm );
	indC.dimension( natm );
	rnm.dimension( natm );

	//------------------------------------------------------------------------------
	//     atmps       array of atom positions
	//     directions  directions of the strands in a sheet
	//     dotcut      cutoff for minimum dotproduct in proper sheets
	//     dstrmin     matrix of (minimum diamer) distances of strands
	//     ind*        vector of indicators (0/1) for type of atom
	//     inddm       vector of indicators (0/1) for strand diamers
	//     lctn        array of residue (CA) positions
	//     locdsm      matrix of locations of (minimum diamer) distances of strands
	//     maxdt       maximum distance between two neighbour strands
	//     mindp       minimum dotproduct between two neighbour strands
	//     natm        number of atoms in the protein
	//     ngbhct      cutoff for strand distances that define neighbours
	//     nloc        number of local strand pairs (one sheet proteins only)
	//     nmamac      name of amino acid (single letter abbrev)
	//     nmatom      name of atom (standard abbrev)
	//     nnloc       number of non-local strand pairs (one sheet proteins only)
	//     nres        number of residues in the protein
	//     nsht        number of sheets in the protein
	//     nstr        number of strands in the protein
	//     order       order of the strands in a sheet
	//     proper      indicator sheet violations (1 no, 0 yes)
	//     rnm         vector of residue numbers of atoms
	//     rubbish     indicator (0 good structure; bigger 0 not)
	//     scstr       vector of secondary structure information
	//                 (1 helix, 2 strand, 3 turn)
	//     slct        array of strand numbers (in sequence) that are
	//                 involved in the sheet under consideration
	//     sequence    sequence of the strands in a sheet
	//     stpppt      vector of positions of strand ends
	//     strdm       array of diamer positions
	//     strdr       array of strand directions (3D)
	//     strlbl      labels of strands, indicating which sheet they belong to
	//     strnm       vector of strand number of residues
	//     strprs      matrix indicating which strands are neighbours (1 yes, 0 no)
	//     strsht      vector with numbers of strands in sheets
	//                 (is of length nstr for simplicity, but actually only the
	//                 first nsh integers are considered, rest is 0)
	//     strtpt      vector of positions of strand starts
	//     indN        vector of 0/1 for each atom in atmps  0=not N, 1=is N
	//     indCA,indC    like indN

	//     rubbish     filter indicator, note that if a structure fails at
	//                 one level, higher level filters are never evaluated
	//$$$  0: passed
	//$$$  1: no strands
	//$$$  2: a strand with more than 2 neighbours
	//$$$  3: a decoy with only one strand
	//$$$  4: a one-stranded sheet
	//$$$  5: bad dotproduct    // no longer a possible result
	//$$$  6: left handed connection between parallel neighbours
	//$$$  7: barrel type sheet
	//------------------------------------------------------------------------------

	nstr = 0;
	nsht = 0;
	maxdt = 0;
	mindp = 0;
	rubbish = 0;
	nloc = 0;
	nnloc = 0;

	ingo_clean_ss(nres,scstr);
	ingo_number_of_strands(nres,scstr,nstr);
	if ( nstr == 0 ) {
		rubbish = 1;
	} else if ( nstr == 1 ) {
		rubbish = 3;
	} else {
		FArray1D_int inddm( nres, 0 );
		FArray1D_int strnm( nres, 0 );
		FArray1D_int directions( nstr );
		FArray1D_int order( nstr );
		FArray1D_int slct( nstr );
		FArray1D_int sequence( nstr );
		FArray1D_int strlbl( nstr, 0 );
		FArray1D_int strsht( nstr, 0 );
		FArray1D_int strtpt( nstr, 0 );
		FArray1D_int stpppt( nstr, 0 );
		FArray2D_int locdsm( nstr, nstr, 0 );
		FArray2D_int strprs( nstr, nstr, 0 );
		FArray2D_float dstrmin( nstr, nstr, 0.0 );
		FArray2D_float lctn( nres, 3, 0.0 );
		FArray2D_float strdm( nres, 3, 0.0 );
		FArray2D_float strdr( nstr, 3, 0.0 );

		ingo_start_stop(nres,scstr,nstr,strnm,strtpt,stpppt);
		ingo_diamers(nres,strnm,natm,atmps,rnm,indC,indN,strdm,inddm);
		ingo_locations(nres,natm,indCA,atmps,lctn);
		ingo_strand_dirs(nres,nstr,strtpt,stpppt,lctn,strdr);
		ingo_strand_dists_min(nres,nstr,inddm,strnm,strdm,dstrmin,locdsm);
		ingo_ident_sheets(nstr,dstrmin,ngbhct,strprs,nsht,strsht,strlbl,
			rubbish);
		//rhiju  Rather than escape, want handedness to always be evaluated...
		//if ( rubbish > 0 ) goto L1234;
		ingo_proper_sheets(nstr,nsht,ngbhct,dstrmin,strdr,strsht,strlbl,
			rubbish,maxdt,mindp,strtpt,stpppt,locdsm,lctn,nres);
		//rhiju  Rather than escape, want handedness to always be evaluated...
		//if ( rubbish > 0 ) goto L1234;

		for ( int j = 1; j <= nsht; ++j ) {
			//js possible bug causing spurious rejections of native structures with
			//js code 7: sends strsht (an array) but expects hm, an integer.
			//js logically it looks like it should be strsh(j), the number of strands
			//js in sheet j.
			//js I am adding index (j) to strsht:
			int const hm = strsht(j);
			for ( int k = 1; k <= hm; ++k ) {
				order(k) = 0;
				sequence(k) = 0;
				slct(k) = 0;
				directions(k) = 0;
			}
			ingo_find_ord(j,hm,nstr,strlbl,dstrmin,ngbhct,order,sequence,slct,
				rubbish);
			//rhiju ingo_find_ord doesn't assign strand order and directions if the sheet
			//rhiju forms a barrel. this might be worth fixing in the future, since we'd like to check
			//rhiju handedness of BAB's in the barrel.
			if ( rubbish > 6 ) break;
			ingo_lnl(hm,order,nloc,nnloc);
			ingo_find_dir(j,hm,nstr,slct,order,strlbl,strdr,directions);
			bool use_whole_helix (false);
			for ( int k = 2; k <= hm; ++k ) {
				ingo_hand(k,hm,stpppt,strtpt,slct(k-1),slct(k),order,strdr,
					nstr,nres,sequence,directions,scstr,lctn,use_whole_helix,rubbish);
			}
			// Following is new; slightly different computation of handedness, where
			// the whole body of the helix (rather than just the regions closest to
			// the strand takeoff points) is used to calculate beta-alpha-beta handedness.
			// For native PDBs, the usual definition is OK, but need a more robust definition
			// for Rosetta decoys (especially those encountered during jumping).
			use_whole_helix = true;
			for ( int k = 2; k <= hm; ++k ) {
				ingo_hand(k,hm,stpppt,strtpt,slct(k-1),slct(k),order,strdr,
					nstr,nres,sequence,directions,scstr,lctn,use_whole_helix,rubbish);
			}
		}
		//L1234:;
	}

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS FUNCTION IDENTIFIES THE START AND END POINTS OF THE
///     STRANDS
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  scstr - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  strnm - [in/out]? -
/// @param  strtpt - [in/out]? -
/// @param  stpppt - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_start_stop(
	int const & nres,
	FArray1A_int scstr,
	int & nstr,
	FArray1A_int strnm,
	FArray1A_int strtpt,
	FArray1A_int stpppt
) const
{
	scstr.dimension( nres );
	strnm.dimension( nres );
	strtpt.dimension( nstr );
	stpppt.dimension( nstr );

	int cnt = 0;
	int start = 0;
	int stopp = 0;
	if ( scstr(1) == 2 ) start = 1;
	int j = 0;
	while ( j < (nres-1) ) {
		++j;
		strnm(j) = 0;
		if ( scstr(j) == 2 ) {
			if ( start == 0 ) start = j;
			strnm(j) = cnt+1;
			if ( scstr(j) != scstr(j+1) ) {
				++cnt;
				stopp = j;
				strtpt(cnt) = start;
				stpppt(cnt) = stopp;
				start = 0;
				stopp = 0;
			}
		}
	}

	if ( scstr(nres) == 2 ) {
		++cnt;
		strtpt(cnt) = start;
		stpppt(cnt) = nres;
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM DETERMINES THE STRAND DIRECTIONS
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  strtpt - [in/out]? -
/// @param  stpppt - [in/out]? -
/// @param  lctn - [in/out]? -
/// @param  strdr - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_strand_dirs(
	int const & nres,
	int & nstr,
	FArray1A_int strtpt,
	FArray1A_int stpppt,
	FArray2A_float lctn,
	FArray2A_float strdr
) const
{
	strtpt.dimension( nstr );
	stpppt.dimension( nstr );
	lctn.dimension( nres, 3 );
	strdr.dimension( nstr, 3 );

	for ( int i = 1; i <= nstr; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			strdr(i,j) = lctn(stpppt(i),j)-lctn(strtpt(i),j);
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM DETERMINES THE MINIMUM STRAND DISTANCES
///
/// @details
///
/// @param  nres - [in/out]? -
/// @param  nstr - [in/out]? -
/// @param  inddm - [in/out]? -
/// @param  strnm - [in/out]? -
/// @param  strdm - [in/out]? -
/// @param  dstrmin - [in/out]? -
/// @param  locdsm - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_strand_dists_min(
	int const & nres,
	int & nstr,
	FArray1A_int inddm,
	FArray1A_int strnm,
	FArray2A_float strdm,
	FArray2A_float dstrmin,
	FArray2A_int locdsm
) const
{
	inddm.dimension( nres );
	strnm.dimension( nres );
	strdm.dimension( nres, 3 );
	dstrmin.dimension( nstr, nstr );
	locdsm.dimension( nstr, nstr );

	for ( int i = 1; i <= nstr; ++i ) {
		for ( int j = 1; j <= nstr; ++j ) {
			locdsm(i,j) = 0;
			dstrmin(i,j) = 10000.0;
		}
	}

	for ( int i = 1; i <= nres; ++i ) {
		if ( inddm(i) > 0 ) {
			float const si1 = strdm(i,1);
			float const si2 = strdm(i,2);
			float const si3 = strdm(i,3);
			for ( int j = 1; j <= nres; ++j ) {
				if ( (inddm(j) > 0) && (strnm(j) != strnm(i)) ) {
					float const sd1 = si1 - strdm(j,1);
					float const sd2 = si2 - strdm(j,2);
					float const sd3 = si3 - strdm(j,3);
					float const dist =
						std::sqrt( ( sd1 * sd1 ) + ( sd2 * sd2 ) + ( sd3 * sd3 ) );
					if ( dist < dstrmin(strnm(i),strnm(j)) ) {
						dstrmin(strnm(i),strnm(j)) = dist;
						locdsm(strnm(i),strnm(j)) = i;
						locdsm(strnm(j),strnm(i)) = j;
					}
				}
			}
		}
	}

	// std::cout << std::endl;
	// for ( int i = 1; i <= nstr; ++i ) {
	//  for ( int j = 1; j <= nstr; ++j ) {
	//   std::cout << F( 9, 2, dstrmin(i,j) );
	//  } std::cout << std::endl;
	// }
	// std::cout << std::endl;
	// for ( int i = 1; i <= nstr; ++i ) {
	//  for ( int j = 1; j <= nstr; ++j ) {
	//   std::cout << I( 5, locdsm(i,j) );
	//  } std::cout << std::endl;
	// }
	// std::cout << std::endl;

}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///     THIS PROGRAM CLEANS THE SECONDARY STRUCTURE OF DECOYS
///
/// @details
///car notes from ingo
///
///$$$ If I understand my notes
///$$$correctly, the big trouble makers were mostly the single strand residue
///$$$that popped up quite frequently. So they are cleaned up. Also, if a strand
///$$$is interrupted by one non-strand residue, that gets turned into a strand
///$$$residue too, so two short strands get turned into one long strand. That
///$$$worked the best, according to my notes.
///
///$$$The secstr cleanup doesn't contain one of the issues that David had
///$$$mentioned: if you have strands like that:
///$$$
///$$$   xxxxxxxx  xxxxxxxx
///$$$   xxxxxxxxxxxxxxxxxx
///$$$
///$$$That is, there's a hole in one strand and it pops up as two strands,
///$$$paired up with another strand. A single residue missing gets cought by the
///$$$secstr cleaner, so I think most cases are taken care of in the beginning.
///$$$There should be something in there as well that catches two or three
///$$$residues missing. But that couldn't occur of course before you identified
///$$$sheets. I didn't see anything in the ident_sheets function (the most
///$$$obvious place where something like that could be).
///
/// @param  nres - [in/out]? -
/// @param  scstr - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
SheetFilter::ingo_clean_ss(
	int const & nres,
	FArray1A_int scstr
) const
{
	scstr.dimension( nres );

	// --- first and last residue strand check
	if ( scstr(nres) == 2 && scstr(nres-1) != 2 ) {
		scstr(nres) = 3;
	}
	if ( scstr(1) == 2 && scstr(2) != 2 ) {
		scstr(1) = 3;
	}

	// --- single residue strand checks
	for ( int j = 2; j <= (nres-1); ++j ) {
		if ( scstr(j-1) == 2 && scstr(j) != 2 && scstr(j+1) == 2 ) {
			scstr(j) = 2;
		}
	}
	for ( int j = 2; j <= (nres-1); ++j ) {
		if ( scstr(j-1) != 2 && scstr(j) == 2 && scstr(j+1) != 2 ) {
			scstr(j) = scstr(j-1);
		}
	}
}

} // filters
} // protocols
