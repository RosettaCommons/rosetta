// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//  CVS information:
//  $Revision: 1.1.2.1 $
//  $Date: 2005/11/07 21:05:35 $
//  $Author: rhiju $
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#ifndef INCLUDED_core_scoring_rna_RNA_BaseDoubletClasses_hh
#define INCLUDED_core_scoring_rna_RNA_BaseDoubletClasses_hh

#include <core/types.hh>
#include <core/chemical/rna/RNA_Util.hh>

// C++ Headers
#include <iomanip>
#include <iostream>
#include <list>

//using core::Size;
//using core::Real;

namespace core {
namespace scoring {
namespace rna {

/////////////////////////////////////////////////////////////////////
// Useful type definitions.
/////////////////////////////////////////////////////////////////////

class Base_pair
{

 public:

	Base_pair():
		res1( 0 ),
		res2( 0 ),
		edge1( 0 ),
		edge2( 0 ),
		orientation( 0 ),
		LW_orientation( 0 ), //Leontis Westhof base-pair orientation (1 = cis; 2 = trans). This is not yet implemented! (PS. 12/26/2011)
		num_hbonds( 0 )
	{
	};

	Base_pair( Size const & res1_input, Size const res2_input,
						 Size const & edge1_input, Size const edge2_input,
						 Size const & orientation_input ):
		res1( res1_input ),
		res2( res2_input ),
		edge1( edge1_input ),
		edge2( edge2_input ),
		orientation( orientation_input ),
		LW_orientation( 0 ),
		num_hbonds( 0 )
	{
	};

  Size res1;
  Size res2;
  Size edge1;
  Size edge2;
  Size orientation; // 1 = antiparallel; 2 = parallel
	Size LW_orientation; // 1 = cis; 2 = trans
	Size num_hbonds;


	void
	print_info( std::ostream & out = std::cout ) const {
		out << "res1 = "  << std::setw( 4 ) << res1;
		out << " res2 = " << std::setw( 4 ) << res2;
		out << " edge1 =	 " << std::setw( 6 ) << core::chemical::rna::get_full_edge_from_num( edge1 );
		out << " edge2 =	 " << std::setw( 6 ) << core::chemical::rna::get_full_edge_from_num( edge2 );
		out << " LW_orient = " << std::setw( 5 ) << core::chemical::rna::get_full_LW_orientation_from_num( LW_orientation );
		out << " orient = " << std::setw( 5 ) << core::chemical::rna::get_full_orientation_from_num( orientation );
		out << " #hbonds= " << std::setw(4) << num_hbonds;
		out << " ";
	}

  friend
    bool operator < ( Base_pair const & lhs, Base_pair const & rhs )
  {
    //There must be a more elegant way to do this...
    if ( lhs.res1 < rhs.res1 ) {
      return true;
		} else if ( lhs.res1 == rhs.res1 ) {
      if ( lhs.res2 < rhs.res2 ) {
				return true;
			} else if ( lhs.res2 == rhs.res2 ) {
				if ( lhs.edge1 < rhs.edge1 ) {
					return true;
				} else if ( lhs.edge1 == rhs.edge1 ) {
					if ( lhs.edge2 < rhs.edge2 ) {
						return true;
					}	else if ( lhs.edge2 == rhs.edge2 ) {
						return ( lhs.orientation < rhs.orientation );
					}
				}
			}
		}
		return false;
  };

  friend
    bool operator == ( Base_pair const & lhs, Base_pair const & rhs )
  {
		return ( lhs.res1 == rhs.res1 &&
						lhs.res2 == rhs.res2 &&
						lhs.edge1 == rhs.edge1 &&
						lhs.edge2 == rhs.edge2 &&
						lhs.orientation == rhs.orientation );
  };


  friend
    std::ostream &
    operator << ( std::ostream & out, Base_pair const & s ){
    out << s.res1 << " " << s.res2 << " " << s.edge1 << " " << s.edge2 << " " << s.orientation;
    return out;
  }
};

typedef std::pair< Real, Base_pair > Energy_base_pair;
typedef std::list < Energy_base_pair > Energy_base_pair_list;

class Base_stack{
 public:

	Base_stack():
		res1( 0 ),
		res2( 0 ),
		orientation( 0 ),
		which_side( 0 )
	{
	};


  Size res1;
  Size res2;
  Size orientation; // 1 = antiparallel; 2 = parallel
  Size which_side;  // 1 = residue 2 is 3' to residue1;  2 = residue 2 is 5' to residue 1

  friend
    bool operator < ( Base_stack const & lhs, Base_stack const & rhs ){
    //There must be a more elegant way to do this...
    if ( lhs.res1 < rhs.res1 ) {
      return true;
		}  else if ( lhs.res1 == rhs.res1 ) {
      if ( lhs.res2 < rhs.res2 ) {
				return true;
			} else if ( lhs.res2 == rhs.res2 ) {
				if ( lhs.orientation < rhs.orientation ) {
					return true;
				}	else if ( lhs.orientation == rhs.orientation ) {
					return ( lhs.which_side < rhs.which_side );
				}
			}
		}
		return false;
  }


  friend
    std::ostream &
    operator << ( std::ostream & out, Base_stack const & s )
    {
      out << s.res1 << " " << s.res2 << " " <<  s.orientation << " " << s.which_side;
      return out;
    }


};

typedef std::pair< Real, Base_stack > Energy_base_stack;
typedef std::list < Energy_base_stack > Energy_base_stack_list;

} //rna
} //scoring
} //core

#endif
