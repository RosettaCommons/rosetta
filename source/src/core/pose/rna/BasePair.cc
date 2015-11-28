// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/rna/BasePair.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/rna/BasePair.hh>
#include <core/pose/rna/leontis_westhof_util.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.rna.BasePair" );

namespace core {
namespace pose {
namespace rna {

	/////////////////////////////////////////////////////////////////////////
  BasePair::BasePair( Size const res1, Size const res2,
		      BaseEdge const edge1, BaseEdge const edge2,
		      BaseDoubletOrientation const orientation ):
    res1_( res1 ),
    res2_( res2 ),
    edge1_( edge1 ),
    edge2_( edge2 ),
    orientation_( orientation ),
    LW_orientation_( ANY_LW_BASE_DOUBLET_ORIENTATION )
  {
    derive_LW_orientation();
  }

	/////////////////////////////////////////////////////////////////////////
  BasePair
  BasePair::flipped() const
  {
    return BasePair( res2_, res1_, edge2_, edge1_, orientation_ );
  }

	/////////////////////////////////////////////////////////////////////////
  void
  BasePair::derive_LW_orientation() {
    if ( ( edge1_ != WATSON_CRICK && edge1_ != HOOGSTEEN && edge1_ != SUGAR ) ||
				 ( edge2_ != WATSON_CRICK && edge2_ != HOOGSTEEN && edge2_ != SUGAR ) ||
				 ( orientation_ != ANTIPARALLEL && orientation_ != PARALLEL ) ) {
      LW_orientation_ = ANY_LW_BASE_DOUBLET_ORIENTATION;
      return;
    }
    LW_orientation_ = get_LW_orientation( edge1_, edge2_, orientation_ );
  }

	/////////////////////////////////////////////////////////////////////////
  void
  BasePair::print_info( std::ostream & out ) const {
    out << "res1 = "  << std::setw( 4 ) << res1_;
    out << " res2 = " << std::setw( 4 ) << res2_;
    out << " edge1 =\t " << std::setw( 6 ) << core::chemical::rna::get_full_edge_from_num( edge1_ );
    out << " edge2 =\t " << std::setw( 6 ) << core::chemical::rna::get_full_edge_from_num( edge2_ );
    out << " LW_orient = " << std::setw( 5 ) << core::chemical::rna::get_full_LW_orientation_from_num( LW_orientation_ );
    out << " orient = " << std::setw( 5 ) << core::chemical::rna::get_full_orientation_from_num( orientation_ );
    out << " ";
  }

  ///////////////////////////////////////////////////////////////////
  bool operator < ( BasePair const & lhs, BasePair const & rhs )
  {
    //There must be a more elegant way to do this...
    if ( lhs.res1_ < rhs.res1_ ) {
      return true;
    } else if ( lhs.res1_ == rhs.res1_ ) {
      if ( lhs.res2_ < rhs.res2_ ) {
	return true;
      } else if ( lhs.res2_ == rhs.res2_ ) {
	if ( lhs.edge1_ < rhs.edge1_ ) {
	  return true;
	} else if ( lhs.edge1_ == rhs.edge1_ ) {
	  if ( lhs.edge2_ < rhs.edge2_ ) {
	    return true;
	  } else if ( lhs.edge2_ == rhs.edge2_ ) {
	    return ( lhs.orientation_ < rhs.orientation_ );
	  }
	}
      }
    }
    return false;
  }


  ///////////////////////////////////////////////////////////////////
  bool operator == ( BasePair const & lhs, BasePair const & rhs )
  {
    return ( lhs.res1_ == rhs.res1_ &&
	     lhs.res2_ == rhs.res2_ &&
	     lhs.edge1_ == rhs.edge1_ &&
	     lhs.edge2_ == rhs.edge2_ &&
	     lhs.orientation_ == rhs.orientation_ );
  }


  ///////////////////////////////////////////////////////////////////
  std::ostream &
  operator << ( std::ostream & out, BasePair const & s ){
    out << s.res1_ << " " << s.res2_ << " " << s.edge1_ << " " << s.edge2_ << " " << s.orientation_;
    return out;
  }

} //rna
} //pose
} //core
