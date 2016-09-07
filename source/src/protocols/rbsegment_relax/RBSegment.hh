// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Frank DiMaio
/// @author Srivatsan Raman

#ifndef INCLUDED_protocols_rbsegment_relax_RBSegment_hh
#define INCLUDED_protocols_rbsegment_relax_RBSegment_hh
#include <protocols/rbsegment_relax/RBSegment.fwd.hh>

#include <core/types.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <iostream>

#ifdef WIN32
#include <functional>
#endif


namespace protocols {
namespace rbsegment_relax {

//////////////////////////////////////////////////////////
/// @brief Enumeration of RB types
/////////////////////////////////////////////////////////
enum RBSegmentType { RB_HELIX=1, RB_SHEET, RB_DEFAULT };


//////////////////////////////////////////////////////////
/// @brief RB residue range
/////////////////////////////////////////////////////////
class RBResidueRange {
public:
	RBResidueRange() = default;

	RBResidueRange( int begin, int end , RBSegmentType type=RB_DEFAULT ) :
		res_first( begin ),
		res_last( end ),
		seg_type( type )
	{}

	inline core::Size length() const { return res_last - res_first + 1; }
	inline core::Size start() const  { return res_first; }
	inline void set_start(core::Size S)  { res_first = S; }
	inline core::Size end() const    { return res_last; }
	inline void set_end(core::Size E)    { res_last = E; }
	inline RBSegmentType type() const { return (seg_type); }
	inline char char_type() const { return( seg_type==RB_HELIX ? 'H' : (seg_type==RB_SHEET ? 'E' : '-') ); }

private:
	core::Size res_first, res_last;
	RBSegmentType seg_type;
};


//////////////////////////////////////////////////////////
/// @brief Rigid-body segments in a protein
/////////////////////////////////////////////////////////
class RBSegment {
public:
	RBSegment(){
		sigAxisR_ = sigAxisT_ = sigOffAxisR_ = sigOffAxisT_ = 0.0;
	}

	/// construct a simple RB Segment
	RBSegment( int seg_begin, int seg_end, RBSegmentType seg_type ) {
		segments_.push_back( RBResidueRange( seg_begin, seg_end, seg_type ) );
		sigAxisR_ = sigAxisT_ = sigOffAxisR_ = sigOffAxisT_ = 0.0;
	}

	/// construct a simple RBSegment from an RB residue range
	RBSegment ( RBResidueRange const &range_in ) {
		segments_.push_back( range_in );
		sigAxisR_ = sigAxisT_ = sigOffAxisR_ = sigOffAxisT_ = 0.0;
	}

	/// construct a simple RB Segment
	RBSegment( int seg_begin, int seg_end, char type ) {
		RBSegmentType seg_type;
		if ( type == 'H' ) {
			seg_type = RB_HELIX;
		} else if ( type == 'E' || type == 'S' ) {
			seg_type = RB_SHEET;
		} else {
			seg_type = RB_DEFAULT;
		}
		segments_.push_back( RBResidueRange( seg_begin, seg_end, seg_type ) );
		sigAxisR_ = sigAxisT_ = sigOffAxisR_ = sigOffAxisT_ = 0.0;
	}

	/// construct a compound RBSegment from a vector of simple RBSegments
	RBSegment ( utility::vector1 < RBSegment > const &segs_in );

	void set_movement( core::Real sigAxisR, core::Real sigAxisT, core::Real sigOffAxisR=0.0, core::Real sigOffAxisT=0.0);
	void get_movement( core::Real &sigAxisR, core::Real &sigAxisT, core::Real &sigOffAxisR, core::Real &sigOffAxisT) const;

	// create a new RBsegment from this one by remapping residue ids
	RBSegment  remap( core::id::SequenceMapping const &mapping ) const;

	inline core::Size nContinuousSegments() const { return segments_.size(); }

	inline bool isEmpty() const     { return ( nContinuousSegments()==0 ); }
	inline bool isSimple() const    { return ( nContinuousSegments()==1 ); }
	inline bool isHelix() const     { return ( isSimple() && segments_[1].type()==RB_HELIX ); }
	inline bool isSheet() const     { return ( isSimple() && segments_[1].type()==RB_SHEET ); }
	inline bool isGenericRB() const { return ( isSimple() && segments_[1].type()==RB_DEFAULT ); }
	inline bool isCompound() const  { return ( nContinuousSegments()>1 ); }
	inline bool initialized() const { return (sigAxisR_==0&&sigAxisT_==0&&sigOffAxisR_==0&&sigOffAxisT_==0); }

	// accessor
	inline core::Real getSigAxisR() const { return ( sigAxisR_ ); }
	inline core::Real getSigAxisT() const { return ( sigAxisT_ ); }
	inline core::Real getSigOffAxisR() const { return ( sigOffAxisR_ ); }
	inline core::Real getSigOffAxisT() const { return ( sigOffAxisT_ ); }

	RBResidueRange & operator[](int i) { return segments_[i]; }
	RBResidueRange const & operator[](int i) const { return segments_[i]; }

private:
	utility::vector1< RBResidueRange > segments_;

	// rotational params
	core::Real sigAxisR_, sigAxisT_, sigOffAxisR_, sigOffAxisT_;
};

/////////////
// used to sort RB structs by start-res
/////////////
class RB_lt : public std::binary_function<double, double, bool> {
public:
	bool operator()(RBSegment x, RBSegment y) {
		if ( x.isEmpty() ) return true;
		else if ( y.isEmpty() ) return false;
		else return (x[1].start() < y[1].start());
	}
};

//////////////////////////////////////////////////////////
/// @brief Parses an RB segment file into a vector of RBsegments
/////////////////////////////////////////////////////////
void read_RBSegment_file(
	utility::vector1< RBSegment > &rbsegs,
	protocols::loops::Loops &loops,
	std::string filename,
	bool autoGenerateLoops=false,
	int nres=0, // only needed if  autoGenerateLoops==true
	utility::vector1< core::Size > cutpts=utility::vector1< core::Size >(0) // only needed if  autoGenerateLoops==true
);

void select_RBsegments(
	utility::vector1< RBSegment > const &rbsegs_in,
	protocols::loops::Loops const &loops_in,
	utility::vector1< RBSegment > &rbsegs_selected,
	protocols::loops::Loops &loops_selected
);


}
}

#endif
