// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file DnaChains.hh
/// @brief a descriptive but lightweight class to keep track of (usually basepaired) dna chains
/// @author ashworth

#ifndef INCLUDED_protocols_dna_DnaChains_hh
#define INCLUDED_protocols_dna_DnaChains_hh

#include <protocols/dna/DnaChains.fwd.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

#include <iosfwd>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace dna {

class DnaPosition {
/// @brief Stores residue index/indices for a DNA position, which could be a basepair or a single-stranded position.  Accomplishes what std::pair<Size,Size> would, but is nicer to use
public: // constructors
	// default empty constructor
	DnaPosition() : top_(0), bottom_(0), paired_(false) {}
	// construct a single-stranded position
	DnaPosition( core::Size i ) : top_(i), bottom_(0), paired_(false) {}
	// construct a double-stranded position (a basepair)
	DnaPosition( core::Size i, core::Size j )	: top_(i), bottom_(j), paired_(true) {}
	~DnaPosition(){}

public: // const methods
	core::Size top() const { return top_; }
	core::Size bottom() const { return bottom_; }
	bool paired() const { return paired_; }
	// there are no 'settors' -- just delete and reconstruct -- it is safer that way

private:
	core::Size top_;
	core::Size bottom_;
	bool paired_;
};

typedef std::map< core::Size, DnaPosition > DnaPositions;

class DnaChains : public utility::pointer::ReferenceCount {
// this class is a light wrapper for DnaPositions (a map which is typedefed above)
public: // constructors
	DnaChains();
	virtual ~DnaChains();
	DnaChains( DnaChains const & other );
	DnaChainsOP clone() const;

	void clear() { positions_.clear(); }
	bool empty() const { return positions_.empty(); }
	core::Size size() const { return positions_.size(); }

	DnaPosition & operator[] ( core::Size resindex ) { return positions_[ resindex ]; }
	DnaPosition const & operator[] ( core::Size resindex ) const
		{ return ( *( positions_.find( resindex ))).second; }

	DnaPositions & positions() { return positions_; }
	DnaPositions const & positions() const { return positions_; }

	DnaPositions::iterator begin() { return positions_.begin(); }
	DnaPositions::iterator end() { return positions_.end(); }
	DnaPositions::const_iterator begin() const { return positions_.begin(); }
	DnaPositions::const_iterator end() const { return positions_.end(); }

	// top strand resid check (looks in map keys)
	bool is_top( core::Size index ) const
	{ if ( positions_.count(index) != 0 ) { return true; } return false; }

	bool contains( core::Size index ) const;
	void print( core::pose::Pose const & pose, std::ostream & os ) const;

private:
	// map of DNA position info, keyed by Pose residue index
	DnaPositions positions_;
};

} // namespace dna
} // namespace protocols

#endif
