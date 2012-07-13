// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// Copyright in the Rosetta software belongs to the developers and their institutions.
// For more information, see www.rosettacommons.org.

/// @file ./src/protocols/fldsgn/topology/Turn.hh
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@u.washington.edu )

#ifndef INCLUDED_protocols_fldsgn_topology_Turn_hh
#define INCLUDED_protocols_fldsgn_topology_Turn_hh

// unit headers
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/Turn.fwd.hh>

// project headers
#include <core/types.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

class Turn : public utility::pointer::ReferenceCount {
public:

	typedef std::string String;

	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef utility::vector1< Size > VecSize;
	typedef utility::vector1< int > VecInt;
	typedef utility::vector1< Real > VecReal;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;
	typedef protocols::fldsgn::topology::StrandPairing StrandPairing;
	

public:// construct/destruct


	/// @brief default constructor
	Turn();

	/// @brief value constructor
	Turn( SS_Info2_COP const ssinfo, StrandPairing const & sp );

	/// @brief copy constructor
	Turn( Turn const & s );

	/// @brief default destructor
	~Turn();

	/// @brief clone this object
	TurnOP clone() {
		return new Turn( *this );
	}

	/// @brief return strand pairing
	friend
	std::ostream & operator<<( std::ostream & out, const Turn &s );


public:


	/// @brief intialize this class
	void initialize( SS_Info2_COP const ssinfo, StrandPairing const & sp );
	

public: //accessors


	/// @brief turn type
	inline String type() const { return type_; }
	
	/// @brief turn length
	inline Size length() const { return length_;	}
		
	/// @brief residue number of the beginning of turn
	inline Size begin() const { return begin_;	}

	/// @brief residue number of the end of turn
	inline Size end() const { return end_;	}

	/// @brief pleat of anchor residue
	inline Size pleat() const { return pleat_;	}
	
	/// @brief string of abego
	inline String abego() const { return abego_; }
	
	
public: //
	
	/// @brief get info of this
	String get_info() const;
	
private: //
	
	/// @brief calc pleat
	Size calc_pleat( SS_Info2_COP const ssinfo, Size const r1, Size const r2 ) const;

private:  // data

	/// @brief turn type
	String type_;
	
	/// @brief turn length
	Size length_;
	
	/// @brief residue number of the beginning of turn
	Size begin_;

	/// @brief residue number of the end of turn
	Size end_;

	/// @brief pleat of anchor residue
	Size pleat_;

	/// @brief string of abego
	String abego_;

};


class TurnSet : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef utility::vector1< Size > VecSize;
	typedef protocols::fldsgn::topology::StrandPairingSet StrandPairingSet;
	typedef protocols::fldsgn::topology::StrandPairingSetCOP StrandPairingSetCOP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;

 
public:// construct/destruct


	/// @brief default constructor
	TurnSet();

	/// @brief value constructor
	TurnSet( Turns const & turns );

	/// @brief value constructor
	TurnSet( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset );

	/// @brief copy constructor
	TurnSet( TurnSet const & s );

	/// @brief return strand pairing
	friend
	std::ostream & operator<<( std::ostream & out, const TurnSet &s );


public:


	void initialize( SS_Info2_COP const ssinfo, StrandPairingSetCOP const spairset );


public: // mutators


	/// @brief
	void push_back( TurnOP const sop );

	/// @brief
	void clear();


public: // accessors

	
	/// @brief return a turn_
	TurnOP turn( Size const s ) const;

	/// @brief return turns_
	inline Turns turns() const { return turns_; }
	
	/// @brief return max lengh of turn
	inline Size max_length_turn() const { return max_length_turn_; }
	
	/// @brief return 
	TurnOP sspair_turn( Size const s1, Size const s2 ) const;

	/// @brief return number of Turns
	inline Size	size() const { return turns_.size(); }



private:  // data
	
	
	/// @brief maximal length of turn to define
	Size max_length_turn_;
	
	/// @brief 2D table of the pointers of turn, which are sorted by the strand number
	utility::vector1< utility::vector1< TurnOP > > map_sspair_turn_;

	/// @brief
	Turns turns_;


};

} // namespace topology
} // namespace fldsgn
} // namespace protocol

#endif
