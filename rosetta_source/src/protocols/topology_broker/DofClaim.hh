// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file TopologyBroker
/// @brief  top-class (Organizer) of the TopologyBroker mechanism
/// @detailed responsibilities:
///           maintains list of ToplogyClaimers
///           maintains DofClaims -- exclusive or non-exclusively markedup dofs like BackboneClaim, IntraResClaim, JumpClaim
///           generates FoldTree, MoveMap, and collects samplers provided by TopologyClaimers
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_DofClaim_hh
#define INCLUDED_protocols_topology_broker_DofClaim_hh


// Unit Headers
#include <protocols/topology_broker/DofClaim.fwd.hh>


// Package Headers
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>


// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>

// ObjexxFCL Headers

// Utility headers
//#include <utility/io/izstream.hh>
//#include <utility/io/ozstream.hh>
//#include <utility/io/util.hh>
//#include <basic/Tracer.hh>
//#include <basic/options/option.hh>
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>

//#include <basic/options/option_macros.hh>

//// C++ headers
//#include <fstream>
#include <string>

// option key includes


namespace protocols {
namespace topology_broker {

/// A better DofClaims class would provide some extracting functions:
/// by owner
/// by type



class DofClaim : public utility::pointer::ReferenceCount {
public:
	typedef core::Size Size;
	enum ClaimType {
		BB,
		//	INTRA, //like side-chain doesn't propagate along FoldTree
		JUMP,
		CUT, //some way so that rigid-chunks can disallow cuts in their regions ?
		SEQUENCE, //I want control over residue type,
		ROOT //I want to set the Fold-Tree root
	};
	enum ClaimRight {
		NEED_TO_KNOW = 1,
		CAN_INIT,
		INIT,
		EXCLUSIVE, //MUST INIT
		REJECTED
	};
	DofClaim( TopologyClaimer* tc, ClaimRight right ) :
		claim_source_( tc ),
		right_( right ),
		approved_( false )
	{};

	virtual DofClaimOP clone() const = 0;
	//virtual?
	virtual Size size() const = 0;
	virtual Size pos( Size i ) const = 0;
	ClaimRight right() const { return right_; };

	TopologyClaimer const* owner() const { return claim_source_; }
	TopologyClaimer* owner() { return claim_source_; }

	virtual void toggle( core::kinematics::MoveMap&, bool /*new_setting*/ ) const {};

	virtual ClaimType type() const = 0;

	bool exclusive() const {
		return right() == DofClaim::EXCLUSIVE;
	}

	virtual std::string str_type() const = 0;
	virtual void show( std::ostream& os ) const;

	bool approved() const {
		return approved_;
	}

	void set_approved() { //this should only be called by the TopologyBroker ... make sure somehow? friend ?
		approved_ = true;
	}
private:
	TopologyClaimer* claim_source_; //NEVER Make this OP --- circularity in smart-pointers   (wanted this reference but there was some kind of problem ... what was it ?
	ClaimRight right_;
	bool approved_; //keep track of this ?
}; //class DofClaim

//this is a bit different then the other claims.
//overlapping sequences don't make sense I think.
// this is merely used to arrange multiple patches so that everybody has a residue number.
// pos defines the starting residue of the patch.
// length the length.
// the initial claim goes out with pos=0 if you allow this sequence to move.
// length specifies the number of residues.
// pos=0 claims will be assigned a position from the broker.
class SequenceClaim : public DofClaim {
	public:
	SequenceClaim(
			TopologyClaimer* tc,
			Size pos,
			Size length,
			std::string const& label,
			ClaimRight right
	)  : DofClaim( tc, right), pos_( pos ), length_( length ), label_( label ) {};

	virtual DofClaimOP clone() const { return new SequenceClaim( *this ); }

	virtual Size size() const { return 2; }

	virtual Size pos( Size i ) const {
		runtime_assert( i <= 2 );
		if ( i == 1 ) return pos_;
		if ( i == 2 ) return length_;
		return 0;
	};

	std::string const& label() const {
		return label_;
	}

	///@brief if you want to have a residue (eg., for a new ligand) you will be given a number...
	void set_offset( Size pos ) {
		pos_ = pos;
	}

	Size offset() const {
		return pos_;
	}

	Size last_residue() const {
		return pos_+length_-1;
	}

	virtual ClaimType type() const {
		return SEQUENCE;
	}

	virtual std::string str_type() const {
		return "SEQUENCE";
	}

	virtual void show( std::ostream& os ) const {
		os << " with label: " << label();
	};

protected:
	Size pos_;
	Size length_;
	std::string label_;
}; //class BBClaim

class BBClaim : public DofClaim {
public:
	BBClaim( TopologyClaimer* tc, Size pos, ClaimRight right = DofClaim::CAN_INIT ) : DofClaim( tc, right), pos_( pos ) {};

	virtual DofClaimOP clone() const { return new BBClaim( *this ); }

	virtual Size size() const { return 1; }

	virtual Size pos( Size /* i */ ) const {
		//		runtime_assert( i == 1 ); warning in release build
		return pos_;
	};

	virtual void toggle( core::kinematics::MoveMap& mm, bool new_setting ) const {
		mm.set_bb( pos_, new_setting );
	}

	virtual ClaimType type() const {
		return BB;
	}

	virtual std::string str_type() const {
		return "BB";
	}

protected:
	Size pos_;
}; //class BBClaim

class JumpClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:
	JumpClaim( TopologyClaimer* tc, core::Size pos1, core::Size pos2, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		pos1_( pos1 ),
		pos2_( pos2 ),
		atom1_( "" ),
		atom2_( "" )
	{}

	JumpClaim( TopologyClaimer* tc, core::Size pos1, core::Size pos2, std::string const& atom1, std::string const& atom2, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		permanent_( false ),
		pos1_( pos1 ),
		pos2_( pos2 ),
		atom1_( atom1 ),
		atom2_( atom2 )
	{}

	virtual DofClaimOP clone() const { return new JumpClaim( *this ); }

	virtual Size size() const { return 2; }

	virtual Size pos( Size i ) const {
		runtime_assert( i <= 2 && i > 0);
		if ( i==1 ) return pos1_;
		if ( i==2 ) return pos2_;
		return 0;
	}

	virtual void toggle( core::kinematics::MoveMap& mm, bool new_setting ) const {
		mm.set_jump( pos1_, pos2_, new_setting );
	}

	virtual ClaimType type() const {
		return JUMP;
	}

	virtual bool remove() const {
		return !permanent_;
	}

	virtual std::string str_type() const {
		return "JUMP";
	}

	std::string const& jump_atom( Size i ) const {
		runtime_assert( i <= 2 && i > 0 );
		if ( i == 1 ) return atom1_;
		if ( i == 2 ) return atom2_;
		return atom1_; //happy compiler, never reached
	}

	void set_jump_atom( core::Size i, std::string const& str ) {
		runtime_assert( i <= 2 && i > 0 );
		if ( i == 1 ) atom1_ = str;
		if ( i == 2 ) atom2_ = str;
	}

private:
	bool permanent_; //true if this jump should still be present after loop-closing
	Size pos1_;
	Size pos2_;
	std::string atom1_;
	std::string atom2_;
}; //JumpClaim


class CutClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:
	CutClaim( TopologyClaimer* tc, core::Size pos1, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		pos1_( pos1 )
	{}

	virtual DofClaimOP clone() const { return new CutClaim( *this ); }

	virtual Size size() const { return 1; }

	virtual Size pos( Size i ) const {
		runtime_assert( i <= 1 && i > 0);
		if ( i==1 ) return pos1_;
		return 0;
	}

	virtual ClaimType type() const {
		return CUT;
	}

	virtual bool remove() const {
		return false; //so far all CutClaims will be physical --> permanent cuts ... !permanent_;
	}

	virtual std::string str_type() const {
		return "CUT";
	}

private:
	//	bool permanent_; //true if this cut should still be present after loop-closing
	Size pos1_;
}; //CutClaim

class RootClaim : public DofClaim {
	//this class could also specify which atoms to use for the jumps --- but I never used this so far.... might be necessary for Zn jumps.
public:
	RootClaim( TopologyClaimer* tc, core::Size pos1, ClaimRight right = DofClaim::CAN_INIT ) :
		DofClaim( tc, right ),
		pos1_( pos1 )
	{}

	virtual DofClaimOP clone() const { return new RootClaim( *this ); }

	virtual Size size() const { return 1; }

	virtual Size pos( Size i ) const {
		runtime_assert( i <= 1 && i > 0);
		if ( i==1 ) return pos1_;
		return 0;
	}

	virtual ClaimType type() const {
		return ROOT;
	}

	virtual bool remove() const {
		return false; //so far all CutClaims will be physical --> permanent cuts ... !permanent_;
	}

	virtual std::string str_type() const {
		return "ROOT";
	}

private:
	//	bool permanent_; //true if this cut should still be present after loop-closing
	Size pos1_;
}; //CutClaim



extern std::ostream& operator<<( std::ostream& os, DofClaim const& );
extern std::ostream& operator<<( std::ostream& os, DofClaims const& );

}
}

#endif
