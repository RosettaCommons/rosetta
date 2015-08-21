// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Oliver Lange


#ifndef INCLUDED_protocols_topology_broker_claims_SequenceClaim_hh
#define INCLUDED_protocols_topology_broker_claims_SequenceClaim_hh


// Unit Headers
#include <protocols/topology_broker/claims/SequenceClaim.fwd.hh>

// Package Headers
#include <protocols/topology_broker/TopologyClaimer.hh>
#include <protocols/topology_broker/TopologyClaimer.fwd.hh>
#include <protocols/topology_broker/claims/DofClaim.hh>

// Project Headers
#include <core/types.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/sequence/AnnotatedSequence.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/exit.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

//// C++ headers
#include <string>
#include <sstream>


// option key includes


namespace protocols {
namespace topology_broker {
namespace claims {

// this is a bit different then the other claims.
// overlapping sequences don't make sense I think.
// this is merely used to arrange multiple patches so that everybody has a residue number.
// pos defines the starting residue of the patch.
// length the length.
// the initial claim goes out with pos=0 if you allow this sequence to move.
// length specifies the number of residues.
// pos=0 claims will be assigned a position from the broker.
class SequenceClaim : public DofClaim {

public:
	SequenceClaim(
		TopologyClaimerAP tc,
		std::string const& annotated_sequence,
		std::string const& label,
		core::Real priority = 0.0
	) : DofClaim( tc, NEED_TO_KNOW ),
		label_( label ),
		annotated_sequence_(annotated_sequence),
		priority_( priority )
	{}

private:
	// core::Size pos_;
	std::string label_;
	core::sequence::AnnotatedSequence annotated_sequence_;
	core::Real priority_; //lower first
	typedef DofClaim Parent;

public:
	virtual DofClaimOP clone() const { return DofClaimOP( new SequenceClaim( *this ) ); }

	core::Size length() const {
		return annotated_sequence_.length();
	}

	core::Real priority() const {
		return priority_;
	}

	core::sequence::AnnotatedSequence const& annotated_sequence() const {
		return annotated_sequence_;
	}

	std::string const& label() const {
		return label_;
	}

	virtual void show( std::ostream& os ) const {
		os << "SequenceClaim (" << label() << ", length="<<length()<<", priority="<<priority_<<") ";
		Parent::show( os );
	};

	virtual std::string str_type() const {
		return "SEQUENCE";
	}

};
// core::Size position() const {
//  return pos_;
// }

/// @brief if you want to have a residue (eg., for a new ligand) you will be given a number...
// void set_offset( core::Size pos ) {
//  pos_ = pos;
// }

// core::Size offset() const {
//  return pos_;
// }

//  core::Size last_residue() const {
//   return pos_+length_-1;
//  }

//    virtual std::string to_string() const {
//        std::ostringstream str_stream;
//        str_stream << "(LegacyRootClaim; owner, " << owner()->type() << "; pos, " << pos_ << ")" ;
//
//        //str_stream << "(SequenceClaim; owner, " << owner()->type() << "; pos, " << pos_;
//            //<< "; length, " << length_ << "; label, " << label_ << ", annotated_sequence, "
//            //<< annotated_sequence_ << ")";
//        return str_stream.str();
//    }


}
}
}
#endif
