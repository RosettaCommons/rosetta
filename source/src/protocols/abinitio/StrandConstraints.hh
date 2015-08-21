// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange

#ifndef INCLUDED_protocols_abinitio_StrandConstraints_hh
#define INCLUDED_protocols_abinitio_StrandConstraints_hh

// Unit Headers
//#include <protocols/abinitio/StrandConstraints.fwd.hh>


// Package Headers
#include <protocols/abinitio/PairingStatistics.hh>

//#include <protocols/abinitio/Templates.hh>
//#include <protocols/abinitio/Template.hh>
//#include <protocols/abinitio/TemplateJumpSetup.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>

//#include <core/scoring/constraints/AtomPairConstraint.fwd.hh>
//#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/scoring/dssp/PairingsList.fwd.hh>
//#include <core/fragment/SecondaryStructure.fwd.hh>
#include <core/scoring/dssp/StrandPairing.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/exit.hh>


//// C++ headers
#include <string>
#include <map>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class AlternativePairings {
public:
	// AlternativePairings();
	bool compatible( core::scoring::dssp::StrandPairing const& pairing ) const;
	bool antiparallel() const { return anti_; };
	// bool add_pairing( core::scoring::dssp::StrandPairing const&, std::string model = "NO_MODEL" );
	bool add_pairing( PairingStatEntry const& );
	void show( std::ostream& ) const;

	//@brief make constraints for the alternative registers of this strand
	//pose only for reference,  constraints are added to cst object
	void build_constraints( core::pose::Pose const& pose, core::scoring::constraints::ConstraintCOPs& cst ) const;

private:
	typedef utility::vector1< PairingStatEntry > Pairings;
	Pairings pairings_;
	bool anti_;
};

class StrandConstraints : public utility::pointer::ReferenceCount {
	typedef utility::vector1< AlternativePairings > FuzzyTopology;

public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~StrandConstraints();
	StrandConstraints( PairingStatistics const& strand_stats_ );
	void add_pairing( core::scoring::dssp::StrandPairing const&, std::string model = "NO_MODEL" );
	void add_pairing( PairingStatEntry const& );

	//const pose: no constraints are added
	core::scoring::constraints::ConstraintCOPs build_constraints( core::pose::Pose const& pose ) const;

	void show( std::ostream& ) const;
private:
	FuzzyTopology fuzzy_topology_;
};

std::ostream& operator<< ( std::ostream& out, StrandConstraints const& st );
std::ostream& operator<< ( std::ostream& out, AlternativePairings const& alt_pairs );


}
}

#endif
