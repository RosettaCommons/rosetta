// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file ./src/protocols/fldsgn/potentials/sspot/NatbiasSheetPotential.cc
/// @brief header file of NatbiasSheetPotential
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#ifndef INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasSheetPotential_hh
#define INCLUDED_protocols_fldsgn_potentials_sspot_NatbiasSheetPotential_hh

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <protocols/fldsgn/topology/SS_Info2.fwd.hh>
#include <protocols/fldsgn/potentials/sspot/NatbiasSheetPotential.fwd.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace fldsgn {
namespace potentials {
namespace sspot {

class NatbiasSheetPotential : public utility::pointer::ReferenceCount {
public:


	typedef core::Size Size;
	typedef core::Real Real;
	typedef core::Vector Vector;
	typedef core::pose::Pose Pose;
	typedef protocols::fldsgn::topology::SS_Info2 SS_Info2;
	typedef protocols::fldsgn::topology::SS_Info2_OP SS_Info2_OP;
	typedef protocols::fldsgn::topology::SS_Info2_COP SS_Info2_COP;


public: // construct/destruct


	/// @brief default constructor
	NatbiasSheetPotential();

	/// @brief value constructor
	NatbiasSheetPotential( SheetingSetOP const hpairset );

	/// @brief copy constructor
	NatbiasSheetPotential( NatbiasSheetPotential const & src );

	/// @brief default destructor
	~NatbiasSheetPotential();


public: // mutator


	void set_hpairset( SheetingSetOP const hpairset );


public: // useful functions


	void show( Pose const & pose ) const;


public: // scoring


	/// @brief score secondary structure
	void score( SS_Info2_COP const ssinfo ) const;


public: // mutator





private: // secondary structure data




};

} // ns sspot
} // ns potentials
} // ns fldsgn
} // ns protocols


#endif
