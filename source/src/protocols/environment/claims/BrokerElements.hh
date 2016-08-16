// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BrokerElements
/// @brief
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_claims_BrokerElements_hh
#define INCLUDED_protocols_environment_claims_BrokerElements_hh

// Unit Headers
// #include <protocols/environment/claims/BrokerElements.fwd.hh>

// Package Headers
#include <protocols/environment/ClientMover.fwd.hh>

// Project Headers
#include <core/environment/LocalPosition.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/TorsionID.hh>

// ObjexxFCL Headers

// Utility headers

//// C++ headers
#include <string>

namespace protocols {
namespace environment {
namespace claims {

// ControlStrengths indicate the level of control over the dof desired by the mover.
// DOES_NOT_CONTROL is the lowest level, used for protocols that exert no control
//   over the claimed dof (perhaps, for example, because they only require that it exists).
// CAN_CONTROL is also a low, uncommon level indicating the mover would like to control
//   (i.e. samples) the dof, but needn't. Fragment Insertion is often reliant on claims
//   of this type, because other movers can take priority (i.e. if a region is fixed).
// MUST_CONTROL asserts a need to sample the degree of freedom (e.g. jumps in docking).
//   This level allows the mover, however, to share the dof with other movers.
// EXCLUSIVE is as MUST_CONTROL, except that no other mover can be granted access to the dof.
enum ControlStrength  {
	DOES_NOT_CONTROL = 0,
	CAN_CONTROL = 1,
	MUST_CONTROL,
	EXCLUSIVE
};

extern std::ostream& operator<<( std::ostream& os, ControlStrength const& cstr );

/// @note at the moment, this only makes virtual residues, but could be modified to hold on to the name of
///       the residue that should be created so that you could (for example) add a carbohydrate or ligand with
///       this system.
struct ResidueElement {
	ResidueElement() : label(), allow_duplicates( false ) {}
	std::string label;
	bool allow_duplicates; //if this label already exists, simply ignore this Element rather than throwing.
	static std::string const type;
};

struct JumpElement {
	JumpElement() : label(), p1(), p2(), atom1(), atom2(),
		force_stub_intra_residue( false ), has_physical_cut(false) {}
	std::string label;
	core::environment::LocalPosition p1;
	core::environment::LocalPosition p2;
	std::string atom1;
	std::string atom2;
	bool force_stub_intra_residue; //see FoldTree::put_jump_stubs_intra_residue
	bool has_physical_cut;
	static std::string const type;
};

struct CutElement {
	CutElement() : p() {}
	core::environment::LocalPosition p;
	static std::string const type;
};

struct CutBiasElement {
	CutBiasElement() : p(), bias( 0.0 ) {}
	core::environment::LocalPosition p;
	core::Real bias;
	static std::string const type;
};

struct DOFElement {
	DOFElement() : id(), c_str( DOES_NOT_CONTROL ), i_str( DOES_NOT_CONTROL ) {}
	core::id::DOF_ID id;
	ControlStrength c_str;
	ControlStrength i_str;
	static std::string const type;
};

typedef utility::vector1< claims::ResidueElement > ResidueElements;
typedef utility::vector1< claims::JumpElement > JumpElements;
typedef utility::vector1< claims::CutElement > CutElements;
typedef utility::vector1< claims::CutBiasElement > CutBiasElements;
typedef utility::vector1< claims::DOFElement > DOFElements;

typedef utility::vector1< claims::ControlStrength > ControlStrengths;
}
}
}

#endif
