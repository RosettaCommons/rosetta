// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/jumping/PairingTemplate
/// @brief header file for ClassicAbinitio protocol
/// @details
///  from converting jumping_pairings.cc of rosetta++ into mini
///
///
///
/// @author Oliver Lange

#ifndef INCLUDED_protocols_jumping_SheetBuilder_hh
#define INCLUDED_protocols_jumping_SheetBuilder_hh

// Unit Headers

// Package Headers
#include <protocols/jumping/SameStrand.fwd.hh>
#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/JumpSetup.hh>
#include <core/fragment/SecondaryStructure.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/kinematics/FoldTree.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray3D.fwd.hh>

#include <utility/vector1.hh>
#include <ObjexxFCL/FArray2A.fwd.hh>


//// C++ headers
//#include <cstdlib>
//#include <string>
//#include <vector>

namespace protocols {
namespace jumping {


/// @brief select jumps to build a given topology
/// @detail this class encapsulates the functionality of choose_random_pairings in jumping_pairings.cc of Rosetta++
class SheetBuilder : public BaseJumpSetup {
public:
	typedef utility::vector1< core::Size > SheetTopology;
	SheetBuilder( core::fragment::SecondaryStructureOP, core::scoring::dssp::PairingsList const&, SheetTopology const& );

	//copy c'stor
	SheetBuilder( SheetBuilder const& );

	//d'stor
	virtual ~SheetBuilder();

	std::string type_name() const {
		return "SheetBuilder";
	}

	virtual
	JumpSample
	create_jump_sample() const;

	JumpSample
	clean_jumps( JumpSample const& js ) const
	{
		std::cerr << "ERROR: JumpSetup::clean_jumps() not implemented" << std::endl;
		return js;
	}

	Size total_residue() const {
		return total_residue_;
	}

protected:
	//default do nothing always use input_sheet_sizes_ as sheet_sizes_.
	virtual SheetTopology create_new_random_topol() const
	{ return sheet_sizes_; };

private:
	bool builder_loop( core::scoring::dssp::PairingsList &jump_pairings ) const;
	void choose_next_pairing( ObjexxFCL::FArray3D_int &, core::Size, core::Size ) const;
	bool check_next_pairing(  ObjexxFCL::FArray3D_int &, core::Size, core::Size ) const;

	bool check_two_pairings(
		ObjexxFCL::FArray1A_int pairing1,
		ObjexxFCL::FArray1A_int pairing2,
		int & common_strands
	) const;

	bool check_sheet_pairings(
		ObjexxFCL::FArray2A_int pairing_list,
		const int last_pairing,
		const bool force_single_sheet
	) const;

	bool check_pairing_intersect( ObjexxFCL::FArray1A_int, ObjexxFCL::FArray1A_int ) const;

	core::Size total_residue_;
	core::scoring::dssp::PairingsList pairings_;
	SameStrandOP same_strand_;
	core::fragment::SecondaryStructureOP secondary_structure_;
	mutable SheetTopology sheet_sizes_;

	bool bForceSingleSheet_;
};

} //protocols
} //jumping

#endif

