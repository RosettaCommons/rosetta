// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Oliver Lange

#ifndef INCLUDED_protocols_abinitio_TemplateJumpSetup_hh
#define INCLUDED_protocols_abinitio_TemplateJumpSetup_hh

// Unit Headers
#include <protocols/abinitio/TemplateJumpSetup.fwd.hh>

// Package Headers
#include <protocols/jumping/JumpSetup.hh>

// Project Headers
#include <core/types.hh>

#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/FragData.fwd.hh>
#include <core/fragment/FrameList.fwd.hh>

#include <core/scoring/dssp/PairingsList.hh>
#include <protocols/jumping/PairingLibrary.hh>
#include <protocols/jumping/JumpSample.hh>

// ObjexxFCL Headers

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

//// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <protocols/abinitio/PairingStatistics.fwd.hh>
#include <protocols/abinitio/Templates.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class TemplateJumpSetup : public jumping::BaseJumpSetup {
public:
	TemplateJumpSetup(
		TemplatesCOP templates,
		core::fragment::SecondaryStructureCOP secstruct,
		PairingStatisticsCOP,
		core::scoring::dssp::PairingList const & helix_pairings
	);

	~TemplateJumpSetup() override;
	std::string type_name() const override {
		return "TemplateJumpSetup";
	}


	jumping::JumpSample create_jump_sample() const override;

	jumping::JumpSample clean_jumps( jumping::JumpSample const& ) const override;

	//' jumping::JumpSample
	// create_jump_sample( std::string ModelID ) const;
	/// @brief returns an ordered FragSet that is compatible with the JumpSample
	/// default: generate jumps from ss-library according to JumpSample

	core::fragment::FragSetOP generate_jump_frags( jumping::JumpSample const&, core::kinematics::MoveMap const& ) const override;

	bool is_helix_jump( core::scoring::dssp::Pairing const& p ) const;

private:
	TemplatesCOP templates_;
	core::fragment::SecondaryStructureCOP secstruct_;

	PairingStatisticsCOP strand_stats_;
	core::scoring::dssp::PairingList helix_pairings_;
	//jumping::BasePairingLibraryOP geometry_library_; ///eventually plug in other sources for geometries .. needs a new interface -- not now
};

class FixTemplateJumpSetup : public TemplateJumpSetup {
public:
	FixTemplateJumpSetup( TemplatesCOP templates, core::fragment::SecondaryStructureCOP secstruct, PairingStatisticsCOP, core::scoring::dssp::PairingList const&,
		jumping::BaseJumpSetupOP jump_def );
	~FixTemplateJumpSetup() override;

	FixTemplateJumpSetup( TemplateJumpSetup const& templ, jumping::BaseJumpSetupOP jump_def );


	jumping::JumpSample create_jump_sample() const override;

private:
	jumping::BaseJumpSetupOP jump_def_;
};


} //abinitio
} //protocols

#endif
