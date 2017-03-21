// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/AbscriptMover.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_AbscriptMover_hh
#define INCLUDED_protocols_abinitio_abscript_AbscriptMover_hh

// Unit Headers
#include <protocols/abinitio/abscript/AbscriptMover.fwd.hh>

// Package headers
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/claims/EnvClaim.hh>

#include <protocols/abinitio/abscript/StageID.hh>
#include <protocols/abinitio/abscript/AbscriptStageMover.hh>
#include <protocols/abinitio/abscript/StagePreparer.fwd.hh>

#include <protocols/moves/MonteCarlo.fwd.hh>
#include <protocols/moves/MoverContainer.fwd.hh>

// Project headers
#include <core/scoring/ScoreFunction.hh>

#include <core/select/residue_selector/ResidueSelector.fwd.hh>

#include <core/pose/Pose.hh>

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

// Utility Headers
#include <utility/vector0.fwd.hh>
#include <utility/tag/Tag.hh>

// C++ Headers
#include <set>
#include <string>

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class AbscriptMover : public protocols::environment::ClientMover {
	typedef environment::claims::EnvClaims EnvClaims;
	typedef environment::ClientMoverOP ClientMoverOP;
	typedef std::set<ClientMoverOP> MoverSet;
	typedef std::map< StageID, MoverSet > IDMoverSetMap;


public:
	AbscriptMover();

	AbscriptMover( AbscriptMover const& );

	void apply( core::pose::Pose& pose ) override;

	static std::string
	stage_complex_namer( std::string );
	// XRW TEMP  virtual std::string get_name() const;

	EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) override;

	void yield_submovers( std::set< ClientMoverOP >& ) const override;

	// the Abscript mover does not make any claims, and should never be given
	// initialization rights
	void initialize( Pose& ) override { utility_exit_with_message("AbscriptMover::initialize() should never be called!"); }

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	moves::MoverOP fresh_instance() const override;

	moves::MoverOP clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	class StageTracker;

	void add_frags( std::string const& small_frags,
		std::string const& large_frags,
		core::select::residue_selector::ResidueSelectorCOP = NULL );

	void parse_subtags( utility::vector0< utility::tag::TagPtr > const&,
		protocols::moves::Movers_map const& );

	void register_submover( protocols::moves::MoverOP, StageIDs const&, core::Real weight );

	void register_preparer( protocols::moves::MoverOP, StageIDs const& );

	std::map< Size, Size > calculate_iterations( core::pose::Pose const& );

	StageIDs parse_stage_id( std::string const& ) const;

	std::map< std::string, StageID > id_map_;
	std::map< StageID, AbscriptStageMoverOP > stage_movers_;
	moves::MonteCarloOP mc_;

	environment::claims::EnvClaims claims_; // list for dynamically made claims to be pushed onto during setup
}; // end AbscriptMover base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_AbscriptMover_hh
