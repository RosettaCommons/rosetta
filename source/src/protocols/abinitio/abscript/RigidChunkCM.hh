// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/RigidChunkCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh
#define INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/RigidChunkCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>
#include <protocols/environment/claims/EnvClaim.hh>

// Package headers
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.tmpl.hh>

#include <core/select/residue_selector/ResidueSelector.hh>

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

// Project headers

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class RigidChunkCM : public protocols::environment::ClientMover {
	typedef ClientMover Parent;
	typedef environment::claims::EnvClaims EnvClaims;

public:
	RigidChunkCM();

	RigidChunkCM( core::select::residue_selector::ResidueSelectorCOP selector,
		core::pose::Pose const& template_pose );

	virtual ~RigidChunkCM() {};

	EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) override;

	// XRW TEMP  virtual std::string get_name() const;

	void sim_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	core::select::residue_selector::ResidueSelectorCOP sim_selector() const;

	void templ_selector( core::select::residue_selector::ResidueSelectorCOP selector );

	core::select::residue_selector::ResidueSelectorCOP templ_selector() const;

	void initialize( Pose& pose ) override;

	void apply( core::pose::Pose& ) override;

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	moves::MoverOP clone() const override;

	loops::Loops select_parts( loops::Loops const& rigid_core,
		core::Size random_grow_loops_by );

	//  void rigid_core( loops::Loops rigid_core );

	//  loops::Loops const& rigid_core() const;

	core::pose::Pose const& templ() const { return *template_; }

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	std::map< core::Size, core::Size > const& sim_origin() const { return sim_origin_; }

	std::map< core::Size, core::Size > const& templ_target() const { return templ_target_; }

protected:
	void passport_updated() override;

private:
	void configure( core::pose::Pose const& in_p,
		utility::vector1< bool > const selection );

	void templ_target( std::map< core::Size, core::Size > const& in ){ debug_assert( templ_target_.empty() ); templ_target_ = in; }

	void sim_origin( std::map< core::Size, core::Size > const& in ){ debug_assert( sim_origin_.empty() ); sim_origin_ = in; }


	EnvClaims claims_;
	//  loops::Loops rigid_core_;
	core::pose::PoseCOP template_;
	core::select::residue_selector::ResidueSelectorCOP sim_selector_;
	core::select::residue_selector::ResidueSelectorCOP templ_selector_;

	std::string xml_name_;

	// configured during claiming
	std::map< core::Size, core::Size > templ_target_;
	std::map< core::Size, core::Size > sim_origin_;

}; // end RigidChunkCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_RigidChunkCM_hh
