// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/FragmentCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abscript_abinitio_FragmentCM_hh
#define INCLUDED_protocols_abscript_abinitio_FragmentCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/FragmentCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>

// Package headers
#include <protocols/simple_moves/FragmentMover.hh>
#include <core/select/residue_selector/ResidueSelector.hh>

// Project headers
#include <basic/datacache/WriteableCacheableMap.fwd.hh>

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class FragmentCM : public protocols::environment::ClientMover {
	typedef protocols::environment::ClientMover Parent;
	typedef environment::claims::EnvClaims EnvClaims;

public:
	FragmentCM();

	FragmentCM( simple_moves::FragmentMoverOP,
		core::select::residue_selector::ResidueSelectorCOP = NULL );

	virtual void set_selector( core::select::residue_selector::ResidueSelectorCOP );

	virtual void set_mover( simple_moves::FragmentMoverOP mover );

	virtual ~FragmentCM();

	void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

	EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) override;

	void initialize( Pose& pose ) override;

	void apply( Pose& pose ) override;

	// XRW TEMP  virtual std::string get_name() const;

	core::select::residue_selector::ResidueSelectorCOP const&
	selector() const { return selector_; }

	bool initialize() const { return bInitialize_; }

	void initialize( bool setting );

	bool yield_cut_bias() const { return bYieldCutBias_; }

	void yield_cut_bias( bool setting );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


protected:
	void passport_updated() override;

	simple_moves::FragmentMoverOP mover() const { return mover_; };

private:
	simple_moves::FragmentMoverOP mover_;
	core::select::residue_selector::ResidueSelectorCOP selector_;
	bool bInitialize_;
	bool bYieldCutBias_;

}; // end FragmentCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_FragmentCM_hh
