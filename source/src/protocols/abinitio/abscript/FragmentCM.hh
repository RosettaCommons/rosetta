// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/FragmentCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abscript_abinitio_FragmentCM_hh
#define INCLUDED_protocols_abscript_abinitio_FragmentCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/FragmentCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>

// Package headers
#include <protocols/simple_moves/FragmentMover.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>

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
		core::pack::task::residue_selector::ResidueSelectorCOP = NULL );

	virtual void set_selector( core::pack::task::residue_selector::ResidueSelectorCOP );

	virtual void set_mover( simple_moves::FragmentMoverOP mover );

	virtual ~FragmentCM();

	virtual void
	parse_my_tag( utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	virtual EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP );

	virtual void initialize( Pose& pose );

	virtual void apply( Pose& pose );

	virtual std::string get_name() const;

	core::pack::task::residue_selector::ResidueSelectorCOP const&
	selector() const { return selector_; }

	bool initialize() const { return bInitialize_; }

	void initialize( bool setting );

	bool yield_cut_bias() const { return bYieldCutBias_; }

	void yield_cut_bias( bool setting );

protected:
	virtual void passport_updated();

	simple_moves::FragmentMoverOP mover() const { return mover_; };

private:
	simple_moves::FragmentMoverOP mover_;
	core::pack::task::residue_selector::ResidueSelectorCOP selector_;
	bool bInitialize_;
	bool bYieldCutBias_;

}; // end FragmentCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_FragmentCM_hh
