// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/abinitio/abscript/StructPerturberCM.hh
/// @author Justin Porter

#ifndef INCLUDED_protocols_abinitio_abscript_StructPerturberCM_hh
#define INCLUDED_protocols_abinitio_abscript_StructPerturberCM_hh

// Unit Headers
#include <protocols/abinitio/abscript/StructPerturberCM.fwd.hh>
#include <protocols/environment/ClientMover.hh>

#ifdef WIN32
#include <basic/datacache/WriteableCacheableMap.hh>
#endif

// Package headers

// Project headers

// C++ Headers

// ObjexxFCL Headers

namespace protocols {
namespace abinitio {
namespace abscript {

class StructPerturberCM : public protocols::environment::ClientMover {
	typedef ClientMover Parent;
	typedef environment::claims::EnvClaims EnvClaims;

public:

	StructPerturberCM();

	StructPerturberCM( std::string const& label,
		core::Real magnitude );

	virtual ~StructPerturberCM() {};

	EnvClaims yield_claims( core::pose::Pose const&,
		basic::datacache::WriteableCacheableMapOP ) override;


	void apply( core::pose::Pose& ) override;

	std::string const& label() const { return label_; }

	void label( std::string const& label ) { label_ = label; }

	core::Real const& magnitude() const { return magnitude_; }

	void magnitude( core::Real const& value ) { magnitude_ = value; }

	void
	parse_my_tag(utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;

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
	core::Real magnitude_;
	std::string label_;

}; // end StructPerturberCM base class

} // abscript
} // abinitio
} // protocols

#endif //INCLUDED_protocols_abinitio_abscript_StructPerturberCM_hh
