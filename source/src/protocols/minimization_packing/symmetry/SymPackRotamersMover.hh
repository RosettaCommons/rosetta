// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre
/// @author Gutted by Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_hh

// Unit headers
#include <protocols/minimization_packing/symmetry/SymPackRotamersMover.fwd.hh>
#include <protocols/minimization_packing/PackRotamersMover.hh>

// Project headers

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

namespace protocols {
namespace minimization_packing {
namespace symmetry {

/// @brief This class is provided for backwards compatibility only
/// The PackRotamersMover should work with symmetry transparently and is preferred.
class SymPackRotamersMover : public protocols::minimization_packing::PackRotamersMover {

public:

	// Constructors same as PackRotamersMover
	using PackRotamersMover::PackRotamersMover;

	~SymPackRotamersMover();

	// Needed to preserve type information
	protocols::moves::MoverOP clone() const override;
	protocols::moves::MoverOP fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:

};

} // symmetry
} // moves
} // protocols

#endif
