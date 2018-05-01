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


#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymMinMover_hh

// Unit headers
#include <protocols/minimization_packing/symmetry/SymMinMover.fwd.hh>
#include <protocols/minimization_packing/MinMover.hh>

#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>

#include <utility/vector1.hh>


// Package headers

// ObjexxFCL Headers

// C++ Headers

// Utility Headers

namespace protocols {
namespace minimization_packing {
namespace symmetry {

///////////////////////////////////////////////////////////////////////////////
/// @brief This class is kept only for backward compatibility --
/// The regular MinMover should now handle symmetric poses transparently and is preferred.
class SymMinMover : public protocols::minimization_packing::MinMover
{
public:

	// Constructors should be inherited from MinMover
	using MinMover::MinMover;

	~SymMinMover();

	// These are kept to keep the correct typing
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

};

} // symmetry
} // moves
} // rosetta
#endif
