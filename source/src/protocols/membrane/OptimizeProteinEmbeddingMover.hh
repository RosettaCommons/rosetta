// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief      Optimizes the protein embedding in the membrane
/// @details Optimizes the protein embedding in the membrane given the smooth
///   high-res score function; transforms the protein into the membrane,
///   optimizes the membrane position (flexible), and uses the optimized
///   embedding to reposition the protein in the membrane
/// @author     JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_hh
#define INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_hh

// Unit Headers
#include <protocols/membrane/OptimizeProteinEmbeddingMover.fwd.hh>
#include <protocols/membrane/OptimizeProteinEmbeddingMoverCreator.hh>
#include <protocols/moves/Mover.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/filters/Filter.fwd.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <utility/tag/Tag.fwd.hh>
#include <basic/datacache/DataMap.fwd.hh>
#include <basic/Tracer.fwd.hh>

namespace protocols {
namespace membrane {

class OptimizeProteinEmbeddingMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default Constructor
	/// @details Defaults: scorefxn = smooth2012
	OptimizeProteinEmbeddingMover();

	/// @brief Destructor
	~OptimizeProteinEmbeddingMover() override;

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/// @brief Parse Rosetta Scripts Options for this Mover
	void parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap &,
		protocols::filters::Filters_map const &,
		protocols::moves::Movers_map const &,
		core::pose::Pose const &
	) override;

	/// @brief Provide xml schema for RosettaScripts compatibility
	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief Get the name of this Mover (OptimizeProteinEmbeddingMover)
	std::string get_name() const override;

	/// @brief Flip the downstream partner in the membrane
	void apply( core::pose::Pose & pose ) override;

private: // methods

	/////////////////////
	/// Setup Methods ///
	/////////////////////

	/// @brief Register Options from Command Line
	void register_options();

	// Tell the compiler that we are not hiding the base
	// function with the parse_my_tag written above
	using protocols::moves::Mover::parse_my_tag;

};

} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_OptimizeProteinEmbeddingMover_hh
