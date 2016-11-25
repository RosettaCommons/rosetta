// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     protocols/membrane/VisualizeEmbeddingMoverCreator.hh
/// @brief      Visualize Embedding normal and center with Virtual Residues
/// @details    Add a set of virtual residues as an additional chain to the
///    membrane pose. This tool is strictly for visualization of
///    the implicit membrane and should not be present in modeling.
///    Last Modified: 11/20/14
/// @author  JKLeman (julia.koehler1982@gmail.com)

#ifndef INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_hh
#define INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_hh

// Unit Headers
#include <protocols/membrane/visualize/VisualizeEmbeddingMover.fwd.hh>

// Package headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <protocols/membrane/geometry/Embedding.hh>

#include <core/conformation/Residue.fwd.hh>

// Project Headers
#include <protocols/moves/Mover.hh>

// Utility Headers
#include <numeric/xyzVector.hh>

namespace protocols {
namespace membrane {
namespace visualize {

/// @brief Adds virtual residues to visualize span embeddings of the membrane protein
class VisualizeEmbeddingMover : public protocols::moves::Mover {

public:

	////////////////////
	/// Constructors ///
	////////////////////

	/// @brief   Default constructor
	/// @details  Gets the embeddings from the pose and topology
	VisualizeEmbeddingMover();

	/// @brief   Constructor from embedding
	/// @details  Visualizes defined embedding
	VisualizeEmbeddingMover( protocols::membrane::geometry::EmbeddingOP embedding );

	/// @brief Copy Constructor for deep copying
	VisualizeEmbeddingMover( VisualizeEmbeddingMover const & src );

	/// @brief Assignment Operator
	VisualizeEmbeddingMover &
	operator=( VisualizeEmbeddingMover const & src );

	/// @brief Destructor
	~VisualizeEmbeddingMover();

	///////////////////////////////
	/// Rosetta Scripts Methods ///
	///////////////////////////////

	/// @brief Create a Clone of this mover
	protocols::moves::MoverOP clone() const override;

	/// @brief Create a Fresh Instance of this Mover
	protocols::moves::MoverOP fresh_instance() const override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

	/// @brief   Apply visualization
	void apply( core::pose::Pose & pose ) override;

	/// @brief   Return the name of this mover
	// XRW TEMP  virtual std::string get_name() const;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


private:

	//////////////////////
	/// Helper Methods ///
	//////////////////////

	/// @brief Register Options with JD2
	void register_options();

	/// @brief Create a Membrane Residue
	/// @details Given a centered position and residue typeset, return
	/// a ResidueOP with the xyz coordinate pos, type EMB, from typeset given
	core::conformation::ResidueOP
	create_embedding_virtual( core::Vector center, core::Vector normal, bool fullatom );

private:

	// Embedding object containing (multiple) EmbeddingDefinition(s)
	protocols::membrane::geometry::EmbeddingOP embeddings_;
};

} // visualize
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_visualize_VisualizeEmbeddingMover_hh
