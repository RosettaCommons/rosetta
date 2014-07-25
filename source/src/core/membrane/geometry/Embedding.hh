// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		core/membrane/geometry/Embedding.hh
///
/// @brief      Methods for Computing Membrane Embeddings
/// @details    Includes methods for computing membrane embeddings from search,
///				sequence, structure, user input, and smaller component embeddings
///				Last Modified: 7/24/14
///
/// @author		Julia Koehler (julia.koehler1982@gmail.com)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_core_membrane_geometry_Embedding_hh
#define INCLUDED_core_membrane_geometry_Embedding_hh

// Unit Headers
#include <core/membrane/geometry/Embedding.fwd.hh>

// Project Headers
#include <core/membrane/geometry/EmbeddingDef.fwd.hh>
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>

// Package Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <cstdlib>

namespace core {
namespace membrane {
namespace geometry {

using namespace core;
using namespace core::pose;
using namespace core::conformation::membrane;

class Embedding : public utility::pointer::ReferenceCount {
    
public: // constructors

	/// @brief	Detault Constructors
	/// @details Construct an empty embedding object
	Embedding();

	/// @brief Custom Cosntructor - From Topology
	/// @details Constructs a membrane embedding object from spanning topology
	/// This is kind of a bogus topology since we don't know about the structure
	Embedding( SpanningTopologyOP topology );

	/// @brief Custom Constructor - from topology & structure
	/// @details Construct Embedding from Structure & Topology
	Embedding( SpanningTopologyOP topology, PoseOP pose );

	/// @brief	Constructs from user-defined embedding
	/// @details COnstruct embedding from topology, pose and user defined embedding
	Embedding( SpanningTopologyOP topology, PoseOP pose, EmbeddingDefOP user_embedding );

	/// @brief	Destructor
	~Embedding();
    
public: // methods

	/// @brief Show this embedding
	/// @details Show the embedding parameters - center & normal
	void show();

	/// @brief Get Embedding associated with a particular transmembrane span
	EmbeddingDefOP span_embedding( Size span_number );

	/// @brief Get a vector of embeddings per transmembrane span
	utility::vector1< EmbeddingDefOP > span_embeddings();

	/// @brief Get the embedding definition for chains
	EmbeddingDefOP chain_embedding();

	/// @brief Compute the total embedding from multiple embeddings
	EmbeddingDefOP sum_of_parts( utility::vector1< EmbeddingDefOP > parts );

	/// @brief Compute Embedding from transmembrane spans
	utility::vector1< EmbeddingDefOP > from_spans( SpanningTopologyOP topology, PoseOP pose );

	/// @brief Compute emmebrane embedding from a single transmembrane span 
	EmbeddingDefOP from_span( Size start, Size end, PoseOP pose );

private: // data

	// embedding per span
	utility::vector1< EmbeddingDefOP > span_embeddings_;

	// embedding per chain
	EmbeddingDefOP chain_embedding_;

}; // Embedding

} // geometry
} // membrane
} // core

#endif // INCLUDED_core_membrane_geometry_Embedding_hh
