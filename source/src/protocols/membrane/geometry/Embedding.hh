// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/geometry/Embedding.hh
///
/// @brief      Methods for Computing Membrane Embeddings
/// @details    Includes methods for computing membrane embeddings from search,
///				sequence, structure, user input, and smaller component embeddings
///				Last Modified: 7/24/14
///
/// @author		Julia Koehler (julia.koehler1982@gmail.com)
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_geometry_Embedding_hh
#define INCLUDED_protocols_membrane_geometry_Embedding_hh

// Unit Headers
#include <protocols/membrane/geometry/Embedding.fwd.hh>

// Project Headers
#include <protocols/membrane/geometry/EmbeddingDef.fwd.hh>
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
#include <iostream>

namespace protocols {
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

	/// @brief	Construction from single EmbeddingDef object
	Embedding( EmbeddingDefOP embedding );

	/// @brief	Constructs bogus object from topology
	Embedding( SpanningTopologyOP topology, Real radius );

	/// @brief Custom Constructor - from topology & structure
	/// @details Construct Embedding from Structure & Topology
	Embedding( SpanningTopologyOP topology, Pose & pose );

	/// @brief Custom Constructor - from topology & structure
	/// @details Construct Embedding from Structure & Topology
	Embedding( SpanningTopologyOP topology, PoseOP pose );

	/// @brief Copy Constructor
	Embedding( Embedding const & Embedding );

	/// @brief Assignment Operator
	Embedding & operator = ( Embedding const & src );

	/// @brief	Destructor
	~Embedding();
    
public: // methods

	// show object
	virtual void show( std::ostream & out=std::cout );
	
	// invert all normals in Embedding object
	void invert();
	
	// number of span embeddings in object
	Size nspans();

	// get span embedding by number
	EmbeddingDefOP embedding( Size span_number );

	// get all span embeddings
	utility::vector1< EmbeddingDefOP > embeddings();

	// get chain embedding
	EmbeddingDefOP total_embed();

	//////////////////////
	/// HELPER METHODS ///
	//////////////////////

	// from TMspans
	utility::vector1< EmbeddingDefOP > from_spans( SpanningTopologyOP topology, PoseOP pose );

	// from TMspans
	utility::vector1< EmbeddingDefOP > from_spans( SpanningTopologyOP topology, Pose & pose );

private: // data

	// embedding per span
	utility::vector1< EmbeddingDefOP > embeddings_;

	// average embedding of all span embeddings
	EmbeddingDefOP total_embed_;

}; // Embedding

} // geometry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_geometry_Embedding_hh
