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

// Unit Headers
#include <protocols/membrane/geometry/Embedding.hh>
#include <protocols/membrane/geometry/util.hh>

// Project Headers
#include <protocols/membrane/geometry/EmbeddingDef.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <numeric/conversions.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <string>
#include <cstdlib>

static thread_local basic::Tracer TR( "protocols.membrane.geometry.Embedding" );

namespace protocols {
namespace membrane {
namespace geometry {

using namespace core; 
using namespace core::conformation::membrane;

////////////////////
/// Constructors ///
////////////////////

/// @brief	Constructs empty object
Embedding::Embedding() :
embeddings_(),
total_embed_(){}


/// @brief	Constructs bogus object from topology
Embedding::Embedding( SpanningTopologyOP topology, Real radius ) :
	embeddings_(),
	total_embed_()
{
	// initialize embedding objects

    TR << "Create a bogus embedding from topology" << std::endl;
	using namespace numeric::conversions;
		
	// initializations
	core::Real increment = 360 / topology->nspans();

	// get position around a circle
	for ( Size i = 1; i <= topology->nspans(); ++i ) {

		core::Real alpha = radians( increment * i );
		core::Real x = radius * cos( alpha );
		core::Real y = radius * sin( alpha );
		
		// For each TMspan, get center
		core::Vector center( x, y, 0 );
		core::Vector normal;
		
		// get normal for each span, respects protein topology
		if ( i % 2 == 1 ){
			normal.assign(0, 0, 1); // odd spans in positive z
		}
		else {
			normal.assign(0, 0, -1); // even spans in negative z
		}

		// create object and push back into vector
		EmbeddingDefOP span_embedding( new EmbeddingDef( center, normal ) );
		embeddings_.push_back( span_embedding );

	}

	// chain embedding is sum of span embeddings
	total_embed_ = average_embeddings( embeddings_ );
}// object from topology


/// @brief	Constructs from topology and pose
Embedding::Embedding( SpanningTopologyOP topology, PoseOP pose ){

    // create span embedding from topology
    embeddings_ = from_spans( topology, pose );
    
    // calculate chain embedding, push back
    total_embed_ = average_embeddings( embeddings_ );
}// from topology and pose

/// @brief	Destructor
Embedding::~Embedding(){}

///////////////
/// Methods ///
///////////////

// show object
void Embedding::show( std::ostream & out ){
    
    out << "Span Embedding: " << std::endl;
    for ( Size i = 1; i <= embeddings_.size(); ++i ){
        out << "\t";
        embeddings_[ i ]->show();
    }
    out << "Chain Embedding: " << std::endl;
	out << "\t";
    total_embed_->show();

}// show

// number of span embeddings in object
Size Embedding::nspans() {
	return embeddings_.size();
}

//////////////////////////////////////////////////////////////////////////////

// get span embedding by number
EmbeddingDefOP Embedding::embedding( Size span_number ){
    return embeddings_[ span_number ];
}

//////////////////////////////////////////////////////////////////////////////

// get all span embeddings
utility::vector1< EmbeddingDefOP > Embedding::embeddings(){
    return embeddings_;
}// span_embeddings

//////////////////////////////////////////////////////////////////////////////

// get chain embedding
EmbeddingDefOP Embedding::total_embed(){
    total_embed_ = average_embeddings( embeddings_ );
    return total_embed_;
}// chain_embedding


// from TMspans
utility::vector1< EmbeddingDefOP > Embedding::from_spans( SpanningTopologyOP topology, PoseOP pose ){

    TR << "Computing membrane embedding from TMspans: " << std::endl;
    topology->show();

    // go through TMspans
    for ( Size i = 1; i <= topology->nspans(); ++i ){

        // residues
        Size start( topology->span( i )->start() );
        Size end( topology->span( i )->end() );
		
		TR << "start: " << start << std::endl;
		TR << "end: " << end << std::endl;
        
        // create embedding object and fill it from span
        EmbeddingDefOP embedding( new EmbeddingDef( pose, start, end ) );
        
        // add embedding to vector
        embeddings_.push_back( embedding );
    }
    return embeddings_;
}// from spans


} // geometry
} // membrane
} // protocols
