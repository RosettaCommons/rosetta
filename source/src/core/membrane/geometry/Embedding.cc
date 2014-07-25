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

// Unit Headers
#include <core/membrane/geometry/Embedding.hh>

// Project Headers
#include <core/membrane/geometry/EmbeddingDef.hh>
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/types.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>

#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <string>
#include <cstdlib>

static basic::Tracer TR( "core.membrane.geometry.Embedding" );

namespace core {
namespace membrane {
namespace geometry {

using namespace core::conformation::membrane;

////////////////////
/// Constructors ///
////////////////////

/// @brief	Constructs empty object
Embedding::Embedding() :
	span_embeddings_(),
	chain_embedding_()
{}


/// @brief	Constructs bogus object from topology
Embedding::Embedding( SpanningTopologyOP topology ){

    TR << "Create a bogus embedding from topology" << std::endl;

    // initialize embedding vectors
    utility::vector1< EmbeddingDefOP > span_embeddings;

    // initializations
    core::Vector chain_com(0, 0, 0);
    core::Real radius(8);
    core::Real increment = 360 / topology->total_spans();

    // get position around a circle
    for ( Size i = 1; i <= topology->total_spans(); ++i ) {
        core::Real alpha = increment * i;
        core::Real x = radius * cos( alpha );
        core::Real y = radius * sin( alpha );

        // For each TMspan, get center and normal
        core::Vector center( x, y, 0 );
        core::Vector normal(0, 0, 1);

        // create object and push back into vector
        EmbeddingDefOP span_embedding = new EmbeddingDef( normal, center );
        span_embeddings_.push_back( span_embedding );
    }

    // chain embedding is sum of span embeddings
    chain_embedding_ = sum_of_parts( span_embeddings_ );
}// object from topology


/// @brief	Constructs from topology and pose
Embedding::Embedding( SpanningTopologyOP topology, PoseOP pose ){

    // create span embedding from topology
    span_embeddings_ = from_spans( topology, pose );

    // calculate chain embedding, push back
    chain_embedding_ = sum_of_parts( span_embeddings_ );
}// from topology and pose


/// @brief	Constructs from user-defined embedding
Embedding::Embedding( SpanningTopologyOP topology, PoseOP pose, EmbeddingDefOP user_embedding ){

    // create span embeddings from topology and pose
    span_embeddings_ = from_spans( topology, pose );

    // calculate chain embedding
    EmbeddingDefOP computed_embedding( sum_of_parts( span_embeddings_ ) );

    // get translation vector from (user - computed)
    core::Vector const translation_center( user_embedding->center() - computed_embedding->center() );
    core::Vector const translation_normal( user_embedding->normal() - computed_embedding->normal() );

    // create embedding object
    EmbeddingDefOP translation = new EmbeddingDef( translation_normal, translation_center );

    // translate span embeddings
    for ( Size i = 1; i <= span_embeddings_.size(); ++i ){
        span_embeddings_[ i ]->translate_by( translation );
    }

    // set chain embedding
    chain_embedding_ = user_embedding;

}// from user-defined embedding

/// @brief	Destructor
Embedding::~Embedding(){}

///////////////
/// Methods ///
///////////////

/// @brief Show this embedding
/// @details Show the embedding parameters - center & normal
void Embedding::show(){
    TR << "Chain Embedding: " << std::endl;
    chain_embedding_->show();

    TR << "Span Embedding: " << std::endl;
    for ( Size i = 1; i <= span_embeddings_.size(); ++i ){
        TR << "			";
        span_embeddings_[ i ]->show();
    }
}// show


/// @brief Get Embedding associated with a particular transmembrane span
EmbeddingDefOP Embedding::span_embedding( Size span_number ){
    return span_embeddings_[ span_number ];
}

/// @brief Get a vector of embeddings per transmembrane span
utility::vector1< EmbeddingDefOP > Embedding::span_embeddings(){
    return span_embeddings_;
}// span_embeddings

/// @brief Get the embedding definition for chains
EmbeddingDefOP Embedding::chain_embedding(){
    chain_embedding_ = sum_of_parts( span_embeddings_ );
    return chain_embedding_;
}// chain_embedding


/// @brief Compute the total embedding from multiple embeddings
EmbeddingDefOP Embedding::sum_of_parts( utility::vector1< EmbeddingDefOP > parts ){

    // Initialize vars
    core::Vector center(0, 0, 0);
    core::Vector normal(0, 0, 0);

    // Compute Resulting Center
    for ( Size i = 1 ; i <= parts.size(); ++i ) {
        center += parts[i]->center();
        normal += parts[i]->normal();
    }
    center /= parts.size();
    normal /= parts.size();

    // Create new embedding setup and return it
    EmbeddingDefOP embedding = new EmbeddingDef( normal, center );
    return embedding;
}// sum of parts


/// @brief Compute Embedding from transmembrane spans
utility::vector1< EmbeddingDefOP > Embedding::from_spans( SpanningTopologyOP topology, PoseOP pose ){

    TR << "Computing membrane embedding from TMspans: " << std::endl;
    topology->show();

    // go through TMspans
    for ( Size i = 1; i <= topology->total_spans(); ++i ){

        // residues
        Size start( topology->span( i )->start() );
        Size end( topology->span( i )->end() );

        // create embedding object and fill it from span
        EmbeddingDefOP embedding = new EmbeddingDef( *from_span( start, end, pose ) );

        // add embedding to vector
        span_embeddings_.push_back( embedding );
    }
    return span_embeddings_;
}// from spans



/// @brief Compute emmebrane embedding from a single transmembrane span
EmbeddingDefOP Embedding::from_span( Size start, Size end, PoseOP pose ){

    TR << "Computing membrane embedding from TMspan " << start << " to " << end << std::endl;

    // get CA atom positions of anchor residues
    core::Vector pos1 = pose->residue( start ).atom( 2 ).xyz();
    core::Vector pos2 = pose->residue( end ).atom( 2 ).xyz();

    // check TMspan
    if ( pos1.z() < 0 && pos2.z() < 0 ){
        throw new EXCN_Illegal_Arguments("Your TMspan does not span the membrane!");
    }
    else if ( pos1.z() > 0 && pos2.z() > 0 ){
        throw new EXCN_Illegal_Arguments("Your TMspan does not span the membrane!");
    }

    // compute center
    core::Vector center = 0.5 * ( pos1 + pos2 );

    // compute normal in direction of increasing z-axis
    core::Vector normal;
    if ( pos1.z() > 0 && pos2.z() < 0 ){
        normal = pos2 - pos1;
    }
    else {
        normal = pos1 - pos2;
    }
    normal.normalize();

    // Create new embedding and return it
    EmbeddingDefOP embedding = new EmbeddingDef( normal, center );
    return embedding;
}// from single span

} // geometry
} // membrane
} // core
