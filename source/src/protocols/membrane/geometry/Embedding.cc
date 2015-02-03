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
#include <core/conformation/membrane/types.hh>
#include <basic/Tracer.hh>

#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

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

/// @brief	Construction from single EmbeddingDef object
Embedding::Embedding( EmbeddingDef const & embedding ) {
	EmbeddingDefOP copy( new EmbeddingDef( embedding ) );
	embeddings_.push_back( copy );
	total_embed_ = copy;
}

////////////////////////////////////////////////////////////////////////////////

/// @brief	Constructs bogus object from topology
Embedding::Embedding( SpanningTopology const & topology, Real radius ) :
	embeddings_()
{
	// initialize embedding objects
	EmbeddingDefOP new_embed( new EmbeddingDef() );
	total_embed_ = new_embed;
	
    TR << "Create a bogus embedding from topology" << std::endl;
	using namespace numeric::conversions;
		
	// initializations
	core::Real increment = 360 / topology.nspans();

	// get position around a circle
	for ( Size i = 1; i <= topology.nspans(); ++i ) {

		core::Real alpha = radians( increment * i );
		core::Real x = radius * cos( alpha );
		core::Real y = radius * sin( alpha );
		
		// For each TMspan, get center
		core::Vector center( x, y, 0 );
		core::Vector normal;
		
		// get normal for each span, respects protein topology
		if ( i % 2 == 1 ){
			normal = mem_normal; // odd spans in positive z
		}
		else {
			normal = mem_normal.negated(); // even spans in negative z
		}

		// create object and push back into vector
		EmbeddingDefOP span_embedding( new EmbeddingDef( center, normal ) );
		embeddings_.push_back( span_embedding );
	}

	// chain embedding is sum of span embeddings
	total_embed_ = average_antiparallel_embeddings( embeddings_ );
    
}// object from topology

////////////////////////////////////////////////////////////////////////////////

/// @brief	Constructs from topology and pose
Embedding::Embedding( SpanningTopology const & topology, Pose const & pose ){

	TR << "Constructing Embedding object from topology and pose" << std::endl;
	
    // create span embedding from topology
    embeddings_ = from_spans( topology, pose );
	
    // calculate chain embedding, push back
    total_embed_ = average_embeddings( embeddings_ );
	
	// make sure the normal of the total embedding shows into positive z direction
	if ( total_embed_->normal().z() < 0 ) {
		this->invert();
	}
	
	this->show();
	
}// from topology and pose

////////////////////////////////////////////////////////////////////////////////

/// @brief Copy Constructor
Embedding::Embedding( Embedding const & src ) :
	utility::pointer::ReferenceCount(),
	embeddings_( src.embeddings_ ),
	total_embed_( src.total_embed_ )
{}

/// @brief Assignment Operator
Embedding & Embedding::operator = ( Embedding const & src ) {

	// Abort self-assignment.
	if ( this == &src ) { return *this; }

	// Deep Copy of the data
	this->embeddings_ = src.embeddings_;
	this->total_embed_ = src.total_embed_;

	return *this;
}

/// @brief	Destructor
Embedding::~Embedding(){}

////////////////////////////////////////////////////////////////////////////////

///////////////
/// Methods ///
///////////////

// show object
void Embedding::show( std::ostream & out ) const {
    
    out << "Span Embedding: " << std::endl;
    for ( Size i = 1; i <= embeddings_.size(); ++i ){
        out << "\t";
        embeddings_[ i ]->show();
    }
    out << "Total Embedding: " << std::endl;
	out << "\t";
    total_embed_->show();

}// show

// invert all normals in Embedding object
void Embedding::invert() {
	
	// go through vector of EmbeddingDefs and invert the normals
	for ( Size i = 1; i <= embeddings_.size(); ++i ) {
		if ( embeddings_[i]->normal().length() > 0 ) {
			embeddings_[i]->invert();
		}
	}

	// invert total embedding
	if ( total_embed_->normal().length() > 0 ) {
		total_embed_->invert();
	}
}

// number of span embeddings in object
Size Embedding::nspans() const {
	return embeddings_.size();
}

//////////////////////////////////////////////////////////////////////////////

// get span embedding by number
EmbeddingDefOP Embedding::embedding( Size span_number ) const {

	if ( span_number <= embeddings_.size() ) {
		return embeddings_[ span_number ];
	}
	else {
		utility_exit_with_message( "embedding of span_number doesn't exist!" );
	}
}

//////////////////////////////////////////////////////////////////////////////

// get all span embeddings
utility::vector1< EmbeddingDefOP > Embedding::embeddings() const {
    return embeddings_;
}// span_embeddings

//////////////////////////////////////////////////////////////////////////////

// get chain embedding
EmbeddingDefOP Embedding::total_embed() const {
    return total_embed_;
}// chain_embedding

// from TMspans, they are all in a similar direction as the first TM span!
utility::vector1< EmbeddingDefOP > Embedding::from_spans( SpanningTopology const & topology, Pose const & pose ){

    TR << "Computing membrane embedding from TMspans: " << std::endl;

	if ( topology.nspans() == 0 ){
		utility_exit_with_message("The topology object is empty!");
	}

	// get first embedding for comparison of inside/out orientation
	Size start1( topology.span( 1 )->start() );
	Size end1( topology.span( 1 )->end() );
	
	// create embedding object and fill it from span for first embedding object in positive z direction
	EmbeddingDef embedding1 = EmbeddingDef( pose, start1, end1 );
	Vector const center1 = embedding1.center();
	Vector const normal1 = embedding1.normal();
	TR << "first span center and normal: " << center1.to_string() << ", " << normal1.to_string() << std::endl;
	
    // go through rest of TMspans
    for ( Size i = 1; i <= topology.nspans(); ++i ){

        // residues
        Size start( topology.span( i )->start() );
        Size end( topology.span( i )->end() );
		
        // create embedding object and fill it from span
        EmbeddingDefOP embedding( new EmbeddingDef( pose, start, end ) );
		Vector center = embedding->center();
		Vector normal = embedding->normal();
		TR << "center: " << center.to_string() << ", normal: " << normal.to_string() << std::endl;

		// calculate points for dihedral angle calculation
		Vector p1 = center1 + normal1;
		Vector p  = embedding->center() + embedding->normal();
		
		// calculate dihedral angle between normals of first object and this one
		Real dih( numeric::dihedral_degrees( p1, center1, center, p ) );
		TR << "dihedral angle to first span: " << dih << std::endl;

		// check if angle of normal is < 100 degrees to first normal
		// if yes, then add to embedding object, if no add inverted vector
		if ( dih > -100 && dih < 100 ) {
			// add embedding to overall embedding
			embeddings_.push_back( embedding );
		}
		else {
			// add inverted embedding to embedding object
			Vector new_normal( normal.negated() );
			EmbeddingDefOP new_embed( new EmbeddingDef( center, new_normal ) );
			embeddings_.push_back( new_embed );
		}

		// TODO: more highres: don't just take a single CA atom for start and end
		// but the COM of the residues i, i-1, i+1

	}
    return embeddings_;
}// from spans

} // geometry
} // membrane
} // protocols
