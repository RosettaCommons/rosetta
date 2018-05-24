// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/sewing/hashing/AlignmentGenerator.cc
/// @brief Determines all possible alignments for an edge given the two SewSegments involved. Also contains the mechanisms for reading an edge file.
/// @author guffysl (guffy@email.unc.edu)

#include <protocols/sewing/hashing/AlignmentGenerator.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
static basic::Tracer TR( "protocols.sewing.hashing.AlignmentGenerator" );


namespace protocols {
namespace sewing {
namespace hashing {

AlignmentGenerator::AlignmentGenerator():
	utility::pointer::ReferenceCount()
{
	//We won't initialize from options because the HasherSettings
	//must come from the edge file we read in
}

AlignmentGenerator::AlignmentGenerator(HasherSettings const  & hasher_settings, SegmentVectorCOP segs  ):
	utility::pointer::ReferenceCount(),
	hasher_settings_( hasher_settings ),
	segment_vector_( segs )
{
}


AlignmentGenerator::AlignmentGenerator(HasherSettings const & hasher_settings, SegmentVectorCOP segs, std::map< core::Size, data_storage::SmartSegmentOP > pdbsegs  ):
	utility::pointer::ReferenceCount(),
	hasher_settings_( hasher_settings ),
	segment_vector_( segs ),
	pdb_segments_( pdbsegs )
{
}


AlignmentGenerator::~AlignmentGenerator(){}

AlignmentGenerator::AlignmentGenerator( AlignmentGenerator const & other ):
	hasher_settings_( other.hasher_settings_ ),
	segment_vector_( other.segment_vector_ ),
	pdb_segments_( other.const_pdb_segments() )
{}



AlignmentGeneratorOP
AlignmentGenerator::clone() const {
	return AlignmentGeneratorOP( new AlignmentGenerator( *this ) );
}

//Getters
HasherSettings
AlignmentGenerator::hasher_settings() const {
	return hasher_settings_;
}
SegmentVectorCOP
AlignmentGenerator::segment_vector() const{
	return segment_vector_;
}

std::map< core::Size, data_storage::SmartSegmentOP > &
AlignmentGenerator::pdb_segments(){
	return pdb_segments_;
}

std::map< core::Size, data_storage::SmartSegmentOP >
AlignmentGenerator::const_pdb_segments() const{
	return pdb_segments_;
}
//Setters

void
AlignmentGenerator::hasher_settings( HasherSettings & hasher_settings ){
	hasher_settings_ = hasher_settings;
}
void
AlignmentGenerator::segment_vector( SegmentVectorCOP segvec )
{
	segment_vector_ = segvec;
}
void
AlignmentGenerator::pdb_segments( std::map< core::Size, data_storage::SmartSegmentOP > &  pdbsegs ){
	pdb_segments_ = pdbsegs;
}


utility::vector1< std::pair< core::Size, core::Size > >
AlignmentGenerator::get_all_alignments( data_storage::SmartSegmentCOP s1, data_storage::SmartSegmentCOP s2 ){
	//Create a new Hasher
	HasherOP alignment_hasher = HasherOP( new Hasher( hasher_settings_, segment_vector_, pdb_segments_ ) );

	//Get iterators for possible basis pairs
	std::pair< utility::vector1< std::pair<core::Size, core::Size > >::const_iterator, utility::vector1< std::pair< core::Size, core::Size > >::const_iterator > iters = alignment_hasher->iterate_over_basis_pairs(s1, s2 );
	//Create our vector that we will return
	utility::vector1< std::pair< core::Size, core::Size > > basis_pairs;
	basis_pairs.clear();
	for ( ; iters.first != iters.second; ++iters.first ) {
		data_storage::BasisPair current_basis = std::make_pair( data_storage::Basis(s1->get_segment_id(), iters.first->first), data_storage::Basis( s2->get_segment_id(), iters.first->second ) );
		if ( alignment_hasher->score_basis_pair( current_basis ) ) {
			basis_pairs.push_back( *iters.first );
		}
	}
	return basis_pairs;
}

} //hashing
} //sewing
} //protocols






