// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///	@file		core/conformation/membrane/SpanningTopology.cc
///
/// @brief      Transmembrane Spans Object
/// @details	Object for storing membrane spanning regions as a vector1 of span objects.
///				Spanning regions can be determined either from an input spanfile, xyz coordinates
///				in the pose, or sequence. Object is constructed from span regions and will do internal
///				checking for validity.
///				Last Modified: 7/20/14
///
/// @author		Julia Koehler (julia.koehler1982@gmail.com)
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

// Project Headers
#include <core/conformation/membrane/Exceptions.hh>

// Package Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>

static thread_local basic::Tracer TR( "core.conformation.membrane.SpanningTopology" );

namespace core {
namespace conformation {
namespace membrane {

using namespace core;

////////////////////
/// Constructors ///
////////////////////

/// @brief	Default Constructor (Private)
/// @details Construct an Empty Spanning Topology Object
SpanningTopology::SpanningTopology() :
	utility::pointer::ReferenceCount(),
	topology_()
{}

/// @brief	Custom Constructor - Transmembrane Spans from Spanfile
/// @details Use transmembrane spans provided to consturct a spanning topology object
SpanningTopology::SpanningTopology(
	std::string spanfile,
	Size total_residues
	) : utility::pointer::ReferenceCount(),
	topology_()
{
    create_from_spanfile( spanfile, total_residues );
	   
} // topology from spanfile

/// @brief	Custom Constructor - Transmembrane Spans from xyz coords
/// @details Use coordinates of residue CA and thickness to determine the spanning regions in the pose
SpanningTopology::SpanningTopology(
	utility::vector1< Real > res_z_coord,
	utility::vector1< Size > chainID, Real thickness
	) : utility::pointer::ReferenceCount(),
	topology_()
{
	create_from_structure( res_z_coord, chainID, thickness);

} // topology from pose and thickness

/// @brief	Copy Constructor
/// @details Create a deep copy of this object copying over all private fields
SpanningTopology::SpanningTopology( SpanningTopology const & src ) :
	utility::pointer::ReferenceCount( src ),
	topology_( src.topology_ )
{}

/// @brief Assignment Operator
/// @details Overload assignemnt operator - required Rosetta method
SpanningTopology &
SpanningTopology::operator=( SpanningTopology const & src ) {
	
	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}
	
	// Otherwise, create a new object
	return *( new SpanningTopology( *this ) );
	
}

/// @brief Destructor
SpanningTopology::~SpanningTopology(){}

///////////////
///	Methods ///
///////////////

/// @brief Show the current spans stored in this SpanningTopology Object
/// @details Generating a String Representation of Spanning Topology Object for debugging purposes
void SpanningTopology::show( std::ostream & output ) const {
    
	// Print Total number of transmembrane Spans
    output << "Total # of TM spans: " << topology_.size() << std::endl;
    
	// print individual spans
    for ( Size i = 1; i <= topology_.size(); ++i ){
        TR << "Span " << i << ": start: " << topology_[ i ]->start();
        TR << ", end: " << topology_[ i ]->end() << std::endl;
    }

} // show

// write spanfile
void SpanningTopology::write_spanfile( std::string output_filename ){

	TR.Debug << "printing spanfile" << std::endl;
	
	// print header
	utility::io::ozstream OUT;
	OUT.open( output_filename );
	OUT << "Rosetta-generated spanfile from SpanningTopology object" << std::endl;
	OUT << topology_.size() << " " << nres_topo_ << std::endl;
	OUT << "antiparallel" << std::endl;
	OUT << "n2c" << std::endl;
	
	// print spans
	for ( Size i = 1; i <= topology_.size(); ++i ){
		OUT << "\t" << topology_[i]->start() << "\t" << topology_[i]->end() << std::endl;
	}
	OUT.close();
	TR << "wrote " << output_filename << std::endl;
	
}// write spanfile


/// @brief Return Spanning Topology
/// @details return spanning topology as a vector1 of transmembrane spans
utility::vector1< SpanOP >
SpanningTopology::get_spans() { return topology_; } // get topology

// add span to end of SpanningTopology object, doesn't reorder
void SpanningTopology::add_span( SpanOP span ) {
	topology_.push_back( span );
}// add span

/// @brief Sort Spans
/// @details  Sort spans in ascending order from begin to end anchor
void SpanningTopology::reorder_spans() {
	std::sort( topology_.begin(), topology_.end() );
}// reorder spans

/// @brief Return Transmembrane Span
/// @details Return transmembrane span by it's index in the spanning topology object
SpanOP
SpanningTopology::span( Size span_number ){ return topology_[ span_number ]; } // get span

/// @brief Get total number of spans
/// @details Return the number of transmembrane spanning regions in this object
Size
SpanningTopology::total_spans() { return topology_.size(); } // total_spans

/// @brief Is the residue in the membrane region?
/// @details Return true if this residue is in a transmembrane span
bool SpanningTopology::in_span( Size resnum ){
    
    // go through spans and check whether residue is in spans
    for ( Size i = 1; i <= total_spans(); ++i ) {
        if ( resnum >= span(i)->start() && resnum <= span(i)->end() ){
            return true;
        }
    }
  
	// if not caught so far, return false
    return false;

} // in_span?

/// @brief Does the span cross the membrane
/// @details Determine if the membrane spanning region crosses the whole membrane
bool
SpanningTopology::spanning( utility::vector1< Real > res_z_coord, SpanOP span ){
	
	// z coordinates of start and end both negative?
	if ( res_z_coord[ span->start() ] < 0 && res_z_coord[ span->end() ] < 0 ) {
		return false;
	}
	
	// z coordinates of start and end both positive?
	if ( res_z_coord[ span->start() ] > 0 && res_z_coord[ span->end() ] > 0 ) {
		return false;
	}
	
	return true;
	
}// spanning?

/// @brief Determine if this Spanning Topology Object is Valid
/// @details Check that spans still span the membrane
bool SpanningTopology::is_valid(){
    
	bool valid( false );
	
    // check all spans
    for ( Size i = 1; i <= topology_.size(); ++i ){
		
		// if any spans invalid, return false
		if ( topology_[i]->is_valid() == false ){
			TR << "Span " << i << " is invalid!" << std::endl;
			return false;
		}
    }
	
    // check if spans are in increasing order
    for ( Size i = 1; i <= topology_.size(); ++i ){
        if ( i > 1 &&
            ( topology_[ i-1 ]->start() > topology_[ i ]->start() ||
             topology_[ i-1 ]->end() > topology_[ i ]->end() ) ){
                show();
                return false;
		}
    }
	
	// if it didn't crash until now, topology is valid
	valid = true;
    return valid;
}// is_valid?


/// @brief Return the number of residues represented by this topology object
Size SpanningTopology::nres_topo(){
	return nres_topo_;
}

//////////////////////
/// Helper Methods ///
//////////////////////

/// @brief Create spanning topology object from spanfile
SpanningTopology
SpanningTopology::create_from_spanfile( std::string spanfile, Size  ){
 
	// Setup vars for reading spanfile using izstream
    std::string line;
    utility::io::izstream stream ( spanfile );
    
    if ( !stream.good() ) {
		TR << "Poor formatting for spanfile - cannot read" << std::endl;
        throw new utility::excn::EXCN_Msg_Exception( "Poor formatting for spanfile - cannot read" );
    }
    
	// Read file Header "TM region prediction for"
    getline( stream, line );
    
		// Read line which includes number of tm spans and total resnum
    getline( stream, line );
    std::istringstream l( line );
    Size total_tmhelix;
    Size total_residues_in_span_file;
    
    l >> total_tmhelix >> total_residues_in_span_file;
	nres_topo_ = total_residues_in_span_file;

		// check number of residues
	//if ( nres != total_residues_in_span_file ){
	//	utility_exit_with_message(
	//			"SpanningTopology: Total_residues in span file " + /utility::to_string(total_residues_in_span_file) +
//				" does not match total_residues in pose " + utility::to_string(nres) );
//	}
    
	// If there is a negative number of total residues, throw exception
    if ( total_residues_in_span_file <= 0 ) {
        utility_exit_with_message( "SpanningTopology: No residues in pose - check file format or data." );
    }
    
	// If there are no tm helices, tell the user to stop using a globular protein
    if ( total_tmhelix <= 0 ) {
        utility_exit_with_message( "SpanningTopology: No TM helices in file - check file format." );
    }
    
	// Read in antiparallel and n2c line
    getline( stream, line );
    getline( stream, line );
    
	// For each line of the file, get spanning region info
    for ( Size i = 1; i <= total_tmhelix; ++i ) {
        
        getline( stream, line );
        std::istringstream l( line );
        Size start, end;
        l >> start >> end;
        
		// add to chain topology
        SpanOP span( new Span( start, end ) );
        topology_.push_back( span );
    }
    
	// Close the izstream
    stream.close();
    stream.clear();
		
	// check created object for validity
	if ( ! is_valid()){
		TR << "SpanningTopology invalid: check your span file!" << std::endl;
		throw utility::excn::EXCN_Msg_Exception( "SpanningTopology invalid: check your span file!" );
//		utility_exit_with_message( "SpanningTopology invalid: check your span file!" );
	}
	return *this;
} // create from spanfile

/// @brief Create Transmembrane SPan OBject from structure
SpanningTopology
SpanningTopology::create_from_structure(
	utility::vector1< Real > res_z_coord,
	utility::vector1< Size > chainID,
	Real thickness )
{
    
	// counter
    Size num_spans( 0 );
	Size start(1);
	Size end(1);
	nres_topo_ = res_z_coord.size();
	
	// cry if vectors are not the same length
	if ( res_z_coord.size() != chainID.size() ){
		utility_exit_with_message( "z_coord and chainID vectors are not the same length!" );
	}
	
	// At each z-coord, get start/end positions
    for ( Size j = 1; j <= res_z_coord.size()-1; ++j ){
		
		TR.Debug << "going through residue " << j << std::endl;
        
		// set start position for a new chain
		if ( j > 1 && chainID[ j ] != chainID[ j-1 ]){
			TR.Debug << "setting new start position: " << j << std::endl;
			start = j;
		}
		
		// get start position: from outside into membrane
        if ( ( chainID[ j ] == chainID[ j+1 ] ) &&
            (( res_z_coord[ j ] <= -thickness && res_z_coord[ j+1 ] > -thickness ) ||
             ( res_z_coord[ j ] >= thickness && res_z_coord[ j+1 ] < thickness )) ){
                
				TR.Debug << "Adding TMspan start at " << j+1 << std::endl;
				start = j+1;
		}
		TR.Debug << "start: " << start << std::endl;
		
		// get end position: from inside to out of the membrane
        if ( ( chainID[ j ] == chainID[ j+1 ] ) &&
            (( res_z_coord[ j ] >= -thickness && res_z_coord[ j+1 ] < -thickness ) ||
             ( res_z_coord[ j ] <= thickness && res_z_coord[ j+1 ] > thickness )) ){
                
				TR.Debug << "Adding TMspan end at " << j << std::endl;
				end = j;
				SpanOP span( new Span( start, end ) );
				
				// if span spans the membrane
				if ( spanning( res_z_coord, span ) ){
					TR << "Adding TMspan " << std::endl;
					span->show();
					topology_.push_back( span );
					++num_spans;
					start = end + 1;
				}
				if ( ! spanning( res_z_coord, span )){
					span->show();
					TR << "...thrown out because it doesn't span the membrane!" << std::endl;
				}
		}
		TR.Debug << "end: " << end << std::endl;
	}
		
	// If no helices predicted, throw an error
	TR.Debug << "Are spans predicted?" << std::endl;
	if ( num_spans == 0 ) {
		TR << "No spans calculated from structure" << std::endl;
		utility_exit_with_message( "No transmembrane helices predicted from pose!" );
	}
	
	// check created object for validity
	TR.Debug << "Are spans valid?" << std::endl;
	if ( ! is_valid()){
		TR << "SpanningTopology invalid" << std::endl;
		utility_exit_with_message( "SpanningTopology invalid: check your span file!" );
	}
	return *this;
	
} // create from structure

} // membrane
} // conformation
} // core
