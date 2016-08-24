// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/conformation/membrane/SpanningTopology.cc
///
/// @brief      Transmembrane Spans Object
/// @details Object for storing membrane spanning regions as a vector1 of span objects.
///    Spanning regions can be determined either from an input spanfile, xyz coordinates
///    in the pose, or sequence. Object is constructed from span regions and will do internal
///    checking for validity.
///    Last Modified: 3/18/15
///
/// @author  Julia Koehler (julia.koehler1982@gmail.com)
/// @author  Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/conformation/membrane/SpanningTopology.hh>
#include <core/conformation/membrane/Span.hh>

// Project Headers
//#include <core/conformation/membrane/Exceptions.hh>
// This class isn't actually used--it's just an indirect way of including
#include <utility/excn/Exceptions.hh>

// Package Headers
#include <core/types.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <numeric/numeric.functions.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <string>
#include <map>

static THREAD_LOCAL basic::Tracer TR( "core.conformation.membrane.SpanningTopology" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace conformation {
namespace membrane {

////////////////////
/// Constructors ///
////////////////////

/// @brief Default Constructor (Private)
/// @details Construct an Empty Spanning Topology Object
SpanningTopology::SpanningTopology() :
	utility::pointer::ReferenceCount(),
	topology_(),
	nres_topo_( 0 ),
	structure_based_( false )
{}

/// @brief Custom Constructor - Transmembrane Spans from Spanfile
/// @details Use transmembrane spans provided to consturct a spanning topology object
SpanningTopology::SpanningTopology(
	std::string spanfile,
	std::map< std::string, core::Size > pdb2pose_map,
	core::Size total_residues
) : utility::pointer::ReferenceCount(),
	topology_(),
	nres_topo_( 0 ),
	structure_based_( false )
{
	create_from_spanfile( spanfile, pdb2pose_map, total_residues );

} // topology from spanfile

/// @brief Custom Constructor - Transmembrane Spans from xyz coords
/// @details Use coordinates of residue CA and thickness to determine the spanning regions in the pose
SpanningTopology::SpanningTopology(
	utility::vector1< core::Real > res_z_coord,
	utility::vector1< core::Size > chainID,
	utility::vector1< char > secstruct,
	core::Real thickness
) : utility::pointer::ReferenceCount(),
	topology_(),
	nres_topo_( 0 ),
	structure_based_( true )
{
	create_from_structure( res_z_coord, chainID, secstruct, thickness);

} // topology from pose and thickness

/// @brief Copy Constructor
/// @details Create a deep copy of this object copying over all private fields
SpanningTopology::SpanningTopology( SpanningTopology const & src ) :
	utility::pointer::ReferenceCount( src ),
	topology_( src.topology_ ),
	nres_topo_( src.nres_topo_ ),
	structure_based_( src.structure_based_ )
{}

/// @brief Assignment Operator
/// @details Overload assignemnt operator - required Rosetta method
SpanningTopology &
SpanningTopology::operator=( SpanningTopology const & src ) {

	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}

	// Deep Copy of the data
	this->topology_ = src.topology_;
	this->nres_topo_ = src.nres_topo_;
	this->structure_based_ = src.structure_based_;

	return *this;
}

/// @brief Destructor
SpanningTopology::~SpanningTopology(){}

///////////////
/// Methods ///
///////////////

/// @brief fill from spanfile
/// @details fill object from spanfile, can be used after creating empty object
void SpanningTopology::fill_from_spanfile( std::string spanfile, std::map< std::string, core::Size > pdb2pose_map, core::Size total_residues ){

	TR << "Filling membrane spanning topology from spanfile " << spanfile << std::endl;
	create_from_spanfile( spanfile, pdb2pose_map, total_residues );

}// fill from spanfile

// fill from structure - can be used after creating empty object
void SpanningTopology::fill_from_structure( utility::vector1< core::Real > res_z_coord,
	utility::vector1< core::Size > chainID,
	utility::vector1< char > secstruct,
	Real thickness ) {
	TR << "Filling membrane spanning topology from structure with thickness " << thickness << std::endl;
	create_from_structure( res_z_coord, chainID, secstruct, thickness );

	this->show( TR );

	TR << "WATCH OUT: Writing spanfile out.span!" << std::endl;
	write_spanfile( "out.span" );

}// fill from structure

/// @brief Generating a String Representation of Spanning Topology Object for debugging purposes
void SpanningTopology::show( std::ostream & ) const {

	if ( nres_topo_ == 0 || topology_.size() == 0 ) {
		utility_exit_with_message("Can't print SpanningTopology. It's empty!");
	}

	// Print Total number of transmembrane Spans
	TR << "Total # of TM spans: " << topology_.size() << std::endl;
	TR << "Number of residues in spanfile: " << nres_topo_ << std::endl;

	// print individual spans
	for ( core::Size i = 1; i <= topology_.size(); ++i ) {
		TR << "Span " << i << ": start: " << topology_[ i ]->start();
		TR << ", end: " << topology_[ i ]->end() << std::endl;
	}

} // show

// concatenate 2nd topology object
SpanningTopology & SpanningTopology::concatenate_topology( SpanningTopology const & topo ){

	// add spans
	for ( core::Size i = 1; i <= topo.nspans(); ++i ) {
		SpanOP span( topo.span( i ) );
		span->shift( nres_topo_ );
		topology_.push_back( span );
	}

	// update nres of topo1 object
	nres_topo_ += topo.nres_topo();

	return *this;
}

//////////////////////////////////////////////////////////////////////////////

// write spanfile
void SpanningTopology::write_spanfile( std::string output_filename ) const {

	TR.Debug << "printing spanfile" << std::endl;

	// print header
	utility::io::ozstream oz;
	oz.open( output_filename );
	oz << "Rosetta-generated spanfile from SpanningTopology object" << std::endl;
	oz << topology_.size() << " " << nres_topo_ << std::endl;
	oz << "antiparallel" << std::endl;
	oz << "n2c" << std::endl;

	// print spans
	for ( core::Size i = 1; i <= topology_.size(); ++i ) {
		oz << "\t" << topology_[i]->start() << "\t" << topology_[i]->end() << std::endl;
	}
	oz.close();
	TR << "wrote " << output_filename << std::endl;

}// write spanfile


/// @brief Return Spanning Topology
/// @details return spanning topology as a vector1 of transmembrane spans
utility::vector1< SpanOP >
SpanningTopology::get_spans() const { return topology_; } // get topology

// add span to end of SpanningTopology object, doesn't reorder
void SpanningTopology::add_span( Span const & span, core::Size offset ){
	SpanOP new_span ( new Span( span ) );
	new_span->shift( offset );
	topology_.push_back( new_span );

	if ( ! structure_based_ ) {
		nres_topo_ += span.end() - span.start() + 1;
	}
}// add span

//////////////////////////////////////////////////////////////////////////////

// add span to end of SpanningTopology object, doesn't reorder
void SpanningTopology::add_span( core::Size start, core::Size end, core::Size offset ){
	SpanOP span( new Span( start+offset, end+offset ) );
	topology_.push_back( span );

	if ( ! structure_based_ ) {
		nres_topo_ += end - start + 1;
	}
}// add span

//////////////////////////////////////////////////////////////////////////////

/// @brief Sort Spans
/// @details  Sort spans in ascending order from begin to end anchor
void SpanningTopology::reorder_spans() {
	std::sort( topology_.begin(), topology_.end() );
}// reorder spans

//////////////////////////////////////////////////////////////////////////////

/// @brief Return Transmembrane Span
/// @details Return transmembrane span by it's index in the spanning topology object
SpanOP SpanningTopology::span( Size span_number ) const {
	return topology_[ span_number ];
} // get span

//////////////////////////////////////////////////////////////////////////////

/// @brief Get total number of spans
/// @details Return the number of transmembrane spanning regions in this object
Size
SpanningTopology::nspans() const { return topology_.size(); } // total_spans

//////////////////////////////////////////////////////////////////////////////

/// @brief Is the residue in the membrane region?
/// @details Return true if this residue is in a transmembrane span
bool SpanningTopology::in_span( core::Size resnum ) const {

	// go through spans and check whether residue is in spans
	for ( core::Size i = 1; i <= nspans(); ++i ) {
		if ( resnum >= span(i)->start() && resnum <= span(i)->end() ) {
			return true;
		}
	}

	// if not caught so far, return false
	return false;

} // in_span?

//////////////////////////////////////////////////////////////////////////////

/// @brief Does the span cross the membrane
/// @details Determine if the membrane spanning region crosses the whole membrane
bool
SpanningTopology::spanning( utility::vector1< core::Real > res_z_coord, Span const & span ) const {

	// z coordinates of start and end both negative?
	if ( res_z_coord[ span.start() ] < 0 && res_z_coord[ span.end() ] < 0 ) {
		return false;
	}

	// z coordinates of start and end both positive?
	if ( res_z_coord[ span.start() ] > 0 && res_z_coord[ span.end() ] > 0 ) {
		return false;
	}

	return true;

}// spanning?

//////////////////////////////////////////////////////////////////////////////

/// @brief Determine if this Spanning Topology Object is Valid
/// @details Check that spans still span the membrane
bool SpanningTopology::is_valid() const {

	bool valid( false );

	// check all spans
	for ( core::Size i = 1; i <= topology_.size(); ++i ) {

		// if any spans invalid, return false
		if ( topology_[i]->is_valid() == false ) {
			TR << "Span " << i << " is invalid!" << std::endl;
			return false;
		}
	}

	// check if spans are in increasing order
	for ( core::Size i = 1; i <= topology_.size(); ++i ) {
		if ( i > 1 &&
				( topology_[ i-1 ]->start() > topology_[ i ]->start() ||
				topology_[ i-1 ]->end() > topology_[ i ]->end() ) ) {
			show( TR );
			return false;
		}
	}

	// if it didn't crash until now, topology is valid
	valid = true;
	return valid;
}// is_valid?

//////////////////////////////////////////////////////////////////////////////

/// @brief Return the number of residues represented by this topology object
Size SpanningTopology::nres_topo() const {
	return nres_topo_;
}

//////////////////////
/// Helper Methods ///
//////////////////////

/// @brief Create spanning topology object from spanfile
SpanningTopology
SpanningTopology::create_from_spanfile( std::string spanfile, std::map< std::string, core::Size > pdb2pose_map, core::Size nres ){

	using namespace utility;

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
	core::Size total_tmhelix;
	core::Size total_residues_in_span_file;

	l >> total_tmhelix >> total_residues_in_span_file;
	nres_topo_ = total_residues_in_span_file;

	TR << "nres: " << nres << std::endl;
	TR << "nres_topo: " << nres_topo() << std::endl;
	TR << "total residues in spanfile: " << total_residues_in_span_file << std::endl;

	// check number of residues: the pose sometimes has an additional membrane residue
	// if ( nres > 0 &&
	//  ( total_residues_in_span_file != nres - 1 ||
	//    total_residues_in_span_file != nres ||
	//    total_residues_in_span_file != nres + 1 )){
	//  utility_exit_with_message(
	//    "SpanningTopology: Total_residues in span file " + utility::to_string(total_residues_in_span_file) +
	//    " does not match total_residues in pose " + utility::to_string(nres) );
	// }

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

	// numbering scheme
	bool pdb_numbering;
	if ( pdb2pose_map.empty() ) {
		pdb_numbering = false;
	} else {
		pdb_numbering = true;
	}

	// For each line of the file, get spanning region info
	for ( core::Size i = 1; i <= total_tmhelix; ++i ) {

		getline( stream, line );
		std::istringstream l( line );

		// tag that is read in
		std::string start_tag, end_tag;
		l >> start_tag >> end_tag;

		// initialize residue numbers in pose numbering that will be used to
		// create a Span object
		core::Size start( 0 );
		core::Size end( 0 );

		// first TM span tag decides whether it is PDB numbering or pose numbering
		// if first character is numeric, then pose numbering
		if ( i == 1 ) {
			if ( isdigit( start_tag[ 0 ] ) ) {
				TR << "Pose is in pose numbering scheme." << std::endl;
				pdb_numbering = false;
			} else {
				TR << "Pose is in PDB numbering scheme." << std::endl;
			}
		}

		// NEW: supports PDB numbering format with chains and insertion codes
		// format: A1. for chain A (single character), residue number 1, no insertion code
		// check for chain
		if ( pdb_numbering ) {

			// check for chain
			if ( ! isalpha( start_tag[ 0 ] )
					|| ! isalpha( end_tag[ 0 ] ) ) {
				utility_exit_with_message( "Cannot read span input. Should either be in pose numbering (renumbered to start with 1 without gaps) without alphabetic characters or in PDB numbering with the format A15C for chain A, residue 15, and insertion code C. If insertion code is empty, use a dot. Seems like your chain is incorrect." );
			}

			// check for insertion code
			if ( isdigit( start_tag[ -1 ] )
					|| isdigit( end_tag[ -1 ] ) ) {
				utility_exit_with_message( "Cannot read span input. Should either be in pose numbering (renumbered to start with 1 without gaps) without alphabetic characters or in PDB numbering with the format A15C for chain A, residue 15, and insertion code C. If insertion code is empty, use a dot. Seems like your insertion code is incorrect." );
			}

			// get pose numbering from map
			start = pdb2pose_map[ start_tag ];
			end = pdb2pose_map[ end_tag ];

		} else {
			// also supports old format in pose numbering
			// format: 1 for residue number 1, no insertion code, chain is whatever
			// it is for that residue
			// if we assume PDB format, pose numbering is also handled correctly because
			// the PDB is already renumbered

			start = from_string( start_tag, core::Size( 0 ) );
			end = from_string( end_tag, core::Size( 0 ) );

		}

		// add to chain topology
		SpanOP span( new Span( start, end ) );
		topology_.push_back( span );
	}

	// Close the izstream
	stream.close();
	stream.clear();

	// check created object for validity
	if ( ! is_valid() ) {
		TR << "SpanningTopology invalid: check your span file!" << std::endl;
		throw utility::excn::EXCN_Msg_Exception( "SpanningTopology invalid: check your span file!" );
		//  utility_exit_with_message( "SpanningTopology invalid: check your span file!" );
	}

	this->show( TR );
	return *this;
} // create from spanfile

//////////////////////////////////////////////////////////////////////////////

/// @brief Create span object from structure
SpanningTopology
SpanningTopology::create_from_structure(
	utility::vector1< core::Real > res_z_coord,
	utility::vector1< core::Size > chainID,
	utility::vector1< char > secstruct,
	Real thickness )
{
	TR << "create topology from structure" << std::endl;

	// variables
	nres_topo_ = res_z_coord.size();
	utility::vector1< bool > spans;

	// cry if vectors are not the same length
	if ( res_z_coord.size() != chainID.size() || secstruct.size() != chainID.size() ) {
		utility_exit_with_message( "Input vectors for create_from_structure are not the same length!" );
	}

	// At each z-coord, get membrane position, then add to spans vector
	for ( core::Size j = 1; j <= res_z_coord.size(); ++j ) {

		TR.Debug << "going through residue " << j << std::endl;

		// take all TM positions that are either helix or strand
		if ( -thickness <= res_z_coord[ j ] && res_z_coord[ j ] <= thickness &&
				secstruct[ j ] != 'L' ) {

			spans.push_back( true );
		} else {
			spans.push_back( false );
		}
	}

	// more variables
	core::Size num_spans( 0 );
	core::Size start(1);
	core::Size end(1);
	core::Size prev_start(0);
	core::Size prev_end(0);
	core::Real dir_prev(1.0);

	// go through newly written vector and identify spans
	for ( core::Size j = 1; j <= spans.size()-1; ++j ) {

		// add start
		if ( spans[ j ] == false && spans[ j+1 ] == true ) {
			TR << "TMspan start at " << j+1 << std::endl;
			start = j+1;
		}

		// add end
		if ( ( spans[ j ] == true && spans[ j+1 ] == false ) ||
				( spans[ j ] == true && j == spans.size()-1 ) ) {

			TR << "TMspan end at " << j << std::endl;
			end = j;

			// compute vector between previous span and this span
			// if the vectors are "antiparallel", keep the loop in between because these
			// are different spans
			if ( prev_start != 0 ) {
				dir_prev = res_z_coord[ prev_end ] - res_z_coord[ prev_start ];
			}
			core::Real dir_this = res_z_coord[ end ] - res_z_coord[ start ];
			bool parallel( false );

			// check if they are "parallel"
			if ( ( dir_prev < 0 && dir_this < 0 ) ||
					( dir_prev > 0 && dir_this > 0 ) ) {
				parallel = true;
			}

			// compute projection: if to short, span is amphipathic and won't be added
			core::Size length = end - start + 1;
			core::Real z_dist = numeric::abs_difference( res_z_coord[ end ], res_z_coord[ start ] );

			// if span is longer than 3 residues, and it's not ampipathic, it will be added
			if ( length >= 3 && z_dist > 5 ) {

				TR << "prev_end: " << prev_end << ", end: " << end << std::endl;

				// last span will be extended if there are only 2 residue in between and
				// the span vectors aren't facing the same direction
				if ( nspans() > 0 && prev_end >= start-4 && parallel == true ) {

					topology_.pop_back();
					TR << "Extending TMspan from " << prev_start << " to " << end << std::endl;
					add_span( prev_start, end );
					prev_start = start;
					prev_end = end;
				} else {

					TR << "Adding TMspan from " << start << " to " << end << std::endl;
					add_span( start, end );
					++num_spans;
					prev_start = start;
					prev_end = end;
				}
			}
		}

		// set new start position for each chain
		if ( j > 1 && chainID[ j ] != chainID[ j-1 ] ) {
			TR.Debug << "setting new start position: " << j << std::endl;
			end = j-1;
			start = j;
		}
	}

	// If no helices predicted, throw a warning
	// we don't want it to quit because MPs with large helical domains also need
	// split spanfiles from the spanfile_from_pdb application
	TR.Debug << "Are spans predicted?" << std::endl;
	if ( num_spans == 0 ) {
		TR << "WARNING: No spans calculated from structure for chain " << chainID[1] << std::endl;
	}

	// check created object for validity
	TR.Debug << "Are spans valid?" << std::endl;
	if ( ! is_valid() ) {
		TR << "SpanningTopology invalid" << std::endl;
		utility_exit_with_message( "SpanningTopology invalid: check your span file!" );
	}

	this->show( TR );
	return *this;

} // create from structure

//////////////////////////////////////////////////////////////////////////////

/// @brief Show Spanning topology
/// @details For PyRosetta!
std::ostream & operator << ( std::ostream & os, SpanningTopology const & spans ) {

	// Print Total number of transmembrane Spans
	os << "Total # of TM spans: " << spans.nspans() << std::endl;

	// print individual spans
	for ( core::Size i = 1; i <= spans.nspans(); ++i ) {
		os << "Span " << i << ": start: " << spans.get_spans()[ i ]->start();
		os << ", end: " << spans.get_spans()[ i ]->end() << std::endl;
	}
	return os;
}

} // membrane
} // conformation
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::conformation::membrane::SpanningTopology::save( Archive & arc ) const {
	arc( CEREAL_NVP( topology_ ) ); // utility::vector1<SpanOP>
	arc( CEREAL_NVP( nres_topo_ ) ); // Size
	arc( CEREAL_NVP( structure_based_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::conformation::membrane::SpanningTopology::load( Archive & arc ) {
	arc( topology_ ); // utility::vector1<SpanOP>
	arc( nres_topo_ ); // Size
	arc( structure_based_ ); //
}

SAVE_AND_LOAD_SERIALIZABLE( core::conformation::membrane::SpanningTopology );
CEREAL_REGISTER_TYPE( core::conformation::membrane::SpanningTopology )

CEREAL_REGISTER_DYNAMIC_INIT( core_conformation_membrane_SpanningTopology )
#endif // SERIALIZATION
