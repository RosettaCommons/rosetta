// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///	@file		core/conformation/membrane/SpanningTopology.hh
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

#ifndef INCLUDED_core_conformation_membrane_SpanningTopology_hh
#define INCLUDED_core_conformation_membrane_SpanningTopology_hh

// Unit headers
#include <core/conformation/membrane/SpanningTopology.fwd.hh>
#include <core/conformation/membrane/Span.fwd.hh>

#ifdef WIN32
	#include <core/conformation/membrane/Span.hh>
#endif


// Package Header
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/xyzVector.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <cstdlib>
#include <string>

namespace core {
namespace conformation {
namespace membrane {

class SpanningTopology : public utility::pointer::ReferenceCount {

public: // constructors

    /// @brief	Default Constructor (Private)
	/// @details Construct an Empty Spanning Topology Object
    SpanningTopology();

    /// @brief	Custom Constructor - Transmembrane Spans from Spanfile
	/// @details Use transmembrane spans provided to consturct a spanning topology object
    SpanningTopology(
		std::string spanfile,
		Size total_residues = 0
	);

    /// @brief	Custom Constructor - Transmembrane Spans from xyz coords
	/// @details Use coordinates of residue CA and thickness to determine the spanning regions in the pose
    SpanningTopology(
		utility::vector1< Real > res_z_coord,
		utility::vector1< Size > chainID,
		Real thickness
		);

	/// @brief	Copy Constructor
	/// @details Create a deep copy of this object copying over all private fields
	SpanningTopology( SpanningTopology const & src );

	/// @brief Assignment Operator
	/// @details Overload assignemnt operator - required Rosetta method
	SpanningTopology &
	operator=( SpanningTopology const & src );

    /// @brief Destructor
    ~SpanningTopology();

public: // methods

	//////////////////////
	// OUTPUT FUNCTIONS //
	//////////////////////
    
	/// @brief  Generate string representation of Spanning Topology Object for debugging purposes.
	virtual void show( std::ostream & output=std::cout ) const;

	// write spanfile
	void write_spanfile( std::string output_filename );
	
	/////////////////////////
	// GETTERS AND SETTERS //
	/////////////////////////
	
    // get topology
    utility::vector1< SpanOP > get();
    
	// get number of spans
    Size nspans();

    // get span by number
    SpanOP span( Size span_number );
	
	// fill from spanfile - can be used after creating empty object
	void fill_from_spanfile( std::string spanfile, Size total_residues = 0 );
	
	// concatenate 2nd topology object
	SpanningTopologyOP concatenate_topology( SpanningTopologyOP topo );
	
	// add span to end of SpanningTopology object, doesn't reorder
	void add_span( SpanOP span, Size offset = 0 );

	// add span to end of SpanningTopology object, doesn't reorder
	void add_span( Size start, Size end, Size offset = 0 );

	// reorder spans, for instance after adding one
	void reorder_spans();
    	
	//////////////////
	// FOR CHECKING //
	//////////////////
	
    // is residue in membrane?
    bool in_span( Size residue );
	
	// does the span cross z=0, i.e. really spanning the membrane?
	bool spanning( utility::vector1< Real > res_z_coord, SpanOP span );

    /// @brief Determine if this Spanning Topology Object is Valid
	/// @details Check that spans still span the membrane
    bool is_valid();

	// return number of residues in spanfile - for checking
	Size nres_topo();

private: // methods

	/// @brief Create spanning topology object from spanfile
	SpanningTopology create_from_spanfile( std::string spanfile, Size nres);

	/// @brief Create Transmembrane SPan OBject from structure
	SpanningTopology create_from_structure( utility::vector1< Real > res_z_coord, utility::vector1< Size > chainID, Real thickness );

private: // data

    // vector of spans
    utility::vector1< SpanOP > topology_;

	// keep track of nres for checks
	Size nres_topo_;

}; // class SpanningTopology

} // membrane
} // conformation
} // core

#endif // INCLUDED_core_conformation_membrane_SpanningTopology_hh
