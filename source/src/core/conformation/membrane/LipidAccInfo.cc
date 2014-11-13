// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		core/conformation/membrane/LipidAccInfo.cc
///
/// @brief      Membrane Lipid Accessibility Data
/// @details    Object for storing per-residue lipid exposed and buried surface
///				area values. Predicted from sequence, transmembrane spans, and psiblast
///				prediction using server called from the run_lips.pl script.
///				Last Modified: 7/7/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

// Unit headers
#include <core/conformation/membrane/LipidAccInfo.hh>

// Package Headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <cmath>

static thread_local basic::Tracer TR( "core.conformation.membrane.LipidAccInfo" );

/// @brief      Membrane Lipid Accessibility Data
/// @details    Stores lipid accessibility data derived from OCTOPUS spanning file
///             and psiblast search using run_lips.pl script

namespace core {
namespace conformation {
namespace membrane {
    
/// @brief Constructor
/// @details Create a blank copy of the lipid accessibility data object
LipidAccInfo::LipidAccInfo() :
	utility::pointer::ReferenceCount()
{}

/// @brief Custom Constructor
/// @brief Construct from user-provided lipid Acc Info File
LipidAccInfo::LipidAccInfo( std::string lipsfile ) {
	
	TR << "Initializing lips exposure info using " << lipsfile << std::endl;
	
	// Initialize local vars
	Size num_of_csts( 0 );
	Size resnum;
	Real exposure;
	
	// Initialize izstream
	std::string line;
	utility::io::izstream stream ( lipsfile );
	
	// Open stream and start reading
	stream.open( lipsfile );
	if ( stream ) {
		
		// Grab the first line
		getline( stream, line );
		getline( stream, line );
		
		// Grab total residues
		int nres = std::atoi( line.c_str() );

		lipid_exposure_.resize( nres, 0.0 );
		lipid_burial_.resize( nres, 0.0 );
		
		while ( !stream.eof() ) {
			
			std::istringstream l( line );
			
			// Add exposure info to vectors
			l >> resnum;
			l >> exposure;
			
			if ( exposure > 0 ) {
				
				lipid_exposure_[ resnum ] = exposure;
				num_of_csts++;
				
			} else {
				
				lipid_burial_[ resnum ] = std::fabs( exposure );
				num_of_csts++;
			}
			
			getline( stream, line );
		}
		
		stream.close();
		stream.clear();

	} else {
		utility_exit_with_message( "Lips data file not found!" );
	}
}

/// @brief Conpy Constructor
/// @details Create a deep copy of this object
LipidAccInfo::LipidAccInfo( LipidAccInfo const & src ) :
	utility::pointer::ReferenceCount()
{
	copy_data(*this, src);
}

/// @brief Assignment Operator
/// @details Create a deep copy of this object, overloading the assignment operator
LipidAccInfo &
LipidAccInfo::operator=( LipidAccInfo const & src ) {
	
	// Abort self-assignment.
	if ( this == &src ) {
		return *this;
	}
	
    // Make a deep copy of data in this object
    copy_data(*this, src);
    return *this;
	
}

/// @brief Destructor
LipidAccInfo::~LipidAccInfo() {}

/// @brief Access Lipid exposed surface area per-residue
utility::vector1< core::Real > LipidAccInfo::lipid_exposure() { return lipid_exposure_; }

/// @details Access Lipid buried surface area per-residue
utility::vector1< core::Real > LipidAccInfo::lipid_burial() { return lipid_burial_; }

/// @brief Copy Data
void
LipidAccInfo::copy_data( LipidAccInfo src, LipidAccInfo copy ) {
	
	src.lipid_exposure_ = copy.lipid_exposure_;
	src.lipid_burial_ = copy.lipid_burial_;
}
    
    
} // membrane
} // conformation
} // core

