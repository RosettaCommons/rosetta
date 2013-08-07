// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RNA de novo fragment assembly Structure Parameters
/// @brief User input parameters for RNA structure inference
/// @detailed
/// @author Rhiju Das


// Unit Headers
#include <protocols/rna/RNA_DataReader.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/RNA_DataInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>

// ObjexxFCL Headers
// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>


namespace protocols {
namespace rna {

/// @details Auto-generated virtual destructor
RNA_DataReader::~RNA_DataReader() {}

static basic::Tracer tr( "protocols.rna.rna_structure_parameters" ) ;

using namespace core;

RNA_DataReader::RNA_DataReader(){
}

//////////////////////////////////////////////////////////////////
void
RNA_DataReader::initialize(
	core::pose::Pose & pose,
	std::string const rna_data_file		 )
{

	if ( rna_data_file.length() == 0 ) return;

	Size const nres( pose.total_residue() );

	backbone_burial_.dimension( nres, false );
	backbone_exposed_.dimension( nres, false );

	read_data_from_file( rna_data_file );

	for (Size i = 1; i <= nres; i++ ) {
		if ( backbone_burial_( i ) && backbone_exposed_( i ) ) {
			utility_exit_with_message(  "Cannot ask for both backbone exposure and burial: "+ObjexxFCL::string_of(i) );
		}
	}

	setup_rna_data( pose );

}


/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_burial_info( std::istringstream & line_stream, ObjexxFCL::FArray1D< bool > & array_to_fill ) {

	Size pos( 0 );

	while ( !line_stream.fail() ) {

		line_stream >> pos;

		if ( pos > array_to_fill.size() ) {
			utility_exit_with_message(  "Problem with data format on BACKBONE_BURIAL line. " );
		}

		array_to_fill( pos ) = true;

	}

}



/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_data_info( std::istringstream & line_stream ) {

	using namespace core::scoring::rna;

	Size pos( 0 );
	char edge( 'X' );
	Real weight( 0.0 );

	while ( !line_stream.fail() ) {

		line_stream >> pos >> edge >> weight;

		if ( edge == 'W' ){
			rna_data_info_.add_datum( RNA_Datum( pos, WATSON_CRICK, weight ) );
		} else if ( edge == 'H' ) {
			rna_data_info_.add_datum( RNA_Datum( pos, HOOGSTEEN, weight ) );
		} else if ( edge == 'S' ) {
			rna_data_info_.add_datum( RNA_Datum( pos, SUGAR, weight ) );
		} else if ( edge == 'X' ) {
			rna_data_info_.add_datum( RNA_Datum( pos, WATSON_CRICK, weight ) );
			rna_data_info_.add_datum( RNA_Datum( pos, HOOGSTEEN, weight ) );
			rna_data_info_.add_datum( RNA_Datum( pos, SUGAR, weight ) );
		} else {
			utility_exit_with_message(  "Problem with data format on DATA line. " );
		}
	}

}



/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_data_from_file( std::string const & filename ) {

	tr << "Reading RNA data file: " << filename << std::endl;

	//	Size a;

	utility::io::izstream data_stream( filename );
	if ( !data_stream ) {
		data_stream.close();
		utility_exit_with_message(  "Data file? " + filename );
	}

	std::string line, tag;

	while ( getline( data_stream, line ) ) {

		std::istringstream line_stream( line );
		line_stream >> tag;

		if (line_stream.fail() ) continue; //Probably a blank line.

		if  ( tag == "EXPOSE" ) {

			read_data_info( line_stream );

		} else if ( tag == "BACKBONE_BURIAL" ) {

			read_burial_info( line_stream, backbone_burial_ );

		} else if ( tag == "BACKBONE_EXPOSE" ) {

			read_burial_info( line_stream, backbone_exposed_ );

		} else if  ( tag[0] == '#' ) {
			continue;

		} else {
			utility_exit_with_message(   "Unrecognized tag in data file: " + tag );
		}
	}
	data_stream.close();

}


/////////////////////////////////////////////////////////////////////
void
RNA_DataReader::setup_rna_data( core::pose::Pose & pose )
{
	rna_data_info_.set_backbone_burial(  backbone_burial_ );
	rna_data_info_.set_backbone_exposed(  backbone_exposed_ );

	core::scoring::rna::RNA_ScoringInfo & rna_scoring_info( core::scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	core::scoring::rna::RNA_DataInfo & rna_data_info( rna_scoring_info.rna_data_info() );
	rna_data_info = rna_data_info_;

	//	rna_data_info.set_backbone_burial( backbone_burial_ );

}

}
}



