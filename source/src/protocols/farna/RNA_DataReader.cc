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
#include <protocols/farna/RNA_DataReader.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <fstream>
#include <iostream>

#include <utility/vector1.hh>

///////////////////////////////////////////////////////////////////////////////////
// This was written a while back, and was meant for experimental testing of
// RNA data input from single-nucleotide-resolution chemical mapping experiments.
//
// In 2014, now upgrading to allow direct input of reactivity values, since we
// finally have solid ways to measure those & normalize precisely.
//
// May soon improve further to allow input of RDAT file format, which comes out
// of HiTRACE/MAPseeker toolkits developed in Das laboratory.
//
// later might be better to index with conventional numbering & chains.
//
// Does this really need pose to initialize? better to just store user-inputted
//  data without cross-checks requiring pose.
//
///////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace farna {

/// @details Auto-generated virtual destructor
RNA_DataReader::~RNA_DataReader() {}

static basic::Tracer TR( "protocols.rna.RNA_DataReader" ) ;

using namespace core;

RNA_DataReader::RNA_DataReader(){
}

//////////////////////////////////////////////////////////////////
void
RNA_DataReader::initialize(
	pose::Pose & pose,
	std::string const rna_data_file		 )
{
	using namespace pose::full_model_info;

	if ( rna_data_file.length() == 0 ) return;
	rna_data_info_ = new core::scoring::rna::data::RNA_DataInfo;

	Size const nres( pose.total_residue() );
	if ( full_model_info_defined( pose ) ) full_model_parameters_ = const_full_model_info( pose ).full_model_parameters();

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

	int pos( 0 );

	while ( !line_stream.fail() ) {

		line_stream >> pos;
		if ( full_model_parameters_ ) pos = full_model_parameters_->conventional_to_full( pos );
		if ( pos > int( array_to_fill.size() ) ) {
			utility_exit_with_message(  "Problem with data format on BACKBONE_BURIAL line. " );
		}

		array_to_fill( pos ) = true;

	}

}



/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_data_info( std::istringstream & line_stream ) {

	using namespace scoring::rna::data;

	int pos( 0 );
	char edge( 'X' );
	Real weight( 0.0 );

	while ( !line_stream.fail() ) {

		line_stream >> pos >> edge >> weight;
		if ( full_model_parameters_ ) pos = full_model_parameters_->conventional_to_full( pos );

		if ( edge == 'W' ){
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::WATSON_CRICK, weight ) );
		} else if ( edge == 'H' ) {
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::HOOGSTEEN, weight ) );
		} else if ( edge == 'S' ) {
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::SUGAR, weight ) );
		} else if ( edge == 'X' ) {
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::WATSON_CRICK, weight ) );
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::HOOGSTEEN, weight ) );
			rna_data_info_->add_datum( RNA_Datum( pos, chemical::rna::SUGAR, weight ) );
		} else {
			utility_exit_with_message(  "Problem with data format on DATA line. " );
		}
	}

}


/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_reactivity_info( std::istringstream & line_stream, scoring::rna::data::RNA_ReactivityType const type ) {

	int pos;
	Real value;
	line_stream >> pos >> value;
	if ( full_model_parameters_ ){
		if ( !full_model_parameters_->has_conventional_residue( pos ) ) return;
		pos = full_model_parameters_->conventional_to_full( pos );
	}
	rna_data_info_->add_reactivity( scoring::rna::data::RNA_Reactivity( pos, type, value ) );

}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_data_from_file( std::string const & filename ) {

	TR << "Reading RNA data file: " << filename << std::endl;

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

		} else if ( tag == "DMS" ) {

			read_reactivity_info( line_stream, scoring::rna::data::DMS );

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
RNA_DataReader::setup_rna_data( pose::Pose & pose )
{
	rna_data_info_->set_backbone_burial(  backbone_burial_ );
	rna_data_info_->set_backbone_exposed(  backbone_exposed_ );

	scoring::rna::RNA_ScoringInfo & rna_scoring_info( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	scoring::rna::data::RNA_DataInfo & rna_data_info( rna_scoring_info.rna_data_info() );
	rna_data_info = *rna_data_info_;

	//	rna_data_info.set_backbone_burial( backbone_burial_ );

}


/////////////////////////////////////////////////////////////////////
core::scoring::rna::data::RNA_DataInfo const &
get_rna_data_info( pose::Pose & pose, std::string const & rna_data_file ) {

	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace core::scoring::rna;
	using namespace core::scoring::rna::data;

	if ( rna_data_file.size() > 0 ) {
		RNA_DataReaderOP rna_data_reader = new RNA_DataReader;
		rna_data_reader->initialize( pose, rna_data_file );
	}

	RNA_DataInfo const & rna_data_info = scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info();

	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ){
		nonconst_rna_scoring_info_from_pose( *(other_pose_list[ n ]) ).rna_data_info() = rna_data_info;
	}

	return rna_data_info;
}



} //farna
} //protocols
