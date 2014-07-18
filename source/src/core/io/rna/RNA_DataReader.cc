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
#include <core/io/rna/RNA_DataReader.hh>

// Package Headers
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/scoring/rna/data/RNA_DataInfo.hh>
#include <core/io/rna/RDAT.hh>

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
// finally have solid ways to measure those & normalize precisely. Currently 3 styles of input
// are accepted. The best, associated with a careful log-odds potential is:
//
//  DMS <pos> <value>
//
// but is only defined for normalized DMS (dimethyl sulfate) reactivity values, at the moment.
// See http://dx.doi.org/10.1021/bi5003426.
//
// Another format requires user to guess bonus/penalty for pairing
// a base edge (Watson/Crick W, Hoogsteen H, Sugar S, or unknown X) during fragment assembly:
//
//  EXPOSE <pos> <edge: W,H,S,or X> <weight>
//
// There is one other style of input:
//
//  BACKBONE_EXPOSE <pos>
//  BACKBONE_BURIAL <pos>
//
// that was tested briefly for incorporation of *OH radical cleavage, but will be deprecated after implementing
// a solid log-odds potential that takes into account quantitative radical reactivity values.
//
// Positions will be read in as 'conventional' numbering, but then converted to full_model numbering
//  when pose is filled, if the pose has full_model_info is defined.
//
// May soon improve further to allow input of RDAT file format, which comes out
// of HiTRACE/MAPseeker toolkits developed in the Das laboratory
//  (http://rmdb.stanford.edu/repository/specs/).
//
//  -- rhiju, 2014.
//
///////////////////////////////////////////////////////////////////////////////////

using namespace core::scoring::rna::data;

namespace core {
namespace io {
namespace rna {

/// @details Auto-generated virtual destructor
RNA_DataReader::~RNA_DataReader() {}

static basic::Tracer TR( "protocols.rna.RNA_DataReader" ) ;

using namespace core;

RNA_DataReader::RNA_DataReader( std::string const rna_data_file )
{
	initialize( rna_data_file );
}

//////////////////////////////////////////////////////////////////
void
RNA_DataReader::initialize(	std::string const rna_data_file )
{

	rna_data_info_with_conventional_numbering_ = new RNA_DataInfo;
	backbone_burial_res_.clear();
	backbone_exposed_res_.clear();

	if ( rna_data_file.length() == 0 ) return;

	read_data_from_file( rna_data_file );

	for (Size i = 1; i <= backbone_burial_res_.size(); i++ ){
		runtime_assert( !backbone_exposed_res_.has_value( backbone_burial_res_[i] ) );
	}

}


//////////////////////////////////////////////////////////////////
// this style of input ('BACKBONE_EXPOSE', 'BACKBONE_BURIAL') was tested briefly for
// incorporation of *OH radical cleavage, but will be deprecated after implementing
// a solid log-odds potential that takes into account quantitative reactivity values, rather
// than boolean guesses.
void
RNA_DataReader::read_backbone_info( std::istringstream & line_stream,
																		utility::vector1< Size > & backbone_res ) {

	int pos( 0 ); // later allow this to be chain:res, using get_resnum_and_chain in string_util.hh
	while ( !line_stream.fail() ) {
		line_stream >> pos;
		backbone_res.push_back( pos );
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
		if ( edge == 'W' ){
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::WATSON_CRICK, weight ) );
		} else if ( edge == 'H' ) {
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::HOOGSTEEN, weight ) );
		} else if ( edge == 'S' ) {
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::SUGAR, weight ) );
		} else if ( edge == 'X' ) {
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::WATSON_CRICK, weight ) );
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::HOOGSTEEN, weight ) );
			rna_data_info_with_conventional_numbering_->add_datum( RNA_Datum( pos, chemical::rna::SUGAR, weight ) );
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
	rna_data_info_with_conventional_numbering_->add_reactivity( scoring::rna::data::RNA_Reactivity( pos, type, value ) );

}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::read_data_from_rdat( std::string const & filename ){
	RDAT rdat( filename );
	get_reactivity_from_rdat( rdat, DMS, "DMS" );
}

/////////////////////////////////////////////////////////////////////////////////////
void
RNA_DataReader::get_reactivity_from_rdat( core::io::rna::RDAT const & rdat,
																				  core::scoring::rna::data::RNA_ReactivityType const & type,
																					std::string const modifier_name ) {
	Size idx( 0 );
	if ( get_tag( rdat.annotations(), "modifier" ) == modifier_name ) {
		idx = 1;
	} else {
		utility::vector1< std::string > modifiers = get_tags( rdat.data_annotations(), "modifier" );
		if ( modifiers.has_value( modifier_name ) )	idx = modifiers.index( modifier_name );
	}
	if ( idx == 0 ) return;

	utility::vector1< Real > reactivity = rdat.reactivity()[ idx ];
	for ( Size n = 1; n <= reactivity.size(); n++ ){
		rna_data_info_with_conventional_numbering_->add_reactivity( scoring::rna::data::RNA_Reactivity( rdat.seqpos()[n], type, reactivity[n] ) );
	}
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

		if ( tag == "RDAT_VERSION" || tag == "VERSION" ) { // RDAT file!
			data_stream.close();
			read_data_from_rdat( filename );
			return;
		}

		if (line_stream.fail() ) continue; //Probably a blank line.

		if  ( tag == "EXPOSE" ) {

			read_data_info( line_stream );

		} else if ( tag == "BACKBONE_BURIAL" ) {

			read_backbone_info( line_stream, backbone_burial_res_ );

		} else if ( tag == "BACKBONE_EXPOSE" ) {

			read_backbone_info( line_stream, backbone_exposed_res_ );

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
ObjexxFCL::FArray1D< bool >
RNA_DataReader::fill_backbone_array( utility::vector1< Size > const & backbone_res,
																		 core::pose::Pose const & pose ) {

	using namespace core::pose::full_model_info;

	Size const nres = pose.total_residue();
	ObjexxFCL::FArray1D< bool > backbone_array( nres, false );

	FullModelParametersCOP full_model_parameters;
	if ( full_model_info_defined( pose ) ) full_model_parameters = const_full_model_info( pose ).full_model_parameters();

	for ( Size n = 1; n <= backbone_res.size(); n++ ){
		Size pos( backbone_res[ n ] );
		if ( full_model_parameters ) {
			if ( !full_model_parameters->has_conventional_residue( pos ) ) continue; // should we give a warning?
			pos = full_model_parameters->conventional_to_full( pos );
		}
		runtime_assert( pos >= 1 && pos <= nres );
		backbone_array( pos ) = true;
	}

	return backbone_array;
}

/////////////////////////////////////////////////////////////////////
void
RNA_DataReader::fill_rna_data_info( pose::Pose & pose )
{
	using namespace core::pose::full_model_info;

	RNA_DataInfoOP rna_data_info_new = new RNA_DataInfo;

	FullModelParametersCOP full_model_parameters;
	if ( full_model_info_defined( pose ) ) full_model_parameters = nonconst_full_model_info( pose ).full_model_parameters();

	rna_data_info_new->set_backbone_burial (   fill_backbone_array( backbone_burial_res_, pose ) );
	rna_data_info_new->set_backbone_exposed(   fill_backbone_array( backbone_exposed_res_, pose ) );

	runtime_assert( rna_data_info_with_conventional_numbering_ != 0 );
	RNA_Data const & rna_data_with_conventional_numbering = rna_data_info_with_conventional_numbering_->rna_data();
	for ( Size n = 1; n <= rna_data_with_conventional_numbering.size(); n++ ){
		RNA_Datum rna_datum = rna_data_with_conventional_numbering[ n ];
		if ( full_model_parameters ){
			if ( !full_model_parameters->has_conventional_residue( rna_datum.position() ) ) continue; // should we give a warning?
			rna_datum.set_position( full_model_parameters->conventional_to_full( rna_datum.position() ) );
		}
		rna_data_info_new->add_datum( rna_datum );
	}

	RNA_Reactivities const & rna_reactivities_with_conventional_numbering = rna_data_info_with_conventional_numbering_->rna_reactivities();
	for ( Size n = 1; n <= rna_reactivities_with_conventional_numbering.size(); n++ ){
		RNA_Reactivity rna_reactivity = rna_reactivities_with_conventional_numbering[ n ];
		if ( full_model_parameters ) {
			if ( !full_model_parameters->has_conventional_residue( rna_reactivity.position() ) ) continue; // should we give a warning?
			rna_reactivity.set_position( full_model_parameters->conventional_to_full( rna_reactivity.position() ) );
		}
		rna_data_info_new->add_reactivity( rna_reactivity );
	}

	scoring::rna::RNA_ScoringInfo & rna_scoring_info( scoring::rna::nonconst_rna_scoring_info_from_pose( pose ) );
	scoring::rna::data::RNA_DataInfo & rna_data_info( rna_scoring_info.rna_data_info() );
	rna_data_info = *rna_data_info_new;

}


/////////////////////////////////////////////////////////////////////
core::scoring::rna::data::RNA_DataInfo const &
get_rna_data_info( pose::Pose & pose, std::string const & rna_data_file,
									 core::scoring::ScoreFunctionOP scorefxn /* = 0 */) {

	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace core::scoring::rna;
	using namespace core::scoring::rna::data;

	if ( rna_data_file.size() > 0 ) {
		RNA_DataReaderOP rna_data_reader = new RNA_DataReader( rna_data_file );
		rna_data_reader->fill_rna_data_info( pose );
	}

	RNA_DataInfo const & rna_data_info = scoring::rna::nonconst_rna_scoring_info_from_pose( pose ).rna_data_info();

	utility::vector1< PoseOP > const & other_pose_list = const_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ){
		nonconst_rna_scoring_info_from_pose( *(other_pose_list[ n ]) ).rna_data_info() = rna_data_info;
	}

	if ( scorefxn != 0 &&
			 scoring::rna::rna_scoring_info_from_pose( pose ).rna_data_info().rna_reactivities().size() > 0 ) {
		scorefxn->set_weight( core::scoring::rna_chem_map, 1.0 );
	}

	return rna_data_info;
}



} //rna
} //io
} //core
