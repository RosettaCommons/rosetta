// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/hbonds/HBondDatabase.cc
/// @brief  Database containing params for HBondEnergy
/// @author John Karanicolas
/// @author Matthew O'Meara


// Unit Headers
#include <core/scoring/hbonds/HBondDatabase.hh>

// Package Headers
#include <core/scoring/hbonds/types.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/polynomial.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/constants.hh>
#include <core/scoring/hbonds/FadeInterval.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <basic/database/schema_generator/DbDataType.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


// Utility Headers
#include <utility/file/file_sys_util.hh>
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>

// External Headers
#include <cppdb/frontend.h>


// C++ headers
#include <cmath>
#include <sstream>
#include <ObjexxFCL/FArray3D.hh>


namespace core {
namespace scoring {
namespace hbonds {

using std::endl;
using std::map;
using std::pair;
using std::string;
using std::stringstream;
using cppdb::result;
using cppdb::statement;
using basic::database::open;
using basic::database::full_name;
using basic::Tracer;
using utility::file::file_exists;
using utility::io::izstream;
using utility::string_split;
using utility::vector1;
using utility::sql_database::sessionOP;

static Tracer tr("core.scoring.hbonds.HBondDatabase");
// Initialize private static data
map< const string, HBondDatabaseCOP > HBondDatabase::initialized_databases_;


HBondDatabase::HBondDatabase():
	initialized_(false),
	params_database_tag_(""),
	HBFadeInterval_lookup_by_name_(),
	HBFadeInterval_lookup_(),
	AHdist_short_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	AHdist_long_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	HBPoly1D_lookup_by_name_(),
	HBPoly1D_lookup_(),
	AHdist_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_short_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_long_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_short_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_long_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	chi_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	don_strength_lookup_(hbdon_MAX, 1.0),
	acc_strength_lookup_(hbacc_MAX, 1.0),
	weight_type_lookup_(HB_EVAL_TYPE_COUNT, hbw_NONE)
{
	HBondOptions hb_options; // default ctor reads options system, initializes default parameters tag from which this database will initialize itself.
	params_database_tag_ = hb_options.params_database_tag();
	initialize();
}

HBondDatabase::HBondDatabase(
	//HBondOptionsCOP hb_options
	std::string const & hbond_params_database_tag
) :
	initialized_(false),
	//hb_options_( hb_options ),
	params_database_tag_( hbond_params_database_tag ),
	HBFadeInterval_lookup_by_name_(),
	HBFadeInterval_lookup_(),
	AHdist_short_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	AHdist_long_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_fade_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	HBPoly1D_lookup_by_name_(),
	HBPoly1D_lookup_(),
	AHdist_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_short_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosBAH_long_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_short_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	cosAHD_long_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	chi_poly_lookup_(HB_EVAL_TYPE_COUNT, NULL),
	don_strength_lookup_(hbdon_MAX, 1.0),
	acc_strength_lookup_(hbacc_MAX, 1.0),
	weight_type_lookup_(HB_EVAL_TYPE_COUNT, hbw_NONE)
{
	initialize();
}

HBondDatabase::HBondDatabase(
	const HBondDatabase & src
) :
	ReferenceCount( src ),
	initialized_(src.initialized_),
	//hb_options_( src.hb_options_ ),
	params_database_tag_( src.params_database_tag_ ),
	HBFadeInterval_lookup_by_name_(src.HBFadeInterval_lookup_by_name_),
	HBFadeInterval_lookup_(src.HBFadeInterval_lookup_),
	AHdist_short_fade_lookup_(src.AHdist_short_fade_lookup_),
	AHdist_long_fade_lookup_(src.AHdist_long_fade_lookup_),
	cosBAH_fade_lookup_(src.cosBAH_fade_lookup_),
	cosAHD_fade_lookup_(src.cosAHD_fade_lookup_),
	HBPoly1D_lookup_by_name_(src.HBPoly1D_lookup_by_name_),
	HBPoly1D_lookup_(src.HBPoly1D_lookup_),
	AHdist_poly_lookup_(src.AHdist_poly_lookup_),
	cosBAH_short_poly_lookup_(src.cosBAH_short_poly_lookup_),
	cosBAH_long_poly_lookup_(src.cosBAH_long_poly_lookup_),
	cosAHD_short_poly_lookup_(src.cosAHD_short_poly_lookup_),
	cosAHD_long_poly_lookup_(src.cosAHD_long_poly_lookup_),
	chi_poly_lookup_(src.chi_poly_lookup_),
	don_strength_lookup_(hbdon_MAX, 1.0),
	acc_strength_lookup_(hbacc_MAX, 1.0),
	weight_type_lookup_(src.weight_type_lookup_)
{}

HBondDatabaseCOP
HBondDatabase::get_database(){

	HBondOptions hb_options; // default ctor queries options system, initializes default hbond params

	map< string const, HBondDatabaseCOP >::const_iterator
		param_db = initialized_databases_.find( hb_options.params_database_tag() );
	if ( param_db == initialized_databases_.end() ) {
		HBondDatabaseCOP newdb( HBondDatabaseOP( new HBondDatabase( hb_options.params_database_tag() ) ) );
		initialized_databases_[ hb_options.params_database_tag() ] = newdb;
		return newdb;
	}
	return param_db->second;
}

HBondDatabaseCOP
HBondDatabase::get_database( string const & tag ){

	map< string const, HBondDatabaseCOP >::const_iterator
		param_db = initialized_databases_.find(tag);
	if ( param_db == initialized_databases_.end() ) {
		HBondDatabaseCOP newdb( HBondDatabaseOP( new HBondDatabase( tag ) ) );
		initialized_databases_[ tag ] = newdb;
		return newdb;
	}
	return param_db->second;
}


HBondDatabase::~HBondDatabase(){}


/// @details initialize hydrogen bond parameters
void
HBondDatabase::initialize()
{
	if ( initialized_ ) {
		tr << "Re-intializing HBond Database when it has already been initialized!";
	}

	initialize_HBPoly1D();
	initialize_HBFadeInterval();
	initialize_HBEval();

	// Note if these aren't defined, then they are silently ignored
	// Currently their use is experimental and are only defined in newer
	// hbond parameter sets.
	initialize_don_strength();
	initialize_acc_strength();

	initialized_ = true;
}

/// @details has the database already been initialized?
bool
HBondDatabase::initialized() const
{
	return initialized_;
}


void
HBondDatabase::initialize_HBFadeInterval()
{
	string HBFadeInterval_fname = "scoring/score_functions/hbonds/" + params_database_tag_ + "/HBFadeIntervals.csv";

	izstream s;
	if ( !open(s, HBFadeInterval_fname) ) {
		stringstream message;
		message << "Unable to open hbond parameter file HBFadeInterval:" << endl;
		message << "'" << HBFadeInterval_fname << "'";
		utility_exit_with_message(message.str());
	}

	Size line_no(0);
	string line;
	vector1<string> tokens;
	Size id;
	string fade_interval_name;
	bool smoothed(false);
	Real min0, fmin, fmax, max0;
	while ( getline( s, line ) ) {
		++line_no;
		tokens = string_split( line, ',');
		Size ntokens = 8;
		if ( tokens.size() != ntokens ) {
			stringstream message;
			message << "FadeInterval definition line does not have the expected number of fields " << endl;
			message << "Expected '" << ntokens << "' tokens but found '" << tokens.size() << "' tokens. " << endl;
			message << line << endl;
			utility_exit_with_message(message.str());
		}

		Size i(1);
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> id; }
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> fade_interval_name; }
		{
			string junction_type;
			stringstream buf; buf << tokens[i]; i++; buf >> junction_type;
			if ( junction_type == "smoothed" ) {
				smoothed = true;
			} else if ( junction_type == "piecewise_linear" ) {
				smoothed = false;
			} else {
				stringstream message;
				message
					<< "On line " << HBFadeInterval_fname << ":" << line_no << ": '" << line << "'"
					<< "the junction_type should be either 'smoothed' or 'piecewise_linear', "
					<< "however it the unrecognized string '" << junction_type << "'." << endl;
				utility_exit_with_message(message.str());
			}
		}
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> min0; }
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> fmin; }
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> fmax; }
		{ stringstream buf;buf.precision(16); buf << tokens[i]; i++; buf >> max0; }

		FadeIntervalOP fade_interval( new FadeInterval(fade_interval_name, min0, fmin, fmax, max0, smoothed) );

		if ( HBFadeInterval_lookup_.size() + 1 != id ) {
			stringstream message;
			message << "The id fields in the HBFadeInterval file '" << HBFadeInterval_fname << "'" << endl;
			message << "are out of order or missing: Expected id: " << HBFadeInterval_lookup_.size() + 1 << " but instead found id: " << id << endl;
			message << "on line: '" << line << "'" << endl;
			utility_exit_with_message(message.str());
		}
		HBFadeInterval_lookup_.push_back(fade_interval);
		HBFadeInterval_lookup_by_name_[fade_interval_name] = fade_interval;
	}
}

/// @details read one dimensional polynomial definition file
// File Format:
//    -fields are space delimited
//    -Columns are: polynomial_name, geometric_dimension, xmin, xmax, root1, root2, degree, c_a, c_b, ..., c_k
void
HBondDatabase::initialize_HBPoly1D()
{

	string HBPoly1D_fname = "scoring/score_functions/hbonds/" + params_database_tag_ + "/HBPoly1D.csv";

	izstream s;
	if ( !open(s, HBPoly1D_fname) ) {
		stringstream message;
		message << "Unable to open hbond parameter file HBPoly1D:" << endl;
		message << HBPoly1D_fname;
		utility_exit_with_message(message.str());
	}

	string line;
	vector1<string> tokens;
	Size id;
	string polynomial_name;
	string geo_dim_name;
	HBGeoDimType geometric_dimension;
	Real xmin, xmax, min_val, max_val, root1, root2;
	Size degree;
	vector1< Real > coefficients_;
	while ( getline( s, line ) ) {
		tokens = string_split( line, ',');
		Size ntokens = 22;
		if ( tokens.size() != ntokens ) {
			stringstream message;
			message << "Polynomial definition line does not have enough fields" << endl;
			message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
			message << line << endl;
			utility_exit_with_message(message.str());
		}

		Size i(1);
		{ stringstream buf; buf << tokens[i]; i++; buf >> id;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> polynomial_name;}
		i++; // classic name field
		{ stringstream buf; buf << tokens[i]; i++; buf >> geo_dim_name;
			geometric_dimension = HBondTypeManager::geo_dim_type_from_name( geo_dim_name ); }
		{ stringstream buf; buf << tokens[i]; i++; buf >> xmin;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> xmax;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> min_val;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> max_val;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> root1;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> root2;}
		{ stringstream buf; buf << tokens[i]; i++; buf >> degree;}

		vector1< Real > coefficients_;
		Real c;

		while ( i <= tokens.size() && i <= 11 + degree ) {
			stringstream buf; buf << tokens[i]; i++; buf >> c;
			coefficients_.push_back(c);
		}

		Polynomial_1dOP p;
		try{
			p = Polynomial_1dOP( new Polynomial_1d(
				polynomial_name,
				geometric_dimension,
				xmin, xmax, min_val, max_val, root1, root2,
				degree,
				coefficients_) );
		} catch ( utility::excn::EXCN_Msg_Exception& excn ) {
			std::stringstream msg;
			msg
				<< "ERROR: Unable to construct the polynomial '" << polynomial_name << "' "
				<< "with the following parameters line:" << std::endl
				<< std::endl
				<< "    " << line << std::endl
				<< std::endl
				<< "because of the following errors:" << std::endl
				<< excn.msg() << std::endl;
			throw( utility::excn::EXCN_Msg_Exception( msg.str() ) );
		}

		if ( HBPoly1D_lookup_.size() + 1 != id ) {
			stringstream message;
			message << "The id fields in the HBPoly1D file '" << HBPoly1D_fname << "'" << endl;
			message << "are out of order or missing: Expected id: " << HBPoly1D_lookup_.size() + 1 << " but instead found id: " << id << endl;
			message << "on line: '" << line << "'" << endl;
			utility_exit_with_message(message.str());
		}
		HBPoly1D_lookup_.push_back(p);
		HBPoly1D_lookup_by_name_[polynomial_name] = p;

	}
}

/// @details read one dimensional polynomial definition file
// File Format:
//    -fields are space delimited
//    -Columns are: HBDonChemType, HBAccChemType, HBSeqSep, AHdist_short_fade_name, AHdist_long_fade_name, cosBAH_fade_name, cosAHD_fade_name, AHDist_poly_name, cosBAH_poly_name, cosAHD_poly_name,

void
HBondDatabase::initialize_HBEval()
{

	string HBEval_fname = "scoring/score_functions/hbonds/" + params_database_tag_ + "/HBEval.csv";

	izstream s;
	open(s, HBEval_fname);
	Size line_number(0);
	vector1<string> tokens;
	string line;
	vector1<bool> initialized_hbe_types(hbe_MAX,false);
	string AHdist_poly_name, cosBAH_short_poly_name, cosBAH_long_poly_name,
		cosAHD_short_poly_name, cosAHD_long_poly_name, chi_poly_name,
		don_chem_type_name, acc_chem_type_name, seq_sep_type_name,
		AHdist_short_fade_name, AHdist_long_fade_name,
		cosBAH_fade_name, cosAHD_fade_name;
	HBDonChemType don_chem_type;
	HBAccChemType acc_chem_type;
	HBSeqSep seq_sep_type;
	Polynomial_1dCOP AHdist_poly, cosBAH_short_poly, cosBAH_long_poly,
		cosAHD_short_poly, cosAHD_long_poly, chi_poly;
	FadeIntervalCOP AHdist_short_fade, AHdist_long_fade, cosBAH_fade, cosAHD_fade;
	string weight_type_name;
	HBondWeightType weight_type;

	initialize_HBEval_lookup();

	while ( getline(s, line) ) {
		++line_number;
		tokens = string_split( line, ',');
		if ( tokens.size() != 15 ) {
			stringstream message;
			message << "HBond evaluation data line does not have enough fields" << endl;
			message << "Expected '" << 15 << "' tokens but found '" << tokens.size() << "' tokens." << endl;
			message << line << endl;
			utility_exit_with_message(message.str());
		}

		Size i(1);
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> don_chem_type_name;
			don_chem_type = HBDonChemType(HBondTypeManager::don_chem_type_from_name(don_chem_type_name));
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> acc_chem_type_name;
			acc_chem_type = HBAccChemType(HBondTypeManager::acc_chem_type_from_name(acc_chem_type_name));
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> seq_sep_type_name;
			seq_sep_type = HBSeqSep(HBondTypeManager::seq_sep_type_from_name(seq_sep_type_name));
		}

		HBEvalType hbe_type = (*HBEval_lookup)(don_chem_type, acc_chem_type, seq_sep_type);
		if ( initialized_hbe_types[hbe_type] ) {
			tr << "Duplicate parameter specification in HBEval.csv:" << endl;
			tr << "  hbe_type: " << hbe_type << " line :" << line_number << endl;
		} else {
			initialized_hbe_types[hbe_type] = true;
		}

		if ( hbe_type > static_cast<HBEvalType>(HB_EVAL_TYPE_COUNT) ) {
			stringstream message;
			message << "hb_eval_type created from" << endl;
			message << "\tdon_chem_type:'"<< don_chem_type_name <<"'" << endl;
			message << "\tacc_chem_type:'"<< acc_chem_type_name <<"'" << endl;
			message << "\tseq_sep_type: '"<< seq_sep_type_name << "'" << endl;
			message << "gives type" << hbe_type << ", which is out of the range valid range (0," << HB_EVAL_TYPE_COUNT << ")" << endl;
			utility_exit_with_message(message.str());
		}

		// NOTE: The rows in the HBEval table contain all combinations of
		// donor chemical types, acceptor chemical types and sequence
		// separation types, however the HBEvalType classification
		// collapses these types for generic types. For example
		// (hbdon_NONE, hbacc_NONE, seq_sep_other) and (hbdon_PBA,
		// hbacc_NONE, seq_sep_other) both map to hbe_NONE. To make sure
		// this collapsing does not cause problems, assert that the
		// fade and plynomial functions for each hbe_type are not assigned
		// two different values.

		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> AHdist_short_fade_name;
			AHdist_short_fade = HBFadeInterval_from_name(AHdist_short_fade_name);
			if ( AHdist_short_fade_lookup_[hbe_type] ) {
				debug_assert(AHdist_short_fade_lookup_[hbe_type] == AHdist_short_fade);
			} else {
				AHdist_short_fade_lookup_[hbe_type] = AHdist_short_fade;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> AHdist_long_fade_name;
			AHdist_long_fade = HBFadeInterval_from_name(AHdist_long_fade_name);
			if ( AHdist_long_fade_lookup_[hbe_type] ) {
				debug_assert(AHdist_long_fade_lookup_[hbe_type] == AHdist_long_fade);
			} else {
				AHdist_long_fade_lookup_[hbe_type] = AHdist_long_fade;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosBAH_fade_name;
			cosBAH_fade = HBFadeInterval_from_name(cosBAH_fade_name);
			if ( cosBAH_fade_lookup_[hbe_type] ) {
				debug_assert(cosBAH_fade_lookup_[hbe_type] == cosBAH_fade);
			} else {
				cosBAH_fade_lookup_[hbe_type] = cosBAH_fade;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosAHD_fade_name;
			cosAHD_fade = HBFadeInterval_from_name(cosAHD_fade_name);
			if ( cosAHD_fade_lookup_[hbe_type] ) {
				debug_assert(cosAHD_fade_lookup_[hbe_type] == cosAHD_fade);
			} else {
				cosAHD_fade_lookup_[hbe_type] = cosAHD_fade;
			}
		}
		++i; // fade for chi dimension
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> AHdist_poly_name;
			AHdist_poly = HBPoly1D_from_name(AHdist_poly_name);
			/// Error handling: make sure we're given a distance polynomial
			if ( AHdist_poly->geometric_dimension() != hbgd_AHdist ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", expected to read a distance polynomial (i.e. geometric_dimension == hbgd_AHdist), but instead, found " + AHdist_poly_name + " of geometric dimension " + utility::to_string(AHdist_poly->geometric_dimension()) );
			}
			if ( AHdist_poly_lookup_[hbe_type] ) {
				debug_assert(AHdist_poly_lookup_[hbe_type] == AHdist_poly);
			} else {
				AHdist_poly_lookup_[hbe_type] = AHdist_poly;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosBAH_short_poly_name;
			cosBAH_short_poly = HBPoly1D_from_name(cosBAH_short_poly_name);
			/// Error handling: make sure we're given a cosBAH polynomial
			if ( cosBAH_short_poly->geometric_dimension() != hbgd_cosBAH ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", expected to read a short-range cosBAH polynomial (i.e. geometric_dimension == hbgd_cosBAH), but instead, found " + cosBAH_short_poly_name + " of geometric dimension " + utility::to_string(cosBAH_short_poly->geometric_dimension()) );
			}
			if ( cosBAH_short_poly_lookup_[hbe_type] ) {
				debug_assert(cosBAH_short_poly_lookup_[hbe_type] == cosBAH_short_poly);
			} else {
				cosBAH_short_poly_lookup_[hbe_type] = cosBAH_short_poly;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosBAH_long_poly_name;
			cosBAH_long_poly = HBPoly1D_from_name(cosBAH_long_poly_name);
			/// Error handling: make sure we're given a cosBAH polynomial
			if ( cosBAH_long_poly->geometric_dimension() != hbgd_cosBAH ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", expected to read a long-range cosBAH polynomial (i.e. geometric_dimension == hbgd_cosBAH), but instead, found " + cosBAH_long_poly_name + " of geometric dimension " + utility::to_string(cosBAH_long_poly->geometric_dimension()) );
			}
			if ( cosBAH_long_poly_lookup_[hbe_type] ) {
				debug_assert(cosBAH_long_poly_lookup_[hbe_type] == cosBAH_long_poly);
			} else {
				cosBAH_long_poly_lookup_[hbe_type] = cosBAH_long_poly;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosAHD_short_poly_name;
			cosAHD_short_poly = HBPoly1D_from_name(cosAHD_short_poly_name);
			/// Error handling: make sure we're given either a cosAHD polynomial or an AHD polynomial
			if ( cosAHD_short_poly->geometric_dimension() != hbgd_cosAHD && cosAHD_short_poly->geometric_dimension() != hbgd_AHD  ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", expected to read a short-range cosAHD or AHD polynomial (i.e. geometric_dimension == hbgd_cosAHD or hbgd_AHD), but instead, found " + cosAHD_short_poly_name + " of geometric dimension " + utility::to_string(cosAHD_short_poly->geometric_dimension()) );
			}
			if ( cosAHD_short_poly_lookup_[hbe_type] ) {
				debug_assert(cosAHD_short_poly_lookup_[hbe_type]);
			} else {
				cosAHD_short_poly_lookup_[hbe_type] = cosAHD_short_poly;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> cosAHD_long_poly_name;
			cosAHD_long_poly = HBPoly1D_from_name(cosAHD_long_poly_name);
			/// Error handling: make sure we're given either a cosAHD polynomial or an AHD polynomial
			if ( cosAHD_long_poly->geometric_dimension() != hbgd_cosAHD && cosAHD_long_poly->geometric_dimension() != hbgd_AHD  ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", expected to read a long-range cosAHD or AHD polynomial (i.e. geometric_dimension == hbgd_cosAHD or hbgd_AHD), but instead, found " + cosAHD_long_poly_name + " of geometric dimension " + utility::to_string(cosAHD_long_poly->geometric_dimension()) );
			}
			/// Error handling: make sure the type of the long and short polynomials are consistent -- no mixing cosAHD and AHD polynomials!
			if ( cosAHD_long_poly->geometric_dimension() != cosAHD_short_poly->geometric_dimension()  ) {
				utility_exit_with_message("When reading HBEval.csv parameters for " + don_chem_type_name + " and " + acc_chem_type_name + ", found that the short- and long-range polynomials for the AHD angle are of different geometric types: one of hbgd_cosAHD and the other of hbgd_AHD.  These types cannot be mixed" );
			}
			if ( cosAHD_long_poly_lookup_[hbe_type] ) {
				debug_assert(cosAHD_long_poly_lookup_[hbe_type] == cosAHD_long_poly);
			} else {
				cosAHD_long_poly_lookup_[hbe_type] = cosAHD_long_poly;
			}
		}
		{
			stringstream buf;
			buf << tokens[i]; i++;
			buf >> weight_type_name;
			weight_type = HBondWeightType(HBondTypeManager::weight_type_from_name(weight_type_name));
			if ( weight_type_lookup_[hbe_type] != hbw_NONE ) {
				debug_assert(weight_type_lookup_[hbe_type] == weight_type);
			} else {
				weight_type_lookup_[hbe_type] = weight_type;
			}
		}
	}
	for ( Size i=1; i <= hbe_MAX; ++i ) {
		if ( !initialized_hbe_types[i] ) {
			tr << "hbe_type: " << i << " is not initialized in HBEval.csv" << endl;
		}
	}
}

void
HBondDatabase::initialize_don_strength() {
	string don_strength_fname =
		"scoring/score_functions/hbonds/" + params_database_tag_ + "/DonStrength.csv";
	izstream s;

	bool open_failed( true );
	if ( file_exists( full_name( don_strength_fname, false /*warn*/ ) ) ) {
		open_failed = false;
		open(s, don_strength_fname, false);
	}

	if ( !open_failed ) {
		Size line_number(0);
		Size expected_n_tokens(2);
		vector1<string> tokens;
		string line;
		string don_type_name;
		Real site_strength;
		HBDonChemType don_chem_type;
		while ( getline(s, line) ) {
			++line_number;
			tokens = string_split(line, ',');
			if ( tokens.size() != expected_n_tokens ) {
				stringstream message;
				message
					<< "BondStrength.csv:" << line_number << " "
					<< "should have " << expected_n_tokens << ", "
					<< "however it has " << tokens.size() << endl;
				message
					<< line << endl;
				utility_exit_with_message(message.str());
			}

			Size i(1);
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> don_type_name;
			}
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> site_strength;
			}

			don_chem_type = HBondTypeManager::don_chem_type_from_name(don_type_name);
			don_strength_lookup_[don_chem_type] = site_strength;
		}
	} // db file present?
	/// now adjust from cmdline
	if ( basic::options::option[ basic::options::OptionKeys::score::hb_don_strength ].user() ) {
		utility::vector1< std::string > const hb_don_strength_tags
			( basic::options::option[ basic::options::OptionKeys::score::hb_don_strength ] );
		for ( std::string const & tag : hb_don_strength_tags ) {
			utility::vector1< std::string > const tokens( string_split( tag, ':' ) );
			runtime_assert( tokens.size() == 2 );
			string don_type_name("tmp");
			Real site_strength(0);
			Size i(1);
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> don_type_name;
			}
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> site_strength;
			}

			HBDonChemType don_chem_type( HBondTypeManager::don_chem_type_from_name(don_type_name) );
			don_strength_lookup_[don_chem_type] = site_strength;
			tr.Trace << "initialize_don_strength from cmdline " << don_type_name << ' ' << don_chem_type << ' ' <<
				site_strength << std::endl;
		}
	}
}

void
HBondDatabase::initialize_acc_strength() {
	string acc_strength_fname =
		"scoring/score_functions/hbonds/" + params_database_tag_ + "/AccStrength.csv";
	izstream s;

	bool open_failed( true );
	if ( file_exists( full_name( acc_strength_fname, false /*warn*/ ) ) ) {
		open_failed = false;
		open(s, acc_strength_fname, false);
	}

	if ( !open_failed ) {
		Size line_number(0);
		Size expected_n_tokens(2);
		vector1<string> tokens;
		string line;
		string acc_type_name;
		Real site_strength;
		HBAccChemType acc_chem_type;
		while ( getline(s, line) ) {
			++line_number;
			tokens = string_split(line, ',');
			if ( tokens.size() != expected_n_tokens ) {
				stringstream message;
				message
					<< "BondStrength.csv:" << line_number << " "
					<< "should have " << expected_n_tokens << ", "
					<< "however it has " << tokens.size() << endl;
				message
					<< line << endl;
				utility_exit_with_message(message.str());
			}

			Size i(1);
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> acc_type_name;
			}
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> site_strength;
			}

			acc_chem_type = HBondTypeManager::acc_chem_type_from_name(acc_type_name);
			acc_strength_lookup_[acc_chem_type] = site_strength;
		}
	} // read from db file if present

	/// now adjust from cmdline
	if ( basic::options::option[ basic::options::OptionKeys::score::hb_acc_strength ].user() ) {
		utility::vector1< std::string > const hb_acc_strength_tags
			( basic::options::option[ basic::options::OptionKeys::score::hb_acc_strength ] );
		for ( std::string const & tag : hb_acc_strength_tags ) {
			utility::vector1< std::string > const tokens( string_split( tag, ':' ) );
			runtime_assert( tokens.size() == 2 );
			string acc_type_name("tmp");
			Real site_strength(0);
			Size i(1);
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> acc_type_name;
			}
			{
				stringstream buf;
				buf << tokens[i]; i++;
				buf >> site_strength;
			}

			HBAccChemType acc_chem_type( HBondTypeManager::acc_chem_type_from_name(acc_type_name) );
			acc_strength_lookup_[acc_chem_type] = site_strength;
			tr.Trace << "initialize_acc_strength from cmdline " << acc_type_name << ' ' << acc_chem_type << ' ' <<
				site_strength << std::endl;
		}
	}
}

FadeIntervalCOP
HBondDatabase::HBFadeInterval_from_name(
	string const & name
) const {
	map< const string, FadeIntervalCOP >::const_iterator it(HBFadeInterval_lookup_by_name_.find(name));
	if ( it == HBFadeInterval_lookup_by_name_.end() ) {
		stringstream message;
		message << "Fade Interval '" << name << "' has not been defined.";
		utility_exit_with_message(message.str());
		return NULL;
	} else {
		return it->second;
	}
}


/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
FadeIntervalCOP
HBondDatabase::AHdist_short_fade_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type <<
			"' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}

	FadeIntervalCOP p(AHdist_short_fade_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No short fade interval for AHdist has been defined for hb eval type '"
			<< hb_eval_type << "'" <<endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
FadeIntervalCOP
HBondDatabase::AHdist_long_fade_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	FadeIntervalCOP p(AHdist_long_fade_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No long fade interval for AHdist has been defined for hb eval type '"
			<< hb_eval_type << "'" <<endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
FadeIntervalCOP
HBondDatabase::cosBAH_fade_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	FadeIntervalCOP p(cosBAH_fade_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No fade interval for cosBAH has been defined for hb eval type '"
			<< hb_eval_type << "'" <<endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
FadeIntervalCOP
HBondDatabase::cosAHD_fade_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	FadeIntervalCOP p(cosAHD_fade_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No fade interval for cosAHD has been defined for hb eval type '"
			<< hb_eval_type << "'" <<endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

Polynomial_1dCOP
HBondDatabase::HBPoly1D_from_name(
	string const & name
) const {
	map< const string, Polynomial_1dCOP >::const_iterator it(HBPoly1D_lookup_by_name_.find(name));
	if ( it == HBPoly1D_lookup_by_name_.end() ) {
		stringstream message;
		message << "1d Polynomial '" << name << "' has not been defined.";
		utility_exit_with_message(message.str());
		return NULL;
	} else {
		return it->second;
	}
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::AHdist_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(AHdist_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No AHdist polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" <<endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::cosBAH_short_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(cosBAH_short_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No cosBAH_short polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;

}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::cosBAH_long_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(cosBAH_long_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No cosBAH_long polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::cosAHD_short_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(cosAHD_short_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No cosAHD_short polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::cosAHD_long_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(cosAHD_long_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No cosAHD_long polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;
}
/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
Polynomial_1dCOP
HBondDatabase::chi_poly_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	Polynomial_1dCOP p(chi_poly_lookup_[hb_eval_type]);

	if ( !p ) {
		stringstream message;
		message << "No chi polynomial has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;
}


Real
HBondDatabase::don_strength(
	HBDonChemType const don_chem_type
) const {
	debug_assert(don_chem_type >= 1 && don_chem_type <= hbdon_MAX );

	return don_strength_lookup_[don_chem_type];
}

Real
HBondDatabase::acc_strength(
	HBAccChemType const acc_chem_type
) const {
	debug_assert(acc_chem_type >= 1 && acc_chem_type <= hbacc_MAX );

	return acc_strength_lookup_[acc_chem_type];
}

/// @details use get_hbond_evaluation_type(...) or HBEval_lookup(...)
///determine hb_eval_type.
HBondWeightType
HBondDatabase::weight_type_lookup(
	Size const hb_eval_type
) const {
	if ( hb_eval_type < 1 || hb_eval_type > HB_EVAL_TYPE_COUNT ) {
		stringstream message;
		message << "HBond eval type '" << hb_eval_type
			<< "' is out side of the valid range (1," << HB_EVAL_TYPE_COUNT << ")";
		utility_exit_with_message(message.str());
	}
	HBondWeightType p(weight_type_lookup_[hb_eval_type]);

	if ( p == hbw_NONE ) {
		stringstream message;
		message << "No weight type has been defined for hb eval type '"
			<< hb_eval_type << "'" << endl;
		utility_exit_with_message(message.str());
	}
	return p;
}

void
HBondDatabase::report_parameter_features_schema_to_db(
	sessionOP db_session
) const {
	write_hbond_fade_interval_table_schema(db_session);
	write_hbond_polynomial_1d_table_schema(db_session);
	write_hbond_evaluation_types_table_schema(db_session);
}

void
HBondDatabase::write_hbond_fade_interval_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column database_tag("database_tag", DbDataTypeOP( new DbText() ));
	Column name("name", DbDataTypeOP( new DbText() ));
	Column junction_type("junction_type", DbDataTypeOP( new DbText() ));
	Column min0("min0", DbDataTypeOP( new DbReal() ));
	Column fmin("fmin", DbDataTypeOP( new DbReal() ));
	Column fmax("fmax", DbDataTypeOP( new DbReal() ));
	Column max0("max0", DbDataTypeOP( new DbReal() ));


	Columns primary_key_columns;
	primary_key_columns.push_back(database_tag);
	primary_key_columns.push_back(name);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("hbond_fade_interval", primary_key);
	table.add_column(junction_type);
	table.add_column(min0);
	table.add_column(fmin);
	table.add_column(fmax);
	table.add_column(max0);

	table.write(db_session);
}

void
HBondDatabase::write_hbond_polynomial_1d_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column database_tag("database_tag", DbDataTypeOP( new DbText() ));
	Column name("name", DbDataTypeOP( new DbText() ));
	Column dimension("dimension", DbDataTypeOP( new DbText() ));
	Column xmin("xmin", DbDataTypeOP( new DbReal() ));
	Column xmax("xmax", DbDataTypeOP( new DbReal() ));
	Column min_val("min_val", DbDataTypeOP( new DbReal() ));
	Column max_val("max_val", DbDataTypeOP( new DbReal() ));
	Column root1("root1", DbDataTypeOP( new DbReal() ));
	Column root2("root2", DbDataTypeOP( new DbReal() ));
	Column degree("degree", DbDataTypeOP( new DbInteger() ));
	Column c_a("c_a", DbDataTypeOP( new DbReal() ));
	Column c_b("c_b", DbDataTypeOP( new DbReal() ));
	Column c_c("c_c", DbDataTypeOP( new DbReal() ));
	Column c_d("c_d", DbDataTypeOP( new DbReal() ));
	Column c_e("c_e", DbDataTypeOP( new DbReal() ));
	Column c_f("c_f", DbDataTypeOP( new DbReal() ));
	Column c_g("c_g", DbDataTypeOP( new DbReal() ));
	Column c_h("c_h", DbDataTypeOP( new DbReal() ));
	Column c_i("c_i", DbDataTypeOP( new DbReal() ));
	Column c_j("c_j", DbDataTypeOP( new DbReal() ));
	Column c_k("c_k", DbDataTypeOP( new DbReal() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(database_tag);
	primary_key_columns.push_back(name);
	PrimaryKey primary_key(primary_key_columns);

	Schema table("hbond_polynomial_1d", primary_key);
	table.add_column(dimension);
	table.add_column(xmin);
	table.add_column(xmax);
	table.add_column(min_val);
	table.add_column(max_val);
	table.add_column(root1);
	table.add_column(root2);
	table.add_column(degree);
	table.add_column(c_a);
	table.add_column(c_b);
	table.add_column(c_c);
	table.add_column(c_d);
	table.add_column(c_e);
	table.add_column(c_f);
	table.add_column(c_g);
	table.add_column(c_h);
	table.add_column(c_i);
	table.add_column(c_j);
	table.add_column(c_k);
	table.write(db_session);
}

void
HBondDatabase::write_hbond_evaluation_types_table_schema(
	sessionOP db_session
) const {
	using namespace basic::database::schema_generator;

	Column database_tag("database_tag", DbDataTypeOP( new DbText() ));
	Column don_chem_type("don_chem_type", DbDataTypeOP( new DbText() ));
	Column acc_chem_type("acc_chem_type", DbDataTypeOP( new DbText() ));
	Column separation("separation", DbDataTypeOP( new DbText() ));
	Column AHdist_short_fade("AHdist_short_fade", DbDataTypeOP( new DbText() ));
	Column AHdist_long_fade("AHdist_long_fade", DbDataTypeOP( new DbText() ));
	Column cosBAH_fade("cosBAH_fade", DbDataTypeOP( new DbText() ));
	Column cosAHD_fade("cosAHD_fade", DbDataTypeOP( new DbText() ));
	Column AHdist("AHdist", DbDataTypeOP( new DbText() ));
	Column cosBAH_short("cosBAH_short", DbDataTypeOP( new DbText() ));
	Column cosBAH_long("cosBAH_long", DbDataTypeOP( new DbText() ));
	Column cosAHD_short("cosAHD_short", DbDataTypeOP( new DbText() ));
	Column cosAHD_long("cosAHD_long", DbDataTypeOP( new DbText() ));
	Column weight_type("weight_type", DbDataTypeOP( new DbText() ));

	Columns primary_key_columns;
	primary_key_columns.push_back(database_tag);
	primary_key_columns.push_back(don_chem_type);
	primary_key_columns.push_back(acc_chem_type);
	primary_key_columns.push_back(separation);
	PrimaryKey primary_key(primary_key_columns);

	Columns foreign_key_columns_AHdist_short_fade;
	foreign_key_columns_AHdist_short_fade.push_back(database_tag);
	foreign_key_columns_AHdist_short_fade.push_back(AHdist_short_fade);
	vector1< std::string > reference_columns_AHdist_short_fade;
	reference_columns_AHdist_short_fade.push_back("database_tag");
	reference_columns_AHdist_short_fade.push_back("name");
	ForeignKey foreign_key_AHdist_short_fade(
		foreign_key_columns_AHdist_short_fade,
		"hbond_fade_interval",
		reference_columns_AHdist_short_fade,
		true);

	Columns foreign_key_columns_AHdist_long_fade;
	foreign_key_columns_AHdist_long_fade.push_back(database_tag);
	foreign_key_columns_AHdist_long_fade.push_back(AHdist_long_fade);
	vector1< std::string > reference_columns_AHdist_long_fade;
	reference_columns_AHdist_long_fade.push_back("database_tag");
	reference_columns_AHdist_long_fade.push_back("name");
	ForeignKey foreign_key_AHdist_long_fade(
		foreign_key_columns_AHdist_long_fade,
		"hbond_fade_interval",
		reference_columns_AHdist_long_fade,
		true);

	Columns foreign_key_columns_cosBAH_fade;
	foreign_key_columns_cosBAH_fade.push_back(database_tag);
	foreign_key_columns_cosBAH_fade.push_back(cosBAH_fade);
	vector1< std::string > reference_columns_cosBAH_fade;
	reference_columns_cosBAH_fade.push_back("database_tag");
	reference_columns_cosBAH_fade.push_back("name");
	ForeignKey foreign_key_cosBAH_fade(
		foreign_key_columns_cosBAH_fade,
		"hbond_fade_interval",
		reference_columns_cosBAH_fade,
		true);

	Columns foreign_key_columns_cosAHD_fade;
	foreign_key_columns_cosAHD_fade.push_back(database_tag);
	foreign_key_columns_cosAHD_fade.push_back(cosAHD_fade);
	vector1< std::string > reference_columns_cosAHD_fade;
	reference_columns_cosAHD_fade.push_back("database_tag");
	reference_columns_cosAHD_fade.push_back("name");
	ForeignKey foreign_key_cosAHD_fade(
		foreign_key_columns_cosAHD_fade,
		"hbond_fade_interval",
		reference_columns_cosAHD_fade,
		true);

	Columns foreign_key_columns_AHdist;
	foreign_key_columns_AHdist.push_back(database_tag);
	foreign_key_columns_AHdist.push_back(AHdist);
	vector1< std::string > reference_columns_AHdist;
	reference_columns_AHdist.push_back("database_tag");
	reference_columns_AHdist.push_back("name");
	ForeignKey foreign_key_AHdist(
		foreign_key_columns_AHdist,
		"hbond_polynomial_1d",
		reference_columns_AHdist,
		true);

	Columns foreign_key_columns_cosBAH_short;
	foreign_key_columns_cosBAH_short.push_back(database_tag);
	foreign_key_columns_cosBAH_short.push_back(cosBAH_short);
	vector1< std::string > reference_columns_cosBAH_short;
	reference_columns_cosBAH_short.push_back("database_tag");
	reference_columns_cosBAH_short.push_back("name");
	ForeignKey foreign_key_cosBAH_short(
		foreign_key_columns_cosBAH_short,
		"hbond_polynomial_1d",
		reference_columns_cosBAH_short,
		true);

	Columns foreign_key_columns_cosBAH_long;
	foreign_key_columns_cosBAH_long.push_back(database_tag);
	foreign_key_columns_cosBAH_long.push_back(cosBAH_long);
	vector1< std::string > reference_columns_cosBAH_long;
	reference_columns_cosBAH_long.push_back("database_tag");
	reference_columns_cosBAH_long.push_back("name");
	ForeignKey foreign_key_cosBAH_long(
		foreign_key_columns_cosBAH_long,
		"hbond_polynomial_1d",
		reference_columns_cosBAH_long,
		true);

	Columns foreign_key_columns_cosAHD_short;
	foreign_key_columns_cosAHD_short.push_back(database_tag);
	foreign_key_columns_cosAHD_short.push_back(cosAHD_short);
	vector1< std::string > reference_columns_cosAHD_short;
	reference_columns_cosAHD_short.push_back("database_tag");
	reference_columns_cosAHD_short.push_back("name");
	ForeignKey foreign_key_cosAHD_short(
		foreign_key_columns_cosAHD_short,
		"hbond_polynomial_1d",
		reference_columns_cosAHD_short,
		true);

	Columns foreign_key_columns_cosAHD_long;
	foreign_key_columns_cosAHD_long.push_back(database_tag);
	foreign_key_columns_cosAHD_long.push_back(cosAHD_long);
	vector1< std::string > reference_columns_cosAHD_long;
	reference_columns_cosAHD_long.push_back("database_tag");
	reference_columns_cosAHD_long.push_back("name");
	ForeignKey foreign_key_cosAHD_long(
		foreign_key_columns_cosAHD_long,
		"hbond_polynomial_1d",
		reference_columns_cosAHD_long,
		true);

	Schema table("hbond_evaluation_types", primary_key);
	table.add_foreign_key(foreign_key_AHdist_short_fade);
	table.add_foreign_key(foreign_key_AHdist_long_fade);
	table.add_foreign_key(foreign_key_cosBAH_fade);
	table.add_foreign_key(foreign_key_cosAHD_fade);
	table.add_foreign_key(foreign_key_AHdist);
	table.add_foreign_key(foreign_key_cosBAH_short);
	table.add_foreign_key(foreign_key_cosBAH_long);
	table.add_foreign_key(foreign_key_cosAHD_short);
	table.add_foreign_key(foreign_key_cosAHD_long);
	table.add_column(AHdist_short_fade);
	table.add_column(AHdist_long_fade);
	table.add_column(cosBAH_fade);
	table.add_column(cosAHD_fade);
	table.add_column(AHdist);
	table.add_column(cosBAH_short);
	table.add_column(cosBAH_long);
	table.add_column(cosAHD_short);
	table.add_column(cosAHD_long);
	table.add_column(weight_type);

	table.write(db_session);
}

Size
HBondDatabase::report_parameter_features(
	sessionOP db_session
) const {

	string const & database_tag(params_database_tag_ );

	std::string select_string = "SELECT * FROM hbond_evaluation_types WHERE database_tag = ?;";
	statement select_statement(basic::database::safely_prepare_statement(select_string,db_session));
	select_statement.bind(1,database_tag);
	result res(basic::database::safely_read_from_database(select_statement));
	if ( res.next() ) return 0;

	//pair<string, FadeIntervalCOP> fade_name_interval;
	std::string hbond_interval_string = "INSERT INTO hbond_fade_interval (database_tag, name, junction_type, min0, fmin, fmax, max0) VALUES (?,?,?,?,?,?,?);";
	statement hbond_interval_statement(basic::database::safely_prepare_statement(hbond_interval_string,db_session));
	for ( auto const & fade_name_interval : HBFadeInterval_lookup_by_name_ ) {
		hbond_interval_statement.bind(1,database_tag);
		hbond_interval_statement.bind(2,fade_name_interval.first);
		hbond_interval_statement.bind(3,(fade_name_interval.second->get_smooth() ? "smooth" : "piecewise_linear"));
		hbond_interval_statement.bind(4,fade_name_interval.second->get_min0());
		hbond_interval_statement.bind(5,fade_name_interval.second->get_fmin());
		hbond_interval_statement.bind(6,fade_name_interval.second->get_fmax());
		hbond_interval_statement.bind(7,fade_name_interval.second->get_max0());
		basic::database::safely_write_to_database(hbond_interval_statement);

	}

	//pair<string, Polynomial_1dCOP> poly_name_fn;
	std::string hbond_polynomial_string = "INSERT INTO hbond_polynomial_1d (database_tag, name, dimension, xmin, xmax, min_val, max_val, root1, root2, degree, c_a, c_b, c_c, c_d, c_e, c_f, c_g, c_h, c_i, c_j, c_k) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement hbond_polynomial_statement(basic::database::safely_prepare_statement(hbond_polynomial_string,db_session));
	for ( auto const & poly_name_fn : HBPoly1D_lookup_by_name_ ) {
		hbond_polynomial_statement.bind(1,database_tag);
		hbond_polynomial_statement.bind(2,poly_name_fn.first);
		hbond_polynomial_statement.bind(3,HBondTypeManager::name_from_geo_dim_type(poly_name_fn.second->geometric_dimension()));
		hbond_polynomial_statement.bind(4,poly_name_fn.second->xmin());
		hbond_polynomial_statement.bind(5,poly_name_fn.second->xmax());
		hbond_polynomial_statement.bind(6,poly_name_fn.second->min_val());
		hbond_polynomial_statement.bind(7,poly_name_fn.second->max_val());
		hbond_polynomial_statement.bind(8,poly_name_fn.second->root1());
		hbond_polynomial_statement.bind(9,poly_name_fn.second->root2());
		hbond_polynomial_statement.bind(10,poly_name_fn.second->degree());
		Size index = 10;
		for ( Size i = 1; i <= poly_name_fn.second->degree(); ++i ) {
			index++;
			hbond_polynomial_statement.bind(index,poly_name_fn.second->coefficients()[i]);
		}
		for ( Size i = 1; i <= 11-poly_name_fn.second->degree(); ++i ) {
			index++;
			hbond_polynomial_statement.bind_null(index);
		}
		basic::database::safely_write_to_database(hbond_polynomial_statement);
	}

	std::string hbond_evaluation_string = "INSERT INTO hbond_evaluation_types (database_tag, don_chem_type, acc_chem_type, separation, AHdist_short_fade, AHdist_long_fade, cosBAH_fade, cosAHD_fade, AHdist, cosBAH_short, cosBAH_long, cosAHD_short, cosAHD_long, weight_type) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?);";
	statement hbond_evaluation_statement(basic::database::safely_prepare_statement(hbond_evaluation_string,db_session));
	for ( Size hbdon=1; hbdon <= hbdon_MAX; ++hbdon ) {
		string const & don_chem_type(HBondTypeManager::name_from_don_chem_type(HBDonChemType(hbdon)));
		for ( Size hbacc=1; hbacc <= hbacc_MAX; ++hbacc ) {
			string const & acc_chem_type(HBondTypeManager::name_from_acc_chem_type(HBAccChemType(hbacc)));
			for ( Size hbseq_sep=1; hbseq_sep <= seq_sep_MAX; ++hbseq_sep ) {
				string const & separation(HBondTypeManager::name_from_seq_sep_type(HBSeqSep(hbseq_sep)));

				HBEvalType const hbe((*HBEval_lookup)(hbdon, hbacc, hbseq_sep));
				if ( hbe == hbe_UNKNOWN ) {
					continue;
				}

				hbond_evaluation_statement.bind(1,database_tag);
				hbond_evaluation_statement.bind(2,don_chem_type);
				hbond_evaluation_statement.bind(3,acc_chem_type);
				hbond_evaluation_statement.bind(4,separation);
				hbond_evaluation_statement.bind(5,AHdist_short_fade_lookup(hbe)->get_name());
				hbond_evaluation_statement.bind(6,AHdist_long_fade_lookup(hbe)->get_name());
				hbond_evaluation_statement.bind(7,cosBAH_fade_lookup(hbe)->get_name());
				hbond_evaluation_statement.bind(8,cosAHD_fade_lookup(hbe)->get_name());
				hbond_evaluation_statement.bind(9,AHdist_poly_lookup(hbe)->name());
				hbond_evaluation_statement.bind(10,cosBAH_short_poly_lookup(hbe)->name());
				hbond_evaluation_statement.bind(11,cosBAH_long_poly_lookup(hbe)->name());
				hbond_evaluation_statement.bind(12,cosAHD_short_poly_lookup(hbe)->name());
				hbond_evaluation_statement.bind(13,cosAHD_long_poly_lookup(hbe)->name());
				hbond_evaluation_statement.bind(14,HBondTypeManager::name_from_weight_type(weight_type_lookup(hbe)));
				basic::database::safely_write_to_database(hbond_evaluation_statement);
			}
		}
	}
	return 0;
}

} // geometric_solvation
} // namespace scoring
} // namespace core

