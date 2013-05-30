// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/hbonds/polynomial.cxxtest.hh
/// @brief  Test hbond polynomial classes
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Test Headers
#include <cxxtest/TestSuite.h>
#include <test/UTracer.hh>
#include <test/core/init_util.hh>

// Package Headers
#include <core/scoring/hbonds/types.hh>
#include <core/types.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondTypeManager.hh>
#include <core/scoring/hbonds/FadeInterval.hh>

// Project Headers
#include <core/scoring/hbonds/polynomial.hh>
#include <basic/database/open.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/string_util.hh>

// C++ Headers
#include <string>

//Auto Headers


using namespace std;
using namespace core;
using namespace core::scoring::hbonds;
using namespace utility::io;
using utility::vector1;

class HBondDatabaseTest : public CxxTest::TestSuite {

public:
  void setUp() {
    core_init();

	}

  void tearDown(){}

	void test_validate_score12_hbond_database(){

		vector1< string > tags;
		tags.push_back( "score12_params" );
		tags.push_back( "helix_hb_06_2009" );
		tags.push_back( "newOH_params" );
		tags.push_back( "His_Phil_fix" );
		tags.push_back( "extended_BAH_params" );
		tags.push_back( "newCHI_params" );
		tags.push_back( "sp2_params" );
		tags.push_back( "sp2_hackelec_params" );
		for( Size i = 1; i <= tags.size(); ++i ){
			// Validate data integrity
			HBondOptionsCOP hb_options(new HBondOptions( tags[i] ));
			string path("scoring/score_functions/hbonds/" + tags[i] );
			validate_HBAccChemType(path);
			validate_HBDonChemType(path);
			validate_HBondWeightType(path);
			validate_HBSeqSep(path);
			validate_HybridizationType(path);
			validate_HBAccHybridization(path);
			validate_HBEval(hb_options);
		}

	}

  void validate_HBAccChemType(string const & path){
    string fname( path + "/HBAccChemType.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		utility::vector1< string > tokens;
		Size line_number(1);

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 4;
			Size id;
			string name; // name_long, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Acceptors definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // id
				stringstream buf;
				buf << tokens[i];
				buf >> id;
				TS_ASSERT_EQUALS(id, line_number);
			}
			++i;

			{ // name
				stringstream buf;
				buf << tokens[i];
				buf >> name;
				HBAccChemType id_chem_type(static_cast<HBAccChemType>(id));
				TS_ASSERT_EQUALS(id_chem_type, HBondTypeManager::acc_chem_type_from_name(name));
				TS_ASSERT_EQUALS(name, HBondTypeManager::name_from_acc_chem_type(id_chem_type));
			}
			++i;

			// name_long
			++i;

			// comment

			++line_number;
		}
	}

  void validate_HBAccHybridization(string const & path){
    string fname( path + "/HBAccHybridization.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		utility::vector1< string > tokens;

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 3;
			string acc_chem_type_name, hybridization_type_name, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Acceptors definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // acc_chem_type
				stringstream buf;
				buf << tokens[i];
				buf >> acc_chem_type_name;
				TS_ASSERT(HBondTypeManager::is_acc_chem_type(acc_chem_type_name));
			}
			++i;

			{ // hybridization
				stringstream buf;
				buf << tokens[i];
				buf >> hybridization_type_name;
				TS_ASSERT(HBondTypeManager::is_hybridization_type(hybridization_type_name));
			}
			++i;

			{ // comment
			}
			// Once we have the Chemical types into the database,
			// Go through the chemical types and make sure everything is ok

		}
	}

  void validate_HBDonChemType(string const & path){
    string fname( path + "/HBDonChemType.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		Size line_number(1);
		utility::vector1< string > tokens;

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 4;
			Size id;
			string name; // name_long, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Donors definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // id
				stringstream buf;
				buf << tokens[i];
				buf >> id;
				TS_ASSERT_EQUALS(id, line_number);
			}
			++i;

			{ // name
				stringstream buf;
				buf << tokens[i];
				buf >> name;
				HBDonChemType id_chem_type(static_cast<HBDonChemType>(id));
				TS_ASSERT_EQUALS(id_chem_type, HBondTypeManager::don_chem_type_from_name(name));
				TS_ASSERT_EQUALS(name, HBondTypeManager::name_from_don_chem_type(id_chem_type));
			}
			++i;

			{ // name_long
			}
			++i;

			{ // comment
			}
			line_number++;
		}
	}

  void validate_HBondWeightType(string const & path){
    string fname( path + "/HBondWeightType.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		Size line_number(1);
		utility::vector1< string > tokens;

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 3;
			Size id;
			string name; // name_long, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Weight Type definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // id
				stringstream buf;
				buf << tokens[i];
				buf >> id;
				TS_ASSERT_EQUALS(id, line_number);
			}
			++i;

			{ // name
				stringstream buf;
				buf << tokens[i];
				buf >> name;
				HBondWeightType id_type(static_cast<HBondWeightType>(id));
				TS_ASSERT_EQUALS(id_type, HBondTypeManager::weight_type_from_name(name));
				TS_ASSERT_EQUALS(name, HBondTypeManager::name_from_weight_type(id_type));
			}
			++i;
			{ // comment
			}
			line_number++;
		}
	}

  void validate_HBEval(HBondOptionsCOP const hb_options){
		HBondDatabaseCOP hb_db( HBondDatabase::get_database( hb_options->params_database_tag()));

		FadeIntervalCOP fi;
		Polynomial_1dCOP p;
		double x, value, deriv;
		x=3;
		for(Size hbe=1; hbe != hbe_MAX; ++hbe){
			fi = hb_db->AHdist_short_fade_lookup(hbe);
			fi->value_deriv(x, value, deriv);
			fi = hb_db->AHdist_long_fade_lookup(hbe);
			fi->value_deriv(x, value, deriv);
			fi = hb_db->cosBAH_fade_lookup(hbe);
			fi->value_deriv(x, value, deriv);
			fi = hb_db->cosAHD_fade_lookup(hbe);
			fi->value_deriv(x, value, deriv);

			p = hb_db->AHdist_poly_lookup(hbe);
			TS_ASSERT_EQUALS(p->geometric_dimension(), hbgd_AHdist);
			(*p)(x, value, deriv);
			p = hb_db->cosBAH_short_poly_lookup(hbe);
			TS_ASSERT_EQUALS(p->geometric_dimension(), hbgd_cosBAH);
			(*p)(x, value, deriv);
			p = hb_db->cosBAH_long_poly_lookup(hbe);
			TS_ASSERT_EQUALS(p->geometric_dimension(), hbgd_cosBAH);
			(*p)(x, value, deriv);
			p = hb_db->cosAHD_short_poly_lookup(hbe);
			TS_ASSERT(p->geometric_dimension() == hbgd_cosAHD || p->geometric_dimension() == hbgd_AHD );
			(*p)(x, value, deriv);
			p = hb_db->cosAHD_long_poly_lookup(hbe);
			TS_ASSERT(p->geometric_dimension() == hbgd_cosAHD || p->geometric_dimension() == hbgd_AHD );
			(*p)(x, value, deriv);
		}
	}

  void validate_HBSeqSep(string const & path){
    string fname( path + "/HBSeqSep.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		Size line_number(1);
		utility::vector1< string > tokens;

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 3;
			Size id;
			string name; // name_long, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Separation Type Type definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // id
				stringstream buf;
				buf << tokens[i];
				buf >> id;
				TS_ASSERT_EQUALS(id, line_number);
			}
			++i;

			{ // name
				stringstream buf;
				buf << tokens[i];
				buf >> name;
				HBSeqSep id_type(static_cast<HBSeqSep>(id));
				TS_ASSERT_EQUALS(id_type, HBondTypeManager::seq_sep_type_from_name(name));
				TS_ASSERT_EQUALS(name, HBondTypeManager::name_from_seq_sep_type(id_type));
			}
			++i;
			{ // comment
			}
			line_number++;
		}
	}

  void validate_HybridizationType(string const & path){
    string fname( path + "/HybridizationType.csv");
		izstream s;
		basic::database::open(s, fname);
		string line;
		Size line_number(1);
		utility::vector1< string > tokens;

		while( getline( s, line ) ) {
			tokens = utility::string_split( line, ',');
			Size ntokens = 3;
			Size id;
			string name; // name_long, comment;
			if (tokens.size() != ntokens){
				stringstream message;
				message << "HBond Separation Type Type definition file does not have correct number of fields" << endl;
				message << "Expected " << ntokens << " tokens but found " << tokens.size() << " tokens. " << endl;
				message << "Path: '" << fname << "'" << endl;
				message << "Line: '" << line << "'" << endl;
				TS_ASSERT(false);
			}

			Size i(1);
			{ // id
				stringstream buf;
				buf << tokens[i];
				buf >> id;
				TS_ASSERT_EQUALS(id, line_number);
			}
			++i;

			{ // name
				stringstream buf;
				buf << tokens[i];
				buf >> name;
				chemical::Hybridization id_type(static_cast<chemical::Hybridization>(id));
				TS_ASSERT_EQUALS(id_type, HBondTypeManager::hybridization_type_from_name(name));
				TS_ASSERT_EQUALS(name, HBondTypeManager::name_from_hybridization_type(id_type));
			}
			++i;
			{ // comment
			}
			line_number++;
		}
	}

};
