// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file utility/options/OptionCollection.cxxtest.hh
/// @brief test suite for options system
/// @author Matthew O'Meara mattjomeara@gmail.com

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/protocols/init_util.hh>

// Project headers
#include <basic/options/option.hh>
#include <basic/options/util.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// For the core_init_w_additional_options call.
#include <test/core/init_util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/options/keys/all.hh>
#include <utility/options/keys/OptionKey.hh>
#include <utility/options/OptionCollection.hh>

// C++ headers
#include <string>


utility::options::BooleanOptionKey       b1( "b1" );
utility::options::BooleanOptionKey       b2( "b2" );
utility::options::BooleanOptionKey       b3( "b3" );
utility::options::BooleanOptionKey       b4( "b4" );
utility::options::BooleanOptionKey       b5( "b5" );
utility::options::IntegerOptionKey       i1( "i1" );
utility::options::IntegerOptionKey       i2( "i2" );
utility::options::IntegerOptionKey       i3( "i3" );
utility::options::IntegerOptionKey       i4( "i4" );
utility::options::IntegerOptionKey       i5( "i5" );
utility::options::RealOptionKey          r1( "r1" );
utility::options::RealOptionKey          r2( "r2" );
utility::options::RealOptionKey          r3( "r3" );
utility::options::RealOptionKey          r4( "r4" );
utility::options::RealOptionKey          r5( "r5" );
utility::options::StringOptionKey        s1( "s1" );
utility::options::StringOptionKey        s2( "s2" );
utility::options::StringOptionKey        s3( "s3" );
utility::options::StringOptionKey        s4( "s4" );
utility::options::StringOptionKey        s5( "s5" );
utility::options::FileOptionKey          f1( "f1" );
utility::options::FileOptionKey          f2( "f2" );
utility::options::FileOptionKey          f3( "f3" );
utility::options::FileOptionKey          f4( "f4" );
utility::options::FileOptionKey          f5( "f5" );
utility::options::PathOptionKey          p1( "p1" );
utility::options::PathOptionKey          p2( "p2" );
utility::options::PathOptionKey          p3( "p3" );
utility::options::PathOptionKey          p4( "p4" );
utility::options::PathOptionKey          p5( "p5" );
utility::options::BooleanVectorOptionKey bv1( "bv1" );
utility::options::BooleanVectorOptionKey bv2( "bv2" );
utility::options::BooleanVectorOptionKey bv3( "bv3" );
utility::options::BooleanVectorOptionKey bv4( "bv4" );
utility::options::BooleanVectorOptionKey bv5( "bv5" );
utility::options::IntegerVectorOptionKey iv1( "iv1" );
utility::options::IntegerVectorOptionKey iv2( "iv2" );
utility::options::IntegerVectorOptionKey iv3( "iv3" );
utility::options::IntegerVectorOptionKey iv4( "iv4" );
utility::options::IntegerVectorOptionKey iv5( "iv5" );
utility::options::RealVectorOptionKey    rv1( "rv1" );
utility::options::RealVectorOptionKey    rv2( "rv2" );
utility::options::RealVectorOptionKey    rv3( "rv3" );
utility::options::RealVectorOptionKey    rv4( "rv4" );
utility::options::RealVectorOptionKey    rv5( "rv5" );
utility::options::StringVectorOptionKey  sv1( "sv1" );
utility::options::StringVectorOptionKey  sv2( "sv2" );
utility::options::StringVectorOptionKey  sv3( "sv3" );
utility::options::StringVectorOptionKey  sv4( "sv4" );
utility::options::StringVectorOptionKey  sv5( "sv5" );
utility::options::FileVectorOptionKey    fv1( "fv1" );
utility::options::FileVectorOptionKey    fv2( "fv2" );
utility::options::FileVectorOptionKey    fv3( "fv3" );
utility::options::FileVectorOptionKey    fv4( "fv4" );
utility::options::FileVectorOptionKey    fv5( "fv5" );
utility::options::PathVectorOptionKey    pv1( "pv1" );
utility::options::PathVectorOptionKey    pv2( "pv2" );
utility::options::PathVectorOptionKey    pv3( "pv3" );
utility::options::PathVectorOptionKey    pv4( "pv4" );
utility::options::PathVectorOptionKey    pv5( "pv5" );
utility::options::ResidueChainVectorOptionKey rcv1( "rcv1" );
utility::options::ResidueChainVectorOptionKey rcv2( "rcv2" );
utility::options::ResidueChainVectorOptionKey rcv3( "rcv3" );
utility::options::ResidueChainVectorOptionKey rcv4( "rcv4" );
utility::options::ResidueChainVectorOptionKey rcv5( "rcv5" );



class OptionsSystemTests : public CxxTest::TestSuite {

public:

	void setUp() {
		basic::options::initialize();
	}

	void test_find_key_cl_good(){
		using namespace basic::options;
		try{
			std::string key(option.find_key_cl("in:file:s", "", true));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "in:file", false));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "", false));
		} catch (...){
			TS_ASSERT(false);
		}

		try{
			std::string key(option.find_key_cl("s", "", true));
		} catch (...){
			TS_ASSERT(false);
		}


		try{
			std::string key(option.find_key_cl("in:file:s", "", false));
		} catch (...){
			TS_ASSERT(false);
		}


	}

	void test_find_key_cl_error() {
		using namespace basic::options;
		try{
			std::string key(option.find_key_cl(":::s", "in:file", true));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			// unique best suffix match
			std::string key(option.find_key_cl("asdfasdfasf:s", "in:file", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			// ignore context when 'top' == true
			std::string key(option.find_key_cl("s", "in:file", true));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("in::file::s", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("-in::file::s-", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("in:file", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

		try{
			std::string key(option.find_key_cl("-in:file:s", "", false));
		} catch (...){
			TS_ASSERT(true);
		}

	}

	void test_initialize_local_option_collection_from_global() {
		using namespace utility::options;
		using namespace basic::options;

		// pattern here;
		// option #1 -- on CL; (scalar) value / (vector) one element value; added to local option collection
		// option #2 -- on CL; (scalar) value / (vector) multiple element value; added to local option collection
		// option #3 -- not on CL; has default value; added to local option collection, active but not set by user()
		// option #4 -- on CL; (scalar) value / (vector) multiple element value; unknown to the local option collection
		// option #5 -- not on CL; added to local option collection, but not active
		option.add( b1,  "" );
		option.add( b2,  "" );
		option.add( b3,  "" ).def( false );
		option.add( b4,  "" );
		option.add( b5,  "" ); // note that boolean options have a default "false" value unless you ask for one that doesn't.
		option.add( i1,  "" );
		option.add( i2,  "" );
		option.add( i3,  "" ).def( 15 );
		option.add( i4,  "" );
		option.add( i5,  "" );
		option.add( r1,  "" );
		option.add( r2,  "" );
		option.add( r3,  "" ).def( 2.625 );
		option.add( r4,  "" );
		option.add( r5,  "" );
		option.add( s1,  "" );
		option.add( s2,  "" );
		option.add( s3,  "" ).def( "youcandoit!" );
		option.add( s4,  "" );
		option.add( s5,  "" );
		option.add( f1,  "" );
		option.add( f2,  "" );
		option.add( f3,  "" ).def( "score12.wts" );
		option.add( f4,  "" );
		option.add( f5,  "" );
		option.add( p1,  "" );
		option.add( p2,  "" );
		option.add( p3,  "" ).def( "/home/andrew/rosetta/all_the_good_stuff/" );
		option.add( p4,  "" );
		option.add( p5,  "" );
		option.add( bv1, "" );
		option.add( bv2, "" );
		option.add( bv3, "" ).def(); // empty vector default
		option.add( bv4, "" );
		option.add( bv5, "" );
		option.add( iv1, "" );
		option.add( iv2, "" );
		option.add( iv3, "" ).def( 12 );
		option.add( iv4, "" );
		option.add( iv5, "" );
		option.add( rv1, "" );
		option.add( rv2, "" );
		option.add( rv3, "" ).def( 12.5 );
		option.add( rv4, "" );
		option.add( rv5, "" );
		option.add( sv1, "" );
		option.add( sv1, "" );
		option.add( sv2, "" );
		option.add( sv3, "" ).def( "failblog" );
		option.add( sv4, "" );
		option.add( sv5, "" );
		option.add( fv1, "" );
		option.add( fv2, "" );
		option.add( fv3, "" ).def( "talaris2014.wts" );
		option.add( fv4, "" );
		option.add( fv5, "" );
		option.add( pv1, "" );
		option.add( pv2, "" );
		option.add( pv3, "" ).def( "/home/sweet/home" );
		option.add( pv4, "" );
		option.add( pv5, "" );

		option.add( rcv1, "" );
		option.add( rcv2, "" );
		utility::vector1< int > rcv3_default( 5 );
		for ( int ii = 1; ii <= 5; ++ii ) rcv3_default[ ii ] = ii;
		option.add( rcv3, "" ).def( rcv3_default );
		option.add( rcv4, "" );
		option.add( rcv5, "" );

		core_init_with_additional_options( "-b1 false -b2 -b4 -i1 1234 -i2 2345 -i4 3456 -r1 1.25 -r2 2.5 -r4 2.75 "
			" -s1 testing -s2 yomomma -s4 haha -f1 protocol.xml -f2 /usr/bin/gcc -f4 native.cst"
			" -p1 /usr/include -p2 /home/andrew/GIT -p4 /home -bv1 true -bv2 true true false -bv4 true false"
			" -iv1 1234 -iv2 123 234 345 456 -iv4 123 111 -rv1 1.25 -rv2 2.125 3.5 -rv4 12.5 13.5 15.5"
			" -sv1 test -sv2 test foo bar baz -sv4 oh hey yeah about that -fv1 des1.resfile -fv2"
			" /path/to/des2.resfile /other/path/to/des3.resfile des4.resfile -fv4 native.pdb"
			" -pv1 /home/andrew/GIT -pv2 /home/andrew/GIT/main1/database /home/andrew/GIT/main2/database -pv4 /usr/bin /usr/lib"
			" -rcv1 A:1 -rcv2 A:5-9 B:11-15" );

		OptionKeyList option_list;
		option_list + b1 + b2 + b3 + b5 + i1 + i2 + i3 + i5 + r1 + r2 + r3 + r5 + s1 + s2 + s3 + s5 +
			f1 + f2 + f3 + f5 + p1 + p2 + p3 + p5 + bv1 + bv2 + bv3 + bv5 +
			iv1 + iv2 + iv3 + iv5 + rv1 + rv2 + rv3 + rv5 + sv1 + sv2 + sv3 + sv5 +
			fv1 + fv2 + fv3 + fv5 + pv1 + pv2 + pv3 + pv5 + rcv1 + rcv2 + rcv3 + rcv5;

		OptionCollectionOP local_options_ptr = read_subset_of_global_option_collection( option_list );
		TS_ASSERT( local_options_ptr );
		if ( ! local_options_ptr ) return;
		OptionCollection const & local_options( *local_options_ptr );

		TS_ASSERT( local_options.has( b1 ));
		TS_ASSERT( local_options[ b1 ].active() );
		TS_ASSERT( local_options[ b1 ].user() );
		TS_ASSERT( local_options[ b1 ] == false );

		TS_ASSERT( local_options.has( b2 ));
		TS_ASSERT( local_options[ b2 ].active() );
		TS_ASSERT( local_options[ b2 ].user() );
		TS_ASSERT( local_options[ b2 ] == true );

		TS_ASSERT( local_options.has( b3 ));
		TS_ASSERT( local_options[ b3 ].active() );
		TS_ASSERT( ! local_options[ b3 ].user() );
		TS_ASSERT( local_options[ b3 ] == false );

		TS_ASSERT( ! local_options.has( b4 ));

		TS_ASSERT( local_options.has( b5 ));
		TS_ASSERT( local_options[ b5 ].active() );
		TS_ASSERT( ! local_options[ b5 ].user() );
		TS_ASSERT( local_options[ b5 ] == false );

		TS_ASSERT( local_options.has( i1 ));
		TS_ASSERT( local_options[ i1 ].active() );
		TS_ASSERT( local_options[ i1 ].user() );
		TS_ASSERT( local_options[ i1 ] == 1234 );

		TS_ASSERT( local_options.has( i2 ));
		TS_ASSERT( local_options[ i2 ].active() );
		TS_ASSERT( local_options[ i2 ].user() );
		TS_ASSERT( local_options[ i2 ] == 2345 );

		TS_ASSERT( local_options.has( i3 ));
		TS_ASSERT( local_options[ i3 ].active() );
		TS_ASSERT( ! local_options[ i3 ].user() );
		TS_ASSERT( local_options[ i3 ] == 15 );

		TS_ASSERT( ! local_options.has( i4 ));

		TS_ASSERT( local_options.has( i5 ));
		TS_ASSERT( ! local_options[ i5 ].active() );
		TS_ASSERT( ! local_options[ i5 ].user() );


		TS_ASSERT( local_options.has( r1 ));
		TS_ASSERT( local_options[ r1 ].active() );
		TS_ASSERT( local_options[ r1 ].user() );
		TS_ASSERT( local_options[ r1 ] == 1.25 );

		TS_ASSERT( local_options.has( r2 ));
		TS_ASSERT( local_options[ r2 ].active() );
		TS_ASSERT( local_options[ r2 ].user() );
		TS_ASSERT( local_options[ r2 ] == 2.5 );

		TS_ASSERT( local_options.has( r3 ));
		TS_ASSERT( local_options[ r3 ].active() );
		TS_ASSERT( ! local_options[ r3 ].user() );
		TS_ASSERT( local_options[ r3 ] == 2.625 );

		TS_ASSERT( ! local_options.has( r4 ));

		TS_ASSERT( local_options.has( r5 ));
		TS_ASSERT( ! local_options[ r5 ].active() );
		TS_ASSERT( ! local_options[ r5 ].user() );


		TS_ASSERT( local_options.has( s1 ));
		TS_ASSERT( local_options[ s1 ].active() );
		TS_ASSERT( local_options[ s1 ].user() );
		TS_ASSERT( local_options[ s1 ]() == "testing" );

		TS_ASSERT( local_options.has( s2 ));
		TS_ASSERT( local_options[ s2 ].active() );
		TS_ASSERT( local_options[ s2 ].user() );
		TS_ASSERT( local_options[ s2 ]() == "yomomma" );

		TS_ASSERT( local_options.has( s3 ));
		TS_ASSERT( local_options[ s3 ].active() );
		TS_ASSERT( ! local_options[ s3 ].user() );
		TS_ASSERT( local_options[ s3 ]() == "youcandoit!" );

		TS_ASSERT( ! local_options.has( s4 ));

		TS_ASSERT( local_options.has( s5 ));
		TS_ASSERT( ! local_options[ s5 ].active() );
		TS_ASSERT( ! local_options[ s5 ].user() );


		TS_ASSERT( local_options.has( f1 ));
		TS_ASSERT( local_options[ f1 ].active() );
		TS_ASSERT( local_options[ f1 ].user() );
		TS_ASSERT( local_options[ f1 ]() == utility::file::FileName( "protocol.xml" ));

		TS_ASSERT( local_options.has( f2 ));
		TS_ASSERT( local_options[ f2 ].active() );
		TS_ASSERT( local_options[ f2 ].user() );
		TS_ASSERT( local_options[ f2 ]() == utility::file::FileName( "/usr/bin/gcc" ));

		TS_ASSERT( local_options.has( f3 ));
		TS_ASSERT( local_options[ f3 ].active() );
		TS_ASSERT( ! local_options[ f3 ].user() );
		TS_ASSERT( local_options[ f3 ]() == utility::file::FileName( "score12.wts" ));

		TS_ASSERT( ! local_options.has( f4 ));

		TS_ASSERT( local_options.has( f5 ));
		TS_ASSERT( ! local_options[ f5 ].active() );
		TS_ASSERT( ! local_options[ f5 ].user() );

		TS_ASSERT( local_options.has( p1 ));
		TS_ASSERT( local_options[ p1 ].active() );
		TS_ASSERT( local_options[ p1 ].user() );
		TS_ASSERT( local_options[ p1 ]() == utility::file::PathName( "/usr/include" ));

		TS_ASSERT( local_options.has( p2 ));
		TS_ASSERT( local_options[ p2 ].active() );
		TS_ASSERT( local_options[ p2 ].user() );
		TS_ASSERT( local_options[ p2 ]() == utility::file::PathName( "/home/andrew/GIT" ));

		TS_ASSERT( local_options.has( p3 ));
		TS_ASSERT( local_options[ p3 ].active() );
		TS_ASSERT( ! local_options[ p3 ].user() );
		TS_ASSERT( local_options[ p3 ]() == utility::file::PathName( "/home/andrew/rosetta/all_the_good_stuff/" ));

		TS_ASSERT( ! local_options.has( p4 ));

		TS_ASSERT( local_options.has( p5 ));
		TS_ASSERT( ! local_options[ p5 ].active() );
		TS_ASSERT( ! local_options[ p5 ].user() );
		TS_ASSERT( local_options.has( bv1 ));
		TS_ASSERT( local_options[ bv1 ].active() );
		TS_ASSERT( local_options[ bv1 ].user() );
		utility::vector1< bool > bv1_expected( 1, true );
		TS_ASSERT( local_options[ bv1 ]() == bv1_expected );

		TS_ASSERT( local_options.has( bv2 ));
		TS_ASSERT( local_options[ bv2 ].active() );
		TS_ASSERT( local_options[ bv2 ].user() );
		utility::vector1< bool > bv2_expected( 3 );
		bv2_expected[ 1 ] = true;
		bv2_expected[ 2 ] = true;
		bv2_expected[ 3 ] = false;
		TS_ASSERT( local_options[ bv2 ]() == bv2_expected );

		TS_ASSERT( local_options.has( bv3 ));
		TS_ASSERT( local_options[ bv3 ].active() );
		TS_ASSERT( ! local_options[ bv3 ].user() );
		utility::vector1< bool > bv3_expected;
		TS_ASSERT( local_options[ bv3 ]() == bv3_expected );

		TS_ASSERT( ! local_options.has( bv4 ));

		TS_ASSERT( local_options.has( bv5 ));
		TS_ASSERT( ! local_options[ bv5 ].active() );
		TS_ASSERT( ! local_options[ bv5 ].user() );
		TS_ASSERT( local_options.has( iv1 ));
		TS_ASSERT( local_options[ iv1 ].active() );
		TS_ASSERT( local_options[ iv1 ].user() );
		utility::vector1< int > iv1_expected( 1, 1234 );
		TS_ASSERT( local_options[ iv1 ]() ==  iv1_expected );

		TS_ASSERT( local_options.has( iv2 ));
		TS_ASSERT( local_options[ iv2 ].active() );
		TS_ASSERT( local_options[ iv2 ].user() );
		utility::vector1< int > iv2_expected( 4 );
		iv2_expected[ 1 ] = 123;
		iv2_expected[ 2 ] = 234;
		iv2_expected[ 3 ] = 345;
		iv2_expected[ 4 ] = 456;
		TS_ASSERT( local_options[ iv2 ]() == iv2_expected );

		TS_ASSERT( local_options.has( iv3 ));
		TS_ASSERT( local_options[ iv3 ].active() );
		TS_ASSERT( ! local_options[ iv3 ].user() );
		utility::vector1< int > iv3_expected( 1, 12 );
		TS_ASSERT( local_options[ iv3 ]() == iv3_expected );

		TS_ASSERT( ! local_options.has( iv4 ));

		TS_ASSERT( local_options.has( iv5 ));
		TS_ASSERT( ! local_options[ iv5 ].active() );
		TS_ASSERT( ! local_options[ iv5 ].user() );

		TS_ASSERT( local_options.has( rv1 ));
		TS_ASSERT( local_options[ rv1 ].active() );
		TS_ASSERT( local_options[ rv1 ].user() );
		utility::vector1< double > rv1_expected( 1, 1.25 );
		TS_ASSERT( local_options[ rv1 ]() == rv1_expected );

		TS_ASSERT( local_options.has( rv2 ));
		TS_ASSERT( local_options[ rv2 ].active() );
		TS_ASSERT( local_options[ rv2 ].user() );
		utility::vector1< double > rv2_expected( 2 );
		rv2_expected[ 1 ] = 2.125;
		rv2_expected[ 2 ] = 3.5;
		TS_ASSERT( local_options[ rv2 ]() == rv2_expected );

		TS_ASSERT( local_options.has( rv3 ));
		TS_ASSERT( local_options[ rv3 ].active() );
		TS_ASSERT( ! local_options[ rv3 ].user() );
		utility::vector1< double > rv3_expected( 1, 12.5 );
		TS_ASSERT( local_options[ rv3 ]() == rv3_expected );

		TS_ASSERT( ! local_options.has( rv4 ));

		TS_ASSERT( local_options.has( rv5 ));
		TS_ASSERT( ! local_options[ rv5 ].active() );
		TS_ASSERT( ! local_options[ rv5 ].user() );

		TS_ASSERT( local_options.has( sv1 ));
		TS_ASSERT( local_options[ sv1 ].active() );
		TS_ASSERT( local_options[ sv1 ].user() );
		utility::vector1< std::string > sv1_expected( 1, "test" );
		TS_ASSERT( local_options[ sv1 ]() ==  sv1_expected );

		TS_ASSERT( local_options.has( sv2 ));
		TS_ASSERT( local_options[ sv2 ].active() );
		TS_ASSERT( local_options[ sv2 ].user() );
		utility::vector1< std::string > sv2_expected( 4 );
		sv2_expected[ 1 ] = "test";
		sv2_expected[ 2 ] = "foo";
		sv2_expected[ 3 ] = "bar";
		sv2_expected[ 4 ] = "baz";
		TS_ASSERT( local_options[ sv2 ]() == sv2_expected );

		TS_ASSERT( local_options.has( sv3 ));
		TS_ASSERT( local_options[ sv3 ].active() );
		TS_ASSERT( ! local_options[ sv3 ].user() );
		utility::vector1< std::string > sv3_expected( 1, "failblog"  );
		TS_ASSERT( local_options[ sv3 ]() == sv3_expected );

		TS_ASSERT( ! local_options.has( sv4 ));

		TS_ASSERT( local_options.has( sv5 ));
		TS_ASSERT( ! local_options[ sv5 ].active() );
		TS_ASSERT( ! local_options[ sv5 ].user() );

		TS_ASSERT( local_options.has( fv1 ));
		TS_ASSERT( local_options[ fv1 ].active() );
		TS_ASSERT( local_options[ fv1 ].user() );
		utility::vector1< utility::file::FileName > fv1_expected( 1, utility::file::FileName( "des1.resfile" ));
		TS_ASSERT( local_options[ fv1 ]() == fv1_expected );

		TS_ASSERT( local_options.has( fv2 ));
		TS_ASSERT( local_options[ fv2 ].active() );
		TS_ASSERT( local_options[ fv2 ].user() );
		utility::vector1< utility::file::FileName > fv2_expected( 3 );
		fv2_expected[ 1 ] = utility::file::FileName( "/path/to/des2.resfile" );
		fv2_expected[ 2 ] = utility::file::FileName( "/other/path/to/des3.resfile" );
		fv2_expected[ 3 ] = utility::file::FileName( "des4.resfile" );
		TS_ASSERT( local_options[ fv2 ]() == fv2_expected );

		TS_ASSERT( local_options.has( fv3 ));
		TS_ASSERT( local_options[ fv3 ].active() );
		TS_ASSERT( ! local_options[ fv3 ].user() );
		utility::vector1< utility::file::FileName > fv3_expected( 1, utility::file::FileName( "talaris2014.wts" ));
		TS_ASSERT( local_options[ fv3 ]() == fv3_expected );

		TS_ASSERT( ! local_options.has( fv4 ));

		TS_ASSERT( local_options.has( fv5 ));
		TS_ASSERT( ! local_options[ fv5 ].active() );
		TS_ASSERT( ! local_options[ fv5 ].user() );

		TS_ASSERT( local_options.has( pv1 ));
		TS_ASSERT( local_options[ pv1 ].active() );
		TS_ASSERT( local_options[ pv1 ].user() );
		utility::vector1< utility::file::PathName > pv1_expected( 1, utility::file::PathName( "/home/andrew/GIT" ));
		TS_ASSERT( local_options[ pv1 ]() == pv1_expected );

		TS_ASSERT( local_options.has( pv2 ));
		TS_ASSERT( local_options[ pv2 ].active() );
		TS_ASSERT( local_options[ pv2 ].user() );
		utility::vector1< utility::file::PathName > pv2_expected( 2 );
		pv2_expected[ 1 ] = "/home/andrew/GIT/main1/database";
		pv2_expected[ 2 ] = "/home/andrew/GIT/main2/database";
		TS_ASSERT( local_options[ pv2 ]() == pv2_expected );

		TS_ASSERT( local_options.has( pv3 ));
		TS_ASSERT( local_options[ pv3 ].active() );
		TS_ASSERT( ! local_options[ pv3 ].user() );
		utility::vector1< utility::file::PathName > pv3_expected( 1, utility::file::PathName( "/home/sweet/home" ));
		TS_ASSERT( local_options[ pv3 ]() == pv3_expected );

		TS_ASSERT( ! local_options.has( pv4 ));

		TS_ASSERT( local_options.has( pv5 ));
		TS_ASSERT( ! local_options[ pv5 ].active() );
		TS_ASSERT( ! local_options[ pv5 ].user() );

		TS_ASSERT( local_options.has( rcv1 ));
		TS_ASSERT( local_options[ rcv1 ].active() );
		TS_ASSERT( local_options[ rcv1 ].user() );
		utility::vector1< int > rcv1_expected( 1, 1 );
		TS_ASSERT( local_options[ rcv1 ]() == rcv1_expected );
		std::pair< utility::vector1<int>, utility::vector1<char> > rcv1_resnum_and_chain = local_options[ rcv1 ].resnum_and_chain();
		utility::vector1< char > rcv1_expected_chain( 1, 'A' );
		TS_ASSERT_EQUALS( rcv1_resnum_and_chain.first,  rcv1_expected );
		TS_ASSERT_EQUALS( rcv1_resnum_and_chain.second, rcv1_expected_chain );

		TS_ASSERT( local_options.has( rcv2 ));
		TS_ASSERT( local_options[ rcv2 ].active() );
		TS_ASSERT( local_options[ rcv2 ].user() );
		utility::vector1< int > rcv2_expected( 10 );
		rcv2_expected[  1 ] = 5;
		rcv2_expected[  2 ] = 6;
		rcv2_expected[  3 ] = 7;
		rcv2_expected[  4 ] = 8;
		rcv2_expected[  5 ] = 9;
		rcv2_expected[  6 ] = 11;
		rcv2_expected[  7 ] = 12;
		rcv2_expected[  8 ] = 13;
		rcv2_expected[  9 ] = 14;
		rcv2_expected[ 10 ] = 15;
		TS_ASSERT( local_options[ rcv2 ]() == rcv2_expected );
		std::pair< utility::vector1<int>, utility::vector1<char> > rcv2_resnum_and_chain = local_options[ rcv2 ].resnum_and_chain();
		utility::vector1< char > rcv2_expected_chain( 10 );
		for ( int ii = 1; ii <= 10; ++ii ) rcv2_expected_chain[ ii ] = ii <= 5 ? 'A' : 'B';
		TS_ASSERT_EQUALS( rcv2_resnum_and_chain.first,  rcv2_expected );
		TS_ASSERT_EQUALS( rcv2_resnum_and_chain.second, rcv2_expected_chain );

		TS_ASSERT( local_options.has( rcv3 ));
		TS_ASSERT( local_options[ rcv3 ].active() );
		TS_ASSERT( ! local_options[ rcv3 ].user() );
		TS_ASSERT( local_options[ rcv3 ]() == rcv3_default );

		TS_ASSERT( ! local_options.has( rcv4 ));

		TS_ASSERT( local_options.has( rcv5 ));
		TS_ASSERT( ! local_options[ rcv5 ].active() );
		TS_ASSERT( ! local_options[ rcv5 ].user() );
	}

	void test_set_value() {
		using namespace utility::options;
		using namespace basic::options;

		option.add( b1,  "" );
		option.add( b3,  "" ).def( true );
		option.add( s1,  "" );
		option.add( s3,  "" ).def( "youcandoit!" );
		option.add( sv1, "" );
		option.add( sv3, "" ).def( "failblog" );
		option.add( rv1, "" );

		TS_ASSERT_EQUALS( option[b1](), false );
		option[b1].set_value( "1" );
		TS_ASSERT_EQUALS( option[b1](), true );
		TS_ASSERT_EQUALS( option[b3](), true );
		option[b3].set_cl_value( "false" );
		TS_ASSERT_EQUALS( option[b3](), false );

		option[s1].set_cl_value( "my spaced string" );
		TS_ASSERT_EQUALS( option[s1](), "my spaced string" );

		option[s3].set_value( "my spaced string" );
		TS_ASSERT_EQUALS( option[s3](), "my spaced string" );

		option[sv3].set_cl_value( "new item" );
		TS_ASSERT_EQUALS( option[sv3].size(), 1 );
		TS_ASSERT_EQUALS( option[sv3][1], "new item" );
		option[sv3].set_cl_value( "additional item" );
		TS_ASSERT_EQUALS( option[sv3].size(), 2 );
		TS_ASSERT_EQUALS( option[sv3][2], "additional item" );
		option[sv3].set_cl_value( "I haz' da\" quot-es" );
		TS_ASSERT_EQUALS( option[sv3].size(), 3 );
		TS_ASSERT_EQUALS( option[sv3][3], "I haz' da\" quot-es" );

		option[sv1].set_value("space separated");
		TS_ASSERT_EQUALS( option[sv1].size(), 2 );
		TS_ASSERT_EQUALS( option[sv1][1], "space" );
		TS_ASSERT_EQUALS( option[sv1][2], "separated" );
		option[sv1].set_value("'quotes will'keep \"things together"); // deliberately missing end quote
		TS_ASSERT_EQUALS( option[sv1].size(), 4 );
		TS_ASSERT_EQUALS( option[sv1][3], "quotes will'keep" );
		TS_ASSERT_EQUALS( option[sv1][4], "things together" );
		option[sv1].set_value("\tescap\\'age it\\'s neat\n  ", /*reset=*/ true);
		TS_ASSERT_EQUALS( option[sv1].size(), 3 ); //also check reset
		TS_ASSERT_EQUALS( option[sv1][1], "escap\\'age" );
		TS_ASSERT_EQUALS( option[sv1][2], "it\\'s" );
		TS_ASSERT_EQUALS( option[sv1][3], "neat" );

		option[rv1].set_cl_value("12.5");
		TS_ASSERT_EQUALS( option[rv1].size(), 1 );
		TS_ASSERT_EQUALS( option[rv1][1], 12.5 );
		option[rv1].set_value("1.0 23.25 7");
		TS_ASSERT_EQUALS( option[rv1].size(), 4 );
		TS_ASSERT_EQUALS( option[rv1][1], 12.5 );
		TS_ASSERT_EQUALS( option[rv1][2],  1.0 );
		TS_ASSERT_EQUALS( option[rv1][3], 23.25 );
		TS_ASSERT_EQUALS( option[rv1][4], 7 );

	}


};
