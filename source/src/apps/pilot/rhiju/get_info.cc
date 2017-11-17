// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <devel/init.hh>
#include <core/options/option.hh>
#include <core/util/Tracer.hh>
#include <core/util/datacache/BasicDataCache.hh>
#include <core/util/datacache/CacheableString.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/io/pose_stream/PoseInputStream.fwd.hh>
#include <core/io/pose_stream/PDBPoseInputStream.hh>
#include <core/io/pose_stream/SilentFilePoseInputStream.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>
#include <protocols/jumping/Dssp.hh>
#include <protocols/viewer/viewers.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <core/options/option_macros.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>
#include <core/options/keys/score.OptionKeys.gen.hh>
#include <numeric/xyzVector.hh>

#include <string>


OPT_KEY( Integer, seq_sep_cutoff )
OPT_KEY( Real, freq_cutoff )
OPT_KEY( Real, max_cst_strength )

using namespace core;
using namespace core::options;
using namespace core::options::OptionKeys;
using ObjexxFCL::FArray2D;
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS


////////////////////////////////////////////////////////////////////////////////////
void
save_ss_info( FArray2D< Real > & ds, pose::Pose & pose ){

	// get DSSP, assign secondary structure
	protocols::jumping::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );

	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( pose.secstruct( i ) == 'H' ){
			ds( i,1 ) += 1.0;
		} else if ( pose.secstruct( i ) == 'E' ){
			ds( i,2 ) += 1.0;
		} else {
			ds( i,3 ) += 1.0;
		}
	}

}

Distance const DIST_CUTOFF = 9.0;
Distance const FADE = 4.0;

////////////////////////////////////////////////////////////////////////////////////
void
save_contact_info( FArray2D< Real > & ds, pose::Pose & pose ){

	using namespace core::id;
	using namespace core::scoring;
	using namespace core::scoring::hbonds;

	for ( Size i = 1; i <= pose.size(); i++ ) {
		Vector const & vec_i = pose.xyz( NamedAtomID( " CA ", i ) );

		for ( Size j = i; j <= pose.size(); j++ ) {
			Vector const & vec_j = pose.xyz( NamedAtomID( " CA ", j ) );

			if ( ( vec_i - vec_j ).length() < DIST_CUTOFF ){
				ds( i, j ) += 1.0;
			}

		}
	}


	//////////////////////////////////////////////////////////
	//H-bonds! Yea!
	static ScoreFunctionOP scorefxn = get_score_function();
	(*scorefxn)( pose );

	FArray2D< bool > hb( pose.size(), pose.size(), false );

	HBondOptionsOP hbond_options( new HBondOptions() );
	hbond_options->use_hb_env_dep( false );
	HBondSetOP hbond_set( new HBondSet( hbond_options ) );

	fill_hbond_set( pose, false /*calc deriv*/, *hbond_set );

	//	std::cout << "NUM HBONDS: " << hbond_set->nhbonds() << std::endl;

	for (Size i = 1; i <= hbond_set->nhbonds(); i++ ) {
		HBond const & hbond( hbond_set->hbond( i ) );

		if ( !hbond.don_hatm_is_protein_backbone() ) continue;
		if ( !hbond.acc_atm_is_protein_backbone() ) continue;

		Size const don_res_num = hbond.don_res();
		Size const don_hatm = hbond.don_hatm();

		Size const acc_res_num = hbond.acc_res();
		Size const acc_atm = hbond.acc_atm();

		hb( don_res_num, acc_res_num ) = true;

	}


	for ( Size i = 1; i <= pose.size(); i++ ) {
		for ( Size j = 1; j < i; j++ ) {
			if ( hb(i,j) || hb(j,i) ) {
				ds( i, j ) += 1.0;
			}
		}
	}


}


////////////////////////////////////////////////////////////////////////////////////
void
normalize_info( FArray2D< Real > & ds, Size const count ){

	assert( count > 0 );

	for ( Size i = 1; i <= ds.size1(); i++ ){
		for ( Size j = 1; j <= ds.size2(); j++ ){
			ds(i,j) /= count;
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////
void
output_ss_info( FArray2D< Real > & ds, 	utility::io::ozstream & out){
	for ( Size i = 1; i <= ds.size1(); i++ ){
		out << "DS " << (i-1) << "  H " << ds(i,1) << "  E " << ds(i,2) << "  L " << ds(i,3) << std::endl;
	}
}

////////////////////////////////////////////////////////////////////////////////////
void
output_native_ss_info( pose::Pose & pose, 	utility::io::ozstream & out){

	// get DSSP, assign secondary structure
	protocols::jumping::Dssp dssp_obj( pose );
	dssp_obj.insert_ss_into_pose( pose );
	out << "NS ";
	for ( Size i = 1; i <= pose.size(); i++ ) {
		out << pose.secstruct( i );
	}
	out << std::endl;

}

////////////////////////////////////////////////////////////////////////////////////
void
output_contact_info( FArray2D< Real > & dc,  utility::io::ozstream & out, std::string const tag = "DC"){
	for ( Size i = 1; i <= dc.size1(); i++ ) {
		for ( Size j = 1; j <= dc.size2(); j++ ) {
			if( dc(i,j) > 0.0 ) {
				out << tag << ' ' << i << ' ' << j << ' ' << dc(i,j) << std::endl;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////////
void
output_csts( FArray2D< Real > & dc,  utility::io::ozstream & out ){

	out << "[ atompairs ]" << std::endl;

	Real const FREQ_CUTOFF = option[ freq_cutoff ]();
	Size const SEQ_SEP_CUTOFF = option[ seq_sep_cutoff ]();
	Real const MAX_CST_STRENGTH = option[ max_cst_strength ]();

	for ( Size i = 1; i <= dc.size1(); i++ ) {
		for ( Size j = 1; j <= dc.size2(); j++ ) {


			if( (dc(i,j) > FREQ_CUTOFF) && ( j > i ) && ( (j - i) > SEQ_SEP_CUTOFF )  ) {

				Real const cst_strength = -1.0 * MAX_CST_STRENGTH * ( dc(i,j) - FREQ_CUTOFF );

				out << " CA " << I(3,i) << " CA " << I(3,j)
						<< "  FADE " << F(8,3, (-1 * FADE) )
						<< " " << F(8,3, (DIST_CUTOFF + FADE))
						<< " " << F(8,3,FADE) << " "
						<< F(8,3,cst_strength) << std::endl;

			}

		}
	}

}

////////////////////////////////////////////////////////////////////////////////////
void
get_info_test(){

	using namespace core::chemical;
	using namespace core::io::pose_stream;
	using namespace core::options;
	using namespace core::options::OptionKeys;

	core::chemical::ResidueTypeSetCAP rsd_set = ChemicalManager::get_instance()->residue_type_set(
		option[ out::file::residue_type_set ]()
	);

	PoseInputStreamOP input;


	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
			);
		} else {
			input = new SilentFilePoseInputStream( option[ in::file::silent ]() );
		}
	} else {
		assert( option[ in::file::s ].user() );
		input =  new PDBPoseInputStream( option[ in::file::s]() );
	}

	if (!input) utility_exit_with_message( "Must specify -in::file::s or -in::file::silent" );

	core::pose::Pose pose;
	Size count( 0 );

	FArray2D< Real > ds, dc;
	Size nres( 0 );
	while ( input->has_another_pose() ) {

		input->fill_pose( pose, *rsd_set );

		if ( count == 0 ){
			nres = pose.size();
			ds.dimension( nres, 3, 0.0 );
			dc.dimension( nres, nres, 0.0 );
			protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 200, 200 );
		}
		if ( pose.size() != nres ) continue;

		save_ss_info( ds, pose );
		save_contact_info( dc, pose );

		count++;
	}

	normalize_info( ds, count );
	normalize_info( dc, count );

	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	output_ss_info( ds, out);
	output_contact_info( dc, out);

	std::string outfile2  = outfile + ".cst";
	utility::io::ozstream out2( outfile2 );
	output_csts( dc, out2 );

	//////////////////////////////////////////////
	// Need to also put in native info
	//////////////////////////////////////////////
	if ( option[ in::file::native ].user() ){

		io::pdb::pose_from_file( pose, *rsd_set, option[ in::file::native ](), core::import_pose::PDB_file);

		output_native_ss_info( pose, out );

		dc = 0.0;
		save_contact_info( dc, pose );
		output_contact_info( dc, out, "NC" );
	}


}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace core::options;

	get_info_test();

	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

	using namespace core::options;

	NEW_OPT( seq_sep_cutoff, "for constraints, minimum sequence separation", 8 );
	NEW_OPT( max_cst_strength, "for constraints, maximum constraint strength", 1.0 );
	NEW_OPT( freq_cutoff, "for constraints, minimum frequency to output constraint", 0.2 );


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

	exit( 0 );

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
