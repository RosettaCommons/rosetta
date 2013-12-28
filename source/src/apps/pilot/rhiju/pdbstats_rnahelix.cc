// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/farna/RNA_ProtocolUtil.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/Stub.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//silly using/typedef
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/conformation/Conformation.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <numeric/xyz.functions.hh>

//Auto using namespaces
namespace ObjexxFCL { namespace format { } } using namespace ObjexxFCL::format; // AUTO USING NS
//Auto using namespaces end

using namespace core;
//using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////
std::string get_WC_atom( core::chemical::AA const & res_type ){
	using namespace core::chemical;
	std::string WC_atom( "" );
	if ( res_type == na_rad ) WC_atom = " N1 ";
	if ( res_type == na_rcy ) WC_atom = " N3 ";
	if ( res_type == na_rgu ) WC_atom = " N1 ";
	if ( res_type == na_ura ) WC_atom = " N3 ";
	return WC_atom;
}

///////////////////////////////////////////////////////////////////////////////
void
rna_helix_pdbstats_from_pose( utility::io::ozstream & out, pose::Pose & pose, Size const count, char const chain, Size & total_residues )
{

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring;
	using namespace protocols::farna;
	using namespace chemical::rna;

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	static bool init( false );

	for (Size i = 1; i <= nres; i++) {
		// get beginning and end points of helices.
	}

	// go through end of each helix, to start of next one.
	// define base-pair coordinate system, 3DNA style.
	// save 3d coordinates of next O5', O3', and next base-pair-centroid
	// save how many residues intervening until next base pair.

	total_residues += res_count;

}


///////////////////////////////////////////////////////////////////////////////
// Following is same as pdbstats.cc.
// Just created a new app to help accelerate compile time.
///////////////////////////////////////////////////////////////////////////////
void
rhiju_pdbstats()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::import_pose;

	//	utility::vector1 < std::string> pdb_files( option[ in::file::s ]() );
	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	std::string const pdb_list( option[ in::file::l ](1) );

	utility::io::izstream instream( pdb_list );
	if (!instream){
		std::cerr  << "Can't find list file " << pdb_list << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		return;
	}

	ResidueTypeSetCAP rsd_set;
	if ( option[rna_stack]() ){
  	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	} else {
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );

	}
	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile );

	//	for (Size i = 1; i <= pdb_files.size(); i++) {
	//	std::string const pdb_file = pdb_files[i];
	std::string pdb_file;
	char chain(' ');
	Size count( 0 );
	pose::Pose pose;

	Size total_residues( 0 );

	std::string line;
	while ( 	getline( instream, line )  ) {
		std::istringstream line_stream( line );

		line_stream >> pdb_file;

		line_stream >> chain;

		if ( line_stream.fail() ) chain = '?';
		//		chain =  pdb_file.at(4) ;
		//		pdb_file = pdb_file.substr(0,4);
		//		lowercase( pdb_file );

		if (chain == '_' ) chain = ' ';

		pose_from_pdb( pose, *rsd_set, file_path + '/' + pdb_file );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << " " << chain << std::endl;

		rna_helix_pdbstats_from_pose( out, pose, count, chain, total_residues );
		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	using namespace basic::options;

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	core::init::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	rhiju_pdbstats();

	exit( 0 );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}
