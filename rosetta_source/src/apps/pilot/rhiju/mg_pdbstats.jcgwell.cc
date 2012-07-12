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
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/import_pose/import_pose.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

#include <ObjexxFCL/format.hh>

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
#include <numeric/xyz.functions.hh>

using namespace core;
using namespace ObjexxFCL::fmt;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;


typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////////
void
mg_pdbstats_from_pose( utility::io::ozstream & out, utility::io::ozstream & out2, pose::Pose & pose, Size const count, Size & total_residues )
{

	using namespace core::conformation;
	using namespace core::chemical;

	static Real const DIST_CUTOFF ( 4.5 );
	static Real const LIGAND_DIST_CUTOFF ( 3.0 );

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	// find all the magnesiums in here.
	for (Size i = 1; i <= nres; i++) {

		if ( pose.residue( i ).name3() != " MG" ) continue;
		Vector xyz_mg = pose.residue( i ).xyz( "MG  " );

		res_count++;
		Size nligands( 0 );

		for ( Size j = 1; j <= nres; j++ ){
			if ( i == j) continue;

			Residue const & rsd_j = pose.residue( j );

			// look through all non-hydrogen atoms
			for ( Size jj = 1; jj <= rsd_j.nheavyatoms(); jj++ ){
				Vector xyz_jj = rsd_j.xyz( jj );

				Real distance = (xyz_jj - xyz_mg).length();
				if ( distance < DIST_CUTOFF ) {

					out2 << I(4, count) << " " <<  I(4, i) << " " << I(4, j) << " " << I(4, jj) << " " << F(8, 3, distance) << " " << I(4, rsd_j.atom_type_set().atom_type_index(  rsd_j.atom_type( jj ).name() ) ) << std::endl;

					if ( distance < LIGAND_DIST_CUTOFF ) {
						nligands++;
					}

				}
			}

		}

		out << I(4, count )	<< ' ' << I(4,i) << ' ' << nligands	<< std::endl;

	} // i

	std::cout << "Processed " << res_count <<  std::endl;

	total_residues += res_count;

}



///////////////////////////////////////////////////////////////////////////////
// Go over the input pdbs in a loop.
///////////////////////////////////////////////////////////////////////////////
void
mg_pdbstats_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	std::string const file_path( option[ in::path::pdb ]( 1 ) );
	std::string const pdb_list( option[ in::file::l ](1) );

	utility::io::izstream instream( pdb_list );
	if (!instream){
		std::cerr  << "Can't find list file " << pdb_list << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		return;
	}

	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile );
	utility::io::ozstream out2( "ligandinfo_"+outfile );

	std::string pdb_file;
	Size count( 0 );

	pose::Pose pose;

	Size total_residues( 0 );

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	std::string line;
	while ( 	getline( instream, line )  ) {
		std::istringstream line_stream( line );

		line_stream >> pdb_file;

		import_pose::pose_from_pdb( pose, *rsd_set,  file_path + '/' + pdb_file );

		pose.dump_pdb( "imported_"+pdb_file );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;

		mg_pdbstats_from_pose( out, out2, pose, count, total_residues );

		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	mg_pdbstats_test();

	exit( 0 );

}
