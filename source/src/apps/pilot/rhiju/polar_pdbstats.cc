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
#include <core/scoring/rna/RNA_AtomVDW.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/import_pose/import_pose.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Stub.hh>
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

OPT_KEY( Boolean, dump )


///////////////////////////////////////////////////////////////////////////////
core::kinematics::Stub
get_phosphate_stub( core::conformation::Residue const & rsd ){

	using namespace core::kinematics;

	Vector const centroid = rsd.xyz( " P  " );
	Vector const a = rsd.xyz( " OP2" );
	Vector const b = rsd.xyz( " OP1" );
	Vector y = (a + b)/2.0 - centroid;
	y.normalize();

	Vector x = b - a;
	Vector z = cross( x, y );
	z.normalize();

	x = cross( y, z );

	Matrix M = Matrix::cols( x, y, z );
	return Stub( M, centroid );

}


///////////////////////////////////////////////////////////////////////////////
void
search_other_atoms( pose::Pose const & pose, Vector const & xyz_polarH, Vector const & xyz_base, Size const count, conformation::Residue const & rsd_i, Size const ii, utility::io::ozstream & out )
{
	using namespace core::conformation;

	static Distance const DIST_CUTOFF ( 8.0 );

	Size const nres = pose.total_residue();
	Size const & i = rsd_i.seqpos();

	for ( Size j = 1; j <= nres; j++ ){
		if ( i == j) continue;

		Residue const & rsd_j = pose.residue( j );

		// look through all atoms
		for ( Size jj = 1; jj <= rsd_j.natoms(); jj++ ){

			Vector xyz_jj = rsd_j.xyz( jj );

			Real distance = (xyz_jj - xyz_polarH).length();

			if ( distance < DIST_CUTOFF ) {

				// let's also get the angles!
				Real theta( 0.0 );
				theta = numeric::conversions::degrees( angle_of( xyz_jj, xyz_polarH, xyz_base ) );

				out << I(4, count) << " " <<  I(4, i) << " " << I( 4, ii ) << " " << rsd_i.name1() << " " << rsd_i.atom_name( ii ) << "       " << I(4, j) << " " << I(4, jj)  << " " << rsd_j.name1() << " " << rsd_j.atom_name(jj ) << " " << F(8, 3, distance) << " " << theta << std::endl;

			}

		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
polar_pdbstats_from_pose( utility::io::ozstream & out1,
													utility::io::ozstream & out2,
													pose::Pose & pose, Size const count, Size & total_residues, std::string const pdb_file)
{

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring::hbonds;
	using namespace core::chemical::rna;

	Size const nres = pose.total_residue();

	// start with polar hydrogens [H-bond donors]
	std::cout << "about to start processing " << std::endl;

	HBondOptions hbond_options; // needed to derive Polar(2+)-acceptor-base angles.

	// set up nucleobase coordinate systems
	// not currently in use, but could put in later.
	utility::vector1< Stub > all_base_stubs;
	for ( Size j = 1; j <= nres; j++ ){
		Residue rsd = pose.residue( j );
		kinematics::Stub base_stub;
		if ( rsd.is_RNA() ){
			Vector centroid = get_rna_base_centroid( rsd, false /*verbose*/ );
			base_stub = Stub( get_rna_base_coordinate_system( rsd, centroid ), centroid );
		}
		all_base_stubs.push_back( base_stub );
	}

	Size res_count = 0;

	for (Size i = 1; i <= nres; i++) {

		Residue const & rsd_i = pose.residue( i );

		// find all the exocyclic amines.
		utility::vector1< Size > const & polarHlist = rsd_i.Hpos_polar();
		for ( Size qq = 1; qq <= polarHlist.size(); qq++ ){

			Size const ii = polarHlist[ qq ];
			Size const base_atom = rsd_i.atom_base( ii );
			if ( rsd_i.atom_name( base_atom )[1] != 'N' ) continue;

			res_count++;
			Vector xyz_polarH = pose.residue( i ).xyz( ii );
			Vector xyz_base   = pose.residue( i ).xyz( base_atom );

			search_other_atoms( pose, xyz_polarH, xyz_base, count, rsd_i, ii, out1 );
		}

		// find all the exocyclic amines.
		utility::vector1< Size > const & accptlist = rsd_i.accpt_pos();
		for ( Size qq = 1; qq <= accptlist.size(); qq++ ){

			Size const ii = accptlist[ qq ];

			res_count++;

			Vector xyz_polar = pose.residue( i ).xyz( ii );

			Vector dummy, xyz_base;
			chemical::Hybridization acc_hybrid( rsd_i.atom_type( ii ).hybridization());
			make_hbBasetoAcc_unitvector(
																	hbond_options,
																	acc_hybrid,
																	rsd_i.atom( ii ).xyz(),
																	rsd_i.xyz( rsd_i.atom_base( ii ) ),
																	rsd_i.xyz( rsd_i.abase2( ii ) ),
																	xyz_base, dummy );

			search_other_atoms( pose, xyz_polar, xyz_base, count, rsd_i, ii, out2 );
		}

	} // i


	std::cout << "Processed " << res_count <<  std::endl;

	total_residues += res_count;

}



///////////////////////////////////////////////////////////////////////////////
// Go over the input pdbs in a loop.
///////////////////////////////////////////////////////////////////////////////
void
polar_pdbstats_test()
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
	utility::io::ozstream out1( "donor_"+outfile );
	utility::io::ozstream out2( "acceptor_"+outfile );

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

		if ( option[dump]() ) pose.dump_pdb( "imported_"+pdb_file );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;

		polar_pdbstats_from_pose( out1, out2, pose, count, total_residues, pdb_file );

		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	NEW_OPT( dump, "dump pdb", false );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	polar_pdbstats_test();

	exit( 0 );

}
