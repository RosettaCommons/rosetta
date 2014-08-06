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
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/rna/RNA_AtomVDW.hh>
#include <core/chemical/rna/util.hh>
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
using namespace ObjexxFCL::format;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::dump_pdb;


typedef  numeric::xyzMatrix< Real > Matrix;


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
mg_pdbstats_from_pose( utility::io::ozstream & out,
											 utility::io::ozstream & out2,
											 utility::io::ozstream & out3,
											 utility::io::ozstream & out4,
											 utility::io::ozstream & out5,
											 utility::io::ozstream & out6,
											 pose::Pose & pose, Size const count, Size & total_residues, std::string const pdb_file)
{

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::kinematics;
	using namespace core::scoring::hbonds;
	using namespace core::chemical::rna;

	static Real const DIST_CUTOFF ( 8.0 );
	static Real const VDW_DIST_CUTOFF ( 4.5 );
	static Real const LIGAND_DIST_CUTOFF ( 3.0 );

	Size const nres = pose.total_residue();

	Size res_count( 0 );

	core::scoring::rna::RNA_AtomVDW rna_atom_vdw;

	std::cout << "about to start processing " << std::endl;

	HBondOptions hbond_options; // needed to derive Mg(2+)-acceptor-base angles.

	// set up nucleobase coordinate systems
	utility::vector1< Stub > all_base_stubs, all_phosphate_stubs;
	for ( Size j = 1; j <= nres; j++ ){
		Residue rsd = pose.residue( j );
		kinematics::Stub base_stub, phosphate_stub;
		if ( rsd.is_RNA() ){
			Vector centroid = get_rna_base_centroid( rsd, false /*verbose*/ );
			base_stub = Stub( get_rna_base_coordinate_system( rsd, centroid ), centroid );

			phosphate_stub = get_phosphate_stub( rsd );
		}
		all_base_stubs.push_back( base_stub );
		all_phosphate_stubs.push_back( phosphate_stub );
	}



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

						//<< I(4, rsd_j.atom_type_set().atom_type_index(  rsd_j.atom_type( jj ).name() ) ) << std::endl;

					// let's also get the angles!
					Real theta( 0.0 );
					if ( rsd_j.heavyatom_is_an_acceptor( jj ) ){
						Vector dummy, xyz_base;
						chemical::Hybridization acc_hybrid( rsd_j.atom_type( jj ).hybridization());
						make_hbBasetoAcc_unitvector(
																				hbond_options,
																				acc_hybrid,
																				rsd_j.atom( jj ).xyz(),
																				rsd_j.xyz( rsd_j.atom_base( jj ) ),
																				rsd_j.xyz( rsd_j.abase2( jj ) ),
																				xyz_base, dummy );
						theta = numeric::conversions::degrees( angle_of( xyz_mg, xyz_jj, xyz_base ) );
						//						if ( rsd_j.atom_name( jj ) == " OP2" && distance < 3.0){
							//							std::cout << "OP2 " << jj << " " << " " << xyz_jj.x() << " " << xyz_base.x() << " " << xyz_mg.x() << " " << theta << std::endl;
						//						}
					}

					out2 << I(4, count) << " " <<  I(4, i) << " " << I(4, j) << " " << I(4, jj) << " " << F(8, 3, distance) << " "  << rsd_j.name1() << " " << rsd_j.atom_name(jj ) << " " << theta << std::endl;

					if ( distance < LIGAND_DIST_CUTOFF ) {
						nligands++;
					}

				}
			}

			// also make sure to look through vdw representatives.
			utility::vector1< std::string > const & vdw_atom_list =rna_atom_vdw.vdw_atom_list( rsd_j.name1() );
			for ( Size m = 1;  m <= vdw_atom_list.size(); m++ ){

				std::string const atom_name = vdw_atom_list[ m ];
				Size const jj = rsd_j.atom_index( atom_name );
				Vector xyz_jj = rsd_j.xyz( jj );

				Real distance = (xyz_jj - xyz_mg).length();
				if ( distance < VDW_DIST_CUTOFF ) {
					out3 << I(4, count) << " " <<  I(4, i) << " " << I(4, j) << " " << I(4, jj) << " " << F(8, 3, distance) << " "  << rsd_j.name1() << " " << rsd_j.atom_name(jj ) << std::endl;
				}

			}

			// where does Mg(2+) lie in coordinate system of phosphate?
			if ( rsd_j.is_RNA()  ){
				Stub const & stub_j = all_phosphate_stubs[j];
				Vector const xyz_mg_local = stub_j.global2local( xyz_mg );
					Real const x = xyz_mg_local.x();
					Real const y = xyz_mg_local.y();
					Real const z = xyz_mg_local.z();
					if ( xyz_mg_local.length() < 12.0 ) {
						out4 << I(4,i) << ' ' << I(4,j) << ' ' << rsd_j.name1() << ' ' << F(8,3,x) << ' ' << F(8,3,y) << ' ' << F(8,3,z) << std::endl;
					}
			}


			// where does Mg(2+) lie in coordinate system of base?
			if ( rsd_j.is_RNA()  ){
				Stub const & stub_j = all_base_stubs[j];
				Vector const xyz_mg_local = stub_j.global2local( xyz_mg );
				Real const x = xyz_mg_local.x();
				Real const y = xyz_mg_local.y();
				Real const z = xyz_mg_local.z();
				if ( std::abs( z ) < 8.0 && ( x*x+y*y ) < 12.0*12.0 ){
					out5 << I(4,i) << ' ' << I(4,j) << ' ' << rsd_j.name1() << ' ' << F(8,3,x) << ' ' << F(8,3,y) << ' ' << F(8,3,z) << std::endl;
				}
			}

			// let's do hydrogen atoms too
			for ( Size jj = rsd_j.nheavyatoms() + 1; jj <= rsd_j.natoms(); jj++ ){

				if ( !rsd_j.atom_is_hydrogen( jj ) ) continue;

				Vector xyz_jj = rsd_j.xyz( jj );
				Real distance = (xyz_jj - xyz_mg).length();
				if ( distance < DIST_CUTOFF ) {

						//<< I(4, rsd_j.atom_type_set().atom_type_index(  rsd_j.atom_type( jj ).name() ) ) << std::endl;

					// let's also get the angles!
					Real theta( 0.0 );
					Vector const & xyz_base = rsd_j.xyz(  rsd_j.atom_base( jj ) );
					theta = numeric::conversions::degrees( angle_of( xyz_mg, xyz_jj, xyz_base ) );

					out6 << I(4, count) << " " <<  I(4, i) << " " << I(4, j) << " " << I(4, jj) << " " << F(8, 3, distance) << " "  << rsd_j.name1() << " " << rsd_j.atom_name(jj ) <<  " " <<  F(8,3,theta) << " " <<   rsd_j.atom_type( jj ).name() << std::endl;

				}
			}



		} // loop over RNA residues.

		// info on total number of close (;inner-sphere') ligands for each Mg(2+)
		out << I(4, count )	<< ' ' << pdb_file << ' ' << I(4,i) << ' ' << nligands	<< std::endl;

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
	utility::io::ozstream out3( "vdwinfo_"+outfile );
	utility::io::ozstream out4( "phosphate3d_"+outfile );
	utility::io::ozstream out5( "base3d_"+outfile );
	utility::io::ozstream out6( "Hinfo_"+outfile );

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

		//pose.dump_pdb( "imported_"+pdb_file );

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;

		mg_pdbstats_from_pose( out, out2, out3, out4, out5, out6, pose, count, total_residues, pdb_file );

		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	std::cout << "FINISHED ==> TOTAL RESIDUES Processed: " << total_residues << std::endl;

}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {


	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);

	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	mg_pdbstats_test();

	exit( 0 );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
