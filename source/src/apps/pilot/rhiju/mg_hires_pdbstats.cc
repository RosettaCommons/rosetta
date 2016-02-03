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
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/import_pose/import_pose.hh>
#include <core/id/AtomID.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <protocols/magnesium/util.hh>
#include <protocols/magnesium/MgOrbitalFrameFinder.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/string.functions.hh>

//silly using/typedef
// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>


using namespace core;
using namespace ObjexxFCL::format;
using namespace basic::options::OptionKeys;

using utility::vector1;
using io::pdb::old_dump_pdb;


typedef  numeric::xyzMatrix< Real > Matrix;

OPT_KEY( Boolean, dump )

///////////////////////////////////////////////////////////////////////////////
//
// Revisiting Mg(2+), including waters for hexahydrates -- output PDB stats
//  to help define a potential. Now use virtual atoms to help track the 'ligand field'
//  which favors octahedral coordination of lone pair donors.
//
//  Rhiju, April 2015
//
///////////////////////////////////////////////////////////////////////////////


///////////////////////////////////////////////////////////////////////////////
Size
mg_hires_pdbstats_from_pose( utility::io::ozstream & out,
pose::Pose & pose, Size const &,
Size & total_residues, std::string const pdb_file)
{

	using namespace core::conformation;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace core::scoring::hbonds;
	using namespace core::chemical::rna;
	using namespace protocols::magnesium;

	Size const nres = pose.total_residue();
	core::pose::PDBInfoCOP pdb_info = pose.pdb_info();
	core::chemical::AtomTypeSet const & atom_type_set = pose.residue_type( 1 ).atom_type_set();

	Size res_count( 0 );
	HBondOptions hbond_options; // needed to derive Mg(2+)-acceptor-base angles.

	MgOrbitalFrameFinder mg_orbital_frame_finder;
	mg_orbital_frame_finder.apply( pose );

	// find all the magnesiums in here.
	for ( Size i = 1; i <= nres; i++ ) {

		Residue const & rsd_i = pose.residue( i );
		if ( rsd_i.name3() != " MG" ) continue;
		Vector xyz_mg = rsd_i.xyz( "MG  " );

		res_count++;

		clock_t start_time = clock();
		utility::vector1< AtomID > const ligands = get_mg_ligands( pose, i, false /* will include everything, not just acceptors*/ );
		clock_t end_time = clock();
		std::cout << "ligand finder finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

		// to determine orbital frame, should only use ligands that really are acceptors.
		utility::vector1< AtomID > const acceptor_ligands = filter_acceptor_ligands( pose, ligands );

		//  determine_mg_orbital_frame( pose, i, acceptor_ligands );

		// output stats...
		for ( Size k = 1; k <= ligands.size(); k++ ) {

			Size const jj = ligands[ k ].atomno();
			Size const j  = ligands[ k ].rsd();
			Residue const & rsd_j = pose.residue( j );

			// measure angle at closest ligand-field atom  (ligand -- VX -- Mg )
			Vector const xyz_ligand = rsd_j.xyz( jj );
			utility::vector1< std::pair< Distance, std::string > > V_distances;
			for ( Size q = 1; q <= 6; q++ ) {
				std::string const V_name = "V"+I(1,q);
				V_distances.push_back( std::make_pair( ( xyz_ligand - rsd_i.xyz( V_name ) ).length(), V_name ) );
			}
			std::sort( V_distances.begin(), V_distances.end() );
			std::string const V_name = V_distances[ 1 ].second /* e.g., "V3" */;
			Vector xyz_V = rsd_i.xyz( V_name );
			Real const theta = numeric::conversions::degrees( angle_of( xyz_V, xyz_mg,  xyz_ligand ) );
			Distance const distance = ( xyz_ligand - xyz_mg ).length();

			// measure angle to base atom of acceptor (should be 120 degrees for sp2)
			Real theta_acceptor( 0.0 );
			if ( rsd_j.heavyatom_is_an_acceptor( jj ) ) {
				Vector dummy, xyz_base;
				chemical::Hybridization acc_hybrid( rsd_j.atom_type( jj ).hybridization());
				make_hbBasetoAcc_unitvector(
					hbond_options,
					acc_hybrid,
					rsd_j.atom( jj ).xyz(),
					rsd_j.xyz( rsd_j.atom_base( jj ) ),
					rsd_j.xyz( rsd_j.abase2( jj ) ),
					xyz_base, dummy );
				theta_acceptor = numeric::conversions::degrees( angle_of( xyz_mg, xyz_ligand, xyz_base ) );
			}
			if ( theta_acceptor == 0.0 ) {
				theta_acceptor = numeric::conversions::degrees( angle_of( xyz_mg, xyz_ligand, rsd_j.xyz( rsd_j.atom_base( jj ) ) ) );
				//theta_acceptor = numeric::conversions::degrees( angle_of( xyz_mg, xyz_ligand, rsd_j.xyz( rsd_j.bonded_neighbor( jj )[ 1 ] ) ) );
			}


			out << I(4, res_count) << " "
				<< I(4, pdb_info->number( i ) ) << " " << rsd_i.name3()
				<< " " << I( 2, ligands.size() ) // coordination number
				<< " " << I( 2, acceptor_ligands.size() ) // coordination number, filtered for acceptors
				<< " " << I( 2, k )
				<< " " << V_name
				<< " " << I(4, pdb_info->number( j ) ) << " " <<  rsd_j.name3()
				<< " " << A(4, rsd_j.atom_name( jj ) )
				<< " " << A(4, rsd_j.atom_type( jj ).name() )
				<< " " << I(4, atom_type_set.atom_type_index( rsd_j.atom_type( jj ).name() ) )
				<< " " << F(8, 3, distance  )
				<< " " << F(8, 3, theta )
				<< " " << F(8, 3, theta_acceptor )
				<< " " << pdb_file
				<<  std::endl;
		}

	} // i




	std::cout << "Processed " << res_count <<  std::endl;

	total_residues += res_count;
	return res_count;
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
	if ( !instream ) {
		std::cerr  << "Can't find list file " << pdb_list << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
		return;
	}

	std::string outfile  = option[ out::file::o ];
	utility::io::ozstream out( outfile ), out_mg_files( utility::replace_in( pdb_list, ".list", ".just_mg.list" ) );

	std::string pdb_file;
	Size count( 0 );

	pose::Pose pose;

	Size total_residues( 0 );

	ResidueTypeSetCOP rsd_set( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD ) );

	std::string line;
	while (  getline( instream, line )  ) {
		std::istringstream line_stream( line );

		line_stream >> pdb_file;

		clock_t start_time = clock();
		import_pose::pose_from_file( pose, *rsd_set,  file_path + '/' + pdb_file , core::import_pose::PDB_file);
		clock_t end_time = clock();
		std::cout << "pdb readin finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

		count++;
		std::cout << "Doing input file " << I(4,count) << " ==> " << pdb_file << std::endl;

		Size const n_mg = mg_hires_pdbstats_from_pose( out, pose, count, total_residues, pdb_file );

		if ( option[ dump ]() ) pose.dump_pdb( "imported_"+pdb_file );
		if ( n_mg > 0 ) out_mg_files << pdb_file << std::endl;

		std::cout << "TOTAL RESIDUES Processed: " << total_residues << std::endl;

	}

	out_mg_files.close();
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

		NEW_OPT( dump, "dump", false );

		devel::init(argc, argv);

		//   -extra_res ~/rosetta_database/chemical/residue_type_sets/fa_standard/residue_types/water/TP3.params

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
