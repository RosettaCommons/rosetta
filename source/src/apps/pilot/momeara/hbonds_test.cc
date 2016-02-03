// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/momeara/hbonds_test.cc
/// @brief this outputs all relavent information about hbonds for one or many pdbs.
/// @author Chris Sheldahl
/// @author Matthew O'Meara


// core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <devel/init.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
//#include <protocols/relax_protocols.hh>
#include <basic/Tracer.hh>
#include <core/scoring/hbonds/HBEvalTuple.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/hbonds/hbonds_geom.hh>
#include <core/scoring/hbonds/types.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AtomType.hh>


// utility headers
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>

// c++ headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/hbonds/HBondDatabase.hh>
#include <core/scoring/hbonds/HBondOptions.hh>


// namespaces
using namespace core;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace hbonds;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::T;
using basic::Error;
using basic::Warning;
using utility::file::FileName;


void dump_hbonds( std::string pdb_filename );


int
main( int argc, char * argv [] )
{
	try{
		// initialize
		devel::init(argc, argv);

		// concatenate -s and -l flags together to get total list of PDB files
		// (This was taken from Ian's early job distributor, thanks Ian)
		std::vector< FileName > pdb_file_names;
		if ( option[ in::file::s ].active() ) {
			pdb_file_names = option[ in::file::s ]().vector(); // make a copy (-s)
		}
		std::vector< FileName > list_file_names;
		if ( option[ in::file::l ].active() ) {
			list_file_names = option[ in::file::l ]().vector(); // make a copy (-l)
		}

		for ( std::vector< FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end(); i != i_end; ++i ) {
			std::string filename( i->name() );
			std::ifstream data( filename.c_str() );
			if ( !data.good() ) {
				utility_exit_with_message( "Unable to open file: " + filename + '\n' );
			}
			std::string line;
			while ( getline(data, line) ) {
				pdb_file_names.push_back( FileName(line) );
			}
			data.close();
		}


		// run dump_hbonds for each name in list
		for ( std::vector< FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end(); i != i_end; ++i ) {
			dump_hbonds( i->name() );
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void dump_hbonds( std::string pdb_filename )
{

	std::string out_filename = pdb_filename + "_hbonds.txt";
	std::ofstream fout( out_filename.c_str() );
	bool allowNonProtein = true;

	// read in pose
	Pose pose;
	core::import_pose::pose_from_file( pose, pdb_filename , core::import_pose::PDB_file);
	HBondDatabaseCOP hb_database( HBondDatabase::get_database());
	HBondOptions hboptions;
	HBondSet set1;

	ScoreFunctionOP scfxn( get_score_function() );
	(*scfxn)(pose);


	pose.update_residue_neighbors();
	fill_hbond_set( pose, false, set1, false );
	//note that the list of HBonds is indexed starting at 1
	for ( Size i = 1; i<= set1.nhbonds(); i++ ) {
		HBond bond = set1.hbond( i );
		//need to access donor and acc Residues as well as ints
		int hatm = bond.don_hatm();
		//int acc_atm = bond.acc_atm();
		HBEvalType type = bond.eval_type();
		int accResNum = bond.acc_res();
		int donResNum = bond.don_res();
		//get acc and donor residues from sequence numbers
		Residue accRes = pose.residue( accResNum );
		Residue donRes = pose.residue( donResNum );
		int const datm( donRes.atom_base( hatm ) );
		Vector const & hatm_xyz( donRes.atom( hatm ).xyz() );
		//donRes.atom() is of type Atom
		Vector const & datm_xyz( donRes.atom( datm ).xyz() );
		int aatm = bond.acc_atm();
		int const base( accRes.atom_base( aatm ) );
		int const base_of_base( accRes.atom_base( base ) );

		int const base2( accRes.abase2( aatm ) );
		Vector const & aatm_xyz = accRes.atom( aatm ).xyz();
		Vector const & base_xyz = accRes.atom( base ).xyz();
		Vector const & base_of_base_xyz = accRes.atom( base_of_base ).xyz();

		Vector const & base2_xyz = accRes.atom( base2 ).xyz();
		int donType = donRes.atom( datm ).type();
		int accType = accRes.atom( aatm ).type();
		std::string donResName = donRes.name3();
		std::string accResName = accRes.name3();
		bool isProtein = true;
		if ( ! donRes.is_protein() || ! accRes.is_protein() ) {
			isProtein = false;
		}

		const AtomType & donAtomType = donRes.atom_type( datm );
		const AtomType & accAtomType = accRes.atom_type( aatm );
		const std::string donElement = donAtomType.element();
		const std::string accElement = accAtomType.element();
		const std::string & donAtomName = donRes.atom_name( datm );
		const std::string & accAtomName = accRes.atom_name( aatm );
		int donElemNum, accElemNum; //usual atomic numbers from periodic table
		if ( ! donElement.compare( "O" ) ) {
			donElemNum = 8;
		} else if ( ! donElement.compare("N") ) {
			donElemNum = 7;
		} else if ( ! donElement.compare("S") ) {
			donElemNum = 16;
		} else {
			donElemNum = -1;
		}

		if ( ! accElement.compare( "O" ) ) {
			accElemNum = 8;
		} else if ( ! accElement.compare("N") ) {
			accElemNum = 7;
		} else if ( ! donElement.compare("S") ) {
			accElemNum = 16;
		} else {
			accElemNum = -1;
		}

		/*
		std::string donor_back;
		if ( bond.don_hatm_is_protein_backbone() )
		{
		donor_back = "donBK";
		}
		else
		{
		donor_back = "donSC";
		}
		std::string acc_back;
		if ( bond.acc_atm_is_protein_backbone() )
		{
		acc_back = "accBK";
		}
		else
		{
		acc_back = "accSC";
		}
		*/

		if ( ( set1.allow_hbond( i ) && isProtein ) || ( set1.allow_hbond( i ) && allowNonProtein ) ) {  //make sure it's not excluded
			fout << type << "\t";
			fout << datm_xyz.x()         << "\t" << datm_xyz.y()         << "\t" << datm_xyz.z()         << "\t";
			fout << hatm_xyz.x()         << "\t" << hatm_xyz.y()         << "\t" << hatm_xyz.z()         << "\t";
			fout << aatm_xyz.x()         << "\t" << aatm_xyz.y()         << "\t" << aatm_xyz.z()         << "\t";
			fout << base_xyz.x()         << "\t" << base_xyz.y()         << "\t" << base_xyz.z()         << "\t";
			fout << base_of_base_xyz.x() << "\t" << base_of_base_xyz.y() << "\t" << base_of_base_xyz.z() << "\t";
			fout << base2_xyz.x()        << "\t" << base2_xyz.y()        << "\t" << base2_xyz.z()        << "\t";
			//fout << donResName << "\t" <<  donResNum << "\t" << donType << "\t" << donor_back << "\t" << donAtomName << "\t" << donElemNum <<"\t";
			fout << donResName << "\t" <<  donResNum << "\t" << donType << "\t" << donAtomName << "\t" << donElemNum << "\t";
			//fout << accResName << "\t" << accResNum << "\t" << accType << "\t" << acc_back << "\t" << accAtomName << "\t" << accElemNum << "\t";
			fout << accResName << "\t" << accResNum << "\t" << accType << "\t" << accAtomName << "\t" << accElemNum << "\t";
			//Deriv deriv = bond.deriv();
			HBondDerivs deriv;
			Real energy = bond.energy();
			hb_energy_deriv(*hb_database, hboptions, HBEvalTuple( datm, donRes, aatm, accRes ), datm_xyz, hatm_xyz, aatm_xyz, base_xyz, base2_xyz, energy, true, deriv );
			Real weight = bond.weight();
			fout << energy << "\t" << weight << "\t";
			fout << deriv.h_deriv.f1()[ 0 ] << "\t" << deriv.h_deriv.f1()[ 1] << "\t" << deriv.h_deriv.f1()[ 2] << "\t";
			fout << deriv.h_deriv.f2()[ 0 ] << "\t" << deriv.h_deriv.f2()[ 1] << "\t" << deriv.h_deriv.f2()[ 2];
			fout << "\t" << pdb_filename;
			fout << "\n";
		} //close if ( set1.allow_hbond( i ) )

	} //close for
	fout.close();
}


