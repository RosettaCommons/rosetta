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
/// @author Chu Wang

// Unit Headers

// Rosetta Headers
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <basic/options/option.hh>
#include <protocols/loops/loops_main.hh>

// ObjexxFCL header
#include <ObjexxFCL/format.hh>

// numeric header
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

// utility header
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
// C++ Headers
#include <iostream>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>



int
main( int argc, char * argv [] )
{
	try {

	using namespace ObjexxFCL::format;

	// initialize option and random number system
	devel::init( argc, argv );
	// get input file names
	utility::vector1< utility::file::FileName > list_file_names;
	utility::vector1< utility::file::FileName > pdb_file_names;
	list_file_names = basic::options::option[ basic::options::OptionKeys::in::file::l ]().vector(); // make a copy (-l)
	for(utility::vector1< utility::file::FileName >::iterator i = list_file_names.begin(), i_end = list_file_names.end();
			i != i_end; ++i) {
		std::string filename( i->name() );
		utility::io::izstream data( filename.c_str() );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open file: " + filename + '\n' );
		}
		std::string line;
		while( getline(data, line) ) {
			pdb_file_names.push_back( utility::file::FileName(line) );
		}
		data.close();
	}
	std::string outname("zinc_stats.txt");
	if ( basic::options::option[ basic::options::OptionKeys::out::file::o ].user() ) {
		outname =basic::options::option[ basic::options::OptionKeys::out::file::o ]();
	}
	std::ofstream outfile(outname.c_str(), std::ios::out | std::ios::app);
	// for each input pdb, do zinc statiscs
	for(utility::vector1< utility::file::FileName >::iterator i = pdb_file_names.begin(), i_end = pdb_file_names.end();
			i != i_end; ++i) {
		std::string base_name = i->base();
		std::cerr << "starting " << i->name() << std::endl;
		core::pose::Pose pose;
		core::import_pose::pose_from_pdb( pose, i->name() );
		core::scoring::ScoreFunctionOP scorefxn = core::scoring::getScoreFunction();
		(*scorefxn)(pose); //scoring for tenA neighbor graph

		core::pose::PDBInfoOP pdb_info = pose.pdb_info();

		core::Size zinc_count(0);
		for ( core::Size j = 1; j <= pose.total_residue(); ++j ) {
			if ( pose.residue(j).name() == "ZN" ) {
				core::Vector const & zn_position( pose.residue(j).xyz("ZN") ); //coord of zinc
				utility::vector1<bool> zn_nbr( pose.total_residue(), false );
				zn_nbr[j] = true;
				protocols::loops::get_tenA_neighbor_residues( pose, zn_nbr );
				utility::vector1<core::Size> ligand_res;
				utility::vector1<core::Size> ligand_atm;
				for ( core::Size k = 1; k <= zn_nbr.size(); ++k ) {
					if ( ! zn_nbr[k]  || k==j ) continue; // not neighbor or zn itself
					core::conformation::Residue const & rsd = pose.residue(k);
					if ( ! rsd.is_protein() ) continue; // not protein residues
					core::Real best_dis = 10000.0;
					core::Size best_atom = 0;
					// go over all sidechain heavy atoms in this residue
					for ( core::Size kk=rsd.first_sidechain_atom(); kk<=rsd.nheavyatoms(); ++kk ) {
						core::Real dis = rsd.xyz(kk).distance( zn_position );
						if ( dis <= best_dis ) {
							best_dis = dis;
							best_atom = kk;
						}
					}
					// if at least atom within 3.0 A of zinc, add it to ligand residue list
					if (best_dis < 3.0) {
						ligand_res.push_back( k );
						ligand_atm.push_back( best_atom );
					}
				} // finish screening all zinc neighbors

				// consider only cases with four zinc ligand residues( structural zinc binding ).
				if ( ligand_res.size() != 4 )  continue;

				zinc_count++;
				core::Size zinc_pdb_res = pdb_info->number(j);
				char zinc_pdb_chain = pdb_info->chain(j);
				if ( zinc_pdb_chain == ' ' ) zinc_pdb_chain = '_';

				for ( core::Size n=1; n<= ligand_res.size(); ++n ) {
					core::Size zres = ligand_res[n];
					core::Size ligand_pdb_res = pdb_info->number(zres);
					char ligand_pdb_chain = pdb_info->chain(zres);
					if ( ligand_pdb_chain == ' ' ) ligand_pdb_chain = '_';

					core::conformation::Residue const & rsd = pose.residue(zres);
					core::Size zatm = ligand_atm[n];
					core::Size zatm_base = rsd.type().atom_base(zatm);
					core::Size zatm_base2 = rsd.type().atom_base(zatm_base);
					std::string cb_atom = ( rsd.name3() == "GLY" ) ? "CA" : "CB";
					core::Real distance = rsd.xyz(zatm).distance(zn_position);
					core::Real cb_zn_distance = rsd.xyz(cb_atom).distance(zn_position);
					core::Real angle = numeric::conversions::degrees( numeric::angle_radians(zn_position, rsd.xyz(zatm), rsd.xyz( zatm_base) ) );
					core::Real dihedral = numeric::conversions::degrees( numeric::dihedral_radians( zn_position, rsd.xyz(zatm), rsd.xyz(zatm_base), rsd.xyz(zatm_base2) ) ) ;
					outfile << A(4, base_name) << " " << I(3, zinc_count) << A(3, " ZN") << " " << I(4,zinc_pdb_res) << " " << A(1,zinc_pdb_chain) << " "
									<< A(3,rsd.name() ) << " " << I(4,ligand_pdb_res) << " " << A(1,ligand_pdb_chain) << " " << A(4,rsd.atom_name(zatm) ) << " "
									<< F(4,2,distance) <<  " " << F(7,2,angle) << " " << F(7,2,dihedral) << " " << F(7,2, cb_zn_distance) << std::endl;
				} // finish this zinc - residue pair
			} // finish this zinc binding site
		} // finiish all residues in this protein
		std::cerr << "finishing " << i->name() << std::endl;
	} // finish this pdb
	outfile.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

	return 0;
} // finish main

