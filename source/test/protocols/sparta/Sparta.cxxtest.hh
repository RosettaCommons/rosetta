// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/sparta/Sparta.hh>
#include <protocols/sparta/ANN.hh>

#include <utility/vector1.hh>
#include <numeric>
#include <ObjexxFCL/string.functions.hh>
//Auto Headers
#include <utility/vector0.hh>


using namespace protocols::sparta;
using namespace ObjexxFCL;

class SpartaTest : public CxxTest::TestSuite {

private:
	Sparta::SpartaLib* lib_instance_;
	Sparta::SpartaLib& lib() { return *lib_instance_; }
	core::pose::Pose test_pose_;
public:

	SpartaTest() : lib_instance_( NULL ) {}

	void setUp() {
		if ( !lib_instance_ ) {
			core_init();
			lib_instance_ = new Sparta::SpartaLib;
			std::string const cs_file ( "protocols/sparta/data_16988.tab" );
			std::string const pdb_file( "protocols/sparta/2kywA.pdb" );
			core::import_pose::pose_from_file( test_pose_, pdb_file , core::import_pose::PDB_file);
			lib_instance_->setup_for_scoring( test_pose_ );
			lib_instance_->getResInfo( false );
		}
	}
	void tearDown() {};

	~SpartaTest() {
		delete lib_instance_;
	}

	void test_BLOSUM62() {
		core::Real const TOLERATED_ERROR( 5e-3 );
		//  for ( Sparta::SpartaLib::BlosumMatrix::iterator it1 = lib().BLOSUM_62.begin(); it1 != lib().BLOSUM_62.end(); ++it1 ) {
		//core::Size ct( 0 );
		//for ( Sparta::SpartaLib::BlosumMatrix::mapped_type::iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2, ++ct ) {
		// std::cout << it1->first << " " << ct << " " << *it2 << std::endl;
		// }
		//}
		//  TS_ASSERT_EQUAL( lib().BLOSUM_62["A"][4]
		#include "asserts_BLOSUM62.cc"
	}

	void test_PDB() {
		PDB& PDB=lib().inPDB;

		//   std::cout << "dumping the ATOMS entry of PDB " << std::endl;
		//   for ( PDB::AtomsMap::const_iterator it1 = PDB.ATOMS.begin(); it1!=PDB.ATOMS.end(); ++it1 ) {
		//    for ( PDB::EntryMap::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ) {
		//     for ( PDB::AtomEntries::const_iterator it3 = it2->second.begin(); it3 != it2->second.end(); ++it3 ) {
		//      std::cout << it1->first << " " << it2->first << " " << it3->first << " " << it3->second << std::endl;
		//     }
		//    }
		//   }
		TS_ASSERT_EQUALS( PDB.ATOMS[1][87]["HN"].atomNum, 1245 );
		TS_ASSERT_EQUALS( PDB.ATOMS[1][87]["HN"].resNum, 87 );
		TS_ASSERT_EQUALS( PDB.ATOMS[1][87]["HN"].resName, "HIS" );
		TS_ASSERT_EQUALS( PDB.ATOMS[1][32]["C"].resName, "SER" );
		TS_ASSERT_EQUALS( PDB.ATOMS[1][32]["C"].resNum, 32 );
		TS_ASSERT_EQUALS( PDB.ATOMS[1][32]["C"].atomNum, 435 );
	}

	//  void test_PDB_Hbonds() {
	//   PDB& inPDB=lib().inPDB;
	//   inPDB.initOrbitalShift();
	//   inPDB.initHBond();
	//   PDB::HBondMap const& hbd( inPDB.HBDistList );
	//   for ( PDB::HBondMap::const_iterator it1 = hbd.begin(); it1 != hbd.end(); ++it1 ) {
	//    for ( PDB::HBondMap::mapped_type::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ) {
	//     std::cout << "hbd-map: " << it1->first << " " << it2->first << " " << it2->second << std::endl;
	//    }
	//   }
	//   PDB::HBondMap const& dho_angles( inPDB.HB_DHO_AngleList );
	//   for ( PDB::HBondMap::const_iterator it1 = dho_angles.begin(); it1 != dho_angles.end(); ++it1 ) {
	//    for ( PDB::HBondMap::mapped_type::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ) {
	//     std::cout << "dho-map: " << it1->first << " " << it2->first << " " << it2->second << std::endl;
	//    }
	//   }
	//   PDB::HBondMap const& hoa_angles( inPDB.HB_HOA_AngleList );
	//   for ( PDB::HBondMap::const_iterator it1 = hoa_angles.begin(); it1 != hoa_angles.end(); ++it1 ) {
	//    for ( PDB::HBondMap::mapped_type::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); ++it2 ) {
	//     std::cout << "hoa-map: " << it1->first << " " << it2->first << " " << it2->second << std::endl;
	//    }
	//   }
	// }

	void test_ANN_PARAMETERS() {
		using namespace core;
		ANN ann;
		ann.init( 113,30,1,9,6,3,lib().TAB_DIR, "HN" );
		//typedef utility::vector0<float> WeightVector;
		core::Real const TOLERATED_ERROR( 5e-3 );
		//   Size ct = 0;
		//   for ( WeightVector::const_iterator it = ann.BI_1.begin(); it != ann.BI_1.end(); ++it, ++ct ) {
		//    std::cout << "BI_1 " << ct << " " << *it << std::endl;
		//   }
		//   ct = 0;
		//   for ( WeightVector::const_iterator it = ann.BI_2.begin(); it != ann.BI_2.end(); ++it, ++ct ) {
		//    std::cout << "BI_2 " << ct << " " << *it << std::endl;
		//   }
		//   ct = 0;
		//   for ( WeightVector::const_iterator it = ann.BI_3.begin(); it != ann.BI_3.end(); ++it, ++ct ) {
		//    std::cout << "BI_3 " << ct << " " << *it << std::endl;
		//   }
		//   ct = 0;
		//   for ( WeightVector::const_iterator it = ann.BI_1.begin(); it != ann.BI_1.end(); ++it, ++ct ) {
		//    std::cout << "BI_1 " << ct << " " << *it << std::endl;
		//   }
		//   typedef boost::unordered_map<int, utility::vector0<float> > WeightMap;
		//   for ( WeightMap::const_iterator it = ann.WI_1.begin(); it != ann.WI_1.end(); ++it ) {
		//    std::cout << "WI_1 " << it->first << " " << it->second[0] << " " <<it->second[1] << " " <<it->second[2] <<  std::endl;
		//   }
		//   for ( WeightMap::const_iterator it = ann.WI_2.begin(); it != ann.WI_2.end(); ++it ) {
		//    std::cout << "WI_2 " << it->first << " " << it->second[0] << " " <<it->second[1] << " " <<it->second[2] <<  std::endl;
		//   }
		//   for ( WeightMap::const_iterator it = ann.WI_3.begin(); it != ann.WI_3.end(); ++it ) {
		//    std::cout << "WI_3 " << it->first << " " << it->second[0] << " " <<it->second[1] << " " <<it->second[2] <<  std::endl;
		//   }
		#include "asserts_ANN_Parameters.cc"

	}
	void test_ANN_INPUT() {
		core::Real const TOLERATED_ERROR( 5e-3 );
		//typedef utility::vector0<float> Values;

		// for ( int resid = lib().r1+1; resid < lib().rN; resid++ ) {
		//    Values const& temp( lib().ANN_IN_MTX[resid] );
		//    core::Size ct( 0 );
		//    std::cout << "DUMP ANN INPUT for residue " << resid << std::endl;
		//    for ( Values::const_iterator it = temp.begin(); it != temp.end(); ++it, ++ct ) {
		//     std::cout << resid << " " << ct << " " << *it << std::endl;
		//    }
		//   }
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][3],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][4],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][5],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][6],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][7],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][8],   0.5000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][9],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][10],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][11],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][12],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][13],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][14],   0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][15],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][16],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][17],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][18],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][19],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][20],  -0.9527, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][21],  -0.3038, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][22],   0.1566, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][23],  -0.9877, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][24],  -0.9484, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][25],   0.3171, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][26],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][27],   0.0111, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][28],  -0.9999, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][29],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][30],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][31],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][32],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][33],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][34],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][35],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][36],   0.8000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][37],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][38],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][39],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][40],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][41],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][42],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][43],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][44],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][45],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][46],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][47],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][48],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][49],   0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][50],  -0.9999, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][51],   0.0150, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][52],   0.7799, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][53],  -0.6259, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][54],  -0.8967, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][55],   0.4426, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][56],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][57],   0.5093, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][58],  -0.8606, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][59],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][60],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][61],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][62],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][63],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][64],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][65],   0.6000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][66],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][67],  -0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][68],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][69],  -0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][70],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][71],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][72],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][73],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][74],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][75],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][76],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][77],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][78],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][79],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][80],  -0.9980, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][81],  -0.0632, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][82],  -0.8297, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][83],  -0.5581, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][84],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][85],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][86],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][87],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][88],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][89],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][90],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][91],   1.9609, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][92],  -0.7884, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][93],  -0.7316, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][94],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][95],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][96],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][97],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][98],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][99],   2.3910, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][100],  -0.8734, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][101],  -0.9260, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][102],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][103],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][104],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][105],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][106],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][107],   2.0971, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][108],  -0.9649, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][109],  -0.5843, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][110],   0.8926, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][111],   0.8765, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[78][112],   0.7871, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][0],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][1],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][2],   0.6000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][3],   0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][4],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][5],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][6],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][7],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][8],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][9],  -0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][10],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][11],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][12],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][13],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][14],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][15],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][16],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][17],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][18],  -0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][19],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][20],  -0.9529, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][21],  -0.3033, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][22],   0.1866, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][23],  -0.9824, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][24],   0.8494, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][25],   0.5278, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][26],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][27],  -0.9673, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][28],   0.2535, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][29],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][30],   0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][31],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][32],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][33],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][34],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][35],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][36],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][37],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][38],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][39],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][40],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][41],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][42],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][43],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][44],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][45],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][46],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][47],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][48],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][49],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][50],  -0.5287, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][51],  -0.8488, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][52],   0.8761, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][53],  -0.4820, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][54],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][55],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][56],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][57],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][58],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][59],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][60],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][61],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][62],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][63],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][64],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][65],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][66],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][67],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][68],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][69],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][70],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][71],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][72],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][73],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][74],  -0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][75],   0.1000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][76],   0.4000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][77],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][78],  -0.3000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][79],  -0.2000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][80],  -0.9625, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][81],  -0.2712, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][82],   0.8319, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][83],  -0.5550, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][84],   0.8756, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][85],   0.4830, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][86],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][87],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][88],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][89],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][90],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][91],   2.2897, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][92],  -0.8294, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][93],  -0.7737, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][94],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][95],   2.0099, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][96],  -0.8568, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][97],  -0.9907, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][98],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][99],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][100],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][101],   0.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][102],   1.0000, TOLERATED_ERROR );
		TS_ASSERT_DELTA( lib().ANN_IN_MTX[5][103],   2.2206, TOLERATED_ERROR );
	}

	void test_PRED_SUM() {
		core::Real const TOLERATED_ERROR( 5e-3 );
		//  lib().setup_for_scoring( test_pose_ );
		GDB Pred_Sum = lib().get_ANN_data( false );
		//dump all entries
		//  for ( Sparta::SpartaLib::AtomNameList::const_iterator it = lib().aN.begin(); it != lib().aN.end(); ++it ) {
		//     std::string const& atom_name( it->second );
		//     std::cout << " dump predicted shifts for " << atom_name << std::endl;
		//     for ( core::Size resid = 1; resid <= 87; ++resid ) {
		//      GDB_Entry temp = Pred_Sum.getEntry("RESID",string_of( resid ),"ATOMNAME",atom_name,1);
		//      float pred_shift = atof( (temp["SHIFT"]).c_str() );
		//      std::cout << " residue " << resid << " " << pred_shift
		//          << " " << temp["HM_SHIFT"]
		//          << " " << temp["EF_SHIFT"]
		//          << " " << temp["RC_SHIFT"]
		//          << " " << temp["SS_SHIFT"]
		//          << " " << temp["SIGMA"]
		//          << std::endl;
		//     }
		//   }
		#include "asserts_PRED_SUM.cc"
	}

	core::Real calc_sparta_score(
		std::string const & pdb_file,
		std::string const & cs_file
	) {
		using core::Size;
		using core::Real;
		using std::string;
		using utility::vector1;
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,pdb_file, core::import_pose::PDB_file);
		Sparta sparta(cs_file);

		Real const sparta_score(sparta.score_pose(pose));
		//std::cout << "score(" << pdb_file << "," << cs_file << ") = " << sparta_score << std::endl;
		return sparta_score;
	}


	void test_scores() {
		using core::Size;
		using core::Real;
		using std::string;
		using utility::vector1;
		using namespace protocols::sparta;
		string const cs_file ( "protocols/sparta/data_16988.tab" );
		core::Real const TOLERATED_ERROR( 5e-3 );

		Sparta sparta(cs_file);

		core::pose::Pose pose( test_pose_ );
		Real const sparta_score( sparta.score_pose(pose) );
		vector1< float > scores( sparta.score_pose_per_residue(pose) );
		float sum_scores( std::accumulate( scores.begin(), scores.end(), 0.0 ) );

		TS_ASSERT_DELTA( sparta_score, 180.3510, TOLERATED_ERROR );
		TS_ASSERT( scores.size() == 87 );
		TS_ASSERT_DELTA( sparta_score, sum_scores/4, TOLERATED_ERROR );
		#include "asserts_scores.cc"
	}

	// failing example from Firas
	void test_sparta_foldit() {
		using core::Size;
		using core::Real;
		using std::string;
		using utility::vector1;
		using namespace protocols::sparta;
		float const TOLERATED_ERROR( 5e-1 );

		string const cs_file ( "protocols/sparta/data_17280.tab" );
		string const pdb_file( "protocols/sparta/solution_0022372568.ir_solution.pdb" );
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,pdb_file, core::import_pose::PDB_file);
		Sparta sparta(cs_file);

		Real const sparta_score(sparta.score_pose(pose));
		TS_ASSERT_DELTA( sparta_score, 702.378, TOLERATED_ERROR );
		//std::cout << "score(" << pdb_file << "," << cs_file << ") = " << sparta_score << std::endl;
	}

};
