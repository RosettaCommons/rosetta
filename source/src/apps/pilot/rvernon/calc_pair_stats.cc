// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file calc_pair_stats.cc
/// @brief
/// @author Robert Vernon

#include <core/types.hh>
#include <devel/init.hh>

#include <core/chemical/AA.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <core/pose/Pose.hh>
// Auto-header: duplicate removed #include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>


#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/robert.OptionKeys.gen.hh>

#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/FileName.hh>

#include <ObjexxFCL/char.functions.hh>
#include <ObjexxFCL/string.functions.hh>

#include <numeric/random/random.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/jumping/StrandPairing.hh>
#include <ObjexxFCL/format.hh>


///////////////////////////////////////////////////////////////////////////////

using std::string;
using core::Size;
using core::Real;
using utility::vector1;
using utility::file::FileName;
using ObjexxFCL::format::A;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
//using namespace basic;
//using namespace core::sequence;
//using namespace basic::options;
//using namespace basic::options::OptionKeys;

using namespace core;
using namespace pose;
using namespace conformation;

int
main( int argc, char* argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	// options, random initialization
	devel::init( argc, argv );

	std::string const pdb_file_list( option[ robert::pairdata_input_pdb_list ]() );
	//std::string const pdb_file_list("list");

	utility::io::izstream infile(pdb_file_list);
	utility::io::ozstream outfile(pdb_file_list+".data");

	std::string pdb_file_location;
	infile >> pdb_file_location;

	bool iterate = true;
	while (iterate) {

		pose::Pose pdb;

		std::cout << "PROCESSING PDB: " << pdb_file_location << std::endl;
		core::import_pose::centroid_pose_from_pdb( pdb, pdb_file_location );
		std::cout << "PROCESSING COMPLETE: " << pdb_file_location << std::endl;

		utility::vector1< pose::PoseOP > pdb_chains(pdb.split_by_chain());


		for ( utility::vector1< core::pose::PoseOP >::iterator it = pdb_chains.begin(), end = pdb_chains.end(); it != end; ++it ) {
			core::scoring::dssp::Dssp dssp(**it);
			dssp.insert_ss_into_pose( **it );


			//(*it)->dump_pdb("HEYO.pdb");

			core::Size n1(1);
			//for ( ResidueOPs::iterator res1( (*it)->res_begin() ); res1 != (*it)->res_end(); ++res1, ++n1 ) {
			for ( Size n1 = 1; n1 <= (*it)->total_residue(); ++n1 ) {
				core::conformation::Residue const & res1( (*it)->residue( n1 ) );
				core::Size n2(1);
				//for ( ResidueOPs::iterator res2( (*it)->res_begin() ); res2 != (*it)->res_end(); ++res2, ++n2 ) {
				for ( Size n2 = 1; n2 <= (*it)->total_residue(); ++n2 ) {
					core::conformation::Residue const & res2( (*it)->residue( n2 ) );
					if ( res1.is_protein() && res2.is_protein()) {


						if ( (n1 != n2) and (res1.xyz("CEN").distance( res2.xyz("CEN")) < 20.0 ) ) {

							Size d2h1(999), d2l1(999), d2e1(999);
							for ( Size ss_dist(1); ss_dist <= (*it)->total_residue(); ss_dist++ ) {
								if (( std::abs(ss_dist - n1) < d2h1) and ((*it)->secstruct(ss_dist) == 'H')) d2h1 = ( std::abs(ss_dist - n1));
								if (( std::abs(ss_dist - n1) < d2l1) and ((*it)->secstruct(ss_dist) == 'L')) d2l1 = ( std::abs(ss_dist - n1));
								if (( std::abs(ss_dist - n1) < d2e1) and ((*it)->secstruct(ss_dist) == 'E')) d2e1 = ( std::abs(ss_dist - n1));
							}

							Size d2h2(999), d2l2(999), d2e2(999);
							for ( Size ss_dist(1); ss_dist <= (*it)->total_residue(); ss_dist++ ) {
								if (( std::abs(ss_dist - n2) < d2h2) and ((*it)->secstruct(ss_dist) == 'H')) d2h2 = ( std::abs(ss_dist - n2));
								if (( std::abs(ss_dist - n2) < d2l2) and ((*it)->secstruct(ss_dist) == 'L')) d2l2 = ( std::abs(ss_dist - n2));
								if (( std::abs(ss_dist - n2) < d2e2) and ((*it)->secstruct(ss_dist) == 'E')) d2e2 = ( std::abs(ss_dist - n2));
							}


							//if ((res1.name3() != "GLY") and (res2.name3() != "GLY")) {

							//float dot_ca_cb = dot_product( (res1.xyz("CA")-res1.xyz("CB")).normalize(), (res2.xyz("CA")-res2.xyz("CB")).normalize() );
							float dot_ca_cen = dot_product( (res1.xyz("CA")-res1.xyz("CEN")).normalize(), (res2.xyz("CA")-res2.xyz("CEN")).normalize() );
							float dot_ca_c = dot_product( (res1.xyz("CA")-res1.xyz("C")).normalize(), (res2.xyz("CA")-res2.xyz("C")).normalize() );


								//outfile << pdb_file_location << "." << n1 << "." << n2 << " " << res1.aa() << " "  << d2l1 << " " << d2h1 << " " << d2e1 << " " << res2.aa() << " " << d2l2 << " " << d2h2 << " " << d2e2 << " " << (*it)->secstruct(n1) << " " << (*it)->secstruct(n2) << " " << std::abs(n2 - n1) << " " << (*it)->total_residue() << " " << res1.xyz("CA").distance( res2.xyz("CA") ) << " " << res1.xyz("CB").distance( res2.xyz("CB") ) << " " << dot_product(res1.xyz("CA")-res1.xyz("CB").normalize(), res2.xyz("CA")-res2.xyz("CB").normalize() ) << " " << angle_of(res1.xyz("CA")-res1.xyz("CB"), res2.xyz("CA")-res2.xyz("CB") ) << " " << dot_product(res1.xyz("C")-res1.xyz("CA"), res2.xyz("C")-res2.xyz("CA") ) << " " << angle_of(res1.xyz("C")-res1.xyz("CA"), res2.xyz("C")-res2.xyz("CA") ) << " " << angle_of( res1.xyz("CB"), res1.xyz("CA"), res2.xyz("CA") ) << " " << angle_of( res2.xyz("CB"), res2.xyz("CA"), res1.xyz("CA") ) << std::endl;

								//outfile << pdb_file_location << "." << n1 << "." << n2 << " " << res1.aa() << " "  << d2l1 << " " << d2h1 << " " << d2e1 << " " << res2.aa() << " " << d2l2 << " " << d2h2 << " " << d2e2 << " " << (*it)->secstruct(n1) << " " << (*it)->secstruct(n2) << " " << std::abs(n2 - n1) << " " << (*it)->total_residue() << " " << res1.xyz("CA").distance( res2.xyz("CA") ) << " " << res1.xyz("CB").distance( res2.xyz("CB") ) << " " << dot_ca_cb << " " << angle_of(res1.xyz("CA")-res1.xyz("CB"), res2.xyz("CA")-res2.xyz("CB") ) << " " << dot_ca_c << " " << angle_of(res1.xyz("CA")-res1.xyz("C"), res2.xyz("CA")-res2.xyz("C") ) << " " << angle_of( res1.xyz("CB"), res1.xyz("CA"), res2.xyz("CA") ) << " " << angle_of( res2.xyz("CB"), res2.xyz("CA"), res1.xyz("CA") ) << std::endl;


							outfile << pdb_file_location << "." << n1 << "." << n2 << " " << res1.aa() << " " << res2.aa() << " " << (*it)->secstruct(n1) << " " << (*it)->secstruct(n2) << " " << std::abs(n2 - n1) << " " << res1.xyz("CEN").distance( res2.xyz("CEN") ) << " " << dot_ca_cen << " " << dot_ca_c << std::endl;

								//} else {


								//float dot_ca_cen = dot_product( (res1.xyz("CA")-res1.xyz("CEN")).normalize(), (res2.xyz("CA")-res2.xyz("CEN")).normalize() );
								//float dot_ca_c = dot_product( (res1.xyz("CA")-res1.xyz("C")).normalize(), (res2.xyz("CA")-res2.xyz("C")).normalize() );
								//float dot_ca_cb = dot_product( (res1.xyz("CA")-res1.xyz("CB")).normalize(), (res2.xyz("CA")-res2.xyz("CB")).normalize() );

								//outfile << pdb_file_location << "." << n1 << "." << n2 << " "  << res1.aa() << " "  << d2l1 << " " << d2h1 << " " << d2e1 << " " << res2.aa() << " " << d2l2 << " " << d2h2 << " " << d2e2 << " " << (*it)->secstruct(n1) << " " << (*it)->secstruct(n2) << " " << std::abs(n2 - n1) << " " << (*it)->total_residue() << " " << res1.xyz("CA").distance( res2.xyz("CA") ) << " " << dot_ca_cb << " " << dot_ca_c << " " << dot_ca_cen << std::endl;
								//}
						}

					}
				}
			}
		}


		std::cout << "DONE!" << std::endl;

		pdb_file_location = "";
		infile >> pdb_file_location;
		if (pdb_file_location == "") iterate = false;
	}

		//std::string name( option[ in::file::vall ]() );
	//FileName fn( name );
	//VallReader vall( name );
} // int main( int argc, char * argv [] )
