// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author jk
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

#include <numeric/xyz.functions.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <utility/io/mpistream.hh>

// Eigen header for SVD support
#include <Eigen/SVD>

// Armadillo can be used, but this requires LAPACK and BLAS (not distributed with Rosetta)
//#include <armadillo/armadillo>

//Protocol Headers
#include <protocols/pockets/PocketGrid.hh>
//#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>


using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace core::pose::datacache;
using namespace core::id;


/*
void test_PCA_by_armadillo() {

using namespace std;
using namespace arma;

int const num_data_points = 3514;
std::ifstream infile;
infile.open("/Users/johnk/temp/data_for_princomp.txt");
if (!infile) {
std::cerr << "Unable to open input file" << std::endl;
exit(1);
}
std::cerr << "Reading input file" << std::endl;

// setup
mat A = randu<mat>(num_data_points,2);
double new_C, new_P;
double avg_C(0.), avg_P(0.);
int j=0;
while (infile >> new_C >> new_P ) {
A(j,0) = new_C;
A(j,1) = new_P;
avg_C += new_C;
avg_P += new_P;
j++;
}
infile.close();
avg_C /= num_data_points;
avg_P /= num_data_points;
for ( int j = 0; j < num_data_points; ++j ) {
A(j,0) -= avg_C;
A(j,1) -= avg_P;
}
// std::cout << "starting A" << std::endl;
// std::cout << A << std::endl;

// PCA code (armadillo)
mat coeff, score;
colvec latent;
princomp(coeff, score, latent, A);
std::cout << "coeff" << std::endl;
std::cout << coeff << std::endl;
// we get the desired result by computing A * coeff, where coeff is:
//   0.2514  -0.9679
//   0.9679   0.2514
// note: these are equivalent (and both are what we want):
// mat transformed = score;
// mat transformed = A * coeff;

// SVD code (armadillo)
mat U;
vec s;
mat V;
svd(U,s,V,A);
std::cout << "V is: " << std::endl;
std::cout << V << std::endl;
std::cout << "U is: " << U.n_cols << " by " << U.n_cols << std::endl;
std::cout << "V is: " << V.n_cols << " by " << V.n_cols << std::endl;
mat transformed = A * V;

// Write output and exit
utility::io::ozstream fout;
fout.open("transformed.txt", std::ios::out);
fout << transformed << std::endl;
fout.close();
fout.clear();

exit(1);
return;
}
*/


void test_PCA_eigen() {

	using namespace std;
	using namespace Eigen;

	int const num_data_points = 3514;

	std::ifstream infile;
	infile.open("/Users/johnk/temp/data_for_princomp.txt");
	if ( !infile ) {
		std::cerr << "Unable to open input file" << std::endl;
		exit(1);
	}
	std::cerr << "Reading input file" << std::endl;

	// setup
	MatrixXf A(num_data_points,2);
	double new_C, new_P;
	double avg_C(0.), avg_P(0.);
	int j=0;
	while ( infile >> new_C >> new_P ) {
		A(j,0) = new_C;
		A(j,1) = new_P;
		avg_C += new_C;
		avg_P += new_P;
		j++;
	}
	infile.close();
	avg_C /= num_data_points;
	avg_P /= num_data_points;
	for ( int j = 0; j < num_data_points; ++j ) {
		A(j,0) -= avg_C;
		A(j,1) -= avg_P;
	}

	// std::cout << "Matrix A is: " << std::endl << A << std::endl;

	JacobiSVD<MatrixXf> svd(A, ComputeFullU | ComputeFullV );
	svd.computeV();
	std::cout << "Matrix V is: " << std::endl << svd.matrixV() << std::endl;
	MatrixXf transformed = A * svd.matrixV();

	// Write output and exit
	utility::io::ozstream fout;
	fout.open("transformed.txt", std::ios::out);
	fout << transformed << std::endl;
	fout.close();
	fout.clear();

	exit(1);
	return;
}


OPT_KEY( String, template_pdb_name )
OPT_KEY( String, template_target_resnum )
OPT_KEY( String, comparison_target_resnum )

static THREAD_LOCAL basic::Tracer TR( "apps.pilot.karen_pocket_compare.main" );

/// General testing code
int main( int argc, char * argv [] ) {

	try{

		NEW_OPT( template_pdb_name, "template pdb", "template.pdb" );
		NEW_OPT( template_target_resnum, "template target residue", "-1" );
		NEW_OPT( comparison_target_resnum, "comparison target residue", "-1" );

		//initializes Rosetta functions
		devel::init(argc, argv);

		TR << "Starting pocket compare" << std::endl;

		//sets pdb as a Rosetta pose
		std::string const template_fname ( option[ template_pdb_name ] );
		pose::Pose template_pose;
		core::import_pose::pose_from_file( template_pose, template_fname , core::import_pose::PDB_file);
		TR << "set template pdb" << std::endl;

		// set the template target residue
		int template_target_residue_number;
		char chain_t = ' ';
		std::string const resid_t = option[template_target_resnum];
		std::size_t fpos_t( resid_t.find(':') );
		if ( fpos_t != std::string::npos ) {
			template_target_residue_number = ObjexxFCL::int_of( resid_t.substr(0,fpos_t) );
			if ( fpos_t != resid_t.size()-1 ) {
				chain_t = resid_t[ fpos_t+1 ];
			}
		} else {
			template_target_residue_number = ObjexxFCL::int_of( resid_t );
		}
		int seqpos_t = 0;
		for ( int j = 1, resnum = template_pose.total_residue(); j <= resnum; ++j ) {
			if ( template_pose.pdb_info()->number(j) == template_target_residue_number ) {
				seqpos_t = j;
				if ( chain_t != ' ' ) {
					if ( template_pose.pdb_info()->chain(j) == chain_t ) {
						seqpos_t = j;
					}
				} else {
					seqpos_t = j;
				}
			}
		}
		if ( seqpos_t == 0 ) {
			std::cout << "ERROR!! Invalid template target residue" << std::endl;
			exit(1);
		}
		TR << "Set template target residue " << resid_t << " (seqpos_t " << seqpos_t << ")" << std::endl;


		//call function to make a grid around a target residue (seqpos_t)
		protocols::pockets::PocketGrid template_pg( template_pose.conformation().residue(seqpos_t) );
		// note: this realigns the pose!!
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.move_pose_to_standard_orie( seqpos_t, template_pose );
		template_pg.autoexpanding_pocket_eval( template_pose.conformation().residue(seqpos_t), template_pose );
		template_pose.dump_pdb( "template_aligned_pose.pdb" );
		template_pg.dumpGridToFile( "template_pocket.pdb" );

		// create pose for comparison pose from pdb
		std::string const comparison_pdb_name ( basic::options::start_file() );
		pose::Pose comparison_pose;
		core::import_pose::pose_from_file( comparison_pose, comparison_pdb_name , core::import_pose::PDB_file);
		TR << "set comparison pdb" << std::endl;

		// set the comparison target residue
		int comparison_target_residue_number;
		char chain_c = ' ';
		std::string const resid_c = option[comparison_target_resnum];
		std::size_t fpos_c( resid_c.find(':') );
		if ( fpos_c != std::string::npos ) {
			comparison_target_residue_number = ObjexxFCL::int_of( resid_c.substr(0,fpos_c) );
			if ( fpos_c != resid_c.size()-1 ) {
				chain_c = resid_c[ fpos_c+1 ];
			}
		} else {
			comparison_target_residue_number = ObjexxFCL::int_of( resid_c );
		}
		int seqpos_c = 0;
		for ( int j = 1, resnum = comparison_pose.total_residue(); j <= resnum; ++j ) {
			if ( comparison_pose.pdb_info()->number(j) == comparison_target_residue_number ) {
				seqpos_c = j;
				if ( chain_c != ' ' ) {
					if ( comparison_pose.pdb_info()->chain(j) == chain_c ) {
						seqpos_c = j;
					}
				} else {
					seqpos_c = j;
				}
			}
		}
		if ( seqpos_c == 0 ) {
			std::cout << "ERROR!! Invalid comparison target residue" << std::endl;
			exit(1);
		}
		TR << "Set comparison target residue " << resid_c << " (seqpos_c " << seqpos_c << ")" << std::endl;

		//call function to make a grid around a target residue (seqpos_c)
		protocols::pockets::PocketGrid comparison_pg( comparison_pose.conformation().residue(seqpos_c) );
		// note: this realigns the pose!!
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.move_pose_to_standard_orie( seqpos_c, comparison_pose );
		comparison_pg.autoexpanding_pocket_eval( comparison_pose.conformation().residue(seqpos_c), comparison_pose );

		core::Real best_dist = template_pg.get_pocket_distance( comparison_pg );
		comparison_pose.dump_pdb( "comparison_aligned_pose.pdb" );
		comparison_pg.dumpGridToFile( "comparison_pocket.pdb" );

		// from PCA, the best match could happen if we rotate pose 180 degrees about the x- or y-axis (or both)
		numeric::xyzVector<core::Real> zero_vec;
		zero_vec.zero();

		// rotate about y-axis, recompute distance
		{
			numeric::xyzMatrix<core::Real> rot_mat = numeric::y_rotation_matrix_degrees(180.);
			comparison_pose.apply_transform_Rx_plus_v(rot_mat, zero_vec);
			protocols::pockets::PocketGrid tmp_pg( comparison_pose.conformation().residue(seqpos_c) );
			tmp_pg.autoexpanding_pocket_eval( comparison_pose.conformation().residue(seqpos_c), comparison_pose );
			core::Real tmp_dist = template_pg.get_pocket_distance( tmp_pg );
			if ( tmp_dist < best_dist ) {
				best_dist = tmp_dist;
				comparison_pose.dump_pdb( "comparison_aligned_pose.pdb" );
				tmp_pg.dumpGridToFile( "comparison_pocket.pdb" );
			}
		}

		// rotate about x-axis, recompute distance
		{
			numeric::xyzMatrix<core::Real> rot_mat = numeric::x_rotation_matrix_degrees(180.);
			comparison_pose.apply_transform_Rx_plus_v(rot_mat, zero_vec);
			protocols::pockets::PocketGrid tmp_pg( comparison_pose.conformation().residue(seqpos_c) );
			tmp_pg.autoexpanding_pocket_eval( comparison_pose.conformation().residue(seqpos_c), comparison_pose );
			core::Real tmp_dist = template_pg.get_pocket_distance( tmp_pg );
			if ( tmp_dist < best_dist ) {
				best_dist = tmp_dist;
				comparison_pose.dump_pdb( "comparison_aligned_pose.pdb" );
				tmp_pg.dumpGridToFile( "comparison_pocket.pdb" );
			}
		}

		// rotate about y-axis, recompute distance
		{
			numeric::xyzMatrix<core::Real> rot_mat = numeric::y_rotation_matrix_degrees(180.);
			comparison_pose.apply_transform_Rx_plus_v(rot_mat, zero_vec);
			protocols::pockets::PocketGrid tmp_pg( comparison_pose.conformation().residue(seqpos_c) );
			tmp_pg.autoexpanding_pocket_eval( comparison_pose.conformation().residue(seqpos_c), comparison_pose );
			core::Real tmp_dist = template_pg.get_pocket_distance( tmp_pg );
			if ( tmp_dist < best_dist ) {
				best_dist = tmp_dist;
				comparison_pose.dump_pdb( "comparison_aligned_pose.pdb" );
				tmp_pg.dumpGridToFile( "comparison_pocket.pdb" );
			}
		}

		// report best distance
		TR << "Distance: " << best_dist << std::endl;

		TR << "Done!" << std::endl;
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}


