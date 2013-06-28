// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Ragul Gowthaman

//GPU enabling is not default
//To test how many threads are fastest for your computer,
//use -gpu:threads 1024 (or other number) on the command line

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <devel/init.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
#include <basic/options/option_macros.hh>

// Utility Headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <core/conformation/Conformation.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <utility/options/StringOption.hh>

using namespace core;
using namespace basic::options;
using namespace std;
using namespace core::scoring;
using namespace core::optimization;
using namespace basic::options::OptionKeys;

OPT_KEY( String, protein )
OPT_KEY( Boolean, adt_grid )
OPT_KEY( Boolean, lig_grid )
//OPT_KEY( Real, gc_x )
//OPT_KEY( Real, gc_y )
//OPT_KEY( Real, gc_z )
//OPT_KEY( Real, gd_x )
//OPT_KEY( Real, gd_y )
//OPT_KEY( Real, gd_z )
OPT_KEY( Real, gs )
OPT_KEY( String, bound_ligand )
OPT_KEY( Integer, add_grid_size )
//OPT_KEY( Integer, add_ext_grid_size )

int main( int argc, char * argv [] ) {
	try{

  NEW_OPT( protein, "protein file name", "protein.pdb" );
	NEW_OPT( adt_grid, "resize grid based on user entered AUTODOCK grid values", false );
	NEW_OPT( lig_grid, "resize grid based on bound ligand", false );
	//  NEW_OPT( gc_x, "gid center : X ", 1.0 );
	//  NEW_OPT( gc_y, "gid center : Y ", 1.0 );
	//  NEW_OPT( gc_z, "gid center : Z ", 1.0 );
	//  NEW_OPT( gd_x, "gid dimension : X ", 20.0 );
	//  NEW_OPT( gd_y, "gid dimension : Y ", 20.0 );
	//  NEW_OPT( gd_z, "gid dimension : Z ", 20.0 );
  NEW_OPT( bound_ligand, "use bound ligand to set the grid for generating eggshell", "bound_ligand.pdb" );
	NEW_OPT( add_grid_size, "add grid dimension along x,y,z axis", 2 );
	//	NEW_OPT( add_ext_grid_size, "add extra grid dimension along x,y,z axis", 2 );

	devel::init(argc, argv);

	std::string const input_protein = option[ protein ];
	//	core::Real const grid_cen_x = option[ gc_x ];
	//	core::Real const grid_cen_y = option[ gc_y ];
	//	core::Real const grid_cen_z = option[ gc_z ];
	//	core::Real grid_dim_x = option[ gd_x ];
	//	core::Real grid_dim_y = option[ gd_y ];
	//	core::Real grid_dim_z = option[ gd_z ];
	std::string const bound_ligand_file = option[ bound_ligand ];
	int add_grid_dim  = option[ add_grid_size ];
	//	int add_ext_grid_dim  = option[ add_ext_grid_size ];

	using namespace basic::options;
	core::Real const spacing = option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();

	protocols::pockets::NonPlaidFingerprint npf;
	pose::Pose protein_pose;
	core::import_pose::pose_from_pdb( protein_pose, input_protein );
	pose::Pose bound_ligand_pose;
	core::import_pose::pose_from_pdb( bound_ligand_pose, bound_ligand_file );
	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = bound_ligand_pose.total_residue(); j <= resnum; ++j ) {
		if (!bound_ligand_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}
	if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}
	numeric::xyzVector<core::Real> input_ligand_CoM(0.);
	conformation::Residue const & curr_rsd = bound_ligand_pose.conformation().residue(lig_res_num);
	core::Real minx(999.), miny(999.), minz(999.), maxx(-999.), maxy(-999.), maxz(-999.);
	for(Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i) {
		if (curr_rsd.atom(i).xyz()(1) > maxx){maxx = curr_rsd.atom(i).xyz()(1);}
		if (curr_rsd.atom(i).xyz()(1) < minx){minx = curr_rsd.atom(i).xyz()(1);}
		if (curr_rsd.atom(i).xyz()(2) > maxy){maxy = curr_rsd.atom(i).xyz()(2);}
		if (curr_rsd.atom(i).xyz()(2) < miny){miny = curr_rsd.atom(i).xyz()(2);}
		if (curr_rsd.atom(i).xyz()(3) > maxz){maxz = curr_rsd.atom(i).xyz()(3);}
		if (curr_rsd.atom(i).xyz()(3) < minz){minz = curr_rsd.atom(i).xyz()(3);}
	}
	core::Real num_pts_in_xdim, num_pts_in_ydim, num_pts_in_zdim;
	core::Real const cen_x = (maxx + minx)/2;
	core::Real const cen_y = (maxy + miny)/2;
	core::Real const cen_z = (maxz + minz)/2;

	int xcntr = (cen_x > 0.0) ? floor(cen_x + 0.5) : ceil(cen_x - 0.5);
	int ycntr = (cen_y > 0.0) ? floor(cen_y + 0.5) : ceil(cen_y - 0.5);
	int zcntr = (cen_z > 0.0) ? floor(cen_z + 0.5) : ceil(cen_z - 0.5);

	num_pts_in_xdim = (std::abs(maxx - minx)) + add_grid_dim;
	num_pts_in_ydim = (std::abs(maxy - miny)) + add_grid_dim;
  num_pts_in_zdim = (std::abs(maxz - minz)) + add_grid_dim;

	//round to nearest integer
	int x_pts = ceil(num_pts_in_xdim);
	int y_pts=  ceil(num_pts_in_ydim);
	int z_pts = ceil(num_pts_in_zdim);

	//round up to nearest even unmber
	x_pts =  (x_pts + 1) & ~1;
	y_pts =  (y_pts + 1) & ~1;
	z_pts =  (z_pts + 1) & ~1;

	std::cout<<xcntr<<" "<< ycntr<<" "<< zcntr <<" "<< x_pts <<" "<< y_pts <<" "<< z_pts <<" "<<spacing<<std::endl;

	utility::io::ozstream VINAfile;
  VINAfile.open("adt_vina_config.txt", std::ios::out);
	VINAfile <<"center_x = "<<xcntr<<"\n";
	VINAfile <<"center_y = "<<ycntr<<"\n";
	VINAfile <<"center_z = "<<zcntr<<"\n";
	VINAfile <<"\n";
	VINAfile <<"size_x = "<<x_pts<<"\n";
	VINAfile <<"size_y = "<<y_pts<<"\n";
	VINAfile <<"size_z = "<<z_pts<<"\n";
	VINAfile.close();
	VINAfile.clear();

	std::cout<< "DONE!"<< std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
        std::cerr << "caught exception " << e.msg() << std::endl;
    }
	return 0;

}
