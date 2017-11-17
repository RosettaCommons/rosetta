// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/io/pdb/pdb_writer.hh>
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
OPT_KEY( String, bound_ligand )
OPT_KEY( Integer, add_grid_size )
OPT_KEY( String, espGrid_file )

int main( int argc, char * argv [] ) {

	try {

		NEW_OPT( protein, "protein file name", "protein.pdb" );
		NEW_OPT( bound_ligand, "use bound ligand to set the grid for generating eggshell", "" );
		NEW_OPT( add_grid_size, "add grid dimension along x,y,z axis", 2 );
		NEW_OPT( espGrid_file, "input electrostatic potential grid file", "" );


		devel::init(argc, argv);

		std::string const resid(option[ OptionKeys::pocket_grid::central_relax_pdb_num ]);
		std::string const input_protein = option[ protein ];
		std::string const bound_ligand_file = option[ bound_ligand ];
		int add_grid_dim = option[ add_grid_size ];
		std::string const input_espGrid_file = option[ espGrid_file ];

		//core::Real const spacing = option[ OptionKeys::pocket_grid::pocket_grid_spacing ]();

		protocols::pockets::NonPlaidFingerprint npf;
		pose::Pose protein_pose;
		core::import_pose::pose_from_file( protein_pose, input_protein , core::import_pose::PDB_file);

		pose::Pose bound_ligand_pose;
		core::import_pose::pose_from_file( bound_ligand_pose, bound_ligand_file , core::import_pose::PDB_file);
		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = bound_ligand_pose.size(); j <= resnum; ++j ) {
			if ( !bound_ligand_pose.residue(j).is_protein() ) {
				lig_res_num = j;
				break;
			}
		}
		if ( lig_res_num == 0 ) {
			std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
			exit(1);
		}

		numeric::xyzVector<core::Real> input_ligand_CoM(0.);
		conformation::Residue const & curr_rsd = bound_ligand_pose.conformation().residue(lig_res_num);
		core::Real minx(999.), miny(999.), minz(999.), maxx(-999.), maxy(-999.), maxz(-999.);
		core::Real lig_net_charge(0.);
		for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
			if ( curr_rsd.atom(i).xyz()(1) > maxx ) { maxx = curr_rsd.atom(i).xyz()(1);}
			if ( curr_rsd.atom(i).xyz()(1) < minx ) { minx = curr_rsd.atom(i).xyz()(1);}
			if ( curr_rsd.atom(i).xyz()(2) > maxy ) { maxy = curr_rsd.atom(i).xyz()(2);}
			if ( curr_rsd.atom(i).xyz()(2) < miny ) { miny = curr_rsd.atom(i).xyz()(2);}
			if ( curr_rsd.atom(i).xyz()(3) > maxz ) { maxz = curr_rsd.atom(i).xyz()(3);}
			if ( curr_rsd.atom(i).xyz()(3) < minz ) { minz = curr_rsd.atom(i).xyz()(3);}
			lig_net_charge += curr_rsd.atomic_charge(i);

		}
		core::Real x_from_grd_cen, y_from_grd_cen, z_from_grd_cen, x_halfwidth, y_halfwidth, z_halfwidth;
		core::Real cen_x = (maxx + minx)/2;
		core::Real cen_y = (maxy + miny)/2;
		core::Real cen_z = (maxz + minz)/2;
		numeric::xyzVector<core::Real> const grid_center(cen_x,cen_y,cen_z);
		x_halfwidth = std::abs(maxx - minx)/2;
		y_halfwidth = std::abs(maxy - miny)/2;
		z_halfwidth = std::abs(maxz - minz)/2;
		x_from_grd_cen = x_halfwidth + add_grid_dim;
		y_from_grd_cen = y_halfwidth + add_grid_dim;
		z_from_grd_cen = z_halfwidth + add_grid_dim;
		protocols::pockets::PocketGrid pg( cen_x, cen_y, cen_z, x_from_grd_cen, y_from_grd_cen, z_from_grd_cen );
		std::list< numeric::xyzVector<core::Real> > surfacePoints_list;
		std::list< numeric::xyzVector<core::Real> > grid_surfacePoints_list;
		surfacePoints_list = pg.get_connolly_surfacePoints(protein_pose);
		grid_surfacePoints_list = pg.get_connolly_surfacePoints_within_grid(surfacePoints_list);
		core::Real surf_esp =  npf.get_surface_esp(grid_surfacePoints_list, input_espGrid_file);

		//std::list< numeric::xyzVector<core::Real> > lig_surfacePoints_list;
		//lig_surfacePoints_list = pg.get_connolly_surfacePoints(bound_ligand_pose);

		//create 'tag' for eggshell output filename
		int pfounddir = input_protein.find_last_of("/\\");
		int pfounddot = input_protein.find_last_of(".");
		std::string protein_name = input_protein.substr((pfounddir+1),(pfounddot-(pfounddir+1)));
		std::string pro_output_pdbname = protein_name+"_pro_surf.pdb";
		std::string lig_output_pdbname = protein_name+"_lig_surf.pdb";
		std::cout<<"surf_esp : "<<protein_name<<" "<<surf_esp<<std::endl;
		std::cout<<"lig_net_charge : "<<protein_name<<" "<<lig_net_charge<<std::endl;

		utility::io::ozstream outPDB_stream;
		outPDB_stream.open(pro_output_pdbname, std::ios::out);
		for ( std::list< numeric::xyzVector<core::Real> >::const_iterator sf = grid_surfacePoints_list.begin(); sf != grid_surfacePoints_list.end(); ++sf ) {
			outPDB_stream<<"HETATM   "<<std::setw(2)<<3<<"  C   CNY A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->z() <<std::endl;
		}
		outPDB_stream.close();
		outPDB_stream.clear();

		/*outPDB_stream.open(lig_output_pdbname, std::ios::out);
		for (std::list< numeric::xyzVector<core::Real> >::const_iterator sf = lig_surfacePoints_list.begin(); sf != lig_surfacePoints_list.end(); ++sf) {
		outPDB_stream<<"HETATM   "<<std::setw(2)<<3<<"  C   CNY A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->x() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<< sf->z() <<std::endl;
		}
		outPDB_stream.close();
		outPDB_stream.clear();
		*/

	}//try
catch (utility::excn::Exception const & e ) {
	std::cerr << "caught exception " << e.msg() << std::endl;
	return -1;
}
	return 0;

}
