// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/pockets/Fingerprint.cc
/// @brief  protocols::pockets::Fingerprint functions
/// @author Ragul Gowthaman

// Protocol Headers
#include <numeric/constants.hh>
#include <protocols/pockets/DarcElectrostatics.hh>
#include <protocols/pockets/PocketGrid.hh>
// AUTO-REMOVED #include <core/init/init.hh>

// Core Headers
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/PDBInfo.hh>
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/AtomType.hh>
#include <core/types.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>

// Utility Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

#include <utility/vector1.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>

using namespace core;
using namespace core::scoring;
using namespace std;


namespace protocols {
namespace pockets {

DarcElectrostaticsBase::DarcElectrostaticsBase () :
	ReferenceCount()
{}



void DelphiElectrostatics::setup_from_DelphiGrid( std::string const & input_filename, Size const & esp_grid_size, core::Real const & esp_grid_spacing, core::Real const & esp_grid_midpoint_x, core::Real const & esp_grid_midpoint_y, core::Real const & esp_grid_midpoint_z) {

	utility::vector1<core::Real> grid_point(4);

	grid_spacing_ = esp_grid_spacing;
	core::Real midpoint_x = esp_grid_midpoint_x;
	core::Real midpoint_y = esp_grid_midpoint_y;
	core::Real midpoint_z = esp_grid_midpoint_z;
	Size height = esp_grid_size;
	Size width = esp_grid_size;
	Size depth = esp_grid_size;

	std::string line;
  ifstream inFile(input_filename.c_str());
  if (!inFile) {
		std::cout<< "Can't open input file " << input_filename << std::endl;
    exit(1);
  }

	//intialize electrostatic potential grid
	std::vector < std::vector < std::vector <core::Real> > > esp_grid;
	esp_grid.resize(height);
	for( Size i=0; i<height; ++i){
		esp_grid[i].resize(width);
		for (Size j =0; j<width; ++j)
			esp_grid[i][j].resize(depth);
	}

	/*
	for( Size i=1; i<=height; ++i){
		for( Size j=1; j=<width; ++j){
			for( Size k=1; k<=depth; ++k){
				inFile >> esp_grid[i][j][k];
			}
		}
	}
	*/

	for( Size i=1; i<=height; ++i){
		//		grid_point.x() = i * grid_spacing_ + ( midpoint_x - ( ((height+1)/2) * grid_spacing_) );
		grid_point[1] = i * grid_spacing_ + midpoint_x;
		for( Size j=1; j<=width; ++j){
			grid_point[2] = j * grid_spacing_ + midpoint_y;// - ( ((width+1)/2) * grid_spacing_) );
			for( Size k=1; k<=depth; ++k){
				grid_point[3] = k * grid_spacing_ + midpoint_z;// - ( ((depth+1)/2) * grid_spacing_) );
				inFile >> grid_point[4];
				esp_grid_point_list_.push_back(grid_point);
			}
		}
	}



	//DUMP ESP_GRID TO A PDB FILE
	utility::io::ozstream outPDB_stream;
  outPDB_stream.open("esp_grid.pdb", std::ios::out);
	core::Size count =0;
	for (std::list< utility::vector1<core::Real> >::const_iterator pd = esp_grid_point_list_.begin(); pd != esp_grid_point_list_.end(); ++pd) {
		++count;
    outPDB_stream<<"HETATM"<<std::setw(5)<<count<<" C   ESP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->at(1)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->at(2)<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->at(3)<<"      "<<std::setw(6)<<std::fixed<<std::setprecision(2)<<pd->at(4)<<std::endl;
  }
  outPDB_stream.close();
  outPDB_stream.clear();
	//END DUMP PDB FILE
	}

	core::Real DelphiElectrostatics::get_electrostatics_energy(core::pose::Pose const & ligand_pose) {
	//trilinear_interpolation

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = ligand_pose.total_residue(); j <= resnum; ++j ) {
		if (!ligand_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}
  if (lig_res_num == 0){
		std::cout<<"Error, no ligand to Calculate Electrostatic energy" << std::endl;
    exit(1);
  }

	//DUMP ligand_GRID TO A PDB FILE
	utility::io::ozstream outPDB_stream;
  outPDB_stream.open("lig_grid.pdb", std::ios::out);


	core::Real x0, y0, z0,
			x1, y1, z1,
			x_d, y_d, z_d,
			c_00, c_10, c_01, c_11, c_0, c_1, c;
	core::Real V_000 = 0.0;
	core::Real V_100 = 0.0;
	core::Real V_010 = 0.0;
	core::Real V_001 = 0.0;
	core::Real V_110 = 0.0;
	core::Real V_101 = 0.0;
	core::Real V_011 = 0.0;
	core::Real V_111 = 0.0;

	core::Real E_energy(0.);
	numeric::xyzVector<core::Real> ligand(0.);
	conformation::Residue const & curr_rsd = ligand_pose.conformation().residue(lig_res_num);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
    ligand.x() = curr_rsd.atom(i).xyz()(1);
    ligand.y() = curr_rsd.atom(i).xyz()(2);
    ligand.z() = curr_rsd.atom(i).xyz()(3);

 		x0 = (static_cast<int>(ligand.x()/grid_spacing_))*grid_spacing_;
		y0 = (static_cast<int>(ligand.y()/grid_spacing_))*grid_spacing_;
		z0 = (static_cast<int>(ligand.z()/grid_spacing_))*grid_spacing_;
		x1 = x0 + grid_spacing_;
		y1 = y0 + grid_spacing_;  // This was originally x1; I'm assuming it was a bug. ~ Labonte
		z1 = z0 + grid_spacing_;  // This was originally x1; I'm assuming it was a bug. ~ Labonte

		x_d = (ligand.x() - x0)/(x1 - x0);
		y_d = (ligand.y() - y0)/(y1 - y0);
		z_d = (ligand.z() - z0)/(z1 - z0);
		for (std::list< utility::vector1<core::Real> >::const_iterator pd = esp_grid_point_list_.begin(); pd != esp_grid_point_list_.end(); ++pd) {
			if (pd->at(1) == x0 && pd->at(2)== y0 && pd->at(3) == z0 ) V_000 = pd->at(4);
			if (pd->at(1) == x1 && pd->at(2)== y0 && pd->at(3) == z0 ) V_100 = pd->at(4);
			if (pd->at(1) == x0 && pd->at(2)== y1 && pd->at(3) == z0 ) V_010 = pd->at(4);
			if (pd->at(1) == x1 && pd->at(2)== y1 && pd->at(3) == z0 ) V_110 = pd->at(4);
			if (pd->at(1) == x0 && pd->at(2)== y0 && pd->at(3) == z1 ) V_001 = pd->at(4);
			if (pd->at(1) == x1 && pd->at(2)== y1 && pd->at(3) == z0 ) V_101 = pd->at(4);
			if (pd->at(1) == x0 && pd->at(2)== y1 && pd->at(3) == z1 ) V_011 = pd->at(4);
			if (pd->at(1) == x1 && pd->at(2)== y1 && pd->at(3) == z1 ) V_111 = pd->at(4);
		}

		c_00 = (V_000 * (1 - x_d)) + (V_100 * x_d);
		c_10 = (V_010 * (1 - x_d)) + (V_110 * x_d);
		c_01 = (V_001 * (1 - x_d)) + (V_101 * x_d);
		c_11 = (V_011 * (1 - x_d)) + (V_111 * x_d);
		c_0 = (c_00 * (1 - y_d)) + (c_10 * y_d);
		c_1 = (c_01 * (1 - y_d)) + (c_11 * y_d);
		c = (c_0 * (1 - z_d)) + (c_1 * z_d);
		E_energy += c * curr_rsd.atomic_charge(i);

		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0+grid_spacing_<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0+grid_spacing_<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0+grid_spacing_<<std::endl;
		outPDB_stream<<"HETATM"<<std::setw(5)<<1<<" C   GRD A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<x0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<y0+grid_spacing_<<std::setw(8)<<std::fixed<<std::setprecision(3)<<z0+grid_spacing_<<std::endl;
	}

	outPDB_stream.close();
	outPDB_stream.clear();
	//END DUMP PDB FILE

		std::cout<< "Delphi Electrostatics Energy: " << E_energy << std::endl;
	return E_energy;

}

} // Pockets
} // protocols
