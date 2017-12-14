// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/pockets/Fingerprint.cc
/// @brief  protocols::pockets::Fingerprint functions
/// @author Ragul Gowthaman

// Protocol Headers
#include <numeric/constants.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
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
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
#include <basic/options/keys/pocket_grid.OptionKeys.gen.hh>
#include <basic/database/open.hh>

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

FingerprintBase::FingerprintBase () :
	ReferenceCount()
{
	origin_.zero();
	CoM_.zero();
}

/// @details Auto-generated virtual destructor
FingerprintBase::~FingerprintBase() = default;

void FingerprintBase::print_to_file(std::string const & output_filename) const {

	utility::io::ozstream out_stream;
	out_stream.open(output_filename, std::ios::out);
	out_stream<<"/DIM/"<<" "<<std::fixed<<std::setprecision(4)<< pocketGrid_dim_.x() << "\t" <<std::fixed<<std::setprecision(4)<< pocketGrid_dim_.y() << "\t"<<std::fixed<<std::setprecision(4)<< pocketGrid_dim_.z() <<std::endl;
	out_stream<<"/MID/"<<" "<<std::fixed<<std::setprecision(4)<< pocketGrid_mid_.x() << "\t" <<std::fixed<<std::setprecision(4)<< pocketGrid_mid_.y() << "\t"<<std::fixed<<std::setprecision(4)<< pocketGrid_mid_.z() <<std::endl;
	out_stream<<"/SPACING/"<<" "<<std::fixed<<std::setprecision(4)<< pocketGrid_spacing_<<std::endl;
	out_stream<<"/COM/"<<" "<<std::fixed<<std::setprecision(4)<< CoM_.x() << "\t" <<std::fixed<<std::setprecision(4)<< CoM_.y() << "\t"<<std::fixed<<std::setprecision(4)<< CoM_.z() <<std::endl;
	core::Size curr_origin_index = 1;
	for ( auto i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori, ++curr_origin_index ) {
		out_stream<<"/ORI/"<< " " <<std::fixed<<std::setprecision(4)<< i_mori->x() << "\t" <<std::fixed<<std::setprecision(4)<< i_mori->y() << "\t"<<std::fixed<<std::setprecision(4)<< i_mori->z() << "\t" << curr_origin_index <<std::endl;
	}
	for ( auto const & fi : triplet_fingerprint_data_ ) {
		out_stream<<"/RAY/"<< " " << fi.phi << "\t" <<fi.psi << "\t"<< fi.rho << "\t"<< fi.ori <<std::endl;
	}

	out_stream.close();
	out_stream.clear();

	return;
}

void FingerprintBase::print_to_pdb(std::string const & output_pdbname) const {
	numeric::xyzVector<core::Real> no_translation(0.);
	print_to_pdb( output_pdbname, no_translation );
}

void FingerprintBase::print_to_pdb(std::string const & output_pdbname, numeric::xyzVector<core::Real> const & translation) const {

	utility::io::ozstream out_stream;
	out_stream.open(output_pdbname, std::ios::out);

	core::Size counter = 1;
	std::string record_name = "HETATM";
	core::Size atom_serial_number = counter;
	std::string gap1 = " ";
	std::string atom_name = "C";
	std::string alternate_location_indicator = " ";
	std::string residue_name = "EGG";
	std::string gap2 = " ";
	std::string chain_identifier = "A";
	core::Size residue_sequence_number = counter;
	std::string code_for_insertion_of_residues = " ";
	std::string gap3 = "   ";


	out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()+translation.z()<<std::endl;
	core::Size curr_origin_index = 1;
	for ( auto i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori, ++curr_origin_index ) {
		//set atom number and residue name for PDB file
		atom_serial_number = counter;
		residue_sequence_number = curr_origin_index;
		residue_name = "OR" + utility::to_string(curr_origin_index);
		//print origin PDB coordinates
		out_stream<<
			std::setw(6)<<record_name<<
			std::setw(5)<<atom_serial_number<<
			std::setw(1)<<gap1<<
			std::setw(4)<<atom_name<<
			std::setw(1)<<alternate_location_indicator<<
			std::setw(3)<<residue_name<<
			std::setw(1)<<gap2<<
			std::setw(1)<<chain_identifier<<
			std::setw(4)<<residue_sequence_number<<
			std::setw(1)<<code_for_insertion_of_residues<<
			std::setw(3)<<gap3<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<i_mori->x()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<i_mori->y()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<i_mori->z()<<
			std::endl;
	}
	for ( auto const & pd : triplet_fingerprint_data_ ) {
		numeric::xyzVector<core::Real> new_coor;
		convert_spherical_coor_triplet_to_cartesian( pd, new_coor );
		new_coor += multi_origin_list_[pd.ori];
		//set atom number and residue name for PDB file
		atom_serial_number = counter;
		residue_sequence_number = pd.ori;
		residue_name = "EG" + utility::to_string(pd.ori);
		//print eggshell PDB coordinates
		out_stream<<
			std::setw(6)<<record_name<<
			std::setw(5)<<atom_serial_number<<
			std::setw(1)<<gap1<<
			std::setw(4)<<atom_name<<
			std::setw(1)<<alternate_location_indicator<<
			std::setw(3)<<residue_name<<
			std::setw(1)<<gap2<<
			std::setw(1)<<chain_identifier<<
			std::setw(4)<<residue_sequence_number<<
			std::setw(1)<<code_for_insertion_of_residues<<
			std::setw(3)<<gap3<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.x()+translation.x()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.y()+translation.y()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.z()+translation.z()<<
			std::endl;
	}
	out_stream.close();
	out_stream.clear();

	return;
}

NonPlaidFingerprint::NonPlaidFingerprint() {
#ifdef USEOPENCL
  gpu_.profiling(0);
  gpu_.Init();
  memset(&gpu_memory_, 0, sizeof(gpu_memory_));
#endif
}

NonPlaidFingerprint::~NonPlaidFingerprint() {
	// All GPU memory is freed automatically when gpu_ is destroyed
}

void NonPlaidFingerprint::setup_from_PlaidFingerprint( PlaidFingerprint const & pfp ) {
	origin_ = pfp.origin();
	triplet_fingerprint_data_.resize(pfp.triplet_fingerprint_data().size());
	std::copy (pfp.triplet_fingerprint_data().begin(),pfp.triplet_fingerprint_data().end(), triplet_fingerprint_data_.begin());
	return;
}

void NonPlaidFingerprint::setup_from_espGrid( std::string const & input_filename, core::pose::Pose const &, bool delphi ) {

	ElectrostaticpotentialGrid eGrid;
	if ( delphi ) {
		eGrid.get_DELPHI_espGrid_values(input_filename);
	} else {
		eGrid.get_ZAP_espGrid_values_with_type(input_filename);
	}
	//eGrid.cap_espGrid();

	//Print esp grid points into a PDB file
	using namespace basic::options;
	if ( option[ OptionKeys::pocket_grid::dump_espGrid ]() ) {
		std::string output_espGrid_pdbname = "espGrid.pdb";
		eGrid.write_espGrid_to_pdb(output_espGrid_pdbname);
	}

	esp_spacing_ = eGrid.espGrid_spacing_;
	esp_mid_ = eGrid.espGrid_mid_;
	esp_dim_ = eGrid.espGrid_dim_;
	typGrid_ = eGrid.typGrid_;
	espGrid_ = eGrid.espGrid_;

#ifdef USEOPENCL
	gpu_setup_espGrid();
#endif

	return;

}

//use same grid for both egg and ext shell
void NonPlaidFingerprint::setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid ) {

	PocketGrid grid_for_extshell = pocket_grid;
	setup_from_PocketGrid(protein_pose, pocket_grid, grid_for_extshell);

	return;
}

void NonPlaidFingerprint::setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell ) {
	EggshellGrid egg_sg(pocket_grid);
	eggshell_list_ = egg_sg.eggshell_coord_list();
	EggshellGrid ext_sg(grid_for_extshell, eggshell_list_);
	extshell_list_ = ext_sg.extra_coord_list();

	//combine eggshell & extra coord list into a single list
	egg_and_ext_list_.clear();
	egg_and_ext_list_ = combine_xyz_lists(eggshell_list_ , extshell_list_);

	//set CoM_ to the eggshell CoM
	CoM_ = egg_sg.eggshell_CoM_;
	pocket_CoM_ = egg_sg.eggshell_CoM_;

	//get PocketGrid dimensions
	pocketGrid_mid_ = pocket_grid.mid_;
	pocketGrid_dim_ = pocket_grid.dim_;
	pocketGrid_spacing_ = pocket_grid.spacing_;

	//set origin
	multi_origin_list_.clear();
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_extrashell_to_set_origin ]() ) {
		set_origin( protein_pose, egg_and_ext_list_);
	} else {
		set_origin( protein_pose, eggshell_list_);
	}

	setup_from_EggshellGrid();

	return;
}

void NonPlaidFingerprint::setup_from_ConnollySurface( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell ) {

	EggshellGrid Eshell( protein_pose, pocket_grid, grid_for_extshell );
	eggshell_list_ = Eshell.eggshell_coord_list();
	extshell_list_ = Eshell.extra_coord_list();

	//combine eggshell & extra coord list into a single list
	egg_and_ext_list_.clear();
	egg_and_ext_list_ = combine_xyz_lists(eggshell_list_ , extshell_list_);

	//set CoM_ to the eggshell CoM
	CoM_ = Eshell.eggshell_CoM_;

	//get PocketGrid dimensions
	pocketGrid_mid_ = pocket_grid.mid_;
	pocketGrid_dim_ = pocket_grid.dim_;
	pocketGrid_spacing_ = pocket_grid.spacing_;

	//set origin
	multi_origin_list_.clear();
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_extrashell_to_set_origin ]() ) {
		set_origin( protein_pose, egg_and_ext_list_);
	} else {
		set_origin( protein_pose, eggshell_list_);
	}

	setup_from_EggshellGrid();

	return;
}

void NonPlaidFingerprint::setup_from_PocketGrid_and_known_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist ) {
	EggshellGrid egg_sg(pocket_grid);
	eggshell_list_ = egg_sg.eggshell_coord_list();
	include_eggshell_points_based_on_known_ligand(known_ligand_pose, trim_dist);
	EggshellGrid ext_sg(grid_for_extshell, eggshell_list_);
	extshell_list_ = ext_sg.extra_coord_list();

	//combine eggshell & extra coord list into a single list
	egg_and_ext_list_.clear();
	egg_and_ext_list_ = combine_xyz_lists(eggshell_list_ , extshell_list_);

	//set CoM_ to the eggshell CoM
	CoM_ = egg_sg.eggshell_CoM_;

	//get PocketGrid dimensions
	pocketGrid_mid_ = pocket_grid.mid_;
	pocketGrid_dim_ = pocket_grid.dim_;
	pocketGrid_spacing_ = pocket_grid.spacing_;

	//set origin
	multi_origin_list_.clear();
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_extrashell_to_set_origin ]() ) {
		set_origin( protein_pose, egg_and_ext_list_);
	} else {
		set_origin( protein_pose, eggshell_list_);
	}

	setup_from_EggshellGrid();

	return;
}

void NonPlaidFingerprint::setup_from_PocketGrid_using_bound_ligand( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid, PocketGrid const & grid_for_extshell, core::pose::Pose const & known_ligand_pose ) {

	EggshellGrid egg_sg(pocket_grid,known_ligand_pose);
	eggshell_list_ = egg_sg.eggshell_coord_list();
	EggshellGrid ext_sg(grid_for_extshell, eggshell_list_);
	extshell_list_ = ext_sg.extra_coord_list();

	//combine eggshell & extra coord list into a single list
	egg_and_ext_list_.clear();
	egg_and_ext_list_ = combine_xyz_lists(eggshell_list_ , extshell_list_);

	//set CoM_ to the eggshell CoM
	CoM_ = egg_sg.eggshell_CoM_;

	//get PocketGrid dimensions
	pocketGrid_mid_ = pocket_grid.mid_;
	pocketGrid_dim_ = pocket_grid.dim_;
	pocketGrid_spacing_ = pocket_grid.spacing_;

	//set origin
	multi_origin_list_.clear();
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_extrashell_to_set_origin ]() ) {
		set_origin( protein_pose, egg_and_ext_list_);
	} else {
		set_origin( protein_pose, eggshell_list_);
	}

	setup_from_EggshellGrid();

	return;
}

void NonPlaidFingerprint::write_eggshell_to_pdb_file( std::string const & output_eggshell_name ) const {

	utility::io::ozstream outPDB_stream;
	outPDB_stream.open(output_eggshell_name, std::ios::out);
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()<<std::endl;
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()<<std::endl;

	Size counter = 0;
	for ( auto const & pd : eggshell_list_ ) {
		++counter;
		std::string record_name = "ATOM  ";
		Size atom_serial_number = counter;
		std::string gap1 = " ";
		std::string atom_name = "C";
		std::string alternate_location_indicator = " ";
		std::string residue_name = "EGG";
		std::string gap2 = " ";
		std::string chain_identifier = "A";
		Size residue_sequence_number = counter;
		std::string code_for_insertion_of_residues = " ";
		std::string gap3 = "   ";

		outPDB_stream<<
			std::setw(6)<<record_name<<
			std::setw(5)<<atom_serial_number<<
			std::setw(1)<<gap1<<
			std::setw(4)<<atom_name<<
			std::setw(1)<<alternate_location_indicator<<
			std::setw(3)<<residue_name<<
			std::setw(1)<<gap2<<
			std::setw(1)<<chain_identifier<<
			std::setw(4)<<residue_sequence_number<<
			std::setw(1)<<code_for_insertion_of_residues<<
			std::setw(3)<<gap3<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.x()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.y()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.z()<<
			std::endl;
	}

	counter = 0;
	for ( auto const & pd : extshell_list_ ) {
		++counter;
		std::string record_name = "ATOM  ";
		Size atom_serial_number = counter;
		std::string gap1 = " ";
		std::string atom_name = "O";
		std::string alternate_location_indicator = " ";
		std::string residue_name = "EXT";
		std::string gap2 = " ";
		std::string chain_identifier = "B";
		Size residue_sequence_number = counter;
		std::string code_for_insertion_of_residues = " ";
		std::string gap3 = "   ";

		outPDB_stream<<
			std::setw(6)<<record_name<<
			std::setw(5)<<atom_serial_number<<
			std::setw(1)<<gap1<<
			std::setw(4)<<atom_name<<
			std::setw(1)<<alternate_location_indicator<<
			std::setw(3)<<residue_name<<
			std::setw(1)<<gap2<<
			std::setw(1)<<chain_identifier<<
			std::setw(4)<<residue_sequence_number<<
			std::setw(1)<<code_for_insertion_of_residues<<
			std::setw(3)<<gap3<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.x()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.y()<<
			std::setw(8)<<std::fixed<<std::setprecision(3)<<pd.z()<<
			std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

}

void NonPlaidFingerprint::set_origin( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell ) {

	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::multiple_origin ]() ) {
		set_multiple_origin( protein_pose, egg_and_extra_shell );
	} else {
		int origin_option = option[ OptionKeys::fingerprint::set_origin ];
		set_origin_from_option_( protein_pose, egg_and_extra_shell, origin_option );
		multi_origin_list_.push_back(origin_);
	}

	return;

}

void NonPlaidFingerprint::choose_origin_with_lowest_eggshell_ruggedness( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell ) {

	//choose best origin by finding the lowest R (ruggedness) value from different origin position
	Size best_origin_option(0);
	core::Real best_R(999.);
	for ( Size set_origin_option = 1; set_origin_option < 4; ++set_origin_option ) {
		core::Real new_R = get_Rvalue( protein_pose, egg_and_extra_shell, set_origin_option );
		//std::cout << "set_origin_option " << set_origin_option << " new_Rvalue " << new_R << std::endl;
		if ( new_R < best_R ) {
			best_R = new_R;
			best_origin_option = set_origin_option;
		}
	}
	//std::cout << "best_origin_option " << best_origin_option << " best_Rvalue " << best_R << std::endl;
	set_origin_from_option_( protein_pose, egg_and_extra_shell, best_origin_option );

	return;

}


core::Real NonPlaidFingerprint::get_Rvalue( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option ) {

	set_origin_from_option_( protein_pose, egg_and_extra_shell, set_origin_option );

	core::Real min_rho_difference(0.);
	core::Size num_points(0);
	for ( auto pt1 = egg_and_extra_shell.begin(); pt1 != egg_and_extra_shell.end(); ++pt1 ) {
		numeric::xyzVector< core::Real > eggshell_point1 = *pt1 - origin_;
		spherical_coor_triplet triplet1 = {0,0,0,0};
		convert_cartesian_to_spherical_coor_triplet( eggshell_point1, triplet1 );
		core::Real min_angle(999.);
		spherical_coor_triplet triplet2 = {0,0,0,0};
		for ( auto pt2 = egg_and_extra_shell.begin(); pt2 != egg_and_extra_shell.end(); ++pt2 ) {
			if ( pt1 == pt2 ) continue;
			core::Real const curr_angle = std::abs(cos_of( *pt1, *pt2 ));
			if ( curr_angle < min_angle ) {
				min_angle = curr_angle;
				numeric::xyzVector< core::Real > eggshell_point2 = *pt2 - origin_;
				convert_cartesian_to_spherical_coor_triplet( eggshell_point2, triplet2 );
			}
		}
		min_rho_difference += std::abs( triplet1.rho - triplet2.rho );
		num_points++;
	}
	return min_rho_difference / num_points;
}

void NonPlaidFingerprint::set_origin_from_option_( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option ) {

	if ( set_origin_option == 0 ) {
		choose_origin_with_lowest_eggshell_ruggedness(protein_pose, egg_and_extra_shell);
		//std::cout<< " ORIGIN1 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else if ( set_origin_option == 1 ) {
		set_origin_away_from_protein_center(protein_pose);
		//std::cout<< " ORIGIN1 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else if ( set_origin_option == 2 ) {
		set_origin_away_from_eggshell(egg_and_extra_shell, protein_pose);
		//std::cout<< " ORIGIN2 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else if ( set_origin_option == 3 ) {
		set_origin_away_from_eggshell_plane(egg_and_extra_shell, protein_pose, set_origin_option);
		//std::cout<< " ORIGIN3 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else if ( set_origin_option == 4 ) {
		set_origin_away_from_eggshell_plane(egg_and_extra_shell, protein_pose, set_origin_option);
		//std::cout<< " ORIGIN4 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else if ( set_origin_option == 5 ) {
		set_origin_from_residue(protein_pose);
		//std::cout<< " ORIGIN5 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
	} else {
		std::cout<<"Error, wrong option to set_origin" << std::endl;
		exit(1);
	}
	return;
}

void NonPlaidFingerprint::set_multiple_origin( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell ) {

	//choose initial origin
	using namespace basic::options;
	int origin_option = option[ OptionKeys::fingerprint::set_origin ];
	set_origin_from_option_( protein_pose, egg_and_extra_shell, origin_option );
	multi_origin_list_.push_back(origin_);

	core::Real angle;
	numeric::xyzVector<core::Real> new_origin;
	numeric::xyzVector<core::Real> tmp_origin;
	numeric::xyzVector<core::Real> random_pt;

	random_pt.x() = 5.;
	random_pt.y() = 0.;
	random_pt.z() = 0.;

	numeric::xyzVector<core::Real> vec1 = origin_ - CoM_;
	numeric::xyzVector<core::Real> vec2 = random_pt - CoM_;
	numeric::xyzVector<core::Real> proj_vector = cross_product( vec1, vec2 );
	numeric::xyzMatrix<core::Real> rot_matrix;

	angle = 45.;
	rot_matrix = numeric::rotation_matrix_degrees( proj_vector, angle );
	tmp_origin = origin_ - CoM_;
	new_origin = numeric::product(rot_matrix, tmp_origin);
	new_origin += CoM_;
	multi_origin_list_.push_back(new_origin);

	angle = -45.;
	rot_matrix = numeric::rotation_matrix_degrees( proj_vector, angle );
	tmp_origin = origin_ - CoM_;
	new_origin = numeric::product(rot_matrix, tmp_origin);
	new_origin += CoM_;
	multi_origin_list_.push_back(new_origin);

	numeric::xyzVector<core::Real> new_proj_vector = cross_product( vec1, proj_vector );

	angle = 45.;
	rot_matrix = numeric::rotation_matrix_degrees( new_proj_vector, angle );
	tmp_origin = origin_ - CoM_;
	new_origin = numeric::product(rot_matrix, tmp_origin);
	new_origin += CoM_;
	multi_origin_list_.push_back(new_origin);

	angle = -45.;
	rot_matrix = numeric::rotation_matrix_degrees( new_proj_vector, angle );
	tmp_origin = origin_ - CoM_;
	new_origin = numeric::product(rot_matrix, tmp_origin);
	new_origin += CoM_;
	multi_origin_list_.push_back(new_origin);

	return;
}

numeric::xyzVector<core::Real> NonPlaidFingerprint::place_origin_point( core::Real const & angle ) {

	//move old-origin to reference point 'origin'
	numeric::xyzVector<core::Real> tmp_origin = origin_ - CoM_;
	//multiply by rotation matrix
	numeric::xyzVector<core::Real> origin_pt;
	numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(angle) );
	numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(angle) );
	numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(angle) );
	numeric::xyzMatrix<core::Real> tot_rot_mat = z_rot_mat * y_rot_mat * x_rot_mat;
	origin_pt = numeric::product(tot_rot_mat, tmp_origin);
	//move it back from reference point 'origin'
	origin_pt += CoM_;

	std::cout<< " ORIGIN "<<angle<<" "<<origin_pt.x()<<" "<<origin_pt.y()<<" "<<origin_pt.z()<<std::endl;

	return origin_pt;

}

void NonPlaidFingerprint::setup_from_EggshellGrid() {

	//find best origin for Eggshell
	//best origin is the shortest point to the eggshell
	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_list_.begin(); pd != eggshell_list_.end(); ++pd ) {
		core::Real best_rho = 9999.;
		core::Size curr_origin_index = 0;
		core::Size best_origin_index = 0;
		spherical_coor_triplet best_triplet = {0,0,0,0};
		for ( utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori ) {
			++curr_origin_index;
			spherical_coor_triplet ray_triplet;
			convert_cartesian_to_spherical_coor_triplet( *pd - *i_mori, ray_triplet );
			if ( ray_triplet.rho < best_rho ) {
				best_rho = ray_triplet.rho;
				best_triplet = ray_triplet;
				best_origin_index = curr_origin_index;
			}
		}
		spherical_coor_triplet tmp_ray = {best_triplet.phi, best_triplet.psi, best_triplet.rho, best_origin_index};
		triplet_fingerprint_data_.push_back(tmp_ray);
	}
	//find best origin for Extrashell
	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator pd = extshell_list_.begin(); pd != extshell_list_.end(); ++pd ) {
		core::Real best_rho = 9999.;
		core::Size curr_origin_index = 0;
		core::Size best_origin_index = 0;
		spherical_coor_triplet best_triplet = {0,0,0,0};
		for ( utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori ) {
			++curr_origin_index;
			spherical_coor_triplet ray_triplet;
			convert_cartesian_to_spherical_coor_triplet( *pd - *i_mori, ray_triplet );
			if ( ray_triplet.rho < best_rho ) {
				best_rho = ray_triplet.rho;
				best_triplet = ray_triplet;
				best_origin_index = curr_origin_index;
			}
		}
		spherical_coor_triplet tmp_ray = {best_triplet.phi, best_triplet.psi, 0.0, best_origin_index};
		triplet_fingerprint_data_.push_back(tmp_ray);
	}
	/*
	using namespace basic::options;
	if (option[ OptionKeys::fingerprint::multiple_origin ]()){

	NonPlaidFingerprint pair_list;
	for (std::list<numeric::xyzVector<core::Real> >::const_iterator i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori) {
	std::cout<< "morigin : " << i_mori->x() <<" "<< i_mori->y() <<" "<< i_mori->z() <<std::endl;
	std::list< spherical_coor_triplet > ray_triplet_list;
	ray_triplet_list.clear();
	//for eggshell
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_list_.begin(); pd != eggshell_list_.end(); ++pd) {
	spherical_coor_triplet best_triplet = {0,0,0,0};
	convert_cartesian_to_spherical_coor_triplet( *pd - *i_mori, best_triplet );
	core::Real best_rho = best_triplet.rho;
	bool best_ori = false;
	for (std::list<numeric::xyzVector<core::Real> >::const_iterator j_mori = multi_origin_list_.begin(); j_mori != multi_origin_list_.end(); ++j_mori) {
	spherical_coor_triplet ray_triplet;
	convert_cartesian_to_spherical_coor_triplet( *pd - *j_mori, ray_triplet );
	if(ray_triplet.rho < best_rho)
	{
	best_ori = false;
	break;
	}
	else
	{
	best_ori = true;
	}
	}
	if(best_ori)
	{
	ray_triplet_list.push_back(best_triplet);
	}
	}
	//for extshell
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = extshell_list_.begin(); pd != extshell_list_.end(); ++pd) {
	spherical_coor_triplet best_ext_triplet;
	convert_cartesian_to_spherical_coor_triplet( *pd - *i_mori, best_ext_triplet );
	core::Real best_ext_rho = best_ext_triplet.rho;
	bool best_ext_ori = false;
	for (std::list<numeric::xyzVector<core::Real> >::const_iterator j_mori = multi_origin_list_.begin(); j_mori != multi_origin_list_.end(); ++j_mori) {
	spherical_coor_triplet ray_triplet;
	convert_cartesian_to_spherical_coor_triplet( *pd - *j_mori, ray_triplet );
	if(ray_triplet.rho < best_ext_rho)
	{
	best_ext_ori = false;
	break;
	}
	else
	{
	best_ext_ori = true;
	}
	}
	if(best_ext_ori)
	{
	spherical_coor_triplet tmp_triplet;
	tmp_triplet = best_ext_triplet;
	tmp_triplet.rho = 0;
	ray_triplet_list.push_back(tmp_triplet);
	}
	}
	pair_list.origin_ = *i_mori;
	pair_list.triplet_fingerprint_data_ = ray_triplet_list;
	npf_list_.push_back(pair_list);
	}
	}
	else {
	//write_eggshell_to_pdb_file("original_eggshell.pdb");
	std::cout<< "sorigin : " << origin_.x() <<" "<< origin_.y() <<" "<< origin_.z() <<std::endl;
	// convert from cartesian eggshell to spherical coors
	std::list< spherical_coor_triplet > rounded_egg_triplet = convert_cart_to_spherical_and_round(eggshell_list_);
	std::list< spherical_coor_triplet > rounded_ext_triplet = convert_cart_to_spherical_and_round(extshell_list_);
	std::list< spherical_coor_triplet > unq_egg_triplet = remove_duplicate_phi_psi(rounded_egg_triplet);
	std::list< spherical_coor_triplet > unq_ext_triplet = remove_duplicate_phi_psi(rounded_ext_triplet);
	//remove ext_tripltet that matches with egg_triplet phi psi angles
	for (std::list<spherical_coor_triplet>::const_iterator aa = unq_egg_triplet.begin(); aa != unq_egg_triplet.end(); ++aa) {
	for (std::list<spherical_coor_triplet>::iterator bb = unq_ext_triplet.begin(); bb != unq_ext_triplet.end();) {
	if( (aa->phi == bb->phi) && (aa->psi == bb->psi) ) {
	bb = unq_ext_triplet.erase(bb);
	}
	else {
	++bb;
	}
	}
	}
	write_eggshell_to_pdb_file("original_eggshell.pdb");
	eggshell_list_.clear();
	eggshell_list_ = convert_spherical_list_to_cartesian_list(unq_egg_triplet);
	extshell_list_.clear();
	extshell_list_ = convert_spherical_list_to_cartesian_list(unq_ext_triplet);
	triplet_fingerprint_data_.insert(triplet_fingerprint_data_.end(), unq_egg_triplet.begin(), unq_egg_triplet.end());
	unq_ext_triplet = set_rho_to_zero(unq_ext_triplet);
	triplet_fingerprint_data_.insert(triplet_fingerprint_data_.end(), unq_ext_triplet.begin(), unq_ext_triplet.end());

	//MAKE COMPATIBLE WITH MULTI ORIGIN
	NonPlaidFingerprint pair_list;
	pair_list.origin_ = origin_;
	pair_list.triplet_fingerprint_data_ = triplet_fingerprint_data_;
	npf_list_.push_back(pair_list);
	//write_eggshell_to_pdb_file("modified_eggshell.pdb");
	}
	*/
	return;
}

#ifdef USEOPENCL
	void NonPlaidFingerprint::gpu_setup( core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, int & num_particles, core::Real const & electrostatics_weight ) {

		// Initialize GPU memory and variables
		if(!gpu_.use()) return;
		if(!gpu_.RegisterProgram("gpu/DARC_VERSION2.cl")) utility_exit_with_message("Failed to load OpenCL program");

		//Get no. of particles and write to GPU memory
		gpu_memory_.num_particles = num_particles;

		// Custom weights
		float weights[4] = {static_cast<float>(missing_point_weight), static_cast<float>(steric_weight), static_cast<float>(extra_point_weight), static_cast<float>(num_particles)};
		if(!gpu_memory_.weights) gpu_memory_.weights = gpu_.AllocateMemory(sizeof(weights));
		gpu_.WriteData(gpu_memory_.weights, weights, sizeof(weights));

		// electrostatic weights
		float elec_weights[4] = {static_cast<float>(electrostatics_weight),static_cast<float>(num_particles),0.0,0.0};
		if(!gpu_memory_.elec_weights) gpu_memory_.elec_weights = gpu_.AllocateMemory(sizeof(elec_weights));
		gpu_.WriteData(gpu_memory_.elec_weights, elec_weights, sizeof(elec_weights));

	}

	void NonPlaidFingerprint::gpu_setup_atomcoords(PlaidFingerprint & pf, int & num_particles ) {
		if(!gpu_.use()) return;

		core::Size const nconformers = pf.compute_ligand_nconformers();

		using namespace basic::options;
		core::Real const radius_scale = option[ OptionKeys::fingerprint::atom_radius_scale ];
		core::Real const atom_buffer = option[ OptionKeys::fingerprint::atom_radius_buffer ];

		//Get no. of atoms present in the input ligand and copy to GPU memory
		//for shape-only calculation we are not including hydrogens
		//for electrosstatic calculations we include hydrogens
		//so we keep two copies of atomcoords with and without hydrogens

		//for shape only calculations with or without hydrogens(option -include_hydrogens)
		core::Size const ligand_natoms_shapecalc = pf.compute_ligand_natoms();
		gpu_memory_.num_ligatoms_shapecalc = ligand_natoms_shapecalc;
		gpu_memory_.num_totatoms_shapecalc = ligand_natoms_shapecalc * num_particles;

		//for electrostatics calculations with hydrogens
		core::Size const ligand_natoms_elstscalc = pf.compute_ligand_natoms_with_hydrogens();
		gpu_memory_.num_ligatoms_elstscalc = ligand_natoms_elstscalc;
		gpu_memory_.num_totatoms_elstscalc = ligand_natoms_elstscalc * num_particles;

		std::vector<basic::gpu::float4> atomcoords_shapecalc;
		std::vector<basic::gpu::float4> atomcoords_elstscalc;
		std::vector<basic::gpu::float4> ligCoM;
		atomcoords_shapecalc.clear();
		atomcoords_elstscalc.clear();
		ligCoM.clear();

		for ( core::Size c = 0; c < nconformers; ++c ) {
			core::pose::PoseOP tmp_pose(new core::pose::Pose(pf.pose(), c+1, c+1));
			core::Size const lig_res_num = pf.compute_ligand_resnum(*tmp_pose);
			core::conformation::ResidueCOP ligand_rsd( core::conformation::ResidueOP( new core::conformation::Residue( tmp_pose->conformation().residue(lig_res_num) ) ) );
			numeric::xyzVector<core::Real> ligand_com(0.);
			//copy atomcoords without hydrogens (included 'radius') for shapeonly calculations
			//aslo copy center of mass of each conformers
			for (Size i = 1; i <= ligand_natoms_shapecalc; ++i) {
				core::Real const this_atom_radius = ( ligand_rsd->atom_type(i).lj_radius() - atom_buffer ) * radius_scale;
				basic::gpu::float4 atomcoord_shapecalc = { static_cast<float>(ligand_rsd->atom(i).xyz()(1)), static_cast<float>(ligand_rsd->atom(i).xyz()(2)), static_cast<float>(ligand_rsd->atom(i).xyz()(3)), static_cast<float>(this_atom_radius) };
				atomcoords_shapecalc.push_back(atomcoord_shapecalc);
				ligand_com.x() += ligand_rsd->atom(i).xyz()(1);
				ligand_com.y() += ligand_rsd->atom(i).xyz()(2);
				ligand_com.z() += ligand_rsd->atom(i).xyz()(3);
			}
			ligand_com /= ligand_natoms_shapecalc;
			basic::gpu::float4 lig_conf_com = { static_cast<float>(ligand_com.x()), static_cast<float>(ligand_com.y()), static_cast<float>(ligand_com.z()), static_cast<float>(ligand_natoms_shapecalc) };
			ligCoM.push_back(lig_conf_com);

			//copy atomcoords with hydrogens (and 'charges') for electrostatic calculations
			for (Size i = 1; i <= ligand_natoms_elstscalc; ++i) {
				basic::gpu::float4 atomcoord_elstscalc = { static_cast<float>(ligand_rsd->atom(i).xyz()(1)), static_cast<float>(ligand_rsd->atom(i).xyz()(2)), static_cast<float>(ligand_rsd->atom(i).xyz()(3)), static_cast<float>(ligand_rsd->atomic_charge(i)) };
			atomcoords_elstscalc.push_back(atomcoord_elstscalc);
			}
		}

		gpu_.AllocateMemoryReuse(gpu_memory_.ligCoM, gpu_memory_.ligCoM_size, sizeof(basic::gpu::float4) * ligCoM.size());
		gpu_.WriteData(gpu_memory_.ligCoM, &ligCoM[0], gpu_memory_.ligCoM_size);

		gpu_.AllocateMemoryReuse(gpu_memory_.atomcoords_shapecalc, gpu_memory_.atomcoords_shapecalc_size, sizeof(basic::gpu::float4) * atomcoords_shapecalc.size());
		gpu_.WriteData(gpu_memory_.atomcoords_shapecalc, &atomcoords_shapecalc[0], gpu_memory_.atomcoords_shapecalc_size);

		gpu_.AllocateMemoryReuse(gpu_memory_.atomcoords_elstscalc, gpu_memory_.atomcoords_elstscalc_size, sizeof(basic::gpu::float4) * atomcoords_elstscalc.size());
		gpu_.WriteData(gpu_memory_.atomcoords_elstscalc, &atomcoords_elstscalc[0], gpu_memory_.atomcoords_elstscalc_size);

		return;

	}

	void NonPlaidFingerprint::gpu_setup_rays() {

		if(!gpu_.use()) return;

		std::vector<basic::gpu::float4> rays;
		std::vector<basic::gpu::float4> RAYorigins;

  for (utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori) {
		basic::gpu::float4 ray_origin = { static_cast<float>(i_mori->x()), static_cast<float>(i_mori->y()), static_cast<float>(i_mori->z()), 0.0 };
			RAYorigins.push_back(ray_origin);
		}

		for (std::list<spherical_coor_triplet>::const_iterator gi = triplet_fingerprint_data_.begin(); gi != triplet_fingerprint_data_.end(); ++gi) {
			basic::gpu::float4 ray = { static_cast<float>(gi->phi), static_cast<float>(gi->psi), static_cast<float>(gi->rho), static_cast<float>(gi->ori) };
			rays.push_back(ray);
		}

    // Allocate GPU memory (reused) and copy ray coordinates over to GPU
    unsigned int rays_size = sizeof(basic::gpu::float4) * rays.size();
    gpu_.AllocateMemoryReuse(gpu_memory_.rays, gpu_memory_.rays_size, rays_size);
    gpu_.WriteData(gpu_memory_.rays, &rays[0], rays_size);
    gpu_memory_.num_rays = rays.size();

    // Allocate GPU memory (reused) and copy ray coordinates over to GPU
    unsigned int RAYorigins_size = sizeof(basic::gpu::float4) * RAYorigins.size();
    gpu_.AllocateMemoryReuse(gpu_memory_.RAYorigins, gpu_memory_.RAYorigins_size, RAYorigins_size);
    gpu_.WriteData(gpu_memory_.RAYorigins, &RAYorigins[0], RAYorigins_size);

	}

	void NonPlaidFingerprint::gpu_setup_espGrid() {

		if(!gpu_.use()) return;

		float grid_dim[4] = {
			static_cast<float>(esp_dim_.x()),
			static_cast<float>(esp_dim_.y()),
			static_cast<float>(esp_dim_.z()),
			static_cast<float>(esp_spacing_)
		};

		float grid_mid[4] = {
			static_cast<float>(esp_mid_.x()),
			static_cast<float>(esp_mid_.y()),
			static_cast<float>(esp_mid_.z()),
			static_cast<float>(esp_spacing_)
		};

		if(!gpu_memory_.grid_dim) gpu_memory_.grid_dim = gpu_.AllocateMemory(sizeof(grid_dim));
		gpu_.WriteData(gpu_memory_.grid_dim, grid_dim, sizeof(grid_dim));

		if(!gpu_memory_.grid_mid) gpu_memory_.grid_mid = gpu_.AllocateMemory(sizeof(grid_mid));
		gpu_.WriteData(gpu_memory_.grid_mid, grid_mid, sizeof(grid_mid));

		//		int num_esp_grid_points = espGrid_.size();
		std::vector<float> gpu_espGrid;
		std::vector<float> gpu_typGrid;

		for (Size ix=0; ix<esp_dim_.x(); ++ix){
			for (Size iy=0; iy<esp_dim_.y(); ++iy){
				for (Size iz=0; iz<esp_dim_.z(); ++iz){
					float esp = espGrid_[ix][iy][iz];
					float typ = typGrid_[ix][iy][iz];
					gpu_espGrid.push_back(esp);
					gpu_typGrid.push_back(typ);
				}
			}
		}

    // Allocate GPU memory (reused) and copy ray coordinates over to GPU
		unsigned int gpu_espGrid_size = sizeof(float) * gpu_espGrid.size();
		gpu_.AllocateMemoryReuse(gpu_memory_.gpu_espGrid, gpu_memory_.gpu_espGrid_size, gpu_espGrid_size);
		gpu_.WriteData(gpu_memory_.gpu_espGrid, &gpu_espGrid[0], gpu_espGrid_size);

		unsigned int gpu_typGrid_size = sizeof(float) * gpu_typGrid.size();
		gpu_.AllocateMemoryReuse(gpu_memory_.gpu_typGrid, gpu_memory_.gpu_typGrid_size, gpu_typGrid_size);
		gpu_.WriteData(gpu_memory_.gpu_typGrid, &gpu_typGrid[0], gpu_typGrid_size);

	}

	int NonPlaidFingerprint::gpu_calculate_particle_scores( core::optimization::ParticleOPs & particles, std::vector<basic::gpu::float4> &particles_rotation_offset, std::vector<basic::gpu::float4> &particles_translation_offset){

		if(!gpu_.use()) return 0;

		// For timing of the rest of the function from here on:
		// basic::gpu::Timer t("gpu_calculate_particle_scores");
		gpu_.AllocateMemoryReuse(gpu_memory_.particles_translation_offset, gpu_memory_.particles_translation_offset_size, sizeof(basic::gpu::float4) * particles_translation_offset.size());
		gpu_.WriteData(gpu_memory_.particles_translation_offset, &particles_translation_offset[0], gpu_memory_.particles_translation_offset_size);
		gpu_.AllocateMemoryReuse(gpu_memory_.particles_rotation_offset, gpu_memory_.particles_rotation_offset_size, sizeof(basic::gpu::float4) * particles_rotation_offset.size());
		gpu_.WriteData(gpu_memory_.particles_rotation_offset, &particles_rotation_offset[0], gpu_memory_.particles_rotation_offset_size);

		gpu_.AllocateMemoryReuse(gpu_memory_.particletest_scores, gpu_memory_.particletest_scores_size, sizeof(float) * gpu_memory_.num_particles);
		gpu_.AllocateMemoryReuse(gpu_memory_.raytest_scores, gpu_memory_.raytest_scores_size, sizeof(float) * gpu_memory_.num_particles * gpu_memory_.num_rays);

		//		std::cout<<"particletest_scores_size : "<<espGrid_.size()<<" "<<gpu_espGrid_size<<" "<<gpu_espGrid.size()<<" "<<gpu_memory_.gpu_espGrid_size<<std::endl;

		// Execute kernels
		if(!gpu_.ExecuteKernel("CALC_RAY_SCORE", gpu_memory_.num_rays, gpu_memory_.num_rays, 64,
													 GPU_DEVMEM, gpu_memory_.rays,
													 GPU_DEVMEM, gpu_memory_.atomcoords_shapecalc,
													 GPU_DEVMEM, gpu_memory_.particles_translation_offset,
													 GPU_DEVMEM, gpu_memory_.particles_rotation_offset,
													 GPU_DEVMEM, gpu_memory_.ligCoM,
													 GPU_DEVMEM, gpu_memory_.weights,
													 GPU_DEVMEM, gpu_memory_.RAYorigins,
													 GPU_DEVMEM, gpu_memory_.raytest_scores,
													 GPU_IN | GPU_INT, gpu_memory_.num_rays,
													 NULL)) {
			std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
			return 0;
		}

		using namespace basic::options;
		if (!option[ OptionKeys::fingerprint::darc_shape_only ].user()) {
			gpu_.AllocateMemoryReuse(gpu_memory_.elec_scores, gpu_memory_.elec_scores_size, sizeof(float) * gpu_memory_.num_particles * gpu_memory_.num_ligatoms_elstscalc);
			if(!gpu_.ExecuteKernel("CALC_ELECTROSTATICS_SCORE", gpu_memory_.num_totatoms_elstscalc, gpu_memory_.num_totatoms_elstscalc, 64,
														 GPU_DEVMEM, gpu_memory_.atomcoords_elstscalc,
														 GPU_DEVMEM, gpu_memory_.particles_translation_offset,
														 GPU_DEVMEM, gpu_memory_.particles_rotation_offset,
														 GPU_DEVMEM, gpu_memory_.ligCoM,
														 GPU_DEVMEM, gpu_memory_.elec_scores,
														 GPU_DEVMEM, gpu_memory_.gpu_espGrid,
														 GPU_DEVMEM, gpu_memory_.gpu_typGrid,
														 GPU_DEVMEM, gpu_memory_.grid_dim,
														 GPU_DEVMEM, gpu_memory_.grid_mid,
														 GPU_IN | GPU_INT, gpu_memory_.num_ligatoms_elstscalc,
														 GPU_IN | GPU_INT, gpu_memory_.num_totatoms_elstscalc,
														 NULL)) {
				std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
				return 0;
			}

			if(!gpu_.ExecuteKernel("GET_DARC_SHAPE_AND_ELECTROSTATICS_SCORE", gpu_memory_.num_particles, gpu_memory_.num_particles, 64,
														 GPU_DEVMEM, gpu_memory_.raytest_scores,
														 GPU_DEVMEM, gpu_memory_.elec_scores,
														 GPU_DEVMEM, gpu_memory_.elec_weights,
														 GPU_DEVMEM, gpu_memory_.particletest_scores,
														 GPU_IN | GPU_INT, gpu_memory_.num_ligatoms_elstscalc,
														 GPU_IN | GPU_INT, gpu_memory_.num_rays,
														 GPU_IN | GPU_INT, gpu_memory_.num_particles,
														 NULL)) {
				std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
				return 0;
			}
		}
		else {
			if(!gpu_.ExecuteKernel("GET_DARC_SHAPE_SCORE", gpu_memory_.num_particles, gpu_memory_.num_particles, 64,
														 GPU_DEVMEM, gpu_memory_.raytest_scores,
														 GPU_DEVMEM, gpu_memory_.particletest_scores,
														 GPU_IN | GPU_INT, gpu_memory_.num_rays,
														 GPU_IN | GPU_INT, gpu_memory_.num_particles,
														 NULL)) {
				std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
				return 0;
			}
		}

		std::vector<float> particletest_scores(gpu_memory_.num_particles);
		gpu_.ReadData(&particletest_scores[0], gpu_memory_.particletest_scores, gpu_memory_.particletest_scores_size);

		Size j =0;
		for (std::vector<float>::const_iterator ci = particletest_scores.begin(); ci != particletest_scores.end(); ++ci) {
			core::Real score = (float)*ci;
			//std::cout<<"Particle score(gpu): "<<score<<std::endl;
			particles[++j]->set_score(score);
		}

		return 1;
	}


#endif

void NonPlaidFingerprint::setup_from_eggshell_pdb_file(std::string const & input_filename) {

	ifstream inFile(input_filename.c_str());
	if ( !inFile ) {
		std::cout<< "Can't open input file " << input_filename << std::endl;
		exit(1);
	}

	std::list< numeric::xyzVector<core::Real> > temp_eggshell_coord_list;
	std::list< numeric::xyzVector<core::Real> > temp_extra_coord_list;
	std::list< numeric::xyzVector<core::Real> > egg_and_extra_coord_list;
	temp_eggshell_coord_list.clear();
	temp_extra_coord_list.clear();
	egg_and_extra_coord_list.clear();

	std::string line;
	std::string name;
	//std::string chainID;
	std::string restype;

	numeric::xyzVector<core::Real> eggshell_CoM;
	core::Size comcounter = 0;
	core::Size oricounter = 0;

	while ( std::getline(inFile, line) ) {
		name = line.substr(0,6);
		//chainID = line.substr(21,1);
		restype = line.substr(17,3);

		numeric::xyzVector<core::Real> pdb_coord;
		std::string Xstring, Ystring, Zstring;

		if ( name=="END" ) {
			break;
		} else if ( (name=="HETATM")&&(restype=="ORI") ) {
			Xstring = line.substr(30,8);
			origin_.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			origin_.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			origin_.z() = atof(Zstring.c_str());
			oricounter++;
		} else if ( (name=="HETATM")&&(restype=="COM") ) {
			Xstring = line.substr(30,8);
			CoM_.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			CoM_.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			CoM_.z() = atof(Zstring.c_str());
			comcounter++;
		} else if ( (name=="HETATM")&&(restype=="EGG") ) {
			Xstring = line.substr(30,8);
			pdb_coord.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			pdb_coord.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			pdb_coord.z() = atof(Zstring.c_str());
			temp_eggshell_coord_list.push_back(pdb_coord);
		} else if ( (name=="HETATM")&&(restype=="EXT") ) {
			Xstring = line.substr(30,8);
			pdb_coord.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			pdb_coord.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			pdb_coord.z() = atof(Zstring.c_str());
			temp_extra_coord_list.push_back(pdb_coord);
		}
	}
	inFile.close();

	if ( (oricounter == 0) || (comcounter == 0) ) {
		std::cout<<"Error, No ORIGIN or CoM value specified in input Eggshell file" << std::endl;
		exit(1);
	} else if ( (oricounter > 1) || (comcounter >1) ) {
		std::cout<<"Error, More than one ORIGIN or CoM value specified in input Eggshell file" << std::endl;
		exit(1);
	}

	// convert from cartesian_coord to spherical coors
	spherical_coor_triplet new_triplet;
	triplet_fingerprint_data_.clear();

	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator pd = temp_eggshell_coord_list.begin(); pd != temp_eggshell_coord_list.end(); ++pd ) {
		convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
		triplet_fingerprint_data_.push_back(new_triplet);
	}

	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator pd = temp_extra_coord_list.begin(); pd != temp_extra_coord_list.end(); ++pd ) {
		convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
		new_triplet.rho = 0.;
		triplet_fingerprint_data_.push_back(new_triplet);
	}

	return;
}

void NonPlaidFingerprint::change_CoM_to_ligandCoM(numeric::xyzVector<core::Real> const & ligandCoM) {

	CoM_ = ligandCoM;
	return;

}
/*
void NonPlaidFingerprint::setup_from_eggshell_triplet_file(std::string const & input_filename) {
//input_file should contain a line starting with "/COM/" for CoM values
//input_file should contain a line starting with "/ORI/" for ORIGIN values
//input_file should contain 3 tab separated columns for the triplet

ifstream inFile(input_filename.c_str());
if (!inFile) {
std::cout<< "Can't open input file " << input_filename << std::endl;
exit(1);
}

std::string lineread;
std::string Line;
std::string Field;

spherical_coor_triplet new_triplet;
std::list< spherical_coor_triplet> ray_triplet_list;
ray_triplet_list.clear();
NonPlaidFingerprint pair_list;
numeric::xyzVector<core::Real> new_origin;
numeric::xyzVector<core::Real> new_CoM;

while (std::getline(inFile, lineread)) {

std::stringstream sss(lineread);
std::string Pock_string_phi, Pock_string_psi, Pock_string_rho;
core::Real Pock_real_phi, Pock_real_psi, Pock_real_rho;

//start new ray_triplet_list if the line has the "TER" field
if (lineread[0] == 'T' && lineread[1] == 'E' && lineread[2] == 'R' ){
pair_list.origin_ = new_origin;
pair_list.CoM_ = new_CoM;
pair_list.triplet_fingerprint_data_ = ray_triplet_list;
npf_list_.push_back(pair_list);
triplet_fingerprint_data_ = ray_triplet_list;
ray_triplet_list.clear();
continue;
}

//parse CoM values from line starting with "/COM/"
if (lineread[0] == '/' && lineread[1] == 'C' && lineread[2] == 'O' && lineread[3] == 'M' && lineread[4] == '/') {
lineread.erase(0,5);
std::stringstream com_line(lineread);
std::getline(com_line, Pock_string_phi, '\t');
new_CoM.x() = atof(Pock_string_phi.c_str());
std::getline(com_line, Pock_string_psi, '\t');
new_CoM.y() = atof(Pock_string_psi.c_str());
std::getline(com_line, Pock_string_rho, '\t');
new_CoM.z() = atof(Pock_string_rho.c_str());
//std::cout<<"setup from triplet file : "<< " " <<CoM_.x()<<" "<<CoM_.y()<<" "<<CoM_.z()<<std::endl;
continue;
}

//parse GRID dimension values from line starting with "/DIM/"
if (lineread[0] == '/' && lineread[1] == 'D' && lineread[2] == 'I' && lineread[3] == 'M' && lineread[4] == '/') {
lineread.erase(0,5);
std::stringstream dim_line(lineread);
std::getline(dim_line, Pock_string_phi, '\t');
pocketGrid_dim_.x() = atof(Pock_string_phi.c_str());
std::getline(dim_line, Pock_string_psi, '\t');
pocketGrid_dim_.y() = atof(Pock_string_psi.c_str());
std::getline(dim_line, Pock_string_rho, '\t');
pocketGrid_dim_.z() = atof(Pock_string_rho.c_str());
//std::cout<<"setup from triplet file : "<< " " <<CoM_.x()<<" "<<CoM_.y()<<" "<<CoM_.z()<<std::endl;
continue;
}

//parse GRID midpoint values from line starting with "/MID/"
if (lineread[0] == '/' && lineread[1] == 'M' && lineread[2] == 'I' && lineread[3] == 'D' && lineread[4] == '/') {
lineread.erase(0,5);
std::stringstream mid_line(lineread);
std::getline(mid_line, Pock_string_phi, '\t');
pocketGrid_mid_.x() = atof(Pock_string_phi.c_str());
std::getline(mid_line, Pock_string_psi, '\t');
pocketGrid_mid_.y() = atof(Pock_string_psi.c_str());
std::getline(mid_line, Pock_string_rho, '\t');
pocketGrid_mid_.z() = atof(Pock_string_rho.c_str());
//std::cout<<"setup from triplet file : "<< " " <<CoM_.x()<<" "<<CoM_.y()<<" "<<CoM_.z()<<std::endl;
continue;
}

//parse GRID spacing value from line starting with "/SPACING/"
if (lineread[0] == '/' && lineread[1] == 'S' && lineread[2] == 'P' && lineread[3] == 'A' && lineread[4] == 'C' && lineread[5] == 'I' && lineread[6] == 'N' && lineread[7] == 'G' && lineread[8] == '/') {
lineread.erase(0,9);
std::stringstream spacing_line(lineread);
std::getline(spacing_line, Pock_string_phi);
pocketGrid_spacing_ = atof(Pock_string_phi.c_str());
//std::cout<<"setup from triplet file : "<< " " <<CoM_.x()<<" "<<CoM_.y()<<" "<<CoM_.z()<<std::endl;
continue;
}

//parse ORIGIN values from line starting with "/ORI/"
if (lineread[0] == '/' && lineread[1] == 'O' && lineread[2] == 'R' && lineread[3] == 'I' && lineread[4] == '/') {
lineread.erase(0,5);
std::stringstream ori_line(lineread);
std::getline(ori_line, Pock_string_phi, '\t');
new_origin.x() = atof(Pock_string_phi.c_str());
std::getline(ori_line, Pock_string_psi, '\t');
new_origin.y() = atof(Pock_string_psi.c_str());
std::getline(ori_line, Pock_string_rho, '\t');
new_origin.z() = atof(Pock_string_rho.c_str());
//std::cout<<"setup from triplet file : "<< " " <<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
continue;
}

std::getline(sss, Pock_string_phi, '\t');
Pock_real_phi = atof(Pock_string_phi.c_str());
std::getline(sss, Pock_string_psi, '\t');
Pock_real_psi = atof(Pock_string_psi.c_str());
std::getline(sss, Pock_string_rho, '\t');
Pock_real_rho = atof(Pock_string_rho.c_str());

new_triplet.phi = Pock_real_phi;
new_triplet.psi = Pock_real_psi;
new_triplet.rho = Pock_real_rho;
ray_triplet_list.push_back(new_triplet);

}
inFile.close();

//set CoM_
CoM_ = new_CoM;

return;
}
*/
void NonPlaidFingerprint::setup_from_eggshell_triplet_file(std::string const & input_filename) {

	//input_file should contain a line starting with "/COM/" for CoM values
	//input_file should contain a line starting with "/ORI/" for ORIGIN values
	//input_file should contain 3 tab separated columns for the triplet

	ifstream inFile(input_filename.c_str());

	if ( !inFile ) {
		std::cout<< "Can't open input file " << input_filename << std::endl;
		exit(1);
	}

	std::string lineread;
	//std::string Line;
	//std::string Field;

	triplet_fingerprint_data_.clear();
	multi_origin_list_.clear();
	while ( std::getline(inFile, lineread) ) {

		//std::stringstream sss(lineread);
		std::string Pock_string_phi, Pock_string_psi, Pock_string_rho, Pock_string_ori;
		//core::Real Pock_real_phi, Pock_real_psi, Pock_real_rho;
		//core::Size Pock_real_ori;

		//parse GRID dimension values from line starting with "/DIM/"
		if ( lineread[0] == '/' && lineread[1] == 'D' && lineread[2] == 'I' && lineread[3] == 'M' && lineread[4] == '/' && lineread[5] == ' ' ) {
			lineread.erase(0,6);
			std::stringstream dim_line(lineread);
			std::getline(dim_line, Pock_string_phi, '\t');
			pocketGrid_dim_.x() = atof(Pock_string_phi.c_str());
			std::getline(dim_line, Pock_string_psi, '\t');
			pocketGrid_dim_.y() = atof(Pock_string_psi.c_str());
			std::getline(dim_line, Pock_string_rho, '\t');
			pocketGrid_dim_.z() = atof(Pock_string_rho.c_str());
			//std::cout<<"setup from triplet file dim : "<< " " <<pocketGrid_dim_.x()<<" "<<pocketGrid_dim_.y()<<" "<<pocketGrid_dim_.z()<<std::endl;
			continue;
		}

		//parse GRID midpoint values from line starting with "/MID/"
		if ( lineread[0] == '/' && lineread[1] == 'M' && lineread[2] == 'I' && lineread[3] == 'D' && lineread[4] == '/' && lineread[5] == ' ' ) {
			lineread.erase(0,6);
			std::stringstream mid_line(lineread);
			std::getline(mid_line, Pock_string_phi, '\t');
			pocketGrid_mid_.x() = atof(Pock_string_phi.c_str());
			std::getline(mid_line, Pock_string_psi, '\t');
			pocketGrid_mid_.y() = atof(Pock_string_psi.c_str());
			std::getline(mid_line, Pock_string_rho, '\t');
			pocketGrid_mid_.z() = atof(Pock_string_rho.c_str());
			//std::cout<<"setup from triplet file mid : "<< " " <<pocketGrid_mid_.x()<<" "<<pocketGrid_mid_.y()<<" "<<pocketGrid_mid_.z()<<std::endl;
			continue;
		}

		//parse GRID spacing value from line starting with "/SPACING/"
		if ( lineread[0] == '/' && lineread[1] == 'S' && lineread[2] == 'P' && lineread[3] == 'A' && lineread[4] == 'C' && lineread[5] == 'I' && lineread[6] == 'N' && lineread[7] == 'G' && lineread[8] == '/' && lineread[9] == ' ' ) {
			lineread.erase(0,10);
			std::stringstream spacing_line(lineread);
			std::getline(spacing_line, Pock_string_phi);
			pocketGrid_spacing_ = atof(Pock_string_phi.c_str());
			//    std::cout<<"setup from triplet file spacing : "<< " " <<pocketGrid_spacing_<<std::endl;
			continue;
		}

		//parse ORIGIN values from line starting with "/ORI/"
		if ( lineread[0] == '/' && lineread[1] == 'O' && lineread[2] == 'R' && lineread[3] == 'I' && lineread[4] == '/' && lineread[5] == ' ' ) {
			lineread.erase(0,6);
			std::stringstream ori_line(lineread);
			std::getline(ori_line, Pock_string_phi, '\t');
			origin_.x() = atof(Pock_string_phi.c_str());
			std::getline(ori_line, Pock_string_psi, '\t');
			origin_.y() = atof(Pock_string_psi.c_str());
			std::getline(ori_line, Pock_string_rho, '\t');
			origin_.z() = atof(Pock_string_rho.c_str());
			//std::getline(ori_line, Pock_string_ori, '\t');
			//core::Size inx = core::Size(atoi(Pock_string_ori.c_str()));
			//std::cout<<"setup from triplet file : "<< " " <<inx<<std::endl;
			multi_origin_list_.push_back(origin_);
			//std::cout<<"setup from triplet file : "<< " " <<inx<<std::endl;
			//std::cout<<"setup from triplet file ori : "<< " " <<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
			continue;
		}
		//parse CoM values from line starting with "/COM/"
		if ( lineread[0] == '/' && lineread[1] == 'C' && lineread[2] == 'O' && lineread[3] == 'M' && lineread[4] == '/' && lineread[5] == ' ' ) {
			lineread.erase(0,6);
			std::stringstream com_line(lineread);
			std::getline(com_line, Pock_string_phi, '\t');
			CoM_.x() = atof(Pock_string_phi.c_str());
			std::getline(com_line, Pock_string_psi, '\t');
			CoM_.y() = atof(Pock_string_psi.c_str());
			std::getline(com_line, Pock_string_rho, '\t');
			CoM_.z() = atof(Pock_string_rho.c_str());
			//std::cout<<"setup from triplet file : "<< " " <<CoM_.x()<<" "<<CoM_.y()<<" "<<CoM_.z()<<std::endl;
			continue;
		}
		if ( lineread[0] == '/' && lineread[1] == 'R' && lineread[2] == 'A' && lineread[3] == 'Y' && lineread[4] == '/' && lineread[5] == ' ' ) {
			lineread.erase(0,6);
			std::stringstream ray_line(lineread);
			spherical_coor_triplet new_triplet;
			std::getline(ray_line, Pock_string_phi, '\t');
			new_triplet.phi = atof(Pock_string_phi.c_str());
			std::getline(ray_line, Pock_string_psi, '\t');
			new_triplet.psi = atof(Pock_string_psi.c_str());
			std::getline(ray_line, Pock_string_rho, '\t');
			new_triplet.rho = atof(Pock_string_rho.c_str());
			std::getline(ray_line, Pock_string_ori, '\t');
			new_triplet.ori = atoi(Pock_string_ori.c_str());
			triplet_fingerprint_data_.push_back(new_triplet);
			//std::cout<<"triplet file : "<< " " <<new_triplet.ori<<std::endl;
		}
	}

	inFile.close();
	num_origins_ = multi_origin_list_.size();

#ifdef USEOPENCL
		gpu_setup_rays();
#endif

	return;
}

void NonPlaidFingerprint::trim_based_on_known_ligand(core::pose::Pose const & known_ligand_pose){

	protocols::pockets::PlaidFingerprint known_pf( known_ligand_pose, *this );
	std::list< spherical_coor_triplet > triplet_trim_data;
	for ( std::list<spherical_coor_triplet>::const_iterator pro = triplet_fingerprint_data_.begin(), lig = known_pf.triplet_fingerprint_data().begin(); ( pro != triplet_fingerprint_data_.end()) && (lig != known_pf.triplet_fingerprint_data().end()); ++pro, ++lig ) {

		//jk note: these are no longer necessarily true, since we alter the Plaid one by shifting up/down by 2*pi
		//these are useful asserts though, it's worth thinking about how to make them valid again...
		//debug_assert( std::abs( pro->phi - lig->phi ) < 0.001 );
		//debug_assert( std::abs( pro->psi - lig->psi ) < 0.001 );

		if ( ( pro->rho > 0.001 ) && ( lig->rho < 0.001 ) ) continue;
		triplet_trim_data.push_back(*pro);

	}
	triplet_fingerprint_data_.clear();
	triplet_fingerprint_data_ = triplet_trim_data;

	//DUMP TRIMED_EGGSHELL TO A PDB FILE
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("trim_eggshell.pdb", std::ios::out);
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()<<std::endl;
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM B   2    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()<<std::endl;
	for ( std::list<spherical_coor_triplet>::const_iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd ) {
		numeric::xyzVector<core::Real> new_coor;
		convert_spherical_coor_triplet_to_cartesian( *pd, new_coor );
		new_coor += origin_;
		outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG C   3    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.z()<<std::endl;
	}
	outPDB_stream.close();
	outPDB_stream.clear();

	return;
}

void NonPlaidFingerprint::include_eggshell_points_based_on_known_ligand( core::pose::Pose const & known_ligand_pose, core::Real const & trim_dist) {
	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = known_ligand_pose.size(); j <= resnum; ++j ) {
		if ( !known_ligand_pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		std::cout<<"Error, no ligand to include_eggshell_points_based_on_known_ligand" << std::endl;
		exit(1);
	}

	core::conformation::Residue const & curr_rsd = known_ligand_pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms = curr_rsd.nheavyatoms();
	numeric::xyzVector<core::Real> lig_atom_coord;
	std::list< numeric::xyzVector<core::Real> > lig_atom_coord_list;
	for ( Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i ) {
		lig_atom_coord.x() = curr_rsd.atom(i).xyz()(1);
		lig_atom_coord.y() = curr_rsd.atom(i).xyz()(2);
		lig_atom_coord.z() = curr_rsd.atom(i).xyz()(3);
		lig_atom_coord_list.push_back(lig_atom_coord);
	}

	std::list< numeric::xyzVector<core::Real> > new_egg_coord_list;
	std::list< numeric::xyzVector<core::Real> > new_ext_coord_list;
	numeric::xyzVector<core::Real> xyz_coord;

	new_egg_coord_list.clear();
	new_ext_coord_list.clear();
	for ( std::list< numeric::xyzVector<core::Real> >::const_iterator aa = eggshell_list_.begin(); aa != eggshell_list_.end(); ++aa ) {
		xyz_coord = *aa;
		bool found = false;
		for ( std::list< numeric::xyzVector<core::Real> >::const_iterator bb = lig_atom_coord_list.begin(); bb != lig_atom_coord_list.end(); ++bb ) {
			if ( xyz_coord.distance(*bb) <= trim_dist ) { found = true;break;}
		}
		if ( found ) {
			new_egg_coord_list.push_back(xyz_coord);
		}
		//  else{
		// new_ext_coord_list.push_back(xyz_coord);
		// }
	}
	eggshell_list_.clear();
	eggshell_list_ = new_egg_coord_list;
	// extshell_list_.clear();
	//extshell_list_ = new_ext_coord_list;

	return;
}

void NonPlaidFingerprint::set_origin_away_from_eggshell_plane( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose, Size const & set_origin_option ) {

	origin_.zero();

	core::Real A11(0.), A12(0.), A13(0.), A22(0.), A23(0.), b1(0.), b2(0.), b3(0.);
	for ( auto const & pd : egg_and_extra_shell ) {
		A11 += pd.x()*pd.x();
		A12 += pd.x()*pd.y();
		A13 += pd.x();
		A22 += pd.y()*pd.y();
		A23 += pd.y();
		b1 += pd.x()*pd.z();
		b2 += pd.y()*pd.z();
		b3 += pd.z();
	}

	numeric::xyzMatrix<core::Real> A;
	A(1,1)=A11;
	A(1,2)=A12;
	A(2,1)=A12;
	A(1,3)=A13;
	A(3,1)=A13;
	A(2,2)=A22;
	A(2,3)=A23;
	A(3,2)=A23;
	A(3,3)=egg_and_extra_shell.size();
	numeric::xyzVector<core::Real> b;
	b(1)=b1;
	b(2)=b2;
	b(3)=b3;

	numeric::xyzVector<core::Real> soln = inverse(A) * b;
	core::Real best_fit_a = soln.x();
	core::Real best_fit_b = soln.y();
	core::Real best_fit_c = soln.z();

	// best fit plane has the equation z = best_fit_a * x + best_fit_b * y + best_fit_c

	// RAGUL TEST
	numeric::xyzVector<core::Real> plane_coord;
	std::list< numeric::xyzVector<core::Real> > plane_coord_list;

	for ( auto const & pdebug : egg_and_extra_shell ) {
		plane_coord.x() = pdebug.x();
		plane_coord.y() = pdebug.y();
		plane_coord.z() = best_fit_a * pdebug.x() + best_fit_b * pdebug.y() + best_fit_c;
		plane_coord_list.push_back(plane_coord);
	}

	//plane_CoM is the same as eggshell_CoM
	numeric::xyzVector<core::Real> plane_CoM(0.);
	plane_CoM = CoM_;

	//cross product of two vectors that lies on the plane gives a vector perpendicular to the plane
	numeric::xyzVector<core::Real> plane_normal(0.);
	plane_normal = cross_product( (plane_coord_list.front() - plane_CoM ), (plane_coord_list.back() - plane_CoM ) );
	// set plane_normal to have length 30 A
	plane_normal.normalize(30.);

	numeric::xyzVector<core::Real> temp_origin1 = plane_CoM + plane_normal;
	numeric::xyzVector<core::Real> temp_origin2 = plane_CoM - plane_normal;
	numeric::xyzVector<core::Real> closest_origin(0.);
	numeric::xyzVector<core::Real> distant_origin(0.);

	//calculate distance b/w protein_CoM and +/- origin and choose the shortest
	numeric::xyzVector<core::Real> protein_CoM = calculate_protein_CoM(protein_pose);
	if ( protein_CoM.distance(temp_origin1) > protein_CoM.distance(temp_origin2) ) {
		closest_origin = temp_origin2;
		distant_origin = temp_origin1;
	} else {
		closest_origin = temp_origin1;
		distant_origin = temp_origin2;
	}

	if ( set_origin_option == 3 ) {
		origin_ = closest_origin;
	} else if ( set_origin_option == 4 ) {
		origin_ = distant_origin;
	}

	//DUMP EGGSHELL_GRID(plane) TO A PDB FILE
	/*
	utility::io::ozstream outPDB_stream;
	outPDB_stream.open("plane.pdb", std::ios::out);
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = plane_coord_list.begin(); pd != plane_coord_list.end(); ++pd) {
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   PLN A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
	}
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM C   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<plane_CoM.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<plane_CoM.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<plane_CoM.z()<<std::endl;
	//outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORA A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin1.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin1.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin1.z()<<std::endl;
	//outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORB B   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin2.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin2.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<temp_origin2.z()<<std::endl;
	outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI C   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()<<std::endl;
	outPDB_stream.close();
	outPDB_stream.clear();
	*/
	//END printing plane.pdb

	return;

}

void NonPlaidFingerprint::set_origin_away_from_eggshell(std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose ) {

	// set origin_ 30 A below the eggshell
	origin_.zero();
	numeric::xyzVector<core::Real> temp_vec(0.);
	for ( auto const & pd : egg_and_extra_shell ) {
		temp_vec += (CoM_ - pd);
	}
	temp_vec.normalize(30.);

	numeric::xyzVector<core::Real> temp_origin1 = CoM_ + temp_vec;
	numeric::xyzVector<core::Real> temp_origin2 = CoM_ - temp_vec;
	numeric::xyzVector<core::Real> closest_origin(0.);
	numeric::xyzVector<core::Real> distant_origin(0.);

	//calculate distance b/w protein_CoM and +/- origin and choose the shortest
	numeric::xyzVector<core::Real> protein_CoM = calculate_protein_CoM(protein_pose);
	if ( protein_CoM.distance(temp_origin1) > protein_CoM.distance(temp_origin2) ) {
		closest_origin = temp_origin2;
		distant_origin = temp_origin1;
	} else {
		closest_origin = temp_origin1;
		distant_origin = temp_origin2;
	}
	origin_ = closest_origin;
	return;
}


void NonPlaidFingerprint::set_origin_away_from_protein_center(core::pose::Pose const & protein_pose) {

	// set origin_ to the protein CoM, then move origin_ 30Angstrom away from protein center
	origin_ = calculate_protein_CoM(protein_pose);
	numeric::xyzVector<core::Real> temp_vec(0.);
	temp_vec = origin_ - CoM_;
	temp_vec.normalize(30.);
	origin_ = temp_vec + CoM_;

	return;
}

void NonPlaidFingerprint::set_origin_from_residue(core::pose::Pose const & protein_pose) {

	// set origin_ to the given residue CoM
	using namespace basic::options;
	std::string const resid(option[ OptionKeys::fingerprint::origin_res_num ]);
	int pdb_number;
	char chain = ' ';
	std::size_t fpos( resid.find(':') );
	if ( fpos != std::string::npos ) {
		pdb_number = ObjexxFCL::int_of( resid.substr(0,fpos) );
		if ( fpos != resid.size()-1 ) {
			chain = resid[ fpos+1 ];
		}
	} else {
		pdb_number = ObjexxFCL::int_of( resid );
	}
	core::Size rosetta_pose_pdb_number = get_pose_resnum(pdb_number, chain, protein_pose);
	core::conformation::Residue const & curr_rsd  = protein_pose.residue(rosetta_pose_pdb_number);
	numeric::xyzVector<core::Real> residue_com(0.);
	for ( Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i ) {
		residue_com.x() += curr_rsd.atom(i).xyz()(1);
		residue_com.y() += curr_rsd.atom(i).xyz()(2);
		residue_com.z() += curr_rsd.atom(i).xyz()(3);
	}
	residue_com /= curr_rsd.nheavyatoms();
	numeric::xyzVector<core::Real> temp_vec(0.);
	temp_vec = residue_com - CoM_;
	temp_vec.normalize(30.);
	origin_ = temp_vec + CoM_;

	return;
}

core::Size NonPlaidFingerprint::get_pose_resnum(int const pdbnum, char const pdbchn, core::pose::Pose const & ps) {

	for ( core::Size j = 1; j <= ps.size(); ++j ) {
		if ( ( ps.pdb_info()->chain(j) == pdbchn ) && (ps.pdb_info()->number(j) == pdbnum) ) {
			return j;
		}
	}
	// residue not found
	std::cout << "ERROR!! Could not find residue" << pdbnum << " and chain " << pdbchn << std::endl;
	exit(1);
}


numeric::xyzVector<core::Real> NonPlaidFingerprint::calculate_protein_CoM(core::pose::Pose const & protein_pose) {

	numeric::xyzVector<core::Real> protein_com(0.);
	core::Size total_atoms(0);
	for ( int j = 1, resnum = protein_pose.size(); j <= resnum; ++j ) {
		core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
		if ( curr_rsd.is_protein() ) {
			for ( Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i ) {
				protein_com.x() += curr_rsd.atom(i).xyz()(1);
				protein_com.y() += curr_rsd.atom(i).xyz()(2);
				protein_com.z() += curr_rsd.atom(i).xyz()(3);
				total_atoms++;
			}
		}
	}
	protein_com /= total_atoms;
	return protein_com;

}

PlaidFingerprint::PlaidFingerprint( core::pose::Pose const & input_pose, FingerprintBase & fp ) :
	FingerprintBase(),
	pose_(input_pose)
{

	//calculate ligand center of mass
	numeric::xyzVector<core::Real> ligand_CoM = calculate_ligand_CoM(pose_);

	//move pose_ such that ligand CoM is at the pocket center of mass
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	pose_.apply_transform_Rx_plus_v(I_mat, fp.CoM() - ligand_CoM);

	core::Size const lig_res_num = compute_ligand_resnum();
	core::conformation::ResidueCOP ligand_rsd( core::conformation::ResidueOP( new core::conformation::Residue( pose_.conformation().residue(lig_res_num) ) ) );

	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	select_conf_and_move_ligand_(fp,no_CoM_offset,0,0,0,0);
	update_rhos_(fp, ligand_rsd);
}

void PlaidFingerprint::apply_rotation_offset_to_pose_( core::pose::Pose & pose, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) const {

	//calculate_ligand_CoM function not working here (const )!!
	//numeric::xyzVector<core::Real> ligand_CoM = calculate_ligand_CoM(pose);

	core::Size lig_res_num = compute_ligand_resnum(pose);
	numeric::xyzVector<core::Real> ligand_CoM(0.);
	core::conformation::Residue const & curr_rsd = pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms = compute_ligand_natoms(pose);
	for ( Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i ) {
		ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
		ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
		ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
	}
	ligand_CoM /= ligand_total_atoms;

	// move pose such that ligand CoM is at the origin
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	pose.apply_transform_Rx_plus_v(I_mat, -ligand_CoM);

	numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_radians(angle1_offset) );
	numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_radians(angle2_offset) );
	numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_radians(angle3_offset) );
	numeric::xyzMatrix<core::Real> tot_rot_mat = z_rot_mat * y_rot_mat * x_rot_mat;
	core::Vector v(0,0,0);
	pose.apply_transform_Rx_plus_v(tot_rot_mat, v);

	// move pose back to original ligand CoM
	pose.apply_transform_Rx_plus_v(I_mat, ligand_CoM);

}

core::Size PlaidFingerprint::compute_ligand_resnum( core::pose::Pose const & pose ) const {
	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = pose.size(); j <= resnum; ++j ) {
		if ( !pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}
	return lig_res_num;
}

core::Size PlaidFingerprint::compute_ligand_natoms( core::pose::Pose const & pose ) const {
	core::Size lig_res_num = compute_ligand_resnum(pose);
	core::conformation::Residue const & curr_rsd = pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms;
	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::include_hydrogens ]() ) {
		ligand_total_atoms = curr_rsd.natoms();
	} else {
		ligand_total_atoms = curr_rsd.nheavyatoms();
	}
	return ligand_total_atoms;
}

core::Size PlaidFingerprint::compute_ligand_natoms_with_hydrogens( core::pose::Pose const & pose ) const {
	core::Size lig_res_num = compute_ligand_resnum(pose);
	core::conformation::Residue const & curr_rsd = pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms;
	ligand_total_atoms = curr_rsd.natoms();
	return ligand_total_atoms;
}

core::Size PlaidFingerprint::compute_ligand_nconformers( core::pose::Pose const & pose ) const {
	return pose.size();
}

core::conformation::ResidueCOP PlaidFingerprint::select_conf_and_move_ligand_( FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer ) {

	// note: conformer passed in is indexed from 0, below is indexed from 1
	core::pose::PoseOP tmp_pose( new core::pose::Pose(pose_, conformer+1, conformer+1) );

	apply_rotation_offset_to_pose_( *tmp_pose, angle1_offset, angle2_offset, angle3_offset );

	// set the fingerprint CoM to be the ligand CoM, then apply offset as needed
	CoM_ = calculate_ligand_CoM(*tmp_pose);

	// jk note: this next step of setting the origins to match is unnecessary as written, because given current implementation CoM_ and fp.CoM() are the same
	// jk is there a case when this isn't true, or is this leftover from an old approach??
	origin_ = fp.origin() + CoM_ - fp.CoM();
	origin_ += CoM_offset;

	num_origins_ = fp.num_origins();
	multi_origin_list_ = fp.multi_origin_list();
	for ( core::Size i=1; i <=fp.multi_origin_list().size(); ++i ) {
		multi_origin_list_[i] += CoM_offset;
	}

	core::Size const lig_res_num = compute_ligand_resnum(*tmp_pose);
	// core::Size const ligand_natoms = compute_ligand_natoms(*tmp_pose);
	core::conformation::ResidueCOP ligand_rsd( core::conformation::ResidueOP( new core::conformation::Residue( tmp_pose->conformation().residue(lig_res_num) ) ) );

	return ligand_rsd;

}
void PlaidFingerprint::update_rhos_(FingerprintBase & fp, core::conformation::ResidueCOP curr_ligand_rsd, bool const update_derivatives ) {

	using namespace basic::options;
	core::Real const radius_scale = option[ OptionKeys::fingerprint::atom_radius_scale ];
	core::Real const atom_buffer = option[ OptionKeys::fingerprint::atom_radius_buffer ];

	triplet_fingerprint_data_.clear();
	core::Size const ligand_natoms = compute_ligand_natoms();

	// compute the max and min possible phi and psi for each atom, store these in 4 arrays
	utility::vector1< utility::vector1<core::Real> > ori_atom_max_phi(num_origins_);
	utility::vector1< utility::vector1<core::Real> > ori_atom_max_psi(num_origins_);
	utility::vector1< utility::vector1<core::Real> > ori_atom_min_phi(num_origins_);
	utility::vector1< utility::vector1<core::Real> > ori_atom_min_psi(num_origins_);
	utility::vector1< utility::vector1< core::Real > > ori_atomX(num_origins_);
	utility::vector1< utility::vector1< core::Real > > ori_atomY(num_origins_);
	utility::vector1< utility::vector1< core::Real > > ori_atomZ(num_origins_);
	utility::vector1< utility::vector1< core::Real > > ori_atom_radius(num_origins_);

	int origin_num = 1;
	for ( utility::vector1<numeric::xyzVector<core::Real> >::const_iterator i_mori = multi_origin_list_.begin(); i_mori != multi_origin_list_.end(); ++i_mori, ++origin_num ) {
		utility::vector1<core::Real> atom_max_phi( ligand_natoms, 0 );
		utility::vector1<core::Real> atom_max_psi( ligand_natoms, 0 );
		utility::vector1<core::Real> atom_min_phi( ligand_natoms, 0 );
		utility::vector1<core::Real> atom_min_psi( ligand_natoms, 0 );
		utility::vector1< core::Real > atomX( ligand_natoms, 0 );
		utility::vector1< core::Real > atomY( ligand_natoms, 0 );
		utility::vector1< core::Real > atomZ( ligand_natoms, 0 );
		utility::vector1<core::Real> atom_radius( ligand_natoms, 0 );

		core::Real max_phi( -999. );
		core::Real max_psi( -999. );
		core::Real min_phi( 999. );
		core::Real min_psi( 999. );

		for ( Size i = 1, i_end = ligand_natoms; i <= i_end; ++i ) {

			numeric::xyzVector<core::Real> this_atomcoors = curr_ligand_rsd->atom(i).xyz() - *i_mori;
			core::Real const this_atom_radius = ( curr_ligand_rsd->atom_type(i).lj_radius() - atom_buffer ) * radius_scale;

			// find the atom center (in spherical coors)
			spherical_coor_triplet atom_center;
			convert_cartesian_to_spherical_coor_triplet( this_atomcoors, atom_center );
			core::Real const tmp_atomx=this_atomcoors.x();
			core::Real const tmp_atomy=this_atomcoors.y();
			core::Real const tmp_atomz=this_atomcoors.z();

			// find the max/min angular phi displacement that will intersect the atom
			// project atom onto the z-axis (ie. set x and y coors to zero) to calculate this for phi
			core::Real curr_max_phi( 999 );
			core::Real curr_min_phi( -999 );
			if ( std::abs(tmp_atomz) > 0.00001 ) {
				core::Real const inside_phi_asin = this_atom_radius / tmp_atomz;
				if ( std::abs(inside_phi_asin) < 1. ) {
					core::Real const max_angular_displacement_phi = std::abs( asin( inside_phi_asin ) );
					curr_max_phi = atom_center.phi + max_angular_displacement_phi;
					curr_min_phi = atom_center.phi - max_angular_displacement_phi;
				}
			}

			// find the max/min angular psi displacement that will intersect the atom
			// project atom onto the x-y plane (ie. set the z coor to zero) to calculate this for psi
			core::Real curr_max_psi( 999 );
			core::Real curr_min_psi( -999 );
			if ( ( std::abs(tmp_atomx) > 0.00001 ) ) { //|| ( std::abs(tmp_atomx) > 0.00001 ) ) {
				core::Real const inside_psi_asin = this_atom_radius / sqrt( (tmp_atomx*tmp_atomx) + (tmp_atomy*tmp_atomy) );
				if ( std::abs(inside_psi_asin) < 1. ) {
					core::Real const max_angular_displacement_psi = std::abs( asin( inside_psi_asin ) );
					curr_max_psi = atom_center.psi + max_angular_displacement_psi;
					curr_min_psi = atom_center.psi - max_angular_displacement_psi;
				}
			}

			if ( curr_max_phi > max_phi ) max_phi = curr_max_phi;
			if ( curr_max_psi > max_psi ) max_psi = curr_max_psi;
			if ( curr_min_phi < min_phi ) min_phi = curr_min_phi;
			if ( curr_min_psi < min_psi ) min_psi = curr_min_psi;

			atom_max_phi.at(i) = curr_max_phi;
			atom_max_psi.at(i) = curr_max_psi;
			atom_min_phi.at(i) = curr_min_phi;
			atom_min_psi.at(i) = curr_min_psi;
			atomX.at(i) = this_atomcoors.x();
			atomY.at(i) = this_atomcoors.y();
			atomZ.at(i) = this_atomcoors.z();
			atom_radius.at(i) = this_atom_radius;
		}
		ori_atom_max_phi[origin_num] = atom_max_phi;
		ori_atom_max_psi[origin_num] = atom_max_psi;
		ori_atom_min_phi[origin_num] = atom_min_phi;
		ori_atom_min_psi[origin_num] = atom_min_psi;
		ori_atomX[origin_num] = atomX;
		ori_atomY[origin_num] = atomY;
		ori_atomZ[origin_num] = atomZ;
		ori_atom_radius[origin_num] = atom_radius;
	}

	derivs_of_ray_distances_.clear();

	//float orig_cpu_total_dist = 0.;
	//Size orig_cpu_total_num = 0;
	//Size orig_cpu_num_evaluations = 0;

	for ( auto const & ni : fp.triplet_fingerprint_data() ) {

		core::Real curr_phi = ni.phi;
		core::Real curr_psi = ni.psi;
		core::Size curr_ori = ni.ori;
		// if the current phi and/or psi is outside the overall max/min, set best_rho to zero and jumpout (ie. ray misses the ligand)
		core::Real best_rho_sq(9999.);
		//core::Size best_intersecting_atom(0);
		for ( Size i = 1, i_end = ligand_natoms; i <= i_end; ++i ) {
			if ( ori_atom_radius[curr_ori][i] < 0.001 ) continue;
			while ( curr_phi < ori_atom_min_phi[curr_ori][i] ) {
				curr_phi += numeric::constants::r::pi_2;
			}

			while ( curr_phi > ori_atom_max_phi[curr_ori][i] ) {
				curr_phi -= numeric::constants::r::pi_2;
			}
			while ( curr_psi < ori_atom_min_psi[curr_ori][i] ) {
				curr_psi += numeric::constants::r::pi_2;
			}
			while ( curr_psi > ori_atom_max_psi[curr_ori][i] ) {
				curr_psi -= numeric::constants::r::pi_2;
			}
			if ( curr_phi < ori_atom_min_phi[curr_ori][i] ) continue;
			if ( curr_psi < ori_atom_min_psi[curr_ori][i] ) continue;
			if ( curr_phi > ori_atom_max_phi[curr_ori][i] ) continue;
			if ( curr_psi > ori_atom_max_psi[curr_ori][i] ) continue;

			core::Real const min_intersect_SQ = Find_Closest_Intersect_SQ(curr_phi,curr_psi,ori_atomX[curr_ori][i],ori_atomY[curr_ori][i],ori_atomZ[curr_ori][i],ori_atom_radius[curr_ori][i]);

			//++orig_cpu_num_evaluations;
			if ( min_intersect_SQ < best_rho_sq ) {
				best_rho_sq = min_intersect_SQ;
				//best_intersecting_atom = i;
			}
		}
		spherical_coor_triplet new_triplet;
		new_triplet.phi = curr_phi;
		new_triplet.psi = curr_psi;
		new_triplet.ori = curr_ori;

		if ( best_rho_sq < 9998. ) {
			new_triplet.rho = sqrt(best_rho_sq);
		} else {
			new_triplet.rho = 0.;
		}
		triplet_fingerprint_data_.push_back(new_triplet);

		if ( update_derivatives ) {
			std::cout<<"DARC derivatives are currently commented out...." << std::endl;
			exit(1);
		}

		/*

		if ( update_derivatives ) {

		ray_distance_derivs new_deriv_set;

		// compute derivatives only if the ray hits both pocket and ligand - otherwise derivative is zero
		if ( (ni->rho > 0.001) && (new_triplet.rho > 0.001) ) {

		//   (x0,y0,z0) = the coordinates of the "origin" point from which all rays originate, called P0
		core::Real const x0 = origin_.x();
		core::Real const y0 = origin_.y();
		core::Real const z0 = origin_.z();

		//   (x2,y2,z2) = the coordinates of the pocket point that intersects the ray whose derivative we are calculating, called P2
		numeric::xyzVector<core::Real> fp_coord;
		convert_spherical_coor_triplet_to_cartesian( *ni, fp_coord );
		fp_coord += origin_;
		core::Real const x2 = fp_coord.x();
		core::Real const y2 = fp_coord.y();
		core::Real const z2 = fp_coord.z();

		//   (x3,y3,z3) = the coordinates of the center of the atom that is intersecting the ray, called P3
		// note: can't use atomX / atomY / atomZ, since these are relative to the origin
		numeric::xyzVector<core::Real> best_atom_coors = curr_rsd.atom(best_intersecting_atom).xyz();
		core::Real const x3 = best_atom_coors.x();
		core::Real const y3 = best_atom_coors.y();
		core::Real const z3 = best_atom_coors.z();

		//   (xc,yc,zc) = the coordinates of the center of mass of the ligand BEFORE it has been moved into the current pose, called Pc
		// note: before moving ligand, the ligand CoM (in pose_) is at fp.CoM() - this happens in the PlaidFingerprint constructor
		core::Real const xc = fp.CoM().x();
		core::Real const yc = fp.CoM().y();
		core::Real const zc = fp.CoM().z();

		//   (xs,ys,zs) = the coordinates of the center of the atom that is intersecting the ray BEFORE it has been moved into the current pose, called Ps
		conformation::Residue const & ligand_res_before_rotation = pose_.conformation().residue(lig_res_num);
		numeric::xyzVector<core::Real> atom_before_rotation = ligand_res_before_rotation.atom(best_intersecting_atom).xyz();
		core::Real const xs = atom_before_rotation.x();
		core::Real const ys = atom_before_rotation.y();
		core::Real const zs = atom_before_rotation.z();

		//   r = the radius of the atom that is intersecting the ray
		core::Real const r = atom_radius.at(best_intersecting_atom);

		new_deriv_set.dDist_dv1 = dD_dv1(x0,x2,x3,y0,y2,y3,z0,z2,z3,r);
		new_deriv_set.dDist_dv2 = dD_dv2(x0,x2,x3,y0,y2,y3,z0,z2,z3,r);
		new_deriv_set.dDist_dv3 = dD_dv3(x0,x2,x3,y0,y2,y3,z0,z2,z3,r);
		new_deriv_set.dDist_dv4 = dD_dv4(x0,x2,x3,y0,y2,y3,z0,z2,z3,r,angle1_offset,angle2_offset,angle3_offset,xc,xs,yc,ys,zc,zs);
		new_deriv_set.dDist_dv5 = dD_dv5(x0,x2,x3,y0,y2,y3,z0,z2,z3,r,angle1_offset,angle2_offset,angle3_offset,xc,xs,yc,ys,zc,zs);
		new_deriv_set.dDist_dv6 = dD_dv6(x0,x2,x3,y0,y2,y3,z0,z2,z3,r,angle1_offset,angle2_offset,angle3_offset,xc,xs,yc,ys,zc,zs);

		} else {

		// set derivatives to zero if there's no intersection with the ligand or pocket
		new_deriv_set.dDist_dv1 = 0.;
		new_deriv_set.dDist_dv2 = 0.;
		new_deriv_set.dDist_dv3 = 0.;
		new_deriv_set.dDist_dv4 = 0.;
		new_deriv_set.dDist_dv5 = 0.;
		new_deriv_set.dDist_dv6 = 0.;

		}

		derivs_of_ray_distances_.push_back(new_deriv_set);

		}
		*/

	}

	// cout << "ORIG_CPU TOTAL DISTANCE: " << orig_cpu_total_dist << std::endl;
	// cout << "ORIG_CPU TOTAL NUM: " << orig_cpu_total_num << std::endl;
	// cout << "ORIG_CPU NUM INTERSECTION EVALUATIONS: " << orig_cpu_num_evaluations << std::endl;

}


core::Real PlaidFingerprint::fp_compare( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) const {

	core::Real Total_score = 0;
	core::Size num_rays = 0;
	core::Real underpack_dist = 0, steric_dist = 0, rays_missing_ligand = 0, rays_missing_pocket = 0;

	using namespace basic::options;
	bool square  = option[ OptionKeys::fingerprint::square_score ]();

	for ( auto pi = fp.triplet_fingerprint_data().begin(), li = triplet_fingerprint_data_.begin(); (pi != fp.triplet_fingerprint_data().end()) &&  (li != triplet_fingerprint_data_.end()); ++pi, ++li ) {

		// jk note: these are no longer necessarily true, since we alter the Plaid one by shifting up/down by 2*pi
		// these are useful asserts though, it's worth thinking about how to make them valid again...
		// debug_assert( std::abs( pi->phi - li->phi ) < 0.001 );
		// debug_assert( std::abs( pi->psi - li->psi ) < 0.001 );
		debug_assert( pi->ori == li->ori );

		if ( (li->rho < 0.001) && (pi->rho < 0.001) ) {
			continue;
		} else if ( li->rho < 0.001 ) {
			rays_missing_ligand++;
			Total_score += missing_point_weight;
		} else if ( pi->rho < 0.001 ) {
			rays_missing_pocket++;
			Total_score += extra_point_weight;
		} else {
			core::Real dist_deviation = std::abs( pi->rho - li->rho );
			if ( square ) {
				if ( li->rho > pi->rho ) dist_deviation *= dist_deviation;
				if ( li->rho < pi->rho ) dist_deviation *= (dist_deviation * steric_weight);
			} else if ( !square ) {
				if ( li->rho > pi->rho ) {
					/*underpack weight = 1*/
					underpack_dist += dist_deviation;
				}
				if ( li->rho < pi->rho ) {
					steric_dist += dist_deviation;
					dist_deviation = (dist_deviation * steric_weight);
				}
			}
			Total_score += dist_deviation;
		}
		num_rays++;
	}

	using namespace basic::options;
	bool get_darc_components = option[ OptionKeys::fingerprint::darc_components ]();
	if ( get_darc_components ) {
		std::cout<<"OPT "<<underpack_dist<<" "<<steric_dist<< " "<<rays_missing_ligand<<" "<<rays_missing_pocket<<" "<<num_rays;
	}

	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::return_zero_darc_score ]() ) {
		return 0.;
	}

	return (Total_score/num_rays);

}


void PlaidFingerprint::fp_compare_deriv( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, core::Real & dE_dx, core::Real & dE_dy, core::Real & dE_dz, core::Real & dE_dv4, core::Real & dE_dv5, core::Real & dE_dv6 ) const {

	dE_dx = 0.; dE_dy = 0.; dE_dz = 0.; dE_dv4 = 0.; dE_dv5 = 0.; dE_dv6 = 0.;

	if ( derivs_of_ray_distances_.size() < 2 ) {
		std::cout<<"Error, fingerprint derivatives have not been computed" << std::endl;
		exit(1);
	}
	debug_assert( derivs_of_ray_distances_.size() == fp.triplet_fingerprint_data().size() );

	core::Real Total_score = 0;
	core::Real Differentiable_score = 0;
	core::Size num_rays = 0;
	auto di = derivs_of_ray_distances_.begin();
	for ( auto pi = fp.triplet_fingerprint_data().begin(), li = triplet_fingerprint_data_.begin(); (pi != fp.triplet_fingerprint_data().end()) && (li != triplet_fingerprint_data_.end()) && (di != derivs_of_ray_distances_.end()); ++pi, ++li, ++di ) {
		debug_assert( std::abs( pi->phi - li->phi ) < 0.001 );
		debug_assert( std::abs( pi->psi - li->psi ) < 0.001 );

		if ( (li->rho < 0.001) && (pi->rho < 0.001) ) {
			continue;
		} else if ( li->rho < 0.001 ) {
			Total_score += missing_point_weight;
		} else if ( pi->rho < 0.001 ) {
			Total_score += extra_point_weight;
		} else {
			core::Real dist_deviation = std::abs( pi->rho - li->rho );
			// derivative is zero except in the case where the ray hits BOTH ligand and pocket
			// (ie. "missing point" and "extra point" scores don't contribute to the derivatives)
			if ( li->rho > pi->rho ) {
				dE_dx += di->dDist_dv1;
				dE_dy += di->dDist_dv2;
				dE_dz += di->dDist_dv3;
				dE_dv4 += di->dDist_dv4;
				dE_dv5 += di->dDist_dv5;
				dE_dv6 += di->dDist_dv6;
			} else {
				dist_deviation *= steric_weight;
				dE_dx += steric_weight * di->dDist_dv1;
				dE_dy += steric_weight * di->dDist_dv2;
				dE_dz += steric_weight * di->dDist_dv3;
				dE_dv4 += steric_weight * di->dDist_dv4;
				dE_dv5 += steric_weight * di->dDist_dv5;
				dE_dv6 += steric_weight * di->dDist_dv6;
			}
			Total_score += dist_deviation;
			Differentiable_score += dist_deviation;
		}

		num_rays++;
	}

	dE_dx /= num_rays;
	dE_dy /= num_rays;
	dE_dz /= num_rays;
	dE_dv4 /= num_rays;
	dE_dv5 /= num_rays;
	dE_dv6 /= num_rays;
	Total_score /= num_rays;
	Differentiable_score /= num_rays;

	std::cout<<"DARC score while computing derivatives: " << Total_score << std::endl;
	// std::cout<<"DARC score while computing derivatives are total: " << Total_score << " and differentiable part " << Differentiable_score << std::endl;
	// std::cout<<"Derivatives are " << dE_dx << " , " << dE_dy << " , " << dE_dz << " , " << dE_dv4 << " , " << dE_dv5 << " , " << dE_dv6 << std::endl;

	return;

}

core::Real PlaidFingerprint::search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) {

	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	return search_random_poses( fp, num_pose_search, optimal_angle1, optimal_angle2, optimal_angle3, missing_point_weight, steric_weight, extra_point_weight, no_CoM_offset);
}


core::Real PlaidFingerprint::search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & CoM_offset ) {

	core::Real best_score = std::numeric_limits<core::Real>::max();

	for ( core::Size j = 0; j < num_pose_search; ++j ) {
		core::Real curr_angle1 = (int) (numeric::random::uniform() *359.999);
		core::Real curr_angle2 = (int) (numeric::random::uniform() *359.999);
		core::Real curr_angle3 = (int) (numeric::random::uniform() *359.999);

		std::cout<< "JK this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to move_ligand_and_update_rhos_ below..." << std::endl;
		exit(1);
		move_ligand_and_update_rhos_( fp, CoM_offset, curr_angle1, curr_angle2, curr_angle3, 0 );
		core::Real curr_score = fp_compare( fp, missing_point_weight, steric_weight, extra_point_weight );
		//std::cout<<"curr_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
		if ( curr_score < best_score ) {
			best_score = curr_score;
			optimal_angle1 = curr_angle1;
			optimal_angle2 = curr_angle2;
			optimal_angle3 = curr_angle3;
			//std::cout<<"best_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
		}
	}
	return best_score;
}

core::Real PlaidFingerprint::find_optimal_rotation( FingerprintBase & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) {
	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	return find_optimal_rotation( fp, angle_increment, optimal_angle1, optimal_angle2, optimal_angle3, missing_point_weight, steric_weight, extra_point_weight, no_CoM_offset);
}

core::Real PlaidFingerprint::find_optimal_rotation( FingerprintBase & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & CoM_offset ) {

	core::Real best_score = std::numeric_limits<core::Real>::max();
	auto num_steps = core::Size ( 360. / angle_increment );

	core::Real curr_angle1=0.;
	for ( core::Size i = 0; i < num_steps; ++i ) {
		core::Real curr_angle2=0.;
		for ( core::Size j = 0; j < num_steps; ++j ) {
			core::Real curr_angle3=0.;

			for ( core::Size k = 0; k < num_steps; ++k ) {
				std::cout<< "JK this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to move_ligand_and_update_rhos_ below..." << std::endl;
				exit(1);
				move_ligand_and_update_rhos_( fp, CoM_offset, curr_angle1, curr_angle2, curr_angle3, 0 );
				core::Real curr_score = fp_compare( fp, missing_point_weight, steric_weight, extra_point_weight );
				//   std::cout<<"curr_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
				if ( curr_score < best_score ) {
					best_score = curr_score;
					optimal_angle1 = curr_angle1;
					optimal_angle2 = curr_angle2;
					optimal_angle3 = curr_angle3;
					//    std::cout<<"best_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
				}
				curr_angle3 += angle_increment;
			}
			curr_angle2 += angle_increment;
		}
		curr_angle1 += angle_increment;
	}

	return best_score;
}

void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) {

	utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	dump_oriented_pose_and_fp_to_pdb(pose_filename, fp_filename, fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, no_CoM_offset );

}

void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, numeric::xyzVector<core::Real> const & CoM_offset ) {

	utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	dump_oriented_pose_and_fp_to_pdb(pose_filename, fp_filename, fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, CoM_offset );

}

void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform ) {

	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	dump_oriented_pose_and_fp_to_pdb(pose_filename, fp_filename, fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, no_CoM_offset );

}


void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset  ) {

	std::cout<< "JK this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to move_ligand_and_update_rhos_ below..." << std::endl;
	exit(1);
	move_ligand_and_update_rhos_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset, 0 );

	core::pose::Pose tmp_pose = pose_;
	apply_rotation_offset_to_pose_( tmp_pose, angle1_offset, angle2_offset, angle3_offset );

	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians( original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians( original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians( original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> best_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	numeric::xyzMatrix<core::Real> inverse_best_mat = numeric::inverse(best_mat);
	core::Vector v(0,0,0);
	tmp_pose.apply_transform_Rx_plus_v(inverse_best_mat, v);

	numeric::xyzVector<core::Real> back_to_FingerprintBase_origin = fp.origin() - origin_;

	print_to_pdb( fp_filename, back_to_FingerprintBase_origin );

	utility::io::ozstream out_stream;
	out_stream.open(pose_filename, std::ios::out);

	out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().z()<<std::endl;
	out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()+back_to_FingerprintBase_origin.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()+back_to_FingerprintBase_origin.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()+back_to_FingerprintBase_origin.z()<<std::endl;

	core::Size lig_res_num = compute_ligand_resnum( tmp_pose );

	core::conformation::Residue const & curr_rsd = tmp_pose.conformation().residue(lig_res_num);
	for ( Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i ) {
		//out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(1) <<std::setw(8)<<std::fixed<<std::setprecision(3)<< curr_rsd.atom(i).xyz()(2) <<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(3) <<std::endl;
		out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(1)+back_to_FingerprintBase_origin.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<< curr_rsd.atom(i).xyz()(2)+back_to_FingerprintBase_origin.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(3)+back_to_FingerprintBase_origin.z() <<std::endl;
	}

	out_stream.close();
	out_stream.clear();

}


core::pose::Pose PlaidFingerprint::get_oriented_pose( FingerprintBase &fp , core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, numeric::xyzVector<core::Real> const & CoM_offset, core::Size const conformer ) {

	utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	core::pose::Pose oriented_pose = get_oriented_pose( fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, CoM_offset, conformer );
	return oriented_pose;
}

core::pose::Pose PlaidFingerprint::get_oriented_pose( FingerprintBase & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset, core::Size const conformer  ) {

	//select_conformer_move_ligand_and_update_rhos
	//note: conformer passed in is indexed from 0, below is indexed from 1
	core::pose::PoseOP tmp_pose( new core::pose::Pose(pose_, conformer+1, conformer+1) );
	apply_rotation_offset_to_pose_( *tmp_pose, angle1_offset, angle2_offset, angle3_offset );
	CoM_ = calculate_ligand_CoM(*tmp_pose);
	origin_ = fp.origin() + CoM_ - fp.CoM();
	origin_ += CoM_offset;
	core::Size const lig_res_num = compute_ligand_resnum(*tmp_pose);
	core::conformation::ResidueCOP ligand_rsd( core::conformation::ResidueOP( new core::conformation::Residue( tmp_pose->conformation().residue(lig_res_num) ) ) );
	bool update_derivatives = false;
	update_rhos_( fp, ligand_rsd, update_derivatives );
	core::pose::Pose new_pose = *tmp_pose;

	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians( original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians( original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians( original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> best_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	numeric::xyzMatrix<core::Real> inverse_best_mat = numeric::inverse(best_mat);
	core::Vector v(0,0,0);
	new_pose.apply_transform_Rx_plus_v(inverse_best_mat, v);

	// move pose such that ligand CoM is at the origin
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	numeric::xyzVector<core::Real> back_to_FingerprintBase_origin = fp.origin() - origin_;
	new_pose.apply_transform_Rx_plus_v(I_mat, back_to_FingerprintBase_origin);

	return new_pose;

}

numeric::xyzVector<core::Real> PlaidFingerprint::calculate_ligand_CoM(core::pose::Pose const & ligand_pose) {

	core::Size lig_res_num = compute_ligand_resnum( ligand_pose );
	numeric::xyzVector<core::Real> ligand_com(0.);
	core::conformation::Residue const & curr_rsd = ligand_pose.conformation().residue(lig_res_num);
	core::Size ligand_total_atoms = compute_ligand_natoms( ligand_pose );
	for ( Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i ) {
		ligand_com.x() += curr_rsd.atom(i).xyz()(1);
		ligand_com.y() += curr_rsd.atom(i).xyz()(2);
		ligand_com.z() += curr_rsd.atom(i).xyz()(3);
	}
	ligand_com /= ligand_total_atoms;

	return ligand_com;
}

core::Real PlaidFingerprint::rmsd(core::pose::Pose const & original_pose, core::pose::Pose const & oriented_pose) {

	core::Size lig_res_num = compute_ligand_resnum( original_pose );

	core::conformation::Residue const & pose1_rsd = original_pose.conformation().residue(lig_res_num);
	core::conformation::Residue const & pose2_rsd = oriented_pose.conformation().residue(lig_res_num);

	core::Real rmsd, dist_sum(0.);
	for ( Size i = 1, i_end = pose1_rsd.nheavyatoms(); i <= i_end; ++i ) {
		core::Real x_dist =  ( (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(i).xyz()(1)) * (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(i).xyz()(1)) );
		core::Real y_dist =  ( (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(i).xyz()(2)) * (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(i).xyz()(2)) );
		core::Real z_dist =  ( (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(i).xyz()(3)) * (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(i).xyz()(3)) );
		dist_sum += x_dist + y_dist + z_dist;
	}
	rmsd = sqrt(dist_sum/pose1_rsd.nheavyatoms());
	return rmsd;

}


// note: not used
//void PlaidFingerprint::move_origin( numeric::xyzVector<core::Real> const & new_origin ){
// numeric::xyzVector<core::Real> fp_coor;
//  for (std::list<spherical_coor_triplet>::iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
//  // find cartesian coors relative to old origin
//  convert_spherical_coor_triplet_to_cartesian( *pd, fp_coor );
//  // move cartesian coors from old origin to new origin
//    fp_coor += origin_ - new_origin;
//  // convert from cartesian coors to polar
//  convert_cartesian_to_spherical_coor_triplet( fp_coor, *pd );
// }
// origin_ = new_origin;
//}


// note: not used, also deprecated because we now express all angles in radians
//void correct_phi_psi( core::Real & phi, core::Real & psi ){
// while ( phi < 0. ) {
//  phi *= -1.;
//  psi += 180.;
// }
// while ( phi > 360. ) {
//  phi -= 360.;
// }
// while ( phi > 180. ) {
//  phi = 360. - phi;
//  psi += 180.;
// }
// while ( psi < -180. ) psi += 360.;
// while ( psi > 180. ) psi -= 360.;
//}


/* DEFINITIONS for the following code:
(x0,y0,z0) = the coordinates of the "origin" point from which all rays originate, called P0
(x2,y2,z2) = the coordinates of the pocket point that intersects the ray whose derivative we are calculating, called P2
(x3,y3,z3) = the coordinates of the center of the atom that is intersecting the ray, called P3
(xc,yc,zc) = the coordinates of the center of mass of the ligand BEFORE it has been moved into the current pose, called Pc
(xs,ys,zs) = the coordinates of the center of the atom that is intersecting the ray BEFORE it has been moved into the current pose, called Ps
r          = the radius of the atom that is intersecting the ray
v1         = variable that determines how much we translate the ligand's center of mass in the x direction, x component of translation vector called V13
v2         = variable that determines how much we translate the ligand's center of mass in the y direction, y component of translation vector called V13
v3         = variable that determines how much we translate the ligand's center of mass in the z direction, z component of translation vector called V13
v4         = variable that determines the angle of rotation about the x axis, called rotation matrix A4
v5         = variable that determines the angle of rotation about the y axis, called rotation matrix A5
v6         = variable that determines the angle of rotation about the z axis, called rotation matrix A6
(a,b,c)    = (x2-x0,y2-y0,z2-z0) is the directional vector for the ray, called Pd = P2 -P0
*/
/* ASSUMPTIONS for the following code:
rotations are applied to the starting vector in the order: x first, y second, z third
the coordinate basis has been chosen such that the NEGATIVE BRANCH of the solution for the intersection of the ray with the sphere is ALWAYS the correct one
P0, P2, PC, PS and r are all constants
we calculate P3 from the following equation:
P3 = A6.A5.A4.(Ps-Pc) + Pc + V13
which implies that v1-v6 are the only variables which influence P3
ALL ANGLES ARE IN RADIANS
*/

/* the function dD_dv1 calculates the value of the partial derivative of the distance, D, with respect to the variable v1 */
double dD_dv1(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r){

	double u ; // paremeter that determines the point of intersection between the ray and the sphere
	double xI; // x coordinate for the intersection point
	double yI; // y coordinate for the intersection point
	double zI; // z coordinate for the intersection point
	double d ; // current distance between PI and P2

	double du_dx3 ; // derivative of u with respect to x3
	double Q ; // constant which multiplies the derivative du_dx3 to obtain dD_dv1

	// calculate directional vector
	const double a = x2-x0;
	const double b = y2-y0;
	const double c = z2-z0;

	// calculate paremter that determines intersection point between the line and the sphere
	u =(-(a*x0) + a*x3 - b*y0 + b*y3 - c*z0 - sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3)))/2.0 + c*z3)/(a*a + b*b + c*c);

	// calculate xI, yI and zI
	xI = x0+a*u ;
	yI = y0+b*u ;
	zI = z0+c*u ;

	// calculate current distance
	d = sqrt((xI-x2)*(xI-x2)+(yI-y2)*(yI-y2)+(zI-z2)*(zI-z2)) ;

	// claculate Q
	Q = (a*(xI-x2)+b*(yI-y2)+c*(zI-z2))/d ;
	//cout << "simple version:\n" << (xI-x2)/d << endl ;

	// calculate the derivative du_dx3
	du_dx3 = (a + (2.0*(b*b*(-x0 + x3) + a*b*(y0 - y3) + c*(c*(-x0 + x3) + a*(z0 - z3))))/sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3))))/(a*a + b*b + c*c) ;

	// cout << "complex version:\n" << Q*du_dx3 << endl ;
	return Q * du_dx3 ;
}

/* the function dD_dv2 calculates the value of the partial derivative of the distance, D, with respect to the variable v2 */
double dD_dv2(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r){

	double u ; // paremeter that determines the point of intersection between the ray and the sphere
	double xI; // x coordinate for the intersection point
	double yI; // y coordinate for the intersection point
	double zI; // z coordinate for the intersection point
	double d ; // current distance between PI and P2

	double du_dy3 ; // derivative of u with respect to y3
	double Q ; // constant which multiplies the derivative du_dy3 to obtain dD_dv2

	// calculate directional vector
	const double a = x2-x0;
	const double b = y2-y0;
	const double c = z2-z0;

	// calculate paremter that determines intersection point between the line and the sphere
	u = (-(a*x0) + a*x3 - b*y0 + b*y3 - c*z0 - sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3)))/2.0 + c*z3)/(a*a + b*b + c*c);

	// calculate xI, yI and zI
	xI = x0+a*u ;
	yI = y0+b*u ;
	zI = z0+c*u ;

	// calculate current distance
	d = sqrt((xI-x2)*(xI-x2)+(yI-y2)*(yI-y2)+(zI-z2)*(zI-z2)) ;

	// claculate Q
	Q = (a*(xI-x2)+b*(yI-y2)+c*(zI-z2))/d ;

	// calculate the derivative du_dx3
	du_dy3 = (b - (2.0*(a*b*(-x0 + x3) + a*a*(y0 - y3) + c*(c*(y0 - y3) + b*(-z0 + z3))))/sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3))))/(a*a + b*b + c*c) ;

	return Q * du_dy3 ;
}

/* the function dD_dv3 calculates the value of the partial derivative of the distance, D, with respect to the variable v3 */
double dD_dv3(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r){

	double u ; // paremeter that determines the point of intersection between the ray and the sphere
	double xI; // x coordinate for the intersection point
	double yI; // y coordinate for the intersection point
	double zI; // z coordinate for the intersection point
	double d ; // current distance between PI and P2

	double du_dz3 ; // derivative of u with respect to z3
	double Q ; // constant which multiplies the derivative du_dz3 to obtain dD_dv3

	// calculate directional vector
	const double a = x2-x0;
	const double b = y2-y0;
	const double c = z2-z0;

	// calculate paremter that determines intersection point between the line and the sphere
	u = (-(a*x0) + a*x3 - b*y0 + b*y3 - c*z0 - sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3)))/2.0 + c*z3)/(a*a + b*b + c*c);

	// calculate xI, yI and zI
	xI = x0+a*u ;
	yI = y0+b*u ;
	zI = z0+c*u ;

	// calculate current distance
	d = sqrt((xI-x2)*(xI-x2)+(yI-y2)*(yI-y2)+(zI-z2)*(zI-z2)) ;

	// claculate Q
	Q = (a*(xI-x2)+b*(yI-y2)+c*(zI-z2))/d ;

	// calculate the derivative du_dx3
	du_dz3 = (c + (2.0*(a*c*(x0 - x3) + a*a*(-z0 + z3) + b*(c*(y0 - y3) + b*(-z0 + z3))))/sqrt(4.0*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3))*(a*(x0 - x3) + b*(y0 - y3) + c*(z0 - z3)) - 4.0*(a*a + b*b + c*c)*(-r*r + (x0 - x3)*(x0-x3) + (y0 - y3)*(y0 - y3) + (z0 - z3)*(z0 - z3))))/(a*a + b*b + c*c);

	return Q * du_dz3 ;
}

/* the function dD_dv4 calculates the value of the partial derivative of the distance, D, with respect to the variable v4 */
double dD_dv4(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r,const double v4,const double v5,const double v6,const double ,const double ,const double yc,const double ys,const double zc,const double zs) {

	double dx3_dv4 ; // derivative of x3 with respect to v4
	double dy3_dv4 ; // derivative of y3 with respect to v4
	double dz3_dv4 ; // derivative of y3 with respect to v4

	// calculate derivatives for P3
	dx3_dv4 = (-zc + zs)*(-(cos(v6)*sin(v4)*sin(v5)) + cos(v4)*sin(v6)) + (-yc + ys)*(cos(v4)*cos(v6)*sin(v5) + sin(v4)*sin(v6));
	dy3_dv4 = (-yc + ys)*(-(cos(v6)*sin(v4)) + cos(v4)*sin(v5)*sin(v6)) + (-zc + zs)*(-(cos(v4)*cos(v6)) - sin(v4)*sin(v5)*sin(v6)) ;
	dz3_dv4 = (-yc + ys)*cos(v4)*cos(v5) - (-zc + zs)*cos(v5)*sin(v4) ;

	return dD_dv1(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dx3_dv4 + dD_dv2(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dy3_dv4 + dD_dv3(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dz3_dv4 ;
}

/* the function dD_dv5 calculates the value of the partial derivative of the distance, D, with respect to the variable v5 */
double dD_dv5(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r,const double v4,const double v5,const double v6,const double xc,const double xs,const double yc,const double ys,const double zc,const double zs) {

	double dx3_dv5 ; // derivative of x3 with respect to v5
	double dy3_dv5 ; // derivative of y3 with respect to v5
	double dz3_dv5 ; // derivative of y3 with respect to v5

	// calculate derivatives for P3
	dx3_dv5 = (-zc + zs)*cos(v4)*cos(v5)*cos(v6) + (-yc + ys)*cos(v5)*cos(v6)*sin(v4) - (-xc + xs)*cos(v6)*sin(v5) ;
	dy3_dv5 = (-zc + zs)*cos(v4)*cos(v5)*sin(v6) + (-yc + ys)*cos(v5)*sin(v4)*sin(v6) - (-xc + xs)*sin(v5)*sin(v6) ;
	dz3_dv5 = (xc - xs)*cos(v5) - (-zc + zs)*cos(v4)*sin(v5) - (-yc + ys)*sin(v4)*sin(v5) ;

	return dD_dv1(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dx3_dv5 + dD_dv2(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dy3_dv5 + dD_dv3(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dz3_dv5 ;
}

/* the function dD_dv6 calculates the value of the partial derivative of the distance, D, with respect to the variable v6 */
double dD_dv6(const double x0,const double x2,const double x3,const double y0,const double y2,const double y3,const double z0,const double z2,const double z3,const double r,const double v4,const double v5,const double v6,const double xc,const double xs,const double yc,const double ys,const double zc,const double zs) {

	double dx3_dv6 ; // derivative of x3 with respect to v6
	double dy3_dv6 ; // derivative of y3 with respect to v6
	double dz3_dv6 ; // derivative of y3 with respect to v6

	// calculate derivatives for P3
	dx3_dv6 = (xc - xs)*cos(v5)*sin(v6) + (-zc + zs)*(cos(v6)*sin(v4) - cos(v4)*sin(v5)*sin(v6)) + (-yc + ys)*(-(cos(v4)*cos(v6)) - sin(v4)*sin(v5)*sin(v6)) ;
	dy3_dv6 = (-xc + xs)*cos(v5)*cos(v6) + (-yc + ys)*(cos(v6)*sin(v4)*sin(v5) - cos(v4)*sin(v6)) + (-zc + zs)*(cos(v4)*cos(v6)*sin(v5) + sin(v4)*sin(v6)) ;
	dz3_dv6 = 0.0 ;

	return dD_dv1(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dx3_dv6 + dD_dv2(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dy3_dv6 + dD_dv3(x0,x2,x3,y0,y2,y3,z0,z2,z3,r)*dz3_dv6 ;
}


core::Real Find_Closest_Intersect_SQ(core::Real const & phiAngle, core::Real const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius){

	// note: phi/psi are in radians

	core::Real const large_dist(999.);
	core::Real dirX,dirY,dirZ;
	//compute (dirX,dirY,dirZ) from phi/psi with random large_dist, relative to Origin
	//Reference: http://paulbourke.net/geometry/sphereline/
	dirX = large_dist*sin(phiAngle)*cos(psiAngle);
	dirY = large_dist*sin(phiAngle)*sin(psiAngle);
	dirZ = large_dist*cos(phiAngle);

	// setup our quadratic equation
	core::Real const a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);
	core::Real const b = 2.0 * ( (dirX*(-atomX)) + (dirY*(-atomY)) + (dirZ*(-atomZ)) );
	core::Real const c = atomX*atomX + atomY*atomY + atomZ*atomZ - (atom_radius * atom_radius);

	// test for intersection
	core::Real const inside_sqrt = ( b * b ) - ( 4. * a * c );

	if ( inside_sqrt > 0. ) {
		//              std::cout << "Line intersects atom\n" << std::endl;
		core::Real const inside = sqrt(inside_sqrt);
		core::Real const mu1 = -(b-inside) / ( 2. * a);
		core::Real const x1 =  mu1 * dirX;
		core::Real const y1 =  mu1 * dirY;
		core::Real const z1 =  mu1 * dirZ;
		core::Real const dist1_sq = x1*x1 + y1*y1 + z1*z1;
		core::Real const mu2 = -(b+inside) / ( 2. * a);
		core::Real const x2 = mu2 * dirX;
		core::Real const y2 = mu2 * dirY;
		core::Real const z2 = mu2 * dirZ;
		core::Real const dist2_sq = x2*x2 + y2*y2 + z2*z2;
		if ( dist2_sq < dist1_sq ) {
			return dist2_sq;
		}
		return dist1_sq;
	}
	//             std::cout <<phiAngle<<" " <<psiAngle<< " " << "No Intersection" << std::endl;
	return 9999.;
}

std::list<numeric::xyzVector<core::Real> > NonPlaidFingerprint::combine_xyz_lists(std::list< numeric::xyzVector<core::Real> > const & xyz_list_1, std::list< numeric::xyzVector<core::Real> > const & xyz_list_2 ) {

	std::list<numeric::xyzVector<core::Real> > combined_list;
	combined_list.clear();
	for ( auto const & pd : xyz_list_1 ) {
		combined_list.push_back(pd);
	}
	for ( auto const & pd : xyz_list_2 ) {
		combined_list.push_back(pd);
	}
	return combined_list;
}


//round triplet values
std::list<spherical_coor_triplet> NonPlaidFingerprint::convert_cart_to_spherical_and_round(std::list< numeric::xyzVector<core::Real> > const & xyz_list) {
	std::list< spherical_coor_triplet > rounded_triplet_list;
	rounded_triplet_list.clear();
	spherical_coor_triplet ray_triplet;
	for ( auto const & pd : xyz_list ) {
		convert_cartesian_to_spherical_coor_triplet( pd - origin_, ray_triplet );
		ray_triplet.phi = (floor(ray_triplet.phi * 100+0.5)/100);
		ray_triplet.psi = (floor(ray_triplet.psi * 100+0.5)/100);
		ray_triplet.rho = (floor(ray_triplet.rho * 100+0.5)/100);
		rounded_triplet_list.push_back(ray_triplet);
	}

	return rounded_triplet_list;
}

//set rho values to zero
std::list<spherical_coor_triplet> NonPlaidFingerprint::set_rho_to_zero(std::list<spherical_coor_triplet> const & spherical_triplet_list) {
	std::list< spherical_coor_triplet > rho_zero_triplet_list;
	rho_zero_triplet_list.clear();
	spherical_coor_triplet ray_triplet;
	for ( auto const & aa : spherical_triplet_list ) {
		ray_triplet = aa;
		ray_triplet.rho = 0;
		rho_zero_triplet_list.push_back(ray_triplet);
	}

	return rho_zero_triplet_list;
}

//remove duplicates from ext_triplet
std::list<spherical_coor_triplet> NonPlaidFingerprint::remove_duplicate_phi_psi(std::list<spherical_coor_triplet> const & rounded_triplet ) {
	std::list< spherical_coor_triplet > temp_triplet = rounded_triplet;
	std::list< spherical_coor_triplet > unique_triplet;
	unique_triplet.clear();

	for ( auto const & aa : rounded_triplet ) {
		bool found = false;
		spherical_coor_triplet best_triplet = aa;
		for ( auto bb = temp_triplet.begin(); bb != temp_triplet.end(); ) {
			if ( (aa.phi == bb->phi) && (aa.psi == bb->psi) ) {
				found = true;
				if ( bb->rho < best_triplet.rho ) { best_triplet = *bb;}
				bb = temp_triplet.erase(bb);
			} else {
				++bb;
			}
		}
		if ( found ) unique_triplet.push_back(best_triplet);
	}
	return unique_triplet;
}

std::list<numeric::xyzVector<core::Real> > NonPlaidFingerprint::convert_spherical_list_to_cartesian_list(std::list<spherical_coor_triplet> const & unique_triplet) {
	std::list<numeric::xyzVector<core::Real> > xyz_list;
	xyz_list.clear();
	for ( auto const & pd : unique_triplet ) {
		numeric::xyzVector<core::Real> new_coor;
		convert_spherical_coor_triplet_to_cartesian( pd, new_coor );
		new_coor += origin_;
		xyz_list.push_back(new_coor);
	}
	return xyz_list;
}

core::Real NonPlaidFingerprint::get_electrostatics_energy(core::pose::Pose const & ligand_pose){

	using namespace basic::options;
	core::Real const esp_wt(option[ OptionKeys::fingerprint::esp_weight ]);

	core::Real E_energy(0.);
	numeric::xyzVector<core::Real> ligand_atom(0.);

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = ligand_pose.size(); j <= resnum; ++j ) {
		if ( !ligand_pose.residue(j).is_protein() ) {
			lig_res_num = j;
			break;
		}
	}
	if ( lig_res_num == 0 ) {
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}
	core::conformation::Residue const & ligand_rsd = ligand_pose.residue(lig_res_num);
	//we include ligand hyfrogen atoms for electrostatics calculations
	Size ligand_total_atoms = ligand_rsd.natoms();

	for ( Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i ) {

		ligand_atom.x() = ligand_rsd.atom(i).xyz()(1);
		ligand_atom.y() = ligand_rsd.atom(i).xyz()(2);
		ligand_atom.z() = ligand_rsd.atom(i).xyz()(3);

		//core::Real lig_atom_esp_energy = get_nearest_neighbour_esp_energy(ligand_atom, ligand_rsd.atomic_charge(i));

		core::Real lig_atom_esp_energy = get_interpolated_esp_energy_with_type(ligand_atom, ligand_rsd.atomic_charge(i));
		//std::cout<<"ligand_atom : "<<ligand_rsd.atom_name(i)<<" " <<ligand_rsd.atomic_charge(i)<<" "<<lig_atom_esp_energy<<std::endl;

		E_energy += lig_atom_esp_energy;

	}

	using namespace basic::options;
	bool get_darc_components = option[ OptionKeys::fingerprint::darc_components ]();
	if ( get_darc_components ) {
		std::cout<<" "<<E_energy<<std::endl;

		//std::string tag_inp = option[ OptionKeys::fingerprint::inp_lig ]();
		//std::string tag_ref = option[ OptionKeys::fingerprint::ref_lig ]();

		//  ligand_pose.dump_pdb(tag_inp);
		//std::system(("/Users/ragul/Desktop/oe_ragul/bin/ragul_get_rmsd -input "+tag_inp+" -reference "+tag_ref+" > rmsd.txt").c_str());

	}

	//ESP weight
	E_energy *= esp_wt;


	using namespace basic::options;
	if ( option[ OptionKeys::fingerprint::return_zero_darc_score ]() ) {
		return 0.;
	}

	return E_energy;

}

core::Real NonPlaidFingerprint::get_interpolated_esp_energy(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge){

	numeric::xyzVector<core::Real> grid_coord;
	convert_cartesian_to_grid(ligand_atom, esp_mid_, esp_dim_, esp_spacing_, grid_coord);

	core::Real X = grid_coord.x();
	core::Real Y = grid_coord.y();
	core::Real Z = grid_coord.z();

	core::Real X1 = std::floor(X+1);
	core::Real Y1 = std::floor(Y+1);
	core::Real Z1 = std::floor(Z+1);

	core::Real X0 = std::ceil(X-1);
	core::Real Y0 = std::ceil(Y-1);
	core::Real Z0 = std::ceil(Z-1);
	core::Real atm_esp_energy;

	//if the ligand atom goes out of the grid
	if ( ((X0 >= esp_dim_.x())||(X0 < 0)) || ((Y0 >= esp_dim_.y())||(Y0 < 0)) || ((Z0 >= esp_dim_.z())||(Z0 < 0)) || ((X1 >= esp_dim_.x())||(X1 < 0)) || ((Y1 >= esp_dim_.y())||(Y1 < 0)) || ((Z1 >= esp_dim_.z())||(Z1 < 0)) ) {
		std::cout<<"\nAtom out of the grid "<<esp_dim_.x()<<" "<<esp_dim_.y()<<" "<<esp_dim_.z()<<" "<<esp_mid_.x()<<" "<<esp_mid_.y()<<" "<<esp_mid_.z()<<" "<<esp_spacing_<<std::endl;
		return 0.;
	} else {
		core::Real Xd = (X-X0)/(X1-X0);
		core::Real Yd = (Y-Y0)/(Y1-Y0);
		core::Real Zd = (Z-Z0)/(Z1-Z0);

		core::Real C00 = espGrid_[Size(X0)][Size(Y0)][Size(Z0)]*(1-Xd) + espGrid_[Size(X1)][Size(Y0)][Size(Z0)]*Xd;
		core::Real C10 = espGrid_[Size(X0)][Size(Y1)][Size(Z0)]*(1-Xd) + espGrid_[Size(X1)][Size(Y1)][Size(Z0)]*Xd;
		core::Real C01 = espGrid_[Size(X0)][Size(Y0)][Size(Z1)]*(1-Xd) + espGrid_[Size(X1)][Size(Y0)][Size(Z1)]*Xd;
		core::Real C11 = espGrid_[Size(X0)][Size(Y1)][Size(Z1)]*(1-Xd) + espGrid_[Size(X1)][Size(Y1)][Size(Z1)]*Xd;
		core::Real C0 = C00*(1-Yd) + C10*Yd;
		core::Real C1 = C01*(1-Yd) + C11*Yd;
		core::Real C = C0*(1-Zd) + C1*Zd;

		atm_esp_energy = C * atom_charge;

		return atm_esp_energy;
	}
}

core::Real NonPlaidFingerprint::get_interpolated_esp_energy_with_type(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge){

	using namespace basic::options;

	numeric::xyzVector<core::Real> grid_coord;
	convert_cartesian_to_grid(ligand_atom, esp_mid_, esp_dim_, esp_spacing_, grid_coord);

	core::Real X = grid_coord.x();
	core::Real Y = grid_coord.y();
	core::Real Z = grid_coord.z();

	core::Real X1 = std::floor(X+1);
	core::Real Y1 = std::floor(Y+1);
	core::Real Z1 = std::floor(Z+1);

	core::Real X0 = std::ceil(X-1);
	core::Real Y0 = std::ceil(Y-1);
	core::Real Z0 = std::ceil(Z-1);
	core::Real atm_esp_energy;

	//std::cout<<"TEST1 : "<<X<<" "<<Y<<" "<<Z<<" "<<X0<<" "<<Y0<<" "<<Z0<<" "<<X1<<" "<<Y1<<" "<<Z1<<std::endl;

	//if the ligand atom goes out of the grid
	if ( ((X0 >= esp_dim_.x())||(X0 < 0)) || ((Y0 >= esp_dim_.y())||(Y0 < 0)) || ((Z0 >= esp_dim_.z())||(Z0 < 0)) || ((X1 >= esp_dim_.x())||(X1 < 0)) || ((Y1 >= esp_dim_.y())||(Y1 < 0)) || ((Z1 >= esp_dim_.z())||(Z1 < 0)) ) {
		//std::cout<<"\nAtom out of the grid "<<esp_dim_.x()<<" "<<esp_dim_.y()<<" "<<esp_dim_.z()<<" "<<esp_mid_.x()<<" "<<esp_mid_.y()<<" "<<esp_mid_.z()<<" "<<esp_spacing_<<std::endl;
		return 100.;
	} else {
		core::Real Xd = (X-X0)/(X1-X0);
		core::Real Yd = (Y-Y0)/(Y1-Y0);
		core::Real Zd = (Z-Z0)/(Z1-Z0);

		core::Real C00 = espGrid_[Size(X0)][Size(Y0)][Size(Z0)]*(1-Xd) + espGrid_[Size(X1)][Size(Y0)][Size(Z0)]*Xd;
		core::Real C10 = espGrid_[Size(X0)][Size(Y1)][Size(Z0)]*(1-Xd) + espGrid_[Size(X1)][Size(Y1)][Size(Z0)]*Xd;
		core::Real C01 = espGrid_[Size(X0)][Size(Y0)][Size(Z1)]*(1-Xd) + espGrid_[Size(X1)][Size(Y0)][Size(Z1)]*Xd;
		core::Real C11 = espGrid_[Size(X0)][Size(Y1)][Size(Z1)]*(1-Xd) + espGrid_[Size(X1)][Size(Y1)][Size(Z1)]*Xd;
		core::Real C0 = C00*(1-Yd) + C10*Yd;
		core::Real C1 = C01*(1-Yd) + C11*Yd;
		core::Real C = C0*(1-Zd) + C1*Zd;
		//std::cout<<"Found_s : "<<C<<" "<<atom_charge<<" "<<C * atom_charge<<std::endl;
		atm_esp_energy = C * atom_charge;

		if ( ( typGrid_[Size(X0)][Size(Y0)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X0)][Size(Y1)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X0)][Size(Y0)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X0)][Size(Y1)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X1)][Size(Y0)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X1)][Size(Y1)][Size(Z0)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X1)][Size(Y0)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ||
				( typGrid_[Size(X1)][Size(Y1)][Size(Z1)] == ElectrostaticpotentialGrid::PROTEIN ) ) {
			// std::cout<<"\nFound_protein"<<std::endl;
			// return zero if the energy is favourable and the atom goes into protein
			if ( atm_esp_energy < 0. ) atm_esp_energy = 0.;
		}

		return atm_esp_energy;

	}
}


core::Real NonPlaidFingerprint::get_nearest_neighbour_esp_energy(numeric::xyzVector<core::Real> const & ligand_atom, core::Real const & atom_charge ){

	//convert to grid points
	core::Real grdX =(ligand_atom.x() - ( esp_mid_.x() - ((static_cast<Real>(esp_dim_.x()+1)/2) * esp_spacing_) ))/esp_spacing_;
	core::Real grdY =(ligand_atom.y() - ( esp_mid_.y() - ((static_cast<Real>(esp_dim_.y()+1)/2) * esp_spacing_) ))/esp_spacing_;
	core::Real grdZ =(ligand_atom.z() - ( esp_mid_.z() - ((static_cast<Real>(esp_dim_.z()+1)/2) * esp_spacing_) ))/esp_spacing_;

	//rounded grid points are the nearest neighbour points
	auto nn_grdX = (core::Size) std::floor(grdX+.5);
	auto nn_grdY = (core::Size) std::floor(grdY+.5);
	auto nn_grdZ = (core::Size) std::floor(grdZ+.5);

	core::Real nn_espEnergy;

	// Return 0. if the ligand atom goes out of the grid
	if ( ((nn_grdX >= esp_dim_.x())||(nn_grdX <= 0)) || ((nn_grdY >= esp_dim_.y())||(nn_grdY <= 0)) || ((nn_grdZ >= esp_dim_.z())||(nn_grdZ <= 0)) ) {
		nn_espEnergy = 0.;
	} else {
		nn_espEnergy = espGrid_[nn_grdX][nn_grdY][nn_grdZ];
	}

	//std::cout<<"ligand_atom : "<<ligand_atom.x()<<" "<<ligand_atom.y()<<" "<<ligand_atom.z()<<std::endl;
	//std::cout<<"grid_points : "<<grdX<<" "<<grdY<<" "<<grdZ<<std::endl;
	//std::cout<<"rounded grid_points : "<<nn_grdX<<" "<<nn_grdY<<" "<<nn_grdZ<<std::endl;
	//std::cout<<"nn_espEnergy : "<<nn_espEnergy<<std::endl;

	return nn_espEnergy * atom_charge;

}

core::Real NonPlaidFingerprint::get_surface_esp( std::list< numeric::xyzVector<core::Real> > const & surfacePoints_list, std::string const & inp_espGrid_fname ){


	ElectrostaticpotentialGrid eGrid;
	eGrid.get_ZAP_espGrid_values(inp_espGrid_fname);
	esp_spacing_ = eGrid.espGrid_spacing_;
	esp_mid_ = eGrid.espGrid_mid_;
	esp_dim_ = eGrid.espGrid_dim_;
	espGrid_ = eGrid.espGrid_;

	core::Real charge(1.);
	core::Real surface_esp(0.);

	for ( auto const & sf : surfacePoints_list ) {
		surface_esp += get_interpolated_esp_energy(sf, charge);
	}

	return surface_esp;

}


} // Pockets
} // protocols
