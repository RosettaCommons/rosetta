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

#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

// Protocol Headers
#include <numeric/constants.hh>
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>
// AUTO-REMOVED #include <core/init.hh>

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

#include <utility/vector1.hh>

//Auto Headers
#include <numeric/random/random.fwd.hh>


// Utility Headers

using namespace core;
using namespace core::scoring;
using namespace std;


namespace protocols {
namespace pockets {

/// @details Auto-generated virtual destructor
FingerprintBase::~FingerprintBase() {}

FingerprintBase::FingerprintBase () :
	ReferenceCount()
{
	origin_.zero();
	CoM_.zero();
}



void FingerprintBase::print_to_file(std::string const & output_filename) const {

	std::filebuf f1;
	f1.open (output_filename.c_str(), std::ios::out);
	std::ostream os1(&f1);
	std::string f1_info;
	std::stringstream f1_tmp;

	f1_tmp<<"//"<<std::fixed<<std::setprecision(2)<< origin_.x() << "\t" <<std::fixed<<std::setprecision(2)<< origin_.y() << "\t"<<std::fixed<<std::setprecision(3)<< origin_.z() <<std::endl;
	f1_info += f1_tmp.str();
	f1_tmp.str(std::string());

	for (std::list<spherical_coor_triplet>::const_iterator fi = triplet_fingerprint_data_.begin(); fi != triplet_fingerprint_data_.end(); ++fi) {
		f1_tmp<<std::fixed<<std::setprecision(2)<< fi->phi << "\t" <<std::fixed<<std::setprecision(2)<<fi->psi << "\t"<<std::fixed<<std::setprecision(3)<< fi->rho <<std::endl;
		f1_info += f1_tmp.str();
		f1_tmp.str(std::string());
	}
	os1<<f1_info;
	f1.close();

	return;
}


void FingerprintBase::print_to_pdb(std::string const & output_pdbname) const {
	numeric::xyzVector<core::Real> no_translation(0.);
	print_to_pdb( output_pdbname, no_translation );
}


void FingerprintBase::print_to_pdb(std::string const & output_pdbname, numeric::xyzVector<core::Real> const & translation) const {

	std::filebuf f2;
	f2.open (output_pdbname.c_str(), std::ios::out);
	std::ostream os2(&f2);
	std::string f2_info;
	std::stringstream f2_tmp;
	f2_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()+translation.z()<<std::endl;
	f2_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()+translation.z()<<std::endl;
	f2_info += f2_tmp.str();
	f2_tmp.str(std::string());
	for (std::list<spherical_coor_triplet>::const_iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
		numeric::xyzVector<core::Real> new_coor;
		convert_spherical_coor_triplet_to_cartesian( *pd, new_coor );

		new_coor += origin_;
		f2_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.z()+translation.z()<<std::endl;
		f2_info += f2_tmp.str();
		f2_tmp.str(std::string());
	}

	os2<<f2_info;
	f2.close();

	return;
}

void NonPlaidFingerprint::setup_from_PlaidFingerprint( PlaidFingerprint const & pfp ) {
	origin_ = pfp.origin();
	triplet_fingerprint_data_.resize(pfp.triplet_fingerprint_data().size());
	std::copy (pfp.triplet_fingerprint_data().begin(),pfp.triplet_fingerprint_data().end(), triplet_fingerprint_data_.begin());
	return;
}


void NonPlaidFingerprint::setup_from_PocketGrid( core::pose::Pose const & protein_pose, PocketGrid const & pocket_grid ) {
	EggshellGrid esg(pocket_grid);
	setup_from_EggshellGrid( protein_pose, esg );
	return;
}


void NonPlaidFingerprint::setup_from_EggshellGrid( core::pose::Pose const & protein_pose, EggshellGrid const & eggshell_grid ) {

	// set CoM_ to the pocket CoM
	CoM_ = eggshell_grid.pocket_CoM();

	// set origin_ to the protein CoM, then	move origin_ 30Angstrom away from protein center
	origin_.zero();
	core::Size total_atoms(0);
	for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
		conformation::Residue const & curr_rsd = protein_pose.residue(j);
		if ( curr_rsd.is_protein() ) {
			for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
				origin_.x() += curr_rsd.atom(i).xyz()(1);
				origin_.y() += curr_rsd.atom(i).xyz()(2);
				origin_.z() += curr_rsd.atom(i).xyz()(3);
				total_atoms++;
			}
		}
	}
	origin_ /= total_atoms;
	numeric::xyzVector<core::Real> temp_vec(0.);
	temp_vec = origin_ - eggshell_grid.pocket_CoM();
	core::Real temp_vec_magnitude = sqrt(temp_vec.x()*temp_vec.x() + temp_vec.y()*temp_vec.y() + temp_vec.z()*temp_vec.z());
	core::Real set_mag = 30/temp_vec_magnitude;
	temp_vec = temp_vec * set_mag;
	origin_ = temp_vec + eggshell_grid.pocket_CoM();

	// convert from cartesian eggshell to spherical coors
	triplet_fingerprint_data_.clear();
	spherical_coor_triplet new_triplet;
	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_grid.eggshell_coord_list().begin(); pd != eggshell_grid.eggshell_coord_list().end(); ++pd) {
		convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
		triplet_fingerprint_data_.push_back(new_triplet);
	}

	for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_grid.extra_coord_list().begin(); pd != eggshell_grid.extra_coord_list().end(); ++pd) {
		convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
		new_triplet.rho = 0.;
		triplet_fingerprint_data_.push_back(new_triplet);
	}

	return;
}

void NonPlaidFingerprint::trim_based_on_known_ligand(core::pose::Pose const & known_ligand_pose){

	protocols::pockets::PlaidFingerprint known_pf( known_ligand_pose, *this );
	std::list< spherical_coor_triplet > triplet_trim_data;
	for (std::list<spherical_coor_triplet>::const_iterator pro = triplet_fingerprint_data_.begin(),
				lig = known_pf.triplet_fingerprint_data().begin();
			pro != triplet_fingerprint_data_.end() && lig != known_pf.triplet_fingerprint_data().end();
			++pro, ++lig) {
		assert( std::abs( pro->phi - lig->phi ) < 0.001 );
		assert( std::abs( pro->psi - lig->psi ) < 0.001 );
		if (lig->rho < 0.001) continue;
		triplet_trim_data.push_back(*pro);
	}
	triplet_fingerprint_data_.clear();
	triplet_fingerprint_data_ = triplet_trim_data;
}

void NonPlaidFingerprint::setup_from_file(std::string const & input_filename) {

	ifstream inFile(input_filename.c_str());

	if (!inFile) {
		std::cout<< "Can't open input file " << input_filename << std::endl;
		exit(1);
	}

	std::string lineread;
	std::string Line;
	std::string Field;

	spherical_coor_triplet new_triplet;
	std::list<spherical_coor_triplet>::iterator it;
	triplet_fingerprint_data_.clear();

	while (std::getline(inFile, lineread)) {

		std::stringstream sss(lineread);
		std::string Pock_string_phi, Pock_string_psi, Pock_string_rho;
		core::Real Pock_real_phi, Pock_real_psi, Pock_real_rho;

		//parse COM values from line starting with "//"
		if (lineread[0] == '/' && lineread[1] == '/') {
			lineread.erase(0,2);
			std::stringstream com_line(lineread);
			std::getline(com_line, Pock_string_phi, '\t');
			origin_.x() = atof(Pock_string_phi.c_str());
			std::getline(com_line, Pock_string_psi, '\t');
			origin_.y() = atof(Pock_string_psi.c_str());
			std::getline(com_line, Pock_string_rho, '\t');
			origin_.z() = atof(Pock_string_rho.c_str());
			//std::cout<<"setupfromfile"<< " " <<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
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
		triplet_fingerprint_data_.push_back(new_triplet);

	}
	inFile.close();
	return;
}


PlaidFingerprint::PlaidFingerprint( core::pose::Pose const & input_pose, FingerprintBase const & fp ) :
		FingerprintBase(),
		pose_(input_pose)
{

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = pose_.total_residue(); j <= resnum; ++j ) {
		if (!pose_.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}

	if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}

	numeric::xyzVector<core::Real> ligand_CoM(0.);
	conformation::Residue const & curr_rsd = pose_.conformation().residue(lig_res_num);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
		ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
		ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
		ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
	}
	ligand_CoM /= curr_rsd.nheavyatoms();

	// move pose_ such that ligand CoM is at the pocket center of mass
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	pose_.apply_transform_Rx_plus_v(I_mat, fp.CoM() - ligand_CoM);

	build_from_pose_(fp);
}

void PlaidFingerprint::apply_rotation_offset_to_pose_( core::pose::Pose & pose, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) const {

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = pose.total_residue(); j <= resnum; ++j ) {
		if (!pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}

	if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}

	numeric::xyzVector<core::Real> ligand_CoM(0.);
	conformation::Residue const & curr_rsd = pose.conformation().residue(lig_res_num);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
		ligand_CoM.x() += curr_rsd.atom(i).xyz()(1);
		ligand_CoM.y() += curr_rsd.atom(i).xyz()(2);
		ligand_CoM.z() += curr_rsd.atom(i).xyz()(3);
	}
	ligand_CoM /= curr_rsd.nheavyatoms();

	// move pose such that ligand CoM is at the origin
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	pose.apply_transform_Rx_plus_v(I_mat, -ligand_CoM);

	numeric::xyzMatrix<core::Real> x_rot_mat( numeric::x_rotation_matrix_degrees(angle1_offset) );
	numeric::xyzMatrix<core::Real> y_rot_mat( numeric::y_rotation_matrix_degrees(angle2_offset) );
	numeric::xyzMatrix<core::Real> z_rot_mat( numeric::z_rotation_matrix_degrees(angle3_offset) );
	numeric::xyzMatrix<core::Real> tot_rot_mat = z_rot_mat * y_rot_mat * x_rot_mat;
	core::Vector v(0,0,0);
	pose.apply_transform_Rx_plus_v(tot_rot_mat, v);

	// move pose back to original ligand CoM
	pose.apply_transform_Rx_plus_v(I_mat, ligand_CoM);

}

void PlaidFingerprint::build_from_pose_(FingerprintBase const & fp){
	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	build_from_pose_( fp, no_CoM_offset, 0., 0., 0. );
}

void PlaidFingerprint::build_from_pose_(FingerprintBase const & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset) {

  using namespace basic::options;
	core::Real const radius_scale = option[ OptionKeys::fingerprint::atom_radius_scale ];
	core::Real const atom_buffer = option[ OptionKeys::fingerprint::atom_radius_buffer ];

	triplet_fingerprint_data_.clear();

	core::pose::Pose tmp_pose = pose_;

	apply_rotation_offset_to_pose_( tmp_pose, angle1_offset, angle2_offset, angle3_offset );

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = tmp_pose.total_residue(); j <= resnum; ++j ) {
		if (!tmp_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}

	if (lig_res_num == 0){
		std::cout<<"Error, no ligand for PlaidFingerprint" << std::endl;
		exit(1);
	}

	CoM_.zero();
	conformation::Residue const & curr_rsd = tmp_pose.conformation().residue(lig_res_num);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
		CoM_.x() += curr_rsd.atom(i).xyz()(1);
		CoM_.y() += curr_rsd.atom(i).xyz()(2);
		CoM_.z() += curr_rsd.atom(i).xyz()(3);
	}
	CoM_ /= curr_rsd.nheavyatoms();

	// set the fingerprint CoM to be the ligand CoM, then apply offset as needed
	origin_ = fp.origin() + CoM_ - fp.CoM();
	origin_ += CoM_offset;

	for (std::list<spherical_coor_triplet>::const_iterator ni = fp.triplet_fingerprint_data().begin(); ni != fp.triplet_fingerprint_data().end(); ++ni) {

		spherical_coor_triplet new_triplet;
		new_triplet.psi = ni->psi;
		new_triplet.phi = ni->phi;
		//core::Real desired_rho = ni->rho;

		// look for best_rho instead of max_rho

		core::Real atomX,atomY,atomZ,atom_radius;
		// NOTE: THIS DEPENDS ON WHETHER ORIGIN IS INSIDE OR OUTSIDE THE PROTEIN!!!
		core::Real best_rho(999.);
		//core::Real best_diff(999.);

		for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
			atomX = 0;  atomY = 0;   atomZ = 0;   atom_radius = 0;
			atomX = curr_rsd.atom(i).xyz()(1)-origin_.x();
			atomY = curr_rsd.atom(i).xyz()(2)-origin_.y();
			atomZ = curr_rsd.atom(i).xyz()(3)-origin_.z();
			//atom radius * "radius_scale" to shrink a little to match with protein surface
			atom_radius = ( curr_rsd.atom_type(i).lj_radius() - atom_buffer ) * radius_scale;
			if ( atom_radius < 0.001 ) continue;

			//			core::Real const intersect = Find_Intersect(ni->phi,ni->psi,atomX,atomY,atomZ,atom_radius,desired_rho);
			//			core::Real const intersect_diff = std::abs(intersect - desired_rho);
			//			if ( intersect_diff < best_diff ) {
			//				best_rho = intersect;
			//				best_diff = intersect_diff;
			//			}

			// for now, set best_rho to be the max rho (instead of best match)
			//			core::Real const max_intersect = Find_Intersect(ni->phi,ni->psi,atomX,atomY,atomZ,atom_radius,49.);
			//			if ( max_intersect > best_rho ) {
			//				best_rho = max_intersect;
			//			}

			// for now, set best_rho to be the max rho (instead of best match)
			// NOTE: THIS DEPENDS ON WHETHER ORIGIN IS INSIDE OR OUTSIDE THE PROTEIN!!!
			core::Real const min_intersect = Find_Intersect(ni->phi,ni->psi,atomX,atomY,atomZ,atom_radius,0.);
			if ( ( std::abs(min_intersect) > 0.00001 ) && ( min_intersect < best_rho ) ) {
				best_rho = min_intersect;
			}

		}

		if ( best_rho > 998. ) {
			best_rho = 0.;
		}

		new_triplet.rho = best_rho;
		triplet_fingerprint_data_.push_back(new_triplet);

	}

}

core::Real PlaidFingerprint::Find_Intersect(core::Real const & phiAngle, core::Real const & psiAngle, core::Real const & atomX, core::Real const & atomY, core::Real const & atomZ, core::Real const & atom_radius, core::Real const & desired_rho){

	core::Real const large_dist(999.);
	core::Real dirX,dirY,dirZ,dot_direction;

	//Define line with two points: (dirX,dirY,dirZ) and Origin (0,0,0)
	//compute (dirX,dirY,dirZ) from phi/psi with randon large_dist, relative to Origin
	//Reference: http://local.wasp.uwa.edu.au/~pbourke/geometry/sphereline/
	dirX = large_dist*sin(phiAngle*(numeric::constants::f::pi_over_180))*cos(psiAngle*(numeric::constants::f::pi_over_180));
	dirY = large_dist*sin(phiAngle*(numeric::constants::f::pi_over_180))*sin(psiAngle*(numeric::constants::f::pi_over_180));
	dirZ = large_dist*cos(phiAngle*(numeric::constants::f::pi_over_180));

	// figure out whether vector points towards atom or away from atom
	dot_direction = (dirX*atomX) + (dirY*atomY) + (dirZ*atomZ);
	if ( dot_direction < 0.0 ){
		//      std::cout <<phiAngle<<" "<<psiAngle<< " Selected direction points away from atom" << std::endl;
		return 0.;
	}

	// setup our quadratic equation
	core::Real const a = (dirX*dirX) + (dirY*dirY) + (dirZ*dirZ);
	core::Real const b = 2.0 * ( (dirX*(-atomX)) + (dirY*(-atomY)) + (dirZ*(-atomZ)) );
	core::Real const c = atomX*atomX + atomY*atomY + atomZ*atomZ - (atom_radius * atom_radius);

	// test for intersection
	core::Real const inside_sqrt = b * b - 4 * a * c;

	if ( std::abs(inside_sqrt) < 0.00001 ) {
		//      std::cout << "Line is tangent to atom\n" << std::endl;
		core::Real const mu = -b / ( 2 * a);
		core::Real const x = mu *dirX;
		core::Real const y = mu *dirY;
		core::Real const z = mu *dirZ;
		core::Real const dist = sqrt( x*x + y*y + z*z );
		return dist;
		//      std::cout <<phiAngle<< " " << psiAngle<< " "<< "Tangent at     " << x << " " << y << " " << z << std::endl;
		//             std::cout << "Distance from Origin is " << dist << std::endl;
	} else if ( inside_sqrt < 0 ) {
		//             std::cout <<phiAngle<<" " <<psiAngle<< " " << "No Intersection" << std::endl;
		return 0.;
	} else {
		//              std::cout << "Line intersects atom\n" << std::endl;
		core::Real const mu1 = -(b-sqrt(inside_sqrt)) / ( 2 * a);
		core::Real const mu2 = -(b+sqrt(inside_sqrt)) / ( 2 * a);
		core::Real const x1 =  mu1 * dirX;
		core::Real const y1 =  mu1 * dirY;
		core::Real const z1 =  mu1 * dirZ;
		core::Real const dist1 = sqrt( x1*x1 + y1*y1 + z1*z1 );
		core::Real const x2 = mu2 * dirX;
		core::Real const y2 = mu2 * dirY;
		core::Real const z2 = mu2 * dirZ;
		core::Real const dist2 = sqrt( x2*x2 + y2*y2 + z2*z2 );

		core::Real const diff1 = std::abs( dist1 - desired_rho );
		core::Real const diff2 = std::abs( dist2 - desired_rho );

		if ( diff2 < diff1 ) {
			return dist2;
		}
		return dist1;
	}
	return 0.;
}


core::Real PlaidFingerprint::fp_compare( FingerprintBase const & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) const {

	core::Real Total_score = 0;
	core::Size num_points = 0;

	for (std::list<spherical_coor_triplet>::const_iterator pi = fp.triplet_fingerprint_data().begin(),
				li = triplet_fingerprint_data_.begin();
			pi != fp.triplet_fingerprint_data().end() && li != triplet_fingerprint_data_.end();
			++pi, ++li) {
		assert( std::abs( pi->phi - li->phi ) < 0.001 );
		assert( std::abs( pi->psi - li->psi ) < 0.001 );

		if ( (li->rho < 0.001) && (pi->rho < 0.001) ) {
			continue;
		} else if (li->rho < 0.001) {
			Total_score += missing_point_weight;
		} else if (pi->rho < 0.001 ) {
			Total_score += extra_point_weight;
		} else {
			core::Real dist_deviation = std::abs( pi->rho - li->rho );
			if (li->rho > pi->rho) dist_deviation *= dist_deviation;
			if (li->rho < pi->rho) dist_deviation *= steric_weight;
			//if ((li->rho < pi->rho) && (dist_deviation > 0.5)) dist_deviation *= steric_weight;

			Total_score += dist_deviation;
		}
		num_points++;
	}
	//	if ( num_points < 25 ) return 999.;
	return (Total_score/num_points);
}


core::Real PlaidFingerprint::search_random_poses( FingerprintBase const & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) {

	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	return search_random_poses( fp, num_pose_search, optimal_angle1, optimal_angle2, optimal_angle3, missing_point_weight, steric_weight, extra_point_weight, no_CoM_offset);
}


core::Real PlaidFingerprint::search_random_poses( FingerprintBase const & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & CoM_offset ) {

	core::Real best_score = std::numeric_limits<core::Real>::max();

	for (core::Size j = 0; j < num_pose_search; ++j ){
		core::Real curr_angle1 = (int) (numeric::random::uniform() *359.999);
		core::Real curr_angle2 = (int) (numeric::random::uniform() *359.999);
		core::Real curr_angle3 = (int) (numeric::random::uniform() *359.999);

		build_from_pose_( fp, CoM_offset, curr_angle1, curr_angle2, curr_angle3 );
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

core::Real PlaidFingerprint::find_optimal_rotation( FingerprintBase const & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) {
	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	return find_optimal_rotation( fp, angle_increment, optimal_angle1, optimal_angle2, optimal_angle3, missing_point_weight, steric_weight, extra_point_weight, no_CoM_offset);
}

core::Real PlaidFingerprint::find_optimal_rotation( FingerprintBase const & fp, core::Real const & angle_increment, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & CoM_offset ) {

	core::Real best_score = std::numeric_limits<core::Real>::max();
	core::Size num_steps = core::Size ( 360. / angle_increment );

	core::Real curr_angle1=0.;
  for (core::Size i = 0; i < num_steps; ++i ){
		core::Real curr_angle2=0.;
		for (core::Size j = 0; j < num_steps; ++j ){
			core::Real curr_angle3=0.;

			for (core::Size k = 0; k < num_steps; ++k ){
				build_from_pose_( fp, CoM_offset, curr_angle1, curr_angle2, curr_angle3 );
				core::Real curr_score = fp_compare( fp, missing_point_weight, steric_weight, extra_point_weight );
				//			std::cout<<"curr_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
				if ( curr_score < best_score ) {
					best_score = curr_score;
					optimal_angle1 = curr_angle1;
					optimal_angle2 = curr_angle2;
					optimal_angle3 = curr_angle3;
					//				std::cout<<"best_score "<<curr_score<< " " << curr_phi << " " <<curr_psi << std::endl;
				}
				curr_angle3 += angle_increment;
			}
			curr_angle2 += angle_increment;
		}
		curr_angle1 += angle_increment;
	}

	return best_score;
}

void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset ) {

	utility::vector1<core::Real> original_pocket_angle_transform(3, 0.);
	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	dump_oriented_pose_and_fp_to_pdb(pose_filename, fp_filename, fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, no_CoM_offset );

}

void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform ) {

	numeric::xyzVector<core::Real> no_CoM_offset(0.);
	dump_oriented_pose_and_fp_to_pdb(pose_filename, fp_filename, fp, angle1_offset, angle2_offset, angle3_offset, original_pocket_angle_transform, no_CoM_offset );

}


void PlaidFingerprint::dump_oriented_pose_and_fp_to_pdb( std::string const & pose_filename, std::string const & fp_filename, FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset  ) {

	build_from_pose_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset);

	core::pose::Pose tmp_pose = pose_;
 	apply_rotation_offset_to_pose_( tmp_pose, angle1_offset, angle2_offset, angle3_offset );

	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_degrees( original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_degrees( original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_degrees( original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> best_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	numeric::xyzMatrix<core::Real> inverse_best_mat = numeric::inverse(best_mat);
	core::Vector v(0,0,0);
  tmp_pose.apply_transform_Rx_plus_v(inverse_best_mat, v);

	numeric::xyzVector<core::Real> back_to_FingerprintBase_origin = fp.origin() - origin_;

	print_to_pdb( fp_filename, back_to_FingerprintBase_origin );

	std::filebuf f3;
	f3.open (pose_filename.c_str(), std::ios::out);
	std::ostream os3(&f3);
	std::string f3_info;
	std::stringstream f3_tmp;
	f3_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<fp.origin().z()<<std::endl;
	f3_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()+back_to_FingerprintBase_origin.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()+back_to_FingerprintBase_origin.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()+back_to_FingerprintBase_origin.z()<<std::endl;

f3_info += f3_tmp.str();
	f3_tmp.str(std::string());
	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = tmp_pose.total_residue(); j <= resnum; ++j ) {
		if (!tmp_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}

	conformation::Residue const & curr_rsd = tmp_pose.conformation().residue(lig_res_num);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
		//f3_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(1) <<std::setw(8)<<std::fixed<<std::setprecision(3)<< curr_rsd.atom(i).xyz()(2) <<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(3) <<std::endl;
		f3_tmp<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(1)+back_to_FingerprintBase_origin.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<< curr_rsd.atom(i).xyz()(2)+back_to_FingerprintBase_origin.y() <<std::setw(8)<<std::fixed<<std::setprecision(3)<<curr_rsd.atom(i).xyz()(3)+back_to_FingerprintBase_origin.z() <<std::endl;
		f3_info += f3_tmp.str();
		f3_tmp.str(std::string());
	}
	os3<<f3_info;
	f3.close();
}

core::pose::Pose PlaidFingerprint::get_oriented_pose( FingerprintBase const & fp, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, utility::vector1<core::Real> const & original_pocket_angle_transform, numeric::xyzVector<core::Real> const & CoM_offset  ) {

	build_from_pose_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset);

	core::pose::Pose tmp_pose = pose_;
 	apply_rotation_offset_to_pose_( tmp_pose, angle1_offset, angle2_offset, angle3_offset );

	numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_degrees( original_pocket_angle_transform[1] ) );
	numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_degrees( original_pocket_angle_transform[2] ) );
	numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_degrees( original_pocket_angle_transform[3] ) );
	numeric::xyzMatrix<core::Real> best_mat = bestz_rot_mat * besty_rot_mat * bestx_rot_mat;
	numeric::xyzMatrix<core::Real> inverse_best_mat = numeric::inverse(best_mat);
	core::Vector v(0,0,0);
  tmp_pose.apply_transform_Rx_plus_v(inverse_best_mat, v);


	// move pose such that ligand CoM is at the origin
	numeric::xyzMatrix<core::Real> I_mat;
	I_mat.to_identity();
	numeric::xyzVector<core::Real> back_to_FingerprintBase_origin = fp.origin() - origin_;
	tmp_pose.apply_transform_Rx_plus_v(I_mat, back_to_FingerprintBase_origin);

	return tmp_pose;

}


core::Real PlaidFingerprint::rmsd(core::pose::Pose const & original_pose, core::pose::Pose const & oriented_pose) {

	core::Size lig_res_num = 0;
	for ( int j = 1, resnum = original_pose.total_residue(); j <= resnum; ++j ) {
		if (!original_pose.residue(j).is_protein()){
			lig_res_num = j;
			break;
		}
	}

	if (lig_res_num == 0){
		std::cout<<"Error, no ligand for RMSD Calculation" << std::endl;
		exit(1);
	}

	conformation::Residue const & pose1_rsd = original_pose.conformation().residue(lig_res_num);
	conformation::Residue const & pose2_rsd = oriented_pose.conformation().residue(lig_res_num);

	core::Real rmsd, dist_sum(0.);
	for(Size i = 1, i_end = pose1_rsd.nheavyatoms(); i <= i_end; ++i) {
		core::Real x_dist =  ( (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(i).xyz()(1)) * (pose1_rsd.atom(i).xyz()(1) - pose2_rsd.atom(i).xyz()(1)) );
		core::Real y_dist =  ( (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(i).xyz()(2)) * (pose1_rsd.atom(i).xyz()(2) - pose2_rsd.atom(i).xyz()(2)) );
		core::Real z_dist =  ( (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(i).xyz()(3)) * (pose1_rsd.atom(i).xyz()(3) - pose2_rsd.atom(i).xyz()(3)) );
		dist_sum += x_dist + y_dist + z_dist;
	}
	rmsd = sqrt(dist_sum/pose1_rsd.nheavyatoms());
	return rmsd;

}


void PlaidFingerprint::move_origin( numeric::xyzVector<core::Real> const & new_origin ){

	numeric::xyzVector<core::Real> fp_coor;
  for (std::list<spherical_coor_triplet>::iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
		// find cartesian coors relative to old origin
		convert_spherical_coor_triplet_to_cartesian( *pd, fp_coor );
		// move cartesian coors from old origin to new origin
    fp_coor += origin_ - new_origin;
		// convert from cartesian coors to polar
		convert_cartesian_to_spherical_coor_triplet( fp_coor, *pd );
	}
	origin_ = new_origin;
}


void correct_phi_psi( core::Real & phi, core::Real & psi ){

	while ( phi < 0. ) {
		phi *= -1.;
		psi += 180.;
	}
	while ( phi > 360. ) {
		phi -= 360.;
	}
	while ( phi > 180. ) {
		phi = 360. - phi;
		psi += 180.;
	}
	while ( psi < -180. ) psi += 360.;
	while ( psi > 180. ) psi -= 360.;
}

} // Pockets
} // protocols
