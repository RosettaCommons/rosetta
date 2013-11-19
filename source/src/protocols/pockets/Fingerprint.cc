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
#include <protocols/pockets/Fingerprint.hh>
#include <protocols/pockets/PocketGrid.hh>

// Core Headers
#include <core/pose/Pose.hh>
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
#include <basic/options/keys/fingerprint.OptionKeys.gen.hh>
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
FingerprintBase::~FingerprintBase() {}

void FingerprintBase::print_to_file(std::string const & output_filename) const {

  utility::io::ozstream out_stream;
  out_stream.open(output_filename, std::ios::out);
  out_stream<<"/ORI/"<<std::fixed<<std::setprecision(2)<< origin_.x() << "\t" <<std::fixed<<std::setprecision(2)<< origin_.y() << "\t"<<std::fixed<<std::setprecision(3)<< origin_.z() <<std::endl;
  out_stream<<"/COM/"<<std::fixed<<std::setprecision(2)<< CoM_.x() << "\t" <<std::fixed<<std::setprecision(2)<< CoM_.y() << "\t"<<std::fixed<<std::setprecision(3)<< CoM_.z() <<std::endl;
  for (std::list<spherical_coor_triplet>::const_iterator fi = triplet_fingerprint_data_.begin(); fi != triplet_fingerprint_data_.end(); ++fi) {
    out_stream<< fi->phi << "\t" <<fi->psi << "\t"<< fi->rho <<std::endl;
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

  out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()+translation.z()<<std::endl;
  out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM X   0    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()+translation.z()<<std::endl;
  for (std::list<spherical_coor_triplet>::const_iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
    numeric::xyzVector<core::Real> new_coor;
    convert_spherical_coor_triplet_to_cartesian( *pd, new_coor );
    new_coor += origin_;
    out_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   MAP A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.x()+translation.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.y()+translation.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.z()+translation.z()<<std::endl;
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

  //set origin
  set_origin( protein_pose, egg_and_ext_list_);

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

  //set origin
  set_origin( protein_pose, egg_and_ext_list_);

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

  //set origin
  set_origin( protein_pose, egg_and_ext_list_);
  setup_from_EggshellGrid();

  return;
}

void NonPlaidFingerprint::write_eggshell_to_pdb_file( std::string const & output_eggshell_name ) const {

  utility::io::ozstream outPDB_stream;
  outPDB_stream.open(output_eggshell_name, std::ios::out);
  outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()<<std::endl;
  outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()<<std::endl;
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = eggshell_list_.begin(); pd != eggshell_list_.end(); ++pd) {
    outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
  }
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = extshell_list_.begin(); pd != extshell_list_.end(); ++pd) {
    outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  O   EXT B   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<pd->z()<<std::endl;
  }
  outPDB_stream.close();
  outPDB_stream.clear();

}

void NonPlaidFingerprint::set_origin( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell ) {
  using namespace basic::options;
  int origin_option = option[ OptionKeys::fingerprint::set_origin ];
  Size best_origin_option(0);
  if (origin_option == 0){
    //choose best origin by finding the lowest R (ruggedness) value from different origin position
    core::Real best_R(999.);
    for(Size set_origin_option = 1; set_origin_option < 4; ++set_origin_option){
      core::Real new_R = get_Rvalue( protein_pose, egg_and_extra_shell, set_origin_option );
      //std::cout << "set_origin_option " << set_origin_option << " new_Rvalue " << new_R << std::endl;
      if ( new_R < best_R ) {
	best_R = new_R;
	best_origin_option = set_origin_option;
      }
    }
    //std::cout << "best_origin_option " << best_origin_option << " best_Rvalue " << best_R << std::endl;
  } else {
    best_origin_option = origin_option;
  }
  set_origin_from_option_( protein_pose, egg_and_extra_shell, best_origin_option );

}

core::Real NonPlaidFingerprint::get_Rvalue( core::pose::Pose const & protein_pose, std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, Size const & set_origin_option ) {

  set_origin_from_option_( protein_pose, egg_and_extra_shell, set_origin_option );

  core::Real min_rho_difference(0.);
  core::Size num_points(0);
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pt1 = egg_and_extra_shell.begin(); pt1 != egg_and_extra_shell.end(); ++pt1) {
    numeric::xyzVector< core::Real > eggshell_point1 = *pt1 - origin_;
    spherical_coor_triplet triplet1;
    convert_cartesian_to_spherical_coor_triplet( eggshell_point1, triplet1 );
    core::Real min_angle(999.);
    spherical_coor_triplet triplet2;
    for (std::list< numeric::xyzVector<core::Real> >::const_iterator pt2 = egg_and_extra_shell.begin(); pt2 != egg_and_extra_shell.end(); ++pt2) {
      if ( pt1 == pt2 ) continue;
      core::Real const curr_angle = std::abs(cos_of( *pt1, *pt2 ));
      if(curr_angle < min_angle){
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

  if (set_origin_option == 1){
    set_origin_away_from_protein_center(protein_pose);
    //std::cout<< " ORIGIN1 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
  } else if (set_origin_option == 2){
    set_origin_away_from_eggshell(egg_and_extra_shell, protein_pose);
    //std::cout<< " ORIGIN2 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
  } else if (set_origin_option == 3){
    set_origin_away_from_eggshell_plane(egg_and_extra_shell, protein_pose, set_origin_option);
    //std::cout<< " ORIGIN3 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
  } else if (set_origin_option == 4){
    set_origin_away_from_eggshell_plane(egg_and_extra_shell, protein_pose, set_origin_option);
    //std::cout<< " ORIGIN4 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
  } else if (set_origin_option == 5){
    set_origin_from_residue(protein_pose);
    //std::cout<< " ORIGIN5 "<<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
  } else {
    std::cout<<"Error, wrong option for set_origin" << std::endl;
    exit(1);
  }
  return;
}

void NonPlaidFingerprint::setup_from_EggshellGrid() {

  // convert from cartesian eggshell to spherical coors
	std::list< spherical_coor_triplet > rounded_egg_triplet = convert_cart_to_spherical_and_round(eggshell_list_);
	std::list< spherical_coor_triplet > rounded_ext_triplet = convert_cart_to_spherical_and_round(extshell_list_);

	std::list< spherical_coor_triplet >	unq_egg_triplet = remove_duplicate_phi_psi(rounded_egg_triplet);
	std::list< spherical_coor_triplet >	unq_ext_triplet = remove_duplicate_phi_psi(rounded_ext_triplet);

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

  //DUMP EGGSHELL TO A PDB FILE
  utility::io::ozstream outPDB_stream;
  outPDB_stream.open("eggshell.pdb", std::ios::out);
  outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   ORI A   1    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<origin_.z()<<std::endl;
  outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   COM B   2    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<CoM_.z()<<std::endl;
  for (std::list<spherical_coor_triplet>::const_iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
    numeric::xyzVector<core::Real> new_coor;
    convert_spherical_coor_triplet_to_cartesian( *pd, new_coor );
    new_coor += origin_;
    outPDB_stream<<"HETATM   "<<std::setw(2)<<1<<"  C   EGG C   3    "<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.x()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.y()<<std::setw(8)<<std::fixed<<std::setprecision(3)<<new_coor.z()<<std::endl;
  }
  outPDB_stream.close();
  outPDB_stream.clear();


#ifdef USEOPENCL
  gpu_setup_rays();
#endif

	return;
}


#ifdef USEOPENCL
void NonPlaidFingerprint::gpu_setup( core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, int & num_particles, PlaidFingerprint & pf ) {

  // Initialize GPU memory and variables
  if(!gpu_.use()) return;

  if(!gpu_.RegisterProgram("gpu/DARC_PSO.cl")) utility_exit_with_message("Failed to load OpenCL program");

  // Custom weights
  float weights[4] = {
    missing_point_weight,
    steric_weight,
    extra_point_weight,
    num_particles
  };

  // Allocate GPU memory for weights (reused), and copy weights to GPU
  if(!gpu_memory_.weights) gpu_memory_.weights = gpu_.AllocateMemory(sizeof(weights));
  gpu_.WriteData(gpu_memory_.weights, weights, sizeof(weights));
  gpu_memory_.num_atoms = pf.compute_ligand_natoms();
  gpu_memory_.num_particles = num_particles;
}

void NonPlaidFingerprint::gpu_setup_rays() {

  if(!gpu_.use()) return;

  std::vector<basic::gpu::float4> rays;

  for (std::list<spherical_coor_triplet>::const_iterator fi = triplet_fingerprint_data_.begin(); fi != triplet_fingerprint_data_.end(); ++fi) {
    basic::gpu::float4 ray = { fi->phi, fi->psi, fi->rho, 0.0 };
    rays.push_back(ray);
  }

  // Allocate GPU memory (reused) and copy ray coordinates over to GPU
  unsigned int rays_size = sizeof(basic::gpu::float4) * rays.size();
  gpu_.AllocateMemoryReuse(gpu_memory_.rays, gpu_memory_.rays_size, rays_size);
  gpu_.WriteData(gpu_memory_.rays, &rays[0], rays_size);
  gpu_memory_.num_rays = rays.size();
}

int NonPlaidFingerprint::gpu_calculate_particle_scores( core::optimization::ParticleOPs & particles, std::vector<basic::gpu::float4> &atoms, std::vector<basic::gpu::float4> &atom_maxmin_phipsi ) {

  if(!gpu_.use()) return 0;

  // For timing of the rest of the function from here on:
  // basic::gpu::Timer t("gpu_calculate_particle_scores");

  // Allocate GPU memory for particles, intermediate results, and final scores (all reused)
  gpu_.AllocateMemoryReuse(gpu_memory_.atoms, gpu_memory_.atoms_size, sizeof(basic::gpu::float4) * atoms.size());
  gpu_.AllocateMemoryReuse(gpu_memory_.atom_maxmin_phipsi, gpu_memory_.atom_maxmin_phipsi_size, sizeof(basic::gpu::float4) * atom_maxmin_phipsi.size());
  gpu_.AllocateMemoryReuse(gpu_memory_.particle_scores, gpu_memory_.particle_scores_size, sizeof(float) * gpu_memory_.num_particles);
  gpu_.AllocateMemoryReuse(gpu_memory_.ray_scores, gpu_memory_.ray_scores_size, sizeof(float) * gpu_memory_.num_particles * gpu_memory_.num_rays);

  // Copy particle coordinates over to GPU
  gpu_.WriteData(gpu_memory_.atoms, &atoms[0], gpu_memory_.atoms_size);
  // gpu_.WriteData(gpu_memory_.atom_maxmin_phipsi, &atom_maxmin_phipsi[0], gpu_memory_.atom_maxmin_phipsi_size);

  // Execute kernels
  if(!gpu_.ExecuteKernel("Check_for_intersection", gpu_memory_.num_rays, gpu_memory_.num_rays, 64,
                           GPU_DEVMEM, gpu_memory_.rays,
                           GPU_DEVMEM, gpu_memory_.atoms,
                           GPU_DEVMEM, gpu_memory_.ray_scores,
                           GPU_IN | GPU_INT, gpu_memory_.num_atoms,
                           GPU_DEVMEM, gpu_memory_.weights,
                           GPU_IN | GPU_INT, gpu_memory_.num_rays,
                           NULL)) {
    std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
    return 0;
  }

  if(!gpu_.ExecuteKernel("Get_scores", gpu_memory_.num_particles, gpu_memory_.num_particles, 64,
                           GPU_DEVMEM, gpu_memory_.ray_scores,
                           GPU_DEVMEM, gpu_memory_.particle_scores,
                           GPU_IN | GPU_INT, gpu_memory_.num_rays,
                           GPU_IN | GPU_INT, gpu_memory_.num_particles,
                           NULL)) {
    std::cout << "Failed to launch kernel: " << gpu_.lastErrorStr() << std::endl;
    return 0;
  }

  std::vector<float> particle_scores(gpu_memory_.num_particles);
  gpu_.ReadData(&particle_scores[0], gpu_memory_.particle_scores, gpu_memory_.particle_scores_size);

  Size j =0;
  for (std::vector<float>::const_iterator ci = particle_scores.begin(); ci != particle_scores.end(); ++ci) {
    core::Real score = (float)*ci;
    particles[++j]->set_score(score);
  }

  return 1;
}
#endif

void NonPlaidFingerprint::setup_from_eggshell_pdb_file(std::string const & input_filename) {

  ifstream inFile(input_filename.c_str());
  if (!inFile) {
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
  std::string chainID;
  std::string restype;

  numeric::xyzVector<core::Real> eggshell_CoM;
  core::Size comcounter = 0;
  core::Size oricounter = 0;

  while (std::getline(inFile, line)) {
    name = line.substr(0,6);
    chainID = line.substr(21,1);
    restype = line.substr(17,3);

    numeric::xyzVector<core::Real> pdb_coord;
    std::string Xstring, Ystring, Zstring;

    if (name=="END") {
			break;
		} else if ((name=="HETATM")&&(restype=="ORI")) {
			Xstring = line.substr(30,8);
			origin_.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			origin_.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			origin_.z() = atof(Zstring.c_str());
			oricounter++;
		} else if ((name=="HETATM")&&(restype=="COM")) {
			Xstring = line.substr(30,8);
			CoM_.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			CoM_.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			CoM_.z() = atof(Zstring.c_str());
			comcounter++;
		} else if ((name=="HETATM")&&(restype=="EGG")) {
			Xstring = line.substr(30,8);
			pdb_coord.x() = atof(Xstring.c_str());
			Ystring = line.substr(38,8);
			pdb_coord.y() = atof(Ystring.c_str());
			Zstring = line.substr(46,8);
			pdb_coord.z() = atof(Zstring.c_str());
			temp_eggshell_coord_list.push_back(pdb_coord);
		} else if ((name=="HETATM")&&(restype=="EXT")) {
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

  if ( (oricounter == 0) || (comcounter == 0) ){
    std::cout<<"Error, No ORIGIN or CoM value specified in input Eggshell file" << std::endl;
    exit(1);
  } else if ( (oricounter > 1) || (comcounter >1) ){
    std::cout<<"Error, More than one ORIGIN or CoM value specified in input Eggshell file" << std::endl;
    exit(1);
  }

  // convert from cartesian_coord to spherical coors
  spherical_coor_triplet new_triplet;
  triplet_fingerprint_data_.clear();

  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = temp_eggshell_coord_list.begin(); pd != temp_eggshell_coord_list.end(); ++pd) {
    convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
    triplet_fingerprint_data_.push_back(new_triplet);
  }

  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = temp_extra_coord_list.begin(); pd != temp_extra_coord_list.end(); ++pd) {
    convert_cartesian_to_spherical_coor_triplet( *pd - origin_, new_triplet );
    new_triplet.rho = 0.;
    triplet_fingerprint_data_.push_back(new_triplet);
  }

  return;
}

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
  triplet_fingerprint_data_.clear();

  while (std::getline(inFile, lineread)) {

    std::stringstream sss(lineread);
    std::string Pock_string_phi, Pock_string_psi, Pock_string_rho;
    core::Real Pock_real_phi, Pock_real_psi, Pock_real_rho;

    //parse ORIGIN values from line starting with "/ORI/"
    if (lineread[0] == '/' && lineread[1] == 'O' && lineread[2] == 'R' && lineread[3] == 'I' && lineread[4] == '/') {
      lineread.erase(0,5);
      std::stringstream ori_line(lineread);
      std::getline(ori_line, Pock_string_phi, '\t');
      origin_.x() = atof(Pock_string_phi.c_str());
      std::getline(ori_line, Pock_string_psi, '\t');
      origin_.y() = atof(Pock_string_psi.c_str());
      std::getline(ori_line, Pock_string_rho, '\t');
      origin_.z() = atof(Pock_string_rho.c_str());
      //std::cout<<"setup from triplet file : "<< " " <<origin_.x()<<" "<<origin_.y()<<" "<<origin_.z()<<std::endl;
      continue;
    }
    //parse CoM values from line starting with "/COM/"
    if (lineread[0] == '/' && lineread[1] == 'C' && lineread[2] == 'O' && lineread[3] == 'M' && lineread[4] == '/') {
      lineread.erase(0,5);
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

#ifdef USEOPENCL
  gpu_setup_rays();
#endif

  return;
}

void NonPlaidFingerprint::trim_based_on_known_ligand(core::pose::Pose const & known_ligand_pose){

  protocols::pockets::PlaidFingerprint known_pf( known_ligand_pose, *this );
  std::list< spherical_coor_triplet > triplet_trim_data;
  for (std::list<spherical_coor_triplet>::const_iterator pro = triplet_fingerprint_data_.begin(), lig = known_pf.triplet_fingerprint_data().begin(); pro != triplet_fingerprint_data_.end(), lig != known_pf.triplet_fingerprint_data().end(); ++pro, ++lig) {

    //jk note: these are no longer necessarily true, since we alter the Plaid one by shifting up/down by 2*pi
    //these are useful asserts though, it's worth thinking about how to make them valid again...
    //assert( std::abs( pro->phi - lig->phi ) < 0.001 );
    //assert( std::abs( pro->psi - lig->psi ) < 0.001 );

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
  for (std::list<spherical_coor_triplet>::const_iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
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
  for ( int j = 1, resnum = known_ligand_pose.total_residue(); j <= resnum; ++j ) {
    if (!known_ligand_pose.residue(j).is_protein()){
      lig_res_num = j;
      break;
    }
  }
  if (lig_res_num == 0){
    std::cout<<"Error, no ligand to include_eggshell_points_based_on_known_ligand" << std::endl;
    exit(1);
  }

  core::conformation::Residue const & curr_rsd = known_ligand_pose.conformation().residue(lig_res_num);
  core::Size ligand_total_atoms = curr_rsd.nheavyatoms();
  numeric::xyzVector<core::Real> lig_atom_coord;
  std::list< numeric::xyzVector<core::Real> > lig_atom_coord_list;
  for(Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i) {
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
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator aa = eggshell_list_.begin(); aa != eggshell_list_.end(); ++aa) {
		xyz_coord = *aa;
		bool found = false;
		for (std::list< numeric::xyzVector<core::Real> >::const_iterator bb = lig_atom_coord_list.begin(); bb != lig_atom_coord_list.end(); ++bb) {
			if(	xyz_coord.distance(*bb) <= trim_dist ) {found = true;break;}
		}
		if (found){
			new_egg_coord_list.push_back(xyz_coord);
		}
		//		else{
		//	new_ext_coord_list.push_back(xyz_coord);
		//	}
	}
	eggshell_list_.clear();
	eggshell_list_ = new_egg_coord_list;
	//	extshell_list_.clear();
	//extshell_list_ = new_ext_coord_list;

  return;
}

void NonPlaidFingerprint::set_origin_away_from_eggshell_plane( std::list< numeric::xyzVector<core::Real> > const & egg_and_extra_shell, core::pose::Pose const & protein_pose, Size const & set_origin_option ) {

  origin_.zero();

  core::Real A11(0.), A12(0.), A13(0.), A22(0.), A23(0.), b1(0.), b2(0.), b3(0.);
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = egg_and_extra_shell.begin(); pd != egg_and_extra_shell.end(); ++pd) {
    A11 += pd->x()*pd->x();
    A12 += pd->x()*pd->y();
    A13 += pd->x();
    A22 += pd->y()*pd->y();
    A23 += pd->y();
    b1 += pd->x()*pd->z();
    b2 += pd->y()*pd->z();
    b3 += pd->z();
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

  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pdebug = egg_and_extra_shell.begin(); pdebug != egg_and_extra_shell.end(); ++pdebug) {
    plane_coord.x() = pdebug->x();
    plane_coord.y() = pdebug->y();
    plane_coord.z() = best_fit_a * pdebug->x() + best_fit_b * pdebug->y() + best_fit_c;
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
  if(	protein_CoM.distance(temp_origin1) > protein_CoM.distance(temp_origin2) ) {
    closest_origin = temp_origin2;
    distant_origin = temp_origin1;
  } else {
    closest_origin = temp_origin1;
    distant_origin = temp_origin2;
  }

  if (set_origin_option == 3){
    origin_ = closest_origin;
  } else if (set_origin_option == 4){
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
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = egg_and_extra_shell.begin(); pd != egg_and_extra_shell.end(); ++pd) {
    temp_vec += (CoM_ - *pd);
  }
  temp_vec.normalize(30.);

  numeric::xyzVector<core::Real> temp_origin1 = CoM_ + temp_vec;
  numeric::xyzVector<core::Real> temp_origin2 = CoM_ - temp_vec;
  numeric::xyzVector<core::Real> closest_origin(0.);
  numeric::xyzVector<core::Real> distant_origin(0.);

  //calculate distance b/w protein_CoM and +/- origin and choose the shortest
  numeric::xyzVector<core::Real> protein_CoM = calculate_protein_CoM(protein_pose);
  if(	protein_CoM.distance(temp_origin1) > protein_CoM.distance(temp_origin2) ) {
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

  // set origin_ to the protein CoM, then	move origin_ 30Angstrom away from protein center
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
	int const ori_resid(option[ OptionKeys::fingerprint::origin_res_num ]);
	//	protein_pose.residue(ori_resid).atoms();
	core::conformation::Residue const & curr_rsd = protein_pose.residue(ori_resid);
  numeric::xyzVector<core::Real> residue_com(0.);
	for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
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

numeric::xyzVector<core::Real> NonPlaidFingerprint::calculate_protein_CoM(core::pose::Pose const & protein_pose) {

  numeric::xyzVector<core::Real> protein_com(0.);
  core::Size total_atoms(0);
  for ( int j = 1, resnum = protein_pose.total_residue(); j <= resnum; ++j ) {
    core::conformation::Residue const & curr_rsd = protein_pose.residue(j);
    if ( curr_rsd.is_protein() ) {
      for(Size i = 1, i_end = curr_rsd.nheavyatoms(); i <= i_end; ++i) {
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

  // move pose_ such that ligand CoM is at the pocket center of mass
  numeric::xyzMatrix<core::Real> I_mat;
  I_mat.to_identity();
  pose_.apply_transform_Rx_plus_v(I_mat, fp.CoM() - ligand_CoM);

  core::Size const lig_res_num = compute_ligand_resnum();
  //	core::Size const ligand_natoms = compute_ligand_natoms();
  core::conformation::ResidueCOP ligand_rsd = new core::conformation::Residue( pose_.conformation().residue(lig_res_num) );

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
  for(Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i) {
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
  return lig_res_num;
}

core::Size PlaidFingerprint::compute_ligand_natoms( core::pose::Pose const & pose ) const {
  core::Size lig_res_num = compute_ligand_resnum(pose);
  core::conformation::Residue const & curr_rsd = pose.conformation().residue(lig_res_num);
  core::Size ligand_total_atoms;
  using namespace basic::options;
  if (option[ OptionKeys::fingerprint::include_hydrogens ]()){
    ligand_total_atoms = curr_rsd.natoms();
  } else {
    ligand_total_atoms = curr_rsd.nheavyatoms();
  }
  return ligand_total_atoms;
}

core::Size PlaidFingerprint::compute_ligand_nconformers( core::pose::Pose const & pose ) const {
  return pose.total_residue();
}

core::conformation::ResidueCOP PlaidFingerprint::select_conf_and_move_ligand_( FingerprintBase & fp, numeric::xyzVector<core::Real> const & CoM_offset, core::Real const & angle1_offset, core::Real const & angle2_offset, core::Real const & angle3_offset, core::Size const & conformer ) {

	// note: conformer passed in is indexed from 0, below is indexed from 1
  core::pose::PoseOP tmp_pose = new core::pose::Pose(pose_, conformer+1, conformer+1);

  apply_rotation_offset_to_pose_( *tmp_pose, angle1_offset, angle2_offset, angle3_offset );

  // set the fingerprint CoM to be the ligand CoM, then apply offset as needed
  CoM_ = calculate_ligand_CoM(*tmp_pose);
  // jk note: this next step of setting the origins to match is unnecessary as written, because given current implementation CoM_ and fp.CoM() are the same
  // jk is there a case when this isn't true, or is this leftover from an old approach??
  origin_ = fp.origin() + CoM_ - fp.CoM();
  origin_ += CoM_offset;

  core::Size const lig_res_num = compute_ligand_resnum(*tmp_pose);
  //	core::Size const ligand_natoms = compute_ligand_natoms(*tmp_pose);
  core::conformation::ResidueCOP ligand_rsd = new core::conformation::Residue( tmp_pose->conformation().residue(lig_res_num) );

  return ligand_rsd;

}

void PlaidFingerprint::update_rhos_(FingerprintBase & fp, core::conformation::ResidueCOP curr_ligand_rsd, bool const update_derivatives ) {

  using namespace basic::options;
  core::Real const radius_scale = option[ OptionKeys::fingerprint::atom_radius_scale ];
  core::Real const atom_buffer = option[ OptionKeys::fingerprint::atom_radius_buffer ];

  triplet_fingerprint_data_.clear();

  core::Size const ligand_natoms = compute_ligand_natoms();

  // compute the max and min possible phi and psi for each atom, store these in 4 arrays
  // also compute the overall max and min possible for the whole ligand
  utility::vector1<core::Real> atom_max_phi( ligand_natoms );
  utility::vector1<core::Real> atom_max_psi( ligand_natoms );
  utility::vector1<core::Real> atom_min_phi( ligand_natoms );
  utility::vector1<core::Real> atom_min_psi( ligand_natoms );

  utility::vector1< core::Real > atomX( ligand_natoms );
  utility::vector1< core::Real > atomY( ligand_natoms );
  utility::vector1< core::Real > atomZ( ligand_natoms );
  utility::vector1<core::Real> atom_radius( ligand_natoms );

  core::Real max_phi( -999. );
  core::Real max_psi( -999. );
  core::Real min_phi( 999. );
  core::Real min_psi( 999. );

  for (Size i = 1, i_end = ligand_natoms; i <= i_end; ++i) {

    numeric::xyzVector<core::Real> this_atomcoors = curr_ligand_rsd->atom(i).xyz() - origin_;
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
    if ( ( std::abs(tmp_atomx) > 0.00001 ) || ( std::abs(tmp_atomx) > 0.00001 ) ) {
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

  derivs_of_ray_distances_.clear();

  //float orig_cpu_total_dist = 0.;
  //Size orig_cpu_total_num = 0;
  //Size orig_cpu_num_evaluations = 0;
  for (std::list<spherical_coor_triplet>::const_iterator ni = fp.triplet_fingerprint_data().begin(); ni != fp.triplet_fingerprint_data().end(); ++ni) {

    core::Real curr_phi = ni->phi;
    core::Real curr_psi = ni->psi;
    // if the current phi and/or psi is outside the overall max/min, set best_rho to zero and jumpout (ie. ray misses the ligand)
    core::Real best_rho_sq(9999.);
    core::Size best_intersecting_atom(0);
    for (Size i = 1, i_end = ligand_natoms; i <= i_end; ++i) {

      if ( atom_radius.at(i) < 0.001 ) continue;

      while ( curr_phi < atom_min_phi.at(i) ) {
	curr_phi += numeric::constants::r::pi_2;
      }
      while ( curr_phi > atom_max_phi.at(i) ) {
	curr_phi -= numeric::constants::r::pi_2;
      }
      while ( curr_psi < atom_min_psi.at(i) ) {
	curr_psi += numeric::constants::r::pi_2;
      }
      while ( curr_psi > atom_max_psi.at(i) ) {
	curr_psi -= numeric::constants::r::pi_2;
      }
      if ( curr_phi < atom_min_phi.at(i) ) continue;
      if ( curr_psi < atom_min_psi.at(i) ) continue;
      if ( curr_phi > atom_max_phi.at(i) ) continue;
      if ( curr_psi > atom_max_psi.at(i) ) continue;

      core::Real const min_intersect_SQ = Find_Closest_Intersect_SQ(curr_phi,curr_psi,atomX.at(i),atomY.at(i),atomZ.at(i),atom_radius.at(i));
      //++orig_cpu_num_evaluations;

      if ( min_intersect_SQ < best_rho_sq ) {
	best_rho_sq = min_intersect_SQ;
	best_intersecting_atom = i;
      }
    }

    spherical_coor_triplet new_triplet;
    new_triplet.phi = curr_phi;
    new_triplet.psi = curr_psi;

    if ( best_rho_sq < 9998. ) {
      new_triplet.rho = sqrt(best_rho_sq);
      //orig_cpu_total_dist += sqrt(best_rho_sq);
      //++orig_cpu_total_num;
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

	//	cout << "ORIG_CPU TOTAL DISTANCE: " << orig_cpu_total_dist << std::endl;
	//	cout << "ORIG_CPU TOTAL NUM: " << orig_cpu_total_num << std::endl;
	//	cout << "ORIG_CPU NUM INTERSECTION EVALUATIONS: " << orig_cpu_num_evaluations << std::endl;

}


core::Real PlaidFingerprint::fp_compare( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) const {

  core::Real Total_score = 0;
  core::Size num_rays = 0;
  using namespace basic::options;
  bool square  = option[ OptionKeys::fingerprint::square_score ]();

  for (std::list<spherical_coor_triplet>::const_iterator pi = fp.triplet_fingerprint_data().begin(), li = triplet_fingerprint_data_.begin(); pi != fp.triplet_fingerprint_data().end(), li != triplet_fingerprint_data_.end(); ++pi, ++li) {

    // jk note: these are no longer necessarily true, since we alter the Plaid one by shifting up/down by 2*pi
    // these are useful asserts though, it's worth thinking about how to make them valid again...
    // assert( std::abs( pi->phi - li->phi ) < 0.001 );
    // assert( std::abs( pi->psi - li->psi ) < 0.001 );

    if ( (li->rho < 0.001) && (pi->rho < 0.001) ) {
      continue;
    } else if (li->rho < 0.001) {
      Total_score += missing_point_weight;
    } else if (pi->rho < 0.001 ) {
      Total_score += extra_point_weight;
    } else {
      core::Real dist_deviation = std::abs( pi->rho - li->rho );
      if (square){
	if (li->rho > pi->rho) dist_deviation *= dist_deviation;
	if (li->rho < pi->rho) dist_deviation *= (dist_deviation * steric_weight);
      }
      else if (!square) {
	if (li->rho > pi->rho) dist_deviation = dist_deviation;
	if (li->rho < pi->rho) dist_deviation = (dist_deviation * steric_weight);
      }
      Total_score += dist_deviation;
    }
    num_rays++;
  }
  //	if ( num_rays < 25 ) return 999.;
  return (Total_score/num_rays);
}

void PlaidFingerprint::fp_compare_deriv( FingerprintBase & fp, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, core::Real & dE_dx, core::Real & dE_dy, core::Real & dE_dz, core::Real & dE_dv4, core::Real & dE_dv5, core::Real & dE_dv6 ) const {

  dE_dx = 0.; dE_dy = 0.; dE_dz = 0.; dE_dv4 = 0.; dE_dv5 = 0.; dE_dv6 = 0.;

  if ( derivs_of_ray_distances_.size() < 2 ) {
    std::cout<<"Error, fingerprint derivatives have not been computed" << std::endl;
    exit(1);
  }
  assert( derivs_of_ray_distances_.size() == fp.triplet_fingerprint_data().size() );

  core::Real Total_score = 0;
  core::Real Differentiable_score = 0;
  core::Size num_rays = 0;
  std::list<ray_distance_derivs>::const_iterator di = derivs_of_ray_distances_.begin();
  for (std::list<spherical_coor_triplet>::const_iterator pi = fp.triplet_fingerprint_data().begin(), li = triplet_fingerprint_data_.begin(); pi != fp.triplet_fingerprint_data().end(), li != triplet_fingerprint_data_.end(), di != derivs_of_ray_distances_.end(); ++pi, ++li, ++di) {
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
      // derivative is zero except in the case where the ray hits BOTH ligand and pocket
      // (ie. "missing point" and "extra point" scores don't contribute to the derivatives)
      if (li->rho > pi->rho) {
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
  //	std::cout<<"DARC score while computing derivatives are total: " << Total_score << " and differentiable part " << Differentiable_score << std::endl;
  //	std::cout<<"Derivatives are " << dE_dx << " , " << dE_dy << " , " << dE_dz << " , " << dE_dv4 << " , " << dE_dv5 << " , " << dE_dv6 << std::endl;

  return;

}

core::Real PlaidFingerprint::search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight ) {

  numeric::xyzVector<core::Real> no_CoM_offset(0.);
  return search_random_poses( fp, num_pose_search, optimal_angle1, optimal_angle2, optimal_angle3, missing_point_weight, steric_weight, extra_point_weight, no_CoM_offset);
}


core::Real PlaidFingerprint::search_random_poses( FingerprintBase & fp, core::Size const & num_pose_search, core::Real & optimal_angle1, core::Real & optimal_angle2, core::Real & optimal_angle3, core::Real const & missing_point_weight, core::Real const & steric_weight, core::Real const & extra_point_weight, numeric::xyzVector<core::Real> const & CoM_offset ) {

  core::Real best_score = std::numeric_limits<core::Real>::max();

  for (core::Size j = 0; j < num_pose_search; ++j ){
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
  core::Size num_steps = core::Size ( 360. / angle_increment );

  core::Real curr_angle1=0.;
  for (core::Size i = 0; i < num_steps; ++i ){
    core::Real curr_angle2=0.;
    for (core::Size j = 0; j < num_steps; ++j ){
      core::Real curr_angle3=0.;

      for (core::Size k = 0; k < num_steps; ++k ){
				std::cout<< "JK this code is not yet conformer-enabled, fix it in the app by removing the zero in the call to move_ligand_and_update_rhos_ below..." << std::endl;
				exit(1);
				move_ligand_and_update_rhos_( fp, CoM_offset, curr_angle1, curr_angle2, curr_angle3, 0 );
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
  for(Size i = 1, i_end = curr_rsd.natoms(); i <= i_end; ++i) {
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

  move_ligand_and_update_rhos_( fp, CoM_offset, angle1_offset, angle2_offset, angle3_offset, conformer );

  core::pose::Pose tmp_pose = pose_;
  apply_rotation_offset_to_pose_( tmp_pose, angle1_offset, angle2_offset, angle3_offset );

  numeric::xyzMatrix<core::Real> bestx_rot_mat( numeric::x_rotation_matrix_radians( original_pocket_angle_transform[1] ) );
  numeric::xyzMatrix<core::Real> besty_rot_mat( numeric::y_rotation_matrix_radians( original_pocket_angle_transform[2] ) );
  numeric::xyzMatrix<core::Real> bestz_rot_mat( numeric::z_rotation_matrix_radians( original_pocket_angle_transform[3] ) );
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

numeric::xyzVector<core::Real> PlaidFingerprint::calculate_ligand_CoM(core::pose::Pose const & ligand_pose) {

  core::Size lig_res_num = compute_ligand_resnum( ligand_pose );
  numeric::xyzVector<core::Real> ligand_com(0.);
  core::conformation::Residue const & curr_rsd = ligand_pose.conformation().residue(lig_res_num);
  core::Size ligand_total_atoms = compute_ligand_natoms( ligand_pose );
  for(Size i = 1, i_end = ligand_total_atoms; i <= i_end; ++i) {
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
  for(Size i = 1, i_end = pose1_rsd.nheavyatoms(); i <= i_end; ++i) {
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
//	numeric::xyzVector<core::Real> fp_coor;
//  for (std::list<spherical_coor_triplet>::iterator pd = triplet_fingerprint_data_.begin(); pd != triplet_fingerprint_data_.end(); ++pd) {
//		// find cartesian coors relative to old origin
//		convert_spherical_coor_triplet_to_cartesian( *pd, fp_coor );
//		// move cartesian coors from old origin to new origin
//    fp_coor += origin_ - new_origin;
//		// convert from cartesian coors to polar
//		convert_cartesian_to_spherical_coor_triplet( fp_coor, *pd );
//	}
//	origin_ = new_origin;
//}



// note: not used, also deprecated because we now express all angles in radians
//void correct_phi_psi( core::Real & phi, core::Real & psi ){
//	while ( phi < 0. ) {
//		phi *= -1.;
//		psi += 180.;
//	}
//	while ( phi > 360. ) {
//		phi -= 360.;
//	}
//	while ( phi > 180. ) {
//		phi = 360. - phi;
//		psi += 180.;
//	}
//	while ( psi < -180. ) psi += 360.;
//	while ( psi > 180. ) psi -= 360.;
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
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = xyz_list_1.begin(); pd != xyz_list_1.end(); ++pd) {
    combined_list.push_back(*pd);
  }
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = xyz_list_2.begin(); pd != xyz_list_2.end(); ++pd) {
    combined_list.push_back(*pd);
  }
  return combined_list;
}


//round triplet values
std::list<spherical_coor_triplet> NonPlaidFingerprint::convert_cart_to_spherical_and_round(std::list< numeric::xyzVector<core::Real> > const & xyz_list) {
	std::list< spherical_coor_triplet > rounded_triplet_list;
	rounded_triplet_list.clear();
	spherical_coor_triplet ray_triplet;
  for (std::list< numeric::xyzVector<core::Real> >::const_iterator pd = xyz_list.begin(); pd != xyz_list.end(); ++pd) {
    convert_cartesian_to_spherical_coor_triplet( *pd - origin_, ray_triplet );
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
	for (std::list<spherical_coor_triplet>::const_iterator aa = spherical_triplet_list.begin(); aa != spherical_triplet_list.end(); ++aa) {
		ray_triplet = *aa;
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

	for (std::list<spherical_coor_triplet>::const_iterator aa = rounded_triplet.begin(); aa != rounded_triplet.end(); ++aa) {
		bool found = false;
		spherical_coor_triplet best_triplet = *aa;
		for (std::list<spherical_coor_triplet>::iterator bb = temp_triplet.begin(); bb != temp_triplet.end();) {
			if( (aa->phi == bb->phi) && (aa->psi == bb->psi) ) {
				found = true;
				if (bb->rho < best_triplet.rho) { best_triplet = *bb;}
				bb = temp_triplet.erase(bb);
			}
			else {
				++bb;
			}
		}
		if(found) unique_triplet.push_back(best_triplet);
	}
	return unique_triplet;
}

std::list<numeric::xyzVector<core::Real> > NonPlaidFingerprint::convert_spherical_list_to_cartesian_list(std::list<spherical_coor_triplet> const & unique_triplet) {
  std::list<numeric::xyzVector<core::Real> > xyz_list;
  xyz_list.clear();
  for (std::list<spherical_coor_triplet>::const_iterator pd = unique_triplet.begin(); pd != unique_triplet.end(); ++pd) {
    numeric::xyzVector<core::Real> new_coor;
    convert_spherical_coor_triplet_to_cartesian( *pd, new_coor );
    new_coor += origin_;
		xyz_list.push_back(new_coor);
  }
  return xyz_list;
}

} // Pockets
} // protocols
