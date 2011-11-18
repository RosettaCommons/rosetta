// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Dinculeotide_Sampler_Util
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_Dinucleotide_Sampler_Util.hh>
//////////////////////////////////
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>

#include <core/pose/util.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
//////////////////////////////////////#include <set>
#include <numeric/conversions.hh>
#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>

#include <ObjexxFCL/string.functions.hh>

#define PI 3.14159265
//#define centroid_bin_size 2
//#define euler_angle_bin_size 45
//#define euler_z_bin_size 0.5

#define centroid_bin_size 1
#define euler_angle_bin_size 20
#define euler_z_bin_size 0.05


using namespace core;

static basic::Tracer TR( "protocols.swa.stepwise_rna_dinucleotide_sampler_util" ) ;

namespace protocols {
namespace swa {
namespace rna {

	Anchor_ribose_stub
	Get_anchor_ribose_stub(conformation::Residue const & rsd, bool verbose){

  	using namespace chemical;

		if(verbose){ //Will implement verbose output later
		}

		Anchor_ribose_stub anchor_ribose_stub;

  	assert( rsd.is_RNA() );

		Vector x,y,z;

		//C4* atom is the origin of the coordinate system
		anchor_ribose_stub.origin=rsd.xyz(" C4*");
		Vector const & origin = anchor_ribose_stub.origin;

		Vector const H4_coord =rsd.xyz(" H4*");
		x = H4_coord - origin;
		x.normalize();

  	Vector const C5_coord =rsd.xyz(" C5*");
  	y = C5_coord - origin; //not orthonormal yet...
		z = cross(x, y);  //Can cross here even though y is not orthogonal to x since the component of y that is parallel to x dissapear when crossed with x since cross(x, x)=0
		z.normalize(); //Choosen H4 and C5 atom specifically so that z will be roughly parallel with the helical axis.

  	y = cross(z, x); //Now y is orthonormal
  	y.normalize(); //not necessary but doesn't hurt.

	 	anchor_ribose_stub.coordinate_matrix=Matrix::cols( x, y, z );
	 	anchor_ribose_stub.invert_coordinate_matrix=inverse(anchor_ribose_stub.coordinate_matrix);

		return anchor_ribose_stub;

	}


	numeric::xyzVector<core::Real> //This is equivalent to Vector
	Convert_base_centroid_to_anchor_ribose_coord_system( numeric::xyzVector<core::Real> const & centroid, Anchor_ribose_stub const & anchor_ribose_stub){

		//Translate origin to C4 atom and then rotate to the anchor_ribose coordinate frame by applying the invert_coordinate matrix

		numeric::xyzVector<core::Real> translated_centroid=centroid - anchor_ribose_stub.origin;

		 return (anchor_ribose_stub.invert_coordinate_matrix * translated_centroid);

//		return product( anchor_ribose_stub.invert_coordinate_matrix, translated_centroid);

	}

	numeric::xyzMatrix< core::Real >
	Convert_base_coordinate_matrix_to_anchor_ribose_coord_system(numeric::xyzMatrix< core::Real > const & base_coordinate_matrix, Anchor_ribose_stub const & anchor_ribose_stub){

		//rotate to the H4_C4_C5 coordinate frame  applying the invert_coordinate matrix

		return (anchor_ribose_stub.invert_coordinate_matrix * base_coordinate_matrix);

	}

	Euler_angles
	Get_euler_angles( numeric::xyzMatrix< core::Real > const & base_coordinate_matrix_in_anchor_ribose_coord_system){

		numeric::xyzMatrix< core::Real > const & M = base_coordinate_matrix_in_anchor_ribose_coord_system;

		Euler_angles euler_angles;
		euler_angles.alpha= atan2(M.xz(),-M.yz()) * 180 / PI ;
		euler_angles.z= M.zz();
		euler_angles.gamma= atan2(M.zx(), M.zy()) * 180 / PI ;  //tan2(y,x)=gamma

/* Defination of rotation matrix at mathworld is actually the inverse of the rotation matrix....also did not correctly input the coefficient of arctan2
		euler_angles.alpha = (-1)* atan2(M.zy(), M.zx()) * 180 / PI; //Is this correct?? Ambuiguity with the minus sign...return in the [-Pi, Pi] range. atan2(y-value, x-value)
		euler_angles.beta  = acos(M.zz()) * 180 / PI;  // rerun in the [0,Pi] range
		euler_angles.z= M.zz(); //Use z instead of beta to make space uniform
		euler_angles.gamma = atan2(M.yz(), M.xz()) * 180 / PI; //return in the [-Pi, Pi] range. atan2(y-value, x-value)
*/

		return euler_angles;
	}


	Base_bin
	Get_base_bin(numeric::xyzVector<core::Real> const & centriod, Euler_angles const & euler_angles){

//		using namespace Bin_size;

		Base_bin base_bin;
		base_bin.centroid_x=static_cast<int>( centriod[0]/static_cast<Real>(centroid_bin_size) ); //Range roughly [-14:14] for dinucleotide case
		base_bin.centroid_y=static_cast<int>( centriod[1]/static_cast<Real>(centroid_bin_size) ); //Range roughly [-14:14] for dinucleotide case
		base_bin.centroid_z=static_cast<int>( centriod[2]/static_cast<Real>(centroid_bin_size) ); //Range roughly [-14:14] for dinucleotide case

		base_bin.euler_alpha=static_cast<int>(euler_angles.alpha/euler_angle_bin_size);
		base_bin.euler_gamma=static_cast<int>(euler_angles.gamma/euler_angle_bin_size);
		base_bin.euler_z=static_cast<int>(euler_angles.z/euler_z_bin_size);


		if(centriod[0]<0) base_bin.centroid_x--;
		if(centriod[1]<0) base_bin.centroid_y--;
		if(centriod[2]<0) base_bin.centroid_z--;
		if(euler_angles.alpha<0) base_bin.euler_alpha--;
		if(euler_angles.gamma<0) base_bin.euler_gamma--;
		if(euler_angles.z<0) base_bin.euler_z--;

		return base_bin;
	}


	int
	DOF_bin_value(std::map<Base_bin , int , compare_base_bin>::const_iterator const & base_bin_it, std::string const & DOF){

		if(DOF=="x"){
			return base_bin_it->first.centroid_x;
		}else if(DOF=="y"){
			return base_bin_it->first.centroid_y;
		}else if(DOF=="z"){
			return base_bin_it->first.centroid_z;
		}else if(DOF=="alpha"){
			return base_bin_it->first.euler_alpha;
		}else if(DOF=="euler_z"){
			return base_bin_it->first.euler_z;
		}else if(DOF=="gamma"){
			return base_bin_it->first.euler_gamma;
		}else{
			utility_exit_with_message( "Invalid DOF= " + DOF);
			exit(1); //prevent compiler warning
		}
	}

	Real
	DOF_bin_size(std::string const & DOF){

		if(DOF=="x" || DOF=="y" || DOF=="z"){
			return centroid_bin_size;
		}else if(DOF=="alpha" || DOF=="gamma"){
			return euler_angle_bin_size;
		}else if(DOF=="euler_z"){
			return euler_z_bin_size;
		}else{
			utility_exit_with_message( "Invalid DOF= " + DOF);
			exit(1); //prevent compiler warning
		}
	}

	void
	Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const & DOF_one, std::string const & DOF_two, std::string const foldername){


		std::map<std::pair<int, int> , int , compare_int_pair> count_density_map;
		std::map<std::pair<int, int> , int , compare_int_pair>::iterator count_density_it;

		std::map<Base_bin , int , compare_base_bin>::const_iterator base_bin_it;

		int total_count=0;
		int total_occupied_bin=0;

		for (base_bin_it=base_bin_map.begin(); base_bin_it!=base_bin_map.end(); base_bin_it++ ){

			total_occupied_bin++;
			total_count=total_count+base_bin_it->second;

			std::pair< int, int > const & DOF_pair=std::make_pair(DOF_bin_value(base_bin_it, DOF_one), DOF_bin_value(base_bin_it, DOF_two));

			count_density_it=count_density_map.find(DOF_pair);

			if(count_density_it==count_density_map.end()){
				count_density_map[DOF_pair]=1;
			}else{
				count_density_it->second++;
			}
		}

		//////////////////////Output data/////////////////////////////////////////////////////////////////////////////
		std::ofstream outfile;
		std::string filename=foldername + "Bin_" + DOF_one + "_" + DOF_two + ".txt";
		outfile.open(filename.c_str());
		Size const spacing=14;

		outfile << std::setw(spacing) << DOF_one;
		outfile << std::setw(spacing) << DOF_two;
		outfile << std::setw(30) << "occupied_bin_count";
		outfile << "\n";

		int DOF_one_bin_max=0;
		int DOF_one_bin_min=0;
		int DOF_two_bin_max=0;
		int DOF_two_bin_min=0;

		for (count_density_it=count_density_map.begin(); count_density_it!=count_density_map.end(); count_density_it++ ){
			int const & DOF_one_bin_value=count_density_it->first.first;
			int const & DOF_two_bin_value=count_density_it->first.second;

			if(DOF_one_bin_value>DOF_one_bin_max) DOF_one_bin_max=DOF_one_bin_value;
			if(DOF_two_bin_value>DOF_two_bin_max) DOF_two_bin_max=DOF_two_bin_value;
			if(DOF_one_bin_value<DOF_one_bin_min) DOF_one_bin_min=DOF_one_bin_value;
			if(DOF_two_bin_value<DOF_two_bin_min) DOF_two_bin_min=DOF_two_bin_value;

		}

		for(int DOF_one_bin_value=(DOF_one_bin_min-5); DOF_one_bin_value<(DOF_one_bin_max+5); DOF_one_bin_value++){
			for(int DOF_two_bin_value=(DOF_two_bin_min-5); DOF_two_bin_value<(DOF_two_bin_max+5); DOF_two_bin_value++){

				Real const DOF_one_value=DOF_one_bin_value*DOF_bin_size(DOF_one);
				Real const DOF_two_value=DOF_two_bin_value*DOF_bin_size(DOF_two);

				int occupied_bin_count;
				std::pair< int, int > const & DOF_pair=std::make_pair(DOF_one_bin_value, DOF_two_bin_value);
				count_density_it=count_density_map.find(DOF_pair);

				if(count_density_it==count_density_map.end()){
					occupied_bin_count=0;
				}else{
					occupied_bin_count=count_density_it->second;
				}

				outfile << std::setw(spacing) << DOF_one_value;
				outfile << std::setw(spacing) << DOF_two_value;
				outfile << std::setw(spacing) << occupied_bin_count;
				outfile << "\n";
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////

		std::cout << std::setw(50) << std::left << "Analysis " + DOF_one + "_" + DOF_two;
		std::cout << " tot_count = " << std::setw(15) << std::left << total_count << " tot_occ= " << std::setw(15) << std::left << total_occupied_bin;
		std::cout << " tot_count/tot_occ_bin= " << std::setw(5) << std::left << (double(total_count)/double(total_occupied_bin)) << std::endl;


	}

	void
	Analyze_base_bin_map(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, std::string const foldername){

		Analyze_base_bin_map(base_bin_map, "x", "y", foldername);
		Analyze_base_bin_map(base_bin_map, "x", "z", foldername);
		Analyze_base_bin_map(base_bin_map, "x", "alpha", foldername);
		Analyze_base_bin_map(base_bin_map, "x", "euler_z", foldername);
		Analyze_base_bin_map(base_bin_map, "x", "gamma", foldername);

		Analyze_base_bin_map(base_bin_map, "y", "z", foldername);
		Analyze_base_bin_map(base_bin_map, "y", "alpha", foldername);
		Analyze_base_bin_map(base_bin_map, "y", "euler_z", foldername);
		Analyze_base_bin_map(base_bin_map, "y", "gamma", foldername);

		Analyze_base_bin_map(base_bin_map, "z", "alpha", foldername);
		Analyze_base_bin_map(base_bin_map, "z", "euler_z", foldername);
		Analyze_base_bin_map(base_bin_map, "z", "gamma", foldername);

		Analyze_base_bin_map(base_bin_map, "alpha", "euler_z", foldername);
		Analyze_base_bin_map(base_bin_map, "alpha", "gamma", foldername);

		Analyze_base_bin_map(base_bin_map, "euler_z", "gamma", foldername);

//		std::cout << "Is_dinucleotide= " << Is_dinucleotide << std::endl;
		std::cout << "centroid_bin_size= " << centroid_bin_size << "  euler_angle_bin_size= " <<  euler_angle_bin_size << "  euler_z_bin_size= " << euler_z_bin_size << std::endl;


	}


	void
	center(pose::Pose & pose, Size const seq_num)
	{

		using namespace core::chemical;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;

		numeric::xyzVector<core::Real> centroid=get_base_centroid(  pose.residue( seq_num ) , true);



		numeric::xyzMatrix< core::Real > base_coordinate_matrix=get_base_coordinate_system( pose.residue( seq_num ) , centroid );

		numeric::xyzMatrix< core::Real > invert_coordinate_matrix= inverse(base_coordinate_matrix);

		Size count=0;
		for( Size at = 1; at <= pose.residue( seq_num ).natoms(); at++){
			std::string const & atom_name=pose.residue( seq_num ).type().atom_name(at);
			std::cout << "atom_name= " << atom_name << std::endl;
//			AtomID id( j, i )
			id::AtomID const id( at, seq_num);

//			AtomID & id= AtomID( at, seq_num )
			pose.set_xyz( id, pose.xyz( id) - centroid );//I think the order here does matter. Translate centroid to origin.
			pose.set_xyz( id, invert_coordinate_matrix * pose.xyz(id)); //I think the order here does matter. Rotate coordinate so that it equal to Roseeta internal reference frame

		}


	}



	void
	Sample_base_test(){

		using namespace core::chemical;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;


		pose::Pose pose;

  	ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		std::string sequence = "a";

		make_pose_from_sequence( pose, sequence, *rsd_set, false /*auto_termini*/);

		pose.dump_pdb( "Single_nucleotide.pdb" );


		center(pose, 1);

		pose.dump_pdb( "Single_nucleotide_after_center.pdb" );

		int const euler_angle_bin_min=-180/euler_angle_bin_size; //Should be -180/euler_angle_bin_size
		int const euler_angle_bin_max=180/euler_angle_bin_size-1;  //Should be 180/euler_angle_bin_size-1

		int const euler_z_bin_min=-1/euler_z_bin_size;
		int const euler_z_bin_max=(1/euler_z_bin_size)-1;

		Real const max_distance=12;
		int const centroid_bin_min=-max_distance/static_cast<Real>(centroid_bin_size);
		int const centroid_bin_max=(max_distance/static_cast<Real>(centroid_bin_size))-1;

		std::cout << "euler_angle_bin min= " << euler_angle_bin_min << " max " << euler_angle_bin_max << std::endl;
		std::cout << "euler_z_bin_min min= " << euler_z_bin_min << " max " << euler_z_bin_max << std::endl;
		std::cout << "centroid_bin min= " << centroid_bin_min<< " max " << centroid_bin_max << std::endl;

		std::map<Base_bin , int , compare_base_bin> base_bin_map;
		std::map<Base_bin , int , compare_base_bin>::const_iterator it;


		SillyCountStruct count_data;

		base_stub reference_res_base;
		reference_res_base.centroid=get_base_centroid(  pose.residue( 1 ) , true);
		reference_res_base.base_coordinate_matrix=get_base_coordinate_system( pose.residue( 1 ) , reference_res_base.centroid);
		utility::vector1 < base_stub > other_residues_base_list;
		other_residues_base_list.push_back(reference_res_base);

		base_stub moving_res_base;
		numeric::xyzVector<Real> & centroid = moving_res_base.centroid;
		numeric::xyzMatrix< core::Real > & base_coordinate_matrix= moving_res_base.base_coordinate_matrix;

		Base_bin base_bin;

		Size total_bin(0), total_screen_bin(0);

		for(base_bin.centroid_x=centroid_bin_min; base_bin.centroid_x<=centroid_bin_max; base_bin.centroid_x++){
			centroid[0]=(base_bin.centroid_x+0.5)*centroid_bin_size;

		for(base_bin.centroid_y=centroid_bin_min; base_bin.centroid_y<=centroid_bin_max; base_bin.centroid_y++){
			centroid[1]=(base_bin.centroid_y+0.5)*centroid_bin_size;

		for(base_bin.centroid_z=centroid_bin_min; base_bin.centroid_z<=centroid_bin_max; base_bin.centroid_z++){
			centroid[2]=(base_bin.centroid_z+0.5)*centroid_bin_size;

		for(base_bin.euler_alpha=euler_angle_bin_min; base_bin.euler_alpha<=euler_angle_bin_max; base_bin.euler_alpha++){
			Real const alpha=(base_bin.euler_alpha+0.5)*euler_angle_bin_size*(PI/180); //convert to radians

		for(base_bin.euler_z=euler_z_bin_min; base_bin.euler_z<=euler_z_bin_max; base_bin.euler_z++){
			Real const euler_z=(base_bin.euler_z+0.5)*euler_z_bin_size;
			Real const beta=acos(euler_z);

		for(base_bin.euler_gamma=euler_angle_bin_min; base_bin.euler_gamma<=euler_angle_bin_max; base_bin.euler_gamma++){ //Potential optimization here, don't really need gamma for base centroid screening
			Real const gamma=(base_bin.euler_gamma+0.5)*euler_angle_bin_size*(PI/180); //convert to radians

			total_bin++;

//Probably could save time if determine x and y and take cross product to determine z
//Should determine using both ways and check for consistency
			base_coordinate_matrix.xx( cos(alpha)*cos(gamma) - sin(alpha)*cos(beta)*sin(gamma));
			base_coordinate_matrix.xy(-cos(alpha)*sin(gamma) - sin(alpha)*cos(beta)*cos(gamma));
			base_coordinate_matrix.xz( sin(alpha)*sin(beta));
			base_coordinate_matrix.yx( sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*sin(gamma));
			base_coordinate_matrix.yy(-sin(alpha)*cos(gamma) + cos(alpha)*cos(beta)*cos(gamma));
			base_coordinate_matrix.yz(-cos(alpha)*sin(beta));
			base_coordinate_matrix.zx( sin(beta) *sin(gamma));
			base_coordinate_matrix.zy( sin(beta) *cos(gamma));
			base_coordinate_matrix.zz( cos(beta));
/*
		xx_( m.xx_ ), xy_( m.xy_ ), xz_( m.xz_ ),
		yx_( m.yx_ ), yy_( m.yy_ ), yz_( m.yz_ ),
		zx_( m.zx_ ), zy_( m.zy_ ), zz_( m.zz_ )
*/

			if(Base_centroid_screening(moving_res_base, other_residues_base_list, count_data)==false) continue;

			it=base_bin_map.find(base_bin);

			if(it==base_bin_map.end()){
				base_bin_map[base_bin]=1;
				total_screen_bin++;
			}else{
				base_bin_map[base_bin]=base_bin_map[base_bin]+1;
			}

		}
		}
		}
		}
		}
		}

		std::string const foldername="test/";
		system(std::string("rm -r " + foldername).c_str());
		system(std::string("mkdir -p " + foldername).c_str());

		Analyze_base_bin_map( base_bin_map, foldername);

		std::cout << "total_bin= " << total_bin << " total_screen_bin= " << total_screen_bin;
		std::cout << "  stack= " << count_data.base_stack_count << " pair= " << count_data.base_pairing_count;

	}

/*
	void
	Analyze_base_bin_map_old(std::map<Base_bin , int , compare_base_bin> const & base_bin_map, bool const Is_dinucleotide){

	  std::cout << "Analyze_base_bin_map function " << std::endl;
//		using namespace Bin_size;

		Real const max_distance= (Is_dinucleotide) ?  22 : 14;


		int const euler_angle_bin_min=-180/euler_angle_bin_size; //Should be -180/euler_angle_bin_size
		int const euler_angle_bin_max=180/euler_angle_bin_size-1;  //Should be 180/euler_angle_bin_size-1

		int const euler_z_bin_min=-1/euler_z_bin_size;
		int const euler_z_bin_max=(1/euler_z_bin_size)-1;

		int const centroid_bin_min=-max_distance/static_cast<Real>(centroid_bin_size);
		int const centroid_bin_max=(max_distance/static_cast<Real>(centroid_bin_size))-1;


		std::cout << "euler_angle_bin min= " << euler_angle_bin_min << " max " << euler_angle_bin_max << std::endl;
		std::cout << "euler_z_bin_min min= " << euler_z_bin_min << " max " << euler_z_bin_max << std::endl;
		std::cout << "centroid_bin min= " << centroid_bin_min<< " max " << centroid_bin_max << std::endl;


		std::ofstream outfile;
		std::string filename="Bin.txt";
		outfile.open(filename.c_str());

		Size const spacing=8;

		std::map<Base_bin , int , compare_base_bin>::const_iterator it;

		Base_bin base_bin;

		int total_count=0;
		int total_bin=0;
		int total_bin_actual=0;
		int total_occupied_bin=0;


		for(base_bin.centroid_x=centroid_bin_min; base_bin.centroid_x<=centroid_bin_max; base_bin.centroid_x++){
			Real const x=base_bin.centroid_x*centroid_bin_size;

		for(base_bin.centroid_y=centroid_bin_min; base_bin.centroid_y<=centroid_bin_max; base_bin.centroid_y++){
			Real const y=base_bin.centroid_y*centroid_bin_size;

			int occupied_bin_count=0;
			int empty_bin_count=0;

		for(base_bin.euler_alpha=euler_angle_bin_min; base_bin.euler_alpha<=euler_angle_bin_max; base_bin.euler_alpha++){
			Real const alpha=base_bin.euler_alpha*euler_angle_bin_size;

		for(base_bin.euler_gamma=euler_angle_bin_min; base_bin.euler_gamma<=euler_angle_bin_max; base_bin.euler_gamma++){
			Real const gamma=base_bin.euler_gamma*euler_angle_bin_size;

		for(base_bin.euler_z=euler_z_bin_min; base_bin.euler_z<=euler_z_bin_max; base_bin.euler_z++){
			Real const euler_z=base_bin.euler_z*euler_z_bin_size;

		for(base_bin.centroid_z=centroid_bin_min; base_bin.centroid_z<=centroid_bin_max; base_bin.centroid_z++){
			Real const z=base_bin.centroid_z*centroid_bin_size;

			int count=0;

			total_bin++;
			it=base_bin_map.find(base_bin);

			if(it!=base_bin_map.end()){
				count=count+it->second;
				total_count=total_count+it->second;
			}


		if(count==0) empty_bin_count++;
		if(count!=0) occupied_bin_count++;

		total_bin_actual++;
		if(count!=0) total_occupied_bin++;

		}
		}
		}
		}

			outfile << std::setw(spacing) << x;
			outfile << std::setw(spacing) << y;
//			outfile << std::setw(spacing) << z;
//			outfile << std::setw(spacing) << alpha;
//			outfile << std::setw(spacing) << gamma;
//			outfile << std::setw(spacing) << euler_z;
			outfile << std::setw(spacing) << empty_bin_count;
			outfile << std::setw(spacing) << occupied_bin_count;
			outfile << "\n";
		}
		}

		std::cout << "Base_bin_map Analysis " << std::endl;
		std::cout << "total_count = " << total_count << " total_bin= " << total_bin;
		std::cout << " total_bin_actual= " << total_bin_actual << "  total_occupied_bin= " << total_occupied_bin;
		std::cout << " occupied_bin fraction= " << (double(total_occupied_bin)/double(total_bin_actual));
		std::cout << " total_count/total_occupied_bin= " << (double(total_count)/double(total_occupied_bin)) << std::endl;
		std::cout << "centroid_bin_size= " << centroid_bin_size << "  euler_angle_bin_size= " <<  euler_angle_bin_size << "  euler_z_bin_size= " << euler_z_bin_size << std::endl;

//		std::ofstream outfile_old;
//		std::string filename_old="Bin_old.txt";
//		outfile_old.open(filename_old.c_str());

		int total_count_old=0;

//		for(int alpha_bin_value=euler_angle_bin_min; alpha_bin_value<=euler_angle_bin_max; alpha_bin_value++){
//		for(int gamma_bin_value=euler_angle_bin_min; gamma_bin_value<=euler_angle_bin_max; gamma_bin_value++){


//			Real const alpha=alpha_bin_value*euler_angle_bin_size;
//			Real const gamma=gamma_bin_value*euler_angle_bin_size;

//			int count=0;

			for (it=base_bin_map.begin(); it!=base_bin_map.end(); it++ ){

				Base_bin const & base_bin=it->first;

				total_count_old=total_count_old+it->second;

//				if(base_bin.euler_alpha==alpha_bin_value && base_bin.euler_gamma==gamma_bin_value) {
//					count=count+it->second;
//				}
			}

//			outfile_old << std::setw(spacing) << alpha;
//			outfile_old << std::setw(spacing) << gamma;
//			outfile_old << std::setw(spacing) << count;
//			outfile_old << "\n";

//		}
//		}
		std::cout << "Total_count_OLD(from base_bin_map)= " << total_count_old << std::endl;

	}




*/

	/*	//Note: Currently testing this function to build bulge...once this work will integrate back with Build_single_nucleotide, Dec 21, 2009 */
/*
	void
	StepWiseRNA_ResidueSampler::Build_dinucleotide( core::pose::Pose &  pose ) {

		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace protocols::rna;
		using namespace core::id;



		std::ofstream outfile;
		SilentFileData silent_file_data;
//		outfile.open(output_filename_.c_str(), std::ios_base::out | std::ios_base::app); //This does not delete existing content
		outfile.open(output_filename_.c_str());

//		outfile.open(output_filename_.c_str(), std::ios_base::out | std::ios_base::app); //This does not delete existing content

		Size const moving_res(  job_parameters_->working_moving_res() ); // corresponds to user input.
		Size const moving_suite(  job_parameters_->working_moving_suite() ); // dofs betweeen this value and value+1 actually move.
		bool const Is_prepend(  job_parameters_->Is_prepend() ); // if true, moving_suite+1 is fixed. Otherwise, moving_suite is fixed.
		bool const Is_internal(  job_parameters_->Is_internal() ); // no cutpoints before or after moving_res.
		Size const actually_moving_res( job_parameters_->actually_moving_res() );
		Size const gap_size( job_parameters_->gap_size());

		Size const  num_nucleotides(  job_parameters_->working_moving_res_list().size() );
		std::cout << " num_nucleotides= " <<  num_nucleotides << std::endl;
		bool const Is_dinucleotide=(num_nucleotides==2);
		std::cout << " Is_dinucleotide= "; Output_boolean(Is_dinucleotide); std::cout << std::endl;


		std::cout << " GAP SIZE " << gap_size << std::endl;
		std::cout << " MOVING RES " << moving_res << std::endl;
		std::cout << " MOVING SUITE " << moving_suite << std::endl;
		std::cout << " PREPEND " << Is_prepend << std::endl;
		std::cout << " INTERNAL " << Is_internal << std::endl;

		////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		Pose const pose_save = pose;
		pose = pose_save; //this recopy is useful for triggering graphics.

		//////////////// Sets up scorefunctions for a bunch of different screens ////////////////////////////////////////
		initialize_scorefunctions();



		//////////////////////////////Setup Bulge residue///////////////////////////////////////////////////////////////////////////

		Size const bulge_moving_suite = (Is_prepend) ? moving_suite+1: moving_suite-1;
		Size const bulge_moving_res= (Is_prepend) ? bulge_moving_suite : bulge_moving_suite+1;

		utility::vector1 <utility::vector1 <Real > > bulge_rotamer_list;

		if(Is_dinucleotide){
			if(bulge_moving_res!=job_parameters_->working_moving_res_list()[2]){
				utility_exit_with_message( "bulge_moving_res!=working_moving_res_list[2]" );
			}

			pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", bulge_moving_res);
//			if(Is_prepend){
				PuckerState second_res_puckerstate = Get_residue_pucker_state( pose, bulge_moving_suite+ 1);
				get_bulge_rotamers(bulge_rotamer_list, ALL, second_res_puckerstate, 1);
				get_bulge_rotamers(bulge_rotamer_list, ALL, second_res_puckerstate, 2);
				get_bulge_rotamers(bulge_rotamer_list, ALL, second_res_puckerstate, 3);
				get_bulge_rotamers(bulge_rotamer_list, ALL, second_res_puckerstate, 4);

//			}else{
//				PuckerState first_res_puckerstate = Get_residue_pucker_state( pose, bulge_moving_suite);
//				get_bulge_rotamers(bulge_rotamer_list, first_res_puckerstate, ALL, 2);
//			}
		}

		/////////////////////////////////////Dinucleotide_Sampler_Util//////////////////////////////////////////

		Size anchor_ribose_res;

		if(Is_dinucleotide){
			anchor_ribose_res= (Is_prepend) ? moving_suite+2 : moving_suite-2;
		}else{
			anchor_ribose_res= (Is_prepend) ? moving_suite+1 : moving_suite-1;
		}
		conformation::Residue const & anchor_ribose_rsd=pose.residue(anchor_ribose_res);
		Anchor_ribose_stub const anchor_ribose_stub=Get_anchor_ribose_stub( anchor_ribose_rsd );

		/////////////////////////////// O2star sampling/virtualization //////////////////////////
		Pose working_pose = pose;
		Add_virtual_O2Star_hydrogen( working_pose );
		working_pose.set_torsion( TorsionID( moving_res, id::CHI, 4 ), 0 );  //This torsion is not sampled. Arbitary set to zero to prevent randomness

		// if o2star_screen, pose 2'-OH will be sampled!
		if ( o2star_screen_ ) {
			if ( use_green_packer_ ){
				initialize_o2star_green_packer( pose );
			} else {
				initialize_o2star_packer_task( pose );
			}
		} else {
			// Otherwise, virtualize the 2-OH.
			pose = working_pose;
		}

		/////////////////////////////////////////Setup Base-stack/Base-pairing screening/////////////////////////////////////////
		utility::vector1 < base_stub > fixed_residues_base_list;
		utility::vector1 < Size > moving_positions;
		Initialize_base_centroid_screening( fixed_residues_base_list, moving_positions, working_pose);

		///////////////////////////////////////Setup chainbreak_screening//////////////////////////////////////////////////////
		//This assumes that the both the harmonic constraint and chemical::CUTPOINT_LOWER/UPPER is setup for the CCD is already set up
 		RNA_LoopCloser rna_loop_closer;

		// Is following really necessary?  I put in the "linear_chainbreak" term in the rna minimizer.
		// I parametrized the the chain_break_screen wrt to to the distance and angular harmonic score. If I remember correctly, screening with the harmonic score was more robust
		// than screening with the linear_chainbreak score. Parin Jan 2, 2009?
		Size const five_prime_chain_break_res = job_parameters_->five_prime_chain_break_res();
		if ( gap_size == 0 )	{
			std::cout << "five_prime_chain_break_res= " << five_prime_chain_break_res << std::endl;
		 	Add_harmonic_chainbreak_constraint(working_pose, five_prime_chain_break_res );
			pose::add_variant_type_to_pose_residue( working_pose, "VIRTUAL_PHOSPHATE", five_prime_chain_break_res+1 );
		}

	  //////////////////////////////////////////Setup Atr_rep_screening/////////////////////////////////////////////////
		Pose base_pose_screen = working_pose;

		//What about the case where the moving_res is internal??? Could use the partition_definition here.?Parin Jan 2, 2009.

		// I think this should work... push apart different parts of the structure so that
		// whatever fa_atr, fa_rep is left is due to "intra-domain" interactions.
		if(Is_dinucleotide){
			base_pose_screen.set_dof( DOF_ID( AtomID( 1, bulge_moving_suite+1 ), core::id::D ), 1.0e4);
		}else{
			base_pose_screen.set_dof( DOF_ID( AtomID( 1, moving_suite+1 ), core::id::D ), 1.0e4);
		}

		base_pose_screen.dump_pdb( "start_extend.pdb" );
		(*atr_rep_screening_scorefxn_)(base_pose_screen);

		EnergyMap const & energy_map=base_pose_screen.energies().total_energies();
		Real base_atr_score = atr_rep_screening_scorefxn_->get_weight(fa_atr) * energy_map[ scoring::fa_atr ];
		Real base_rep_score = atr_rep_screening_scorefxn_->get_weight(fa_rep) * energy_map[ scoring::fa_rep ];
		std::cout << "base_rep= " << base_rep_score << " base_atr= " << base_atr_score << std::endl;



		/////////////////////////////////////////////////////////////////////
		// Need to generalize this -- perhaps return
		//  a vector of TorsionID's... and a
		//  vector of vectors of numbers.
		////////////////////////////////
		utility::vector1< utility::vector1 <utility::vector1 <Real > > > rotamers_groups;

		// Note: Is_prepend should be generalized to include two more possibilities:
		//   we are creating a dinucleotide from scratch --> sample sugar/chi for both moving_res and
		if ( Is_internal  ) {
				PuckerState first_res_puckerstate = Get_residue_pucker_state( pose, moving_res );
				std::cout << "Previous res pucker state is " <<  first_res_puckerstate << std::endl;
				PuckerState second_res_puckerstate = Get_residue_pucker_state( pose, moving_res + 1);
				std::cout << "Second res pucker state is " <<  second_res_puckerstate << std::endl;
				get_rotamers( rotamers_groups, first_res_puckerstate, second_res_puckerstate );
		} else {
			if ( Is_prepend  ) {
				PuckerState second_res_puckerstate = Get_residue_pucker_state( pose, moving_res + 1);
				std::cout << "Second res pucker state is " <<  second_res_puckerstate << std::endl;
				get_rotamers( rotamers_groups, ALL, second_res_puckerstate );
			} else {
				PuckerState first_res_puckerstate = Get_residue_pucker_state( pose, moving_res - 1 );
				std::cout << "Previous res pucker state is " <<  first_res_puckerstate << std::endl;
				get_rotamers( rotamers_groups, first_res_puckerstate , ALL );
			}
		}

		utility::vector1< pose_data_struct2 > pose_data_list;

		//In my old code, chain_break_screening_pose doesn't have a three_prime_chain_break when the gap_size == 0.??
		pose::Pose screening_pose=working_pose; // screening_pose has virtual_phosphate at the three_prime_chain_break even when gap_size == 0
		pose::Pose chain_break_screening_pose=working_pose;
		//Does CCD work if the Virtual phosphate variant type exist? Parin Jan 2, 2009

		Size const total_rotamer_group = rotamers_groups.size();
		Size const rotamer_group_min( 1 );
		Size const rotamer_group_max = total_rotamer_group;


		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		// MAIN LOOP --> rotamer sampling.
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		/////////////////////////////////////////////////////////////////////////////////////////////////////////
		clock_t const time_start( clock() );

//		Size const bulge_rotamer_list_size= (Is_dinucleotide) ? 36 : 1;
			Size const bulge_rotamer_list_size=  (Is_dinucleotide) ? bulge_rotamer_list.size() : 1;
		std::cout << "bulge_rotamer_list.size()= " << bulge_rotamer_list.size() << std::endl;

		std::map<Base_bin , int , compare_base_bin> base_bin_map;
		std::map<Base_bin , int , compare_base_bin>::const_iterator it;

		Real max_centroid_x=0;
		Real max_centroid_y=0;
		Real max_centroid_z=0;
		Real min_centroid_x=0;
		Real min_centroid_y=0;
		Real min_centroid_z=0;
		Real max_centroid_r=0;

		Real max_bin_centroid_x=0;
		Real max_bin_centroid_y=0;
		Real max_bin_centroid_z=0;
		Real min_bin_centroid_x=0;
		Real min_bin_centroid_y=0;
		Real min_bin_centroid_z=0;
		int occupied_bin_count=0;


		Size const spacing=12;

		outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << "name";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "x";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "y";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "z";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "x_bin";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "y_bin";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left << "z_bin";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(2) << std::left << "alpha";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << "euler_z";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << "gamma";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(2) << std::left << "eu_a_bin";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << "eu_z_bin";
		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << "eu_g_bin";
		outfile << "\n";

		std::ofstream count_outfile;
		count_outfile.open("count.txt");
		count_outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << "tot_count";
		count_outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << "occ_bin";
		count_outfile << "\n";

		for(Size bulge_rotamer_ID=1; bulge_rotamer_ID<=bulge_rotamer_list_size ; bulge_rotamer_ID++){
			if(fast_ && count_data_.both_count==100) break;

			std::cout << "bulge_rotamer_ID= " << bulge_rotamer_ID << std::endl;

			if(Is_dinucleotide){
				apply_rotamer( screening_pose, bulge_rotamer_list[bulge_rotamer_ID], bulge_moving_suite);
				apply_rotamer( pose, bulge_rotamer_list[bulge_rotamer_ID], bulge_moving_suite);
			}

		for(Size group_rotamer= rotamer_group_min; group_rotamer<= rotamer_group_max; group_rotamer++) {
			if(fast_ && count_data_.both_count==100) break;

			if( group_rotamer > total_rotamer_group ) {
				std::cout << "group_rotamer > total_rotamer_group(108), break" << std::endl;
				break; //important just for first residue rebuilt
			}

			for(Size subgroup_rotamer=1; subgroup_rotamer<=rotamers_groups[ group_rotamer ].size(); subgroup_rotamer++){
				if(fast_ && count_data_.both_count==100) break;

				std::string const tag = create_tag("U", bulge_rotamer_ID, group_rotamer, subgroup_rotamer, "");

				utility::vector1<Real> const & current_rotamer=rotamers_groups[group_rotamer][subgroup_rotamer]; //If not const & take 3.2 seconds for 10 bulge groups, with const &-> zero time

				apply_rotamer( screening_pose, current_rotamer, moving_suite ); //86.42 seconds for 10 bulges group...not consistent with old results

				count_data_.tot_rotamer_count++;

				if((count_data_.tot_rotamer_count%10000)==0){
					count_outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << count_data_.tot_rotamer_count;
					count_outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << occupied_bin_count;
					count_outfile << "\n";
				}

//				Real blah_num=moving_rsd.xyz(1)[0]+moving_rsd.xyz(1)[1]; //424.25 seconds just to compute this for 10 bulges group...must involve actual torsion rotation here!!
//				Real blah_num_2=moving_rsd.xyz(1)[1]+moving_rsdxyz(1)[2];
//				Real blah_num_3=moving_rsd.xyz(1)[0]+moving_rsd.xyz(1)[2];

				conformation::Residue const & moving_rsd=screening_pose.residue(moving_res);
				Vector const & base_centroid=get_base_centroid( moving_rsd);
				Matrix const & base_coordinate_matrix=get_base_coordinate_system( moving_rsd, base_centroid);

				Vector const & base_centroid_anchor_ribose=Convert_base_centroid_to_anchor_ribose_coord_system( base_centroid, anchor_ribose_stub);
				Matrix const & base_coordinate_matrix_anchor_ribose =	Convert_base_coordinate_matrix_to_anchor_ribose_coord_system(base_coordinate_matrix, anchor_ribose_stub);
				Euler_angles const & euler_angles= Get_euler_angles(base_coordinate_matrix_anchor_ribose);

				Base_bin const & base_bin=Get_base_bin(base_centroid_anchor_ribose, euler_angles);

				Real const base_centroid_radius_anchor_ribose=sqrt(base_centroid_anchor_ribose[0]*base_centroid_anchor_ribose[0]+
																											 		 base_centroid_anchor_ribose[1]*base_centroid_anchor_ribose[1]+
																													 base_centroid_anchor_ribose[2]*base_centroid_anchor_ribose[2]);

				if(max_centroid_x<base_centroid_anchor_ribose[0]) max_centroid_x=base_centroid_anchor_ribose[0];
				if(max_centroid_y<base_centroid_anchor_ribose[1]) max_centroid_y=base_centroid_anchor_ribose[1];
				if(max_centroid_z<base_centroid_anchor_ribose[2]) max_centroid_z=base_centroid_anchor_ribose[2];
				if(min_centroid_x>base_centroid_anchor_ribose[0]) min_centroid_x=base_centroid_anchor_ribose[0];
				if(min_centroid_y>base_centroid_anchor_ribose[1]) min_centroid_y=base_centroid_anchor_ribose[1];
				if(min_centroid_z>base_centroid_anchor_ribose[2]) min_centroid_z=base_centroid_anchor_ribose[2];
				if(max_centroid_r<base_centroid_radius_anchor_ribose) max_centroid_r=base_centroid_radius_anchor_ribose;


				if(max_bin_centroid_x<base_bin.centroid_x) max_bin_centroid_x=base_bin.centroid_x;
				if(max_bin_centroid_y<base_bin.centroid_y) max_bin_centroid_y=base_bin.centroid_y;
				if(max_bin_centroid_z<base_bin.centroid_z) max_bin_centroid_z=base_bin.centroid_z;
				if(min_bin_centroid_x>base_bin.centroid_x) min_bin_centroid_x=base_bin.centroid_x;
				if(min_bin_centroid_y>base_bin.centroid_y) min_bin_centroid_y=base_bin.centroid_y;
				if(min_bin_centroid_z>base_bin.centroid_z) min_bin_centroid_z=base_bin.centroid_z;


				it=base_bin_map.find(base_bin);

				if(it==base_bin_map.end()){
					base_bin_map[base_bin]=1;
					occupied_bin_count++;
				}else{
					base_bin_map[base_bin]=base_bin_map[base_bin]+1;
				}

//				std::cout << "base_centroid x= " << base_centroid[0] << "y= " << base_centroid[1] << "z= " << base_centroid[2];
//				std::cout << " base_centroid_anchor_ribose x= " << base_centroid_anchor_ribose[0] << "y= " << base_centroid_anchor_ribose[1] << "z= " << base_centroid_anchor_ribose[2];
//				std::cout << " euler_angles: alpha: " << euler_angles.alpha << " beta= " <<  << " gamma= " << euler_angles.gamma << std::endl;

/////////////////////////////////////////Mod this part out//////////////////////////////////////////////////////////////////
				if (native_rmsd_screen_ && get_native_pose()){
					// Following does not look right -- what if pose and native_pose are not superimposed corectly? -- rhiju
					if( suite_rmsd(*get_native_pose(), screening_pose, actually_moving_res, Is_prepend)>(native_screen_rmsd_cutoff_)) continue;
					count_data_.rmsd_count++;
					if(verbose_) std::cout << "rmsd_count = " << count_data_.rmsd_count << " total count= " << count_data_.tot_rotamer_count << std::endl;

				}

				/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				if(base_centroid_screening_ && !Is_internal){ //Warning, this does not work if the nucleotide is a bulge
					if ( !Base_centroid_screening( screening_pose, moving_positions, fixed_residues_base_list) ) continue; //484.79 seconds for 10 bulges groups, I think this cause the actual torsion rotation!
				}

				///////////////Chain_break_screening -- special case /////////////////
				if ( gap_size == 1 ){
					// Gap size = 1. This is currently just a distance cut to make sure i-1 is within spitting distance of i+1.
					if ( !Check_chain_closable(screening_pose, five_prime_chain_break_res) ) continue;
					//					if (verbose_ ) std::cout << "Passed gap=1 distance cut" << std::endl;
				}

				//////////////////////////////////////////////////////////////////////////////////////////
				///////////////Chain_break_screening or Full_atom_van_der_Waals_screening/////////////////
				//////////////////////////////////////////////////////////////////////////////////////////
				// This is also where the main pose gets a copy of the torsions.
				//////////////////////////////////////////////////////////////////////////////////////////

				if ( gap_size == 0 ){ //Close chain_break
					//Problem is that chain_break backbone torsions is currently to initialized to value after CCD of last pose July 27, 2009
					//Need to fix this!! But this maybe beneficial since chain break torsion value of the different sampling pose might be similar??

					// Shouldn't we still do a vdw screen for clashes?
					if ( !Full_atom_van_der_Waals_screening( screening_pose, base_rep_score, base_atr_score, false ) ) continue;

					apply_rotamer(chain_break_screening_pose, current_rotamer, moving_suite );
					if( ! Chain_break_screening(screening_pose, chain_break_screening_pose, chainbreak_scorefxn_, rna_loop_closer, base_rep_score)) continue;

					// Need to be very careful here -- do CCD torsions ever overlap with pose torsions?
					//There is no overlap since I took out the delta torsion from RNA_loopcloser... Parin Jan 2, 2009.
					apply_rotamer( pose, current_rotamer, moving_suite );
					Copy_CCD_torsions(pose, chain_break_screening_pose);

				} else {

					if ( !Full_atom_van_der_Waals_screening( screening_pose, base_rep_score, base_atr_score) ) continue;
					apply_rotamer( pose, current_rotamer, moving_suite );

				}

		    outfile << std::setw(20) << std::fixed << std::setprecision(3) << std::left << tag;
		    outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_centroid_anchor_ribose[0];
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_centroid_anchor_ribose[1];
				outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_centroid_anchor_ribose[2];
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_bin.centroid_x;
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_bin.centroid_y;
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_bin.centroid_z;
		  	outfile << std::setw(spacing) << std::fixed << std::setprecision(2) << std::left  << euler_angles.alpha;
		 		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << euler_angles.z;
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << euler_angles.gamma;
		  	outfile << std::setw(spacing) << std::fixed << std::setprecision(2) << std::left  << base_bin.euler_alpha;
		 		outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_bin.euler_z;
			  outfile << std::setw(spacing) << std::fixed << std::setprecision(3) << std::left  << base_bin.euler_gamma;
			  outfile << "\n";

				if(fast_) screening_pose.dump_pdb( tag + ".pdb" );

      	////////////////Add pose to pose_data_list if pose have good score////////////////////////////////////////////

				if ( o2star_screen_ ) sample_o2star_hydrogen( pose );

				Real current_score=Pose_selection_by_full_score(pose_data_list, pose, tag, sampling_scorefxn_);

				if(verbose_){
					std::cout << tag <<  std::endl;
					Output_data(silent_file_data, silent_file_, tag, true, pose, get_native_pose(), actually_moving_res, Is_prepend);
		  	}
////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}//rotamer sub_group for loop
		}//rotamer group for loop

			if(bulge_rotamer_ID==1 || bulge_rotamer_ID==5 || bulge_rotamer_ID==10 || bulge_rotamer_ID==20 || bulge_rotamer_ID==36
				|| bulge_rotamer_ID==48 || bulge_rotamer_ID==60 || bulge_rotamer_ID==72 || bulge_rotamer_ID==84 || bulge_rotamer_ID==96
				|| bulge_rotamer_ID==108 || bulge_rotamer_ID==120 || bulge_rotamer_ID==132 || bulge_rotamer_ID==144){

				std::string const foldername="density_map/" + ObjexxFCL::string_of(bulge_rotamer_ID) + "/";
				system(std::string("rm -r " + foldername).c_str());
				system(std::string("mkdir -p " + foldername).c_str());

//				std::cout << "bulge_rotamer_ID= " << bulge_rotamer_ID << std::endl;
				clock_t const time_start_analyze( clock() );
				Analyze_base_bin_map(base_bin_map, foldername);
				clock_t const time_end_analyze( clock() );

				std::cout << "time in Analyze_base_bin_map: " << static_cast<Real>( time_end_analyze - time_start_analyze ) / CLOCKS_PER_SEC << std::endl;
			}
		}// Bulge rotamer

		Output_title_text("Final sort and clustering");
		std::sort(pose_data_list.begin(), pose_data_list.end(), sort_criteria);
		cluster_pose_data_list(pose_data_list);
		if( pose_data_list.size()>num_pose_kept_ ) pose_data_list.erase(pose_data_list.begin()+num_pose_kept_, pose_data_list.end());
		std::cout<< "after erasing.. pose_data_list= " << pose_data_list.size() << std::endl;

		pose_data_list_=pose_data_list;

		pose = pose_save;
		outfile.close();

		std::cout << "FINAL COUNTS" << std::endl;
		std::cout << "stack= " << count_data_.base_stack_count << " pair= " << count_data_.base_pairing_count;
		std::cout << " atr= " << count_data_.good_atr_rotamer_count;
		std::cout << " rep= " << count_data_.good_rep_rotamer_count;
		std::cout << " both= " << count_data_.both_count;
		std::cout << " rmsd= " << count_data_.rmsd_count << " tot= " << count_data_.tot_rotamer_count << std::endl;

		std::cout << "max_centroid_x= " << max_centroid_x;
		std::cout << " max_centroid_y= " << max_centroid_y;
		std::cout << " max_centroid_z= " << max_centroid_z;
		std::cout << " min_centroid_x= " << min_centroid_x;
		std::cout << " min_centroid_y= " << min_centroid_y;
		std::cout << " min_centroid_z= " << min_centroid_z;
		std::cout << " max_centroid_r= " << max_centroid_r << std::endl;

		std::cout << "max_bin_centroid_x= " << max_bin_centroid_x;
		std::cout << " max_bin_centroid_y= " << max_bin_centroid_y;
		std::cout << " max_bin_centroid_z= " << max_bin_centroid_z;
		std::cout << " min_bin_centroid_x= " << min_bin_centroid_x;
		std::cout << " min_bin_centroid_y= " << min_bin_centroid_y;
		std::cout << " min_bin_centroid_z= " << min_bin_centroid_z << std::endl;

		std::cout << "bulge_rotamer_list_size= " << bulge_rotamer_list_size << std::endl;


		std::cout << "Total time in StepWiseRNA_ResidueSampler: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

	}
*/

}
}
}
