// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNAResidueSampler
/// @brief Not particularly fancy, just minimizes a list of poses.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/StepWiseUtil.hh>
//////////////////////////////////
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
//////////////////////////////////////
#include <numeric/conversions.hh>

#include <iostream>
#include <fstream>
#include <sstream>
#include <ObjexxFCL/format.hh>
#include <set>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

	//Obsolete...to be converted to Rhiju version
 	Size
 	get_matching_atom_name(std::string const & atom_name, conformation::Residue const & rsd){

 		for ( Size atomno=1; atomno<= rsd.nheavyatoms(); ++atomno ) {
 			if(rsd.type().atom_name(atomno)==atom_name){
 				return atomno;
 			}
 		}

 		//Might not find atom_name since atom_name might not actually be a heavy_atom (this happens since any virtual type atom is consider heavy_atom even if if it is a hydrogen atom)
 		return 9999;
 	//		utility_exit_with_message( "Could not find atom_name= " + atom_name + " in rsd");
 	//		exit(1);

 	}

	//Obsolete...to be converted to Rhiju version
 	void
 	setup_suite_atom_id_map(conformation::Residue const & rsd_1,
 													conformation::Residue const & rsd_2,
 													Size const res_num,
 													std::set<std::string> const & special_atom_set,
 													std::string mode,
 													id::AtomID_Map< id::AtomID > & atom_ID_map){

 		if(name_from_aa( rsd_1.aa() )!=name_from_aa( rsd_2.aa())){
 			utility_exit_with_message( "rsd_1.aa()!= rsd_2.aa()");
 		}


 		for ( Size atomno_1=1; atomno_1<= rsd_1.nheavyatoms(); ++atomno_1 ) {

 			std::string const atom_name_1=rsd_1.type().atom_name(atomno_1);

 			if(mode=="exclude_special_atoms"){
 				if(special_atom_set.find(atom_name_1)!=special_atom_set.end()) continue;
 			}else if(mode=="only_special_atoms"){
 				if(special_atom_set.find(atom_name_1)==special_atom_set.end()) continue;
 			}else{
 				utility_exit_with_message( "Invalid mode= " + mode );
 			}

 			Size const atomno_2=get_matching_atom_name(atom_name_1, rsd_2);

 			if(atomno_2==9999) continue;

 			//Check
 			std::string const atom_name_2=rsd_2.type().atom_name(atomno_2);
 			if(atom_name_1!=atom_name_2){
 				utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1= " + atom_name_1 + " atom_name_2= " + atom_name_2);
 			}

 			if(rsd_1.atom_type(atomno_1).name()=="VIRT") continue; //Check for virtual atoms
 			if(rsd_2.atom_type(atomno_2).name()=="VIRT") continue; //Check for virtual atoms

 			id::AtomID const id1( atomno_1, res_num);
 			id::AtomID const id2( atomno_2, res_num);
 			atom_ID_map.set( id1, id2 );

 		}

 	}

 	//Virtual types mess up numbering...This will make sure that numbering is correct.
 	void
 	setup_suite_atom_id_map(pose::Pose const & pose_1, pose::Pose const & pose_2, Size const base_res, bool Is_prepend, id::AtomID_Map< id::AtomID > & atom_ID_map){


 		std::set<std::string> special_atom_set;
 		std::set<std::string>::iterator it;

 		if(Is_prepend){
 			special_atom_set.insert(std::string(" P  "));
 			special_atom_set.insert(std::string(" O1P"));
 			special_atom_set.insert(std::string(" O2P"));
 			special_atom_set.insert(std::string(" O5*"));

 		}



 		conformation::Residue const & base_rsd_1=pose_1.residue(base_res);
 		conformation::Residue const & base_rsd_2=pose_2.residue(base_res);

 		setup_suite_atom_id_map(base_rsd_1, base_rsd_2, base_res, special_atom_set, "exclude_special_atoms", atom_ID_map);

 		if(Is_prepend){

 			Size const upper_suite_res=base_res+1;

 			conformation::Residue const & base_rsd_1=pose_1.residue(upper_suite_res);
 			conformation::Residue const & base_rsd_2=pose_2.residue(upper_suite_res);

 			setup_suite_atom_id_map(base_rsd_1, base_rsd_2, upper_suite_res, special_atom_set, "only_special_atoms", atom_ID_map);

 		}
 	}

	void
	output_pair_size(std::pair<Size, Size> const & pair_size){
		std::cout << "( " << pair_size.first << ", " << pair_size.second << " ) ";
	}

	void
	output_pair_size_vector(utility::vector1 <std::pair<Size, Size> > const & pair_size_vector, std::string const & output_string){
		std::cout << std::setw(60) << std::left << output_string << " :";
		for(Size n=1; n<=pair_size_vector.size(); n++){
			output_pair_size(pair_size_vector[n]);
		}
		std::cout << std::endl;
	}

	bool
	pair_sort_citeria(std::pair<Size, Size> pair_one, std::pair<Size, Size> pair_two){
		return (pair_one.first < pair_two.second);
	}

	void
	Sort_pair_list(utility::vector1< std::pair<Size, Size> > pair_list){
		sort(pair_list.begin(), pair_list.end(), pair_sort_citeria);
	}


	bool
	seq_num_sort_citeria(core::Size seq_num_1, core::Size seq_num_2){
		return (seq_num_1 < seq_num_2);
	}

	void
	sort_seq_num_list(utility::vector1<core::Size> & seq_num_list) {
		sort(seq_num_list.begin(), seq_num_list.end(), seq_num_sort_citeria);
	}

	void
	Output_seq_num_list(std::string const tag, utility::vector1<core::Size> const & seq_num_list, core::Size const spacing){

		using namespace ObjexxFCL;
		using namespace ObjexxFCL::fmt;

		std::cout <<  std::setw(spacing) << tag;

		utility::vector1<core::Size> sorted_seq_num_list=seq_num_list;
		sort_seq_num_list(sorted_seq_num_list);

		Size seq_num=1;
		for(Size n=1; n<=sorted_seq_num_list.size(); n++){

			while(seq_num<sorted_seq_num_list[n]){
				std::cout << A(4," ");
				seq_num++;
			}
			std::cout << I(4, sorted_seq_num_list[n]);
			seq_num++;
		}

		std::cout << std::endl;

	}


	Size
	Get_five_prime_chain_break(pose::Pose const & pose){

		if(!Is_close_chain_break(pose)){
			std::cout << "In StepWiseRNA_ResidueSampler::Get_five_prime_chain_break, the function should not be called! "<< std::endl;
			exit(1);
		};

		for(Size seq_num = 1; seq_num <= pose.total_residue(); seq_num++){

			if ( !pose.residue( seq_num   ).has_variant_type( chemical::CUTPOINT_LOWER )  ) continue;
			if ( !pose.residue( seq_num+1 ).has_variant_type( chemical::CUTPOINT_UPPER )  ) continue;

			return seq_num;
		}

		std::cout << "In StepWiseRNA_ResidueSampler::Get_five_prime_chain_break, cannot find five_prime_chain_break" << std::endl;
		exit(1);
	}

	void
	Add_harmonic_chainbreak_constraint(pose::Pose & pose, Size const five_prime_res){

		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::id;

		using numeric::conversions::radians;

		Size three_prime_res=five_prime_res+1;


		//    From RAD.param file
		//	  ICOOR_INTERNAL  UPPER -175.907669   60.206192    1.607146   O3*   C3*   C4*   , Upper is P1
		//    ICOOR_INTERNAL  LOWER  -64.027359   71.027062    1.593103   P     O5*   C5*   , Lower is O3'
		//    Bug that bond distance is not the same. Rhiju suggest using 1.593103

		//Original (Amber?) parameter 1.608, 119.8, 103.4

		Real const O3_P_distance( 1.593 ); //amber=1.608
		Real const O3_angle( 119.8 ); // 180-60.206192
		Real const  P_angle( 108.97 ); // Quite off from original (original==Amber's parameter??) ,180-71.027062

		Real const distance_stddev( 0.0659 ); // amber is 0.0659
		Real const angle_stddev_degrees_P( 8.54 ); // amber is 8.54 (P angle), 5.73 (O3 angle)
		Real const angle_stddev_degrees_O3( 5.73 );

		ConstraintSetOP cst_set( pose.constraint_set()->clone() );
		assert( cst_set ); //if ( !cst_set ) cst_set = new ConstraintSet();


		FuncOP const distance_func( new HarmonicFunc( O3_P_distance, distance_stddev ) );
		FuncOP const O3_angle_func( new HarmonicFunc( radians( O3_angle ), radians( angle_stddev_degrees_P ) ) );
		FuncOP const  P_angle_func( new HarmonicFunc( radians(  P_angle ), radians( angle_stddev_degrees_O3 ) ) );



		Residue const & rsd1( pose.residue(five_prime_res) );
		Residue const & rsd2( pose.residue(three_prime_res) );

		AtomID const C3_id( rsd1.atom_index( "C3*" ), five_prime_res);
		AtomID const O3_id( rsd1.atom_index( "O3*" ), five_prime_res);
		AtomID const  P_id( rsd2.atom_index( "P"   ), three_prime_res);
		AtomID const O5_id( rsd2.atom_index( "O5*" ), three_prime_res);

		// distance from O3* to P
		cst_set->add_constraint( new AtomPairConstraint( O3_id, P_id, distance_func ) );

		// angle at O3*
		cst_set->add_constraint( new AngleConstraint( C3_id, O3_id, P_id, O3_angle_func ) );

		// angle at P
		cst_set->add_constraint( new AngleConstraint( O3_id, P_id, O5_id,  P_angle_func ) );

		pose.constraint_set( cst_set );

		//std::cout << "pose_cst_set" << std::endl;
		//cst_set->show(std::cout);

	}



	void
	Output_title_text(std::string const title){

		std::cout << std::endl;

		Size title_length=title.size();
		Size char_per_line=80;
		Size dash_length=char_per_line-title_length;

		for(Size i=1; i<=dash_length/2; i++){
			std::cout << "-";
		}

		std::cout << title;

		for(Size i=1; i<=dash_length/2; i++){
			std::cout << "-";
		}

		std::cout << std::endl;

	}


	void
	Output_fold_tree_info(kinematics::FoldTree const & fold_tree, std::string const pose_name){

		std::cout << "fold tree of " << pose_name << ": " << std::endl;
		for(int i=1; i<=fold_tree.num_cutpoint(); i++){
			std::cout << std::setw(30) << "jump_point_num= " << i;
			std::cout << "   cutpoint= " << fold_tree.cutpoint(i);
			std::cout << "   5' jump_point= " << fold_tree.jump_point( 1, i ) << "," << fold_tree.upstream_atom(i);
			std::cout << "   3' jump_point= " << fold_tree.jump_point( 2, i ) << "," << fold_tree.downstream_atom(i) << std::endl;
		}
	}

	/////////////////////////////////////////////////////////////
	void
	Output_fold_tree_info(pose::Pose const & pose, std::string pose_name){
		Output_fold_tree_info(pose.fold_tree(), pose_name);
	}

	/////////////////////////////////////////////////////////////
	void
	Output_pose_data_list(utility::vector1 <pose_data_struct2> const & pose_data_list, std::string const & silent_file , bool const write_score_only){
		using namespace core::io::silent;

		core::io::silent::SilentFileData silent_file_data;

		for ( Size n = 1; n <= pose_data_list.size(); n++ ) {
			BinaryRNASilentStruct s( *(pose_data_list[n].pose_OP), pose_data_list[n].tag );
			silent_file_data.write_silent_struct( s, silent_file, write_score_only ) ;
		}

	}


	/////////////////////////////////////////////////////////////
	void
	Output_data(core::io::silent::SilentFileData& silent_file_data, std::string const & silent_file, std::string const & tag, bool const write_score_only, pose::Pose const & pose, core::pose::PoseCOP native_poseCOP, Size const moving_base_residue, bool const Is_prepend){

		using namespace core::io::silent;
		using namespace core::scoring;
		using namespace core::pose;

		BinaryRNASilentStruct s( pose, tag ); //Does this take a long time to create?
		if ( native_poseCOP ) {
			s.add_energy( "all_rms", rms_at_corresponding_heavy_atoms( pose, *native_poseCOP ) );
			// Following does not look right -- what if pose and native_pose are not superimposed correctly? -- rhiju.
			// Yes, this assume that pose and native_pose are correctly syperimposed.
			// I added a function in Pose_Setup to make sure this happens. Parin Jan 28, 2010

			s.add_energy( "rmsd", suite_rmsd( pose, *native_poseCOP, moving_base_residue, Is_prepend));
		}

		if ( silent_file.size() > 0 ) silent_file_data.write_silent_struct(s, silent_file, write_score_only);
		silent_file_data.add_structure( s );
	}


	void
	output_rotamer(utility::vector1 <Real > & rotamer){

		Size const spacing=10;

		std::cout <<  std::setw(18) << "Torsions= ";
		std::cout << std::setw(spacing)<< rotamer[1] << " ";
		std::cout << std::setw(spacing)<< rotamer[2] << " ";
		std::cout << std::setw(spacing)<< rotamer[3] << " ";
		std::cout << std::setw(spacing)<< rotamer[4] << " ";
		std::cout << std::setw(spacing)<< rotamer[5] << " ";
		std::cout << std::setw(spacing)<< rotamer[6] << " ";
		std::cout << std::setw(spacing)<< rotamer[7] << " ";
		std::cout << std::setw(spacing)<< rotamer[8] << " ";
		std::cout << std::setw(spacing)<< rotamer[9] << " ";
		std::cout << std::setw(spacing)<< rotamer[10] << " ";
		std::cout << std::setw(spacing)<< rotamer[11] << " ";
		std::cout << std::setw(spacing)<< rotamer[12] << " ";
		std::cout << std::setw(spacing)<< rotamer[13] << " ";
		std::cout << std::endl;

	}

	///////////////////////////////////////////////////////////////////////////
	// Following needs to move to StepWiseRNA_RotamerGenerator.cc!!!
	///////////////////////////////////////////////////////////////////////////
/*
		void
	get_bulge_rotamers( utility::vector1< utility::vector1 <Real> >& rotamer_list, PuckerState const & pucker1, PuckerState const & pucker2, Size const set_num=1) {

		std::cout << "set_num= " << set_num << std::endl;

		scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;

		Real bin_size=20;

		utility::vector1< Real > torsion_samples;
		torsion_samples.push_back( 0.0 );
	//	torsion_samples.push_back( -1.0 );
	//	torsion_samples.push_back( 1.0 );

		if(pucker1==ALL && pucker2==ALL){
			std::cout << "Cannot sample both pucker" << std::endl;
			exit(1);
		}

		Real delta, nu2, nu1, chi;

		Real alpha2, beta2, gamma2, epsilon1, zeta1;
		Real delta1, nu2_1, nu1_1, chi_1;
		Real delta2, nu2_2, nu1_2, chi_2;

		scoring::rna::Gaussian_parameter gp( 0.0, 0.0, 0.0);

//		output_rotamer_title();

		for (int d = 1; d <= 2; d++ ) {

			if (d == 1) {
				delta = rna_fitted_torsion_info.gaussian_parameter_set_delta_north()[1].center;
				chi = rna_fitted_torsion_info.gaussian_parameter_set_chi_north()[1].center;
				nu2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_north()[1].center;
				nu1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_north()[1].center;

			}	else {
				delta = rna_fitted_torsion_info.gaussian_parameter_set_delta_south()[1].center;
				chi = rna_fitted_torsion_info.gaussian_parameter_set_chi_south()[1].center;
				nu2 = rna_fitted_torsion_info.gaussian_parameter_set_nu2_south()[1].center;
			nu1 = rna_fitted_torsion_info.gaussian_parameter_set_nu1_south()[1].center;
			}

			if(pucker1==ALL){
				delta1=delta, nu2_1=nu2, nu1_1=nu1, chi_1=chi;
				delta2=0.0, nu2_2=0.0, nu1_2=0.0, chi_2=0.0;
			}else{
				delta1=0.0, nu2_1=0.0, nu1_1=0.0, chi_1=0.0;
				delta2=delta, nu2_2=nu2, nu1_2=nu1, chi_2=chi;
			}

				//epsilon depends on delta of same residue...
				for (Size e1 = 1; e1 <= 1; e1++ ) {
					for (Size a2 = 1; a2 <= 3; a2++ ) {
					  //zeta depends on alpha of next residue...
						for (Size z1 = 1; z1 <= 2; z1++ ) {
							for (Size b2 = 1; b2 <= 1; b2++ ) {
								for (Size g2 = 1; g2 <= 3; g2++ ) {





	//More Rotamers
	/////////////////////////////////////////////////

								for ( Size e1_std = 1; e1_std <= torsion_samples.size(); e1_std++ ) {

									if(pucker1==0){//Sample sugar pucker 1
										if (d == 1) 		 epsilon1 = -150.17 +	torsion_samples[ e1_std]*bin_size; //North , Rhiju's mean
										else         		 epsilon1 = -98.45 +	torsion_samples[ e1_std]*bin_size; //South , Rhiju's mean
									}else{
										if (pucker1==1 ) epsilon1 = -150.17 +	torsion_samples[ e1_std]*bin_size; //North , Rhiju's mean
										else         		 epsilon1 = -98.45 +	torsion_samples[ e1_std]*bin_size; //South , Rhiju's mean
									}


				for ( Size a2_std = 1; a2_std <= torsion_samples.size(); a2_std++ ) {
						alpha2 = rna_fitted_torsion_info.gaussian_parameter_set_alpha()[a2].center + torsion_samples[a2_std]*bin_size;

						if(set_num==2 || set_num==4) {
							alpha2=alpha2+60;
							if(alpha2>180) alpha2=alpha2-360;
						}

						for ( Size z1_std = 1; z1_std <= torsion_samples.size(); z1_std++ ) {
								if (a2 == 1)    zeta1 = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_minus()[z1].center + torsion_samples[z1_std]*bin_size;
								else if (a2==2) zeta1 = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_sc_plus()[z1].center + torsion_samples[z1_std]*bin_size;
								else            zeta1 = rna_fitted_torsion_info.gaussian_parameter_set_zeta_alpha_ap()[z1].center + torsion_samples[z1_std]*bin_size;

								for ( Size b2_std = 1; b2_std <= torsion_samples.size(); b2_std++ ) {

										beta2 = rna_fitted_torsion_info.gaussian_parameter_set_beta()[b2].center + torsion_samples[b2_std]*bin_size;

										for ( Size g2_std = 1; g2_std <= torsion_samples.size(); g2_std++ ) {
												gamma2 = rna_fitted_torsion_info.gaussian_parameter_set_gamma()[g2].center + torsion_samples[g2_std]*bin_size;

												if(set_num==3 || set_num==4) {
													gamma2=gamma2+60;
													if(gamma2>180) gamma2=gamma2-360;
												}

													utility::vector1 < Real >  backbone_rotamer; //Would be better to make this a class, RNA_SuiteToSuiteRotamer

													backbone_rotamer.push_back( delta1 );    //1
													backbone_rotamer.push_back( chi_1 );     //2
													backbone_rotamer.push_back( nu2_1 );     //3
													backbone_rotamer.push_back( nu1_1 );     //4
													backbone_rotamer.push_back( epsilon1 );  //5
													backbone_rotamer.push_back( zeta1 );     //6
													backbone_rotamer.push_back( alpha2 );    //7
													backbone_rotamer.push_back( beta2 );     //8
													backbone_rotamer.push_back( gamma2 );    //9
													backbone_rotamer.push_back( delta2 );    //10
													backbone_rotamer.push_back( chi_2 );     //11
													backbone_rotamer.push_back( nu2_2 );     //12
													backbone_rotamer.push_back( nu1_2 );     //13


													output_rotamer(backbone_rotamer);
													rotamer_list.push_back( backbone_rotamer );

											  } // gamma2_samples
										  } // beta2_samples
									  } // zeta1_samples
								  } //  alpha2_samples
							  } // epsilon1_samples

								} // gamma2
							} // beta2
						} // zeta1
					} // alpha2
				} // epsilon1
			} // delta2
	}
*/


	void
	Add_virtual_O2Star_hydrogen( core::pose::Pose & pose){

	  for (core::Size i = 1; i <= pose.total_residue(); i++){
	    pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", i);
	  }
	}

	bool
	Remove_virtual_O2Star_hydrogen(pose::Pose & pose){
		bool did_something( false );
		for(Size i=1; i<=pose.total_residue(); i++){
			if ( pose.residue_type( i ).has_variant_type( "VIRTUAL_O2STAR_HYDROGEN" ) ){
				pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", i);
				did_something = true;
			}
		}
		return true;
	}



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//This now works for any rebuild residue.
	Real
	suite_rmsd(pose::Pose const & pose1, pose::Pose const & pose2, Size const & moving_res_num, bool const prepend_res){

		Size atom_count=0;
  	Real sum_sd=0;

		suite_square_deviation(pose1, pose2, prepend_res, moving_res_num, atom_count, sum_sd, false /*verbose*/);

		sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);

		return (std::max(0.01, rmsd));

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	rmsd_over_residue_list(pose::Pose const & pose1, pose::Pose const & pose2, utility::vector1 <Residue_info> const & residue_list, std::map< core::Size, core::Size > const & res_map, std::map< core::Size, bool > const & Is_prepend_map, bool const verbose) {


		if(verbose){
			std::cout << "In residue_list_rmsd function" << std::endl;;
			std::cout <<  std::setw(40)  << "rebuild_residue_list: "; Output_residue_list(residue_list);
		}

  	Size atom_count=0;
  	Real sum_sd=0;

  	for(Size i=1; i<=residue_list.size(); i++){

  		Size const full_seq_num= residue_list[i].seq_num;
 			Size const seq_num=res_map.find(full_seq_num)->second;
			bool Is_prepend=Is_prepend_map.find(full_seq_num)->second;

			if(verbose) std::cout << "Full_seq_num= " << full_seq_num << " partial_seq_num= " << seq_num << std::endl;

			//add atom in the suites to atom_count
			//add sd of each atom to sum_sd
			suite_square_deviation(pose1, pose2, Is_prepend, seq_num, atom_count, sum_sd, verbose);

		}

  	sum_sd=sum_sd/(atom_count);
  	Real rmsd=sqrt(sum_sd);


		if(verbose) std::cout << "sum_sd= " << sum_sd << " atom_count= " << atom_count << " rmsd= " << rmsd << std::endl;

  	return (std::max(0.01, rmsd));
}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Print_heavy_atoms(Size suite_num, pose::Pose const & pose1, pose::Pose const & pose2){

		using namespace conformation;
		Size num_atoms;

		num_atoms=std::max(pose1.residue(suite_num).nheavyatoms(), pose2.residue(suite_num).nheavyatoms());

		std::cout << "num_atoms: " << num_atoms << std::endl;

		for(Size n=1; n<= num_atoms;  n++){

			std::cout << " atom num = " <<  n;
			std::cout << "  atom_name of the pose1 " <<  pose1.residue(suite_num).atom_name(n);
			std::cout << "  atom_name of the pose2 " <<  pose2.residue(suite_num).atom_name(n) << std::endl;

		}
	}
//

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Size
	Get_num_side_chain_atom_from_res_name(chemical::AA const & res_aa, bool const verbose){

		Size num_side_chain_atom;

		if(name_from_aa(res_aa)=="RAD") {
			if(verbose) std::cout << "name_from_aa: RAD" << std::endl;
			num_side_chain_atom=11;
		} else if(name_from_aa(res_aa)=="RCY") {
			if(verbose) std::cout << "name_from_aa: RCY" << std::endl;
			num_side_chain_atom=9;
		} else if(name_from_aa(res_aa)=="RGU") {
			if(verbose) std::cout << "name_from_aa: RGU" << std::endl;
			num_side_chain_atom=12;
		} else if(name_from_aa(res_aa)=="URA") {
			if(verbose) std::cout << "name_from_aa: URA" << std::endl;
			num_side_chain_atom=9;
		} else {
			std::cout << "Error, cannot identify residue type" << std::endl;
			num_side_chain_atom=0;
			exit (1);
		}

		return num_side_chain_atom;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Real
	atom_square_deviation(conformation::Residue const & rsd_1, conformation::Residue const & rsd_2, Size const & atomno_1, Size const & atomno_2, bool verbose){



		//////////////Might take this out later if it significantly slow down the code///////////////
		std::string const & atom_name_1=rsd_1.type().atom_name(atomno_1);
		std::string const & atom_name_2=rsd_2.type().atom_name(atomno_2);

		if(atom_name_1!=atom_name_2){
			utility_exit_with_message( "atom_name_1 != atom_name_2, atom_name_1= " + atom_name_1 + " atom_name_2= " + atom_name_2);
		}

//This can be turned on for debugging, but by default is mod out since string comparison might be slow Parin Jan 28, 2009
//		if(rsd_1.atom_type(atomno_1).name()=="VIRT"  || rsd_2.atom_type(atomno_2).name()=="VIRT") {
//			utility_exit_with_message( "rsd_1.atom_type(n).name()==\"VIRT\"  || rsd_2.atom_type(n).name()==\"VIRT\" =TRUE!");
//		}
		////////////////////////////////////////////////////////////////////////////////////////////////

  	Distance const dist_squared = (rsd_1.xyz( atomno_1) - rsd_2.xyz( atomno_2)).length_squared();


		if(verbose){
			std::cout << " atom_name of the atom1= " << atom_name_1 << " " << rsd_1.seqpos();
			std::cout << " atom_name of the atom2= " << atom_name_2 << " " << rsd_2.seqpos();
			std::cout << " Dist_squared= " << dist_squared << std::endl;
		}

		return dist_squared;

	}

//	Could get the atom_num using the atom_name but I don't think this will be fast Size atom_P = rsd.atom_index( " P  " ); Jan 28, 2009 Parin S.

	void
	suite_square_deviation(pose::Pose const & pose1, pose::Pose const & pose2,

												 bool const & prepend_res, Size const & moving_res_num,
												 Size& atom_count, Real& sum_sd, bool verbose){

  		chemical::AA const & res_aa =  pose1.residue(moving_res_num).aa();
  		chemical::AA const & res_aa2 =  pose2.residue(moving_res_num).aa();

  		Size first_sidechain_atom1=pose1.residue(moving_res_num).first_sidechain_atom();
  		Size first_sidechain_atom2=pose2.residue(moving_res_num).first_sidechain_atom();

			Size num_side_chain_atom=Get_num_side_chain_atom_from_res_name(res_aa, verbose);

 			if(false && verbose){
  			std::cout << " residue type1= " <<  res_aa << " residue type2= " <<  res_aa2;
  			std::cout << " 1st_side_atom1= " <<  first_sidechain_atom1;
  			std::cout << " 1st_side_atom2= " <<  first_sidechain_atom2;
  			std::cout << " nheavyatoms1= " <<  pose1.residue(moving_res_num).nheavyatoms();
  			std::cout << " nheavyatoms2= " <<  pose2.residue(moving_res_num).nheavyatoms() << std::endl;
				std::cout << " num_side_chain_atom = " <<  num_side_chain_atom << std::endl;
				Print_heavy_atoms(moving_res_num,pose1, pose2);
			}

			if ( verbose ) {
				std::cout << " MOVING_RES: " << moving_res_num << "    PREPEND? " << prepend_res << std::endl;
			}

			Size const num_heavy_backbone_atoms=11; //RNA contain 11 heavy backbone atoms.

			for( Size atomno = 1; atomno <= num_heavy_backbone_atoms; atomno++){

				//atom 1-4 are " P  "), " O1P", " O2P" and " O5*"
				Size const res_num= (prepend_res && atomno <= 4) ? moving_res_num + 1: moving_res_num;

				atom_count++;

				conformation::Residue const & rsd_1=pose1.residue(res_num);
				conformation::Residue const & rsd_2=pose2.residue(res_num);

				sum_sd=sum_sd+atom_square_deviation(rsd_1, rsd_2, atomno, atomno, verbose);
  		}

  		//Need to use num_side_chain_atom from pose1 since a silly bug in Rosetta miscalculate num_heavy_atom by considering
			//the virtaul O2star hydrogen to be heavy_atom when it is set to virtual in the current_pose_screen
  		for( Size n = 1; n<= num_side_chain_atom; n++){ //INCLUDE the O2star oxygen
  			atom_count++;

				Size const atomno_1=(n-1)+first_sidechain_atom1;
				Size const atomno_2=(n-1)+first_sidechain_atom2;

				conformation::Residue const & rsd_1=pose1.residue(moving_res_num);
				conformation::Residue const & rsd_2=pose2.residue(moving_res_num);

				sum_sd=sum_sd+atom_square_deviation(rsd_1, rsd_2, atomno_1, atomno_2, verbose);
 		 	}


			// O1P<-->O2P check on phosphate positions closest to fixed side
			if(verbose && false){
				Size const res_num= (prepend_res) ? moving_res_num + 1: moving_res_num;
				Distance dist_squared = (pose1.residue(res_num).xyz( 2 ) - pose2.residue(res_num).xyz( 3 ) ).length_squared();
				std::cout << " atom_name of the pose1= " << pose1.residue(res_num).atom_name(2)  << " " << res_num;
  			std::cout << " atom_name of the pose2= " << pose2.residue(res_num).atom_name(3)  << " " << res_num;
				std::cout << " Switch Phosphate1= " << " dist_squared= " << dist_squared << std::endl;
				dist_squared = (pose1.residue(res_num).xyz( 3 ) - pose2.residue(res_num).xyz( 2 ) ).length_squared();
				std::cout << "atom_name of the pose1= " << pose1.residue(res_num).atom_name(3)  << " " << res_num;
  			std::cout << " atom_name of the pose2= " << pose2.residue(res_num).atom_name(2)  << " " << res_num;
				std::cout << " Switch Phosphate2= " << " dist_squared= " << dist_squared << std::endl;
			}


	}

	/////////////////////////////////////////////////////////////////////////////////////////////////
	//This function is called when the chain_break gap_size = 0 or 1 Parin S, Jan 28, 2010.
	bool
	Check_chain_closable(pose::Pose const & pose, Size const five_prime_chain_break_res, Size const gap_size ){

		if ( gap_size > 1 ) return true;

		//Use to call Calculate_dist_to_close_chain function, but that requires rebuild_unit_struct
		Size three_prime_chain_break_res= five_prime_chain_break_res + 1;

		Distance dist_to_close_chain=(pose.residue(three_prime_chain_break_res).xyz("C5*") - pose.residue(five_prime_chain_break_res).xyz("O3*") ).length();

		if ( gap_size == 1 ) {
			static Distance const cutoff_distance = 11.0138;
			return ( dist_to_close_chain < cutoff_distance );
		}

		assert( gap_size == 0 );

		static Distance const cutoff_distance_min( 2.0 );
		static Distance const cutoff_distance_max( 4.627 );

		//basically cannot close chain if the C5_O3_distance is either too short or too long.
		if( (dist_to_close_chain > cutoff_distance_max) ||
				dist_to_close_chain < cutoff_distance_min ) return false;

		return true;

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	Freeze_sugar_torsions(core::kinematics::MoveMap& mm, Size const total_residue){

	 using namespace core::id;

	 std::cout << "Freeze pose sugar torsions, total_residue=" << total_residue << std::endl;

	 for(Size i=1; i<=total_residue; i++){

			mm.set( TorsionID( i  , id::BB,  4 ), false ); //delta_i
			mm.set( TorsionID( i  , id::CHI, 2 ), false ); //nu2_i
			mm.set( TorsionID( i  , id::CHI, 3 ), false );	//nu1_i

	 }
	}


	void
	o2star_minimize(pose::Pose& pose, core::scoring::ScoreFunctionOP const & packer_scorefxn){

		//TR << "Repacking 2'-OH ... " << std::endl;

		pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
		task->initialize_from_command_line();

		for (Size i = 1; i <= pose.total_residue(); i++) {
			if ( !pose.residue(i).is_RNA() ) continue;
			task->nonconst_residue_task(i).and_extrachi_cutoff( 0 );
			task->nonconst_residue_task(i).or_ex4( true ); // might as well sample a bit extra.
			task->nonconst_residue_task(i).or_include_current( true );
			// How about bump check?
		}

		//	TR << "Orienting 2' hydroxyls..." << std::endl;
		//		pack::pack_rotamers( pose, *packer_scorefxn, task);

		pack::rotamer_trials( pose, *packer_scorefxn, task);
	}

	utility::vector1<std::string>
	Tokenize_with_input_delimiters(std::string const & str, std::string const delimiters)
	{
	  using namespace std;

	  utility::vector1<std::string> tokens;

    // Skip delimiters at beginning.
    string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    string::size_type pos     = str.find_first_of(delimiters, lastPos);

    while (string::npos != pos || string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }

		return tokens;
	}

	/////////////////////////////////////////////////////////////////////////////////////
	//When a CUTPOINT_UPPER is added to 3' chain_break residue, the EXISTENCE of the CUTPOINT_UPPER atoms means that the alpha torsion which previously DOES NOT exist due to the chain_break now exist. The alpha value is automatically defined to the A-form value by Rosetta. However Rosetta does not automatically adjust the O1P and O2P atom position to account for this fact. So it is important that the O1P and O2P atoms position are correctly set to be consistent with A-form alpha torsion before the CUTPOINT_UPPER IS ADDED Parin Jan 2, 2009
	void
	Correctly_position_cutpoint_phosphate_torsions(pose::Pose & current_pose, Size five_prime_chainbreak,  bool verbose /*=false*/){

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::id;
		using namespace core::io::pdb;

		static const ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		if(verbose) Output_title_text("Correctly_position_cutpoint_phosphate_torsions function");

		chemical::AA res_aa = aa_from_name( "RAD" );
		ResidueOP new_rsd = conformation::ResidueFactory::create_residue( *(rsd_set->aa_map( res_aa )[1]) ) ;
		if(verbose) std::cout << "  res_aa: " << res_aa << std::endl;

		Size three_prime_chainbreak=five_prime_chainbreak+1;

		if(verbose){
			dump_pdb(current_pose, "Before_prepending_dummy_nucleotide");
			std::cout <<  std::setw(50) << "Before_prepending_dummy_nucleotide";
			//print_backbone_torsions(current_pose, five_prime_chainbreak);
		}

		current_pose.prepend_polymer_residue_before_seqpos( *new_rsd, three_prime_chainbreak, true);
		scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;


		if(verbose){
			dump_pdb(current_pose, "Before_setting_torsion_to_A_form.pdb");
			std::cout << std::setw(50) << "Before_setting_torsion_to_A_form";
			//print_backbone_torsions(current_pose, five_prime_chainbreak+1);
		}
		//These are the initial value of virtual upper and lower cutpoint atom. Actaully only the alpha (id::BB, 1) is important here since it is need to correctly orient O1P and O2P atom
		//Where the hell did I get these numbers from value...by appending with ideal geometry and look at the initalized value? Oct 13, 2009
		current_pose.set_torsion( TorsionID( five_prime_chainbreak+1, id::BB, 5 ), -151.943 );
		current_pose.set_torsion( TorsionID( five_prime_chainbreak+1, id::BB, 6 ), -76.4185 );
		current_pose.set_torsion( TorsionID( three_prime_chainbreak+1, id::BB, 1 ), -64.0274 );

		if(verbose){
			dump_pdb(current_pose, "After_setting_torsion_to_A_form.pdb");
			std::cout << std::setw(50) << "After_setting_torsion_to_A_form"; //print_backbone_torsions(current_pose, five_prime_chainbreak+1);
		}

		current_pose.delete_polymer_residue(five_prime_chainbreak+1);

		if(verbose){
			dump_pdb(current_pose, "After_deleting_dummy_nucleotide");
			std::cout << std::setw(50) << "After_deleting_dummy_nucleotide"; //print_backbone_torsions(current_pose, five_prime_chainbreak);
		}
	}

	////////////////////////////////////////////////////////////////////////
	void
	apply_rotamer( pose::Pose & pose,
								 utility::vector1< Torsion_Info > const & rotamer_list)
	{
		for ( Size i = 1; i <= rotamer_list.size(); i++ ) {
			pose.set_torsion( rotamer_list[ i ].id , rotamer_list[i].value);
		}
	}


	//This is the version currently used for floating base sampling....slightly different from the old version of Base_centroid_screening which appear to StepWiseRNA_Sampling
	//Should probably integrate this with Rhiju's class Jan 28, 2010. ***ALERT***RHIJU pointed out that the screening condition is slight different in his new class.

	bool
	Base_centroid_screening( core::pose::Pose const & pose,
													 utility::vector1< Size > moving_positions,
													 utility::vector1 < base_stub > const & other_residues_base_list,
													 SillyCountStruct & count_data){

		// static to prevent any overhead of remaking this base.
		static base_stub moving_res_base;

		bool stack_base( false ), base_pairing( false );

		for ( Size m = 1; m <= moving_positions.size(); m++ ) {

			core::conformation::Residue const & residue_object( pose.residue( moving_positions[ m ] ) );

			Size const moving_pos( moving_positions[ m ] );
			moving_res_base.centroid=get_base_centroid( residue_object );
			moving_res_base.base_coordinate_matrix=get_base_coordinate_system( residue_object, moving_res_base.centroid );

			if(Base_centroid_screening(moving_res_base, other_residues_base_list, count_data)==true) return true;

		}

		return false;
	}

	bool
	Base_centroid_screening(base_stub const & moving_res_base, utility::vector1 < base_stub > const & other_residues_base_list, SillyCountStruct & count_data){

		bool stack_base = false;
		for(Size i=1; i<=other_residues_base_list.size(); i++){
			base_stub const & base_info=other_residues_base_list[i];
			numeric::xyzVector<Real> const other_z_vector=base_info.base_coordinate_matrix.col_z();
			numeric::xyzVector<Real> rebuild_z_vector=moving_res_base.base_coordinate_matrix.col_z();

			numeric::xyzVector<Real> centroid_diff;
			subtract( moving_res_base.centroid, base_info.centroid, centroid_diff);
			Real centroid_distance=centroid_diff.length();
			if(centroid_distance>6.364) continue;

			Real base_z_offset_one=std::abs(dot( centroid_diff, other_z_vector));
			Real base_z_offset_two=std::abs(dot( centroid_diff, rebuild_z_vector));

			if((base_z_offset_one>4.5 || base_z_offset_one<2.5) && (base_z_offset_two>4.5 || base_z_offset_two<2.5)) continue;

			Real base_axis_one=base_z_offset_one/centroid_distance;
			Real base_axis_two=base_z_offset_two/centroid_distance;

			if(base_axis_one<0.707 && base_axis_two<0.707) continue;

			Real base_planarity=std::abs(dot( other_z_vector, rebuild_z_vector));

			if(base_planarity<0.707) continue;

			//If reach this point means success!
			stack_base=true;
			break;
		}



		bool base_pairing=false;
		for(Size i=1; i<=other_residues_base_list.size(); i++){
			base_stub const & base_info=other_residues_base_list[i];
			numeric::xyzVector<Real> const other_z_vector=base_info.base_coordinate_matrix.col_z();
			numeric::xyzVector<Real> rebuild_z_vector=moving_res_base.base_coordinate_matrix.col_z();

			numeric::xyzVector<Real> centroid_diff;
			subtract( moving_res_base.centroid, base_info.centroid, centroid_diff);

			Real centroid_distance=centroid_diff.length();
			if(centroid_distance<5 || centroid_distance>12) continue;

			Real base_z_offset_one=std::abs(dot( centroid_diff, other_z_vector));
			Real base_z_offset_two=std::abs(dot( centroid_diff, rebuild_z_vector));

			if(base_z_offset_one> 3 && base_z_offset_two>3) continue;

			Real base_axis_one=base_z_offset_one/centroid_distance;
			Real base_axis_two=base_z_offset_two/centroid_distance;
			if(base_axis_one> 0.5 && base_axis_two> 0.5) continue; //This is a stronger condition compare to baze_z_off_set check

			Real base_planarity=std::abs(dot( rebuild_z_vector, other_z_vector));
			if(base_planarity< 0.866 ) continue;

			numeric::xyzVector<Real> centroid_diff_parallel_one=dot( centroid_diff, other_z_vector)*other_z_vector; //messed with this on Jan 16, 2009 Parin S.
			numeric::xyzVector<Real> centroid_diff_perpendicular_one= centroid_diff-centroid_diff_parallel_one;
			Real rho_one=centroid_diff_perpendicular_one.length(); //length along xy plane
			numeric::xyzVector<Real> centroid_diff_parallel_two=dot( centroid_diff, rebuild_z_vector)*rebuild_z_vector;
			numeric::xyzVector<Real> centroid_diff_perpendicular_two= centroid_diff-centroid_diff_parallel_two;
			Real rho_two=centroid_diff_perpendicular_two.length();

			if((rho_one<5 || rho_one>10) && (rho_two<5 || rho_two>10)) continue;

			//If reach this point means success!
			base_pairing=true;
			break;
		}

		if ( stack_base ) {
			count_data.base_stack_count++;
//			return true; // found an interaction!
		}

		if ( base_pairing ){
			count_data.base_pairing_count++; //under_count, since if base_stack then the loop will already have short_circuited.
///			return true;
		}

		if ( !base_pairing && !stack_base ) return false; //This moving_res does not base_stack or base_pair with the base of any residue in the other_residues_base_list
		return true;

	}


	numeric::xyzVector<core::Real>
	get_base_centroid( conformation::Residue const & rsd , bool verbose){
		using namespace chemical;

  	assert( rsd.is_RNA() );

  	Vector centroid( 0.0 );
  	Size numatoms = 0;

		if(verbose)  std::cout << "Base atoms" << std::endl;
  	for ( Size i=rsd.first_sidechain_atom()+1; i<= rsd.nheavyatoms(); ++i ) { //rsd.first_sidechain_atom()+1 to not include the O2star oxygen.

		if(verbose) std::cout << "atom " << i  << " " << 	"name= " << rsd.type().atom_name(i) << " type= " << rsd.atom_type(i).name()  << " " << rsd.atom_type_index(i) << " " << rsd.atomic_charge(i);


		if(rsd.atom_type(i).name()=="VIRT"){
			if(verbose) std::cout << "  Virtual type: Ignore! " << std::endl;
			continue;
		}
		if(verbose) std::cout << std::endl;

    centroid += rsd.xyz(i);
    numatoms++;
  	}
  	centroid /= static_cast< Real >( numatoms );
  	return centroid;
	}


	// Should not copy code -- need to get rid of this or unify with RNA_CentroidInfo.cc
	///Stolen this from RNA_CentroidInfo.cc
	Matrix
	get_base_coordinate_system( conformation::Residue const & rsd, Vector const & centroid ){
  	using namespace chemical;
 	 Size res_type = rsd.aa();

  	assert( rsd.is_RNA() );

  	Vector x,y,z;

  	// Make an axis pointing from base centroid to Watson-Crick edge.
  	std::string WC_atom;
  	if ( res_type == na_rad ) WC_atom = " N1 ";
  	if ( res_type == na_rcy ) WC_atom = " N3 ";
  	if ( res_type == na_rgu ) WC_atom = " N1 ";
  	if ( res_type == na_ura ) WC_atom = " N3 ";

  	Vector const WC_coord (rsd.xyz( WC_atom ) );
  	x = WC_coord - centroid;
  	x.normalize();

  	// Make a perpendicular axis pointing from centroid towards
  	// Hoogstein edge (e.g., major groove in a double helix).
  	std::string H_atom;
  	if ( res_type == na_rad ) H_atom = "N7";
  	if ( res_type == na_rcy ) H_atom = "C5";
  	if ( res_type == na_rgu ) H_atom = "N7";
  	if ( res_type == na_ura ) H_atom = "C5";

    Vector const H_coord (rsd.xyz( H_atom ) );
  	y = H_coord - centroid; //not orthonormal yet...
  	z = cross(x, y);
  	z.normalize(); // Should poSize roughly 5' to 3' if in a double helix.

  	y = cross(z, x);
  	y.normalize(); //not necessary but doesn't hurt.


	 	Matrix M=Matrix::cols( x, y, z );
		return M;

	}

	//////////////////////////////////////////////////////////////////////////////////////
	PuckerState
	Get_residue_pucker_state(core::pose::Pose const & pose, Size const seq_num){

		using namespace core::scoring;
		using namespace core::scoring::rna;

		static scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info;
		Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );

		//		std::cout << "  DELTA_CUTOFF angle=" << DELTA_CUTOFF;

		conformation::Residue const & rsd(pose.residue(seq_num));
		Real const & delta( rsd.mainchain_torsion( DELTA ) );

		//		std::cout << "  delta angle=" << delta << std::endl;

		if (numeric::principal_angle_degrees( delta ) <= DELTA_CUTOFF) {
			return NORTH;
		} else {
			return SOUTH;
		}
	}


}
}
}
