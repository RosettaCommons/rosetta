// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @file
/// @brief


// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <utility/excn/Exceptions.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/import_pose/import_pose.hh>
//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
///////////////////////////////////////////////////
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh> //for EnergyMap
#include <core/scoring/EnergyMap.fwd.hh> //for EnergyMap
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/rms_util.tmpl.hh> 

#include <core/pack/pack_rotamers.hh> 
#include <core/pack/rotamer_trials.hh> 
#include <core/pack/task/PackerTask.hh> 
#include <core/pack/task/TaskFactory.hh> 

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

//////////////////////////////////////////////////////////
#include <protocols/idealize/idealize.hh>
#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/swa/rna/StepWiseRNA_Util.hh> 
#include <protocols/rna/RNA_ProtocolUtil.hh> 
#include <protocols/rna/RNA_BasePairClassifier.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/StepWiseClusterer.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGeneratorWrapper.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGeneratorWrapper.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseSugarRotamer.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseSugarRotamer.fwd.hh>
#include <protocols/rna/RNA_LoopCloser.hh>
#include <protocols/rna/RNA_LoopCloser.fwd.hh>

#include <core/scoring/rna/RNA_BaseDoubletClasses.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>  //Test
#include <cctype>
#include <iomanip>
#include <map>
#include <cstdlib>
#include <ctime>

//Added by Parin
//#include <core/scoring/ScoreType.hh>
#include <list>
#include <stdio.h>
#include <math.h>


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;
using namespace protocols::swa::rna;

typedef  numeric::xyzMatrix< Real > Matrix;



OPT_KEY( IntegerVector, delete_res )
OPT_KEY( Real, surrounding_radius)
OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, input_res )
OPT_KEY( IntegerVector, minimize_res )
OPT_KEY( StringVector, fold_tree_strings)
OPT_KEY( Boolean, fast )
OPT_KEY( String, 	algorithm)
OPT_KEY( StringVector, alignment_res_pairs)
OPT_KEY( StringVector, alignment_res_pairs_two)
OPT_KEY( StringVector, RMSD_res_pairs)
OPT_KEY( Real, alignment_RMSD_CUTOFF)
OPT_KEY( Boolean, graphic )
OPT_KEY( String, 	tag_name)
OPT_KEY( String, 	output_silent_file)
OPT_KEY( IntegerVector, lower_to_full_res)
OPT_KEY( IntegerVector, upper_to_full_res)
OPT_KEY( IntegerVector, helical_ends_to_full_res)
OPT_KEY( IntegerVector, cutpoint_closed)
OPT_KEY( String, filter_filename)
OPT_KEY( IntegerVector, virtual_res)
OPT_KEY( IntegerVector, virtual_ribose)
OPT_KEY( IntegerVector, native_virtual_res)
OPT_KEY( IntegerVector, native_alignment_res)
OPT_KEY( IntegerVector, rmsd_res)
OPT_KEY( String, 	helical_ends)
OPT_KEY( Integer, user_num_nstruct_per_node)
OPT_KEY( Integer, user_JOB_ID)
OPT_KEY( Integer, user_JOB_ID_MOD_CUTOFF)
OPT_KEY( Boolean, USER_BIOX_SUBMIT)
OPT_KEY( Boolean, user_skip_minimize)
OPT_KEY( Boolean, user_extra_minimize_rounds)
OPT_KEY( Boolean, align_only_over_base_atoms )
OPT_KEY( IntegerVector, additional_slice_res)
OPT_KEY( Boolean, INCLUDE_EDGE_PHOSPHATE )
OPT_KEY( Boolean, double_count_base_pair )
OPT_KEY( Real, rotamer_cluster_rmsd)
OPT_KEY( String, dinucleotide_sequence)
OPT_KEY( Boolean, cluster_rotamers_optimize_screening )
OPT_KEY( Boolean, two_stage_rotamer_clustering )
OPT_KEY( Boolean, quick_test)
OPT_KEY( Boolean, cluster_rotamer_sparse_output)
OPT_KEY( Integer, cluster_rotamer_bin_size)
OPT_KEY( Boolean, cluster_rotamer_replusion_screen)
OPT_KEY( Boolean, cluster_rotamer_VDW_rep_screening_slow_check)
OPT_KEY( Boolean, idl_close_chainbreaks)			  //FOR rna_idealize_test()
OPT_KEY( Real, atom_pair_constraint_weight )  //FOR rna_idealize_test()
OPT_KEY( Real, coordinate_constraint_weight ) //FOR rna_idealize_test()
OPT_KEY( StringVector, list_of_virtual_res )
OPT_KEY( RealVector, list_of_energy ) 
OPT_KEY( String, native_tag_name ) 
OPT_KEY( StringVector, decoy_tag_name ) 
OPT_KEY( Boolean, dump ) 


//////////////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP 
create_scorefxn(){ //Copy from rna_swa_test.cc on Oct 11, 2011

	using namespace core::scoring;


	std::string score_weight_file;

	Size num_score_weight_file=0;

	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file= option[ basic::options::OptionKeys::score::weights ]();
		std::cout << "User passed in score:weight option: " << score_weight_file << std::endl;
		num_score_weight_file++;
	}	


	if(num_score_weight_file==0){
		score_weight_file="rna_loop_hires_07232011.wts";
		std::cout << "Using default score_weight_file=" << score_weight_file << std::endl;
	}

	if(num_score_weight_file>1){
		std::cout << "num_score_weight_file (inputted by user)=" << num_score_weight_file << std::endl;
		utility_exit_with_message("num_score_weight_file>1");
	}
	
	core::scoring::ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( score_weight_file );

	std::cout << "---------score function weights----------" << std::endl;
	scorefxn->show(std::cout);
	std::cout << "-----------------------------------------" << std::endl;


	return scorefxn; 
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
multiple_variant_type_test(){
	
	  using namespace core::pose;
	  using namespace core::chemical;
	  using namespace core::kinematics;
	  using namespace core::scoring;
		using namespace protocols::swa::rna;

		ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		std::string pdb_tag= option[ in::file::s ]()[1] ;

		Size seq_num=5;

		pose::Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set,  pdb_tag );
		std::string output_pose_name="";

		/////////////////////////////////////////////////////////////////////////////////////////////////

		std::cout << "ADD VIRTUAL_RNA_RESIDUE_UPPER" << std::endl;
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_UPPER", seq_num);
		output_pose_name+="iAU";
		dump_pdb( pose, output_pose_name+ ".pdb" );	

		Output_boolean(" VIRTUAL_PHOSPHATE= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_PHOSPHATE" ), TR );
		Output_boolean(" VIRTUAL_RNA_RESIDUE_UPPER= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ), TR );
		std::cout << std::endl;
	
		/////////////////////////////////////////////////////////////////////////////////////////////////

		std::cout << "ADD VIRTUAL_PHOSPHATE" << std::endl;
		pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", seq_num);
		output_pose_name+="iAP";
		dump_pdb( pose, output_pose_name+ ".pdb" );	
	
		
		Output_boolean(" VIRTUAL_PHOSPHATE= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_PHOSPHATE" ), TR );
		Output_boolean(" VIRTUAL_RNA_RESIDUE_UPPER= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ), TR );
		std::cout << std::endl;

	
		/////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout << "REMOVE VIRTUAL_PHOSPHATE " << std::endl;
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", seq_num);
		output_pose_name+="iRP";
		dump_pdb( pose, output_pose_name+ ".pdb" );	

		Output_boolean(" VIRTUAL_PHOSPHATE= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_PHOSPHATE" ), TR );
		Output_boolean(" VIRTUAL_RNA_RESIDUE_UPPER= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ), TR );
		std::cout << std::endl;

		std::cout << "REMOVE VIRTUAL_RNA_RESIDUE_UPPER " << std::endl;
		pose::remove_variant_type_from_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_UPPER", seq_num);
		output_pose_name+="iRU";
		dump_pdb( pose, output_pose_name+ ".pdb" );		

		Output_boolean(" VIRTUAL_PHOSPHATE= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_PHOSPHATE" ), TR );
		Output_boolean(" VIRTUAL_RNA_RESIDUE_UPPER= ", pose.residue( seq_num  ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ), TR );
		std::cout << std::endl;

		/////////////////////////////////////////////////////////////////////////////////////////////////
	

}


void
align_pdbs_function(pose::Pose const static_pose,
									utility::vector1< pose_data_struct2 > & moving_pose_data_list, 
									utility::vector1< std::string > const & alignment_res_pair_list,
									Real const alignment_RMSD_cutoff){

	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::rna;
	using namespace protocols::swa::rna;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	std::cout << "alignment_RMSD_cutoff= " << alignment_RMSD_cutoff << std::endl;


	if(alignment_res_pair_list.size()==0){
		utility_exit_with_message( "alignment_res_pair_list.size()==0" );
	}

	utility::vector1< core::Size > static_pdb_align_res;
	utility::vector1< core::Size > moving_pdb_align_res;

	for(Size n=1; n<=alignment_res_pair_list.size(); n++){

		utility::vector1< std::string > const alignment_res_pair=Tokenize(alignment_res_pair_list[n], "-");
		if(alignment_res_pair.size()!=2){
			 utility_exit_with_message( "alignment_res_pair.size()!=2, alignment_res_pair_list[n]= " + alignment_res_pair_list[n] );
		}

		static_pdb_align_res.push_back(string_to_int(alignment_res_pair[1]));
		moving_pdb_align_res.push_back(string_to_int(alignment_res_pair[2]));
	}

	if(static_pdb_align_res.size()!=moving_pdb_align_res.size()){
		utility_exit_with_message( "static_pdb_align_res.size()!=moving_pdb_align_res.size()" );
	}

	std::cout << "static_pdb_align_res to moving_pdb_align_res:" << std::endl;
	for(Size ii=1; ii<=static_pdb_align_res.size(); ii++){
		std::cout << static_pdb_align_res[ii] << " ---> " << moving_pdb_align_res[ii] << std::endl;
	}


	bool output_pdb=true;
	bool output_silent_file=true;

	for(Size n=1; n<=moving_pose_data_list.size(); n++){

		pose::Pose & moving_pose=*(moving_pose_data_list[n].pose_OP);
		std::string const moving_pdb_tag=moving_pose_data_list[n].tag;

		////////////////////////////////////create the alignment map////////////////////////////////////////////////////////////////////
		id::AtomID_Map< id::AtomID > atom_ID_map; 
		pose::initialize_atomid_map( atom_ID_map, moving_pose, id::BOGUS_ATOM_ID );

		std::string const static_sequence=static_pose.sequence();
		std::string const moving_sequence=moving_pose.sequence();


		for(Size ii=1; ii<=static_pdb_align_res.size(); ii++){

			Size const static_seq_num=static_pdb_align_res[ii];
			Size const moving_seq_num=moving_pdb_align_res[ii];

			if( (static_sequence[static_seq_num-1])!=(moving_sequence[moving_seq_num-1]) ){
				std::cout << "static_seq_num= " << static_seq_num << " static_sequence= " << static_sequence << " static_sequence[static_seq_num-1]= " << static_sequence[static_seq_num-1] << std::endl;
				std::cout << "moving_sequence= " << moving_seq_num << " moving_sequence= " << moving_sequence << " moving_sequence[moving_seq_num-1]= " << moving_sequence[moving_seq_num-1] << std::endl;
				utility_exit_with_message( "(static_sequence[static_seq_num-1])!=(moving_sequence[moving_seq_num-1])" );				
			}

			setup_suite_atom_id_map(moving_pose.residue(moving_seq_num), static_pose.residue(static_seq_num),  atom_ID_map);

		}
		core::scoring::superimpose_pose(moving_pose, static_pose, atom_ID_map);		

		Size total_atom_count=0;
		Real total_sum_sd=0.0;

		for(Size ii=1; ii<=static_pdb_align_res.size(); ii++){

			Size const moving_seq_num=moving_pdb_align_res[ii];
			Size const static_seq_num=static_pdb_align_res[ii];

			//bool verbose= (ii==1) ? true : false;
			bool verbose= false;	

			Size atom_count=0;
			Real sum_sd=0.0;

			base_atoms_square_deviation(moving_pose, static_pose, moving_seq_num, static_seq_num, atom_count, sum_sd, verbose , false /*ignore_virtual_atom*/);

			sum_sd=sum_sd/(atom_count);
			Real rmsd=sqrt(sum_sd);

			if(atom_count==0) rmsd=99.99; //This is different from suite_rmsd function..!!


			if(rmsd>alignment_RMSD_cutoff ){ //change on Sept 26, 2010..problem arise when use this in non-long-loop mode...
				std::cout << "rmsd= " << rmsd  << " is greater than " << alignment_RMSD_cutoff << " Angstrom between res " << moving_seq_num << " of moving_pose and res " << static_seq_num << " of static_pose" << std::endl;
				std::cout << "moving_pdb_tag= " << moving_pdb_tag << std::endl;
				utility_exit_with_message( "rmsd>alignment_RMSD_cutoff!"); 
			}

			total_atom_count+=atom_count;
			total_sum_sd+=sum_sd;

		}

		total_sum_sd=total_sum_sd/(total_atom_count);
		Real all_base_rmsd=sqrt(total_sum_sd);

	}

}

/////////////////////////////////////////////////////////////////////////////////////////////////
utility::vector1< pose_data_struct2 >
convert_silent_file_to_pose_data_list(std::string const silent_file){

	using namespace core::pose;
	using namespace ObjexxFCL;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::scoring;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );


	core::import_pose::pose_stream::SilentFilePoseInputStreamOP silent_file_stream;

	silent_file_stream= new core::import_pose::pose_stream::SilentFilePoseInputStream();
	silent_file_stream->set_order_by_energy( true );
	silent_file_stream->set_record_source( true );

	utility::vector1< std::string > singleton_list;
	singleton_list.push_back( silent_file );
	silent_file_stream->filenames( singleton_list ); //triggers read in of files, too.

	utility::vector1< pose_data_struct2 > pose_data_list;

	while ( silent_file_stream->has_another_pose() ) {

		core::io::silent::SilentStructOP const silent_struct( silent_file_stream->next_struct() );

		PoseOP pose_op( new Pose );

		silent_struct->fill_pose( *pose_op,*(rsd_set)  ); 

		Real score( 0.0 );
		getPoseExtraScores( *pose_op, "score", score );

		std::string const & tag( silent_struct->decoy_tag() );

		if(protocols::swa::rna::check_for_messed_up_structure((*pose_op), tag ) ) {
			utility_exit_with_message("tag= " + tag  + " is messed up!");
		}

		pose_data_struct2 pose_data;
		pose_data.pose_OP=pose_op;
		pose_data.score=score;
		pose_data.tag=tag;

		pose_data_list.push_back(pose_data);

	}

	return pose_data_list;

}

class Combine_Tags_Info{

	public:

	Combine_Tags_Info():
		lower_tag( "" ),
		upper_tag( "" ),
		combine_score( 999999999999.99 ),
		overlap_rmsd( 999999999999.99 )
	{
	}

	~Combine_Tags_Info(){};

	public:

		std::string lower_tag;
		std::string upper_tag;
		core::Real combine_score;
		core::Real overlap_rmsd;
};

void
copy_virtual_variant_type(pose::Pose & full_pose, pose::Pose const & start_pose_with_variant, utility::vector1< core::Size > const & input_res_map){

	using namespace core::chemical;

	for( Size n = 1; n <= start_pose_with_variant.total_residue(); n++  ) {
		if(input_res_map[n]==0){
			continue;
		}

		if(start_pose_with_variant.residue(n).has_variant_type("VIRTUAL_RNA_RESIDUE")){
	
			if( (n+1)>start_pose_with_variant.total_residue()){ //Check in range
				std::cout << "(n+1)= " << (n+1)  << std::endl;
				utility_exit_with_message( "(n+1)>start_pose_with_variant.total_residue()!" );
			} 

			if(start_pose_with_variant.residue(n+1).has_variant_type("VIRTUAL_RNA_RESIDUE_UPPER")==false){
				std::cout << "n= " << n << std::endl;
				utility_exit_with_message("res n has_variant_type VIRTUAL_RNA_RESIDUE but res n+1 does not have variant_type VIRTUAL_RNA_RESIDUE_UPPER");
			}

			apply_virtual_rna_residue_variant_type( full_pose, input_res_map[n] , false /*apply_check*/ ) ;
		} 

		if(start_pose_with_variant.residue(n).has_variant_type("VIRTUAL_RIBOSE")){
			add_variant_type_to_pose_residue( full_pose, "VIRTUAL_RIBOSE", input_res_map[n] );

		} 
	}
}

/*Commented out on Jan 16, 2012 after implemented a more general version of this function in StepWiseRNA_Util.cc.
void
remove_all_variant_types(pose::Pose & pose){

	using namespace core::pose;



	for ( Size n = 1; n <= pose.total_residue(); n++  ) {
		remove_variant_type_from_pose_residue( pose, "VIRTUAL_PHOSPHATE", n );
		remove_variant_type_from_pose_residue( pose, "VIRTUAL_O2STAR_HYDROGEN", n );
		remove_variant_type_from_pose_residue( pose, "CUTPOINT_LOWER", n );
		remove_variant_type_from_pose_residue( pose, "CUTPOINT_UPPER", n );
		remove_variant_type_from_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", n );
		remove_variant_type_from_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_UPPER", n );
		remove_variant_type_from_pose_residue( pose, "BULGE", n );	
		remove_variant_type_from_pose_residue( pose, "VIRTUAL_RIBOSE", n );	
	}

}
*/

void
copy_DOFS_local(pose::Pose & pose, pose::Pose const start_pose, utility::vector1< core::Size >  const & input_res){

	using namespace core::chemical;

	std::map< core::Size, core::Size > res_map;  //This is map from sub numbering to input_res numbering..
	for ( Size n = 1; n <= input_res.size(); n++ ) {
		if(input_res[n]==0) continue;
		res_map[ input_res[n] ] = n;
	}

	//Does this work for the "overlap residue" case?? If there is a overlap residue, then order of input_res will manner...Parin Jan 2, 2010.
	//copy_dofs( pose, start_pose, res_map, true /*copy_dofs_for_junction_residues*/ );
	copy_dofs_match_atom_names( pose, start_pose, res_map, false /*backbone_only*/, false /*ignore_virtual*/); //Dec 28, 2011

}

void
hermann_phase_two_minimize(){

	using namespace core::pose;
	using namespace ObjexxFCL;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::scoring;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;
	using namespace core::optimization;
	//using namespace protocols::rna;
	using namespace core::id;
	using namespace core::conformation;
	using namespace core::scoring::constraints;

	clock_t const time_start( clock() ); 

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	
	protocols::rna::RNA_LoopCloser rna_loop_closer;

	bool const biox_submit=option[ USER_BIOX_SUBMIT ]();

	bool output_silentfile_every_round=true;
	bool output_copy_DOF=true;

	if(biox_submit){
		output_copy_DOF=false; 
		output_silentfile_every_round=false;
	}

	bool const skip_minimize=option[ user_skip_minimize ]();
	bool const extra_minimize_rounds=option[ user_extra_minimize_rounds] ();
	Size const num_extra_rounds=20;


	Output_boolean("biox_submit= ", biox_submit, TR ); std::cout << std::endl;
	Output_boolean("output_silentfile_every_round= ", output_silentfile_every_round, TR ); std::cout << std::endl;
	Output_boolean("output_copy_DOF= ", output_copy_DOF, TR ); std::cout << std::endl;
	Output_boolean("skip_minimize= ", skip_minimize, TR ); std::cout << std::endl;
	Output_boolean("extra_minimize_rounds= ", extra_minimize_rounds, TR ); std::cout << std::endl;


	bool copy_DOF=false;

	pose::Pose helical_ends_pose;

	utility::vector1< core::Size > helical_ends_to_full_map;	
	if( option[ helical_ends ].user() ){
 
		if ( !option[ helical_ends_to_full_res ].user() ) utility_exit_with_message( "User must supply helical_ends_to_full_res in copy_DOF mode!" );
		helical_ends_to_full_map=option[ helical_ends_to_full_res ]();

		copy_DOF=true;
		std::string const helical_ends_pdb=option[ helical_ends ]();
		import_pose::pose_from_pdb( helical_ends_pose, *rsd_set, helical_ends_pdb );
	}



	if ( !option[ lower_to_full_res ].user() ) utility_exit_with_message( "User must supply lower_to_full_res!" );
	if ( !option[ upper_to_full_res ].user() ) utility_exit_with_message( "User must supply upper_to_full_res!" );
	if ( !option[ cutpoint_closed ].user() ) utility_exit_with_message( "User must supply cutpoint_closed!" );
	if ( !option[ filter_filename ].user() ) utility_exit_with_message( "User must supply filter_filename" );
	if ( !option[ minimize_res ].user() ) utility_exit_with_message( "User must supply minimize_res!" );
	if ( !option[ fold_tree_strings ].user() ) utility_exit_with_message( "User must supply fold_tree_strings!" );
	if ( !option[ out::file::silent ].user() ) utility_exit_with_message( "User must supply out::file::silent!" );


	if ( !option[ in::file::native ].user() ) utility_exit_with_message( "User must supply in::file::native!" ); //This is the static pose.
	if ( !option[ native_alignment_res ].user() ) utility_exit_with_message( "User must supply in::file::native!" ); //This is the static pose.
	if ( !option[ rmsd_res ].user() ) utility_exit_with_message( "User must supply in::file::native!" ); //This is the static pose.

	//For parallelizing the jobs
	if ( !option[ user_num_nstruct_per_node ].user() ) utility_exit_with_message( "User must supply user_num_nstruct_per_node!" ); //
	if ( !option[ user_JOB_ID ].user() ) utility_exit_with_message( "User must supply user_JOB_ID!" ); //start from 0
	if ( !option[ user_JOB_ID_MOD_CUTOFF ].user() ) utility_exit_with_message( "User must supply user_JOB_ID_MOD_CUTOFF!" ); //start from 0
	if ( !option[ in::file::fasta  ].user() ) utility_exit_with_message( "User must supply in::file::fasta !" );


	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const full_sequence = fasta_sequence->sequence();
	core::Size const total_res=full_sequence.length();


	utility::vector1< Size > const native_alignment_res_list	=option[ native_alignment_res ]();
	utility::vector1< Size > const native_virtual_res_list		=option[ native_virtual_res ]();
	utility::vector1< Size > const rmsd_res_list						=option[ rmsd_res ]();

	if(native_alignment_res_list.size()==0) utility_exit_with_message( "native_alignment_res_list=.size()==0" );
	if(rmsd_res_list.size()==0) utility_exit_with_message( "rmsd_res_list.size()==0" );


	SilentFileData silent_file_data;
	std::string const silent_outfile=option[ out::file::silent ]();


	std::string const native_filename = option[ in::file::native ]();

	pose::Pose native_pose_ACT;
	import_pose::pose_from_pdb( native_pose_ACT, *rsd_set, native_filename );

	if(native_pose_ACT.sequence()!=full_sequence){
		std::cout << "native_pose_ACT.sequence()= " << native_pose_ACT.sequence() << std::endl;
		std::cout << "full_sequence= " << full_sequence << std::endl;
		utility_exit_with_message( "native_pose_ACT.sequence()!=full_sequence!" );
		
	}

	utility::vector1< core::Size > const cutpoint_closed_list= option[ cutpoint_closed  ]();	

	utility::vector1< core::Size > const minimize_res_list= option[ minimize_res ]();
	
	Output_seq_num_list("minimize_res_list= ", minimize_res_list, TR, 30 );

	utility::vector1< std::string > const	fold_tree_string_list= option[ fold_tree_strings ]();	
	std::string const filtered_tag_file=option[ filter_filename  ]();



	//Copy DOF? (For now, take care of this with python..)

	////////////////////////create score functions//////////////////////////////////////////////////////
	utility::vector1< core::scoring::ScoreFunctionOP > scorefxn_list;
	scorefxn_list.clear();

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn(); 

	scorefxn->set_weight( angle_constraint, 5.0 );
	scorefxn->set_weight( atom_pair_constraint, 5.0 );

	core::scoring::ScoreFunctionOP close_CB_scorefxn = scorefxn->clone();

	close_CB_scorefxn->set_weight( linear_chainbreak, 0.0); //no point having linear_CB since the torsion is not minimized anyways. 
	close_CB_scorefxn->set_weight( coordinate_constraint, 0.00 );



	scorefxn_list.push_back(close_CB_scorefxn);

	if(extra_minimize_rounds){
		for(Size n=1; n<=num_extra_rounds; n++){
			core::scoring::ScoreFunctionOP extra_round_scorefxn= scorefxn->clone();
			Real coord_const_weight=(n-1)*0.02;

			/*
			if(n==1){
				coord_const_weight=0.0;
			}else if(n==2){
				coord_const_weight=0.01;
			}else if(n==3){
				coord_const_weight=0.03;
			}else{
				coord_const_weight=0.05;
			}
			*/
		
			extra_round_scorefxn->set_weight( coordinate_constraint, coord_const_weight );
			scorefxn_list.push_back(extra_round_scorefxn);
		}

		core::scoring::ScoreFunctionOP final_scorefxn = scorefxn->clone();
		final_scorefxn->set_weight( coordinate_constraint, num_extra_rounds*0.02 );
		scorefxn_list.push_back(final_scorefxn);

	}else{

		core::scoring::ScoreFunctionOP final_scorefxn = scorefxn->clone();
		final_scorefxn->set_weight( coordinate_constraint, 0.10 );
		//final_scorefxn->set_weight( coordinate_constraint, 0.40 ); #Minimization not smooth with this weight + only 2 rounds.
		scorefxn_list.push_back(final_scorefxn);

	}



	std::cout << "total_min_rounds= " << scorefxn_list.size() << " num_extra_rounds= " << num_extra_rounds << std::endl;

	////////////////////////////////////Setup fold_tree...////////////////////////////////////////////////////////////////////
	Output_title_text("Setup FOLD_TREE", TR );
	core::kinematics::FoldTree fold_tree( total_res );

	for(Size ii=1; ii<=fold_tree_string_list.size(); ii++){

		utility::vector1< std::string > const fold_tree_string=Tokenize(fold_tree_string_list[ii], "-");
		if(fold_tree_string.size()!=3){
			 utility_exit_with_message( "fold_tree_string.size()!=3, fold_tree_string_list[ii]= " + fold_tree_string_list[ii] );
		}
		Size const five_prime_seq_num=string_to_int(fold_tree_string[1]);
		Size const three_prime_seq_num=string_to_int(fold_tree_string[2]);
		Size const cut_point=string_to_int(fold_tree_string[3]);

		std::cout << "five_prime_seq_num= " << five_prime_seq_num << " three_prime_seq_num= " << three_prime_seq_num << " cut_point= " << cut_point << std::endl;

		fold_tree.new_jump( five_prime_seq_num, three_prime_seq_num, cut_point );
	}		
	Output_fold_tree_info(fold_tree, "fold_tree derived from user inputted fold_tree_strings");


	////////////////////////setup coresponding movemap_list//////////////////////////////////////////////////////
	core::kinematics::MoveMap fixed_mm;

	fixed_mm.set_bb( false );
	fixed_mm.set_chi( false );
	fixed_mm.set_jump( false );

	utility::vector1< core::kinematics::MoveMap> mm_list;

	for (Size n = 1; n <= fold_tree.num_jump(); n++ ){
		Size const jump_pos1( fold_tree.upstream_jump_residue( n ) );
		Size const jump_pos2( fold_tree.downstream_jump_residue( n ) );

		if( ((jump_pos1+1)==(jump_pos2)) || ((jump_pos1==1) && (jump_pos2==total_res)) ){
			fixed_mm.set_jump( n, false );
			std::cout << "jump_pos1= " << jump_pos1 << " jump_pos2= " << jump_pos2 << " mm.jump= false"; ;  std::cout << std::endl;
		}else{
			fixed_mm.set_jump( n, true );
			std::cout << "jump_pos1= " << jump_pos1 << " jump_pos2= " << jump_pos2 << " mm.jump= true"; ;  std::cout << std::endl;
		}
	}
	
	core::kinematics::MoveMap mm_first=fixed_mm;
	core::kinematics::MoveMap mm_final=fixed_mm;

	for(Size n=1; n<=minimize_res_list.size(); n++){
		Size const seq_num=minimize_res_list[n];
		mm_final.set_bb(seq_num, true);
		mm_final.set_chi(seq_num, true );
	}	
	

	mm_list.push_back(mm_first);

	if(extra_minimize_rounds){
		for(Size n=1; n<=num_extra_rounds; n++){
			mm_list.push_back(mm_final);
		}
	}
	
	mm_list.push_back(mm_final);

	if(scorefxn_list.size()!=mm_list.size()){
		utility_exit_with_message( "scorefxn_list.size()"+string_of(scorefxn_list.size())+"!mm_list.size()("+string_of(mm_list.size())+")" );
	}
	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	Output_title_text("IMPORT_FILTERED_FILENAME", TR );

	Size const num_nstruct_per_node=option[ user_num_nstruct_per_node ]();
	Size const JOB_ID_MOD_CUTOFF=option[ user_JOB_ID_MOD_CUTOFF ]();
	Size JOB_ID=							option[ user_JOB_ID ]();

	std::cout << "num_nstruct_per_node=" << num_nstruct_per_node << std::endl;
	std::cout << "JOB_ID_MOD_CUTOFF= " << JOB_ID_MOD_CUTOFF << std::endl;
	std::cout << "USER_JOB_ID=" << JOB_ID << std::endl;

	if(JOB_ID>=JOB_ID_MOD_CUTOFF){
		JOB_ID=JOB_ID-JOB_ID_MOD_CUTOFF;
	}
							
	std::cout << "ACT_JOB_ID=" << JOB_ID << std::endl;							

	
	Size const start_nstruct=(num_nstruct_per_node*JOB_ID)+1;
	Size const end_nstruct=(num_nstruct_per_node*(JOB_ID+1));

	std::cout << "start_nstruct=" << start_nstruct << std::endl;
	std::cout << "end_nstruct=" << end_nstruct << std::endl;


	//extract LOW RMSD tag from filter_filename.	
	std::ifstream infile;	
 	infile.open(filtered_tag_file.c_str());

	if (infile.fail()){
	 utility_exit_with_message("Error! \"" + filtered_tag_file + "\" could not be opened!");
	}else{
		std::cout << "Open \"" << filtered_tag_file << "\" successful!" << std::endl;
	}

	std::string line_string;

	utility::vector1< Combine_Tags_Info > filtered_tag_pair_info_list;
	filtered_tag_pair_info_list.clear();

	bool first_line=true;

	Size curr_nstruct=0;

	while(getline(infile, line_string) ){	

		utility::vector1< std::string > const line_list=Tokenize(line_string," \t\n\f\v"); //Oct 19, 2010..now filterer_outfile contain other terms.

		if(first_line){
			if(line_list[1]!="lower_tag")  utility_exit_with_message( "line_list[1]!=\"lower_tag");
			if(line_list[2]!="upper_tag")  utility_exit_with_message( "line_list[2]!=\"upper_tag");
			if(line_list[3]!="sum_score")  utility_exit_with_message( "line_list[3]!=\"sum_score");
			if(line_list[4]!="base_rmsd")  utility_exit_with_message( "line_list[4]!=\"base_rmsd");                                                              
			first_line=false;
			continue;
		}

		curr_nstruct++;
		if(curr_nstruct<start_nstruct) continue;
		if(curr_nstruct>end_nstruct) continue;

		Combine_Tags_Info combine_tag_info;
		combine_tag_info.lower_tag=line_list[1];
		combine_tag_info.upper_tag=line_list[2];
		combine_tag_info.combine_score=string_to_real(line_list[3]);
		combine_tag_info.overlap_rmsd=string_to_real(line_list[4]);

		std::cout << line_string << " curr_nstruct= " << curr_nstruct << std::endl;

		filtered_tag_pair_info_list.push_back(combine_tag_info);
	}

	if(filtered_tag_pair_info_list.size()!=num_nstruct_per_node){
		utility_exit_with_message( "filtered_tag_pair_info_list.size()=("+ string_of(filtered_tag_pair_info_list.size()) +")!=("+ string_of(num_nstruct_per_node) +")=num_nstruct_per_node" );
	}

	infile.close();

	
	//////////////////////////////////Import lower_to_full_map and upper_to_full_map///////////////////////////////////////////////////

	utility::vector1< core::Size > lower_to_full_map=option[ lower_to_full_res ]();
	utility::vector1< core::Size > upper_to_full_map=option[ upper_to_full_res ]();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////


	utility::vector1< std::string > silent_files_list = option[ in::file::silent ]();

	if(silent_files_list.size()!=2) utility_exit_with_message( "silent_files_list.size()!=2" );
	
	std::string const upper_silent_file=silent_files_list[2];


	utility::vector1< pose_data_struct2 > lower_pose_data_list=convert_silent_file_to_pose_data_list(silent_files_list[1]);
	utility::vector1< pose_data_struct2 > upper_pose_data_list=convert_silent_file_to_pose_data_list(silent_files_list[2]);

	bool verbose=true;

	pose::Pose full_pose;

	if(option[ graphic ]()) protocols::viewer::add_conformation_viewer( full_pose.conformation(), "minimize_hermann_duplex", 400, 400 );

	Size num_tag_pair_found=0;

	for(Size lower_ID=1; lower_ID<=lower_pose_data_list.size(); lower_ID++){	
	for(Size upper_ID=1; upper_ID<=upper_pose_data_list.size(); upper_ID++){

		pose_data_struct2 const & lower_pose_data=lower_pose_data_list[lower_ID];
		pose_data_struct2 const & upper_pose_data=upper_pose_data_list[upper_ID];

		//Add virtual variant type
		pose::Pose const & lower_pose=*(lower_pose_data.pose_OP);
		pose::Pose const & upper_pose=*(upper_pose_data.pose_OP);


		//find low RMSD match
		Size found_tag_pair=0;
		std::string full_pose_tag="";
		std::string import_pose_tag="";

		Combine_Tags_Info FOUND_combine_tag_info;

		for(Size tag_pair_ID=1; tag_pair_ID<=filtered_tag_pair_info_list.size(); tag_pair_ID++){
			Combine_Tags_Info const & combine_tag_info=filtered_tag_pair_info_list[tag_pair_ID];
			if( (combine_tag_info.lower_tag==("aligned_lower_"+lower_pose_data.tag) )  && (combine_tag_info.upper_tag==("aligned_upper_"+upper_pose_data.tag) ) ){
				std::cout << "Found tag_pair: lower_tag= " << combine_tag_info.lower_tag << " upper_tag= " << combine_tag_info.upper_tag << std::endl;
				full_pose_tag="S_lower_"+ lower_pose_data.tag + "_" + "upper_" + upper_pose_data.tag;
				import_pose_tag="COMBINE_PDB/"+ combine_tag_info.lower_tag + "_" + combine_tag_info.upper_tag + ".pdb" ;

				found_tag_pair++;
				FOUND_combine_tag_info=	combine_tag_info;			
			}
		}

		if(found_tag_pair>1){
			utility_exit_with_message( "found_tag_pair ("+ string_of(found_tag_pair)+ ")>1" );
		}

		if(found_tag_pair==0) continue;
		
		num_tag_pair_found++;
		std::cout << "SO_FAR: num_tag_pair_found= " << num_tag_pair_found << " time_taken=" << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC <<std::endl;
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
		if(copy_DOF){
			Output_title_text("COPY DOFS", TR );

			pose::Pose copy_DOF_pose;

			make_pose_from_sequence( copy_DOF_pose, full_sequence, *rsd_set, false /*auto_termini*/);

			full_pose=copy_DOF_pose; //This ensures that everything in full_pose is initialized from scratch
	
			full_pose.fold_tree( fold_tree ); 

			pose::Pose lower_pose_without_variant_types=	lower_pose;
			pose::Pose upper_pose_without_variant_types=	upper_pose;
			pose::Pose helical_ends_without_variant_types= helical_ends_pose;

			utility_exit_with_message("CODE CONTAIN BROKEN CODE THAT HAVE NOT YET BEEN FIXED");

			//Commented out on Jan 12, 2012...This used to called the version of remove_all_variant_types() in parin_test.cc which has been commented out and replaced with the one StepWsieRNA_Util.cc
			//If want to use this code, then need to ensure that the version of remove_all_variant_types() in StepWsieRNA_Util.cc reproduce the results of the local remove_all_variant_types() function
			//remove_all_variant_types(lower_pose_without_variant_types); 
			//remove_all_variant_types(upper_pose_without_variant_types); 
			//remove_all_variant_types(helical_ends_without_variant_types); 
			//////////////////////////////////////////////////////////////////////////////////

			copy_DOFS_local(full_pose, helical_ends_without_variant_types, helical_ends_to_full_map);
		
			if(verbose && output_copy_DOF) dump_pdb(full_pose, "copy_DOF_1_" + full_pose_tag + ".pdb");

			copy_DOFS_local(full_pose, lower_pose_without_variant_types, lower_to_full_map);

			if(verbose && output_copy_DOF) dump_pdb(full_pose, "copy_DOF_2_" + full_pose_tag + ".pdb");


			copy_DOFS_local(full_pose, upper_pose_without_variant_types, upper_to_full_map);

			if(verbose && output_copy_DOF) dump_pdb(full_pose, "copy_DOF_3_" + full_pose_tag + ".pdb");

			//continue;
		}else{
			Output_title_text("IMPORT PDB", TR );
			//import PDB? CAN JUST COPY DOF? MIGHT NEED the alignment helix...

			std::cout << "importing " << import_pose_tag << std::endl;

			pose::Pose import_pose;
			import_pose::pose_from_pdb( import_pose, *rsd_set, import_pose_tag );

			full_pose=import_pose; //This ensures that everything in full_pose is initialized from scratch

			for ( Size n = 1; n <= full_pose.total_residue(); n++  ) {
				remove_variant_type_from_pose_residue( full_pose, "VIRTUAL_PHOSPHATE", n );
				remove_variant_type_from_pose_residue( full_pose, "VIRTUAL_O2STAR_HYDROGEN", n );
				remove_variant_type_from_pose_residue( full_pose, "CUTPOINT_LOWER", n );
				remove_variant_type_from_pose_residue( full_pose, "CUTPOINT_UPPER", n );
				remove_variant_type_from_pose_residue( full_pose, "VIRTUAL_RNA_RESIDUE", n );
				remove_variant_type_from_pose_residue( full_pose, "VIRTUAL_RNA_RESIDUE_UPPER", n );
				remove_variant_type_from_pose_residue( full_pose, "BULGE", n );	
				remove_variant_type_from_pose_residue( full_pose, "VIRTUAL_RIBOSE", n );	
			}

			full_pose.fold_tree( fold_tree );

		}

		protocols::rna::assert_phosphate_nomenclature_matches_mini(full_pose);

		if(full_pose.total_residue()!=total_res){
			utility_exit_with_message( "full_pose.total_residue()(" + string_of(full_pose.total_residue())+")!=total_res("+string_of(total_res)+")" );
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if(lower_pose.total_residue()!=lower_to_full_map.size()){
			utility_exit_with_message( "lower_pose.total_residue()(" + string_of(lower_pose.total_residue())+")!=lower_to_full_map.size()("+string_of(lower_to_full_map.size())+")" );
		}

		if(upper_pose.total_residue()!=upper_to_full_map.size()){
			utility_exit_with_message( "upper_pose.total_residue()(" + string_of(upper_pose.total_residue())+")!=upper_to_full_map.size()("+string_of(upper_to_full_map.size())+")" );
		}
		

		copy_virtual_variant_type(full_pose, lower_pose, lower_to_full_map);
		copy_virtual_variant_type(full_pose, upper_pose, upper_to_full_map);

		////////////////////////////////////////////////////////////////////////////////////////////////////////////
		
		//Setup cutpoint upper and lower variant type.


		for(Size n=1; n<=cutpoint_closed_list.size(); n++){

			Size const cutpoint_closed=cutpoint_closed_list[n]; 	

			if ( cutpoint_closed == 0 ) utility_exit_with_message( "cutpoint_closed == 0" );


			std::cout << "Applying cutpoint variants to " << cutpoint_closed << std::endl;

			Correctly_position_cutpoint_phosphate_torsions( full_pose, cutpoint_closed, verbose /*verbose*/ );

			pose::add_variant_type_to_pose_residue( full_pose, chemical::CUTPOINT_LOWER, cutpoint_closed   );
			pose::add_variant_type_to_pose_residue( full_pose, chemical::CUTPOINT_UPPER, cutpoint_closed+1 );

			Add_harmonic_chainbreak_constraint(full_pose, cutpoint_closed );

		}
		dump_pdb( full_pose, "full_pose_after_add_cutpoint_variant.pdb");	

		//std::cout << "Dec 26 Check point 1 "<< std::endl;	

		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.00000025);
		bool const use_nblist( true );
		MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
		options.nblist_auto_update( true );


		///////////////////////////////////////////////////Add cst to fixed res pos./////////////////////////////////////////////////////////////////////////

		Real const coord_sdev( 2.0 );
		Size const my_anchor( 1 ); //anchor atom on first residue? Basically this residue should be fixed in space?

		ConstraintSetOP cst_set = full_pose.constraint_set()->clone();

		for ( Size i=1; i<= full_pose.total_residue();  i++ ) {
		
			if(Contain_seq_num(i, minimize_res_list)) continue;

			Residue const & i_rsd( full_pose.residue(i) );

			for ( Size ii = 1; ii<= i_rsd.nheavyatoms(); ++ii ) {

				cst_set->add_constraint( new CoordinateConstraint( AtomID(ii,i), AtomID(1,my_anchor), i_rsd.xyz(ii), new HarmonicFunc( 0.0, coord_sdev ) ) );
			}
		}

		full_pose.constraint_set( cst_set );

		///////////////////////////////////////////////////Minimize/////////////////////////////////////////////////////////////////////////
		//Size const total_rounds=cutpoint_closed_list.size(); BUG!! fix on Dec 28, 2010
		Size const total_rounds=scorefxn_list.size();

		bool SKIP_POSE_DUE_TO_ABNORMALLY_BAD_ENERGY=false;

		for(Size round_ID=1; round_ID<=total_rounds; round_ID++){
			if(SKIP_POSE_DUE_TO_ABNORMALLY_BAD_ENERGY==true) continue; //Need this to prevent minimization error message.
		
			core::kinematics::MoveMap const mm=mm_list[round_ID];

			if(skip_minimize==false){

				//if( (scorefxn_list[round_ID]->get_weight(coordinate_constraint)) > 0.0001){

					//if(verbose) std::cout << "scorefxn_list[" << round_ID << "] coordinate_constraint >0.0001, do_minimize_with_constraints!" << std::endl;
				
					//minimize_with_constraints(full_pose, mm, (scorefxn_list[round_ID]), options);
				//}else{
					minimizer.run( full_pose, mm, *(scorefxn_list[round_ID]), options );		
				//}

				o2star_minimize(full_pose, scorefxn_list[round_ID], get_surrounding_O2star_hydrogen(full_pose, minimize_res_list, verbose) );

				for(Size cc_ID=1; cc_ID<=cutpoint_closed_list.size(); cc_ID++){

					Size const cutpoint_closed=cutpoint_closed_list[cc_ID]; 	
		
					Real mean_dist_err = rna_loop_closer.apply( full_pose, cutpoint_closed );
					std::cout << "mean_dist_err for cutpoint_closed ("<< cutpoint_closed <<", round_ID= " << round_ID << " ) = " <<  mean_dist_err << std::endl;
				}

				o2star_minimize(full_pose, scorefxn_list[round_ID], get_surrounding_O2star_hydrogen(full_pose, minimize_res_list, verbose) ); 
				minimizer.run( full_pose, mm, *(scorefxn_list[round_ID]), options );

			}			

			std::string full_pose_tag_MOD=full_pose_tag;

			if(round_ID==scorefxn_list.size()){
				full_pose_tag_MOD[0]='M';
			}else{
				full_pose_tag_MOD.erase(0,1);
				full_pose_tag_MOD=string_of(round_ID)+full_pose_tag_MOD;
			}

			//if(verbose){
				std::cout << "full_pose_tag_MOD= " << full_pose_tag_MOD << std::endl;
				scorefxn_list[round_ID]->show( std::cout, full_pose );
				Real const current_score=(*scorefxn_list[round_ID])(full_pose);
				std::cout << "current_score= " << current_score << std::endl;
				if(current_score>10000.00){
					std::cout << "SKIP_POSE_DUE_TO_ABNORMALLY_BAD_ENERGY: " << full_pose_tag << std::endl;
					SKIP_POSE_DUE_TO_ABNORMALLY_BAD_ENERGY=true;
				}
			//}


			(*scorefxn_list[total_rounds])(full_pose); //Ouput the pose with the final round score function.

			if(output_silentfile_every_round || (round_ID==total_rounds) || (round_ID==1) || (round_ID==2) ){

				std::string act_silent_outfile=silent_outfile;
				if(round_ID!=total_rounds){
					act_silent_outfile+="_round_ID_" + string_of(round_ID);
				}

				BinaryRNASilentStruct s( full_pose, full_pose_tag_MOD );

				if(verbose){
						protocols::swa::rna::Output_seq_num_list("native_virtual_res_list=", native_virtual_res_list, TR );
						protocols::swa::rna::Output_seq_num_list("native_alignment_res_list=",native_alignment_res_list, TR );
						protocols::swa::rna::Output_seq_num_list("rmsd_res_list=",rmsd_res_list, TR );
				}


				pose::Pose native_pose=native_pose_ACT; //HARD COPY

				for(Size n=1; n<=native_virtual_res_list.size(); n++){
					protocols::swa::rna::apply_virtual_rna_residue_variant_type(native_pose, native_virtual_res_list[n], true /*apply_check*/);
				}

				//std::map< core::Size, core::Size > full_to_sub;
				//full_to_sub.clear();	
	
				if(full_pose.total_residue()!=native_pose.total_residue()){
					std::cout << "full_pose.total_residue()= " << full_pose.total_residue() << std::endl;
					std::cout << "native_pose.total_residue()= " << native_pose.total_residue() << std::endl;
					utility_exit_with_message( "full_pose.total_residue()!=native_pose.total_residue()" );
				}

				//for(Size n=1; n<=full_pose.total_residue(); n++){
				//	full_to_sub[n]=n;
				//}

				protocols::swa::rna::align_poses(native_pose, "native", full_pose, full_pose_tag_MOD, native_alignment_res_list);

				s.add_energy( "Full_L_rmsd", full_length_rmsd_over_residue_list(full_pose, native_pose, rmsd_res_list, native_pose.sequence(), verbose, false) );
				s.add_energy( "Full_V_L_rms", full_length_rmsd_over_residue_list(full_pose, native_pose, rmsd_res_list, native_pose.sequence(), verbose, true) );
				s.add_energy( "COMBINE_SCORE" , FOUND_combine_tag_info.combine_score);
				s.add_energy( "OVERLAP_RMSD" , FOUND_combine_tag_info.overlap_rmsd);

				silent_file_data.write_silent_struct(s, act_silent_outfile, false); 

			}

		}


		verbose=false;
	}
	}

	std::cout << "filtered_tag_pair_info_list.size()= " << filtered_tag_pair_info_list.size() << " num_tag_pair_found= " << num_tag_pair_found << std::endl;

	if(filtered_tag_pair_info_list.size()!=num_tag_pair_found){
		utility_exit_with_message( "filtered_tag_pair_info_list.size()!=num_tag_pair_found!" );
	}

	std::cout << "Total time in ROSETTA hermann_phase_two_minimize(): " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;
	std::cout << "ROSETTA hermann_phase_two_minimize EXIT" << std::endl;

}

void
hermann_phase_two(){

	using namespace core::pose;
	using namespace ObjexxFCL;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::scoring;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );


	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if ( !option[ in::file::native ].user() ) utility_exit_with_message( "User must supply in::file::native!" ); //This is the static pose.
	if ( !option[ alignment_res_pairs ].user() ) utility_exit_with_message( "User must supply alignment_res_pairs!" );
	if ( !option[ alignment_res_pairs_two ].user() ) utility_exit_with_message( "User must supply alignment_res_pairs_two!" );
	if ( !option[ alignment_RMSD_CUTOFF ].user() ) utility_exit_with_message( "User must supply alignment_RMSD_CUTOFF!" );

	if ( !option[ RMSD_res_pairs ].user() ) utility_exit_with_message( "User must supply RMSD_res_pairs!" );



	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > const RMSD_res_pair_list = option[ RMSD_res_pairs ]();	

	utility::vector1< core::Size > lower_pdb_RMSD_res;
	utility::vector1< core::Size > upper_pdb_RMSD_res;
	
	lower_pdb_RMSD_res.clear();
	upper_pdb_RMSD_res.clear();

	for(Size n=1; n<=RMSD_res_pair_list.size(); n++){

		utility::vector1< std::string > const RMSD_res_pair=Tokenize(RMSD_res_pair_list[n], "-");
		if(RMSD_res_pair.size()!=2){
			 utility_exit_with_message( "RMSD_res_pair.size()!=2, RMSD_res_pair_list[n]= " + RMSD_res_pair_list[n] );
		}

		lower_pdb_RMSD_res.push_back(string_to_int(RMSD_res_pair[1]));
		upper_pdb_RMSD_res.push_back(string_to_int(RMSD_res_pair[2]));
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	if(lower_pdb_RMSD_res.size()!=upper_pdb_RMSD_res.size()){
		utility_exit_with_message( "lower_pdb_RMSD_res.size()!=upper_pdb_RMSD_res.size()" );
	}

	std::cout << "lower_pdb_RMSD_res to upper_pdb_RMSD_res:" << std::endl;
	for(Size ii=1; ii<=lower_pdb_RMSD_res.size(); ii++){
		std::cout << lower_pdb_RMSD_res[ii] << " ---> " << upper_pdb_RMSD_res[ii] << std::endl;
	}


	std::string const static_pdb_tag = option[ in::file::native ]();

	utility::vector1< std::string > const lower_alignment_res_pair_list = option[ alignment_res_pairs ]();	
	utility::vector1< std::string > const upper_alignment_res_pair_list = option[ alignment_res_pairs_two ]();	

	Real const alignment_RMSD_cutoff=option[ alignment_RMSD_CUTOFF ]();	 


	pose::Pose static_pose;
	import_pose::pose_from_pdb( static_pose, *rsd_set, static_pdb_tag );


	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< std::string > silent_files_list = option[ in::file::silent ]();

	if(silent_files_list.size()!=2) utility_exit_with_message( "silent_files_list.size()!=2" );
	
	std::string const upper_silent_file=silent_files_list[2];


	////////////////////////////////////////////////////////////////////////////////////////////////////////////

	utility::vector1< pose_data_struct2 > lower_pose_data_list=convert_silent_file_to_pose_data_list(silent_files_list[1]);
	utility::vector1< pose_data_struct2 > upper_pose_data_list=convert_silent_file_to_pose_data_list(silent_files_list[2]);


	Output_title_text("align LOWER SWA pdbs", TR );
	align_pdbs_function(static_pose, lower_pose_data_list, lower_alignment_res_pair_list, 	alignment_RMSD_cutoff);

	Output_title_text("align UPPER SWA pdbs", TR );
	align_pdbs_function(static_pose, upper_pose_data_list, upper_alignment_res_pair_list, 	alignment_RMSD_cutoff);


	std::string const align_lower_SWA_folder="align_lower_SWA/";
	std::string const align_upper_SWA_folder="align_upper_SWA/";

	for(Size n=1; n<=lower_pose_data_list.size(); n++){
		pose_data_struct2 & pose_data=lower_pose_data_list[n];
		pose_data.tag= "aligned_lower_" +pose_data.tag;
	
		dump_pdb( *(pose_data.pose_OP), align_lower_SWA_folder + pose_data.tag + ".pdb");	
	
	}

	for(Size n=1; n<=upper_pose_data_list.size(); n++){
		pose_data_struct2 & pose_data=upper_pose_data_list[n];
		pose_data.tag= "aligned_upper_" +pose_data.tag;	

		dump_pdb( *(pose_data.pose_OP), align_upper_SWA_folder + pose_data.tag + ".pdb");	
	
	}

	/////////////////////////////////////////////////////////////////////////////////////


	std::string outfile_name="overlap_base_rmsd.txt";

	std::ofstream outfile;
	outfile.open(outfile_name.c_str()); //Opening the file with this command removes all prior content..

	outfile << std::setw(50) << std::left << "lower_tag" ; 
	outfile << std::setw(50) << std::left << "upper_tag" ; 
	outfile << std::setw(15) << std::left << "sum_score" ;
	outfile << std::setw(15) << std::left << "base_rmsd";
	outfile << "\n"; 
	////////////////////////Calculate the overlap helix RMSD///////////////////////////////////


	for(Size lower_ID=1; lower_ID<=lower_pose_data_list.size(); lower_ID++){	
	for(Size upper_ID=1; upper_ID<=upper_pose_data_list.size(); upper_ID++){

		pose_data_struct2 & lower_pose_data=lower_pose_data_list[lower_ID];
		pose_data_struct2 & upper_pose_data=upper_pose_data_list[upper_ID];

		Size total_atom_count=0;
		Real total_sum_sd=0.0;

		for(Size ii=1; ii<=lower_pdb_RMSD_res.size(); ii++){

			Size const lower_seq_num=lower_pdb_RMSD_res[ii];
			Size const upper_seq_num=upper_pdb_RMSD_res[ii];

			//bool verbose= (ii==1) ? true : false;
			bool verbose = false;

			Size atom_count=0;
			Real sum_sd=0.0;

			base_atoms_square_deviation( *(lower_pose_data.pose_OP), *(upper_pose_data.pose_OP), lower_seq_num, upper_seq_num, atom_count, sum_sd, verbose , false /*ignore_virtual_atom*/);


			if(atom_count==0){
				std::cout << "atom_count=0 for lower_seq_num= " << lower_seq_num << " and upper_seq_num= " << upper_seq_num << std::endl;
				std::cout << "lower_pdb_tag= " << lower_pose_data.tag <<  "upper_pdb_tag= " << upper_pose_data.tag << std::endl;
			}


			total_atom_count+=atom_count;
			total_sum_sd+=sum_sd;

		}


		total_sum_sd=total_sum_sd/(total_atom_count);
		Real const all_base_rmsd=sqrt(total_sum_sd);

		Real const sum_score = lower_pose_data.score + upper_pose_data.score;

		//std::string output_tag= lower_pose_data.tag + "_" + upper_pose_data.tag;

		outfile << std::setw(50) << std::left << lower_pose_data.tag; 
		outfile << std::setw(50) << std::left << upper_pose_data.tag; 
		outfile << std::setw(15) << std::left << sum_score ;
		outfile << std::setw(15) << std::left << all_base_rmsd;
		outfile << "\n"; 


	}
	}

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void
journal_club_syn_chi(){

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::kinematics;

		ResidueTypeSetCAP rsd_set;
		rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

		std::string sequence="g";
		pose::Pose pose;
		make_pose_from_sequence( pose, sequence, *rsd_set, false /*auto_termini*/);

		pose.set_torsion( TorsionID( 1, CHI, 1 ) , 79.430);	

		dump_pdb( pose, "anti_chi.pdb" );	

		pose.set_torsion( TorsionID( 1, CHI, 1 ) , -51.00);	

		dump_pdb( pose, "syn_chi_51.pdb" );

		pose.set_torsion( TorsionID( 1, CHI, 1 ) , -71.00);	

		dump_pdb( pose, "syn_chi_71.pdb" );

		pose.set_torsion( TorsionID( 1, CHI, 1 ) , -31.00);	

		dump_pdb( pose, "syn_chi_31.pdb" );

}


/*usage*/
// ~/src/mini/bin/parin_test.macosgccrelease -algorithm calculate_theoretical_RNA_length -database ~/minirosetta_database > output.txt


void
calculate_theoretical_RNA_length(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::kinematics;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	//Create a trinucleotide
	std::string sequence="auu";
	pose::Pose pose;
	make_pose_from_sequence( pose, sequence, *rsd_set, false /*auto_termini*/);

	protocols::viewer::add_conformation_viewer( pose.conformation(), "test", 400, 400 );
						
	std::string start_atom=" O3'"; //" C5'", " C1'";
	core::Size start_res=2;
	std::string end_atom=" O3'"; //" O3'" , "P";
	core::Size end_res=1;

	Size const bin_size=5;
//	Size num_bin=360/bin_size;

	Real max_dist=0;
	Real mean_dist=0;	
	Size count=0;

	Size const min_angle=0;	
	Size const max_angle=360;
	//prepend
//	for(Size a3=min_angle; a3<=max_angle; a3+=bin_size){
//	for(Size z2=min_angle; z2<=max_angle; z2+=bin_size){
//	for(Size e2=min_angle; e2<=max_angle; e2+=bin_size){
	for(Size g2=min_angle; g2<=max_angle; g2+=bin_size){
	for(Size b2=min_angle; b2<=max_angle; b2+=bin_size){
	for(Size a2=min_angle; a2<=max_angle; a2+=bin_size){


//		pose.set_torsion( TorsionID( 3, BB, 1 ) , a3);		
//		pose.set_torsion( TorsionID( 2, BB, 6 ) , z2);			
//		pose.set_torsion( TorsionID( 2, BB, 5 ) , e2);		
		pose.set_torsion( TorsionID( 2, BB, 3 ) , g2);		
		pose.set_torsion( TorsionID( 2, BB, 2 ) , b2);			
		pose.set_torsion( TorsionID( 2, BB, 1 ) , a2);		

		for(Size d2=1; d2<=2; d2++){
			Real d2_angle = (d2==1) ? 85.418:152.467;  
			pose.set_torsion( TorsionID( 2, BB, 4 ) , d2_angle);


			Real distance=(pose.residue(start_res).xyz(start_atom)- pose.residue(end_res).xyz(end_atom)).length();
			if(max_dist<distance) max_dist=distance;

			mean_dist+=distance;
			count++;

			std::cout << "distance= " << distance;

//			std::cout << " a3_angle= " << a3 << " z2_angle= " << z2 << " e2_angle= " << e2;
			std::cout << " g2_angle= " << g2 << " b2_angle= " << b2 << " a2_angle= " << a2; 
			std::cout << " d2= " << d2_angle << std::endl;
		}
	}
	}
	}
//	}
//	}
//	}

	mean_dist=(mean_dist/count);

	std::cout << "max_distance= " << max_dist;
	std::cout << " mean_dist= " << mean_dist;

}

/*
void
calculate_bulge_length_distribution(){

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::kinematics;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	//Create a trinucleotide
	std::string sequence="auu";
	pose::Pose pose;
	make_pose_from_sequence( pose, sequence, *rsd_set, false );


	protocols::viewer::add_conformation_viewer( pose.conformation(), "test", 400, 400 );

	std::string start_atom_1=" C5'"; //" C1'";
	std::string start_atom_2=" C4'"; //" C1'";
	core::Size start_res=3;

	std::string end_atom_1=" O3'"; //"P";
	std::string end_atom_2=" C3'"; //"P";
	core::Size end_res=1;


	Size const bin_size=10;

	Size const min_angle=0;	
	Size const max_angle=360;




//	for(Size g3=min_angle; g3<=max_angle; g3+=bin_size){

	for(Size b3=min_angle; b3<=max_angle; b3+=bin_size){
	for(Size a3=min_angle; a3<=max_angle; a3+=bin_size){
	for(Size z2=min_angle; z2<=max_angle; z2+=bin_size){

	std::cout << "setup_delta_rotamer_generator" << std::endl;
	BaseState const bulge_base_state = (is_purine( (*pose_data_list[1].pose_OP).residue( FB_job_params.bulge_res ) ) ) ? BOTH: ANTI;
	StepWiseRNA_BaseSugarRotamerOP bulge_base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( bulge_base_state, ALL, rna_fitted_torsion_info);


//for(Size e2=min_angle; e2<=max_angle; e2+=bin_size){

	for(Size g2=min_angle; g2<=max_angle; g2+=bin_size){
	for(Size b2=min_angle; b2<=max_angle; b2+=bin_size){
	for(Size a2=min_angle; a2<=max_angle; a2+=bin_size){
//	for(Size z1=min_angle; z1<=max_angle; z1+=bin_size){


}
*/

void
calculate_theoretical_RNA_length_with_bond_angle_dependence(){

		using namespace core::chemical;
		using namespace core::conformation;
		using namespace core::scoring;
		using namespace core::pose;
		using namespace core::io::silent;
		using namespace core::id;
		using namespace core::kinematics;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	//Create a trinucleotide
	std::string sequence="auu";
	pose::Pose pose;
	make_pose_from_sequence( pose, sequence, *rsd_set, false /*auto_termini*/);

	protocols::viewer::add_conformation_viewer( pose.conformation(), "test", 400, 400 );

	std::string start_atom_1=" C5'"; //" C1'";
	std::string start_atom_2=" C4'"; //" C1'";
	core::Size start_res=3;

	std::string end_atom_1=" O3'"; //"P";
	std::string end_atom_2=" C3'"; //"P";
	core::Size end_res=2;

	Size bin_size=10;
	
	Size min_angle=0;	
	Size max_angle=360;

	std::ofstream outfile;
	outfile.open("distance.txt");

	Size const spacing=8;

	Real max_dist=0;
	Real min_dist=999999;

	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "b3";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "a3";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "z2";
//	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "e2";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "C5_O3";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "C4_C3";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "dot";
	outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << "dot_c";
	outfile << "\n";

	for(Size b3=min_angle; b3<=max_angle; b3+=5){
	for(Size a3=min_angle; a3<=max_angle; a3+=5){
	for(Size z2=min_angle; z2<=max_angle; z2+=5){
//	for(Size e2=min_angle; e2<=max_angle; e2+=10){

		pose.set_torsion( TorsionID( 3, BB, 2 ) , b3);		
		pose.set_torsion( TorsionID( 3, BB, 1 ) , a3);			
		pose.set_torsion( TorsionID( 2, BB, 6 ) , z2);		
//		pose.set_torsion( TorsionID( 2, BB, 5 ) , e2);		

			Real C5_O3_distance=(pose.residue(start_res).xyz(start_atom_1)- pose.residue(end_res).xyz(end_atom_1)).length();
			Real C4_C3_distance=(pose.residue(start_res).xyz(start_atom_2)- pose.residue(end_res).xyz(end_atom_2)).length();

			if(max_dist<C4_C3_distance) max_dist=C4_C3_distance;
			if(min_dist>C4_C3_distance) min_dist=C4_C3_distance;

			numeric::xyzVector<Real> start_vector=pose.residue(start_res).xyz(start_atom_2)-pose.residue(start_res).xyz(start_atom_1);
			numeric::xyzVector<Real> end_vector=pose.residue(end_res).xyz(end_atom_1)-pose.residue(end_res).xyz(end_atom_2);

			start_vector.normalize();
			end_vector.normalize();

			Real dot_product=dot( start_vector, end_vector);
			Real dot_product_check=dot( end_vector, start_vector);

			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << b3;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << a3;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << z2;
//			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << e2;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << C5_O3_distance;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << C4_C3_distance;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << dot_product;
			outfile << std::setw(spacing) << std::fixed << std::setprecision(3)  << std::left << dot_product_check;



			outfile << "\n";

//	}
	}
	}
	}

	outfile.close();

	std::cout << "max_distance= " << max_dist << " min_distance= " << min_dist;


/*
  delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=   85.418        chi_2=   59.430         nu2_2=   37.364         nu1_2=   92.204
        delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=   85.418        chi_2=   79.430         nu2_2=   37.364         nu1_2=   92.204
        delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=   85.418        chi_2=   99.430         nu2_2=   37.364         nu1_2=   92.204
        delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=  152.467        chi_2=   96.600         nu2_2=  -38.085         nu1_2=  150.133
        delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=  152.467        chi_2=  116.600         nu2_2=  -38.085         nu1_2=  150.133
        delta1=   85.418        chi_1=   59.430         nu2_1=   37.364         nu1_1=   92.204         delta2=  152.467        chi_2=  136.600         nu2_2=  -38.085         nu1_2=  150.133


			add_torsion_id(TorsionID( moving_suite_    , BB, 4 ));  //delta1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 1 )); //chi_1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 2 )); //nu2_1
			add_torsion_id(TorsionID( moving_suite_    , CHI, 3 )); //nu2_2

			lower_base_state = (sample_syn_chi1_) ? BOTH: ANTI;
			if(Is_bulge_) lower_base_state= NONE ;
		}

		add_torsion_id( TorsionID( moving_suite_    , BB, 5 ) ); // epsilon1
		add_torsion_id( TorsionID( moving_suite_    , BB, 6 ) ); // zeta1
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 1 ) ); // alpha2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 2 ) ); // beta2
		add_torsion_id( TorsionID( moving_suite_ + 1, BB, 3 ) ); // gamma2

			add_torsion_id( TorsionID( moving_suite_ + 1, BB, 4 ) );  //delta1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 1 ) ); //chi_1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 2 ) ); //nu2_1
			add_torsion_id( TorsionID( moving_suite_ + 1, CHI, 3 ) ); //nu2_2

			upper_base_state = (sample_syn_chi2_) ? BOTH: ANTI;
			if(Is_bulge_) upper_base_state = NONE ;
		}
*/

}


/////////////////////////////////////////////////////////////////////

/* In PROGRESS!
void
extract_clash_list(){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::swa::rna;
	using namespace protocols::rna;
	using namespace scoring::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::id; 

	clock_t const time_start( clock() ); 	

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "Must supply in::file::s!" );

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	std::string const pdb_file=   option[ in::file::s ]()[1];	

	std::cout << "importing " << pdb_file << std::endl;
	pose::Pose pose;

	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );


	utility::vector1< core::Size > sample_res_list= option[ sample_res ]();

	sort_seq_num_list(sample_res_list);

	Output_seq_num_list("sample_res_list= ", sample_res_list, TR );

	for(Size seq_num_1=1; seq_num_1<=sample_res_list.size(); seq_num_1++){

	for(Size seq_num_2=1; seq_num_2<=pose.total_residue(); seq_num_2++){

	
		conformation::Residue const & sample_rsd= pose.residue(sample_res_ID);

		conformation::Residue const & other_rsd= pose.residue(other_res_ID);

			for(Size atomno_1=1; atomno_1<=sample_rsd.natoms(); atomno_1++){ //include hydrogen atoms

			std::string const & atom_name_1=rsd_1.type().atom_name(atomno_1);
			if(atom_name_1=="HO2'") continue;
			if(rsd_1.atom_type(atomno_1).name()=="VIRT") continue; 


			bool const Is_3_prime_phosphate_1=Is_three_prime_phosphate_atom(atom_name_1); //This is just the O3' atom
			bool const Is_5_prime_phosphate_1=Is_five_prime_phosphate_atom(atom_name_1);  //This is O5', OP2, OP1 and P 

			bool const Is_base_sugar_1=(Is_3_prime_phosphate_1==false && Is_5_prime_phosphate_1==false);

			for(Size atomno_2=1; atomno_2<=other_rsd.natoms(); atomno_2++){ //include hydrogen atoms

				if(seq_num_1==seq_num_2){
					if(atomno_1>=atomno_2) continue; 
				}

				std::string const & atom_name_2=rsd_2.type().atom_name(atomno_2);
				if(atom_name_2=="HO2'") continue;
				//if(rsd_2.atom_type(atomno_2).name()=="VIRT") continue; 

				bool const Is_3_prime_phosphate_2=Is_three_prime_phosphate_atom(atom_name_2);
				bool const Is_5_prime_phosphate_2=Is_five_prime_phosphate_atom(atom_name_2);

				bool const Is_base_sugar_2=(Is_3_prime_phosphate_2==false && Is_5_prime_phosphate_2==false);

				if(seq_num_1==seq_num_2){ //same nucleotide
					if( (Is_base_sugar_1 ) && (Is_base_sugar_2 ) ) continue; //Ignore clash between same base-sugar group

					if( Is_3_prime_phosphate_1 && Is_3_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group
					if( Is_5_prime_phosphate_1 && Is_5_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group

					if( Is_3_prime_phosphate_1 && Is_base_sugar_2) continue; //O3' atom is actually part of the 3' sugar
					if( Is_3_prime_phosphate_2 && Is_base_sugar_1) continue; //O3' atom is actually part of the 3' sugar

				}

				if(seq_num_1==1 && seq_num_2==2){ 
					if( Is_3_prime_phosphate_1 && Is_5_prime_phosphate_2) continue; //Ignore clash between same phosphate group
				} 

				if(seq_num_2==1 && seq_num_1==2){ 
					if( Is_3_prime_phosphate_2 && Is_5_prime_phosphate_1) continue; //Ignore clash between same phosphate group
				}

				if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_1, seq_num_1, atom_name_2, seq_num_2)) continue;
				if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_2, seq_num_2, atom_name_1, seq_num_1)) continue;

				if(atom_name_1==atom_name_2 && seq_num_1==seq_num_2){
					std::cout << "atom_name_1=" << atom_name_1 << std::endl;
					std::cout << "seq_num_1="   << seq_num_1   << std::endl;
					std::cout << "atom_name_2=" << atom_name_2 << std::endl;
					std::cout << "seq_num_2="   << seq_num_2   << std::endl;
					utility_exit_with_message("atom_name_1==atom_name_2 && seq_num_1==seq_num_2");
				}

				res_2_atom_ID_list.push_back(std::make_pair(seq_num_2, atomno_2));

			}//atomno_2
		}//seq_num_2

		VDW_rep_atom_map_list.push_back(std::make_pair(res_1_atom_ID, res_2_atom_ID_list));
	}//atomno_1
	}//seq_num_1



	}



}
*/
/////////////////////////////////////////////////////////////////////


void
extract_hydrogen_bonds_statistic(){

	utility_exit_with_message("Source code required to run this function have not yet been committed to TRUNK!");

/*
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::swa::rna;
	using namespace protocols::rna;
	using namespace scoring::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::id; 


	clock_t const time_start( clock() ); 	

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "Must supply in::file::s!" );

	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must supply sample_res!" );


	bool const include_edge_phosphate= option[ INCLUDE_EDGE_PHOSPHATE ]();

	bool const double_count_BP= option[ double_count_base_pair]();

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	std::string const pdb_file=   option[ in::file::s ]()[1];	

	std::cout << "importing " << pdb_file << std::endl;
	pose::Pose pose;

	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );


	utility::vector1< core::Size > input_sample_res_list= option[ sample_res ]();

	sort_seq_num_list(input_sample_res_list);

	utility::vector1< core::Size > sample_res_list=input_sample_res_list;

	Output_boolean("include_edge_phosphate= ", include_edge_phosphate, TR ); std::cout << std::endl;


	if(include_edge_phosphate){ //Assume continuous_range for now.
		Size const min_res=sample_res_list[1];
		Size const max_res=sample_res_list[sample_res_list.size()];

		if(min_res<=1){
			utility_exit_with_message("min_res=" + string_of(min_res) + "<=1!");
		}

		if(max_res>=pose.total_residue()){
			utility_exit_with_message("max_res=" + string_of( max_res ) + "<=" + string_of( pose.total_residue() ) + "=pose.total_residue()!");
		}

		sample_res_list.push_back(min_res-1);
		sample_res_list.push_back(max_res+1);

		sort_seq_num_list(sample_res_list);

	}

	Output_seq_num_list("input_sample_res_list= ", input_sample_res_list, TR );
	Output_seq_num_list("sample_res_list= ", sample_res_list, TR );

	utility::vector1< Hydrogen_Bond_Info > hydrogen_bond_info_list;

	bool verbose=true;
	
	for(Size sample_res_ID=1; sample_res_ID<=sample_res_list.size(); sample_res_ID++){

		bool Is_five_prime_edge=false;
		bool Is_three_prime_edge=false;


		if(include_edge_phosphate){
			if(sample_res_ID==1 && include_edge_phosphate==true) Is_five_prime_edge=true;
			if(sample_res_ID==sample_res_list.size() && include_edge_phosphate==true) Is_three_prime_edge=true;
		}

		core::Size const sample_res=sample_res_list[sample_res_ID];

		core::conformation::Residue const & sample_rsd = pose.residue(sample_res);

		for(Size seq_num=1; seq_num<=pose.total_residue(); seq_num++){

			core::conformation::Residue const & surrounding_rsd=pose.residue(seq_num);

			get_hbonds_1way(hydrogen_bond_info_list, sample_rsd, surrounding_rsd, sample_res, Is_five_prime_edge, Is_three_prime_edge, verbose); //sample_rsd as donor

			get_hbonds_1way(hydrogen_bond_info_list, surrounding_rsd, sample_rsd, sample_res, Is_five_prime_edge, Is_three_prime_edge, verbose); //sample_rsd as acceptor

			if(verbose==true) verbose=false;

		}

	}


	std::cout << "--------------------------Classify Base Pairs (RHIJU METHOD)-------------------------- " << std::endl;

	utility::vector1< core::scoring::rna::Base_pair> const RHIJU_base_pair_list = classify_base_pairs_parin(pose, input_sample_res_list);

	for(Size n=1; n<=RHIJU_base_pair_list.size(); n++){

		std::cout << "base_pair # " << n << " : "; RHIJU_base_pair_list[n].print_info(); std::cout << " (RHIJU METHOD) " << std::endl;
	}

	std::cout << "--------------------------NUM ClASSIFED BASE PAIR (RHIJU METHOD)= " << RHIJU_base_pair_list.size()<< "--------------------------" << std::endl; 

	std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
	std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------ " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------ " << std::endl;
	std::cout << "------------------------------------------------------------------------------------------------------------------------------------------------------------ " << std::endl;
	


	std::cout << "--------------------------Classify Base Pairs-------------------------- " << std::endl;

	utility::vector1< core::scoring::rna::Base_pair> base_pair_list=classify_base_pairs_strict(pose, input_sample_res_list, hydrogen_bond_info_list, false);

	

	if(double_count_BP==false){
		utility::vector1< core::scoring::rna::Base_pair> double_count_base_pair_list=base_pair_list;
		base_pair_list.clear();

		for(Size n=1; n<=double_count_base_pair_list.size(); n++){
			core::scoring::rna::Base_pair const base_pair=double_count_base_pair_list[n];
		
			if(base_pair.res1==base_pair.res2){
				std::cout << "base_pair.res1=" << base_pair.res1 << "base_pair.res2=" << base_pair.res2 << std::endl;
				utility_exit_with_message("base_pair.res1==base_pair.res2");
			}	

			if(base_pair.res1>base_pair.res2) continue;
			base_pair_list.push_back(base_pair);

		}
	}


	std::string const outfile_name="base_pair_list.txt";

	std::ofstream outfile;
	outfile.open(outfile_name.c_str()); //Opening the file with this command removes all prior content..

	for(Size n=1; n<=base_pair_list.size(); n++){
		std::cout << "base_pair # " << n << " : "; base_pair_list[n].print_info(); std::cout << std::endl;
		outfile << "base_pair # " << n << " : "; base_pair_list[n].print_info( outfile ); outfile << "\n"; 
	}

	outfile.close();


	Size num_classify_WC_BP=0;
	Size num_classify_non_WC_BP=0;

	for ( Size n = 1; n <= base_pair_list.size(); n++ ) {

		core::scoring::rna::Base_pair const base_pair = base_pair_list[ n ];

		Residue const & rsd_1= pose.residue( base_pair.res1 );
		Residue const & rsd_2= pose.residue( base_pair.res2 );

		if ( ( base_pair.edge1 == WATSON_CRICK && base_pair.edge2 == WATSON_CRICK && base_pair.orientation == 1 )  && core::scoring::rna::possibly_canonical( rsd_1.aa(), rsd_2.aa() ) ){
			num_classify_WC_BP++;
		} else {
			num_classify_non_WC_BP++;
		}
	}

	Size const total_classified_BP=base_pair_list.size();


	std::cout << "--------------------------NUM ClASSIFED BASE PAIR= " << total_classified_BP << "--------------------------" << std::endl; std::cout << std::endl;

	
	std::cout << "--------------------------Potential Hydrogen Bonds-------------------------- " << std::endl;

	Hydrogen_Bond_Count hbond_count;

	for(Size n=1; n<=hydrogen_bond_info_list.size(); n++){

		Hydrogen_Bond_Info & hbond=hydrogen_bond_info_list[n];

		std::string const hbond_type=get_hbond_type(hbond);

		if(hbond_type=="BASE-BASE" || hbond_type=="BASE-BASE" ) hbond_count.base_base++;
		
		if(hbond_type=="BASE-2*OH" || hbond_type=="2*OH-BASE" ) hbond_count.base_O2star++;

		if(hbond_type=="BASE-O4'*" || hbond_type=="O4'*-BASE" ) hbond_count.base_O4star++;

		if(hbond_type=="BASE-PHOS" || hbond_type=="PHOS-BASE" ) hbond_count.base_phos++;

		if(hbond_type=="2*OH-2*OH" || hbond_type=="2*OH-2*OH" ) hbond_count.O2star_O2star++;

		if(hbond_type=="2*OH-O4'*" || hbond_type=="O4'*-2*OH" ) hbond_count.O2star_O4star++;

		if(hbond_type=="2*OH-PHOS" || hbond_type=="PHOS-2*OH" ) hbond_count.O2star_phos++;


		if(std::abs(int(hbond.acceptor_res-hbond.donor_res))==1) hbond.comment+="I and I+1 case, ";
		
		std::cout << "#" << std::setw(5) << std::left << n << " "; 

		print_hydrogen_bond_info(hbond);
	}



	std::cout << "--------------------------Total_hbonds= " << hydrogen_bond_info_list.size() << "--------------------------" << std::endl;

	if(( hbond_count.base_base			+	hbond_count.base_O2star 		+ hbond_count.base_O4star		+ hbond_count.base_phos+
			 hbond_count.O2star_O2star	+	hbond_count.O2star_O4star  +	hbond_count.O2star_phos )!= hydrogen_bond_info_list.size()){

		utility_exit_with_message("hbond_count INCONSISTENCY!");

	}


	std::cout << "hbond_count.base_base= "     << hbond_count.base_base << std::endl;
	std::cout << "hbond_count.base_O2star= "   << hbond_count.base_O2star << std::endl;
	std::cout << "hbond_count.base_O4star= "   << hbond_count.base_O4star << std::endl;
	std::cout << "hbond_count.base_phos= "     << hbond_count.base_phos << std::endl;

	std::cout << "hbond_count.O2star_O2star= " << hbond_count.O2star_O2star << std::endl;
	std::cout << "hbond_count.O2star_O4star= " << hbond_count.O2star_O4star << std::endl;
	std::cout << "hbond_count.O2star_phos= "   << hbond_count.O2star_phos << std::endl;



	Size const total_sample_res=input_sample_res_list.size();

	std::cout << std::endl;
	std::cout << std::setw(40) << "total_sample_res= " << total_sample_res << std::endl;
	std::cout << std::setw(40) << "total_classified_BP= " << std::setw(5) << (total_classified_BP) << " per nucleotide= " << std::setw(6) << float(total_classified_BP)/float(total_sample_res) << std::endl;
	std::cout << std::setw(40) << "total_classified_WC_BP= " << std::setw(5) << (num_classify_WC_BP) << " per nucleotide= " << std::setw(6) << float(num_classify_WC_BP)/float(total_sample_res) << std::endl;
	std::cout << std::setw(40) << "total_classified_non_WC_BP= " << std::setw(5) << (num_classify_non_WC_BP) << " per nucleotide= " << std::setw(6) << float(num_classify_non_WC_BP)/float(total_sample_res) << std::endl;
	std::cout << std::setw(40) << "total_hbonds= " << std::setw(5) << hydrogen_bond_info_list.size()  << " per nucleotide= " << std::setw(6) << float(hydrogen_bond_info_list.size())/float(total_sample_res) << std::endl;
	std::cout << std::endl;


	std::cout << "--------------------------Successfully RAN extract_hydrogen_bonds_statistic in " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds--------------------------" << std::endl;
*/
}

//////////////////////////////////////////////////////////////////////////////////////


void
test_function(){

	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::swa::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::scoring;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	SilentFileData silent_file_data;

	std::string const pdb_file=   option[ in::file::s ]()[1];	
	std::string pose_name;

  size_t found=pdb_file.rfind(".pdb");


  if(found!=std::string::npos){
		pose_name = pdb_file.substr(found+1);
	} else { 
		pose_name=pdb_file;
	}

	std::string const silent_file=pose_name + ".out";


	std::cout << "importing " << pdb_file << std::endl;
	pose::Pose pose;

	import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );

	protocols::rna::make_phosphate_nomenclature_matches_mini( pose);


	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RIBOSE", 18 );
	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 18 );

//	pose::add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE", 6 );
	apply_virtual_rna_residue_variant_type( pose, 19 , false /*apply_check*/) ;


	BinaryRNASilentStruct s( pose, pose_name );
	silent_file_data.write_silent_struct(s, silent_file, false); 



//~/src/mini/bin/parin_test.macosgccrelease -algorithm test_function -s start.pdb -database ~/minirosetta_database 
}

/*
void
silent_struct_slice(){

	using namespace ObjexxFCL;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::rna;

	// First read in any information on pdb read in from silent files.
	// Assume one to one correspondence between number of tags and number of silent_file
	if ( option[ in::file::silent ].user() ==false) utility_exit_with_message("option[ in::file::silent ].user() ==false");
	if ( option[ in::file::tags ].user() ==false) utility_exit_with_message("option[ in::file::tags ].user() ==false");
	if ( option[input_res].user() ==false) utility_exit_with_message("option[input_res].user() ==false)");


	utility::vector1< Size >  const delete_res_list =	option[delete_res]()[1];

	std::string silent_files_in = option[ in::file::silent ]()[1];
	std::string input_tags = option[ in::file::tags ]()[1];  


	pose::Pose pose;
	import_pose_from_silent_file(pose, silent_files_in_[ silent_file_num ], input_tags_[silent_file_num] );

	for( Size seq_num = pose.total_residue(); seq_num >= 1 ; seq_num--){
 		if(Contain_seq_num(seq_num, delete_res_list)){
			pose.conformation().delete_residue_slow(seq_num );	
		}
	}
	


}
*/

//~/src/mini/bin/parin_test.macosgccrelease -algorithm get_pose_energy_breakdown -s start.pdb -database ~/minirosetta_database > output.txt
//~/src/mini/bin/parin_test.macosgccrelease -algorithm get_pose_energy_breakdown  -in:file:silent_struct_type binary_rna -in:file:silent region_0_2_sample.cluster.out  -tags S_$(Process)   -database ~/minirosetta_database > output.txt


void
get_pose_energy_breakdown(){


	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::kinematics;
	using namespace core::optimization;
	using namespace core::io::silent;
	using namespace protocols::swa::rna;
	using namespace core::conformation;
	using namespace ObjexxFCL;
	using namespace core::io::silent;
	using namespace core::pose;

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	SilentFileData silent_file_data;

	pose::Pose pose;

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn(); 

	if(option[ in::file::s].user()){

		std::string const pdb_file= option[ in::file::s ]()[1];

		std::cout << "importing pose from pdb_file: " << pdb_file << std::endl;

		import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );
	
		protocols::rna::make_phosphate_nomenclature_matches_mini( pose);

	}else if( (option[ in::file::silent ].user()) && (option[ in::file::tags ].user() ) ){ 
		std::string const silent_file=option[ in::file::silent]()[1];
		std::string const input_tag=option[ in::file::tags ]()[1];

		std::cout << "importing pose : " << input_tag << " from silent_file: " << silent_file <<std::endl;

		import_pose_from_silent_file(pose,  silent_file ,  input_tag);

	}else{
		utility_exit_with_message( "Must pdb file or silent_struct!" );
	}

	(*scorefxn)(pose);

	pose.energies().show(std::cout);

}

void
minimize_pdb(){

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
	using namespace protocols::swa::rna;
	using namespace ObjexxFCL;
	using namespace core::optimization;


	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	core::scoring::ScoreFunctionOP scorefxn=create_scorefxn();

	if ( !option[ in::file::s ].user() ) utility_exit_with_message( "User must supply in::file::s!" );
	if ( !option[ minimize_res ].user() ) utility_exit_with_message( "User must supply minimize_res!" );
	if ( !option[ fold_tree_strings].user() ) utility_exit_with_message( "User must supply fold_tree_strings!" );


	utility::vector1< std::string > const	pdb_tag_list= option[ in::file::s ]();	
	utility::vector1< Size > const	minimize_res_list= option[ minimize_res ]();	
	utility::vector1< std::string > const	fold_tree_string_list= option[ fold_tree_strings ]();	

  AtomTreeMinimizer minimizer;
  float const dummy_tol( 0.00000025);
  bool const use_nblist( true );
  MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
  options.nblist_auto_update( true );

	std::string minimize_res_string="";

	core::kinematics::MoveMap mm;

	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	Output_seq_num_list("minimize_res_list= ", minimize_res_list, TR, 30 );

	for(Size n=1; n<=minimize_res_list.size(); n++){
		Size const seq_num=minimize_res_list[n];
		mm.set_bb(seq_num, true);
		mm.set_chi(seq_num, true );

		minimize_res_string+= string_of(seq_num) + "_";
	}		

	for(Size n=1; n<=pdb_tag_list.size(); n++){

		std::string const pdb_tag=pdb_tag_list[n];
		std::string const output_pdb_tag= "minimize_res_" + minimize_res_string + get_tag_from_pdb_filename(pdb_tag);

		pose::Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set, pdb_tag );

		if(option[ graphic ]()) protocols::viewer::add_conformation_viewer( pose.conformation(), pdb_tag, 400, 400 );


		kinematics::FoldTree const start_fold_tree=pose.fold_tree();
		if(start_fold_tree.is_simple_tree()==false){
			Output_fold_tree_info(pose, 	pdb_tag); std::cout << std::endl;	
			utility_exit_with_message("Error: imported_pose does not have a simple fold_tree");
		}	
		Size const nres=pose.total_residue();

		core::kinematics::FoldTree fold_tree( nres );

		for(Size ii=1; ii<=fold_tree_string_list.size(); ii++){

			utility::vector1< std::string > const fold_tree_string=Tokenize(fold_tree_string_list[ii], "-");
			if(fold_tree_string.size()!=3){
				 utility_exit_with_message( "fold_tree_string.size()!=3, fold_tree_string_list[ii]= " + fold_tree_string_list[ii] );
			}
			Size const five_prime_seq_num=string_to_int(fold_tree_string[1]);
			Size const three_prime_seq_num=string_to_int(fold_tree_string[2]);
			Size const cut_point=string_to_int(fold_tree_string[3]);

			std::cout << "five_prime_seq_num= " << five_prime_seq_num << " three_prime_seq_num= " << three_prime_seq_num << " cut_point= " << cut_point << std::endl;

			fold_tree.new_jump( five_prime_seq_num, three_prime_seq_num, cut_point );
		}		


		move_jump_atom_to_base(fold_tree, pose.sequence());

		pose.fold_tree( fold_tree);
		
		std::cout << "----------------------------------------------------------" << std::endl;
		std::cout << "FINAL_fold_Tree:" << std::endl;
		Output_fold_tree_info(pose, 	pdb_tag); std::cout << std::endl;	
		std::cout << "----------------------------------------------------------" << std::endl;


		minimizer.run( pose, mm, *(scorefxn), options );

		dump_pdb( pose, output_pdb_tag + ".pdb");	

	}

}


std::string //silly function to convert to real to string
hack_create_torsion_value_string(core::Real const & torsion_value){

		using namespace ObjexxFCL;

		std::string torsion_string="";

		core::Real	const principal_torsion=numeric::principal_angle_degrees( torsion_value);

		Size const principal_torsion_SIZE=std::abs(principal_torsion+0.00001); //0.00001 is to prevent random ambiguity if the torsion decimal value is exactly .0000 Oct 12, 2010


		if(principal_torsion>0){
			torsion_string="p" + lead_zero_string_of(principal_torsion_SIZE, 3);
		}else{
			torsion_string="n" + lead_zero_string_of(principal_torsion_SIZE, 3);
		}

		return torsion_string;
}



std::string //silly function used for appending the rotamer value to the tag
hack_create_rotamer_string( core::pose::Pose const & pose, bool const Is_prepend, Size const moving_res ){

	std::string rotamer_tag="";

	conformation::Residue const & five_prime_rsd= (Is_prepend) ? pose.residue(moving_res): pose.residue(moving_res-1);
	conformation::Residue const & three_prime_rsd= (Is_prepend) ?  pose.residue(moving_res+1) : pose.residue(moving_res);


	rotamer_tag.append("_E" + hack_create_torsion_value_string(five_prime_rsd.mainchain_torsion( 5  ) ) );
	rotamer_tag.append("_Z" + hack_create_torsion_value_string(five_prime_rsd.mainchain_torsion( 6  ) ) );
	rotamer_tag.append("_A" + hack_create_torsion_value_string(three_prime_rsd.mainchain_torsion( 1 ) ) );
	rotamer_tag.append("_B" + hack_create_torsion_value_string(three_prime_rsd.mainchain_torsion( 2 ) ) );
	rotamer_tag.append("_G" + hack_create_torsion_value_string(three_prime_rsd.mainchain_torsion( 3 ) ) );


	if(Is_prepend){
		rotamer_tag.append("_D" + hack_create_torsion_value_string(five_prime_rsd.mainchain_torsion( 4 ) ) );
		rotamer_tag.append("_C" + hack_create_torsion_value_string(five_prime_rsd.chi(  1) ) );

	}else{
		rotamer_tag.append("_D" + hack_create_torsion_value_string(three_prime_rsd.mainchain_torsion( 4) ) );
		rotamer_tag.append("_C" + hack_create_torsion_value_string(three_prime_rsd.chi( 1) ) );
	}

	return rotamer_tag;

}

bool
Is_new_cluster_center_rotamer(pose::Pose const & pose, 
														pose::Pose & cluster_center_pose, 
														Size const rotamer_count,
														Size const sample_res,
														bool const Is_prepend, 
														Real const cluster_rmsd, 
														utility::vector1< utility::vector1< Torsion_Info > > const & cluster_center_rotamer_list){

	bool const verbose=false;

	for(Size cluster_ID=cluster_center_rotamer_list.size(); cluster_ID>=1; cluster_ID--){

		apply_rotamer( cluster_center_pose, cluster_center_rotamer_list[cluster_ID] );

		Real const rmsd=suite_rmsd(pose, cluster_center_pose, sample_res, Is_prepend, true); 			

		if(rmsd < cluster_rmsd) {
			if(verbose) std::cout << "rotamer_count= "<< rotamer_count << " is part of cluster center " << cluster_ID << ". RMSD= " << rmsd << std::endl;
			return false;
		}  

	}
	
	return true;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
get_residue_xyz_list(pose::Pose const & pose, Size const sample_res, bool const Is_prepend, utility::vector1< numeric::xyzVector<Real> > & xyz_list){

		xyz_list.clear();

		Size const num_heavy_backbone_atoms=11; //RNA contain 11 heavy backbone atoms.

		conformation::Residue const & base_rsd= pose.residue(sample_res);

		conformation::Residue const & phosphate_rsd= (Is_prepend) ? pose.residue(sample_res+1) : pose.residue(sample_res);

		if(num_heavy_backbone_atoms!=(base_rsd.first_sidechain_atom()-1) ){
			std::cout << "num_heavy_backbone_atoms= " << num_heavy_backbone_atoms << std::endl;
			std::cout << "base_rsd.first_sidechain_atom()= " << base_rsd.first_sidechain_atom() << std::endl;
		 	utility_exit_with_message("num_heavy_backbone_atoms!=(base_rsd.first_sidechain_atom()-1)" );
		}

		//Basically want the Base atoms at the top of the list since these atoms will are furthest from anchor and hence largest variation.
		for ( Size atomno=base_rsd.first_sidechain_atom()+1; atomno<= base_rsd.nheavyatoms(); atomno++ ) { //rsd.first_sidechain_atom()+1 to not include the O2star oxygen.

			if(base_rsd.atom_type(atomno).name()=="VIRT"){
				std::cout << "base_rsd.atom_type(atomno).name()==\"VIRT\"!, atomno= " << atomno << std::endl;
				 utility_exit_with_message("base_rsd.atom_type(atomno).name()==\"VIRT\"!");
			}


			xyz_list.push_back(base_rsd.xyz(atomno));

		}


		for ( Size atomno=1; atomno<= base_rsd.first_sidechain_atom(); atomno++ ){ //rsd.first_sidechain_atom() is the O2star oxygen.

			if(Is_prepend && (atomno<=4) ){

				if(phosphate_rsd.atom_type(atomno).name()=="VIRT"){
					std::cout << "phosphate_rsd.atom_type(atomno).name()==\"VIRT\"!, atomno= " << atomno << std::endl;
					std::cout << "atomno= " << atomno << std::endl;
					utility_exit_with_message("phosphate_rsd.atom_type(atomno).name()==\"VIRT\"");
				}

				xyz_list.push_back(phosphate_rsd.xyz(atomno));


			}else{

				if(base_rsd.atom_type(atomno).name()=="VIRT"){
					std::cout << "base_rsd.atom_type(atomno).name()==\"VIRT\"!, atomno= " << atomno << std::endl;
					utility_exit_with_message("base_rsd.atom_type(atomno).name()==\"VIRT\"!");
				}

				xyz_list.push_back(base_rsd.xyz(atomno));

			}

		}
 
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool
Is_new_cluster_center_xyz(pose::Pose const & pose, 
												 Size const rotamer_count, 
												 Size const sample_res, 
												 bool const Is_prepend, 
												 Real const cluster_rmsd,  
												 utility::vector1< utility::vector1< numeric::xyzVector<Real> > >	& cluster_center_xyz_list){

	bool const verbose=false;

	utility::vector1< numeric::xyzVector<Real> > curr_xyz_list;	

	get_residue_xyz_list( pose, sample_res, Is_prepend, curr_xyz_list);

	Real const CUTOFF_sum_square_deviation=(cluster_rmsd*cluster_rmsd)*(curr_xyz_list.size());


	for(Size cluster_ID=cluster_center_xyz_list.size(); cluster_ID>=1; cluster_ID--){ //more likely to be cluster a cluster member of a neighbor rotamer.

		//check
		if(cluster_center_xyz_list[cluster_ID].size()!=curr_xyz_list.size()){
			std::cout << "cluster_center_xyz_list[cluster_ID].size()!=curr_xyz_list.size()" << std::endl;
			std::cout << "cluster_ID=" << cluster_ID << std::endl;
			std::cout << "cluster_center_xyz_list[cluster_ID].size()" << cluster_center_xyz_list[cluster_ID].size() << std::endl;
			std::cout << "curr_xyz_list.size()" << curr_xyz_list.size() << std::endl;
			utility_exit_with_message("cluster_center_xyz_list[cluster_ID].size()!=curr_xyz_list.size()");
		}

		Real sum_square_deviation=0;

		for(Size atomno=1; atomno<=	curr_xyz_list.size(); atomno++){

			sum_square_deviation+=( curr_xyz_list[atomno]-cluster_center_xyz_list[cluster_ID][atomno] ).length_squared();

			if(sum_square_deviation>=CUTOFF_sum_square_deviation) break;


			if(atomno==curr_xyz_list.size()){
				if(verbose){
				 	std::cout << "rotamer_count= "<< rotamer_count << " is part of cluster center " << cluster_ID;
					std::cout << ". num_atoms= " << curr_xyz_list.size() << " sum_square_deviation= " << sum_square_deviation;
					std::cout << " CUTOFF_sum_square_deviation= " << CUTOFF_sum_square_deviation << std::endl;
				}
				return false;
			}

		}
	
	}		

	cluster_center_xyz_list.push_back(curr_xyz_list);
	return true;
}


bool
Is_cluster_member(utility::vector1< numeric::xyzVector<Real> > const & cluster_center_xyzs, 
								 utility::vector1< numeric::xyzVector<Real> > const & curr_xyz_list, 
								 Real const CUTOFF_sum_square_deviation){


	if(cluster_center_xyzs.size()!=curr_xyz_list.size()){
		std::cout << "cluster_center_xyzs.size()!=curr_xyz_list.size()" << std::endl;
		std::cout << "cluster_center_xyzs.size()= " << cluster_center_xyzs.size() << std::endl;
		std::cout << "curr_xyz_list.size()=" << curr_xyz_list.size() << std::endl;
		utility_exit_with_message("cluster_center_xyzs.size()!=curr_xyz_list.size()");
	}

	Real sum_square_deviation=0;

	for(Size atomno=1; atomno<=	curr_xyz_list.size(); atomno++){

		sum_square_deviation+=( curr_xyz_list[atomno]-cluster_center_xyzs[atomno] ).length_squared();

		if(sum_square_deviation>=CUTOFF_sum_square_deviation) return false;

	}

	return true;

}

bool
Is_new_cluster_center_second_stage(pose::Pose const & pose, 
												 				Size const rotamer_count, 
												 				Size const sample_res, 
												 				bool const Is_prepend, 
												 				Real const cluster_rmsd, 
												 				Real const large_cluster_rmsd, 
																utility::vector1< utility::vector1<Size> > & large_to_small_cluster_center_map_list, 
												 				utility::vector1< utility::vector1< numeric::xyzVector<Real> > >	& small_cluster_center_xyz_list,
												 				utility::vector1< utility::vector1< numeric::xyzVector<Real> > >	const & large_cluster_center_xyz_list){


	utility::vector1< numeric::xyzVector<Real> > curr_xyz_list;	

	get_residue_xyz_list( pose, sample_res, Is_prepend, curr_xyz_list);

	Real const first_stage_CUTOFF_sum_square_deviation=((cluster_rmsd+large_cluster_rmsd)*(cluster_rmsd+large_cluster_rmsd))*(curr_xyz_list.size());
	Real const second_stage_CUTOFF_sum_square_deviation=((cluster_rmsd)*(cluster_rmsd))*(curr_xyz_list.size());

	for(Size large_cluster_ID=large_cluster_center_xyz_list.size(); large_cluster_ID>=1; large_cluster_ID--){ //more likely to be cluster a cluster member of a neighbor rotamer.

		//Use triangular identity here.
		if(Is_cluster_member(large_cluster_center_xyz_list[large_cluster_ID], curr_xyz_list, first_stage_CUTOFF_sum_square_deviation)==false) continue;

		utility::vector1<Size> const & large_to_small_cluster_center_map=large_to_small_cluster_center_map_list[large_cluster_ID];

		for(Size n=1; n<=large_to_small_cluster_center_map.size(); n++){
			Size const small_cluster_ID=large_to_small_cluster_center_map[n];
			if(Is_cluster_member(small_cluster_center_xyz_list[small_cluster_ID], curr_xyz_list, second_stage_CUTOFF_sum_square_deviation)){
				 return false;
			}			
		}

	}

	//OK if reach this point, then a new cluster.

	small_cluster_center_xyz_list.push_back(curr_xyz_list);
	bool found_large_to_small_map=false;

	Real const large_CUTOFF_sum_square_deviation=((large_cluster_rmsd)*(large_cluster_rmsd))*(curr_xyz_list.size());

	for(Size large_cluster_ID=large_cluster_center_xyz_list.size(); large_cluster_ID>=1; large_cluster_ID--){ //more likely to be cluster a cluster member of a neighbor rotamer.

		Real sum_square_deviation=0;

		for(Size atomno=1; atomno<=curr_xyz_list.size(); atomno++){

			sum_square_deviation+=( curr_xyz_list[atomno]-large_cluster_center_xyz_list[large_cluster_ID][atomno] ).length_squared();

		} 

		if(sum_square_deviation<large_CUTOFF_sum_square_deviation){	
			large_to_small_cluster_center_map_list[large_cluster_ID].push_back(small_cluster_center_xyz_list.size() );
			found_large_to_small_map=true;
			break;
		}	
	}

	if(found_large_to_small_map==false){
		std::cout << "found_large_to_small_map==false!" << std::endl;
		std::cout << "rotamer_count= " << rotamer_count << std::endl;
		utility_exit_with_message("found_large_to_small_map==false!");
	}
	

	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////





//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Ignore atoms that are within 3-bonds
bool
Is_bonded_neighbor_atoms_at_phosphate_interface(std::string const & atom_name_1, Size const seq_num_1, std::string const & atom_name_2, Size const seq_num_2){



	if(Is_O3star_atom(atom_name_1)){

		if(seq_num_1==seq_num_2){ 

			//Take care of this outside...basically O3star is considered as part of the sugar as well.

		}else if(seq_num_1==(seq_num_2-1)){

			if(Is_C5star_atom(atom_name_2)) return true; //3 bond

		}
	}


	if(Is_P_atom(atom_name_1)){

		if(seq_num_1==(seq_num_2+1)){ 

			if(Is_C3star_atom(atom_name_2)) return true; //2 bond

			if(Is_H3star_atom(atom_name_2)) return true; //3 bond

			if(Is_C4star_atom(atom_name_2)) return true; //3 bond

			if(Is_C2star_atom(atom_name_2)) return true; //3 bond

		}else if(seq_num_1==seq_num_2){ 

			if(Is_C5star_atom(atom_name_2)) return true; //2 bond

			if(Is_1H5star_atom(atom_name_2)) return true; //3 bond

			if(Is_2H5star_atom(atom_name_2)) return true; //3 bond

			if(Is_C4star_atom(atom_name_2)) return true; //3 bond

		}
	}


	if( (Is_OP2_atom(atom_name_1) || Is_OP1_atom(atom_name_1) ) ){

		if(seq_num_1==(seq_num_2+1)){ 

			if(Is_C3star_atom(atom_name_2)) return true; //3 bond

		}else if(seq_num_1==seq_num_2){ 

			if(Is_C5star_atom(atom_name_2)) return true; //3 bond

		}
	}

	if(Is_O5star_atom(atom_name_1)){

		if(seq_num_1==(seq_num_2+1)){ 

			if(Is_C3star_atom(atom_name_2)) return true; //3 bond

		}else if(seq_num_1==seq_num_2){ 

			if(Is_C5star_atom(atom_name_2)) return true; //1 bond

			if(Is_1H5star_atom(atom_name_2)) return true; //2 bond

			if(Is_2H5star_atom(atom_name_2)) return true; //2 bond

			if(Is_C4star_atom(atom_name_2)) return true; //2 bond

			if(Is_H4star_atom(atom_name_2)) return true; //3 bond

			if(Is_C3star_atom(atom_name_2)) return true; //3 bond

			if(Is_O4star_atom(atom_name_2)) return true; //3 bond

		}
	}

	return false;
}



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
setup_VDW_rep_atom_map_list(pose::Pose const & pose, utility::vector1< std::pair< std::pair<Size,Size>, utility::vector1< std::pair<Size,Size> > > > & VDW_rep_atom_map_list){

	for(Size seq_num_1=1; seq_num_1<=2; seq_num_1++){

	conformation::Residue const & rsd_1= pose.residue(seq_num_1);


	for(Size atomno_1=1; atomno_1<=rsd_1.natoms(); atomno_1++){ //include hydrogen atoms

		std::string const & atom_name_1=rsd_1.type().atom_name(atomno_1);
		if(atom_name_1=="HO2'") continue;
		if(rsd_1.atom_type(atomno_1).name()=="VIRT") continue; 


		bool const Is_3_prime_phosphate_1=Is_three_prime_phosphate_atom(atom_name_1); //This is just the O3' atom
		bool const Is_5_prime_phosphate_1=Is_five_prime_phosphate_atom(atom_name_1);  //This is O5', OP2, OP1 and P 

		bool const Is_base_sugar_1=(Is_3_prime_phosphate_1==false && Is_5_prime_phosphate_1==false);

		////////////////////////////////////////////////////////////////////////////////////
		std::pair<Size,Size> const res_1_atom_ID = std::make_pair(seq_num_1, atomno_1);
		utility::vector1< std::pair<Size,Size> > res_2_atom_ID_list;
		////////////////////////////////////////////////////////////////////////////////////		

		for(Size seq_num_2=seq_num_1; seq_num_2<=2; seq_num_2++){
			conformation::Residue const & rsd_2= pose.residue(seq_num_2);

			for(Size atomno_2=1; atomno_2<=rsd_2.natoms(); atomno_2++){ //include hydrogen atoms

				if(seq_num_1==seq_num_2){
					if(atomno_1>=atomno_2) continue; 
				}

				std::string const & atom_name_2=rsd_2.type().atom_name(atomno_2);
				if(atom_name_2=="HO2'") continue;
				if(rsd_2.atom_type(atomno_2).name()=="VIRT") continue; 

				bool const Is_3_prime_phosphate_2=Is_three_prime_phosphate_atom(atom_name_2);
				bool const Is_5_prime_phosphate_2=Is_five_prime_phosphate_atom(atom_name_2);

				bool const Is_base_sugar_2=(Is_3_prime_phosphate_2==false && Is_5_prime_phosphate_2==false);

				if(seq_num_1==seq_num_2){ //same nucleotide
					if( (Is_base_sugar_1 ) && (Is_base_sugar_2 ) ) continue; //Ignore clash between same base-sugar group

					if( Is_3_prime_phosphate_1 && Is_3_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group
					if( Is_5_prime_phosphate_1 && Is_5_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group

					if( Is_3_prime_phosphate_1 && Is_base_sugar_2) continue; //O3' atom is actually part of the 3' sugar
					if( Is_3_prime_phosphate_2 && Is_base_sugar_1) continue; //O3' atom is actually part of the 3' sugar

				}

				if(seq_num_1==1 && seq_num_2==2){ 
					if( Is_3_prime_phosphate_1 && Is_5_prime_phosphate_2) continue; //Ignore clash between same phosphate group
				} 

				if(seq_num_2==1 && seq_num_1==2){ 
					if( Is_3_prime_phosphate_2 && Is_5_prime_phosphate_1) continue; //Ignore clash between same phosphate group
				}

				if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_1, seq_num_1, atom_name_2, seq_num_2)) continue;
				if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_2, seq_num_2, atom_name_1, seq_num_1)) continue;

				if(atom_name_1==atom_name_2 && seq_num_1==seq_num_2){
					std::cout << "atom_name_1=" << atom_name_1 << std::endl;
					std::cout << "seq_num_1="   << seq_num_1   << std::endl;
					std::cout << "atom_name_2=" << atom_name_2 << std::endl;
					std::cout << "seq_num_2="   << seq_num_2   << std::endl;
					utility_exit_with_message("atom_name_1==atom_name_2 && seq_num_1==seq_num_2");
				}

				res_2_atom_ID_list.push_back(std::make_pair(seq_num_2, atomno_2));

			}//atomno_2
		}//seq_num_2

		VDW_rep_atom_map_list.push_back(std::make_pair(res_1_atom_ID, res_2_atom_ID_list));
	}//atomno_1
	}//seq_num_1

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool
pass_VDW_replusion_screen_fast(pose::Pose const & pose, Real const VDW_overlap_dist_cutoff, Size const num_atom_VDW_clash_cutoff,
														 utility::vector1< std::pair< std::pair<Size,Size>, utility::vector1< std::pair<Size,Size> > > > const & VDW_rep_atom_map_list){

	bool verbose=false;

	Size num_atom_clash_so_far=0;

	for(Size n=1; n<=VDW_rep_atom_map_list.size(); n++){
		
		Size const seq_num_1=VDW_rep_atom_map_list[n].first.first;
		Size const atomno_1=VDW_rep_atom_map_list[n].first.second;

		for(Size ii=1; ii<=VDW_rep_atom_map_list[n].second.size(); ii++){

			Size const seq_num_2=VDW_rep_atom_map_list[n].second[ii].first;
			Size const atomno_2=VDW_rep_atom_map_list[n].second[ii].second;


			Real const VDW_radius_ONE=pose.residue(seq_num_1).atom_type(atomno_1).lj_radius();
			Real const VDW_radius_TWO=pose.residue(seq_num_2).atom_type(atomno_2).lj_radius();

			Real const cutoff_sum_VDW_radius=VDW_radius_ONE+VDW_radius_TWO-VDW_overlap_dist_cutoff;
	
			if( ( pose.residue(seq_num_1).xyz(atomno_1)-pose.residue(seq_num_2).xyz(atomno_2) ).length_squared()<(cutoff_sum_VDW_radius*cutoff_sum_VDW_radius) ){

				num_atom_clash_so_far++;
				if(verbose){
					std::cout << "FAST_VDW CLASH between: " << seq_num_1 << pose.residue(seq_num_1).type().atom_name(atomno_1) << " and " << seq_num_2 << pose.residue(seq_num_2).type().atom_name(atomno_2);
					std::cout << " VDW_radius_ONE= " << VDW_radius_ONE << " VDW_radius_TWO= " << VDW_radius_TWO; 
					std::cout << " VDW_overlap_dist_cutoff= " << VDW_overlap_dist_cutoff << " cutoff_sum_VDW_radius= " << cutoff_sum_VDW_radius;
					std::cout << " atom_atom_dist= " << ( pose.residue(seq_num_1).xyz(atomno_1)-pose.residue(seq_num_2).xyz(atomno_2) ).length();
					std::cout << " num_atom_clash= " << num_atom_clash_so_far << "/" << num_atom_VDW_clash_cutoff << std::endl;
				} 				
				if(num_atom_clash_so_far>=num_atom_VDW_clash_cutoff) return false;
			}

		}

	}

	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool
pass_VDW_replusion_screen_slow(pose::Pose const & pose, Real const VDW_overlap_dist_cutoff, Size const num_atom_VDW_clash_cutoff){

	bool verbose=true;
	Size num_atom_clash_so_far=0;

	for(Size seq_num_1=1; seq_num_1<=2; seq_num_1++){
	for(Size seq_num_2=seq_num_1; seq_num_2<=2; seq_num_2++){

	conformation::Residue const & rsd_1= pose.residue(seq_num_1);
	conformation::Residue const & rsd_2= pose.residue(seq_num_2);


	for(Size atomno_1=1; atomno_1<=rsd_1.natoms(); atomno_1++){ //include hydrogen atoms

		std::string const & atom_name_1=rsd_1.type().atom_name(atomno_1);
		if(atom_name_1=="HO2'") continue;
		if(rsd_1.atom_type(atomno_1).name()=="VIRT") continue; 


		bool const Is_3_prime_phosphate_1=Is_three_prime_phosphate_atom(atom_name_1); //This is just the O3' atom
		bool const Is_5_prime_phosphate_1=Is_five_prime_phosphate_atom(atom_name_1);  //This is O5', OP2, OP1 and P 

		bool const Is_base_sugar_1=(Is_3_prime_phosphate_1==false && Is_5_prime_phosphate_1==false);

		for(Size atomno_2=1; atomno_2<=rsd_2.natoms(); atomno_2++){ //include hydrogen atoms

			if(seq_num_1==seq_num_2){
				if(atomno_1>=atomno_2) continue; 
			}

			std::string const & atom_name_2=rsd_2.type().atom_name(atomno_2);
			if(atom_name_2=="HO2'") continue;
			if(rsd_2.atom_type(atomno_2).name()=="VIRT") continue; 

			bool const Is_3_prime_phosphate_2=Is_three_prime_phosphate_atom(atom_name_2);
			bool const Is_5_prime_phosphate_2=Is_five_prime_phosphate_atom(atom_name_2);

			bool const Is_base_sugar_2=(Is_3_prime_phosphate_2==false && Is_5_prime_phosphate_2==false);

			if(seq_num_1==seq_num_2){ //same nucleotide
				if( (Is_base_sugar_1 ) && (Is_base_sugar_2 ) ) continue; //Ignore clash between same base-sugar group

				if( Is_3_prime_phosphate_1 && Is_3_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group
				if( Is_5_prime_phosphate_1 && Is_5_prime_phosphate_2 ) continue; //Ignore clash between same phosphate group

				if( Is_3_prime_phosphate_1 && Is_base_sugar_2) continue; //O3' atom is actually part of the 3' sugar
				if( Is_3_prime_phosphate_2 && Is_base_sugar_1) continue; //O3' atom is actually part of the 3' sugar

			}

			if(seq_num_1==1 && seq_num_2==2){ 
				if( Is_3_prime_phosphate_1 && Is_5_prime_phosphate_2) continue; //Ignore clash between same phosphate group
			} 

			if(seq_num_2==1 && seq_num_1==2){ 
				if( Is_3_prime_phosphate_2 && Is_5_prime_phosphate_1) continue; //Ignore clash between same phosphate group
			}

			if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_1, seq_num_1, atom_name_2, seq_num_2)) continue;
			if(Is_bonded_neighbor_atoms_at_phosphate_interface(atom_name_2, seq_num_2, atom_name_1, seq_num_1)) continue;

			if(atom_name_1==atom_name_2 && seq_num_1==seq_num_2){
				std::cout << "atom_name_1=" << atom_name_1 << std::endl;
				std::cout << "seq_num_1="   << seq_num_1   << std::endl;
				std::cout << "atom_name_2=" << atom_name_2 << std::endl;
				std::cout << "seq_num_2="   << seq_num_2   << std::endl;
				utility_exit_with_message("atom_name_1==atom_name_2 && seq_num_1==seq_num_2");
			}

			Real const VDW_radius_ONE=rsd_1.atom_type(atomno_1).lj_radius();
			Real const VDW_radius_TWO=rsd_2.atom_type(atomno_2).lj_radius();

			Real const cutoff_sum_VDW_radius=VDW_radius_ONE+VDW_radius_TWO-VDW_overlap_dist_cutoff;
	
			if( ( rsd_1.xyz(atomno_1)-rsd_2.xyz(atomno_2) ).length_squared()<(cutoff_sum_VDW_radius*cutoff_sum_VDW_radius) ){

				num_atom_clash_so_far++;
				if(verbose){
					std::cout << "SLOW_VDW CLASH between: " << seq_num_1 << rsd_1.type().atom_name(atomno_1) << " and " << seq_num_2 << rsd_2.type().atom_name(atomno_2);
					std::cout << " VDW_radius_ONE= " << VDW_radius_ONE << " VDW_radius_TWO= " << VDW_radius_TWO; 
					std::cout << " VDW_overlap_dist_cutoff= " << VDW_overlap_dist_cutoff << " cutoff_sum_VDW_radius= " << cutoff_sum_VDW_radius;
					std::cout << " atom_atom_dist= " << ( rsd_1.xyz(atomno_1)-rsd_2.xyz(atomno_2) ).length();
					std::cout << " num_atom_clash= " << num_atom_clash_so_far << "/" << num_atom_VDW_clash_cutoff << std::endl;
				} 				
				if(num_atom_clash_so_far>=num_atom_VDW_clash_cutoff) return false;
			}
		}
	}
	}
	}

	return true;

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
cluster_rotamers(bool const second_stage,
								bool const optimize_screening,
								Real const cluster_rmsd,
								Real const large_cluster_rmsd, 
								utility::vector1< utility::vector1< numeric::xyzVector<Real> > > const & input_cluster_center_xyz_list, 
								utility::vector1< utility::vector1< numeric::xyzVector<Real> > > & cluster_center_xyz_list){

	if(second_stage){
		Output_title_text("Enter Rosetta cluster_rotamers function SECOND_STAGE", TR );
	}else{
		Output_title_text("Enter Rosetta cluster_rotamers function FIRST_STAGE", TR );
	}

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::io::silent;
	using namespace core::id;
	using namespace core::kinematics;
	using namespace ObjexxFCL;
	using namespace ObjexxFCL::fmt;

	clock_t const time_start( clock() );

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	//Create a dinucleotide

	if ( !option[ dinucleotide_sequence ].user() ) utility_exit_with_message( "User must supply dinucleotide_sequence!" );

	std::string const sequence=option[ dinucleotide_sequence ](); //"aa";

	if(sequence.size()!=2) utility_exit_with_message( "sequence.size()!=2, sequence=" + sequence);

	bool const Is_prepend=false;
	Size const sample_res=(Is_prepend) ? 1 : 2;
	Size const sample_suite= (Is_prepend) ? sample_res : sample_res-1;
	Size const anchor_res=(Is_prepend) ? 2 : 1;
	bool const replusion_screen=option[ cluster_rotamer_replusion_screen ]();
	bool const VDW_rep_screening_slow_check=option[ cluster_rotamer_VDW_rep_screening_slow_check ] ();

	std::string const silent_file = option[output_silent_file](); 
	bool create_rotamer_silent_file= (silent_file!="") ? true : false;

	bool const two_stage_clustering =option[ two_stage_rotamer_clustering ]();
	bool const QUICK =option[ quick_test ]();
	bool const sparse_output= option[ cluster_rotamer_sparse_output]();
	Size const bin_size= option[ cluster_rotamer_bin_size]();
	Real const VDW_overlap_dist_cutoff=0.5;
	Size const num_atom_VDW_clash_cutoff=4;  //(0.8,3), (0.8,1), (0.5,1)
	
	if(two_stage_clustering==true && create_rotamer_silent_file==true){
		utility_exit_with_message("two_stage_clustering==true && create_rotamer_silent_file==true");
	}

	std::cout << "------------------------------------------------------" << std::endl;
	std::cout << "sequence=" << sequence << std::endl;
	Output_boolean("Is_prepend= ", Is_prepend, TR ); std::cout << std::endl;
	std::cout << "sample_res= " << sample_res << std::endl;
	std::cout << "sample_suite= " << sample_suite << std::endl;
	std::cout << "bin_size= " << bin_size << std::endl;
	std::cout << "cluster_rmsd= " << cluster_rmsd << std::endl;
	std::cout << "large_cluster_rmsd= " << large_cluster_rmsd << std::endl;
	std::cout << "VDW_overlap_dist_cutoff= " << VDW_overlap_dist_cutoff << std::endl;
	std::cout << "num_atom_VDW_clash_cutoff= " << num_atom_VDW_clash_cutoff << std::endl;
	Output_boolean("replusion_screen= ", replusion_screen, TR ); std::cout << std::endl;
	std::cout << "silent_file= " << silent_file << std::endl;
	Output_boolean("QUICK= ", QUICK, TR ); std::cout << std::endl;
	Output_boolean("sparse_output= ", sparse_output, TR ); std::cout << std::endl;
	Output_boolean("create_rotamer_silent_file= ", create_rotamer_silent_file, TR ); std::cout << std::endl;
	Output_boolean("VDW_rep_screening_slow_check= ", VDW_rep_screening_slow_check, TR ); std::cout << std::endl;
	Output_boolean("optimize_screening= ", optimize_screening ); std::cout << std::endl;




	std::cout << "------------------------------------------------------" << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////
	pose::Pose pose;
	make_pose_from_sequence( pose, sequence, *rsd_set, false /*auto_termini*/);

	if(option[ graphic ]()) protocols::viewer::add_conformation_viewer( pose.conformation(), "pose", 400, 400 );


	dump_pdb( pose, "before_set_to_A_form_" + sequence + ".pdb");	

	set_nucleotide_to_A_form( pose, 1);
	set_nucleotide_to_A_form( pose, 2);

	dump_pdb( pose, "after_set_to_A_form_" + sequence + ".pdb");	
	
	dump_pdb( pose, "before_virtualize_" + sequence + ".pdb");	

	
	//add_variant_type_to_pose_residue( pose, "VIRTUAL_RNA_RESIDUE_EXCLUDE_PHOSPHATE", anchor_res );	
	add_variant_type_to_pose_residue( pose, "VIRTUAL_PHOSPHATE", 1 );	
	
	dump_pdb( pose, "after_virtualize_" + sequence + ".pdb");	
	
	////////////////////////////////////////////////////////////////////////////////////////////////////
	core::kinematics::FoldTree fold_tree=pose.fold_tree();
	if(fold_tree.is_simple_tree()==false) utility_exit_with_message( "starting fold_tree is not a simple tree!" );
	
	core::kinematics::FoldTree rerooted_fold_tree = fold_tree;
	bool reorder_went_OK = rerooted_fold_tree.reorder( anchor_res );
	if( !reorder_went_OK) utility_exit_with_message( "!reorder_went_OK" );
	pose.fold_tree( rerooted_fold_tree); 		
	////////////////////////////////////////////////////////////////////////////////////////////////////


	utility::vector1< core::Size > working_moving_suite_list;
	working_moving_suite_list.push_back(sample_suite);

	bool const sample_sugar_and_base1= (Is_prepend) ? true : false;
	bool const sample_sugar_and_base2= (Is_prepend) ? false : true;


	StepWiseRNA_RotamerGeneratorWrapperOP rotamer_generator = new StepWiseRNA_RotamerGeneratorWrapper( pose,
																																																	working_moving_suite_list,
																																																	sample_sugar_and_base1,
																																																	sample_sugar_and_base2);

	/*
		README_SETUP.write( "command+= '-sampler_extra_epsilon_rotamer true '\n\n" )
		README_SETUP.write( "command+= '-sampler_extra_beta_rotamer false '\n\n" ) #updated the non_extra rotamer mode to include 80 and -80
		README_SETUP.write( "command+= '-sampler_extra_chi_rotamer false '\n\n" ) #change to false on Oct 21, 2010...speed up code.
	*/

	rotamer_generator->set_fast( false );
	rotamer_generator->set_include_syn_chi(true);
	rotamer_generator->set_extra_epsilon(true);
	rotamer_generator->set_extra_beta(false);
	rotamer_generator->set_extra_anti_chi(false);
	rotamer_generator->set_extra_syn_chi(false);
	rotamer_generator->set_exclude_alpha_beta_gamma_sampling(false);
	rotamer_generator->set_allow_syn_pyrimidine(false);
	rotamer_generator->set_bin_size(bin_size); 
	rotamer_generator->initialize_rotamer_generator_list();

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	core::scoring::ScoreFunctionOP scorefxn = new ScoreFunction;

	//scorefxn_->set_weight( fa_atr  , 0.23 ); //standard weight
	//scorefxn_->set_weight( fa_rep  , 0.12 ); //standard weight

	scorefxn->set_weight( fa_rep  , 1.0 ); 
	scorefxn->set_weight( fa_atr  , 1.0 ); 
	scorefxn->set_weight( fa_intra_rep  , 1.0 ); 

	(*scorefxn)(pose);

	EnergyMap const & start_energy_map = pose.energies().total_energies();

	Real const start_rep_score = scorefxn->get_weight(fa_rep) * start_energy_map[scoring::fa_rep];
	Real const start_atr_score =scorefxn->get_weight(fa_atr) * start_energy_map[scoring::fa_atr];
	Real const start_intra_rep_score = scorefxn->get_weight(fa_intra_rep) * start_energy_map[scoring::fa_intra_rep];

	std::cout << "start_rep_score= " << start_rep_score << std::endl;
	std::cout << "start_atr_score= " << start_atr_score << std::endl;
	std::cout << "start_intra_rep_score= " << start_intra_rep_score << std::endl;

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	SilentFileData silent_file_data;
	
	utility::vector1< utility::vector1< Torsion_Info > > cluster_center_rotamer_list;
	cluster_center_rotamer_list.clear();
	
	utility::vector1< utility::vector1< numeric::xyzVector<Real> > > check_cluster_center_xyz_list;
	check_cluster_center_xyz_list.clear();

	utility::vector1< utility::vector1< Size > > large_to_small_cluster_center_map_list;

	if(second_stage){
		utility::vector1< Size > empty_size_list;
		empty_size_list.clear();
		large_to_small_cluster_center_map_list.clear();
		large_to_small_cluster_center_map_list.assign(input_cluster_center_xyz_list.size(), empty_size_list);	

		if(large_to_small_cluster_center_map_list.size()!=input_cluster_center_xyz_list.size()){
			std::cout << "large_to_small_cluster_center_map_list.size()!=input_cluster_center_xyz_list.size()" << std::endl;
			std::cout << "large_to_small_cluster_center_map_list.size()= " << large_to_small_cluster_center_map_list.size() << std::endl;
			std::cout << "input_luster_center_rotamer_list.size()= " << input_cluster_center_xyz_list.size() << std::endl;
			utility_exit_with_message("large_to_small_cluster_center_map_list.size()!=input_cluster_center_xyz_list.size()");
		}
	}


	Size rotamer_count=0;
	Size cluster_count=0;
	Size pass_VDW_rep_screen_count=0;

	pose::Pose cluster_center_pose=pose;	


/////////////////////////////////
/*
	bool const extra_chi_rotamer=false;
	BaseState base_state = (core::scoring::rna::is_purine( pose.residue( sample_res ) )) ? BOTH: ANTI;	

	PuckerState pucker_state=ALL;
	core::scoring::rna::RNA_FittedTorsionInfo rna_fitted_torsion_info;

	StepWiseRNA_BaseSugarRotamerOP base_sugar_rotamer = new StepWiseRNA_BaseSugarRotamer( base_state, pucker_state, rna_fitted_torsion_info);
	base_sugar_rotamer->set_extra_chi(extra_chi_rotamer);

	Size count=0;

	while(base_sugar_rotamer->get_next_rotamer()){
		count++;

		std::cout << " 	delta1= " <<  F(8, 3, base_sugar_rotamer->delta()) << " 	chi_1= " <<  F(8, 3, base_sugar_rotamer->chi());
		std::cout << " 	nu2_1= " <<  F(8, 3, base_sugar_rotamer->nu2())    << " 	nu1_1= " <<  F(8, 3, base_sugar_rotamer->nu1());
		std::cout << std::endl;

		pose::Pose pose_at_origin = pose; //This make sure that it is not possible to link different rsd in the generated list to the same pose Apr 10, 2010 Parin

		pose.set_torsion( TorsionID( sample_res , id::BB, 4 ) , base_sugar_rotamer->delta());
		pose.set_torsion( TorsionID( sample_res , id::CHI, 1 ) , base_sugar_rotamer->chi());
		pose.set_torsion( TorsionID( sample_res , id::CHI, 2 ) , base_sugar_rotamer->nu2());
		pose.set_torsion( TorsionID( sample_res , id::CHI, 3 ) , base_sugar_rotamer->nu1());

		utility::vector1< Torsion_Info > current_rotamer; //empty

//	}
*/
////////////////////////////////
	utility::vector1< std::pair< std::pair<Size,Size> ,utility::vector1< std::pair<Size,Size> > > > VDW_rep_atom_map_list;
	VDW_rep_atom_map_list.clear();

	if(replusion_screen) setup_VDW_rep_atom_map_list(pose, VDW_rep_atom_map_list);
	

	while( rotamer_generator->has_another_rotamer() ){
		utility::vector1< Torsion_Info > const current_rotamer = rotamer_generator->get_next_rotamer();
		apply_rotamer( pose, current_rotamer);

		rotamer_count++;

		if(QUICK){
			if(rotamer_count>10000) break;
		}

		if(replusion_screen){
			bool const pass_VDW_rep_screen=pass_VDW_replusion_screen_fast(pose, VDW_overlap_dist_cutoff, num_atom_VDW_clash_cutoff, VDW_rep_atom_map_list);


			if(VDW_rep_screening_slow_check){
				bool const pass_VDW_rep_screen_check=pass_VDW_replusion_screen_slow(pose, VDW_overlap_dist_cutoff, num_atom_VDW_clash_cutoff);

				if(pass_VDW_rep_screen!=pass_VDW_rep_screen_check){
					std::cout << "pass_VDW_rep_screen!=pass_VDW_rep_screen_check" << rotamer_count << std::endl;
					Output_boolean("pass_VDW_rep_screen=", pass_VDW_rep_screen, TR ); std::cout << std::endl;
					Output_boolean("pass_VDW_rep_screen=", pass_VDW_rep_screen_check, TR ); std::cout << std::endl;
					utility_exit_with_message("pass_VDW_rep_screen!=pass_VDW_rep_screen_check");
				}	

			}

			if(pass_VDW_rep_screen==false) continue;
			pass_VDW_rep_screen_count++;
		}

		if(create_rotamer_silent_file){
			(*scorefxn)(pose);


			std::string const tag= "S" + hack_create_rotamer_string( pose, Is_prepend, sample_res ) + "_" + lead_zero_string_of(rotamer_count, 8); 
		
			BinaryRNASilentStruct s( pose, tag ); //Does this take a long time to create?
			silent_file_data.write_silent_struct(s, silent_file, true);

			dump_pdb( pose, tag+".pdb");	

		}else{

			bool pass_screen;

			if(second_stage){
				pass_screen=Is_new_cluster_center_second_stage(pose, rotamer_count, sample_res, Is_prepend, cluster_rmsd, large_cluster_rmsd, 
																				 large_to_small_cluster_center_map_list, cluster_center_xyz_list, input_cluster_center_xyz_list);													
			}else{
				pass_screen=Is_new_cluster_center_xyz(pose, rotamer_count, sample_res, Is_prepend, cluster_rmsd, cluster_center_xyz_list);
			}


			if(optimize_screening==false){

				bool pass_screen_check;		

				if(second_stage){
					pass_screen_check=Is_new_cluster_center_xyz(pose, rotamer_count, sample_res, Is_prepend, cluster_rmsd, check_cluster_center_xyz_list);
				}else{
					pass_screen_check=Is_new_cluster_center_rotamer(pose, cluster_center_pose, rotamer_count, sample_res, Is_prepend, cluster_rmsd, cluster_center_rotamer_list);
				}

				if(pass_screen_check) cluster_center_rotamer_list.push_back(current_rotamer);

				if(pass_screen!=pass_screen_check){
					std::cout << "pass_screen!=pass_screen_check for rotamer_count= " << rotamer_count << std::endl;
					Output_boolean("pass_screen=", pass_screen, TR ); std::cout << std::endl;
					Output_boolean("pass_screen_check=", pass_screen_check, TR ); std::cout << std::endl;
					utility_exit_with_message("pass_screen!=pass_screen_check for rotamer_count");
				}				
			}


			if(pass_screen){
				cluster_count++;

				bool output_line=true;
				if(second_stage && sparse_output && ((cluster_count % 100)!=0) ) output_line=false;

				if(output_line){
					if(second_stage){
						std::cout << "SECOND_STAGE_new_cluster:";
					}else{
						std::cout << "FIRST_STAGE_new_cluster:";
					}
					std::cout << " rotamer_count= " << rotamer_count << " cluster_count= " << cluster_count << " time_taken= " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC;

					if(replusion_screen) std::cout << " pass_VDW_rep_screen_count= " << pass_VDW_rep_screen_count;
					
					std::cout << std::endl;

				}
			}
		}
	}

	std::cout << "------------------------------------------------------" << std::endl;
	if(second_stage){
		std::cout << "FINAL OUTPUT (SECOND_TAGE): " << std::endl;
	}else{
		std::cout << "FINAL OUTPUT (FIRST_STAGE): " << std::endl;
	}

	std::cout << "rotamer_cluster_centers=" << cluster_count << std::endl;
	std::cout << "pass_VDW_rep_screen_count= " << pass_VDW_rep_screen_count << std::endl;
	std::cout << "total_rotamer_count=" << rotamer_count << std::endl;
	std::cout << "------------------------------------------------------" << std::endl;

	std::cout << "Total time_taken:" <<  static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC <<std::endl;

	if(second_stage){
		Output_title_text("Exit Rosetta cluster_rotamers function SECOND_STAGE", TR );
	}else{
		Output_title_text("Exit Rosetta cluster_rotamers function FIRST_STAGE", TR );
	}

}

///////////////////////////////////////////////////////////////


void
cluster_rotamers_wrapper(){

	Output_title_text("Enter Rosetta cluster_rotamers WRAPPER ", TR );


	utility::vector1< utility::vector1< numeric::xyzVector<Real> > > first_stage_cluster_center_xyz_list;
	first_stage_cluster_center_xyz_list.clear();

	utility::vector1< utility::vector1< numeric::xyzVector<Real> > > second_stage_cluster_center_xyz_list;
	second_stage_cluster_center_xyz_list.clear();

	utility::vector1< utility::vector1< numeric::xyzVector<Real> > > empty_list;
	empty_list.clear();

	if ( !option[ rotamer_cluster_rmsd ].user() ) utility_exit_with_message( "User must supply rotamer_cluster_rmsd!" );

	bool const two_stage_clustering =option[ two_stage_rotamer_clustering ]();

	Real const large_cluster_rmsd = 2.0;

	Real const cluster_rmsd= option[ rotamer_cluster_rmsd ]();

	bool const optimize_screening= option[cluster_rotamers_optimize_screening]();


	if(two_stage_clustering){

			cluster_rotamers(false, true, large_cluster_rmsd, large_cluster_rmsd,  empty_list, first_stage_cluster_center_xyz_list );

			cluster_rotamers(true, optimize_screening, cluster_rmsd, large_cluster_rmsd, first_stage_cluster_center_xyz_list, second_stage_cluster_center_xyz_list );

	}else{
			cluster_rotamers(false, optimize_screening, cluster_rmsd, large_cluster_rmsd, empty_list, first_stage_cluster_center_xyz_list );
	}

	Output_title_text("Exit Rosetta cluster_rotamers WRAPPER ", TR );


}


///////////////////////////Stole from src/apps/pilot/rhiju/rna_test.cc June 11, 2011//////////////////////////////////
void
rna_idealize_test() {

	pose::Pose pose;
	utility::vector1 <std::string> pdb_files ( option[ in::file::s ]() );

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	bool const close_chainbreaks = option[ idl_close_chainbreaks ];

	for (Size n = 1; n <= pdb_files.size(); n++ ){

		std::string const pdb_file = pdb_files[n];

		pose::Pose pose;
		import_pose::pose_from_pdb( pose, *rsd_set, pdb_file );
		/////////////////////////////////////////
		protocols::rna::make_phosphate_nomenclature_matches_mini( pose );
		/////////////////////////////////////////

		if (!close_chainbreaks) protocols::rna::figure_out_reasonable_rna_fold_tree( pose );

		pose::Pose const start_pose( pose );

		pose = start_pose;

		protocols::idealize::IdealizeMover idealizer;

		// set some options
		if ( option[ coordinate_constraint_weight ].user() ) {
			idealizer.coordinate_constraint_weight( option[ coordinate_constraint_weight ] ) ;
		}
		if ( option[ atom_pair_constraint_weight ].user() ) {
			idealizer.atom_pair_constraint_weight( option[ atom_pair_constraint_weight ] );
		}
		idealizer.fast( false /* option[ fast ] */ );

		idealizer.apply( pose );

		// confirm that nothing changes:
		// actually something *does* change!
		idealizer.apply( pose );

		//test this.
		//		pose::Pose refold_pose;
		//		std::string refold_sequence = pose.sequence();
		//		refold_sequence.erase( refold_sequence.size()-1 );
		//		std::cout << "ABOUT TO MAKE POSE FROM SEQUENCE " << refold_sequence << std::endl;
		//		core::chemical::make_pose_from_sequence( refold_pose, refold_sequence,	*rsd_set );
		//		std::cout << "HEY! " << pose.total_residue() << " " << refold_pose.total_residue() << std::endl;
		//		refold_pose.dump_pdb( "extended.pdb" );
		//		for (Size i = 1; i <= refold_pose.total_residue(); i++ ) copy_rna_torsions( i, i, refold_pose, pose );
		//		refold_pose.dump_pdb( "refold.pdb" );

		pose.dump_pdb( "idealize_"+pdb_file );

	}
}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	std::string algorithm_input = option[algorithm];

//	if(option[ parin_favorite_output ]()){
//		system(std::string("mkdir pose/").c_str()); //Output all the poses generated by the code in here. Parin S Jan 28, 2010
//	}

	if (algorithm_input=="multiple_variant_type_test"){
		multiple_variant_type_test();
	} else if (algorithm_input=="calculate_theoretical_RNA_length"){
		calculate_theoretical_RNA_length();
	} else if (algorithm_input=="journal_club_syn_chi"){
		journal_club_syn_chi();
	}	else if (algorithm_input=="test_function"){
	  test_function();
	} else if (algorithm_input=="minimize_pdb"){
		minimize_pdb();
	}	else if (algorithm_input=="get_pose_energy_breakdown"){
	  get_pose_energy_breakdown();
	}	else if (algorithm_input=="hermann_phase_two"){
	  hermann_phase_two();
	}	else if (algorithm_input=="hermann_phase_two_minimize"){
	  hermann_phase_two_minimize();
	}	else if (algorithm_input=="extract_hydrogen_bonds_statistic"){
	  extract_hydrogen_bonds_statistic();
	}	else if (algorithm_input=="cluster_rotamers"){
	  cluster_rotamers_wrapper();
	} else if (algorithm_input=="rna_idealize" ){ //Stole from src/apps/pilot/rhiju/rna_test.cc
	  rna_idealize_test();
	} else {
		std::cout << "Error no algorithm selected" << std::endl;
	}

	protocols::viewer::clear_conformation_viewers();
  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;
	utility::vector1< Real > blank_real_vector;


	NEW_OPT( delete_res, "delete_res", blank_size_vector); 
	NEW_OPT( sample_res, "sample_res", blank_size_vector ); 
	NEW_OPT( minimize_res, "minimize_res", blank_size_vector ); 
	NEW_OPT( fold_tree_strings, "fold_tree_strings", blank_string_vector );  //five_prime_seq_num-three_prime_seq_num-cut_point 
	NEW_OPT( input_res, "Residues already present in starting file", blank_size_vector );
	NEW_OPT( fast, "quick runthrough for debugging", false );
	NEW_OPT( algorithm, "Specify algorithm to execute", "");
	NEW_OPT( alignment_res_pairs, "alignment_res_pairs", blank_string_vector );  //1-3 4-5 ,res 1 of static to res 3 of moving...res 4 of static to res 5 of moving.
	NEW_OPT( alignment_res_pairs_two, "alignment_res_pairs_two", blank_string_vector );  //1-3 4-5 ,res 1 of static to res 3 of moving...res 4 of static to res 5 of moving.
	NEW_OPT( RMSD_res_pairs, "RMSD_res_pairs", blank_string_vector );  //1-3 4-5 ,RMDS res 1 of first pdb to res 3 of second pdb..res 4 of first pdb to res 5 of second pdb.
	NEW_OPT( alignment_RMSD_CUTOFF, "alignment_RMSD_CUTOFF", 0.0 );  
	NEW_OPT( graphic, "Turn graphic on/off", true); //May 5, 2010 
	NEW_OPT( tag_name, "tag_name", "BLAH_TAG"); //Nov 13, 2010 
	NEW_OPT( output_silent_file, "output_silent_file", "");

	NEW_OPT( lower_to_full_res, "lower_to_full_res", blank_size_vector);
	NEW_OPT( upper_to_full_res, "upper_to_full_res", blank_size_vector);
	NEW_OPT( helical_ends_to_full_res, "helical_ends_to_full_res", blank_size_vector);
	NEW_OPT( cutpoint_closed, "cutpoint_closed", blank_size_vector );
	NEW_OPT( filter_filename, "filter_filename", "");
	NEW_OPT( virtual_res , " virtual_res ", blank_size_vector );
	NEW_OPT( virtual_ribose , " virtual_ribose ", blank_size_vector );
	NEW_OPT( native_virtual_res , " native_virtual_res (use for loop_rmsd calculation)", blank_size_vector );
	NEW_OPT( native_alignment_res , " native_alignment_res (use for loop_rmsd calculation) ", blank_size_vector );
	NEW_OPT( rmsd_res, "optional: residues that will be use to calculate rmsd (use for loop_rmsd calculation)", blank_size_vector );
	NEW_OPT( helical_ends, "helical_ends", "" );
	NEW_OPT( user_num_nstruct_per_node, "num_nstruct_per_node", 0 );
	NEW_OPT( user_JOB_ID, "JOB_ID", 0 );
	NEW_OPT( user_JOB_ID_MOD_CUTOFF, "user_JOB_ID_MOD_CUTOFF", 0 );
	NEW_OPT( USER_BIOX_SUBMIT, "USER_BIOX_SUBMIT", false );
	NEW_OPT( user_skip_minimize, "skip_minimize", false );
	NEW_OPT( user_extra_minimize_rounds, "user_extra_minimize_rounds", false );
	NEW_OPT( align_only_over_base_atoms , "align_only_over_base_atoms", true);
	NEW_OPT( additional_slice_res , "additional_slice_res", blank_size_vector);
	NEW_OPT( INCLUDE_EDGE_PHOSPHATE , "INCLUDE_EDGE_PHOSPHATE", true);
	NEW_OPT( double_count_base_pair , "double_count_base_pair", true);
	NEW_OPT( rotamer_cluster_rmsd , "cluster_rmsd, for cluster_rotamers function", 0.0);
	NEW_OPT( dinucleotide_sequence , "dinucleotide_sequence, for cluster_rotamers function", "");
	NEW_OPT( cluster_rotamers_optimize_screening , "cluster_rotamers_optimize_screening, for cluster_rotamers function", true);
	NEW_OPT( two_stage_rotamer_clustering  , "two_stage_clustering , for cluster_rotamers function", true);
	NEW_OPT( quick_test  , "quick_test , for cluster_rotamers function", false);
	NEW_OPT( cluster_rotamer_sparse_output  , "cluster_rotamer_sparse_output , for cluster_rotamers function", false);
	NEW_OPT( cluster_rotamer_bin_size  , "cluster_rotamer_bin_size , for cluster_rotamers function", 20);
	NEW_OPT( cluster_rotamer_replusion_screen  , "cluster_rotamer_replusion_screen  , for cluster_rotamers function", true);
	NEW_OPT( cluster_rotamer_VDW_rep_screening_slow_check , "cluster_rotamer_VDW_rep_screening_slow_check  , for cluster_rotamers function", true);
	NEW_OPT( idl_close_chainbreaks, "RNA idealize, close chain breaks", false ); //rna_idealize_test()
	NEW_OPT( coordinate_constraint_weight, "coordinate constraint weight", 0.0 ); //rna_idealize_test()
	NEW_OPT( atom_pair_constraint_weight, "atompair constraint weight", 0.0 ); //rna_idealize_test()
	NEW_OPT( list_of_virtual_res, " list of virtual_res of each corresponding imported pdb", blank_string_vector);
	NEW_OPT( list_of_energy, " list of energy of each corresponding imported pdb", blank_real_vector);
	NEW_OPT( native_tag_name, "native tag from a silent_file", "" );
	NEW_OPT( decoy_tag_name, "decoy tag from a silent_file", blank_string_vector);
	NEW_OPT( dump, "dump pdb", false);





  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////

  protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl; 
	} 

}



