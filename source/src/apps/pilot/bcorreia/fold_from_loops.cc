// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief


// libRosetta headers

#include <core/types.hh>


#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSelector.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/conformation/ResidueFactory.hh>

#include <core/fragment/FragmentIO.hh>

#include <protocols/abinitio/AbrelaxApplication.hh>
#include <protocols/abinitio/ClassicAbinitio.hh>
#include <protocols/abinitio/FoldConstraints.hh>
#include <protocols/checkpoint/CheckPointer.hh>

#include <protocols/loops/SlidingWindowLoopClosure.hh>
#include <protocols/loops/SlidingWindowLoopClosure.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Exceptions.hh>

#include <protocols/jumping/util.hh>

#include <protocols/Protocol.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
// Auto-header: duplicate removed #include <protocols/evaluation/RmsdEvaluator.hh>


#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>


//constraints
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>


#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/BinarySilentStruct.hh>


#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// Auto-header: duplicate removed #include <basic/Tracer.hh>
#include <utility/string_util.hh>

#include <apps/pilot/bcorreia/fold_from_loops.hh>

// C++ headers
#include <vector>
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>


//Design headers
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>

//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <protocols/relax/ClassicRelax.hh>
#include <core/util/SwitchResidueTypeSet.hh>


using namespace core;
using namespace protocols;
using namespace kinematics;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static THREAD_LOCAL basic::Tracer TR( "bcorreia_fold_from_loops" );


//brief define_ut_point takes a loop file and a pose figures out how many and where the cutpoints should be based on secondary structure and the
// middle point between loops
std::vector<Size> define_cut_points( protocols::loops::Loops & loops , core::pose::Pose & nat_pose);

// verifies if a given residue is within the range of a loop or not
bool is_loop (protocols::loops::Loops & loops , Size & residue);

// Verifies if it is within the loop or a neighbor range

bool is_loop_neighbor( protocols::loops::Loops & loops, Size & residue, Size & range);


// generates fold trees for topologies with several loops

void fold_tree_generator( protocols::loops::Loops & loops , std::vector<Size> & cutpoints, core::pose::Pose & pose, kinematics::FoldTree & f);

//defines the move map and extends the chain to be folded

void define_movemap_extending_chain( core::kinematics::MoveMapOP & movemap, core::pose::Pose & pose, utility::vector1< protocols::loops::Loop > & loops);

//wrapper to read frag files

void get_fragments( core::fragment::FragSetOP & fragset_large_, core::fragment::FragSetOP & fragset_small_ );


//defines map and extends the full chain

void extending_chain( core::kinematics::MoveMapOP & movemap, core::pose::Pose & pose );

//copying the side chains form de initial loops

void copying_side_chains( core::pose::Pose & nat_pose, core::pose::Pose & fold_pose, protocols::loops::Loops & loops, core::kinematics::MoveMapOP & movemap );


//Count neighbors


//setting up mutation vector

void exclude_loop_residues( core::pose::Pose & pose, utility::vector1< bool > & residues_to_mutate, utility::vector1< bool > & allowed_aas, core::pack::task::PackerTaskOP & task, protocols::loops::Loops & loops);

void refresh_cutpoints( core::pose::Pose & pose, std::vector<Size> & cut_points);


void CA_cst_generator(core::pose::Pose & pose, scoring::constraints::ConstraintSetOP & cst, protocols::loops::Loops & loops, std::vector<Size> & cut_points );


void CA_cst_generator(core::pose::Pose & pose, scoring::constraints::ConstraintSetOP & cst, protocols::loops::Loops & loops );


void new_pose_generator(core::pose::Pose & target_loops, core::pose::Pose & nat_prot ,protocols::loops::Loops & loops);


void copying_side_chains_swap_loop (core::pose::Pose & swap_loops, core::pose::Pose & fold_pose, protocols::loops::Loops & loops, core::kinematics::MoveMapOP & movemap);


///////////////////////////////////////////////////////////////////////////////////////


std::vector<Size> define_cut_points (

		protocols::loops::Loops & loops,
		core::pose::Pose & nat_pose
){

	TR <<"Several Loops Defined "<<std::endl;

	std::vector<Size> loop_edges; //vector that will keep the end points of the vector

	for (Size i = 1; i <= loops.size() ; ++i) {
			loop_edges.push_back(loops[i].start());
			loop_edges.push_back(loops[i].stop());
		}

	std::vector<Size> loop_dividers; //middle residues from the segments that connect the loops
	for (Size i= 2; i < loop_edges.size()-1 ; i += 2  ){ //this for loop is supposed to do an offset start in 2
		int cut = loop_edges[i-1]+((loop_edges[i]-loop_edges[i-1])/2);
		loop_dividers.push_back( cut );
		TR << "Loop edges "<< i << cut << std::endl;
	}


	//getting the secondary structure
	char ss;
	char ss_up;
	char ss_down;

	// defining cut points both out of the loop ranges
	std::vector<Size> cut_points;
	for (Size i=0; i < loop_dividers.size(); ++i ){
		int cutpoint = loop_dividers[i];
		ss = nat_pose.secstruct( cutpoint );

		if (ss == 'L'){
			cut_points.push_back(cutpoint);
			TR<<"Found cut on the middle point between loops: "<< cutpoint <<std::endl;
		}

		else{
			for (Size j=1; j <= 5; ++j ){


				Size cut_point_upset = cutpoint + j;
				Size cut_point_downset = cutpoint - j;

				ss_up = nat_pose.secstruct(cut_point_upset);
				ss_down= nat_pose.secstruct(cut_point_downset);
				bool is_loop_up = is_loop(loops, cut_point_upset);
				bool is_loop_down = is_loop(loops, cut_point_downset);

				if ( ss_up == 'L' && !is_loop_up ){
					cut_points.push_back(cut_point_upset);
					TR<<"Found cut upstream the middle point between loops: "<< cut_point_upset <<std::endl;
					}

				else if (ss_down == 'L' && !is_loop_down){
					cut_points.push_back(cut_point_downset);
					TR<<"Found cut downstream on the middle point between loops: "<<  cut_point_downset <<std::endl;
					}

				else{
					cut_points.push_back(cutpoint);
					TR<<"Found cut on the middle point between loops after screening the vicinity: "<<  cutpoint <<std::endl;
					}

				}

			}

	}

	return cut_points;
}


bool is_loop (
		protocols::loops::Loops & loops,
		Size & residue
){

	bool loop_range = false;
	for (Size i = 1; i <= loops.size() ; ++i) {

		if (loops[i].start() <= residue && loops[i].stop() >= residue){
			loop_range = true;
		}
		else {loop_range = false; }

	}
	return loop_range;

}

bool is_loop_neighbor( protocols::loops::Loops & loops,
		Size & residue,
		Size & range
){
	bool loop_range = false;

	for (Size i = 1; i <= loops.size() ; ++i){

		if (loops[i].start() - range <= residue && loops[i].stop() + range >= residue){
			loop_range =true;

		}
		else{ loop_range = false;}

	}

	return loop_range;
}


void fold_tree_generator(
		protocols::loops::Loops & loops ,
		std::vector<Size> & cutpoints,
		core::pose::Pose & pose,
		kinematics::FoldTree & f
){
	f.add_edge( 1, pose.size(), Edge::PEPTIDE );

	for (Size i=1;  i < loops.size() ; ++i){
		TR<<"LOOP "<<  i <<std::endl;
		TR<<"Cut point"<< cutpoints[i-1]<<std::endl;
		Size cutpoint= cutpoints[i-1];
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_LOWER, cutpoint); // residue on the pose has to be assigned as a cut
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_UPPER, cutpoint+1);
		Size loop1_midpoint = ((loops[1].stop()-loops[1].start())/2) + loops[1].start();
		TR<<"loop 1 mid point"<< loop1_midpoint<<std::endl;
		Size variable_midpoint = ((loops[i+1].stop()-loops[i+1].start())/2) + loops[i+1].start();
		TR<<"Variable mid_point"<< variable_midpoint <<std::endl;
		f.new_jump( loop1_midpoint, variable_midpoint, cutpoint );
	}

}


void define_movemap_extending_chain(
		core::kinematics::MoveMapOP & movemap,
		core::pose::Pose & pose,
		protocols::loops::Loops & loops
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	movemap->set_bb(false);

	Size loop_offset(0);
	for ( Size pos = 1; pos <= pose.size(); pos++ ) {
		bool res_is_loop = is_loop_neighbor( loops, pos, loop_offset);
		if ( res_is_loop ) {
			movemap->set_chi(pos, false);
			continue;
		}
		else {
			pose.set_phi( pos, -150 );
			pose.set_psi( pos, 150);
			pose.set_omega( pos, 180 );
			movemap->set_bb(pos, true); // allowing these residues to move
		}
	}

	if ( option[ OptionKeys::fold_from_loops::loop_mov_nterm ].user() ) { // this particular option will only work properly on the context of one loop
		Size n_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_nterm ];

		for (Size i = 1; i <= loops.size() ; ++i) {

				Size loop_n_pos = loops[i].start();

				for ( Size pos_1 = loop_n_pos; pos_1 < loop_n_pos + n_movable_loop_res; ++pos_1 ){

					TR<<"Movable Residue inside the loop (nterm) "<< pos_1 <<std::endl;
					movemap->set_bb(pos_1, true);

					}
			}

		}


	if ( option[ OptionKeys::fold_from_loops::loop_mov_cterm ].user() ) {

		Size c_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_cterm ];

		for (Size i = 1; i <= loops.size() ; ++i) {

			Size loop_c_pos = loops[i].stop();

				for ( Size pos_2 = loop_c_pos; pos_2 < loop_c_pos + c_movable_loop_res; ++pos_2 ){

						TR<<"Movable Residue inside the loop (cterm) "<< pos_2 <<std::endl;
						movemap->set_bb(pos_2, true);

					}

			}


		}


}


void get_fragments(
		core::fragment::FragSetOP & fragset_large_,
		core::fragment::FragSetOP & fragset_small_
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::fragment;

	std::string frag_large_file, frag_small_file;
		if (option[ in::file::fragA ].user()) frag_large_file  = option[ in::file::fragA ]();
		else                                  frag_large_file  = option[ in::file::frag9 ]();
		if (option[ in::file::fragB ].user()) frag_small_file  = option[ in::file::fragB ]();
		else 								  frag_small_file  = option[ in::file::frag3 ]();

		fragset_large_ = FragmentIO(option[ OptionKeys::abinitio::number_9mer_frags ] ).read_data( frag_large_file );
		fragset_small_ = FragmentIO(option[ OptionKeys::abinitio::number_3mer_frags ] ).read_data( frag_small_file );


}


void extending_chain(
		core::kinematics::MoveMapOP & movemap,
		core::pose::Pose & pose

){

	movemap->set_bb(true);

	for ( Size pos = 1; pos <= pose.size(); pos++ ) {
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );

		}
}

void copying_side_chains(
		core::pose::Pose & nat_pose,
		core::pose::Pose & fold_pose,
		protocols::loops::Loops & loops,
		core::kinematics::MoveMapOP & movemap

){

	for ( Size pos = 1; pos <= fold_pose.size(); pos++ ) {
			bool res_is_loop = is_loop(loops, pos);
			if ( res_is_loop ) {
				fold_pose.replace_residue(pos,nat_pose.residue( pos ), true);
				movemap->set_chi(pos, false);
			}
	}
}


void exclude_loop_residues( core::pose::Pose & pose,
							utility::vector1< bool > & residues_to_mutate,
							utility::vector1< bool > & allowed_aas ,
							core::pack::task::PackerTaskOP & task,
							protocols::loops::Loops & loops
){

	for ( Size pos = 1; pos <= pose.size(); pos++){
		bool res_is_loop = is_loop( loops, pos);
		if (!res_is_loop){
			residues_to_mutate[pos] = true;
			task->nonconst_residue_task(pos).restrict_absent_canonical_aas( allowed_aas );
		}

	}

}


void refresh_cutpoints( core::pose::Pose & pose,
						std::vector<Size> & cut_points){

	for ( std::vector<Size>::iterator it = cut_points.begin(); it != cut_points.end(); ++it )
	{

		TR << "Refresh cut_point "<< *it <<std::endl;
		Size residue = *it;
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_LOWER,residue ); // residue on the pose has to be assigned as a cut
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_UPPER, residue + 1);

	}


}


bool is_cut( std::vector<Size> & cut_points,
			Size & residue){
	bool res_cut = false;
	for ( std::vector<Size>::iterator it = cut_points.begin(); it != cut_points.end(); ++it )
		{
			if (*it == residue){ res_cut = true; }

		}
	return res_cut;
}

void CA_cst_generator(core::pose::Pose & pose,
					scoring::constraints::ConstraintSetOP & cst,
					protocols::loops::Loops & loops,
					std::vector<Size> & cut_points
){

using namespace scoring::constraints;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;

	Real sd = option[ OptionKeys::fold_from_loops::ca_csts_dev ];

	TR << "Constrains Standard Deviation - "<< sd <<std::endl;

	for ( Size pos = 1; pos <= pose.size(); ++pos ) {
				for (Size pos_2 = 1; pos_2 <=pose.size(); ++pos_2   ){

					bool res_is_loop = is_loop(loops, pos);
					bool res2_is_loop = is_loop(loops,pos_2);
					bool cut_point_pos = is_cut(cut_points, pos );
					bool cut_point_pos_2 = is_cut(cut_points, pos_2 );

					//Size offset(5);
					//bool res_is_loop_neighbor = is_loop_neighbor(loops,pos, offset);
					//bool res2_is_loop_neighbor = is_loop_neighbor(loops,pos_2, offset);

					Size seq_sep = 0;

					if (pos > pos_2){
						seq_sep = pos - pos_2;
						}
					else { seq_sep = pos_2 - pos; }


					if ( seq_sep >=6  ){

					if ( !res_is_loop && !res2_is_loop ) {

						if ( !cut_point_pos && !cut_point_pos_2){
							if (pos != pos_2) {

								core::conformation::Residue res_pos = pose.residue(pos);
								core::conformation::Residue res_pos_2 = pose.residue(pos_2);

								Real const d( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
								TR <<pos_2 <<" "<< pos << " "<<d <<std::endl;


								cst->add_constraint( new AtomPairConstraint ( AtomID(res_pos.atom_index("CA"),pos), AtomID(res_pos_2.atom_index("CA"),pos_2), new HarmonicFunc( d, sd ) ) );

							}

						}
					}

					}
				}

	}

}


void CA_cst_generator(core::pose::Pose & pose,
					scoring::constraints::ConstraintSetOP & cst,
					protocols::loops::Loops & loops
){

using namespace scoring::constraints;
using namespace id;
using namespace basic::options;
using namespace basic::options::OptionKeys;

	Real sd = option[ OptionKeys::fold_from_loops::ca_csts_dev ];

	TR << "Constrains Standard Deviation - "<< sd <<std::endl;


	for ( Size pos = 1; pos <= pose.size(); ++pos ) {
				for (Size pos_2 = 1; pos_2 <=pose.size(); ++pos_2   ){

					bool res_is_loop = is_loop(loops, pos);
					bool res2_is_loop = is_loop(loops,pos_2);

					Size offset(5);
					bool res_is_loop_neighbor = is_loop_neighbor(loops,pos, offset);
					bool res2_is_loop_neighbor = is_loop_neighbor(loops,pos_2, offset);

					Size seq_sep = 0;

					if (pos > pos_2){
						seq_sep = pos - pos_2;

					}
					else { seq_sep = pos_2 - pos; }


					if ( seq_sep >= 6  ){

					if ( !res_is_loop && !res2_is_loop ) {

							if ( !res_is_loop_neighbor && !res2_is_loop_neighbor) {

								core::conformation::Residue res_pos = pose.residue(pos);
								core::conformation::Residue res_pos_2 = pose.residue(pos_2);

								Real const d( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
								TR <<pos_2 <<" "<< pos << " "<<d <<std::endl;


								cst->add_constraint( new AtomPairConstraint ( AtomID(res_pos.atom_index("CA"),pos), AtomID(res_pos_2.atom_index("CA"),pos_2), new HarmonicFunc( d, sd ) ) );

							}


					}

					}
				}

	}

}


void new_pose_generator(core::pose::Pose & target_loops, core::pose::Pose & nat_prot ,protocols::loops::Loops & loops){

	using namespace core::conformation;

	core::pose::Pose centroid_target_loops;


	core::util::switch_to_residue_type_set( target_loops , core::chemical::CENTROID );

	core::chemical::ResidueTypeSet const & rsd_set( target_loops.residue(1).residue_type_set() );

	std::string nat_seq = nat_prot.sequence();


	centroid_target_loops = target_loops;

	if (loops.size() == 1){


		Size nsegment = loops[1].start()-1;
		Size csegment = nat_prot.size() - loops[1].stop();

		TR << "NSEGMENT "<< nsegment <<std::endl;
		TR << "CSEGMENT "<< csegment <<std::endl;
		TR << "NAT SEQ  " << nat_seq <<std::endl;


		for (Size k=1; k <= nsegment ; ++k ) {


			const char aa = nat_seq[ loops[1].start() - k - 1]; // The minus 1 is here to compensate the indexing strating in zero
			Size residue = loops[1].start() - k;


			TR << "RES AA   " << residue << aa <<std::endl;


			core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueTypeSelector().set_name1( aa ).exclude_variants().select( rsd_set )[1] );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, 1, true);


			target_loops.set_omega( 1, 180.0 );


		}


		for ( Size j = 0 ; j < csegment ; ++j  ){

			const char aa = nat_seq[ loops[1].stop() + j  ];
			Size residue = loops[1].stop() + j;


			core::chemical::ResidueTypeCOP new_rsd_type( core::chemical::ResidueTypeSelector().set_name1( aa ).exclude_variants().select( rsd_set )[1] );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, residue , true );


			target_loops.set_omega( residue , 180.0 );


		}


	}


	TR << "Last Residue    " << target_loops.size() <<std::endl;

	core::pose::add_variant_type_to_pose_residue(target_loops, chemical::LOWER_TERMINUS, 1 );
	core::pose::add_variant_type_to_pose_residue(target_loops, chemical::UPPER_TERMINUS, target_loops.size());

}


void copying_side_chains_swap_loop (
		core::pose::Pose & swap_loops,
		core::pose::Pose & fold_pose,
		protocols::loops::Loops & loops,
		core::kinematics::MoveMapOP & movemap

){

	using namespace core::conformation;
	Size offsetres = 0;

	// Cleanin up the pose - > this step has to be performed before
	core::pose::remove_lower_terminus_type_from_pose_residue(swap_loops, 1 );
	core::pose::remove_upper_terminus_type_from_pose_residue(swap_loops, swap_loops.size());


	movemap->set_chi(true);


	for ( Size pos = 1; pos <= fold_pose.size(); ++pos ) {
				bool res_is_loop = is_loop(loops, pos);

				if ( res_is_loop ) {

					++offsetres ;


					//core::chemical::ResidueTypeSet const & rsd_set( swap_loops.residue( offsetres ).residue_type_set() );

					//core::chemical::ResidueType const & rsd_type( rsd_set.name_map( swap_loops.residue( offsetres ).name3() ) );


					//replace_pose_residue_copying_existing_coordinates( fold_pose, pos, rsd_type ); //this function is not really keeping the same rotamers
																									// potential fix is to idealize and fill missing atoms in our
																									// input loops and then copy the coordinates


					fold_pose.replace_residue( pos, swap_loops.residue( offsetres ), true);

					TR << "After Swap loop Residue  " << offsetres << std::endl;

					//Residue const & old_rsd( fold.residue( seqpos ) );

					//copy_residue_coordinates_and_rebuild_missing_atoms( swap_loops.residue( offsetres ), old_rsd , fold_pose.conformation() );


					movemap->set_chi(pos, false);


				}
	}

}


////////////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char* argv [] )
{
	try {

	protocols::abinitio::ClassicAbinitio::register_options();
	protocols::abinitio::AbrelaxApplication::register_options();
	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::fragment;
	using namespace protocols;
	using namespace constraints;


	pose::Pose nat_pose;

	pose::Pose extended_pose;


	core::import_pose::pose_from_file( nat_pose, basic::options::start_file() , core::import_pose::PDB_file);

	extended_pose = nat_pose; //making working copy

	protocols::checkpoint::CheckPointer sliding_checkpoint("closing");


	//Reading loops
	protocols::loops::Loops lr_loops_in;
	if ( option[ OptionKeys::loops::loop_file ].user() ){
		std::string filename( protocols::loops::get_loop_file_name() );
		lr_loops_in.read_loop_file( filename );
	}


	ConstraintSetOP ca_cst ( new ConstraintSet() );


	//Defining a fold tree
	kinematics::FoldTree f;
	f.clear();

	std::vector<Size> cut_points;

	if (lr_loops_in.size() > 1){

		cut_points= define_cut_points (lr_loops_in, nat_pose);
		fold_tree_generator(lr_loops_in, cut_points, extended_pose, f);
		TR << "Pose fold tree " << f << std::endl;


		if (option [OptionKeys::fold_from_loops::native_ca_cst ].user() ){


			CA_cst_generator( nat_pose, ca_cst ,lr_loops_in, cut_points );

		}

	}

	else{


		TR<< "Only one loop defined no cutpoints necessary "<<std::endl;
		f.add_edge( 1, extended_pose.size(), Edge::PEPTIDE );
		TR << "Pose fold tree " << f << std::endl;

		if (option [OptionKeys::fold_from_loops::native_ca_cst ].user() ){

				CA_cst_generator( nat_pose, ca_cst ,lr_loops_in );

			}


	}


	extended_pose.fold_tree( f ); // apply the fold tree

	core::util::switch_to_residue_type_set( extended_pose, core::chemical::CENTROID );


	if ( option[out::pdb].user()){
		extended_pose.dump_pdb("init_centroid.pdb");
		}


	core::scoring::ScoreFunctionOP scorefxn_centr = core::scoring::ScoreFunctionFactory::create_score_function( "score3" ); // defining scoring function

	scorefxn_centr->set_weight( scoring::linear_chainbreak, 1);


	core::io::silent::SilentFileData sfd_cent;
	core::io::silent::SilentFileData sfd_fa;
	core::io::silent::SilentFileData sfd_des;
	std::string pdb_silent_file_cent( option[ out::file::silent ]() );
	std::string pdb_silent_file_fa(option[ out::file::silent ]()+"_fa");
	std::string pdb_silent_file_des(option [ out::file::silent]()+"_des");


	Real loop_tag=0;

	if ( option[ out::file::silent ].user() ){


		Real native_rmsd= core::scoring::rmsd_with_super(nat_pose, extended_pose, is_protein_CA); // calculating RMSD on CA


		pose::setPoseExtraScore( extended_pose, "rms",  native_rmsd );
		(*scorefxn_centr)(extended_pose); // score of the native structure


				core::io::silent::BinarySilentStructOP ss_init ( new core::io::silent::BinarySilentStruct(extended_pose,"NATIVE" ));


			sfd_cent.add_structure(ss_init);

	}


	core::kinematics::MoveMapOP  movemap = new core::kinematics::MoveMap(); // defining move map


	if (option[ OptionKeys::loops::loop_file ].user()){

		define_movemap_extending_chain( movemap, extended_pose, lr_loops_in); // extending the chain and setting the psi-phis moves
	}

	else{

		extending_chain( movemap, extended_pose );

		TR << "No LOOPS FILE defined -> Classic Abinitio  " << std::endl;
	}


	if ( option[out::pdb].user()){
			extended_pose.dump_pdb("extended_pose.pdb");
		}


	FragSetOP fragset_large_;
	FragSetOP fragset_small_;

	get_fragments(fragset_large_, fragset_small_ ); // reading_fragments


	core::pose::Pose target_loops;
	core::pose::Pose nat_target_loops;
	core::pose::Pose extended_target_loops;


	if (option [OptionKeys::fold_from_loops::swap_loops ].user()){


		std::string swap_loops = option [OptionKeys::fold_from_loops::swap_loops ]().name();

		core::import_pose::pose_from_file( target_loops , swap_loops , core::import_pose::PDB_file);

		nat_target_loops = target_loops; // keep copy to recover side chains later


		new_pose_generator( target_loops, nat_pose , lr_loops_in ); //takes the target_loops an builds a pose with with the new loop and the correct lengths


		extended_pose = target_loops;

		extended_pose.fold_tree( f );

	}


	if (option [OptionKeys::fold_from_loops::native_ca_cst ].user() ){
		extended_pose.constraint_set( ca_cst ); //cst stuff
	}


	if (option[ OptionKeys::in::file::psipred_ss2 ].user() ) {


		protocols::loops::set_secstruct_from_psipred_ss2(extended_pose);


		TR << "Pose Secondary Structure Set by Psipred " << extended_pose.secstruct() << std::endl;


	}


	core::scoring::ScoreFunctionOP scorefxn = core::scoring::ScoreFunctionFactory::create_score_function( "cen_std" ,"score4L" ); //setting up scoring function for loop closure

	scorefxn->set_weight( scoring::linear_chainbreak, 1); //scoring function to close the loops has a chainbreak term ... makes sense


	pose::Pose fold_pose; //defining a pose o be folded

	pose::Pose fold_pose_relax;


	Size nstruct = option[ out::nstruct ];

	Size i = 0;


	load_checkpoint( i );


	for (; i < nstruct ; ++i  ){

		fold_pose = extended_pose; //passing the extended pose

		std::string outfile_extended = option[ out::prefix ]() + "_" + right_string_of(i,3,'0') ; //defining name for the extended pose

		protocols::abinitio::ClassicAbinitioOP abinitio = new protocols::abinitio::ClassicAbinitio( fragset_small_, fragset_large_ ,movemap );

		if (option [OptionKeys::fold_from_loops::native_ca_cst].user() ){

		abinitio = new protocols::abinitio::FoldConstraints( fragset_small_, fragset_large_ ,movemap );

			}


		abinitio->init(fold_pose);
		abinitio->apply( fold_pose );


		TR << "Pose Secondary Structure After abinitio " << fold_pose.secstruct() << std::endl;


		Real rmsd_to_native =  core::scoring::rmsd_with_super(fold_pose, nat_pose, is_protein_CA);


		//Insert code to save the abinitio decoys here

		pose::setPoseExtraScore( fold_pose, "rms",   rmsd_to_native );
		(*scorefxn_centr)(fold_pose); // score of the native structure


		core::io::silent::BinarySilentStructOP ss_cent ( new core::io::silent::BinarySilentStruct(fold_pose, outfile_extended ));


		sfd_cent.add_structure(ss_cent);


		if ( rmsd_to_native < 5.0) {


			if ( option[ OptionKeys::loops::ccd_closure ].user()  ){ // to close loops

				Real chain_break_dist = 0; // chain_break eval

				chain_break_dist = fold_pose.energies().total_energies()[ scoring::linear_chainbreak ];


				TR << "Chain Break -  " << chain_break_dist << std::endl;


				protocols::loops::SlidingWindowLoopClosureOP closure_protocol = new protocols::loops::SlidingWindowLoopClosure( fragset_small_, scorefxn, movemap );


				closure_protocol->scored_frag_cycle_ratio( 0.2 );
				closure_protocol->short_frag_cycle_ratio( 0.1 );

			// TODO: fix the new checkpointer issue
				loop_tag = 1;
				try{
					jumping::close_chainbreaks( closure_protocol, fold_pose, sliding_checkpoint ,"sliding", f );
				} catch ( protocols::loops::EXCN_Loop_not_closed& excn ) {
					loop_tag = 0;
				}

				TR << "LOOPS_TAG " << loop_tag << std::endl;


			}


		std::string outfilename = option[ out::prefix ]() + "_" + right_string_of(i,3,'0') + ".pdb.gz"; //building name


		core::util::switch_to_residue_type_set( fold_pose, chemical::FA_STANDARD ); //switching the pose for full atom


		core::scoring::ScoreFunctionOP scorefxn_fa( get_score_function() );
		scorefxn_fa->set_weight( scoring::chainbreak, 1 ); // in order to pipe this scoring function  maybe I have to refresh the cutpoints
		scorefxn_fa->set_weight(scoring::overlap_chainbreak, 1 );


		if (lr_loops_in.size() > 1 ){

				refresh_cutpoints(fold_pose, cut_points);

		}


		//Copying Side Chains from loops

		if ( option [OptionKeys::fold_from_loops::swap_loops ].user() ){

			copying_side_chains_swap_loop( nat_target_loops, fold_pose, lr_loops_in, movemap );


		}

		else {

			copying_side_chains( nat_pose, fold_pose, lr_loops_in, movemap);
		}

		//relaxing structure


		fold_pose_relax = fold_pose;


		(*scorefxn_fa)(fold_pose_relax);


		relax::ClassicRelax relax_protocol( scorefxn_fa, movemap);

		//Implementing relax cycles


		relax_protocol.apply( fold_pose_relax );


		std::string outfilename_rlx = outfile_extended + "_nat_rlx";

		std::string outfilename_des = outfile_extended + "_des";

		std::string outfilename_rlx_1 = outfile_extended +"_des_rlx_1";

		std::string outfilename_rlx_2 = outfile_extended +"_des_rlx_2";


		core::io::silent::BinarySilentStructOP ss_fa_nat_rlx ( new core::io::silent::BinarySilentStruct( fold_pose_relax, outfilename_rlx ));


		sfd_fa.add_structure( ss_fa_nat_rlx );


		if ( option [OptionKeys::fold_from_loops::add_relax_cycles ].user() ){

			Size n_relax = option [OptionKeys::fold_from_loops::add_relax_cycles ] ;

			for (i=1; i <= n_relax ; ++i ){

				std::string outfilename_n_rlx = outfilename_rlx + "_" + right_string_of(i,1,'0') ;

				relax_protocol.apply( fold_pose_relax );

				core::io::silent::BinarySilentStructOP ss_fa_nat_rlx_cycle ( new core::io::silent::BinarySilentStruct( fold_pose_relax, outfilename_n_rlx ));

				sfd_fa.add_structure( ss_fa_nat_rlx_cycle );
			}

		}


		//Design step


		if (option[ OptionKeys::loops::loop_file ].user()){

			utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, true );

			utility::vector1< bool > residues_to_mutate( fold_pose.size(), false );

			core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( fold_pose ));


			exclude_loop_residues(fold_pose, residues_to_mutate, allowed_aas, task, lr_loops_in );

			task->restrict_to_residues( residues_to_mutate );

			task->initialize_extra_rotamer_flags_from_command_line();


			pack::pack_rotamers( fold_pose, *scorefxn_fa , task );


			//save the designed decoy

			core::io::silent::BinarySilentStructOP ss_fa_des ( new core::io::silent::BinarySilentStruct( fold_pose, outfilename_des ));


			sfd_des.add_structure( ss_fa_des );


			relax_protocol.apply(fold_pose);


			Real final_rmsd  = core::scoring::rmsd_with_super(nat_pose, fold_pose, is_protein_CA); // calculating RMSD on CA

			pose::setPoseExtraScore( fold_pose, "rms",  final_rmsd );


			(*scorefxn_fa)(fold_pose); // score of the fold_pose


			core::io::silent::BinarySilentStructOP ss_fa_rlx ( new core::io::silent::BinarySilentStruct( fold_pose, outfilename_rlx_1 ));


			sfd_fa.add_structure( ss_fa_rlx );


			pack::pack_rotamers( fold_pose, *scorefxn_fa , task );


			relax_protocol.apply( fold_pose );


			Real final_rmsd_rlx  = core::scoring::rmsd_with_super(nat_pose, fold_pose, is_protein_CA); // calculating RMSD on CA

			pose::setPoseExtraScore( fold_pose, "rms",  final_rmsd_rlx );

			core::io::silent::BinarySilentStructOP ss_fa_rlx_2 ( new core::io::silent::BinarySilentStruct( fold_pose, outfilename_rlx_2 ));

			sfd_fa.add_structure( ss_fa_rlx_2 );

		}


		sfd_fa.write_all( pdb_silent_file_fa );
		sfd_des.write_all(pdb_silent_file_des );


		sfd_fa.clear();
		sfd_des.clear();


		if ( option[out::pdb].user()){
			utility::io::ozstream outfile( outfilename ); // outputing gzip pdbs
			core::io::pdb::dump_pdb( fold_pose, outfile );
			outfile.close();
		}


	}


	write_checkpoint( i );

	}


	sfd_cent.write_all(pdb_silent_file_cent);

	TR.flush();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;

}

