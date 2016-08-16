// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/fold_from_loops/FoldFromLoops_functions.cc
/// @brief  May 1, 2009
/// @author bcorreia


//Unite headers
#include <devel/fold_from_loops/FoldFromLoops_functions.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/conformation/ResidueFactory.hh>


#include <core/kinematics/FoldTree.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <core/pose/Pose.hh>

#include <utility/vector1.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

//constraints
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>

//options stuff

#include <basic/options/option.hh>
#include <basic/options/keys/fold_from_loops.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


//design headers

#include <core/pack/task/ResfileReader.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>


#include <vector>
#include <iostream>
#include <string>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector0.hh>


namespace devel {
namespace fold_from_loops {

using namespace core;
using namespace kinematics;
using namespace scoring::constraints;
using namespace core::fragment;


static THREAD_LOCAL basic::Tracer TR( "FoldFromLoops_function" );


bool is_loop (
	protocols::loops::Loops & loops,
	Size & residue
){

	bool loop_range = false;
	for ( Size i = 1; i <= loops.size() ; ++i ) {

		if ( loops[i].start() <= residue && loops[i].stop() >= residue ) {
			loop_range = true;
		} else { loop_range = false; }

	}
	return loop_range;

}


std::vector<Size> define_cut_points (

	protocols::loops::Loops & loops,
	core::pose::Pose & nat_pose
){


	TR <<"Several Loops Defined "<<std::endl;

	std::vector<Size> loop_edges; //vector that will keep the end points of the vector

	for ( core::Size i = 1; i <= loops.size() ; ++i ) {
		loop_edges.push_back(loops[i].start());
		loop_edges.push_back(loops[i].stop());
	}

	std::vector<Size> loop_dividers; //middle residues from the segments that connect the loops
	for ( Size i= 2; i < loop_edges.size()-1 ; i += 2  ) { //this for loop is supposed to do an offset start in 2
		int cut = loop_edges[i-1]+((loop_edges[i]-loop_edges[i-1])/2);
		loop_dividers.push_back( cut );
		TR << "Loop edges "<< i << cut << std::endl;
	}


	//getting the secondary structure
	char ss_up;
	char ss_down;

	// defining cut points both out of the loop ranges
	std::vector<Size> cut_points;
	for ( Size i=0; i < loop_dividers.size(); ++i ) {
		int cutpoint = loop_dividers[i];
		char ss = nat_pose.secstruct( cutpoint );

		if ( ss == 'L' ) {
			cut_points.push_back(cutpoint);
			TR<<"Found cut on the middle point between loops: "<< cutpoint <<std::endl;
		} else {
			for ( Size j=1; j <= 5; ++j ) {


				Size cut_point_upset = cutpoint + j;
				Size cut_point_downset = cutpoint - j;

				ss_up = nat_pose.secstruct(cut_point_upset);
				ss_down= nat_pose.secstruct(cut_point_downset);
				bool is_loop_up = is_loop(loops, cut_point_upset);
				bool is_loop_down = is_loop(loops, cut_point_downset);

				if ( ss_up == 'L' && !is_loop_up ) {
					cut_points.push_back(cut_point_upset);
					TR<<"Found cut upstream the middle point between loops: "<< cut_point_upset <<std::endl;
				} else if ( ss_down == 'L' && !is_loop_down ) {
					cut_points.push_back(cut_point_downset);
					TR<<"Found cut downstream on the middle point between loops: "<<  cut_point_downset <<std::endl;
				} else {
					cut_points.push_back(cutpoint);
					TR<<"Found cut on the middle point between loops after screening the vicinity: "<<  cutpoint <<std::endl;
				}

			}

		}

	}

	return cut_points;
}


void fold_tree_cutpoints_generator(
	protocols::loops::Loops & loops ,
	std::vector<Size> & cutpoints,
	core::pose::Pose & pose,
	kinematics::FoldTree & f
){

	using namespace core;
	using namespace kinematics;


	if ( loops.size() > 1 ) {

		cutpoints = define_cut_points( loops, pose );

		f.add_edge( 1, pose.total_residue(), Edge::PEPTIDE );


		for ( Size i=1;  i < loops.size() ; ++i ) {
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

		TR << "Fold Tree of the initial scaff " << f << std::endl;

	} else {

		TR<< "Only one loop defined no cutpoints necessary "<<std::endl;
		f.add_edge( 1, pose.total_residue(), Edge::PEPTIDE );
		TR << "Pose fold tree " << f << std::endl;

	}


}


bool is_cut( std::vector<Size> & cut_points,
	Size & residue){
	bool res_cut = false;
	for ( std::vector<Size>::iterator it = cut_points.begin(); it != cut_points.end(); ++it ) {
		if ( *it == residue ) { res_cut = true; }

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

	protocols::loops::Loops clear_loops_seg;

	if ( option[ OptionKeys::fold_from_loops::clear_csts].user() ) {

		std::string  clear_cst_file = option[ OptionKeys::fold_from_loops::clear_csts]();

		clear_loops_seg = protocols::loops::Loops(clear_cst_file);

	}


	for ( Size pos = 1; pos <= pose.total_residue(); ++pos ) {
		for ( Size pos_2 = 1; pos_2 <=pose.total_residue(); ++pos_2   ) {

			bool res_is_loop = loops.is_loop_residue( pos );
			bool res2_is_loop = loops.is_loop_residue( pos_2 );
			bool cut_point_pos = false;
			bool cut_point_pos_2 = false;

			bool res_cst_free = false;


			if ( option[ OptionKeys::fold_from_loops::clear_csts].user() ) {

				if ( clear_loops_seg.is_loop_residue(pos) || clear_loops_seg.is_loop_residue(pos_2) ) {

					res_cst_free = true;
				}
			}


			if ( !cut_points.empty() ) {

				cut_point_pos = is_cut(cut_points, pos );
				cut_point_pos_2 = is_cut(cut_points, pos_2 );


			}


			Size seq_sep = 0;

			if ( pos > pos_2 ) {
				seq_sep = pos - pos_2;
			} else { seq_sep = 0; } // This is here to avoid repetition of the csts


			if ( seq_sep >=6  ) {

				if ( !res_is_loop && !res2_is_loop ) {

					if ( !cut_point_pos && !cut_point_pos_2 ) {
						if ( pos != pos_2 ) {
							if ( !res_cst_free ) {
								core::conformation::Residue res_pos = pose.residue(pos);
								core::conformation::Residue res_pos_2 = pose.residue(pos_2);

								Real const d( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
								TR <<pos_2 <<" "<< pos << " "<<d <<std::endl;


								cst->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint ( AtomID(res_pos.atom_index("CA"),pos), AtomID(res_pos_2.atom_index("CA"),pos_2), core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( d, sd ) ) ) ) ) );
							}
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

	//define here a place were to read the free csts segments
	// create the option

	protocols::loops::Loops clear_loops_seg;

	if ( option[ OptionKeys::fold_from_loops::clear_csts].user() ) {

		std::string  clear_cst_file = option[ OptionKeys::fold_from_loops::clear_csts]();

		clear_loops_seg = protocols::loops::Loops(clear_cst_file);

	}


	for ( Size pos = 1; pos <= pose.total_residue(); ++pos ) {
		for ( Size pos_2 = 1; pos_2 <=pose.total_residue(); ++pos_2   ) {

			bool res_is_loop = loops.is_loop_residue( pos );
			bool res2_is_loop = loops.is_loop_residue( pos_2 );

			int offset(5);
			bool res_is_loop_neighbor = loops.is_loop_residue(pos, offset);
			bool res2_is_loop_neighbor = loops.is_loop_residue(pos_2, offset);

			bool res_cst_free = false;
			if ( option [OptionKeys::fold_from_loops::add_cst_loop].user() ) {
				res_is_loop = false;
				res2_is_loop = false;
				res_is_loop_neighbor = false;
				res2_is_loop_neighbor = false;
			}

			if ( option[ OptionKeys::fold_from_loops::clear_csts].user() ) {

				if ( clear_loops_seg.is_loop_residue(pos) || clear_loops_seg.is_loop_residue(pos_2) ) {

					res_cst_free = true;
				}
			}


			Size seq_sep = 0;

			if ( pos > pos_2 ) {
				seq_sep = pos - pos_2;

			} else {
				seq_sep = pos_2 - pos;
			}


			if ( seq_sep >= 6  ) {

				if ( !res_is_loop && !res2_is_loop ) {

					if ( !res_is_loop_neighbor && !res2_is_loop_neighbor ) {

						if ( !res_cst_free ) {
							core::conformation::Residue res_pos = pose.residue(pos);
							core::conformation::Residue res_pos_2 = pose.residue(pos_2);

							Real const d( res_pos.xyz( res_pos.atom_index("CA") ).distance( res_pos_2.xyz( res_pos_2.atom_index("CA") )));
							TR <<pos_2 <<" "<< pos << " "<<d <<std::endl;


							cst->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint ( AtomID(res_pos.atom_index("CA"),pos), AtomID(res_pos_2.atom_index("CA"),pos_2), core::scoring::func::FuncOP( new core::scoring::func::HarmonicFunc( d, sd ) ) ) ) ) );

						}
					}
				}
			}
		}
	}


}


bool is_loop_neighbor( protocols::loops::Loops & loops,
	Size & residue,
	Size & range
){
	bool loop_range = false;

	for ( Size i = 0; i <= loops.size() ; ++i ) {

		if ( loops[i].start() - range <= residue && loops[i].stop() + range >= residue ) {
			loop_range =true;

		} else { loop_range = false; }

	}

	return loop_range;
}


void define_movemap_extending_chain(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose,
	protocols::loops::Loops & loops
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	movemap->set_bb(false);


	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		if ( loops.is_loop_residue( pos ) ) {
			movemap->set_chi(pos, false);
			continue;
		} else {
			pose.set_phi( pos, -150 );
			pose.set_psi( pos, 150);
			pose.set_omega( pos, 180 );
			movemap->set_bb(pos, true); // allowing these residues to move
		}
	}

	if ( option[ OptionKeys::fold_from_loops::loop_mov_nterm ].user() ) { // this particular option will only work properly on the context of one loop
		Size n_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_nterm ];

		for ( Size i = 1; i <= loops.size() ; ++i ) {

			Size loop_n_pos = loops[i].start();

			for ( Size pos_1 = loop_n_pos; pos_1 < loop_n_pos + n_movable_loop_res; ++pos_1 ) {

				TR<<"Movable Residue inside the loop (nterm) "<< pos_1 <<std::endl;
				movemap->set_bb(pos_1, true);

			}
		}

	}


	if ( option[ OptionKeys::fold_from_loops::loop_mov_cterm ].user() ) {

		Size c_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_cterm ];

		for ( Size i = 1; i <= loops.size() ; ++i ) {

			Size loop_c_pos = loops[i].stop();

			for ( Size pos_2 = loop_c_pos; pos_2 < loop_c_pos + c_movable_loop_res; ++pos_2 ) {

				TR<<"Movable Residue inside the loop (cterm) "<< pos_2 <<std::endl;
				movemap->set_bb(pos_2, true);

			}

		}


	}


}


void define_movemap(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose,
	protocols::loops::Loops & loops
){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	movemap->set_bb(false);


	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {


		if ( loops.is_loop_residue( pos ) ) {
			movemap->set_chi(pos, false);
			TR<< "residues that don't move" << pos << std::endl;
			continue;
		} else {
			movemap->set_bb(pos, true); // allowing these residues to move
		}
	}

	if ( option[ OptionKeys::fold_from_loops::loop_mov_nterm ].user() ) { // this particular option will only work properly on the context of one loop
		Size n_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_nterm ];

		for ( Size i = 1; i <= loops.size() ; ++i ) {

			Size loop_n_pos = loops[i].start();

			for ( Size pos_1 = loop_n_pos; pos_1 < loop_n_pos + n_movable_loop_res; ++pos_1 ) {

				TR<<"Movable Residue inside the loop (nterm) "<< pos_1 <<std::endl;
				movemap->set_bb(pos_1, true);

			}
		}

	}


	if ( option[ OptionKeys::fold_from_loops::loop_mov_cterm ].user() ) {

		Size c_movable_loop_res = option[OptionKeys::fold_from_loops::loop_mov_cterm ];

		for ( Size i = 1; i <= loops.size() ; ++i ) {

			Size loop_c_pos = loops[i].stop();

			for ( Size pos_2 = loop_c_pos; pos_2 < loop_c_pos + c_movable_loop_res; ++pos_2 ) {

				TR<<"Movable Residue inside the loop (cterm) "<< pos_2 <<std::endl;
				movemap->set_bb(pos_2, true);

			}

		}


	}


}


void extending_chain(
	core::kinematics::MoveMapOP & movemap,
	core::pose::Pose & pose

){

	movemap->set_bb(true);

	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		pose.set_phi( pos, -150 );
		pose.set_psi( pos, 150);
		pose.set_omega( pos, 180 );

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
	if ( option[ in::file::fragA ].user() ) frag_large_file  = option[ in::file::fragA ]();
	else                                  frag_large_file  = option[ in::file::frag9 ]();
	if ( option[ in::file::fragB ].user() ) frag_small_file  = option[ in::file::fragB ]();
	else           frag_small_file  = option[ in::file::frag3 ]();

	fragset_large_ = FragmentIO(option[ OptionKeys::abinitio::number_9mer_frags ] ).read_data( frag_large_file );
	fragset_small_ = FragmentIO(option[ OptionKeys::abinitio::number_3mer_frags ] ).read_data( frag_small_file );


}


void new_pose_generator(
	core::pose::Pose & target_loops,
	core::pose::Pose & nat_prot,
	protocols::loops::Loops & loops,
	std::vector<Size> & cut_points
){

	using namespace core::conformation;

	core::pose::Pose centroid_target_loops;


	core::util::switch_to_residue_type_set( target_loops , core::chemical::CENTROID );

	core::chemical::ResidueTypeSetCOP rsd_set( target_loops.residue(1).residue_type_set() );

	std::string nat_seq = nat_prot.sequence();


	centroid_target_loops = target_loops;

	if ( loops.size() == 1 ) {


		Size nsegment = loops[1].start()-1;
		Size csegment = nat_prot.total_residue() - loops[1].stop();

		TR << "NSEGMENT "<< nsegment <<std::endl;
		TR << "CSEGMENT "<< csegment <<std::endl;
		TR << "NAT SEQ  " << nat_seq <<std::endl;


		//set up a jump between the two loops


		for ( Size k=1; k <= nsegment ; ++k ) {


			const char aa = nat_seq[ loops[1].start() - k - 1]; // The minus 1 is here to compensate the indexing starting in zero
			Size residue = loops[1].start() - k;


			TR << "RES AA   " << residue << aa <<std::endl;


			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, 1, true);


			target_loops.set_omega( 1, 180.0 );


		}


		for ( Size j = 0 ; j < csegment ; ++j  ) {

			const char aa = nat_seq[ loops[1].stop() + j  ];
			Size residue = loops[1].stop() + j;



			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, residue , true );


			target_loops.set_omega( residue , 180.0 );


		}


	}


	//here is were we going to add a function to deal with the case of 2 loops


	if ( loops.size() == 2 ) {


		FoldTree loop_f;

		loop_f.add_edge( 1, target_loops.total_residue(), Edge::PEPTIDE );

		loop_f.new_jump( 1,target_loops.total_residue()-1, 1+(loops[1].stop()-loops[1].start()) );


		target_loops.fold_tree( loop_f );


		Size n1segment = loops[1].start()-1;
		Size c1segment = cut_points[0]-loops[1].stop();


		Size n2segment = loops[2].start()-cut_points[0];
		Size c2segment = nat_prot.total_residue() - loops[2].stop();


		TR << "NSEGMENT1 "<< n1segment <<std::endl;
		TR << "CSEGMENT1 "<< c1segment <<std::endl;

		TR << "NSEGMENT2 "<< n2segment <<std::endl;
		TR << "CSEGMENT2 "<< c2segment <<std::endl;

		TR << "CUT "<< cut_points[0] <<std::endl;

		TR << "NAT SEQ  " << nat_seq <<std::endl;


		for ( Size k=1; k <= n1segment ; ++k ) {


			const char aa = nat_seq[ loops[1].start() - k - 1]; // The minus 1 is here to compensate the indexing starting in zero
			Size residue = loops[1].start() - k;


			TR << "RES AA N1  " << residue << aa <<std::endl;


			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );



			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, 1, true);


			target_loops.set_omega( 1, 180.0 );


		}

		//target_loops.dump_pdb("addedN1element2loops.pdb");

		for ( Size j = 0 ; j < c1segment ; ++j  ) {

			const char aa = nat_seq[ loops[1].stop() + j  ];
			Size residue = loops[1].stop() + j ;


			TR << "RES AA C1  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, residue , true );


			target_loops.set_omega( residue , 180.0 );


		}


		Size current_tot_res = target_loops.total_residue();

		TR << " Target Loops total residues before the second loop inserts  " << target_loops.total_residue() << std::endl;


		for ( Size l = 1; l < n2segment ; ++l  ) {


			const char aa = nat_seq[ loops[2].start() - l - 1  ];// need to check this one
			Size residue =   current_tot_res - (loops[2].stop()-loops[2].start()) ; // we need to check which one we start


			TR << "RES AA N2  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}


		for ( Size m = 0; m < c2segment ; ++m  ) {

			const char aa = nat_seq[ loops[2].stop() + m  ];
			Size residue = loops[2].stop() + m ;


			TR << "RES AA C2  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}


	}


	if ( loops.size() == 3 ) {

		FoldTree loop_f;

		loop_f.add_edge( 1, target_loops.total_residue(), Edge::PEPTIDE );

		//loop_f.new_jump( 1,target_loops.total_residue()-1, 1+(loops[3].stop()-loops[3].start()) );

		loop_f.new_jump(1, target_loops.total_residue()-1,  target_loops.total_residue()-loops[3].size());

		TR <<"inside loop_3 build "<< std::endl;

		//>loop_f.new_jump( 1,  loops[1].size()+2 , loops[1].size()+1 );

		loop_f.new_jump( 1,  loops[1].size()+2 , loops[1].size() );


		//loop_f.new_jump(1, loops[1].stop()-1, 1+(loops[1].stop()-loops[1].start()));


		target_loops.fold_tree( loop_f );


		Size n1segment = loops[1].start()-1;
		Size c1segment = cut_points[0]-loops[1].stop();


		Size n2segment = loops[2].start()-cut_points[0];
		Size c2segment = cut_points[1] - loops[2].stop();


		Size n3segment = loops[3].start()-cut_points[1];
		Size c3segment = nat_prot.total_residue() - loops[3].stop();


		TR << "NSEGMENT1 "<< n1segment <<std::endl;
		TR << "CSEGMENT1 "<< c1segment <<std::endl;

		TR << "NSEGMENT2 "<< n2segment <<std::endl;
		TR << "CSEGMENT2 "<< c2segment <<std::endl;

		TR << "NSEGMENT3 "<< n3segment <<std::endl;

		TR << "CSEGMENT3 "<< c3segment <<std::endl;

		TR << "CUT1 "<< cut_points[0] <<std::endl;
		TR << "CUT2 "<< cut_points[1] <<std::endl;

		TR << "NAT SEQ  " << nat_seq <<std::endl;


		for ( Size k=1; k <= n1segment ; ++k ) {


			const char aa = nat_seq[ loops[1].start() - k - 1]; // The minus 1 is here to compensate the indexing starting in zero
			Size residue = loops[1].start() - k;


			TR << "RES AA N1  " << residue << aa <<std::endl;


			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, 1, true);


			target_loops.set_omega( 1, 180.0 );

		}

		//target_loops.dump_pdb("n1.pdb");

		for ( Size j = 0 ; j < c1segment ; ++j  ) {

			const char aa = nat_seq[ loops[1].stop() + j  ];
			Size residue = loops[1].stop() + j ;


			TR << "RES AA C1  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos( *new_rsd, residue , true );


			target_loops.set_omega( residue , 180.0 );


		}

		//target_loops.dump_pdb("c1.pdb");

		Size current_tot_res = target_loops.total_residue();

		TR << " Target Loops total residues before the second loop inserts  " << target_loops.total_residue() << std::endl;


		for ( Size l = 1; l < n2segment ; ++l  ) {


			const char aa = nat_seq[ loops[2].start() - l - 1  ];// need to check this one
			Size residue =   current_tot_res - (loops[2].size()+loops[3].size()-1) ; // we need to check which one we start


			TR << "RES AA N2  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}

		//target_loops.dump_pdb("n2.pdb");

		for ( Size m = 0; m < c2segment ; ++m  ) {

			const char aa = nat_seq[ loops[2].stop() + m  ];
			Size residue = loops[2].stop() + m ;


			TR << "RES AA C2  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}

		//target_loops.dump_pdb("c2.pdb");
		current_tot_res = target_loops.total_residue();

		TR << " Target Loops total residues before the third loop inserts  " << target_loops.total_residue() << std::endl;


		for ( Size n = 1; n < n3segment ; ++n  ) {


			const char aa = nat_seq[ loops[3].start() - n - 1  ];// need to check this one
			Size residue =   current_tot_res - (loops[3].stop()-loops[3].start()) ; // we need to check which one we start


			TR << "RES AA N3  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_prepend_polymer_residue_before_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}

		//target_loops.dump_pdb("n3.pdb");

		for ( Size o = 0; o < c3segment ; ++o  ) {

			const char aa = nat_seq[ loops[3].stop() + o  ];
			Size residue = loops[3].stop() + o ;


			TR << "RES AA C3  " << residue << aa <<std::endl;

			core::chemical::ResidueTypeCOP new_rsd_type( rsd_set->get_representative_type_name1( aa ) );


			core::conformation::ResidueOP new_rsd( core::conformation::ResidueFactory::create_residue( *new_rsd_type ) );


			target_loops.conformation().safely_append_polymer_residue_after_seqpos(*new_rsd, residue, true);


			target_loops.set_omega( residue , 180.0 );


		}

		//target_loops.dump_pdb("c3.pdb");

	}


	TR << "Last Residue    " << target_loops.total_residue() <<std::endl;

	core::pose::add_variant_type_to_pose_residue(target_loops, chemical::LOWER_TERMINUS_VARIANT, 1 );
	core::pose::add_variant_type_to_pose_residue(target_loops, chemical::UPPER_TERMINUS_VARIANT, target_loops.total_residue());

}


void refresh_cutpoints( core::pose::Pose & pose,
	std::vector<Size> & cut_points){

	for ( std::vector<Size>::iterator it = cut_points.begin(); it != cut_points.end(); ++it ) {

		TR << "Refresh cut_point "<< *it <<std::endl;
		Size residue = *it;
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_LOWER,residue ); // residue on the pose has to be assigned as a cut
		core::pose::add_variant_type_to_pose_residue(pose, chemical::CUTPOINT_UPPER, residue + 1);

	}


}


void copying_side_chains_swap_loop (
	core::pose::Pose & swap_loops,
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops,
	core::kinematics::MoveMapOP & movemap

){

	using namespace core::conformation;
	Size offsetres = 0;

	// Cleaning up the pose - > this step has to be performed before
	core::pose::remove_lower_terminus_type_from_pose_residue(swap_loops, 1 );
	core::pose::remove_upper_terminus_type_from_pose_residue(swap_loops, swap_loops.total_residue());


	movemap->set_chi(true);


	for ( Size pos = 1; pos <= fold_pose.total_residue(); ++pos ) {


		if ( loops.is_loop_residue( pos ) ) {

			++offsetres ;


			fold_pose.replace_residue( pos, swap_loops.residue( offsetres ), true);

			TR << "After Swap loop Residue  " << offsetres << std::endl;

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

	for ( Size pos = 1; pos <= pose.total_residue(); pos++ ) {
		if ( !loops.is_loop_residue( pos ) ) {
			residues_to_mutate[pos] = true;
			task->nonconst_residue_task(pos).restrict_absent_canonical_aas( allowed_aas );
		}

	}

}


void design_excluding_swap_loops (
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops_in,
	core::scoring::ScoreFunctionOP & scorefxn_fa

){


	using namespace basic::options;
	using namespace basic::options::OptionKeys;


	// PDBInfo setup to set the chains
	core::pose::PDBInfoOP pdb_info( new core::pose::PDBInfo( fold_pose ) );

	pdb_info->set_chains('A');

	fold_pose.pdb_info( pdb_info );


	core::pack::task::PackerTaskOP task( core::pack::task::TaskFactory::create_packer_task( fold_pose ));


	if ( option[ OptionKeys::packing::resfile ].user() ) {
		core::pack::task::parse_resfile(fold_pose, *task);
	} else {

		utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, true );

		utility::vector1< bool > residues_to_mutate( fold_pose.total_residue(), false );


		exclude_loop_residues(fold_pose, residues_to_mutate, allowed_aas, task, loops_in );


		if ( option[ OptionKeys::fold_from_loops::res_design_bs ].user() ) {

			for ( Size ex = 1; ex <= option[ OptionKeys::fold_from_loops::res_design_bs ]().size(); ex ++ ) {

				Size res_number = option[ OptionKeys::fold_from_loops::res_design_bs ]()[ex];
				residues_to_mutate[res_number] = true;

			}
		}


		task->restrict_to_residues( residues_to_mutate );

	}


	task->initialize_extra_rotamer_flags_from_command_line();


	pack::pack_rotamers( fold_pose, *scorefxn_fa , task );


}


void copying_side_chains(
	core::pose::Pose & nat_pose,
	core::pose::Pose & fold_pose,
	protocols::loops::Loops & loops,
	core::kinematics::MoveMapOP & movemap

){

	for ( Size pos = 1; pos <= fold_pose.total_residue(); pos++ ) {

		if ( loops.is_loop_residue( pos ) ) {
			fold_pose.replace_residue(pos,nat_pose.residue( pos ), true);
			movemap->set_chi(pos, false);
		}
	}

}


}
}
