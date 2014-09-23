// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c)University of Washington UW TechTransfer, email:license@u.washington.edu.

/// @file     protocols/antibody_legacy/CDRH3Modeler.cc
/// @brief    models CDR H3 loop using loop modeling
/// @detailed
/// @author   Aroop Sircar (aroopsircar@yahoo.com)

// Unit headers
#include <protocols/antibody_legacy/CDRH3Modeler.hh>

// Rosetta Headers
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/chemical/VariantType.hh>
#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/FragData.hh>
#include <core/fragment/FragSet.hh>
#include <core/id/types.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>

#include <numeric/numeric.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer TR( "protocols.antibody.CDRH3Modeler" );

namespace protocols {
namespace antibody_legacy {

using namespace core;
using namespace protocols::moves;

CDRH3Modeler::CDRH3Modeler(
    utility::vector1< fragment::FragSetOP > cdr_h3_frags
) : moves::Mover( "CDRH3Modeler" ) {
	cdr_h3_frags_ = cdr_h3_frags;
	set_default();
} // CDRH3Modeler default constructor

// CDRH3Modeler default destructor
CDRH3Modeler::~CDRH3Modeler() {}

void CDRH3Modeler::set_default() {
	benchmark_ = false;
	do_h3_modeling_ = false;
	base_ = 1;
	c_ter_stem_ = 3;
	cen_cst_ = 10.0;
	high_cst_ = 100.0; // if changed here, please change at the end of
	// AntibodyModeler as well

	lowres_scorefxn_ = scoring::ScoreFunctionFactory::
	                   create_score_function( "cen_std", "score4L" );
	lowres_scorefxn_->set_weight( scoring::chainbreak, 10./3. );
	// adding constraints
	lowres_scorefxn_->set_weight( scoring::atom_pair_constraint, cen_cst_ );

	highres_scorefxn_ = scoring::get_score_function();
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
	// adding constraints
	highres_scorefxn_->set_weight( scoring::atom_pair_constraint, high_cst_ );

	apply_centroid_mode_ = false;
	apply_fullatom_mode_ = false;

	current_loop_is_H3_ = true;
	H3_filter_ = true;
	antibody_build_ = true;
	antibody_refine_ = true;
	min_base_relax_ = false;
	h3_random_cut_ = false;
	decoy_loop_cutpoint_ = 0;
	snug_fit_ = true;
	loops_flag_ = true;
	docking_local_refine_ = true;
	dle_flag_ = true;
	refine_input_loop_ = true;
	h3_flank_ = 2;
	flank_relax_ = true;
	freeze_h3_ = true;
	is_camelid_ = false;
	base_ = 1;
	// size of loop above which 9mer frags are used
	cutoff_9_ = 16; // default 16
	// size of loop above which 3mer frags are used
	cutoff_3_ = 6; // default 6

	TR << "H3M Finished Setting Defaults" << std::endl;

	return;
} // CDRH3Modeler set_default

void CDRH3Modeler::set_lowres_score_func(
    scoring::ScoreFunctionOP lowres_scorefxn
) {
	lowres_scorefxn_ = lowres_scorefxn;
} // set_lowres_score_func

void CDRH3Modeler::set_highres_score_func(
    scoring::ScoreFunctionOP highres_scorefxn
) {
	highres_scorefxn_ = highres_scorefxn;
} // set_highres_score_func

void CDRH3Modeler::apply( pose::Pose & pose_in ) {
	if( !do_h3_modeling_ )
		return;

	TR << "H3M Applying CDR H3 modeler" << std::endl;

	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::moves;

	start_pose_ = pose_in;
	antibody_in_.set_Fv( pose_in, is_camelid_ );
	setup_packer_task( pose_in );
	pose::Pose start_pose = pose_in;

	if( is_camelid_ && !antibody_in_.extended_ && !antibody_in_.kinked_ )
		c_ter_stem_ = 0;

	Size framework_loop_begin( antibody_in_.cdrh_[3][1] );
	Size frmrk_loop_end_plus_one( antibody_in_.cdrh_[3][2] + 1 );
	//Size framework_loop_size = ( frmrk_loop_end_plus_one -
	//														framework_loop_begin ) + 1;
	Size cutpoint = framework_loop_begin + 1;
	loops::Loop cdr_h3( framework_loop_begin, frmrk_loop_end_plus_one,
	                    cutpoint,	0, true );
	simple_one_loop_fold_tree( antibody_in_.Fv, cdr_h3 );

	// switching to centroid mode
	simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

	// Building centroid mode loop
	if( apply_centroid_mode_ ) {
		to_centroid.apply( antibody_in_.Fv );
		build_centroid_loop();
		if( is_camelid_ )
			loop_centroid_relax( antibody_in_.Fv, antibody_in_.cdrh_[1][1],
			                     antibody_in_.cdrh_[1][2] );
		to_full_atom.apply( antibody_in_.Fv );

		utility::vector1<bool> allow_chi_copy( antibody_in_.Fv.total_residue(),
		                                       true );
		for( Size ii = antibody_in_.cdrh_[3][1];
		        ii <= ( antibody_in_.cdrh_[3][2] + 1 ); ii++ )
			allow_chi_copy[ii] = false;
		//recover sidechains from starting structures
		protocols::simple_moves::ReturnSidechainMover recover_sidechains(
		    start_pose_, allow_chi_copy );
		recover_sidechains.apply( antibody_in_.Fv );

		// Packer
		protocols::simple_moves::PackRotamersMoverOP packer;
		packer = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
		packer->task_factory(tf_);
		packer->apply( antibody_in_.Fv );
	}

	if( apply_fullatom_mode_ ) {
		build_fullatom_loop();
		if( !benchmark_ ) {
			Size repack_cycles(1);
			if( antibody_refine_ && !snug_fit_ )
				repack_cycles = 3;
			protocols::simple_moves::PackRotamersMoverOP packer;
			packer = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
			packer->task_factory(tf_);
			packer->nloop( repack_cycles );
			packer->apply( antibody_in_.Fv );
		}
	}

	// Minimize CDR H2 loop if this is a camelid
	if( is_camelid_ ) {
		bool store_current_loop = current_loop_is_H3_;
		bool store_H3_filter = H3_filter_;
		current_loop_is_H3_ = false;
		H3_filter_ = false;
		bool closed_cutpoints( false );

		Antibody starting_antibody;
		starting_antibody = antibody_in_;

		while( !closed_cutpoints) {
			antibody_in_ = starting_antibody;
			loop_fa_relax( antibody_in_.Fv, antibody_in_.cdrh_[1][1],
			               antibody_in_.cdrh_[1][2]  );
			closed_cutpoints = cutpoints_separation();
		} // while( ( cut_separation > 1.9 )

		// Restoring variables to initial state
		current_loop_is_H3_ = store_current_loop;
		H3_filter_ = store_H3_filter;
	}

	pose_in = antibody_in_.Fv;

	TR << "H3M Finished applying CDR H3 modeler" << std::endl;

	return;
} // CDRH3Modeler::apply()

std::string
CDRH3Modeler::get_name() const {
	return "CDRH3Modeler";
}


void CDRH3Modeler::build_centroid_loop() {
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::moves;

	if( !apply_centroid_mode_ )
		return;

	TR <<  "H3M Modeling Centroid CDR H3 loop" << std::endl;

	Size frmrk_loop_end_plus_one( antibody_in_.cdrh_[3][2] + 1 );
	Size framework_loop_size = ( frmrk_loop_end_plus_one -
	                             antibody_in_.cdrh_[3][1] ) + 1;
	Size cutpoint = antibody_in_.cdrh_[3][1] + 1;
	loops::Loop cdr_h3( antibody_in_.cdrh_[3][1], frmrk_loop_end_plus_one,
	                    cutpoint,	0, true );
	simple_one_loop_fold_tree( antibody_in_.Fv, cdr_h3 );

	// silly hack to make extended loops to work
	loops::LoopsOP cdr_h3_loop_list = new loops::Loops();
	cdr_h3_loop_list->add_loop( cdr_h3 );
	/* Commented out by BDW with JX's consent
	loops::loop_mover::LoopMoverOP my_loop_move =  new loops::loop_mover::LoopMover( cdr_h3_loop_list );
	my_loop_move->set_extended_torsions( antibody_in_.Fv, cdr_h3 );
	my_loop_move->apply( antibody_in_.Fv );
	*/

	Size unaligned_cdr_loop_begin(0)/*, unaligned_cdr_loop_end(0)*/;
	core::import_pose::pose_from_pdb( template_pose_, "hfr.pdb" );
	std::string template_name = "h3";
	Antibody hfr_template( template_pose_, template_name );
	unaligned_cdr_loop_begin = hfr_template.current_start;
	//unaligned_cdr_loop_end = hfr_template.current_end;  // set but never used ~Labonte

	antibody_in_.Fv.set_psi( antibody_in_.cdrh_[3][1] - 1,
	                         template_pose_.psi( unaligned_cdr_loop_begin - 1 ) );
	antibody_in_.Fv.set_omega(antibody_in_.cdrh_[3][1] - 1,
	                          template_pose_.omega( unaligned_cdr_loop_begin - 1 ) );

	Size modified_framework_loop_end = frmrk_loop_end_plus_one - c_ter_stem_;
	loops::Loop trimmed_cdr_h3( antibody_in_.cdrh_[3][1],
	                            modified_framework_loop_end, cutpoint, 0, true );

	Antibody starting_antibody;
	starting_antibody = antibody_in_;
	bool closed_cutpoints( false );

	while( !closed_cutpoints) {
		antibody_in_ = starting_antibody;
		if( framework_loop_size > 5 )
			antibody_modeling_insert_ter();
		scored_frag_close( antibody_in_.Fv, trimmed_cdr_h3 );
		if( trimmed_cdr_h3.size() > cutoff_9_  ) { // aroop_temp default cutoff_9_
			Size saved_cutoff_9 = cutoff_9_;
			cutoff_9_ = 100; // never going to reach
			scored_frag_close( antibody_in_.Fv, trimmed_cdr_h3 );
			cutoff_9_ = saved_cutoff_9; // restoring
		}
		closed_cutpoints = cutpoints_separation();
	} // while( ( cut_separation > 1.9 )

	TR <<  "H3M Finished Modeling Centroid CDR H3 loop" << std::endl;

	return;

} // build_centroid_loop

void CDRH3Modeler::build_fullatom_loop() {
	using namespace core::pose;
	using namespace core::scoring;
	using namespace protocols::moves;

	if( !apply_fullatom_mode_ )
		return;

	TR <<  "H3M Modeling Fullatom CDR H3 loop" << std::endl;

	Antibody starting_antibody;
	starting_antibody = antibody_in_;
	bool closed_cutpoints( false );

	while( !closed_cutpoints) {
		antibody_in_ = starting_antibody;
		loop_fa_relax( antibody_in_.Fv, antibody_in_.cdrh_[3][1],
		               antibody_in_.cdrh_[3][2] + base_ );
		closed_cutpoints = cutpoints_separation();
	} // while( ( cut_separation > 1.9 )

	TR <<  "H3M Finished modeling Fullatom CDR H3 loop" << std::endl;

	return;
} // build_fullatom_loop

void CDRH3Modeler::store_H3_cter_fragment(
    utility::vector1< fragment::FragData > & base_library_in
) {
	H3_base_library = base_library_in;
	return;
} // store_H3_cter_fragment

void simple_one_loop_fold_tree(
    pose::Pose & pose_in,
    loops::Loop const & loop
) {
	using namespace kinematics;

	TR <<  "H3M Setting up simple one loop fold tree" << std::endl;

	//setup fold tree for this loop
	FoldTree f;
	f.clear();
	Size nres = pose_in.total_residue();
	Size jumppoint1 = loop.start() - 1;
	Size jumppoint2 = loop.stop() + 1;

	if( jumppoint1 < 1 )   jumppoint1 = 1;
	if( jumppoint2 > nres) jumppoint2 = nres;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, loop.cut(), Edge::PEPTIDE );
	f.add_edge( loop.cut() + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR <<  "H3M Finished setting up simple one loop fold tree" << std::endl;

	return;
} // simple_one_loop_fold_tree

void simple_fold_tree(
    pose::Pose & pose_in,
    Size jumppoint1,
    Size cutpoint,
    Size jumppoint2
) {
	using namespace kinematics;

	TR <<  "H3M Setting up simple fold tree" << std::endl;

	//setup fold tree for this loop
	FoldTree f;
	f.clear();
	Size nres = pose_in.total_residue();

	if( jumppoint1 < 1 )   jumppoint1 = 1;
	if( jumppoint2 > nres) jumppoint2 = nres;

	f.add_edge( 1, jumppoint1, Edge::PEPTIDE );
	f.add_edge( jumppoint1, cutpoint, Edge::PEPTIDE );
	f.add_edge( cutpoint + 1, jumppoint2, Edge::PEPTIDE );
	f.add_edge( jumppoint2, nres, Edge::PEPTIDE );
	f.add_edge( jumppoint1, jumppoint2, 1 );
	f.reorder( 1 );

	pose_in.fold_tree( f );

	TR <<  "H3M Finished setting up simple fold tree" << std::endl;

	return;
} // simple_fold_tree

void read_H3_cter_fragment(
    Antibody & antibody_in,
    utility::vector1< fragment::FragData > & H3_base_library,
    bool is_camelid
) {
	using namespace fragment;

	TR <<  "H3M Reading CDR H3 C-ter Fragments" << std::endl;

	bool is_kinked( antibody_in.kinked_ );
	bool is_extended( antibody_in.extended_ );

	// extract single letter aa codes for the chopped loop residues
	Size cdr_h3_size = ( antibody_in.cdrh_[3][2] -
	                     antibody_in.cdrh_[3][1] ) + 1;
	utility::vector1< char > aa_1name;
	for( Size ii = antibody_in.cdrh_[3][1] - 2;
	        ii <= ( antibody_in.cdrh_[3][1] - 2 ) + cdr_h3_size + 3; ++ii )
		aa_1name.push_back( antibody_in.Fv_sequence_[ii] );

	// used only when no length & kink match are found
	utility::vector1< FragData > H3_base_library_seq_kink;

	// used only when no (length & kink) or (length & seq) are found
	utility::vector1< FragData > H3_base_library_kink;

	std::string H3_ter_library_filename;
	// file is read in from where other contraints are supposed to exist
	if( is_camelid )
		H3_ter_library_filename = "camelid_H3_CTERM";
	else
		H3_ter_library_filename = "H3_CTERM";

	// Read the file defined by command line option
	utility::io::izstream H3_ter_library_stream( H3_ter_library_filename );

	// Check to see if file exists
	if ( !H3_ter_library_stream ) {
		TR << "[Error]: Could not open H3 base library file: "
		   << H3_ter_library_filename << std::endl
		   << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
		std::exit( EXIT_FAILURE );
	}

	std::string pdb_name;
	std::string res_no;
	char res_name;
	Real phi(0.0);
	Real psi(0.0);
	Real omega(0.0);
	Size H3_length(0);
	Real resolution(0.0);
	std::string base_type;

	Size pdb_H3_length = cdr_h3_size;
	Size h3_base_frag_size( is_camelid ? 6 : 4 ); // changed from 4:6
	bool end_not_reached(true);
	while(end_not_reached) {
		bool seq_match( true );
		bool kink_match( false );

		FragData f;
		f.set_valid( true );

		for ( Size i = 1; i <= h3_base_frag_size; ++i ) {
			H3_ter_library_stream >> pdb_name
			                      >> res_no
			                      >> res_name
			                      >> omega
			                      >> phi
			                      >> psi
			                      >> H3_length
			                      >> resolution
			                      >> base_type
			                      >> std::skipws;
			if ( H3_ter_library_stream.eof() ) {
				end_not_reached = false;
				break;
			}
			if( res_name != aa_1name[aa_1name.size() - 5 + i] )
				seq_match = false;

			BBTorsionSRFDOP res_torsions( new BBTorsionSRFD( 3, 'L', res_name ) ); // 3 protein torsions
			// ugly numbers 1-3, but pose.set_phi also uses explicit numbers
			res_torsions->set_torsion   ( 1, phi   );
			res_torsions->set_torsion   ( 2, psi   );
			res_torsions->set_torsion   ( 3, omega );
			res_torsions->set_secstruct ( 'L' );
			f.add_residue( res_torsions );
		}
		if( !is_camelid ) {
			if( is_kinked && base_type == "KINK" )
				kink_match = true;
			else if( is_extended && base_type == "EXTENDED" )
				kink_match = true;
			else if( !is_kinked && !is_extended && base_type == "NEUTRAL" )
				kink_match = true;
		} else {
			if( is_extended && base_type == "EXTENDED" )
				kink_match = true;
			else if( ( is_kinked && base_type == "KINK" ) ||
			         ( is_kinked && base_type == "EXTENDED" ) )
				kink_match = true;
			else if( !is_kinked && !is_extended )
				kink_match = true;
		}
		if( is_camelid && end_not_reached && kink_match ) {
			H3_base_library.push_back( f );
		} else if( end_not_reached && ( H3_length == pdb_H3_length )
		           && kink_match ) {
			H3_base_library.push_back( f );
		}
		if( end_not_reached && seq_match && kink_match ) {
			H3_base_library_seq_kink.push_back( f );
		}
		if( end_not_reached && kink_match  ) {
			H3_base_library_kink.push_back( f );
		}
	}

	H3_ter_library_stream.close();
	H3_ter_library_stream.clear();

	// if no match found based on sequence and kink match criterion
	// then choose based on size and kink match criterion
	// if still no match, then choose based only on kink
	if( H3_base_library.size() == 0 ) {
		H3_base_library = H3_base_library_seq_kink;
	}
	if( H3_base_library.size() == 0 ) {
		H3_base_library = H3_base_library_kink;
	}

	TR <<  "H3M Finished reading CDR H3 C-ter Fragments" << std::endl;

	return;
}

void CDRH3Modeler::antibody_modeling_insert_ter() {

	TR <<  "H3M Inserting CDR H3 C-ter Fragments" << std::endl;

	// Storing initial fold tree
	kinematics::FoldTree const input_tree( antibody_in_.Fv.fold_tree() );

	Size loop_begin(0), loop_end(0), cutpoint(0), random_H3_ter(0);
	//utility::vector1< fragment::FragData >::const_iterator H3_ter;
	fragment::FragData f;

	loop_begin = antibody_in_.cdrh_[3][1];
	cutpoint = antibody_in_.cdrh_[3][1] + 1;
	random_H3_ter = numeric::random::rg().random_range( 1, H3_base_library.size() );
	//H3_ter = H3_base_library.begin();

	loop_end = antibody_in_.cdrh_[3][2] + 1;

	loops::Loop cdr_h3( loop_begin, loop_end, cutpoint,	0, true );
	simple_one_loop_fold_tree( antibody_in_.Fv, cdr_h3 );

	// choosing a base randomly
	//H3_ter = H3_ter + random_H3_ter;
	f = H3_base_library[ random_H3_ter ];

	pose::Pose start_pose = antibody_in_.Fv;
	//inserting base dihedrals
	Size cter_insertion_pos( is_camelid_ ? 4 : 2 );
	if( (antibody_in_.cdrh_[3][2] - cter_insertion_pos) <=
	        antibody_in_.cdrh_[3][1] ) {
		TR << "H3 LOOP IS TOO SHORT: CAN NOT USE N-TERM INFORMATION" << std::endl;
	} else {
		// H3_ter->apply(...);
		f.apply( antibody_in_.Fv, antibody_in_.cdrh_[3][2] -
		         cter_insertion_pos, antibody_in_.cdrh_[3][2] + 1 );
	}

	// Restoring pose fold tree
	antibody_in_.Fv.fold_tree( input_tree );

	TR <<  "H3M Finished Inserting CDR H3 C-ter Fragments" << std::endl;

	return;
} // antibody_modeling_insert_ter

bool CDRH3Modeler::cutpoints_separation() {

	bool closed_cutpoints = true;

	for( loops::Loops::const_iterator it=antibody_in_.all_cdr_loops.begin(),
	        it_end=antibody_in_.all_cdr_loops.end(),
	        it_next; it != it_end; ++it ) {
		Size cutpoint   = it->cut();
		Real separation = 10.00; // an unlikely high number
		separation = cutpoint_separation( antibody_in_.Fv, cutpoint );

		if( separation > 1.9 ) {
			closed_cutpoints = false;
			break;
		}
	}
	return( closed_cutpoints );
} // cutpoints_separation

Real CDRH3Modeler::cutpoint_separation(
    pose::Pose & pose_in,
    Size cutpoint ) {

	Size const N ( 1 ); // N atom
	Size const C ( 3 ); // C atom

	// Coordinates of the C atom of cutpoint res and N atom of res cutpoint+1
	numeric::xyzVector_float peptide_C(pose_in.residue( cutpoint ).xyz( C )),
	        peptide_N( pose_in.residue( cutpoint + 1 ).xyz( N ) );
	Real cutpoint_separation=peptide_C.distance(peptide_N);

	return( cutpoint_separation );
} // cutpoint_separation

///////////////////////////////////////////////////////////////////////////
/// @begin scored_frag_close
///
/// @brief builds a loop from fragments file.
///
/// @detailed Loop is built by a monte carlo simulation using fragments
///           from a fragment files. CCD moves are used to close loops
///           with gaps at cutpoint.H3_check is enforced if H3_filter flag
///           is set in command line. Loop building results in many files
///           containing differnt conformations of the same loop in
///           phi-psi-omega angles. Parallel processing is utilized.
///
/// @param[in] weight_map: in this case its a centroid weight
///            pose_in: loop to be built on this template provided
///            loop_begin/loop_end: loop termini definition
///            frag_size: 3-mer or 9-mer
///            frag_offset:agreement in frag file numbering & pose numberng
///            cycles1: max cycles to be spent building loops
///            cycles2: # of fragment swaps for each loop(depends on size)
///            do_ccd_moves: should ccd moves be used to close gaps
///
/// @global_read benchmark_
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @authors Aroop 02/04/2010
///
/// @last_modified 02/04/2010
///////////////////////////////////////////////////////////////////////////
void CDRH3Modeler::scored_frag_close (
    pose::Pose & pose_in,
    loops::Loop const trimmed_cdr_h3 ) {
	using namespace fragment;
	using namespace protocols;
	using namespace protocols::simple_moves;
	using namespace protocols::loops;
	using loop_closure::ccd::CcdMover;
	using loop_closure::ccd::CcdMoverOP;
	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	TR <<  "H3M Fragments based centroid CDR H3 loop building" << std::endl;

	if( trimmed_cdr_h3.size() <= 2)
		utility_exit_with_message("Loop too small for modeling");

	// set cutpoint variants for correct chainbreak scoring
	if( !pose_in.residue( trimmed_cdr_h3.cut() ).is_upper_terminus() ) {
		if( !pose_in.residue( trimmed_cdr_h3.cut() ).has_variant_type(
		            chemical::CUTPOINT_LOWER))
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_LOWER,
			        trimmed_cdr_h3.cut() );
		if( !pose_in.residue( trimmed_cdr_h3.cut() + 1 ).has_variant_type(
		            chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_UPPER,
			        trimmed_cdr_h3.cut() + 1 );
	}


	Size cycles1(10);
	// aroop_temp default 25 * loop size
	Size cycles2(25 * trimmed_cdr_h3.size() );

	// params
	Real const ccd_threshold( 0.1);
	Size h3_attempts(0);
	Real h3_fraction = 0.75; // 75% of loops are required to be H3's
	Real current_h3_prob = numeric::random::rg().uniform();;
	bool H3_found_ever(false);
	bool loop_found(false);
	Size total_cycles(0);
	Size frag_size(0);
	FragSetOP frags_to_use;
	{
		if( trimmed_cdr_h3.size() > cutoff_9_ ) {
			frags_to_use = cdr_h3_frags_[1]->empty_clone();
			frags_to_use = cdr_h3_frags_[1];
			frag_size = 9;
		} else {
			frags_to_use = cdr_h3_frags_[2]->empty_clone();
			frags_to_use = cdr_h3_frags_[2];
			frag_size = 3;
		}
	}

	// Storing Fold Tree
	kinematics::FoldTree old_fold_tree = pose_in.fold_tree();
	// New Fold Tree
	simple_one_loop_fold_tree( pose_in, trimmed_cdr_h3 );

	//setting MoveMap
	kinematics::MoveMapOP cdrh3_map;
	cdrh3_map = new kinematics::MoveMap();
	cdrh3_map->clear();
	cdrh3_map->set_chi( true );
	cdrh3_map->set_bb( false );
	for( Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop(); ii++ ) {
		cdrh3_map->set_bb( ii, true );
	}
	cdrh3_map->set_jump( 1, false );

	// setup monte_carlo
	Real temp( 2.0);
	MonteCarloOP mc, outer_mc;
	mc = new moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
	outer_mc = new moves::MonteCarlo( pose_in, *lowres_scorefxn_, temp );
	Size buffer( (is_camelid_ && antibody_in_.extended_) ? 2 : 0 );
	while( !loop_found && ( total_cycles++ < cycles1) ) {
		// insert random fragments over the whole loop
		for(Size ii = trimmed_cdr_h3.start(); ii<=trimmed_cdr_h3.stop()
		        - ( buffer + (frag_size - 1 ) ); ii++ ) {
			ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use,
			        cdrh3_map);
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			cfm->define_start_window( ii );
			cfm->apply( pose_in );
		}
		if( total_cycles == 1 ) {
			mc->reset( pose_in );
		}

		Size local_h3_attempts(0);
		for ( Size c2 = 1; c2 <= cycles2; ++c2 ) {
			// apply a random fragment
			ClassicFragmentMoverOP cfm = new ClassicFragmentMover( frags_to_use,
			        cdrh3_map);
			cfm->set_check_ss( false );
			cfm->enable_end_bias_check( false );
			cfm->apply( pose_in );

			bool H3_found_current(false);
			if( current_loop_is_H3_ && H3_filter_ &&
			        ( local_h3_attempts++ < (50 * cycles2) ) ) {
				H3_found_current = CDR_H3_filter(pose_in,antibody_in_.cdrh_[3][1],
				                                 ( antibody_in_.cdrh_[3][2] - antibody_in_.cdrh_[3][1] ) + 1 );
				if( !H3_found_ever && !H3_found_current) {
					--c2;
					mc->boltzmann( pose_in );
				} else if( !H3_found_ever && H3_found_current ) {
					H3_found_ever = true;
					mc->reset( pose_in );
				} else if( H3_found_ever && !H3_found_current ) {
					--c2;
					continue;
				} else if( H3_found_ever && H3_found_current )
					mc->boltzmann( pose_in );
			} else {
				mc->boltzmann( pose_in );
			}

			if ( (c2 > cycles2/2 && numeric::random::rg().uniform() * cycles2 < c2) ||
			        ( trimmed_cdr_h3.size() <= 5) ) {
				// in 2nd half of simulation, start trying to close the loop:
				CCDLoopClosureMoverOP ccd_moves = new CCDLoopClosureMover( trimmed_cdr_h3, cdrh3_map );
				RepeatMoverOP ccd_cycle;
				if( trimmed_cdr_h3.size() <= 5 ) {
					ccd_cycle = new RepeatMover(ccd_moves,500*trimmed_cdr_h3.size());
					ccd_cycle->apply( pose_in );
				} else {
					ccd_cycle = new RepeatMover(ccd_moves, 10*trimmed_cdr_h3.size());
					ccd_cycle->apply( pose_in );
				}
				mc->boltzmann( pose_in );
			}
		}

		mc->recover_low( pose_in );
		CCDLoopClosureMoverOP ccd_closure = new CCDLoopClosureMover(
		    trimmed_cdr_h3, cdrh3_map );
		ccd_closure->tolerance( ccd_threshold );
		ccd_closure->max_cycles( 500 );
		ccd_closure->apply( pose_in );

		if( total_cycles == 1 )
			outer_mc->reset( pose_in );

		if ( ccd_closure->deviation() <= ccd_threshold ) {
			// CDR-H3 filter for antibody mode
			// introduce enough diversity
			outer_mc->boltzmann( pose_in );
			if( current_loop_is_H3_ && H3_filter_ &&
			        (current_h3_prob < h3_fraction) && (h3_attempts++<50) )
				if( !CDR_H3_filter(pose_in, antibody_in_.cdrh_[3][1],
				                   ( antibody_in_.cdrh_[3][2] - antibody_in_.cdrh_[3][1] ) + 1) )
					continue;
			loop_found = true;
		} else if( H3_filter_ ) {
			h3_attempts++;
		}
	}
	outer_mc->recover_low( pose_in );

	// Restoring Fold Tree
	pose_in.fold_tree( old_fold_tree );

	TR <<  "H3M Finished Fragments based centroid CDR H3 loop building"
	   << std::endl;

	return;
} // scored_frag_close

///////////////////////////////////////////////////////////////////////////
/// @begin CDR_H3_filter
///
/// @brief tests if a loop has H3 like base charachteristics
///
/// @detailed Uses the Shirai rules to find out if the dihedral angle
///           formed by CA atoms of residues n-2,n-1,n and n+1 conform to a
///           kinked/extended structure in accordance with the sequence. If
///           there is a match, a true value is returned
///
/// @param[in] pose: full actual protein
///            loop_begin: seq numbered loop begin corresponding to pose
///            size: size of loop to compute loop_end
///
/// @global_read reads -command line flag -base stored in dle_ns
///              to determine to do the complete H3 filter check or just do
///              a prediction of the H3 base type based on the
///              aforementioned dihedral angle
///
/// @global_write
///
/// @remarks
///
/// @references Structural classification of CDR-H3 in antibodies
///             Hiroki Shirai, Akinori Kidera, Haruki Nakamura
///             FEBS Letters 399 (1996) 1-8
///
/// @authors Aroop 02/04/2010
///
/// @last_modified 02/04/2010
///////////////////////////////////////////////////////////////////////////
bool CDRH3Modeler::CDR_H3_filter( const pose::Pose & pose_in,
		Size const loop_begin,
		Size const size,
		char const light_chain )
{
	TR <<  "H3M Checking Kink/Extended CDR H3 Base Angle" << std::endl;

	if( !H3_filter_ || is_camelid_ )
		return( true );

	// Values read from plot in reference paper. Fig 1 on Page 3
	// Values adjusted to match data from antibody training set
	Real const kink_lower_bound = -10.00; // Shirai: 0
	Real const kink_upper_bound = 70.00; // Shirai: 70
	Real const extended_lower_bound = 125.00; // Shirai: ~180
	Real const extended_upper_bound = 185.00; // Shirai: ~180

	// Hydrogen Bond maximum value is 3.9 Angstroms - not used
	//	Real const h_bond(3.9);
	// Salt Bridge maximum value is 2.0 Angstroms - not used
	//	Real const s_bridge(4.0);

	// chop out the loop
	pose::Pose h3_loop( pose_in, loop_begin - 2, loop_begin + size + 1 );

	//bool is_kinked( false );
	//bool is_extended( false );
	bool is_H3( false );

	// extract 3 letter residue codes for the chopped loop
	std::vector <std::string> aa_name; // loop residue 3 letter codes
	for(Size ii = 1; ii <= size + 3; ii++)
		aa_name.push_back( h3_loop.residue(ii).name3() );

	Size const CA(2);   // CA atom position in full_coord array
	// base dihedral angle to determine kinked/extended conformation
	Real base_dihedral( numeric::dihedral_degrees(
	                        h3_loop.residue( aa_name.size() ).xyz( CA ),
	                        h3_loop.residue( aa_name.size() - 1).xyz( CA ),
	                        h3_loop.residue( aa_name.size() - 2).xyz( CA ),
	                        h3_loop.residue( aa_name.size() - 3).xyz( CA ) ) );

	// std::cout << "Base Dihedral: " << base_dihedral << std::endl;

	// setting up pseudo-periodic range used in extended base computation
	if( base_dihedral < kink_lower_bound )
		base_dihedral = base_dihedral + 360.00;


	// Rule 1a for standard kink
	if ((aa_name[aa_name.size()-3] != "ASP") &&
	        (aa_name[aa_name.size()-1] == "TRP"))	{
		if( (base_dihedral > kink_lower_bound) &&
		        (base_dihedral < kink_upper_bound)) {
			// std::cout << "KINK Found" << std::endl; // aroop_temp remove
			//is_kinked = true;  // set but never used ~Labonte
			is_H3 = true;
		}
	}

	// Rule 1b for standard extended form
	if ( ( aa_name[ aa_name.size() - 3 ] == "ASP" ) &&
	        ( ( aa_name[1] != "LYS" ) && ( aa_name[1] != "ARG" ) )&&
	        ( is_H3 != true ) ) {
		if( ( base_dihedral > extended_lower_bound ) &&
		        ( base_dihedral < extended_upper_bound) ) {
			// std::cout << "EXTENDED Found" << std::endl; // aroop_temp remove
			//is_extended = true;  // set but never used ~Labonte
			is_H3 = true;
		}

		if(!is_H3) {
			// Rule 1b extension for special kinked form
			bool is_basic(false); // Special basic residue exception flag
			for(Size ii = 2; ii <= Size(aa_name.size() - 5); ii++) {
				if( aa_name[ii] == "ARG" || aa_name[ii] == "LYS" ) {
					is_basic = true;
					break;
				}
			}

			if(!is_basic) {
				Size rosetta_number_of_L49 = pose_in.pdb_info()->pdb2pose(
				                                 light_chain, 49 );
				std::string let3_code_L49 =
				    pose_in.residue( rosetta_number_of_L49 ).name3();
				if( let3_code_L49 == "ARG" || let3_code_L49 == "LYS")
					is_basic = true;
			}
			if( is_basic && ( base_dihedral > kink_lower_bound ) &&
			        ( base_dihedral < kink_upper_bound ) ) {
				// aroop_temp remove
				// std::cout << "KINK (special 1b) Found" << std::endl;
				//is_kinked = true;  // set but never used ~Labonte
				is_H3 = true;
			}
		}
	}

	// Rule 1c for kinked form with salt bridge
	if ( ( aa_name[ aa_name.size() - 3 ] == "ASP") &&
	        ( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
	        ( (aa_name[0] != "LYS" ) && ( aa_name[0] != "ARG" ) ) &&
	        ( is_H3 != true) ) {
		if( (base_dihedral > kink_lower_bound ) &&
		        (base_dihedral < kink_upper_bound ) ) {
			// aroop_temp remove
			// std::cout << "KINK (w sb) Found" << std::endl;
			//is_kinked = true;  // set but never used ~Labonte
			is_H3 = true;
		}
		if(!is_H3) {
			bool is_basic(false); // Special basic residue exception flag
			Size rosetta_number_of_L46 = pose_in.pdb_info()->pdb2pose(
			                                 light_chain, 46 );
			std::string let3_code_L46 =
			    pose_in.residue( rosetta_number_of_L46 ).name3();
			if( let3_code_L46 == "ARG" || let3_code_L46 == "LYS")
				is_basic = true;
			if( is_basic && (base_dihedral > extended_lower_bound ) &&
			        ( base_dihedral < extended_upper_bound ) ) {
				// aroop_temp remove
				// std::cout << "EXTENDED (special 1c) Found" << std::endl;
				//is_extended = true;  // set but never used ~Labonte
				is_H3 = true;
			}
		}
	}

	// Rule 1d for extened form with salt bridge
	if ( ( aa_name[ aa_name.size() - 3 ] == "ASP") &&
	        ( ( aa_name[1] == "LYS") || ( aa_name[1] == "ARG" ) ) &&
	        ( ( aa_name[0] == "LYS") || ( aa_name[0] == "ARG") ) &&
	        ( is_H3 != true ) ) {
		if( ( base_dihedral > extended_lower_bound ) &&
		        ( base_dihedral < extended_upper_bound ) ) {
			// aroop_temp remove
			// std::cout << "EXTENDED (w sb) Found" << std::endl;
			//is_extended = true;  // set but never used ~Labonte
			is_H3 = true;
		}
	}

	TR <<  "H3M Finished Checking Kink/Extended CDR H3 Base Angle: "
	   << is_H3 << std::endl;

	return is_H3;
} // CDR_H3_filter

///////////////////////////////////////////////////////////////////////////
/// @begin loop_fa_relax
///
/// @brief actually relaxes the region specified
///
/// @detailed This is all done in high resolution.Hence there are no rigid
///           body moves relative to the docking partners. Only small moves
///           are carried out here to see if there are better fits.
///           Repacking is carried out extensively after each move.
///
/// @param[in] pose, loop begin position, loop end position
///
/// @global_read none
///
/// @global_write none
///
/// @remarks
///
/// @references
///
/// @authors Aroop 02/04/2010
///
/// @last_modified 02/04/2010
///////////////////////////////////////////////////////////////////////////
void CDRH3Modeler::loop_fa_relax( pose::Pose & pose_in, Size const loop_begin, Size const loop_end )
{
	using namespace protocols;
	using namespace protocols::simple_moves;
	using namespace protocols::loops;
	using namespace protocols::moves;
	using namespace protocols::toolbox::task_operations;
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	TR << "H3M Relaxing CDR H3 Loop" << std::endl;

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	//setting MoveMap
	kinematics::MoveMapOP cdrh3_map;
	cdrh3_map = new kinematics::MoveMap();
	cdrh3_map->clear();
	cdrh3_map->set_chi( false );
	cdrh3_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
	for( Size ii = loop_begin; ii <= loop_end; ii++ )
		allow_bb_move[ ii ] = true;
	cdrh3_map->set_bb( allow_bb_move );
	cdrh3_map->set_jump( 1, false );


	// minimize_set_local_min( false, 0 );// all non-move-list rsds minimized

	Size loop_size = ( loop_end - loop_begin ) + 1;
	Size cutpoint = loop_begin + int(loop_size/2);
	if( current_loop_is_H3_ ) {
		if( (antibody_build_ || antibody_refine_ ) &&
		        !min_base_relax_ && !h3_random_cut_ &&
		        (decoy_loop_cutpoint_ != 0))
			cutpoint = decoy_loop_cutpoint_;
		//else if( h3_random_cut_ )
		//	cutpoint = dle_choose_random_cutpoint(loop_begin, loop_end);
	} else
		cutpoint = loop_begin + Size( loop_size / 2);
	/*
	if( snug_fit_ && loops_flag_ && docking_local_refine_ &&
			dle_flag_ ) {
		one_loop = dle_ns::dle_loops;
	}
	else */
	loops::Loop one_loop( loop_begin, loop_end,	cutpoint,	0, false );

	// sets up stuff to use rosetta's fullatom energy function
	//initialize_fullatom();
	// maximum number of rotamers to allow (buried, surface)
	//set_rot_limit( 45, 27 );
	// maximum sum for the dunbrak prob (buried,surf)
	//set_perc_limit( 0.9999, 0.9999 );
	// include extra rotamers in chi1/chi2/chi1aro
	//design::active_rotamer_options.set_ex12( true, true, true);
	// checking if old rosetta full atom flag is on

	if ( !pose_in.is_fullatom() )
		utility_exit_with_message("Fullatom poses only");

	ChangeFoldTreeMoverOP one_loop_fold_tree;
	ChangeFoldTreeMoverOP with_flank_fold_tree;
	simple_fold_tree( pose_in, loop_begin - 1, cutpoint, loop_end + 1 );
	one_loop_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );
	with_flank_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );

	//////////////////
	// setup fold_tree
	utility::vector1< bool> flank_allow_bb_move( allow_bb_move  );
	if( current_loop_is_H3_  && flank_relax_ && freeze_h3_) {
		simple_fold_tree( pose_in, loop_begin - h3_flank_ - 1, cutpoint,
		                  loop_end + h3_flank_ + 1 );
		with_flank_fold_tree = new ChangeFoldTreeMover( pose_in.fold_tree() );
		for( Size i = 1; i <= pose_in.total_residue(); i++ )
			if( (i >= (loop_begin - h3_flank_)) && (i <= (loop_end + h3_flank_)))
				flank_allow_bb_move[i] = true;
	} else
		one_loop_fold_tree->apply( pose_in );

	// set cutpoint variants for correct chainbreak scoring
	if( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
		if( !pose_in.residue( cutpoint ).has_variant_type(
		            chemical::CUTPOINT_LOWER))
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_LOWER,
			        cutpoint );
		if( !pose_in.residue( cutpoint + 1 ).has_variant_type(
		            chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_UPPER,
			        cutpoint + 1 );
	}



	utility::vector1< bool> allow_repack( pose_in.total_residue(), false );
	select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
	                      allow_repack);
	cdrh3_map->set_chi( allow_repack );

	protocols::simple_moves::PackRotamersMoverOP loop_repack=new protocols::simple_moves::PackRotamersMover(highres_scorefxn_);
	setup_packer_task( start_pose_ );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new RestrictToInterface( allow_repack ) );
	loop_repack->task_factory(tf_);
	// loop_repack->apply( pose_in );

	Real min_tolerance = 0.001;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
	bool nb_list = true;
	protocols::simple_moves::MinMoverOP loop_min_mover = new protocols::simple_moves::MinMover( cdrh3_map,
	        highres_scorefxn_, min_type, min_tolerance, nb_list );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	Size inner_cycles( loop_size );
	Size outer_cycles( 1 );
	if( antibody_refine_ || refine_input_loop_ )
		outer_cycles = 5;
	if( antibody_refine_ && snug_fit_ )
		outer_cycles = 2;
	if( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles = 1;
		outer_cycles = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
	protocols::simple_moves::BackboneMoverOP small_mover = new protocols::simple_moves::SmallMover( cdrh3_map,
	        high_move_temp,
	        n_small_moves );
	protocols::simple_moves::BackboneMoverOP shear_mover = new protocols::simple_moves::ShearMover( cdrh3_map,
	        high_move_temp,
	        n_small_moves );
	if( min_base_relax_ ) {
		small_mover->angle_max( 'H', 0.5 );
		small_mover->angle_max( 'E', 0.5 );
		small_mover->angle_max( 'L', 1.0 );
		shear_mover->angle_max( 'H', 0.5 );
		shear_mover->angle_max( 'E', 0.5 );
		shear_mover->angle_max( 'L', 1.0 );
	} else {
		small_mover->angle_max( 'H', 2.0 );
		small_mover->angle_max( 'E', 5.0 );
		small_mover->angle_max( 'L', 6.0 );
		shear_mover->angle_max( 'H', 2.0 );
		shear_mover->angle_max( 'E', 5.0 );
		shear_mover->angle_max( 'L', 6.0 );
	}

	CCDLoopClosureMoverOP ccd_moves = new CCDLoopClosureMover( one_loop, cdrh3_map );
	RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);

	SequenceMoverOP wiggle_cdr_h3( new SequenceMover() );
	wiggle_cdr_h3->add_mover( one_loop_fold_tree );
	wiggle_cdr_h3->add_mover( small_mover );
	wiggle_cdr_h3->add_mover( shear_mover );
	wiggle_cdr_h3->add_mover( ccd_cycle );
	wiggle_cdr_h3->add_mover( with_flank_fold_tree );


	cdrh3_map->set_bb( flank_allow_bb_move );
	with_flank_fold_tree->apply( pose_in );
	loop_min_mover->movemap( cdrh3_map );
	loop_min_mover->apply( pose_in );
	cdrh3_map->set_bb( allow_bb_move );

	// rotamer trials
	select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
	                      allow_repack);
	cdrh3_map->set_chi( allow_repack );
	setup_packer_task( start_pose_ );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( new RestrictToInterface( allow_repack ) );
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial = new protocols::simple_moves::RotamerTrialsMover(
	    highres_scorefxn_, tf_ );

	pack_rottrial->apply( pose_in );


	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	MonteCarloOP mc;
	mc = new moves::MonteCarlo( pose_in, *highres_scorefxn_, temperature );
	mc->reset( pose_in ); // monte carlo reset

	bool relaxed_H3_found_ever( false );
	if( H3_filter_)
		relaxed_H3_found_ever =CDR_H3_filter(pose_in,antibody_in_.cdrh_[3][1],
		                                     (antibody_in_.cdrh_[3][2] - antibody_in_.cdrh_[3][1]) + 1 );

	// outer cycle
	for(Size i = 1; i <= outer_cycles; i++) {
		mc->recover_low( pose_in );
		Size h3_attempts(0); // number of H3 checks after refinement

		// inner cycle
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc->set_temperature( temperature );
			wiggle_cdr_h3->apply( pose_in );
			cdrh3_map->set_bb( flank_allow_bb_move );
			loop_min_mover->movemap( cdrh3_map );
			loop_min_mover->apply( pose_in );
			cdrh3_map->set_bb( allow_bb_move );

			// rotamer trials
			select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
			                      allow_repack);
			cdrh3_map->set_chi( allow_repack );
			setup_packer_task( start_pose_ );
			( *highres_scorefxn_ )( pose_in );
			tf_->push_back( new RestrictToInterface( allow_repack ) );
			protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial = new protocols::simple_moves::RotamerTrialsMover(
			    highres_scorefxn_, tf_ );
			pack_rottrial->apply( pose_in );

			bool relaxed_H3_found_current(false);
			// H3 filter check
			if(H3_filter_ && (h3_attempts <= inner_cycles)) {
				h3_attempts++;
				relaxed_H3_found_current = CDR_H3_filter(pose_in,
				                           antibody_in_.cdrh_[3][1], ( antibody_in_.cdrh_[3][2] -
				                                   antibody_in_.cdrh_[3][1]) + 1 );

				if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
					mc->boltzmann( pose_in );
				} else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
					relaxed_H3_found_ever = true;
					mc->reset( pose_in );
				} else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
					--j;
					continue;
				} else if( relaxed_H3_found_ever && relaxed_H3_found_current )
					mc->boltzmann( pose_in );
			} else {
				if( H3_filter_ ) {
					bool relaxed_H3_found_current(false);
					relaxed_H3_found_current = CDR_H3_filter(pose_in,
					                           antibody_in_.cdrh_[3][1], ( antibody_in_.cdrh_[3][2] -
					                                   antibody_in_.cdrh_[3][1]) + 1 );
					if( !relaxed_H3_found_ever && !relaxed_H3_found_current) {
						mc->boltzmann( pose_in );
					} else if( !relaxed_H3_found_ever && relaxed_H3_found_current ) {
						relaxed_H3_found_ever = true;
						mc->reset( pose_in );
					} else if( relaxed_H3_found_ever && !relaxed_H3_found_current ) {
						mc->recover_low( pose_in );
					} else if( relaxed_H3_found_ever && relaxed_H3_found_current )
						mc->boltzmann( pose_in );
				} else
					mc->boltzmann( pose_in );
			}

			if ( numeric::mod(j,Size(20))==0 || j==inner_cycles ) {
				// repack trial
				loop_repack = new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ );
				setup_packer_task( start_pose_ );
				( *highres_scorefxn_ )( pose_in );
				tf_->push_back( new RestrictToInterface( allow_repack ) );
				loop_repack->task_factory( tf_ );
				loop_repack->apply( pose_in );
				mc->boltzmann( pose_in );
			}
		} // inner cycles
	} // outer cycles
	mc->recover_low( pose_in );

	// minimize
	if( !benchmark_ ) {
		cdrh3_map->set_bb( flank_allow_bb_move );
		with_flank_fold_tree->apply( pose_in );
		loop_min_mover->movemap( cdrh3_map );
		loop_min_mover->apply( pose_in );
		cdrh3_map->set_bb( allow_bb_move );
	}

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	TR << "H3M Finished Relaxing CDR H3 Loop" << std::endl;

	return;
} //  CDRH3Modeler::loop_fa_relax

///////////////////////////////////////////////////////////////////////////
/// @begin loop_centroid_relax
///
/// @brief actually relaxes the region specified
///
/// @detailed This is all done in low resolution. Intention was to give
///           camelid CDR H1 a larger perturbation.
///
/// @param[in] pose, loop begin position, loop end position
///
/// @global_read none
///
/// @global_write none
///
/// @remarks
///
/// @references
///
/// @authors Aroop 05/07/2010
///
/// @last_modified 05/07/2010
///////////////////////////////////////////////////////////////////////////
void
CDRH3Modeler::loop_centroid_relax( pose::Pose & pose_in, Size const loop_begin, Size const loop_end )
{
	using namespace protocols;
	using namespace protocols::simple_moves;
	using namespace protocols::loops;
	using namespace protocols::moves;
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;
	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	TR << "H3M Centroid Relaxing Loop" << std::endl;

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	//setting MoveMap
	kinematics::MoveMapOP loop_map;
	loop_map = new kinematics::MoveMap();
	loop_map->clear();
	loop_map->set_chi( false );
	loop_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.total_residue(), false );
	for( Size ii = loop_begin; ii <= loop_end; ii++ )
		allow_bb_move[ ii ] = true;
	loop_map->set_bb( allow_bb_move );
	loop_map->set_jump( 1, false );


	Size loop_size = ( loop_end - loop_begin ) + 1;
	Size cutpoint = loop_begin + Size(loop_size/2);

	loops::Loop one_loop( loop_begin, loop_end,	cutpoint,	0, false );
	simple_one_loop_fold_tree( pose_in, one_loop );

	// set cutpoint variants for correct chainbreak scoring
	if( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
		if( !pose_in.residue( cutpoint ).has_variant_type(
		            chemical::CUTPOINT_LOWER))
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_LOWER,
			        cutpoint );
		if( !pose_in.residue( cutpoint + 1 ).has_variant_type(
		            chemical::CUTPOINT_UPPER ) )
			core::pose::add_variant_type_to_pose_residue( pose_in,
			        chemical::CUTPOINT_UPPER,
			        cutpoint + 1 );
	}



	Real min_tolerance = 0.001;
	if( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "dfpmin_armijo_nonmonotone" );
	bool nb_list = true;
	protocols::simple_moves::MinMoverOP loop_min_mover = new protocols::simple_moves::MinMover( loop_map,
	        lowres_scorefxn_, min_type, min_tolerance, nb_list );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	Size inner_cycles( loop_size );
	Size outer_cycles( 1 );
	if( antibody_refine_ || refine_input_loop_ )
		outer_cycles = 5;
	if( antibody_refine_ && snug_fit_ )
		outer_cycles = 2;
	if( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles = 1;
		outer_cycles = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
	protocols::simple_moves::BackboneMoverOP small_mover = new protocols::simple_moves::SmallMover( loop_map,
	        high_move_temp,
	        n_small_moves );
	protocols::simple_moves::BackboneMoverOP shear_mover = new protocols::simple_moves::ShearMover( loop_map,
	        high_move_temp,
	        n_small_moves );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );

	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CCDLoopClosureMoverOP ccd_moves = new CCDLoopClosureMover( one_loop, loop_map );
	RepeatMoverOP ccd_cycle = new RepeatMover(ccd_moves, n_small_moves);

	SequenceMoverOP wiggle_cdr_h3( new SequenceMover() );
	wiggle_cdr_h3->add_mover( small_mover );
	wiggle_cdr_h3->add_mover( shear_mover );
	wiggle_cdr_h3->add_mover( ccd_cycle );


	loop_min_mover->apply( pose_in );

	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	MonteCarloOP mc;
	mc = new moves::MonteCarlo( pose_in, *lowres_scorefxn_, temperature );
	mc->reset( pose_in ); // monte carlo reset

	// outer cycle
	for(Size i = 1; i <= outer_cycles; i++) {
		mc->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc->set_temperature( temperature );
			wiggle_cdr_h3->apply( pose_in );
			loop_min_mover->apply( pose_in );

			mc->boltzmann( pose_in );

		} // inner cycles
	} // outer cycles
	mc->recover_low( pose_in );

	// minimize
	if( !benchmark_ )
		loop_min_mover->apply( pose_in );

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	TR << "H3M Finished Centroid Relaxing Loop" << std::endl;

	return;
} // loop_centroid_relax

void
CDRH3Modeler::setup_packer_task(
    pose::Pose & pose_in ) {
	using namespace pack::task;
	using namespace pack::task::operation;

	if( init_task_factory_ ) {
		tf_ = new TaskFactory( *init_task_factory_ );
		TR << "CDRH3Modeler Reinitializing Packer Task" << std::endl;
		return;
	} else
		tf_ = new TaskFactory;

	TR << "CDRH3Modeler Setting Up Packer Task" << std::endl;

	tf_->push_back( new OperateOnCertainResidues( ResLvlTaskOperationOP( new PreventRepackingRLT ), ResFilterOP( new ResidueLacksProperty("PROTEIN") ) ) );
	tf_->push_back( new InitializeFromCommandline );
	tf_->push_back( new IncludeCurrent );
	tf_->push_back( new RestrictToRepacking );
	tf_->push_back( new NoRepackDisulfides );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot =
	    new pack::rotamer_set::UnboundRotamersOperation();
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation =
	    new operation::AppendRotamerSet( unboundrot );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	init_task_factory_ = tf_;

	TR << "CDRH3Modeler Done: Setting Up Packer Task" << std::endl;

} // setup_packer_task



}  // namespace antibody
}  // namespace protocols
