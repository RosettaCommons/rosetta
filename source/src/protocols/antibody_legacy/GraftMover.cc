// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file moves/GraftMover.cc
/// @brief grafts a cdr onto the template of an antibody framework
/// @details
/// @author Aroop Sircar (aroopsircar@yahoo.com)

// Unit headers
#include <protocols/antibody_legacy/GraftMover.hh>

// Package headers
#include <protocols/antibody_legacy/CDRH3Modeler.hh>

// Rosetta Headers
#include <core/chemical/ChemicalManager.fwd.hh>

#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/types.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>
#include <core/pack/task/operation/OperateOnCertainResidues.hh>
#include <core/pack/task/operation/ResFilters.hh>
#include <core/pack/task/operation/ResLvlTaskOperations.hh>
#include <protocols/toolbox/task_operations/RestrictToInterface.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
#include <basic/Tracer.hh>

#include <protocols/antibody_legacy/AntibodyClass.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
//#include <protocols/loops/LoopMover.fwd.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/RotamerTrialsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <utility/exit.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/variant_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.GraftMover" );
static THREAD_LOCAL basic::Tracer TRO( "protocols.antibody.GraftOneMover" );

namespace protocols {
namespace antibody_legacy {

using namespace core;
using namespace moves;

GraftMover::GraftMover():Mover( "GraftMover" ) {
	set_default();
} // GraftMover default constructor

// GraftMover default destructor
GraftMover::~GraftMover() = default;

void GraftMover::set_default() {
	TR <<  "Setting up default settings" << std::endl;
	graft_l1_ = false;
	graft_l2_ = false;
	graft_l3_ = false;
	graft_h1_ = false;
	graft_h2_ = false;
	graft_h3_ = false;
	camelid_ = false;
} // GraftMover set_default

void GraftMover::apply( pose::Pose & pose_in ) {
	TR <<  "Grafting designated CDRs" << std::endl;

	Antibody antibody_in( pose_in, camelid_ );
	SequenceMoverOP graft_sequence( new SequenceMover() );
	Size nres = pose_in.size();

	// Storing secondary structure
	utility::vector1<char> secondary_struct_storage;
	for ( Size i = 1; i <= nres; i++ ) {
		secondary_struct_storage.push_back( pose_in.secstruct( i ) );
	}

	if ( !camelid_ ) {
		if ( graft_l1_ ) {
			GraftOneMoverOP graftone_l1( new GraftOneMover(
				antibody_in.cdrl_[1][1],antibody_in.cdrl_[1][2],"l1" ) );
			graftone_l1->enable_benchmark_mode( benchmark_ );
			graft_sequence->add_mover( graftone_l1 );
		}
		if ( graft_l2_ ) {
			GraftOneMoverOP graftone_l2( new GraftOneMover(
				antibody_in.cdrl_[2][1],antibody_in.cdrl_[2][2],"l2" ) );
			graftone_l2->enable_benchmark_mode( benchmark_ );
			graft_sequence->add_mover( graftone_l2 );
		}
		if ( graft_l3_ ) {
			GraftOneMoverOP graftone_l3( new GraftOneMover(
				antibody_in.cdrl_[3][1],antibody_in.cdrl_[3][2],"l3" ) );
			graftone_l3->enable_benchmark_mode( benchmark_ );
			graft_sequence->add_mover( graftone_l3 );
		}
	}
	if ( graft_h1_ ) {
		GraftOneMoverOP graftone_h1( new GraftOneMover(
			antibody_in.cdrh_[1][1],antibody_in.cdrh_[1][2],"h1" ) );
		graftone_h1->enable_benchmark_mode( benchmark_ );
		graft_sequence->add_mover( graftone_h1 );
	}
	if ( graft_h2_ ) {
		GraftOneMoverOP graftone_h2( new GraftOneMover(
			antibody_in.cdrh_[2][1],antibody_in.cdrh_[2][2],"h2" ) );
		graftone_h2->enable_benchmark_mode( benchmark_ );
		graft_sequence->add_mover( graftone_h2 );
	}
	if ( graft_h3_ ) {
		GraftOneMoverOP graftone_h3( new GraftOneMover(
			antibody_in.cdrh_[3][1],antibody_in.cdrh_[3][2],"h3" ) );
		graftone_h3->enable_benchmark_mode( benchmark_ );
		graft_sequence->add_mover( graftone_h3 );
	}

	graft_sequence->apply(pose_in);

	if ( !graft_h3_ ) {
		TR << "Extending CDR H3" << std::endl;

		Size framework_loop_begin( antibody_in.cdrh_[3][1] - 1 );
		Size frmrk_loop_end_plus_one( antibody_in.cdrh_[3][2] + 1 );
		Size cutpoint = framework_loop_begin + 1;
		loops::Loop cdr_h3( framework_loop_begin, frmrk_loop_end_plus_one,
			cutpoint, 0, false );
		simple_one_loop_fold_tree( pose_in, cdr_h3);

		// silly hack to make extended loops work
		loops::LoopsOP loop_list( new loops::Loops() );
		loop_list->add_loop( cdr_h3 );
		/* Commented out by BDW with JX's consent
		loops::loop_mover::LoopMoverOP my_loop_move = new loops::loop_mover::LoopMover( loop_list );
		my_loop_move->set_extended_torsions( pose_in, cdr_h3 );
		*/
	}

	if ( graft_l1_ || graft_l2_ || graft_l3_ ||
			graft_h1_ || graft_h2_ || graft_h3_ ) {
		// disallowing bogus sidechains while repacking using include_current
		utility::vector1<bool> allow_repack( nres, false );
		for ( Size i = 1; i <= nres; i++ ) {
			if ( pose_in.secstruct(i) == 'X' ) {
				allow_repack[i] = true;
			}
			if ( !graft_h2_ && ( i >= antibody_in.cdrh_[2][1] ) &&
					( i <= antibody_in.cdrh_[2][2] ) ) {
				allow_repack[i] = true;
			}
			if ( !graft_h3_ && ( i >= antibody_in.cdrh_[3][1] ) &&
					( i <= antibody_in.cdrh_[3][2] ) ) {
				allow_repack[i] = true;
			}
		}

		// Recover secondary structures
		for ( Size i = 1; i <= nres; i++ ) {
			pose_in.set_secstruct( i, secondary_struct_storage[ i ] );
		}
		// saving sidechains for use of include_current
		pose::Pose saved_sidechains( pose_in );

		// generating centroids for residues devoid of sidechains
		simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
		simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );
		to_centroid.apply( pose_in );
		to_full_atom.apply( pose_in );
		//recover sidechains from starting structures
		protocols::simple_moves::ReturnSidechainMover recover_sidechains(
			saved_sidechains );
		recover_sidechains.apply( pose_in );
		bool include_current( false );
		// High resolution scores
		scoring::ScoreFunctionOP antibody_score( scoring::get_score_function_legacy( core::scoring::PRE_TALARIS_2013_STANDARD_WTS ) );
		set_packer_default( pose_in, antibody_score, include_current );
		packer_->apply( pose_in );

		// Retain coordinates for which full atomic coordinates were
		// available in the input structure. Even for template residues
		// which match corresponding target residues, keep the template
		// sidechians. For the residues which do not have any sidechain
		// coordinates, grab them from a repacked version of the pose
		// generated without using include_current
		for ( Size i = 1; i <= nres; i++ ) {
			conformation::Residue const & original_rsd(
				saved_sidechains.residue( i ) );
			conformation::Residue const & packed_rsd(
				pose_in.residue( i ) );
			if ( !allow_repack[i] ) {
				for ( Size j=1; j <= packed_rsd.natoms(); ++j ) {
					std::string const & atom_name( packed_rsd.atom_name(j) );
					if ( original_rsd.type().has( atom_name ) ) {
						pose_in.set_xyz( id::AtomID( packed_rsd.atom_index(
							atom_name),i), original_rsd.xyz( atom_name ) );
					}
				}
			}
		}
		include_current = true;
		set_packer_default( pose_in, antibody_score, include_current );
		packer_->apply( pose_in );

		SequenceMoverOP close_grafted_loops( new SequenceMover() );
		if ( !camelid_ ) {
			if ( graft_l1_ ) {
				CloseOneMoverOP closeone_l1( new CloseOneMover(
					antibody_in.cdrl_[1][1],antibody_in.cdrl_[1][2]) );
				closeone_l1->enable_benchmark_mode( benchmark_ );
				close_grafted_loops->add_mover( closeone_l1 );
			}
			if ( graft_l2_ ) {
				CloseOneMoverOP closeone_l2( new CloseOneMover(
					antibody_in.cdrl_[2][1],antibody_in.cdrl_[2][2]) );
				closeone_l2->enable_benchmark_mode( benchmark_ );
				close_grafted_loops->add_mover( closeone_l2 );
			}
			if ( graft_l3_ ) {
				CloseOneMoverOP closeone_l3( new CloseOneMover(
					antibody_in.cdrl_[3][1],antibody_in.cdrl_[3][2]) );
				closeone_l3->enable_benchmark_mode( benchmark_ );
				close_grafted_loops->add_mover( closeone_l3 );
			}
		}
		if ( graft_h1_ ) {
			CloseOneMoverOP closeone_h1( new CloseOneMover(
				antibody_in.cdrh_[1][1],antibody_in.cdrh_[1][2]) );
			closeone_h1->enable_benchmark_mode( benchmark_ );
			close_grafted_loops->add_mover( closeone_h1 );
		}
		if ( graft_h2_ ) {
			CloseOneMoverOP closeone_h2( new CloseOneMover(
				antibody_in.cdrh_[2][1],antibody_in.cdrh_[2][2]) );
			closeone_h2->enable_benchmark_mode( benchmark_ );
			close_grafted_loops->add_mover( closeone_h2 );
		}
		if ( graft_h3_ ) {
			CloseOneMoverOP closeone_h3( new CloseOneMover(
				antibody_in.cdrh_[3][1],antibody_in.cdrh_[3][2]) );
			closeone_h3->enable_benchmark_mode( benchmark_ );
			close_grafted_loops->add_mover( closeone_h3 );
		}
		close_grafted_loops->apply(pose_in);

		for ( Size i = 1; i <= nres; i++ ) {
			if ( pose_in.secstruct(i) == 'Y' ) {
				allow_repack[i] = true;
			}
		}

		to_centroid.apply( pose_in );
		to_full_atom.apply( pose_in );
		recover_sidechains.apply( pose_in );
		set_packer_default( pose_in, antibody_score, include_current );
		packer_->apply( pose_in );

		// Retain coordinates for which full atomic coordinates were
		// available in the input structure. Even for template residues
		// which match corresponding target residues, keep the template
		// sidechians. For the residues which do not have any sidechain
		// coordinates, grab them from a repacked version of the pose
		// generated without using include_current
		for ( Size i = 1; i <= nres; i++ ) {
			conformation::Residue const & original_rsd(
				saved_sidechains.residue( i ) );
			conformation::Residue const & packed_rsd(
				pose_in.residue( i ) );
			if ( !allow_repack[i] ) {
				for ( Size j=1; j <= packed_rsd.natoms(); ++j ) {
					std::string const & atom_name( packed_rsd.atom_name(j) );
					if ( original_rsd.type().has( atom_name ) ) {
						pose_in.set_xyz( id::AtomID( packed_rsd.atom_index(
							atom_name),i), original_rsd.xyz( atom_name ) );
					}
				}
			}
		}
		include_current = true;
		set_packer_default( pose_in, antibody_score, include_current );
		packer_->apply( pose_in ); // packing to account for loop closure

		for ( Size i = 1; i <= nres; i++ ) {
			allow_repack[i] = true;
		}
		include_current = true;
		set_packer_default( pose_in, antibody_score, include_current );
		packer_->apply( pose_in ); // packing all to relieve clashes

		// relax optimized CDR grafted regions
		relax_optimized_CDR_grafts( pose_in );


		// Recover secondary structures
		for ( Size i = 1; i <= nres; i++ ) {
			pose_in.set_secstruct( i, secondary_struct_storage[ i ] );
		}

	}

	// align Fv to native.Fv
	pose::Pose native_pose;
	if ( get_native_pose() ) {
		native_pose = *get_native_pose();
	} else {
		native_pose = pose_in;
	}
	Antibody native_ab( native_pose, camelid_ );
	Antibody antibody_out( pose_in, camelid_ );

	antibody_out.align_to_native( native_ab );
	pose_in = antibody_out.Fv;
} // GraftMover::apply()

std::string
GraftMover::get_name() const {
	return "GraftMover";
}

void GraftMover::set_packer_default(
	pose::Pose & pose_in,
	scoring::ScoreFunctionOP scorefxn,
	bool include_current) {

	//set up packer
	pack::task::PackerTaskOP task;
	task = pack::task::TaskFactory::create_packer_task( pose_in );
	task->restrict_to_repacking();
	task->or_include_current( include_current );
	GraftMover::packer_ = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( scorefxn, task ) );

} // GraftMover set_packer_default

GraftOneMover::GraftOneMover(
	Size query_start,
	Size query_end,
	std::string template_name ) : Mover( "GraftOneMover" ) {
	query_start_ = query_start;
	query_end_ = query_end;
	template_name_ = template_name;
	set_default();
} // GraftOneMover default constructor

void GraftOneMover::set_default() {
	TRO << "Reading in template: " << template_name_ << ".pdb "
		<< std::endl;
	core::import_pose::pose_from_file( template_pose_, template_name_ + ".pdb" , core::import_pose::PDB_file);
	Antibody antibody( template_pose_, template_name_ );
	template_start_ = antibody.current_start;
	template_end_ = antibody.current_end;
} // GraftOneMover::set_default

std::string
GraftOneMover::get_name() const {
	return "GraftOneMover";
}
void GraftOneMover::apply( pose::Pose & pose_in ) {
	Size const nres( pose_in.size() ); // Total residues
	Size query_size = ( query_end_ - query_start_ ) + 1;

	pose::Pose truncated_pose;
	truncated_pose = pose_in;

	// create a sub pose with  5 flanking residues on either side of CDR loop
	truncated_pose.conformation().delete_residue_range_slow(
		query_end_ + 5, nres);
	truncated_pose.conformation().delete_residue_range_slow(
		1, query_start_ - 5);
	truncated_pose.conformation().delete_residue_range_slow(
		5, ( query_size - 1 ) + 5 );

	// create atom map for superimposing 2 flanking resiudes
	id::AtomID_Map< id::AtomID > atom_map;
	core::pose::initialize_atomid_map( atom_map, template_pose_,
		id::AtomID::BOGUS_ATOM_ID() );

	Size flank = 2; // default 2
	for ( Size start_stem = 4 - (flank-1); start_stem <= 4; start_stem++ ) {
		for ( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, start_stem );
			id::AtomID const id2( j, start_stem );
			atom_map[ id1 ] = id2;
		}
	}
	for ( Size end_stem = 4 + query_size + 1; end_stem <= 4+query_size+flank;
			end_stem++ ) {
		for ( Size j=1; j <= 4; j++ ) {
			id::AtomID const id1( j, end_stem );
			id::AtomID const id2( j, end_stem - query_size );
			atom_map[ id1 ] = id2;
		}
	}
	scoring::superimpose_pose( template_pose_, truncated_pose,
		atom_map );


	// copy coordinates of properly oriented template to framework
	for ( Size i = query_start_; i <= query_end_; ++i ) {
		conformation::Residue const & orientable_rsd(pose_in.residue(i));
		conformation::Residue const & template_rsd( template_pose_.
			residue( i - (query_start_ - 5) ) );
		// keeping track of missing rotamers
		if ( template_rsd.name() != orientable_rsd.name() ) {
			pose_in.set_secstruct( i, 'X' );
		}
		for ( Size j=1; j <= template_rsd.natoms(); ++j ) {
			std::string const & atom_name( template_rsd.atom_name(j) );
			if ( orientable_rsd.type().has( atom_name ) ) {
				pose_in.set_xyz( id::AtomID( orientable_rsd.atom_index(
					atom_name),i), template_rsd.xyz(
					atom_name ) );
			}
			// hack for generating coordinates of amide hydrogen for target
			// residues whose template was a Proline residue (Proline residues
			// do not have the amide hydrogen found ubiquitously in all other
			// amino acids)
			if ( template_rsd.name() == "PRO" && orientable_rsd.name() != "PRO"
					&& atom_name == " CD " ) {
				id::AtomID_Mask missing( false );
				// dimension the missing-atom mask
				core::pose::initialize_atomid_map( missing, pose_in );
				missing[ id::AtomID(
					pose_in.residue_type(i).atom_index("H"), i ) ] = true;
				pose_in.conformation().fill_missing_atoms( missing );
			}
		}
	}
} // GraftOneMover::apply

void GraftMover::relax_optimized_CDR_grafts( pose::Pose & pose_in ) {
	Size loop_begin(0), loop_end(0);
	bool detect_flag( false );
	for ( Size ii = 1; ii <= pose_in.size(); ii++ ) {
		if ( (pose_in.secstruct(ii) == 'Y') && !detect_flag ) {
			loop_begin = ii;
			detect_flag = true;
		}
		if ( (pose_in.secstruct(ii) != 'Y') && detect_flag ) {
			loop_end = ii - 1;
			detect_flag = false;
		}
		if ( (detect_flag == false) && (loop_begin != 0) && (loop_end != 0 ) ) {
			LoopRlxMoverOP rlx_one_loop( new LoopRlxMover( loop_begin, loop_end) );
			rlx_one_loop->enable_benchmark_mode( benchmark_ );
			rlx_one_loop->apply( pose_in );
			loop_begin = 0;
			loop_end = 0;
		}
	} // for ii <= nres
} // GraftMover::relax_optimized_CDR_grafts


CloseOneMover::CloseOneMover(
	Size query_start,
	Size query_end
)
: Mover( "CloseOneMover" ) {
	loop_start_ = query_start;
	loop_end_ = query_end;
	set_default();
} // CloseOneMover default constructor

void CloseOneMover::set_default() {
	allowed_separation_ = 1.9;
	flanking_residues_ = 5; // default 5;
} // CloseOneMover::set_default

void CloseOneMover::apply( pose::Pose & pose_in ) {

	Size const N ( 1 ); // N atom
	Size const C ( 3 ); // C atom

	// Coordinates of the C and N atoms at stem
	numeric::xyzVector_float peptide_C, peptide_N;
	// N-terminal
	peptide_C = pose_in.residue( loop_start_ - 1 ).xyz( C );
	peptide_N = pose_in.residue( loop_start_ ).xyz( N );
	Real nter_separation=peptide_C.distance(peptide_N);
	// C-terminal
	peptide_C = pose_in.residue( loop_end_ ).xyz( C );
	peptide_N = pose_in.residue( loop_end_ + 1 ).xyz( N );
	Real cter_separation=peptide_C.distance(peptide_N);

	// switching to centroid mode
	simple_moves::SwitchResidueTypeSetMover to_centroid( chemical::CENTROID );
	simple_moves::SwitchResidueTypeSetMover to_full_atom( chemical::FA_STANDARD );

	bool repack_flag( false );
	if ( ( nter_separation > allowed_separation_ ) ||
			( cter_separation > allowed_separation_ ) ) {
		to_centroid.apply( pose_in );
		repack_flag = true;
	}

	// saving sidechains for use of include_current
	pose::Pose saved_sc( pose_in );
	//recover sidechains from starting structures
	protocols::simple_moves::ReturnSidechainMover recover_sidechains( saved_sc );

	bool nter( false );
	if ( nter_separation > allowed_separation_ ) {
		nter = true;
		close_one_loop_stem( pose_in, loop_start_, nter );
	}

	if ( cter_separation > allowed_separation_ ) {
		nter = false;
		close_one_loop_stem( pose_in, loop_end_, nter );
	}

	Real separation = 0.00;
	for ( Size ii = loop_start_; ii <= loop_end_; ii++ ) {
		peptide_C = pose_in.residue( ii ).xyz( C );
		peptide_N = pose_in.residue( ii + 1 ).xyz( N );
		separation=peptide_C.distance(peptide_N);
		if ( separation > allowed_separation_ ) {
			repack_flag = true;
			Size cutpoint = ii;
			close_one_loop_stem( pose_in, loop_start_, loop_end_, cutpoint );
		}
	} // for

	if ( repack_flag ) {
		to_full_atom.apply( pose_in );
		recover_sidechains.apply( pose_in );
	}

	return;
} // CloseOneMover::apply

std::string
CloseOneMover::get_name() const {
	return "CloseOneMover";
}

void CloseOneMover::close_one_loop_stem (
	pose::Pose & pose_in,
	Size cutpoint_in,
	bool nter
) {
	using namespace protocols;
	using namespace protocols::loops;

	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;


	// storing starting fold tree
	//kinematics::FoldTree tree_in( pose_in.fold_tree() );

	Size loop_flex_begin, loop_flex_end;
	if ( nter ) {
		loop_flex_begin = cutpoint_in;
		loop_flex_end = cutpoint_in + ( flanking_residues_ - 1 );
	} else {
		loop_flex_begin = cutpoint_in - ( flanking_residues_ - 1 );
		loop_flex_end = cutpoint_in;
	}

	//setting MoveMap
	kinematics::MoveMapOP loop_map;
	loop_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	loop_map->clear();
	loop_map->set_chi( false );
	loop_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.size(), false );
	for ( Size ii = loop_flex_begin; ii <= loop_flex_end; ii++ ) {
		allow_bb_move[ ii ] = true;
		pose_in.set_secstruct( ii, 'Y' );
	}
	loop_map->set_bb( allow_bb_move );
	loop_map->set_jump( 1, false );

	Size loop_begin, loop_end, cutpoint;
	if ( nter ) {
		loop_begin = loop_flex_begin - 2;
		loop_end = loop_flex_end;
		cutpoint = loop_flex_begin - 1;
	} else {
		loop_begin = loop_flex_begin;
		loop_end = loop_flex_end + 2;
		cutpoint = loop_flex_end;
	}

	close_one_loop_stem_helper( loop_begin, loop_end, cutpoint, pose_in, loop_map );

} // CloseOneMover::close_one_loop_stem

void CloseOneMover::close_one_loop_stem (
	pose::Pose & pose_in,
	Size loop_begin,
	Size loop_end,
	Size cutpoint
) {
	using namespace protocols;
	using namespace protocols::loops;

	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;


	// storing starting fold tree
	//kinematics::FoldTree tree_in( pose_in.fold_tree() );

	//setting MoveMap
	kinematics::MoveMapOP loop_map;
	loop_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	loop_map->clear();
	loop_map->set_chi( false );
	loop_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.size(), false );
	for ( Size ii = loop_begin; ii <= loop_end; ii++ ) {
		allow_bb_move[ ii ] = true;
		pose_in.set_secstruct( ii, 'Y' );
	}
	loop_map->set_bb( allow_bb_move );
	loop_map->set_jump( 1, false );

	close_one_loop_stem_helper( loop_begin, loop_end, cutpoint, pose_in, loop_map );

} // CloseOneMover::close_one_loop_stem

void
CloseOneMover::close_one_loop_stem_helper(
	Size loop_begin, Size loop_end, Size cutpoint, Pose & pose_in, kinematics::MoveMapOP loop_map
) {

	using namespace protocols;
	using namespace protocols::loops;

	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	loops::Loop one_loop( loop_begin, loop_end, cutpoint, 0, false );
	simple_one_loop_fold_tree( pose_in, one_loop );

	// set cutpoint variants for correct chainbreak scoring
	if ( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
		if ( !pose_in.residue( cutpoint ).has_variant_type(
				chemical::CUTPOINT_LOWER) ) {
			core::pose::add_variant_type_to_pose_residue( pose_in,
				chemical::CUTPOINT_LOWER,
				cutpoint );
		}
		if ( !pose_in.residue( cutpoint + 1 ).has_variant_type(
				chemical::CUTPOINT_UPPER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose_in,
				chemical::CUTPOINT_UPPER,
				cutpoint + 1 );
		}
	}

	scoring::ScoreFunctionOP lowres_scorefxn;
	lowres_scorefxn = scoring::ScoreFunctionFactory::
		create_score_function( "cen_std", "score4L" );
	lowres_scorefxn->set_weight( scoring::chainbreak, 10. / 3. );

	Real min_tolerance = 0.001;
	if ( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "lbfgs_armijo_nonmonotone" );
	bool nb_list = true;
	protocols::simple_moves::MinMoverOP loop_min_mover(
		new protocols::simple_moves::MinMover( loop_map,
		lowres_scorefxn, min_type, min_tolerance, nb_list ) );

	// more params
	Size loop_size = ( loop_end - loop_begin ) + 1;
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	Size inner_cycles( loop_size );
	Size outer_cycles( 2 );
	if ( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles = 1;
		outer_cycles = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
	protocols::simple_moves::BackboneMoverOP small_mover(
		new protocols::simple_moves::SmallMover( loop_map,
		high_move_temp,
		n_small_moves ) );
	protocols::simple_moves::BackboneMoverOP shear_mover(
		new protocols::simple_moves::ShearMover( loop_map,
		high_move_temp,
		n_small_moves ) );

	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CCDLoopClosureMoverOP ccd_moves( new CCDLoopClosureMover( one_loop, loop_map ) );
	RepeatMoverOP ccd_cycle( new RepeatMover(ccd_moves, n_small_moves) );

	SequenceMoverOP wiggle_cdr( new SequenceMover() );
	wiggle_cdr->add_mover( small_mover );
	wiggle_cdr->add_mover( shear_mover );
	wiggle_cdr->add_mover( ccd_cycle );


	loop_min_mover->apply( pose_in );

	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *lowres_scorefxn, temperature ) );
	mc->reset( pose_in ); // monte carlo reset

	// outer cycle
	for ( Size i = 1; i <= outer_cycles; i++ ) {
		mc->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc->set_temperature( temperature );
			wiggle_cdr->apply( pose_in );
			loop_min_mover->apply( pose_in );

			mc->boltzmann( pose_in );

		} // inner cycles
	} // outer cycles
	mc->recover_low( pose_in );

	// minimize
	if ( !benchmark_ ) {
		loop_min_mover->apply( pose_in );
	}

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	return;
}

LoopRlxMover::LoopRlxMover(
	Size query_start,
	Size query_end
) : Mover( "LoopRlxMover" ) {
	loop_start_ = query_start;
	loop_end_ = query_end;
	set_default();
} // LoopRlxMover default constructor

void LoopRlxMover::set_default() {
	highres_scorefxn_ = scoring::get_score_function();
	highres_scorefxn_->set_weight( scoring::chainbreak, 1.0 );
	highres_scorefxn_->set_weight( scoring::overlap_chainbreak, 10./3. );
	return;
} // LoopRlxMover::set_default


///////////////////////////////////////////////////////////////////////////
///
/// @brief relaxes the region specified
///
/// @details This is all done in high resolution.Hence there are no rigid
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
/// @author Aroop 05/12/2010
///
///////////////////////////////////////////////////////////////////////////
void LoopRlxMover::apply( pose::Pose & pose_in ) {
	using namespace protocols;
	using namespace protocols::loops;
	using namespace protocols::moves;
	using namespace pack;
	using namespace pack::task;
	using namespace pack::task::operation;

	using loop_closure::ccd::CCDLoopClosureMover;
	using loop_closure::ccd::CCDLoopClosureMoverOP;

	TR << "LoopRlxMover: Apply" << std::endl;

	// storing starting fold tree
	kinematics::FoldTree tree_in( pose_in.fold_tree() );

	//setting MoveMap
	kinematics::MoveMapOP loop_map;
	loop_map = kinematics::MoveMapOP( new kinematics::MoveMap() );
	loop_map->clear();
	loop_map->set_chi( false );
	loop_map->set_bb( false );
	utility::vector1< bool> allow_bb_move( pose_in.size(), false );
	for ( Size ii = loop_start_; ii <= loop_end_; ii++ ) {
		allow_bb_move[ ii ] = true;
	}
	loop_map->set_bb( allow_bb_move );
	loop_map->set_jump( 1, false );


	Size loop_size = ( loop_end_ - loop_start_ ) + 1;
	Size cutpoint = loop_start_ + int(loop_size/2);

	loops::Loop one_loop( loop_start_, loop_end_, cutpoint, 0, false );
	simple_one_loop_fold_tree( pose_in, one_loop );

	if ( !pose_in.is_fullatom() ) {
		utility_exit_with_message("Fullatom poses only");
	}

	// set cutpoint variants for correct chainbreak scoring
	if ( !pose_in.residue( cutpoint ).is_upper_terminus() ) {
		if ( !pose_in.residue( cutpoint ).has_variant_type(
				chemical::CUTPOINT_LOWER) ) {
			core::pose::add_variant_type_to_pose_residue( pose_in,
				chemical::CUTPOINT_LOWER,
				cutpoint );
		}
		if ( !pose_in.residue( cutpoint + 1 ).has_variant_type(
				chemical::CUTPOINT_UPPER ) ) {
			core::pose::add_variant_type_to_pose_residue( pose_in,
				chemical::CUTPOINT_UPPER,
				cutpoint + 1 );
		}
	}

	utility::vector1< bool> allow_repack( pose_in.size(), false );
	select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
		allow_repack);
	loop_map->set_chi( allow_repack );

	protocols::simple_moves::PackRotamersMoverOP loop_repack( new protocols::simple_moves::PackRotamersMover(highres_scorefxn_) );
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) ) );
	loop_repack->task_factory(tf_);
	loop_repack->apply( pose_in );

	Real min_tolerance = 0.001;
	if ( benchmark_ ) min_tolerance = 1.0;
	std::string min_type = std::string( "lbfgs_armijo_nonmonotone" );
	bool nb_list = true;
	protocols::simple_moves::MinMoverOP loop_min_mover( new protocols::simple_moves::MinMover( loop_map,
		highres_scorefxn_, min_type, min_tolerance, nb_list ) );

	// more params
	Size n_small_moves ( numeric::max(Size(5), Size(loop_size/2)) );
	Size inner_cycles( loop_size );
	Size outer_cycles( 2 );
	if ( benchmark_ ) {
		n_small_moves = 1;
		inner_cycles = 1;
		outer_cycles = 1;
	}

	Real high_move_temp = 2.00;
	// minimize amplitude of moves if correct parameter is set
	protocols::simple_moves::BackboneMoverOP small_mover( new protocols::simple_moves::SmallMover( loop_map,
		high_move_temp,
		n_small_moves ) );
	protocols::simple_moves::BackboneMoverOP shear_mover( new protocols::simple_moves::ShearMover( loop_map,
		high_move_temp,
		n_small_moves ) );
	small_mover->angle_max( 'H', 2.0 );
	small_mover->angle_max( 'E', 5.0 );
	small_mover->angle_max( 'L', 6.0 );
	shear_mover->angle_max( 'H', 2.0 );
	shear_mover->angle_max( 'E', 5.0 );
	shear_mover->angle_max( 'L', 6.0 );

	CCDLoopClosureMoverOP ccd_moves( new CCDLoopClosureMover( one_loop, loop_map ) );
	RepeatMoverOP ccd_cycle( new RepeatMover(ccd_moves, n_small_moves) );

	SequenceMoverOP wiggle_loop( new SequenceMover() );
	wiggle_loop->add_mover( small_mover );
	wiggle_loop->add_mover( shear_mover );
	wiggle_loop->add_mover( ccd_cycle );

	// rotamer trials
	select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
		allow_repack);
	loop_map->set_chi( allow_repack );
	setup_packer_task( pose_in );
	( *highres_scorefxn_ )( pose_in );
	tf_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) ) );
	protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
		highres_scorefxn_, tf_ ) );

	pack_rottrial->apply( pose_in );


	Real const init_temp( 2.0 );
	Real const last_temp( 0.5 );
	Real const gamma = std::pow( (last_temp/init_temp), (1.0/inner_cycles));
	Real temperature = init_temp;

	MonteCarloOP mc;
	mc = MonteCarloOP( new moves::MonteCarlo( pose_in, *highres_scorefxn_, temperature ) );
	mc->reset( pose_in ); // monte carlo reset

	// outer cycle
	for ( Size i = 1; i <= outer_cycles; i++ ) {
		mc->recover_low( pose_in );

		// inner cycle
		for ( Size j = 1; j <= inner_cycles; j++ ) {
			temperature *= gamma;
			mc->set_temperature( temperature );
			wiggle_loop->apply( pose_in );
			loop_min_mover->apply( pose_in );

			// rotamer trials
			select_loop_residues( pose_in, one_loop, true /*include_neighbors*/,
				allow_repack);
			loop_map->set_chi( allow_repack );
			setup_packer_task( pose_in );
			( *highres_scorefxn_ )( pose_in );
			tf_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) ) );
			protocols::simple_moves::RotamerTrialsMoverOP pack_rottrial( new protocols::simple_moves::RotamerTrialsMover(
				highres_scorefxn_, tf_ ) );
			pack_rottrial->apply( pose_in );

			mc->boltzmann( pose_in );


			if ( numeric::mod(j,Size(20))==0 || j==inner_cycles ) {
				// repack trial
				loop_repack = protocols::simple_moves::PackRotamersMoverOP( new protocols::simple_moves::PackRotamersMover( highres_scorefxn_ ) );
				setup_packer_task( pose_in );
				( *highres_scorefxn_ )( pose_in );
				tf_->push_back( TaskOperationCOP( new protocols::toolbox::task_operations::RestrictToInterface( allow_repack ) ) );
				loop_repack->task_factory( tf_ );
				loop_repack->apply( pose_in );
				mc->boltzmann( pose_in );
			}
		} // inner cycles
	} // outer cycles
	mc->recover_low( pose_in );

	// minimize
	if ( !benchmark_ ) {
		loop_min_mover->apply( pose_in );
	}

	// Restoring pose stuff
	pose_in.fold_tree( tree_in ); // Tree

	TR << "LoopRlxMover: Finished Apply" << std::endl;

	return;
} // LoopRlxMover::apply

std::string
LoopRlxMover::get_name() const {
	return "LoopRlxMover";
}

void
LoopRlxMover::setup_packer_task(
	pose::Pose & pose_in
) {
	using namespace pack::task;
	using namespace pack::task::operation;

	if ( init_task_factory_ ) {
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory( *init_task_factory_ ) );
		TR << "LoopRlxMover Reinitializing Packer Task" << std::endl;
		return;
	} else {
		tf_ = core::pack::task::TaskFactoryOP( new TaskFactory );
	}

	TR << "LoopRlxMover Setting Up Packer Task" << std::endl;

	tf_->push_back( TaskOperationCOP( new OperateOnCertainResidues( ResLvlTaskOperationOP( new PreventRepackingRLT ), ResFilterOP( new ResidueLacksProperty("PROTEIN") ) ) ) );
	tf_->push_back( TaskOperationCOP( new InitializeFromCommandline ) );
	tf_->push_back( TaskOperationCOP( new IncludeCurrent ) );
	tf_->push_back( TaskOperationCOP( new RestrictToRepacking ) );
	tf_->push_back( TaskOperationCOP( new NoRepackDisulfides ) );

	// incorporating Ian's UnboundRotamer operation.
	// note that nothing happens if unboundrot option is inactive!
	pack::rotamer_set::UnboundRotamersOperationOP unboundrot( new pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot->initialize_from_command_line();
	operation::AppendRotamerSetOP unboundrot_operation( new operation::AppendRotamerSet( unboundrot ) );
	tf_->push_back( unboundrot_operation );
	// adds scoring bonuses for the "unbound" rotamers, if any
	core::pack::dunbrack::load_unboundrot( pose_in );

	init_task_factory_ = tf_;

	TR << "LoopRlxMover Done: Setting Up Packer Task" << std::endl;

} // LoopRlxMover::setup_packer_task


}  // namespace moves
}  // namespace protocols
