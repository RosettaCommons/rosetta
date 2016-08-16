// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief quick demo program for cm_main
/// @author James Thompson

#include <devel/init.hh>

#include <protocols/relax/MiniRelax.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/ccd_closure.hh>
#include <protocols/comparative_modeling/util.hh>
#include <protocols/comparative_modeling/ThreadingMover.hh>
#include <protocols/idealize/IdealizeMover.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/chemical/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/sequence/util.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/pose/Pose.hh>

#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/import_pose/pose_stream/util.hh>
#include <core/import_pose/pose_stream/MetaPoseInputStream.hh>

#include <core/scoring/ScoreType.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <string>
#include <map>

#include <utility/excn/Exceptions.hh>

void ccd_close_loop(
	core::pose::Pose & pose,
	protocols::loops::Loop & loop
) {
	using core::Real;
	using core::Size;
	using core::kinematics::MoveMap;

	int  const ccd_cycles( 100 ); // num of cycles of ccd_moves
	Real const ccd_tol( 0.01 ); // criterion for a closed loop
	bool const rama_check( false );
	Real const max_rama_score_increase( 2.0 ); // dummy number when rama_check is false
	Real const max_total_delta_helix( 10.0 ); // max overall angle changes for a helical residue
	Real const max_total_delta_strand( 50.0 ); // ... for a residue in strand
	Real const max_total_delta_loop( 75.0 ); // ... for a residue in loop

	core::kinematics::MoveMap mm;
	mm.set_bb(false);
	mm.set_chi(false);
	for ( Size ii = loop.start(); ii <= loop.stop(); ++ii ) {
		mm.set_bb( ii, true );
	}

	// output for ccd_closure
	Real forward_deviation, backward_deviation; // actually loop closure msd, both dirs
	Real torsion_delta, rama_delta; // actually torsion and rama score changes, averaged by loop_size
	protocols::loops::fast_ccd_loop_closure(
		pose, mm, loop.start(), loop.stop(), loop.cut(),
		ccd_cycles, ccd_tol, rama_check, max_rama_score_increase,
		max_total_delta_helix, max_total_delta_strand, max_total_delta_loop,
		forward_deviation, backward_deviation, torsion_delta, rama_delta
	);
}

class FastThreadingMover : public protocols::comparative_modeling::ThreadingMover {

public:
	FastThreadingMover(
		core::sequence::SequenceAlignment aln,
		core::pose::Pose templ
	) : ThreadingMover( aln, templ ) {}

	void apply( core::pose::Pose & query_pose ) {
		// don't rebuild loops
		build_loops( true );
		repack_query( true );
		std::cout << "DEBUG(sequence1) = " << query_pose.sequence() << std::endl;
		ThreadingMover::apply( query_pose );

		//using protocols::loops::Loops;

		//Loops loops = protocols::comparative_modeling::loops_from_alignment(
		//	query_pose.total_residue(), alignment(), min_loop_size()
		//);

		//// foreach loop
		//using core::Size;
		//for ( Size loop_idx = 1; loop_idx <= loops.num_loop(); ++loop_idx ) {
		//	// do a quick ccd closure on the loop residues so long as the loop is not
		//	// at a terminus.
		//	protocols::loops::Loop loop = loops[loop_idx];
		//	if ( loop.start() != 1 ) {
		//		ccd_close_loop(query_pose,loop);
		//	}
		//}
		std::cout << "DEBUG(sequence2) = " << query_pose.sequence() << std::endl;
		query_pose.dump_pdb("after_loop_model.pdb");

		// add coordinate constraints over aligned residues
		using core::conformation::ResidueFactory;
		query_pose.append_residue_by_jump(
			*ResidueFactory::create_residue( query_pose.residue(1).residue_type_set().name_map( "VRT" ) ),
			static_cast< Size > (query_pose.total_residue() / 2)
		);
		Size nres = query_pose.total_residue();

		core::id::SequenceMapping map = alignment().sequence_mapping(1,2);

		using core::Real;
		using core::id::AtomID;
		Real const coord_sdev( 2.0 ); // from idealize -- maybe too small
		for ( Size idx = 1; idx <= nres - 1; ++idx ) {
			if ( map[idx] != 0 ) {
				std::cout << "adding constraints for residue " << idx << std::endl;
				using namespace core::conformation;
				Residue const & nat_i_rsd( query_pose.residue(idx) );
				for ( Size ii = 1; ii <= nat_i_rsd.last_backbone_atom(); ++ii ) {
					using namespace core::scoring::constraints;
					query_pose.add_constraint(
						new CoordinateConstraint(
							AtomID(ii,idx), AtomID(1,nres), nat_i_rsd.xyz(ii),
							new HarmonicFunc(0.0,coord_sdev)
						)
					);
				} // ii
			} // if residue is aligned
		} // nres
		using namespace core::scoring;
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		scorefxn->set_weight( coordinate_constraint, 0.5 );

		// perform a fast repack/minimization on the pose. Maybe do
		// this in a loop.
		using ObjexxFCL::string_of;
		std::cout << "fold_tree pre-min: " << query_pose.fold_tree();

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb(true);
		mm->set_chi(true);
		mm->set_jump( query_pose.num_jump(), true );
		protocols::relax::MiniRelax minimizer( scorefxn );
		minimizer.set_movemap(mm);
		for ( Size ii = 1; ii <= 3; ++ii ) {
			minimizer.apply(query_pose);
			query_pose.dump_pdb( "after_min_" + string_of(ii) + ".pdb" );;
		}

		(*scorefxn)(query_pose);

		// idealize structure. maybe do this in a loop to make sure
		// that all residues are idealized? Or, maybe only idealize
		// around chainbreaks?
		//protocols::idealize::IdealizeMover idealizer;
		//idealizer.apply(query_pose);
		//query_pose.dump_pdb("after_idealization.pdb");
	}
};


int
main( int argc, char * argv [] ) {
	try {

	// options, random initialization
	devel::init( argc, argv );
	core::import_pose::pose_stream::MetaPoseInputStream input
		= core::import_pose::pose_stream::streams_from_cmd_line();

	using std::map;
	using std::string;
	using core::pose::Pose;
	using utility::vector1;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::io::silent;
	using namespace protocols::comparative_modeling;
	ResidueTypeSetCAP rsd_set = rsd_set_from_cmd_line();

	vector1< std::string > aln_fns = option[ in::file::alignment ]();
	vector1< SequenceAlignment > alns = core::sequence::read_aln(
		option[ cm::aln_format ](), aln_fns.front()
	);
	std::map< string, Pose > poses = poses_from_cmd_line(
		option[ in::file::template_pdb ]()
	);

	core::io::silent::SilentFileData sfd;

	while( input.has_another_pose() ) {
		core::pose::Pose pose;
		input.fill_pose( pose, *rsd_set );
		std::cout << "sequence = " << pose.sequence() << std::endl;

		for ( vector1< SequenceAlignment >::iterator aln_it = alns.begin(),
				end = alns.end();
				aln_it != end; ++aln_it
		) {
			string const template_id( aln_it->sequence(2)->id().substr(0,5) );
			map< string, Pose >::iterator pose_it = poses.find( template_id );
			if ( pose_it == poses.end() ) {
				string msg( "Error: can't find pose (id = "
					+ template_id + ")"
				);
				//utility_exit_with_message(msg);
				std::cerr << msg << std::endl;
			} else {
				Pose template_pose = pose_it->second;
				FastThreadingMover mover( *aln_it, template_pose );
				mover.apply( pose );
				SilentStructOP ss_out( SilentStructFactory::get_instance()->get_silent_struct_out() );
				ss_out->fill_struct(pose);
				sfd.write_silent_struct( *ss_out, option[ out::file::silent ]() );
			}
		}
	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
