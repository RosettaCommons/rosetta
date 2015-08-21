// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief

// libRosetta headers

#include <numeric/xyzVector.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>

#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/simple_moves/SaneMinMover.hh>

#include <core/scoring/rms_util.hh>

#include <devel/init.hh>
#include <iostream>
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/string.functions.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <utility/excn/Exceptions.hh>


class SetupMover : public protocols::moves::Mover {
public:

	SetupMover() : Mover("SetupMover") {}

	virtual std::string get_name() const {
		return "SetupMover";
	}

	void apply( core::pose::Pose & pose ) {
		using core::Size;
		using core::Real;
		using std::string;
		using core::id::AtomID;
		using numeric::xyzVector;

		// select residues 18 - 41
		// translate them ~20A away
		Size const start(18);
		Size const stop(41);
		//Real const translation(20);  // unused ~Labonte
		for ( Size ii = start; ii <= stop; ++ii ) {
			for ( Size jj = 1; jj <= pose.residue(ii).natoms(); ++jj ) {
				AtomID id(jj,ii);
				xyzVector< Real > coord = pose.xyz(id);
				coord = coord + 20;
				pose.set_xyz( id, coord );
			}
		}

		// set up a new fold-tree:
		using core::kinematics::FoldTree;
		FoldTree ft;

		ObjexxFCL::FArray2D_int ft_jumps(2,2);
		ObjexxFCL::FArray1D_int ft_cuts(2);
		ft_jumps(1,1) = 1;
		ft_jumps(2,1) = 30;
		ft_jumps(1,2) = 1;
		ft_jumps(2,2) = 50;
		ft_cuts(1)    = start-1;
		ft_cuts(2)    = stop;

		bool status = ft.tree_from_jumps_and_cuts(
			pose.total_residue(), // nres_in
			2,                    // num_jump_in
			ft_jumps,             // jump_point
			ft_cuts,              // cuts
			1                     // root
		);
		if ( !status ) {
			utility_exit_with_message("failed to build fold tree from cuts and jumps");
		}
		std::cout << "ft = " << ft << std::endl;
		pose.fold_tree(ft);

		// add cutpoint variants to 17,18 and 41,42
		const core::kinematics::FoldTree& tree(pose.fold_tree());
		for ( Size ii = 1; ii <= pose.total_residue() - 1; ++ii ) {
			using core::pose::add_variant_type_to_pose_residue;
			//std::cout << "considering " << ii << std::endl;
			if ( tree.is_cutpoint(ii) ) {
				add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_LOWER, ii);
				add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_UPPER, ii+1);
				std::cout << "adding variant to " << ii << "," << ii+1 << std::endl;
			}
		}
	} // apply
};

void dump_pose(
	std::string const & tag, core::pose::Pose & pose, core::pose::Pose const & orig_pose,
	core::scoring::ScoreFunctionOP scorefxn
) {
	std::string const fn_out( tag + ".pdb" );
	pose.dump_pdb(fn_out);
	core::Real score = (*scorefxn)(pose);
	core::Real rmsd = core::scoring::CA_rmsd( orig_pose, pose );
	core::Real gdtmm = core::scoring::CA_gdtmm( orig_pose, pose );
	std::cout << " " << score << " " << rmsd << " " << gdtmm << std::endl;
}

int
main( int argc, char * argv [] ) {
	try {

		using namespace protocols;
		using namespace protocols::jobdist;
		using namespace protocols::moves;
		using namespace basic::options;
		using namespace core::scoring::constraints;
		using namespace basic::options::OptionKeys;

		using core::pose::Pose;
		using core::Real;
		using core::Size;

		devel::init(argc, argv);

		core::pose::Pose pose;
		core::import_pose::pose_from_pdb(pose, option[ in::file::s ]()[1]);
		core::pose::Pose const orig_pose(pose);

		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() );
		//scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
		scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
		//scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );
		//scorefxn->set_weight( core::scoring::distance_chainbreak, 1.0 );

		core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
		mm->set_bb (false);
		mm->set_chi(false);
		mm->set_jump(1,true);
		mm->set_jump(2,false);

		SetupMover setup;
		pose.dump_pdb("pre_setup.pdb");
		scorefxn->show(pose);
		setup.apply(pose);
		pose.dump_pdb("post_setup.pdb");
		scorefxn->show(pose);
		core::optimization::MinimizerOptionsOP min_options( new core::optimization::MinimizerOptions( "dfpmin", 1e-2, true ) );

		//runtime_assert( option[ run::nblist_autoupdate ]() );

		if ( !option[ james::debug ]() ) {
			// case 1: 100 rounds of 20 iterations each.
			min_options->max_iter(20);
			for ( Size ii = 1; ii <= 100; ++ii ) {
				protocols::simple_moves::SaneMinMover min_mover( mm, scorefxn, min_options );
				min_mover.apply(pose);

				std::string tag( "case_1." + ObjexxFCL::string_of(ii) );
				dump_pose(tag,pose,orig_pose,scorefxn);
			}
		} else {
			// case 2: 1 round of 2000 iterations
			min_options->max_iter(2000);
			protocols::simple_moves::SaneMinMover min_mover( mm, scorefxn, min_options );
			std::string tag( "case_2.final" );
			min_mover.apply(pose);
			dump_pose(tag,pose,orig_pose,scorefxn);
		}
		scorefxn->show(pose);

		protocols::simple_moves::SaneMinMover min_mover( mm, scorefxn, min_options );
		min_mover.apply(pose);
		std::string const fn_out( "post_min_final.pdb" );
		pose.dump_pdb(fn_out);

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
