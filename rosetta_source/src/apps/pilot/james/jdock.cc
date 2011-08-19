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
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/scoring/constraints/util.hh>
#include <core/import_pose/import_pose.hh>

#include <protocols/jobdist/not_universal_main.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/MinMover.hh>
#include <protocols/moves/CompositionMover.hh>
#include <protocols/moves/ConstraintSetMover.hh>

#include <devel/init.hh>
#include <iostream>
#include <string>

#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option.hh>

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
		Real const translation(20);
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
		std::cout << "ft = " << ft << std::endl;
 		if (!status)
			utility_exit_with_message("failed to build fold tree from cuts and jumps");
		std::cout << "ft = " << ft << std::endl;
		pose.fold_tree(ft);

		// add cutpoint variants to 17,18 and 41,42
		const core::kinematics::FoldTree& tree(pose.fold_tree());
		for ( Size ii = 1; ii <= pose.total_residue() - 1; ++ii ) {
			using core::pose::add_variant_type_to_pose_residue;
			std::cout << "considering " << ii << std::endl;
			if ( tree.is_cutpoint(ii) ) {
				add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_LOWER, ii);
				add_variant_type_to_pose_residue(pose, core::chemical::CUTPOINT_UPPER, ii+1);
				std::cout << "adding variant to " << ii << "," << ii+1 << std::endl;
			}
		}
	} // apply
};

int
main( int argc, char * argv [] ) {
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace basic::options;
	using namespace core::scoring::constraints;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

	core::pose::Pose pose;
	core::import_pose::pose_from_pdb(pose, option[ in::file::s ]()[1]);

	//pose.constraint_set( new ConstraintSet );
	//using core::Size;
	// add CA-CA constraints from i => i+1
	//for ( Size ii = 1; ii <= pose.total_residue() - 1; ++ii ) {
	//	pose.constraint_set().add_constraint(
	//		new AtomPairConstraint( AtomID(2,ii), AtomID(2,ii+1) ), HarmonicFunc( 3.8, 2 )
	//	);
	//}

	core::scoring::ScoreFunctionOP scorefxn( core::scoring::getScoreFunction() );
	scorefxn->set_weight( core::scoring::chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::linear_chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::overlap_chainbreak, 1.0 );
	scorefxn->set_weight( core::scoring::distance_chainbreak, 1.0 );

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

	pose.dump_pdb("pre_min.pdb");
	scorefxn->show(pose);
	for ( core::Size ii = 1; ii <= 100; ++ii ) {
		std::cout << "iteration " << ii << std::endl;
		protocols::moves::MinMover min_mover( mm, scorefxn, "dfpmin", 1e-5, false );
		min_mover.apply(pose);
		pose.dump_pdb("post_min.pdb");
		scorefxn->show(pose);
	}
}
