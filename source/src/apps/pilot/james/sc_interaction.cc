// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson <tex@u.washington.edu>

// libRosetta headers


#include <core/types.hh>

#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/rms_util.hh>

#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>

#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/util.hh>

#include <core/id/AtomID.hh>

#include <core/sequence/Sequence.hh>
#include <core/id/SequenceMapping.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/PDBSilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <basic/basic.hh>
#include <basic/Tracer.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/loop_build/LoopBuild.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.hh>
#include <protocols/abinitio/MaxSeqSepConstraintSet.fwd.hh>
#include <protocols/comparative_modeling/ConstraintRemodelMover.hh>

#include <protocols/jobdist/standard_mains.hh>
#include <protocols/jobdist/not_universal_main.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>


// option key includes

#include <basic/options/keys/james.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

class InteractionDistMinimizer : public protocols::moves::Mover {
public:
	InteractionDistMinimizer( core::Real const dist, core::Real const sdev )
		: dist_( dist ), sdev_( sdev )
	{}

	// 1. For each pair of residues, define a constraint that
	// attempts to minimize the ends of the side chains to
	// within 2.0A
	// 2. Report the minimum distance between the side chain
	// atoms of the two residues as a comment inside the pose.
	void apply( core::pose::Pose & pose ) {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::scoring;
		using namespace core::scoring::constraints;
		using namespace core::chemical;
		using namespace core::optimization;
		using core::id::AtomID;
		using ObjexxFCL::format::F;

		ScoreFunctionOP scorefxn( get_score_function() );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );

		core::kinematics::MoveMap mm;
		mm.set_chi( true );
		std::string const min_type("dfpmin_armijo_nonmonotone");

		for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		for ( Size jj = ii+1; jj <= pose.total_residue(); ++jj ) {
			std::cout << ii << "," << jj << std::endl;
			core::pose::Pose pose_copy = pose;
			ResidueType const & res_type_ii( pose_copy.residue_type(ii) );
			ResidueType const & res_type_jj( pose_copy.residue_type(jj) );

			Size const nbr_atom_ii( res_type_ii.nbr_atom() );
			Size const nbr_atom_jj( res_type_jj.nbr_atom() );

			AtomID id1( nbr_atom_ii, ii );
			AtomID id2( nbr_atom_jj, jj );

			ConstraintSetOP cstset( new core::scoring::constraints::ConstraintSet );
			FuncOP func( new HarmonicFunc( dist_, sdev_ ) );
			ConstraintOP cst( new AtomPairConstraint( id1, id2, func ) );
			cstset->clear();
			cstset->add_constraint( cst );

			pose_copy.constraint_set( cstset );
      AtomTreeMinimizer().run(
				 pose_copy, mm, (*scorefxn), MinimizerOptions(min_type,0.001,true)
      );

			core::Real const dist( pose_copy.residue(ii).xyz(nbr_atom_ii).distance(
				pose_copy.residue(jj).xyz(nbr_atom_jj)
			) );
			core::Real const dist_cutoff(
				res_type_ii.nbr_radius() + res_type_jj.nbr_radius()
			);
			if ( dist <= dist_cutoff ) {
				std::string key(
					"INTERACTION_DIST_" + string_of(ii) + "_" + string_of(jj)
				);
				add_comment( pose, key, F( 8, 3, dist ) );
			}
		} // ii
		} // jj
	} // apply

private:
	core::Real const dist_;
	core::Real const sdev_;
}; // InteractionDistMinimizer

void*
my_main( void* ) {
	using core::Real;
	Real const sdev( 0.5 );
	Real const min_dist( 2.0 ); // try to get atoms to interact within 2.0A
	InteractionDistMinimizer my_mover( min_dist, sdev );
	protocols::jobdist::not_universal_main( my_mover );
	return 0;
}

int
main( int argc, char* argv [] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );
	protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // int main
