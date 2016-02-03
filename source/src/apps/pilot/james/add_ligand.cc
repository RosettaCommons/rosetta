// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief
/// @author James Thompson

#include <core/types.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <basic/options/util.hh>
#include <basic/options/after_opts.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/rms_util.hh>
#include <core/pose/util.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/relax_protocols.hh>
#include <protocols/frags/TorsionFragment.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <core/io/silent/BinarySilentStruct.hh>

#include <protocols/viewer/visualize.hh>

#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <ctime>

//silly using/typedef
//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

using basic::T;
using basic::Error;
using basic::Warning;

int
main( int argc, char* argv[] ) {
	try {

	basic::Tracer tr( "james.add_calcium" );
	// options, random initialization
	devel::init( argc, argv );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::Size;
	using utility::vector1;

	// add the ligand:

	core::pose::Pose query_pose, ligand_pose;

	std::string query_pose_file  = option[ in::file::s ]()[1];
	std::string ligand_pose_file = option[ in::file::template_pdb ]()[1];

	core::import_pose::pose_from_file( query_pose , query_pose_file , core::import_pose::PDB_file);
	core::import_pose::pose_from_file( ligand_pose, ligand_pose_file , core::import_pose::PDB_file);

	tr.Error << ligand_pose.fold_tree() << std::endl;

	// configuration stuff
	core::Size  const query_pose_anchor ( 10 ); // GLU10 in query
	core::Size  const ligand_pose_anchor( 12 ); // is like GLU12 in template
	std::string const anchor_atom_name("CA");
	std::string const ligand_atom_name("CA");

	vector1< core::Size > ligand_indices;
	ligand_indices.push_back( 320 );
	ligand_indices.push_back( 321 );
	ligand_indices.push_back( 322 );

	core::kinematics::FoldTree new_fold_tree;

	// peptide edges
	new_fold_tree.add_edge(
		1,
		ligand_pose_anchor,
		core::kinematics::Edge::PEPTIDE
	);
	new_fold_tree.add_edge(
		ligand_pose_anchor,
		ligand_pose.total_residue() - 3,
		core::kinematics::Edge::PEPTIDE
	);
	tr.Error << "adding ligand residues to fold-tree" << std::endl;

	// edges from anchor to ligand residues
	for ( Size jj = 1; jj <= ligand_indices.size(); ++jj ) {
		tr.Error << "adding " << jj << std::endl;
		core::kinematics::Edge out_edge(
			ligand_pose_anchor, // start
			ligand_indices[jj],
			(int)jj, // label
			std::string(anchor_atom_name), // start_atom
			std::string(ligand_atom_name), // stop_atom
			false // bKeepStubInResidue
		);
		tr.Error << out_edge << std::endl;
		new_fold_tree.add_edge( out_edge );
	}

	tr.Error << ligand_pose.fold_tree();
	ligand_pose.fold_tree( new_fold_tree );

	for ( Size jj = 1; jj <= ligand_indices.size(); ++jj ) {
		query_pose.append_residue_by_jump(
			ligand_pose.residue( ligand_indices[jj] ),
			query_pose_anchor,
			anchor_atom_name,
			ligand_atom_name
		);
		query_pose.set_jump( jj, ligand_pose.jump(jj) );
	}

	// remodel loops

	query_pose.dump_pdb( "pose_with_ligand.pdb" );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;


} // main
