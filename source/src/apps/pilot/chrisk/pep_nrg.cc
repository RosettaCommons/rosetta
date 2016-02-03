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


// libRosetta headers
 #include <protocols/frags/VallData.hh>
 #include <protocols/frags/TorsionFragment.hh>

 #include <core/scoring/dna/setup.hh>
 #include <core/scoring/dna/base_geometry.hh>
 #include <core/scoring/dna/BasePartner.hh>
 #include <core/scoring/GenBornPotential.hh>
 #include <core/scoring/LREnergyContainer.hh>
 #include <core/scoring/methods/Methods.hh>

 #include <protocols/simple_moves/BackboneMover.hh>
 #include <protocols/simple_moves/MinMover.hh>
 #include <protocols/moves/MonteCarlo.hh>
 #include <protocols/moves/Mover.hh>
 #include <protocols/moves/MoverContainer.hh>
 #include <protocols/moves/OutputMovers.hh>
 #include <protocols/rigid/RigidBodyMover.hh>
 // #include <protocols/moves/rigid_body_moves.hh>
 #include <protocols/moves/TrialMover.hh>
 #include <protocols/simple_moves/PackRotamersMover.hh>
 #include <protocols/simple_moves/RotamerTrialsMover.hh>
 #include <protocols/moves/RepeatMover.hh>

 #include <protocols/loops/ccd_closure.hh>
 #include <protocols/loops/loops_main.hh>

 #include <protocols/viewer/viewers.hh>

 #include <core/types.hh>

 #include <core/scoring/sasa.hh>

// #include <basic/prof.hh> // profiling
// #include <basic/CacheableData.hh> // profiling

 #include <core/id/SequenceMapping.hh>

 #include <core/chemical/AtomTypeSet.hh>
 #include <core/chemical/MMAtomTypeSet.hh>

 #include <core/chemical/AA.hh>
 #include <core/conformation/Residue.hh>
 #include <core/conformation/ResidueMatcher.hh>
 #include <core/pack/rotamer_set/RotamerCouplings.hh>
 #include <core/chemical/ResidueTypeSet.hh>
 #include <core/chemical/ResidueTypeSelector.hh>
#include <core/conformation/ResidueFactory.hh>
 #include <core/chemical/VariantType.hh>

 #include <core/chemical/ChemicalManager.hh>

 #include <core/scoring/etable/Etable.hh>
 #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
 #include <core/scoring/Ramachandran.hh>
 #include <core/pack/dunbrack/RotamerLibrary.hh>
 #include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
 #include <core/scoring/hbonds/HBondSet.hh>
 #include <core/scoring/hbonds/hbonds.hh>
 #include <core/scoring/etable/count_pair/CountPairFunction.hh>

 #include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/TaskOperation.hh>

 #include <core/kinematics/FoldTree.hh>
 #include <protocols/viewer/visualize.hh>
#include <core/kinematics/MoveMap.hh>
 #include <core/kinematics/util.hh>
 #include <core/id/AtomID_Map.hh>
// #include <core/id/AtomID_Map.Pose.hh>

 #include <core/mm/MMTorsionLibrary.hh>
 #include <core/mm/MMTorsionLibrary.fwd.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBPoseMap.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/options/util.hh>//option.hh>
// #include <basic/options/after_opts.hh>

 #include <basic/basic.hh>

 #include <basic/database/open.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

 #include <numeric/xyzVector.hh>
 #include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// //REMOVE LATER!
// #include <utility/io/izstream.hh>


// // C++ headers
 #include <cstdlib>
 #include <fstream>
 #include <iostream>
 #include <string>

//silly using/typedef
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/pep_spec.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

 using basic::T;
 using basic::Error;
 using basic::Warning;


using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace ObjexxFCL::format;
//using namespace protocols;

using utility::vector1;
using std::string;
using io::pdb::old_dump_pdb;

using namespace basic::options;
using namespace basic::options::OptionKeys;

using namespace pose;
using namespace chemical;
using namespace conformation;
using namespace scoring;
using namespace optimization;
using namespace kinematics;
using namespace protocols::moves;
using namespace id;
using namespace protocols::frags;
using utility::file::FileName;

void
RunPepSpec()
{

	vector1< Pose > poses;
	vector1< double > nrgs;
	vector1< std::string > pdb_filenames;
	std::string pdb_list_filename( option[ pep_spec::pdb_list ] );
	Pose pose;
	{
		std::ifstream pdb_list_data( pdb_list_filename.c_str() );
		if ( !pdb_list_data.good() ) {
			utility_exit_with_message( "Unable to open file: " + pdb_list_filename + '\n' );
		}
		std::string pdb_list_line;
		getline( pdb_list_data, pdb_list_line, '\t' );
		std::string filename( pdb_list_line );
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
	}

	Size pep_anchor_in( option[ pep_spec::pep_anchor ] );
	std::string const pep_chain_in( option[ pep_spec::pep_chain ] );
	Size pep_anchor( pose.pdb_info()->pdb2pose( pep_chain_in[0], pep_anchor_in ) );
	Size pep_chain( pose.chain( pep_anchor ) );
	Size pep_begin( pose.conformation().chain_begin( pep_chain ) );
	Size pep_end( pose.conformation().chain_end( pep_chain ) );

	Size prot_chain;
	for( Size i = 1; i <= pose.conformation().num_chains(); ++i ){
		if( i == pep_chain ) continue;
		if( !( pose.residue( pose.conformation().chain_begin( i ) ).is_protein() ) ) continue;
		else{
			prot_chain = i;
			break;
		}
	}
	Size prot_begin( pose.conformation().chain_begin( prot_chain ) );
	Size prot_end( pose.conformation().chain_end( prot_chain ) );
	Size prot_anchor( prot_begin );

	core::scoring::ScoreFunctionOP scorefxn(  ScoreFunctionFactory::create_score_function( option[ pep_spec::wts ] ) );
	ScoreTypes score_types( scorefxn->get_nonzero_weighted_scoretypes() );
	std::string pairvec_file_name_str( option[ out::file::o ]+".nrgkey" );
	char const *pairvec_file_name = pairvec_file_name_str.c_str();
	std::fstream pairvec_file( pairvec_file_name, std::ios::out );
	for( Size i_seqpos = pep_begin; i_seqpos <= pep_end; ++i_seqpos ){
		for ( ScoreTypes::const_iterator score_type( score_types.begin() ); score_type != score_types.end(); ++score_type ){
			pairvec_file << string_of( i_seqpos ) << "\t" << name_from_score_type( *score_type ) << "\n";
		}
	}

	std::string pairmat_file_name_str( option[ out::file::o ]+".nrgmat" );
	char const *pairmat_file_name = pairmat_file_name_str.c_str();
	std::fstream pairmat_file( pairmat_file_name, std::ios::out );

	std::ifstream pdb_list_data( pdb_list_filename.c_str() );
	std::string pdb_list_line;
	while( !getline( pdb_list_data, pdb_list_line, '\t' ).eof() ) {
		std::string filename( pdb_list_line );
		pose::Pose pose;
		core::import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
		( *scorefxn )( pose );
		std::string data;
		getline( pdb_list_data, data );

		for( Size i_seqpos = pep_begin; i_seqpos <= pep_end; ++i_seqpos ){
			for ( ScoreTypes::const_iterator score_type( score_types.begin() ); score_type != score_types.end(); ++score_type ){
				Real nrg( pose.energies().residue_total_energies( i_seqpos )[ *score_type ] * scorefxn->weights()[ *score_type ] );
				pairmat_file << string_of( nrg ) << "\t";
			}
		}
		pairmat_file << "\n";
	}
	pdb_list_data.close();
}

int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init(argc, argv);
		RunPepSpec();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
