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
#include <devel/dna/protocols.hh>
#include <devel/dna/util.hh>
//#include <devel/dna/util.hh>
#include <protocols/frags/VallData.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/frags/TorsionFragment.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/rigid_body_moves.hh>

#include <core/scoring/LREnergyContainer.hh>
#include <core/scoring/dna/setup.hh>
#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/etable/Etable.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/elec/FA_ElecEnergy.hh>
//#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/count_pair/CountPairAll.hh>
#include <core/scoring/etable/count_pair/CountPairFunction.hh>
#include <core/scoring/etable/count_pair/CountPairFactory.hh>
//#include <core/scoring/etable/count_pair/CountPair1BC4.hh>

#include <core/types.hh>

//#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/chemical/VariantType.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/util.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_set/RotamerCouplings.hh>
#include <core/pack/rotamer_set/WaterPackingInfo.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/visualize.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/id/AtomID_Map.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>

#include <basic/prof.hh> // profiling
#include <basic/basic.hh>
#include <core/id/SequenceMapping.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


#include <basic/Tracer.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <set>
#include <cstdlib>
#include <sstream>
#include <math.h>

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>





////////////////////////////////////////////////
// danger USING ////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace id;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;

/// stochastic -- choose cutpoints, may update pose foldtree
/// @note  mapping goes from the source pdb to the target sequence, and has already been extended and applied, ie
/// the sequence of pose is the target sequence and pose.total_residue() == mapping.size2()
/// @note  CURRENTLY SKIPS TERMINAL LOOPS FOR HISTORICAL REASONS!!!
void
setup_loops_from_mapping(
	pose::Pose & pose,
	id::SequenceMapping const & mapping,
	protocols::loops::Loops & loops,
	bool const start_extended
)
{
	//assert( mapping.size2() == pose.total_residue() ); // TOO STRONG
	assert( mapping.size2() <= pose.total_residue() ); // in case pose has extra stuff at the end

	id::SequenceMapping m( mapping );
	m.reverse();


	bool in_loop( false );
	Size loop_begin( 0 );

	for ( Size i=1; i<= m.size1(); ++i ) {
		if ( m[i] == 0 ) {
			if ( !in_loop ) {
				loop_begin = i;
				in_loop = true;
			}
		} else {
			if ( in_loop ) {
				Size const loop_end( i-1 );
				Size const loop_size( loop_end - loop_begin + 1 );

				// detect and skip a terminal loop:
				if ( !pose.residue( loop_begin ).is_protein() || !pose.residue( loop_end ).is_protein() ) {
					// a non-protein "loop
					utility_exit_with_message( "nonprotein loop?" );
				} else if ( loop_begin == 1 || loop_end == pose.total_residue() ||
						 !pose.residue( loop_begin-1 ).is_protein() ||
						 !pose.residue( loop_end+1   ).is_protein() ||
						 pose.chain( loop_begin-1 ) != pose.chain( loop_end+1 ) ) {
					// need to add logic for terminal loops
					//std::cout << "add logic for terminal loops!" << std::endl;
					loops.add_loop( loop_begin, loop_end, loop_end, 0.0, 0, start_extended );
				} else {
					Size const cutpoint = loop_begin-1 + static_cast< Size >( numeric::random::uniform() * (loop_size+1) );
					loops.add_loop( loop_begin, loop_end, cutpoint, 0 /*offset*/, start_extended );
					// remodel the foldtree to reflect this new cutpoint
					protocols::loops::set_loop_cutpoint_in_pose_fold_tree( cutpoint, pose, loop_begin, loop_end );
				}
				in_loop = false;
			}
		}
	}
}


typedef utility::vector1< std::pair< std::string, Size > > SS_Quotas;
static SS_Quotas const empty_ss_quotas;

void
setup_loops_fragment_libraries_with_ss(
	std::string const target_seq,
	protocols::loops::Loops const & loops,
	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > & frag_libs,
	std::map< Size, bool > & frag_libs_init,
	SS_Quotas const & ss_quotas,
	std::map< Size, Real > const & ss_weight
																			 )
{
	using namespace basic::options;


	Real const seq_wt( 1.0 );

	// total size of sequence
	int const nres( target_seq.size() );

	// read in the whole vall, only once
	static protocols::frags::VallData vall( option[ OptionKeys::loops::vall_file ]().name() );

	// for each frag_size library
	for ( std::map<Size, bool>::iterator it = frag_libs_init.begin(),
					it_end = frag_libs_init.end(); it != it_end; ++it ) {
		// skip if it has been initialized
		if ( it->second ) continue;
		it->second = true;
		int const frag_size( it->first );
		protocols::frags::TorsionFragmentLibrary & lib( *(frag_libs.find(frag_size)->second) );
		lib.resize( nres - frag_size + 1 );
		assert( ss_weight.count( frag_size ) );
		Real const ss_wt( ss_weight.find( frag_size )->second );

		// for each ss
		for ( SS_Quotas::const_iterator ssq= ss_quotas.begin(); ssq != ss_quotas.end(); ++ssq ) {
			std::string const & secstruct( ssq->first );
			Size const nfrags( ssq->second );

			// for each loop
			for ( protocols::loops::Loops::const_iterator it2 = loops.begin(), it2_end = loops.end(); it2 != it2_end; ++it2 ) {
				int const loop_begin( it2->start() ), loop_end( it2->stop() );
				assert( loop_end <= secstruct.size() );

				// for each residue window
				for ( int i = loop_begin; i <= loop_end-frag_size+1; ++i ) {
					std::string const frag_seq( target_seq.substr( i-1, frag_size ) );
					std::string const frag_ss (  secstruct.substr( i-1, frag_size ) );
					vall.get_frags( nfrags, frag_seq, frag_ss, seq_wt, ss_wt, false, false, true, lib[i] );
				} // for each residue window
			}// for each loop
		} // for each ss
	}// for each fragment library/size
}

/**
	 fragment setup:
	 1. pick 3mer library from vall
	 2. derive 1mer libary from 3mers
**/


void
setup_frags_from_vall(
	pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	//utility::vector1<int> const & frag_sizes,
	std::map< Size, protocols::frags::TorsionFragmentLibraryOP > & frag_libs,
	std::map< Size, bool > & frag_libs_init,
	SS_Quotas const & ss_quotas = empty_ss_quotas,
	bool const pick_6mers = false
)
{
	frag_libs.clear();
	frag_libs_init.clear();

	std::map< Size, Real > ss_weight;
	ss_weight[ 6 ] = 2.0;
	ss_weight[ 3 ] = 0.5;


	using namespace protocols::frags;
	if ( pick_6mers ) { // 6mers
		Size const frag_size( 6 );
		protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
		frag_libs_init.insert( std::make_pair( frag_size , false ) );
		frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
	}
	{ // 3mers
		Size const frag_size( 3 );
		protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
		frag_libs_init.insert( std::make_pair( frag_size , false ) );
		frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );
	}

	{ // now pick from vall
		if ( ss_quotas.empty() ) {
			protocols::loops::setup_loops_fragment_libraries( pose.sequence(), loops, frag_libs, frag_libs_init );
		} else {
			setup_loops_fragment_libraries_with_ss( pose.sequence(), loops, frag_libs, frag_libs_init, ss_quotas, ss_weight );
		}
	}

	{ // 1mers
		// gets 1mers from 3mers
		Size const frag_size = 1;
		Size const prev_size = 3;
		protocols::frags::TorsionFragmentLibraryOP frag_lib_op( new protocols::frags::TorsionFragmentLibrary );
		frag_libs_init.insert( std::make_pair( frag_size , false ) );
		frag_libs.insert( std::make_pair(frag_size, frag_lib_op) );

		protocols::frags::TorsionFragmentLibraryOP    prev_lib_op( frag_libs.find( prev_size)->second );

		frag_libs_init[ frag_size ] = frag_lib_op->derive_from_src_lib( frag_size, prev_size, prev_lib_op );
	}

	// confirm success
	assert( frag_libs_init[1] && frag_libs_init[3] );
}

void
full_protein_repack( pose::Pose & pose, ScoreFunction const & scorefxn )
{
	Size const nloop( 20 );

	pack::task::PackerTaskOP task( pack::task::TaskFactory::create_packer_task( pose ));
	task->initialize_from_command_line().or_include_current( true );

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( pose.residue(ii).is_protein() ) {
			task->nonconst_residue_task( ii ).restrict_to_repacking();
		} else {
			task->nonconst_residue_task( ii ).prevent_repacking();
		}
	}

	//Real const score_before_packing( scorefxn( pose ) );

	pack::pack_rotamers_loop( pose, scorefxn, task, nloop );

	//TR.Trace << "packing scores: " << score_before_packing << ' '<< scorefxn(pose) << std::endl;
}


