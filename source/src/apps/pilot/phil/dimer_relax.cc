// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/** @details  Protein-DNA relax protocol for protein dimers.

		Current Features:
		  ** No internal protein flexibility
			** No DNA backbone flexibility
			** DNA sidechain sampling
			** protein sidechain sampling
			** rigid-body moves of both monomers via jumps from DNA (system rooted in DNA)

 **/


// libRosetta headers
#include <devel/dna/ProteinDNA_Relax.hh>

#include <protocols/viewer/viewers.hh>

#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/excn/Exceptions.hh>


#include <core/conformation/Residue.hh>


#include <core/kinematics/MoveMap.hh>


#include <core/pose/Pose.hh>

#include <basic/options/util.hh>
#include <basic/options/keys/dna.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>


#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>


#include <numeric/random/random.hh>

#include <ObjexxFCL/string.functions.hh>


// // C++ headers

//silly using/typedef


#include <basic/Tracer.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/scoring/constraints/Constraint.hh>


// option key includes


using basic::T;
using basic::Error;
using basic::Warning;


using namespace core;
//using namespace protocols;

using utility::vector1;
using std::string;


static THREAD_LOCAL basic::Tracer tt( "demo.phil.dimer_relax", basic::t_trace );
static THREAD_LOCAL basic::Tracer td( "demo.phil.dimer_relax", basic::t_debug );
static THREAD_LOCAL basic::Tracer ti( "demo.phil.dimer_relax", basic::t_info );
static THREAD_LOCAL basic::Tracer tw( "demo.phil.dimer_relax", basic::t_warning );


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// danger USING /////////////////////////////////////////////////////////////////////////////////////////
using namespace core;
using namespace protocols;
using namespace pose;
using namespace conformation;
using namespace chemical;
using namespace scoring;
using namespace basic::options;
using namespace optimization;
using namespace id;
using core::import_pose::pose_from_file;
namespace OK = OptionKeys;
using utility::vector1;
using std::string;
/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////


/**
void
setup_dimer_relax_pose( Pose & pose )
{
	Pose const start_pose( pose );

	// foldtree: dont currently need dna-flexibility foldtree
	//
	// choose jump points for DNA->protein jumps


	// requirements: >=4 chains total, 1st 2 protein and 1st 2 DNA taken as simulation chains
	//
	Size const nchains( pose.conformation().num_chains() );
	assert( nchains >= 4 );
	Size pchain1(0), pchain2(0), dchain1(0), dchain2(0);
	for ( Size i=1; i<= nchains; ++i ) {
		Residue const & rsd( pose.residue( pose.conformation().chain_begin(i) ) );
		if ( rsd.is_protein() ) {
			if ( !pchain1 ) pchain1 = i;
			else if ( !pchain2 ) pchain2 = i;
		} else if ( rsd.is_DNA() ) {
			if ( !dchain1 ) dchain1 = i;
			else if ( !dchain2 ) dchain2 = i;
		}
	}


	// make a new pose
	pose.clear();
	vector1< Size > chains;
	chains.push_back( dchain1 ); chains.push_back( dchain2 );
	chains.push_back( pchain1 ); chains.push_back( pchain2 );
	append_chains_from_pose( pose, start_pose, chains ); // in pose/util.cc

	// look for anchor points to the protein chains
	find_protein_DNA_anchor_point( setup_chain_mask( pose, 3 ), prot_root1, dna_anchor1 );
	find_protein_DNA_anchor_point( setup_chain_mask( pose, 4 ), prot_root2, dna_anchor2 );

	FoldTree f( pose.size() );
	f.new_jump( dna_anchor1, prot_root1, pose.conformation().chain_end(2) );
	f.new_jump( dna_anchor2, prot_root2, pose.conformation().chain_end(3) );
	f.new_jump( dna_anchor1, dna_anchor2, pose.conformation().chain_end(1) );
	f.reorder( dna_anchor1 );
	pose.fold_tree(f);

}

**/

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* )
{

	Pose pose;
	core::import_pose::pose_from_file( pose, start_file() , core::import_pose::PDB_file);

	ScoreFunction scorefxn;

	//scorefxn.energy_method_options().exclude_DNA_DNA( false );
	// Safer:
	methods::EnergyMethodOptions options( scorefxn.energy_method_options() );
	options.exclude_DNA_DNA( false );
	scorefxn.set_energy_method_options( options );

	scorefxn.add_weights_from_file( option[ OK::dna::specificity::score_function ] );

	// hard-coded!
	vector1< Size > moving_jumps;
	moving_jumps.push_back(1);
	moving_jumps.push_back(2);

	{ // sanity check
		for ( Size i=1; i<= moving_jumps.size(); ++i ) {
			if ( ! ( pose.residue( pose.fold_tree().  upstream_jump_residue( moving_jumps[i] ) ).is_DNA() &&
							 pose.residue( pose.fold_tree().downstream_jump_residue( moving_jumps[i] ) ).is_protein() ) ) {
				std::cout << pose.fold_tree() << std::endl;
				std::cout << "moving jump: " << i << std::endl;
				utility_exit_with_message( "bad foldtree for dimer_relax!" );
			}
		}
	}

	Pose const start_pose( pose );

	// for graphics:
	protocols::viewer::add_conformation_viewer( pose.conformation(), "dimer_relax" );

	for ( int n=1; n<= option[ OK::out::nstruct ]; ++n ) {
		pose = start_pose;
		devel::dna::ProteinDNA_Relax relax( scorefxn, moving_jumps );
		relax.rot_mag( 1.5 ).trans_mag( 0.25 ); // half the default values
		relax.apply( pose );
		std::string const filename
			( option[ OK::out::output_tag ] + "dimer_relax" + lead_zero_string_of( n, 4 ) + ".pdb" );
		Real const final_score( scorefxn( pose ) );
		std::cout << "final_score: " << filename << ' ' << final_score << ' ' <<
			scoring::all_atom_rmsd( pose, start_pose ) << ' ' <<
			pose.energies().total_energies().weighted_string_of( scorefxn.weights() ) << std::endl;

		pose.dump_scored_pdb( filename, scorefxn );
	}

// 	Pose pose;
// 	core::import_pose::pose_from_file( pose, start_file() , core::import_pose::PDB_file);

// 	pose.dump_pdb("test.pdb");

// 	for ( Size i=1; i<= pose.conformation().num_chains(); ++i ) {
// 		std::cout << "chain bounds: " << i << ' ' <<
// 			pose.conformation().chain_begin(i) << ' ' << pose.conformation().chain_end(i) << std::endl;
// 	}

// 	exit(0);
	return 0;
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try{
// 	{ // add any new options
// 		using namespace utility::options;
// 		BooleanOptionKey const myopt = BooleanOptionKey("phil:dof_constraint_weight");
// 		option.add(myopt, "Does nothing, nothing at all");
// 	}

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
