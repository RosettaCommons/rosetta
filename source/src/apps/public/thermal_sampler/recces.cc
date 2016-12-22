// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Legacy version of recces app, to be merged with thermal_sampler by December 2016

// libRosetta headers
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/thermal_sampling/RECCES_Mover.hh>
#include <protocols/stepwise/modeler/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/viewer/viewers.hh>

#include <utility/io/ozstream.hh>

// C++ headers
#include <string>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Exception handling
#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace core::pose;
using namespace protocols;
using namespace protocols::stepwise;
using namespace protocols::thermal_sampling;
using namespace basic::options::OptionKeys;
using namespace basic::options;

OPT_KEY( String, seq1 )
OPT_KEY( String, seq2 )
OPT_KEY( Integer, n_cycle )
OPT_KEY( Real, a_form_range )
OPT_KEY( RealVector, temps )
OPT_KEY( RealVector, st_weights )
OPT_KEY( String, out_prefix )
OPT_KEY( Boolean, save_score_terms )
OPT_KEY( Boolean, dump_pdb )
OPT_KEY( Integer, n_intermediate_dump )

//////////////////////////////////////////////////////////////////////////////
PoseOP pose_setup(
	std::string const & seq1,
	std::string const & seq2,
	Size const len1
) {
	protocols::stepwise::modeler::rna::helix::RNA_HelixAssembler assembler;
	assembler.use_phenix_geo( true );
	PoseOP pose( assembler.build_init_pose( seq1, seq2 ) );
	add_variant_type_to_pose_residue( *pose, chemical::VIRTUAL_PHOSPHATE, 1 );
	if ( seq1 != "" && seq2 != "" ) {
		add_variant_type_to_pose_residue(
			*pose, chemical::VIRTUAL_PHOSPHATE, len1 + 1 );
	}
	return pose;
}

//////////////////////////////////////////////////////////////////////////////
void
MC_run() {
	using namespace protocols::stepwise::sampler::rna;
	using namespace protocols::moves;
	using namespace scoring;

	RECCES_Mover recces_mover( option[ temps ](),
														 option[ st_weights ]() );
	if ( option[ score::weights ].user() ) recces_mover.set_scorefxn( get_score_function() );
	recces_mover.set_n_cycle( option[n_cycle]() );
	recces_mover.set_dump_pdb( option[ dump_pdb ]() );
	recces_mover.set_n_dump( option[ n_intermediate_dump ]() );
	recces_mover.set_save_scores( option[ save_score_terms ]() );
	recces_mover.set_a_form_range( option[ a_form_range ]() );
	recces_mover.set_out_prefix( option[ out_prefix ]() );

	/////////////////////////////////////////////////
	// Figure out bp and dangling residues, needed
	//  to initialize sampler inside recces_mover.
	//  Later, should be able to figure this out
	//  *inside* recces_mover based on pose sequence
	//  and fold-tree.
	std::string const & seq1_( option[seq1]() );
	std::string const & seq2_( option[seq2]() );
	Size const len1( get_sequence_len( seq1_ ) );
	Size const len2( get_sequence_len( seq2_ ) );
	utility::vector1< Size > bp_rsd, dangling_rsd;
	Size const n_bp( std::min( len1, len2 ) );
	Size const total_len( len1 + len2 );
	for ( Size i = 1; i <= total_len; ++i ) {
		if ( i > n_bp && i <= total_len - n_bp ) {
			dangling_rsd.push_back( i );
		} else {
			bp_rsd.push_back( i );
		}
	}
	recces_mover.set_bp_rsd( bp_rsd );
	recces_mover.set_dangling_rsd( dangling_rsd );
	/////////////////////////////////////////////////

	// Later replace this with FullModelInfoSetupFromCommandLine:
	Pose pose( *pose_setup( option[seq1](), option[seq2](), get_sequence_len( option[seq1]() ) ) );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	recces_mover.apply( pose );

}

//////////////////////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	MC_run();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
//////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	using namespace core;
	utility::vector1< Real > null_real_vector;
	NEW_OPT( seq1, "sequence 1 to model, 3' to 5' ", "" );
	NEW_OPT( seq2, "sequence 2 to model, 3' to 5' ", "" );
	NEW_OPT( n_cycle, "cycle number for random sampling", 0 );
	NEW_OPT( a_form_range, "Range of sampling near A-form for duplexes.", 60.0 );
	NEW_OPT( temps, "Simulated tempering temperatures", null_real_vector );
	NEW_OPT( st_weights, "Simulated tempering weights", null_real_vector );
	NEW_OPT( out_prefix, "prefix for the out file", "recces" );
	NEW_OPT( save_score_terms,
		"Save scores and individual score terms"
		" of all sampled conformers", false );
	NEW_OPT( dump_pdb, "Dump pdb files", false );
	NEW_OPT( n_intermediate_dump,
		"Number of intermediate conformations to be dumped", 0 );

	try {
		devel::init ( argc, argv );
		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

