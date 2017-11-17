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
/// @author Rhiju Das

// libRosetta headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

//RNA stuff.
#include <protocols/rna/denovo/RNA_DeNovoProtocol.hh>
#include <protocols/rna/denovo/setup/RNA_DeNovoSetup.hh>
#include <protocols/rna/denovo/options/RNA_DeNovoProtocolOptions.hh>
//#include <protocols/moves/PyMOLMover.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>

#include <core/pose/annotated_sequence.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/sequence/Sequence.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using namespace protocols::rna::denovo::options;
using namespace protocols::rna::denovo::setup;
using utility::vector1;

///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{
	using namespace basic::options;
	if ( option[ OptionKeys::stepwise::superimpose_over_all ].user() ) {
		std::cout << "The use of -superimpose_over_all is deprecated. The behavior in question now defaults to TRUE and is turned off by providing a particular residue that is part of an anchoring input domain as -alignment_anchor_res." << std::endl;
	}

	using namespace core::pose;
	using namespace protocols::rna::denovo;

	RNA_DeNovoSetupOP rna_de_novo_setup( new RNA_DeNovoSetup );
	rna_de_novo_setup->initialize_inputs_from_options( basic::options::option );
	rna_de_novo_setup->initialize_from_command_line();
	Pose & pose = *( rna_de_novo_setup->pose() );

	protocols::rna::denovo::RNA_DeNovoProtocol rna_de_novo_protocol( rna_de_novo_setup->options(),
		rna_de_novo_setup->rna_params() );
	rna_de_novo_protocol.set_native_pose( rna_de_novo_setup->native_pose() );
	rna_de_novo_protocol.set_refine_pose_list( rna_de_novo_setup->refine_pose_list() );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	// protocols::moves::AddPyMOLObserver( pose, false, 0.01);

	rna_de_novo_protocol.apply( pose );

}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_denovo_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence>  [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::native );
		option.add_relevant( in::file::input_res );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( out::overwrite );
		option.add_relevant( score::weights );
		option.add_relevant( score::set_weights );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::sequence );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::secstruct );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::secstruct_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize_rna );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rounds );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::fixed_stems );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::obligate_pair );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::obligate_pair_explicit );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::relax_rna );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::simple_relax );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::ignore_secstruct );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::lores_scorefxn );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::set_lores_weights);
		option.add_relevant( basic::options::OptionKeys::rna::denovo::cycles );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::temperature );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::jump_change_frequency );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::close_loops );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::close_loops_after_each_move );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::heat );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::staged_constraints );
		// option.add_relevant( basic::options::OptionKeys::rna::denovo::jump_library_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::params_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_lores_base_pairs_early );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure_halfway );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::filter_chain_closure_distance );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::output_filters );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::autofilter );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::no_filters );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::vall_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::use_1jj2_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rna_lores_chainbreak_weight );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::rna_lores_linear_chainbreak_weight );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allow_bulge  );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allowed_bulge_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::allow_consecutive_bulges );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::move_first_rigid_body );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::root_at_first_rigid_body );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::suppress_bp_constraint );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::output_res_num );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::offset );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::tag );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::refine_silent_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::refine_native );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::bps_moves );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::dump );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::output_lores_silent_file );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::out::binary_output );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::minimizer_use_coordinate_constraints );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::extra_minimize_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::extra_minimize_chi_res );
		option.add_relevant( basic::options::OptionKeys::rna::denovo::minimize::minimize_bps );
		option.add_relevant( full_model::cyclize );
		option.add_relevant( full_model::cutpoint_closed );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( full_model::rna::block_stack_above_res );
		option.add_relevant( full_model::rna::block_stack_below_res );
		option.add_relevant( basic::options::OptionKeys::rna::vary_geometry );
		option.add_relevant( basic::options::OptionKeys::rna::data_file );

		option.add_relevant( constraints::cst_file );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" ); // for chemical mapping.
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_SIDE_CHAIN" ); // for proteins for vdw screen

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

