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
#include <core/init/init.hh>
#include <utility/vector1.hh>

//RNA stuff.
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/setup/RNA_DeNovoSetup.hh>
#include <protocols/farna/options/RNA_DeNovoProtocolOptions.hh>
//#include <protocols/moves/PyMolMover.hh>

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
#include <core/sequence/Sequence.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

vector1<pose::PoseOP> get_refine_pose_list( std::string const & input_silent_file,
	std::pair< utility::vector1< int >, utility::vector1< char > > const & output_res_and_chain,
	core::chemical::ResidueTypeSetCOP rsd_set );

///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{

	using namespace core::pose;
	using namespace protocols::farna;

	RNA_DeNovoSetupOP rna_de_novo_setup( new RNA_DeNovoSetup );
	rna_de_novo_setup->initialize_from_command_line();
	Pose & pose = *( rna_de_novo_setup->pose() );

	protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol( rna_de_novo_setup->options(),
		rna_de_novo_setup->rna_params() );
	rna_de_novo_protocol.set_native_pose( rna_de_novo_setup->native_pose() );
	rna_de_novo_protocol.set_refine_pose_list( rna_de_novo_setup->refine_pose_list() );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );
	// protocols::moves::AddPyMolObserver( pose, false, 0.01);

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
		option.add_relevant( score::weights );
		option.add_relevant( rna::farna::sequence );
		option.add_relevant( rna::farna::secstruct );
		option.add_relevant( rna::farna::secstruct_file );
		option.add_relevant( rna::farna::minimize_rna );
		option.add_relevant( rna::farna::fixed_stems );
		option.add_relevant( rna::farna::obligate_pair );
		option.add_relevant( rna::farna::obligate_pair_explicit );
		// option.add_relevant( rna::farna::relax_rna );
		// option.add_relevant( rna::farna::simple_relax );
		// option.add_relevant( rna::farna::ignore_secstruct );
		option.add_relevant( rna::farna::lores_scorefxn );
		option.add_relevant( rna::farna::cycles );
		option.add_relevant( rna::farna::temperature );
		option.add_relevant( rna::farna::jump_change_frequency );
		option.add_relevant( rna::farna::close_loops );
		option.add_relevant( rna::farna::close_loops_after_each_move );
		option.add_relevant( rna::farna::output_lores_silent_file );
		option.add_relevant( rna::farna::heat );
		option.add_relevant( rna::farna::dump );
		option.add_relevant( rna::farna::staged_constraints );
		// option.add_relevant( rna::farna::jump_library_file );
		option.add_relevant( rna::farna::params_file );
		option.add_relevant( rna::farna::filter_lores_base_pairs );
		option.add_relevant( rna::farna::filter_lores_base_pairs_early );
		option.add_relevant( rna::farna::filter_chain_closure );
		option.add_relevant( rna::farna::filter_chain_closure_halfway );
		option.add_relevant( rna::farna::filter_chain_closure_distance );
		option.add_relevant( basic::options::OptionKeys::rna::farna::vall_torsions );
		option.add_relevant( rna::farna::use_1jj2_torsions );
		option.add_relevant( rna::farna::rna_lores_chainbreak_weight );
		option.add_relevant( rna::farna::rna_lores_linear_chainbreak_weight );
		option.add_relevant( rna::farna::allow_bulge  );
		option.add_relevant( rna::farna::allowed_bulge_res );
		option.add_relevant( rna::farna::allow_consecutive_bulges );
		option.add_relevant( rna::farna::binary_output );
		option.add_relevant( rna::farna::move_first_rigid_body );
		option.add_relevant( rna::farna::root_at_first_rigid_body );
		option.add_relevant( rna::farna::suppress_bp_constraint );
		option.add_relevant( rna::farna::output_filters );
		option.add_relevant( rna::farna::autofilter );
		option.add_relevant( rna::farna::output_res_num );
		option.add_relevant( rna::farna::offset );
		option.add_relevant( rna::farna::tag );
		option.add_relevant( rna::farna::refine_silent_file );
		option.add_relevant( rna::farna::refine_native );
		option.add_relevant( rna::farna::bps_moves );
		option.add_relevant( rna::farna::minimize::minimizer_use_coordinate_constraints );
		option.add_relevant( rna::farna::minimize::extra_minimize_res );
		option.add_relevant( rna::farna::minimize::extra_minimize_chi_res );
		option.add_relevant( rna::farna::minimize::minimize_bps );
		option.add_relevant( full_model::cutpoint_closed );
		option.add_relevant( full_model::cutpoint_open );
		option.add_relevant( rna::vary_geometry );
		option.add_relevant( basic::options::OptionKeys::rna::data_file );

		option.add_relevant( constraints::cst_file );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_BASE" ); // for chemical mapping.
		option[ OptionKeys::chemical::patch_selectors ].push_back( "VIRTUAL_SIDE_CHAIN" ); // for proteins for vdw screen

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

