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
/// @author Rhiju Das

// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/util.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <basic/options/option.hh>
#include <basic/database/open.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/rna/RNA_DataReader.hh>
#include <utility/vector1.hh>
#include <ObjexxFCL/string.functions.hh>

//RNA stuff.
#include <protocols/farna/RNA_DeNovoProtocol.hh>
#include <protocols/farna/RNA_DeNovoProtocolOptions.hh>
#include <protocols/farna/util.hh>
//#include <protocols/moves/PyMolMover.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/sequence/Sequence.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>


using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

vector1<pose::PoseOP> get_refine_pose_list( std::string const & input_silent_file,
																						vector1< Size > const & output_res,
																						core::chemical::ResidueTypeSetCOP rsd_set );

///////////////////////////////////////////////////////////////////////////////
void
rna_denovo_test()
{

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace protocols::farna;

	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	// Some initialization
	/////////////////////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////////////////////
	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	pose::PoseOP native_pose_OP( new pose::Pose );
	pose::Pose & native_pose = *native_pose_OP;

	pose::Pose extended_pose;

	std::string const in_path = option[ in::path::path ]()[1];

	bool native_exists = false;
	//Read in native if it exists.
	if ( option[ in::file::native ].user() ) {
		native_exists = true;
		//Read in native if it exists.
		std::string native_pdb_file  = option[ in::file::native ];
		core::import_pose::pose_from_pdb( native_pose, *rsd_set, in_path + native_pdb_file );
	} else {
		runtime_assert( !option[ OptionKeys::rna::farna::refine_native ]() );
	}

	//Prepare starting structure from scratch --> read from fasta.
	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( in_path + fasta_file )[1];
	if ( option[ OptionKeys::rna::farna::refine_native ]() ) {
		extended_pose = native_pose;
	} else {
		core::pose::make_pose_from_sequence( extended_pose, fasta_sequence->sequence(), *rsd_set );
	}
	set_output_res_num( extended_pose, option[ OptionKeys::rna::farna::output_res_num ]() );

	// Silent file input for fine refinement
	vector1<pose::PoseOP> refine_pose_list = get_refine_pose_list( option[ OptionKeys::rna::farna::refine_silent_file ](), option[ OptionKeys::rna::farna::output_res_num ](), rsd_set );

	//Score these suckers.
	pose::Pose pose( extended_pose );
	ScoreFunctionOP scorefxn = ScoreFunctionFactory::create_score_function( RNA_LORES_WTS );

	RNA_DeNovoProtocolOptionsOP options( new RNA_DeNovoProtocolOptions);
	options->initialize_from_command_line();
	if ( option[ OptionKeys::rna::farna::refine_native ]() ) options->set_refine_pose( true );

	protocols::farna::RNA_DeNovoProtocol rna_de_novo_protocol( options );
	if ( native_exists ) rna_de_novo_protocol.set_native_pose( native_pose_OP );
	rna_de_novo_protocol.set_refine_pose_list( refine_pose_list );

	if ( option[ OptionKeys::constraints::cst_file ].user() ) {
		ConstraintSetOP cst_set = ConstraintIO::get_instance()->read_constraints( option[ OptionKeys::constraints::cst_file ](1), ConstraintSetOP( new ConstraintSet ), pose );
		pose.constraint_set( cst_set );
		for ( Size i = 1; i <= refine_pose_list.size(); ++i ) refine_pose_list[i]->constraint_set( cst_set );
	}

	if ( option[ OptionKeys::rna::farna::data_file].user() ) {
		core::io::rna::RNA_DataReader rna_data_reader( in_path + option[ OptionKeys::rna::farna::data_file ]  );
		rna_data_reader.fill_rna_data_info( pose );
	}

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );

	// protocols::moves::AddPyMolObserver( pose, false, 0.01);

	rna_de_novo_protocol.apply( pose );

}

///////////////////////////////////////////////////////////////
vector1<pose::PoseOP>
get_refine_pose_list( std::string const & input_silent_file,
											vector1< Size > const & output_res,
											core::chemical::ResidueTypeSetCOP rsd_set )
{

	vector1<pose::PoseOP> refine_pose_list;
	if ( input_silent_file.size() > 0 ) {
		core::import_pose::pose_stream::SilentFilePoseInputStream input( input_silent_file );
		input.set_order_by_energy( true );
		while ( input.has_another_pose() ) {
			pose::PoseOP new_pose( new pose::Pose );
			input.fill_pose( *new_pose, *rsd_set );
			protocols::farna::set_output_res_num( *new_pose, output_res );
			refine_pose_list.push_back( new_pose );
		}
	}
	return refine_pose_list;
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
		option.add_relevant( in::file::input_res );
		option.add_relevant( in::file::native );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( score::weights );
		option.add_relevant( constraints::cst_file );
		option.add_relevant( rna::farna::minimize_rna );
		option.add_relevant( rna::farna::relax_rna );
		option.add_relevant( rna::farna::simple_relax );
		option.add_relevant( rna::farna::ignore_secstruct );
		option.add_relevant( rna::farna::lores_scorefxn );
		option.add_relevant( rna::farna::filter_lores_base_pairs );
		option.add_relevant( rna::farna::filter_lores_base_pairs_early );
		option.add_relevant( rna::farna::filter_chain_closure );
		option.add_relevant( rna::farna::filter_chain_closure_halfway );
		option.add_relevant( rna::farna::filter_chain_closure_distance );
		option.add_relevant( rna::farna::cycles );
		option.add_relevant( rna::farna::temperature );
		option.add_relevant( rna::farna::jump_change_frequency );
		option.add_relevant( rna::farna::close_loops );
		option.add_relevant( rna::farna::close_loops_after_each_move );
		option.add_relevant( rna::farna::output_lores_silent_file );
		option.add_relevant( rna::farna::heat );
		option.add_relevant( rna::farna::dump );
		option.add_relevant( rna::farna::staged_constraints );
		option.add_relevant( rna::farna::jump_library_file );
		option.add_relevant( rna::farna::params_file );
		option.add_relevant( rna::farna::use_1jj2_torsions );
		option.add_relevant( rna::farna::rna_lores_chainbreak_weight );
		option.add_relevant( rna::farna::rna_lores_linear_chainbreak_weight );
		option.add_relevant( rna::farna::allow_bulge  );
		option.add_relevant( rna::farna::allowed_bulge_res );
		option.add_relevant( rna::farna::extra_minimize_res );
		option.add_relevant( rna::farna::extra_minimize_chi_res );
		option.add_relevant( rna::farna::allow_consecutive_bulges );
		option.add_relevant( rna::farna::binary_output );
		option.add_relevant( rna::farna::move_first_rigid_body );
		option.add_relevant( rna::farna::root_at_first_rigid_body );
		option.add_relevant( rna::farna::suppress_bp_constraint );
		option.add_relevant( rna::farna::output_filters );
		option.add_relevant( rna::farna::autofilter );
		option.add_relevant( rna::farna::output_res_num );
		option.add_relevant( rna::farna::refine_silent_file );
		option.add_relevant( rna::farna::refine_native );
		option.add_relevant( rna::farna::bps_moves );
		option.add_relevant( rna::farna::minimizer_use_coordinate_constraints );

		option.add_relevant( constraints::cst_file );
		option.add_relevant( basic::options::OptionKeys::rna::farna::vary_geometry );
		option.add_relevant( basic::options::OptionKeys::rna::farna::vall_torsions );
		option.add_relevant( basic::options::OptionKeys::rna::farna::data_file );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		core::init::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

