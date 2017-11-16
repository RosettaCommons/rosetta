// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file rna_motif.cc
/// @brief check if we can annotate U-turns, A-minor motifs, etc. in an RNA

// libRosetta headers
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Motif.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/sequence/util.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>

// option key includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <utility/io/ozstream.hh>

#include <ObjexxFCL/string.functions.hh>

static basic::Tracer TR( "rna_motif" );

using namespace basic::options::OptionKeys;
using namespace basic::options;
using namespace utility;
using namespace core;
using namespace core::scoring::rna;
using core::pose::ChainSegID;

OPT_KEY( StringVector, chains )

///////////////////////////////////////////////////////////////////////////////////////////////////////
// @details
//
//  Identify RNA motifs in a model, and output coloring scripts for Pymol.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////
void
rna_motif_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::sequence;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::silent;
	using namespace core::pose::rna;

	ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				) );
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	pose::Pose full_pose, pose;
	RNA_LowResolutionPotential potential;
	ScoreFunctionOP hires_scorefxn( new ScoreFunction );
	hires_scorefxn->set_weight( hbond_sc, 1 );
	ScoreFunctionOP denovo_scorefxn( ScoreFunctionFactory::create_score_function( RNA_LORES_WTS ) );
	std::string const pymol_command_file( "color_motifs.pml" );
	utility::io::ozstream pymol_out( pymol_command_file );
	while ( input->has_another_pose() ) {

		// get the pose.
		input->fill_pose( full_pose, *rsd_set );

		// look for protein, HOH, ions, etc. contacting RNA chains of interest.
		utility::vector1< ChainSegID > chain_segids = figure_out_rna_chains( full_pose, option[ chains ]() );
		std::string const out_ligand_file = tag_from_pose( full_pose )+ ".ligands.txt";
		utility::io::ozstream out_ligand( out_ligand_file );
		(*hires_scorefxn)( full_pose );
		output_ligands( out_ligand, full_pose, chain_segids );
		TR << "Created: " << out_ligand_file << std::endl;

		// extract *only RNA* chains of interest
		// get ready for RNA motif identification.
		pose = extract_rna_chains( full_pose, chain_segids );

		(*denovo_scorefxn)( pose );
		std::cout << tag_from_pose( pose ) << std::endl;
		RNA_Motifs const rna_motifs = get_rna_motifs( pose, potential,
			rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info() );
		output_rna_motifs_detailed( pose, rna_motifs );
		std::cout << std::endl;

		std::string const out_bp_file = tag_from_pose( pose )+ ".base_pairs.txt";
		utility::io::ozstream out_bp( out_bp_file );
		output_base_pairs( out_bp, rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info().base_pair_list(), pose );
		TR << "Created: " << out_bp_file << std::endl;
		out_bp.close();

		std::string const out_stack_file = tag_from_pose( pose )+ ".stacks.txt";
		utility::io::ozstream out_stack( out_stack_file );
		output_base_stacks( out_stack, rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info().base_stack_list(), pose );
		TR << "Created: " << out_stack_file << std::endl;
		out_stack.close();

		std::string const out_othercontacts_file = tag_from_pose( pose ) + ".other_contacts.txt";
		utility::io::ozstream out_othercontacts( out_othercontacts_file );
		output_other_contacts( out_othercontacts, pose );
		TR << "Created: " << out_othercontacts_file << std::endl;
		out_stack.close();

		std::string const out_stem_file = tag_from_pose( pose )+ ".stems.txt";
		utility::io::ozstream out_stem( out_stem_file );
		output_stems( out_stem, rna_motifs, pose );
		TR << "Created: " << out_stem_file << std::endl;
		out_stem.close();

		std::string const out_motif_file = tag_from_pose( pose ) + ".motifs.txt";
		utility::io::ozstream out_motif( out_motif_file );
		output_rna_motifs( out_motif, pose, rna_motifs );
		TR << "Created: " << out_motif_file << std::endl;

		std::string const out_fasta_file = tag_from_pose( pose ) + ".fasta";
		output_fasta_file( out_fasta_file, pose );
		TR << "Created: " << out_fasta_file << std::endl;

		// output motifs to Pymol.
		output_motifs_to_pymol( pymol_out, pose, rna_motifs );

	}
	pymol_out.close();
	TR << "Created: " << pymol_command_file << std::endl;
}

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	rna_motif_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		NEW_OPT( chains, "list of chains or chain:segid over which to assess RNA motifs", utility::vector1< std::string >() );
		option.add_relevant( in::file::s );
		option.add_relevant( in::file::silent );

		devel::init(argc, argv);
		protocols::viewer::viewer_main( my_main );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

