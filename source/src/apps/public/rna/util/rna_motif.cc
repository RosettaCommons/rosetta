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
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rna/RNA_Motif.hh>
#include <core/scoring/rna/RNA_LowResolutionPotential.hh>
#include <core/scoring/rna/RNA_FilteredBaseBaseInfo.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/pose/Pose.hh>
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

static THREAD_LOCAL basic::Tracer TR( "rna_motif" );

using namespace basic::options::OptionKeys;
using namespace basic::options;
using namespace utility;
using namespace core;
using namespace core::scoring::rna;

///////////////////////////////////////////////////////////////////////////////////////////////////////
// @details
//
//  Identify RNA motifs in a model, and output coloring scripts for Pymol.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////

std::map< RNA_MotifType, std::string > motif_color =
{ {U_TURN, "lightblue"},
{UA_HANDLE, "marine" },
{T_LOOP, "tv_blue" },
{INTERCALATED_T_LOOP, "deepblue" },
{LOOP_E_SUBMOTIF, "salmon" },
{BULGED_G, "red" },
{GNRA_TETRALOOP, "ruby"},  // GNRA is already pretty favorable.
{STRICT_WC_STACKED_PAIR, "gray20"},
{WC_STACKED_PAIR, "gray50"},
{A_MINOR, "gold"},
{PLATFORM, "sand"},
{TL_RECEPTOR, "limon"},
{TETRALOOP_TL_RECEPTOR, "orange"}
};

// @brief super-simple helper function for PyMOL commands.
// @details for speed, could also define PyMOL boolean property and then color things at end based on property.
void
output_motifs_to_pymol( std::ostream & out,
	pose::Pose const & pose,
	RNA_Motifs const & rna_motifs ) {
	std::string tag( tag_from_pose( pose ) );
	tag = replace_in( tag, ".pdb", "" );
	for ( auto const & motif : rna_motifs ) {
		for ( auto const & res : motif ) {
			out << "color " << motif_color[ motif.type() ] << ", " << tag << " and chain " << pose.pdb_info()->chain( res ) << " and resi " << pose.pdb_info()->number( res ) << std::endl;

		}
	}
}

void
rna_motif_test()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::scoring;
	using namespace core::scoring::rna;
	using namespace core::import_pose::pose_stream;
	using namespace core::io::silent;

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

	pose::Pose pose;
	RNA_LowResolutionPotential potential;
	ScoreFunctionOP denovo_scorefxn( ScoreFunctionFactory::create_score_function( RNA_LORES_WTS ) );
	std::string const pymol_command_file( "color_motifs.pml" );
	utility::io::ozstream pymol_out( pymol_command_file );
	while ( input->has_another_pose() ) {
		input->fill_pose( pose, *rsd_set );
		(*denovo_scorefxn)( pose );
		std::cout << tag_from_pose( pose ) << std::endl;
		RNA_Motifs const rna_motifs = get_rna_motifs( pose, potential,
			rna_scoring_info_from_pose( pose ).rna_filtered_base_base_info() );
		output_rna_motifs( pose, rna_motifs );
		std::cout << std::endl;
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
		devel::init(argc, argv);

		protocols::viewer::viewer_main( my_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
