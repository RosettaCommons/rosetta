// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/setup_util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/recces/setup_util.hh>
#include <protocols/recces/options/RECCES_Options.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/rna/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/stepwise/modeler/rna/helix/RNA_HelixAssembler.hh>
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/farna/secstruct/RNA_SecStructLegacyInfo.hh>
#include <basic/Tracer.hh>

#include <utility/io/ozstream.hh>

static basic::Tracer TR( "protocols.recces.setup_util" );

using namespace core;
using namespace core::chemical;
using namespace core::pose;

namespace protocols {
namespace recces {


//////////////////////////////////////////////////////////////////////////////
PoseOP
recces_pose_setup( options::RECCES_Options const & options )
{
	if ( options.legacy_turner_mode() ) {

		TR << TR.Green << "Assuming RECCES Turner mode, due to specification of -seq1" << std::endl;
		runtime_assert( options.seq1().size() > 0 );
		return pose_setup_turner( options.seq1(), options.seq2() );

	} else {

		runtime_assert( options.infile().size()  > 0 );
		PoseOP pose = pose_setup_from_file( options );

		return pose;
	}
}

//////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
pose_setup_turner(
	std::string const & seq1,
	std::string const & seq2)
{
	using namespace core::chemical;
	using namespace core::pose::rna;
	using namespace protocols::farna::secstruct;

	protocols::stepwise::modeler::rna::helix::RNA_HelixAssembler assembler;
	assembler.use_phenix_geo( true );
	PoseOP pose( assembler.build_init_pose( seq1, seq2 ) );
	virtualize_5prime_phosphates( *pose );

	// Need to specify which nucleotides get 'A-form' range, which will
	//  also determine how large the Gaussian steps are during sampling.
	std::string secstruct_legacy;
	Size const len1( seq1.size() );
	Size const len2( seq2.size() );
	Size const n_bp( std::min( len1, len2 ) );
	Size const total_len( len1 + len2 );
	for ( Size i = 1; i <= total_len; ++i ) {
		if ( i > n_bp && i <= total_len - n_bp ) secstruct_legacy.push_back( 'X' );
		else secstruct_legacy.push_back( 'H' );
	}
	set_rna_secstruct_legacy( *pose, secstruct_legacy );

	return pose;
}


	//////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
pose_setup_from_file( options::RECCES_Options const & options )
{

	// PoseOP pose = stepwise::setup::get_pdb_and_cleanup( options.infile() );
	// // needs to accept command-line input -- need to guess sequence info if fasta not specified!
	// // also -- force sample_res all on if not specified from command-line. [different from stepwise default behavior]
 	// stepwise::setup::fill_full_model_info_from_command_line( *pose );
	// stepwise::modeler::fix_up_residue_type_variants( *pose ); //virtualizes phosphates, etc.

	using namespace core::chemical::rna;
	// Pose setup
	PoseOP pose( new Pose );
	// should be starter base pair (e.g., "cg.pdb") -- only 2 residues!
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	import_pose::pose_from_file( *pose, *rsd_set_op, options.infile(), core::import_pose::PDB_file);

	TR << "Annotated sequence of pose from " << options.infile() << ": " << pose->annotated_sequence() << std::endl;

	// replace this with fold_tree, full_model_info definition, + cleanup_variants
	if ( pose->size() == 2 ) {
	 	kinematics::FoldTree f( 2 );
	 	f.new_jump( 1, 2, 1 );
	 	f.set_jump_atoms( 1, default_jump_atom(pose->residue_type(1)), default_jump_atom(pose->residue_type(2) )) ;
	 	pose->fold_tree( f );
	 	for ( Size n = 1; n <= 2; n++ ) {
	 		pose::add_variant_type_to_pose_residue( *pose, VIRTUAL_PHOSPHATE, n );
	 		pose::add_variant_type_to_pose_residue( *pose, VIRTUAL_RIBOSE, n );
	 		if ( options.block_stack() ) {
	 			pose::add_variant_type_to_pose_residue( *pose, BLOCK_STACK_ABOVE, n );
	 			pose::add_variant_type_to_pose_residue( *pose, BLOCK_STACK_BELOW, n );
	 		}
	 	}
	}

	return pose;
}

} //recces
} //protocols
