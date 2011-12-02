// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @detailed
///	  Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
#include <protocols/evaluation/StructuralSimilarityEvaluator.hh>
#include <protocols/evaluation/RDC_Evaluator.hh>
#include <protocols/evaluation/ScoreEvaluator.hh>
#include <protocols/evaluation/JScoreEvaluator.hh>
#include <protocols/evaluation/ContactMapEvaluator.hh>
#include <protocols/evaluation/JumpEvaluator.hh>
#include <protocols/evaluation/CamShiftEvaluator.hh>
#include <protocols/evaluation/PalesEvaluator.hh>
#include <protocols/evaluation/PredictedBurialEvaluator.hh>
#include <protocols/loops/Loops.hh>
#include <core/scoring/constraints/util.hh>
#include <protocols/moves/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Project Headers
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/SequenceAlignment.fwd.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/string.functions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>

#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/io/silent/SilentStructFactory.hh>
//// C++ headers

// due to template function
#include <core/io/silent/SilentStruct.hh>


// option key includes
#include <basic/options/option_macros.hh>
#include <basic/options/keys/evaluation.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/fragment/SecondaryStructure.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/evaluation/util.hh>
#include <utility/vector0.hh>

#ifdef WIN32
	#include <core/scoring/constraints/Constraint.hh>
#endif


static basic::Tracer tr("protocols.evalution");

namespace protocols {
namespace evaluation {
using namespace core;

static bool options_registered_=false;

void register_options() {
	using namespace basic::options;
	if ( options_registered_ ) return;
	options_registered_ = true;
	OPT( evaluation::rmsd );
	OPT( evaluation::gdtmm );
	OPT( evaluation::rdc );
	OPT( evaluation::pool );
	OPT( evaluation::constraints );
	OPT( in::file::native );
	OPT( evaluation::chemical_shifts );
}

void invert_include_residues( Size nres, core::scoring::ResidueSelectionVector const& include_list, core::scoring::ResidueSelectionVector& exclude_list ) {

	exclude_list.clear();

	for ( Size ir = 1; ir <= nres; ++ir ) {
		bool include_residue = false;
		for ( Size ex = 1; ex <= include_list.size(); ex ++ ) {
			if ( include_list[ex] == ir ) {
				include_residue = true;
				break;
			}
		}

		if ( !include_residue ) {
			exclude_list.push_back( ir );
		}
	} // for ( Size ir = 1; ir <= native_pose.total_residue(); ++ir )
}

//@detail find residues that don't have missing density
void find_existing_residues(  core::pose::PoseCOP pose, std::string tag, core::scoring::ResidueSelection& selection ) {
	for ( Size pos = 1; pos <= pose->total_residue(); pos++ ) {
		if ( pose->residue_type( pos ).is_protein() && pose->residue_type( pos ).has("CA") ) {
			numeric::xyzVector< core::Real> ca_pos = pose->residue( pos ).atom("CA").xyz();
			bool good ( true );
			for ( Size j=1; j<= pose->residue( pos ).natoms(); ++j ) {
				if ( ( ca_pos - pose->residue( pos ).atom(j).xyz() ).length() > 20 ) {
					good = false;
				}
			}
			if ( good ) selection.push_back( pos );
		}
	}
	if ( tr.Trace.visible() ) {
		tr.Trace << "selection of residues for rmsd of " << tag << std::endl;
		for ( std::list< core::Size >::const_iterator it = selection.begin(), eit = selection.end();
					it != eit; ++it ) {
			tr.Trace << " " << *it;
		}
		tr.Trace << std::endl;
	}
}

void evaluate_pose( core::pose::Pose& pose, PoseEvaluator& eval, std::ostream& os ) {
		//		ProteinSilentStruct pss;
		io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
		pss->fill_struct( pose, "eval" );
		eval.apply( pose, "eval", *pss );
		os << "\n";
		pss->print_score_header( os );
		os << "\n";
		pss->print_scores( os );
		os << std::endl;
}


void define_scorable_core_from_secondary_structure(
   core::fragment::SecondaryStructure const& ss_def,
	 protocols::loops::Loops& scored_core )
{
	using namespace core;
	using namespace basic::options;
	//	Size const max_loop_size( 3 );
	//	Size const max_short_helix( 5 );
	Size const max_loop_size( option[ OptionKeys::evaluation::score_sscore_maxloop ]() );
	Size const max_short_helix( option[ OptionKeys::evaluation::score_sscore_short_helix ]() );

	//find residues that are part of a short helix -- less than or equal to 5 residues
	utility::vector1< bool > short_helix( ss_def.total_residue(), false );

	//selection of loop definitions...
	//these loops define regions that are scored. Add all loops that are 4 residues or longer.
	//subsequently add also helices that have fewer than 6 residues if they terminated a long loop (>=4)
	loops::Loops unscored_loops;

	for ( Size pos=1; pos <= ss_def.total_residue(); pos++ ) {

		//detect loops
		if ( ss_def.loop_fraction( pos ) > 0.1 ) {
			//go to end of loop
			Size lpos = 1;
			for ( ; ( lpos+pos <= ss_def.total_residue() ) && ( ss_def.loop_fraction( pos+lpos ) > 0.1); ++lpos ) {}
			if ( lpos > max_loop_size ) { //this loop has 4 or more residues
				unscored_loops.add_loop( pos, pos+lpos-1 );
			}
			pos+=lpos-1;
		} // have found a loop

		// look for short helices and store in short_helix
		if ( ss_def.helix_fraction( pos ) > 0.1 ) {
			Size hpos = 1;
			for ( ; ( hpos+pos <= ss_def.total_residue() ) && ( ss_def.helix_fraction( pos+hpos ) > 0.1); ++hpos ) {}
			if ( hpos <= max_short_helix   ) { //this helix has 5 or fewer residues
				for ( Size ipos = 0; ipos < hpos; ++ipos ) {
					short_helix[ pos+ipos] = true;
				}
			}
		}

		//finished parsing secondary structure definition
	}

	//elongate loops if they are terminated by a short helix
	loops::Loops removed_short_helices( unscored_loops );
	for ( loops::Loops::const_iterator it=unscored_loops.begin(); it != unscored_loops.end(); ++it ) {
		Size npos( it->stop() + 1 );
		while ( short_helix[ npos ] ) {
			removed_short_helices.add_loop( npos-1, npos );
			npos++;
		}
	}

	scored_core =	removed_short_helices.invert( ss_def.total_residue() );
}


}
}
