// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file relax_initialization_protocols
/// @brief initialization protocols for relax
/// @details
///   Contains currently: Classic Abinitio
///
///
/// @author Oliver Lange


// Unit Headers

// Package Headers
#include <protocols/evaluation/PoseEvaluator.fwd.hh>
#include <protocols/evaluation/PoseEvaluator.hh>
#include <core/scoring/constraints/util.hh>
//#include <protocols/canonical_sampling/mc_convergence_checks/Pool_ConvergenceCheck.hh>

// Project Headers
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/pose/Pose.hh>

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

#include <core/fragment/SecondaryStructure.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/evaluation/util.hh>
#include <utility/vector0.hh>

#ifdef WIN32
#include <core/scoring/constraints/Constraint.hh>
#endif


static THREAD_LOCAL basic::Tracer tr( "protocols.evaluation" );

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
	//  ProteinSilentStruct pss;
	io::silent::SilentStructOP pss = io::silent::SilentStructFactory::get_instance()->get_silent_struct_out();
	pss->fill_struct( pose, "eval" );
	eval.apply( pose, "eval", *pss );
	os << "\n";
	pss->print_score_header( os );
	os << "\n";
	pss->print_scores( os );
	os << std::endl;
}


}
}
