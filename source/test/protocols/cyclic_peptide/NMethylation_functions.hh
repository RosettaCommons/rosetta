// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  protocols/cyclic_peptide/NMethylation_functions.cxxtest.hh
/// @brief  Code for N-methylation unit tests (shared between beta_nov15 and regular tests).
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

#ifndef INCLUDED_test_protocols_cyclic_peptide_NMethylation_functions_HH
#define INCLUDED_test_protocols_cyclic_peptide_NMethylation_functions_HH

// Test headers

// Project Headers
#include <protocols/cyclic_peptide/DeclareBond.hh>
#include <protocols/cyclic_peptide/PeptideStubMover.hh>
#include <protocols/simple_moves/MutateResidue.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/PackerTask_.hh>

// Core Headers
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

// Protocol Headers
#include <basic/Tracer.hh>

namespace test {
namespace protocols {
namespace cyclic_peptide {

static basic::Tracer TR("protocols.cyclic_peptide.NMethylation_functions");

class NMethylationTests_functions {

public:

	/// @brief Just build an N-methylated peptide, manipulate it a bit, and do a FastRelax.
	///
	void common_test_linear_nmethyl_peptide() {
		using namespace ::protocols::cyclic_peptide;
		using namespace ::protocols::simple_moves;
		using namespace ::protocols::relax;

		core::pose::PoseOP mypose( new core::pose::Pose );
		core::scoring::ScoreFunctionOP scorefxn( core::scoring::get_score_function() ); //Get whatever the current scorefunction is.
		//scorefxn->set_weight( core::scoring::fa_dun, 0.0 ); //DELETE ME

		PeptideStubMoverOP maker( new PeptideStubMover );
		maker->set_reset_mode(true);
		maker->add_residue("Append", "ALA", 0, false, "", 15, 0, "");
		maker->apply(*mypose);

		MutateResidueOP mutres5( new MutateResidue( 5, "TRP:N_Methylation" ) );
		mutres5->set_update_polymer_dependent( true );
		mutres5->apply(*mypose);
		MutateResidueOP mutres6( new MutateResidue( 6, "TRP:N_Methylation" ) );
		mutres6->set_update_polymer_dependent( true );
		mutres6->apply(*mypose);
		MutateResidueOP mutres7( new MutateResidue( 7, "TRP:N_Methylation" ) );
		mutres7->set_update_polymer_dependent( true );
		mutres7->apply(*mypose);

		DeclareBondOP termini( new DeclareBond );
		termini->set( 2, "C", 3, "N", true ); //Note: if you change this to 1 and 2 (instead of 2 and 3), you reliably recreate that "cannot normalize xyzvector of length zero" error that we really need to fix...
		termini->apply( *mypose );

		for ( core::Size i=1, imax=mypose->total_residue(); i<=imax; ++i ) {
			if ( i>1 ) mypose->set_phi( i, -60 );
			if ( i<imax ) {
				mypose->set_psi( i, 60 );
				mypose->set_omega( i, 180 );
			}
		}
		mypose->update_residue_neighbors();

		//mypose->dump_pdb( "vtemp_nmethyl_pre_pack.pdb" ); //DELETE ME

		core::pose::PoseOP mypose2( mypose->clone() );

		core::pack::task::PackerTaskOP task( new core::pack::task::PackerTask_(*mypose2) );
		task->restrict_to_repacking();
		PackRotamersMoverOP pack( new PackRotamersMover( scorefxn, task ) );
		pack->apply(*mypose2);
		//mypose2->dump_pdb( "vtemp_nmethyl_post_pack.pdb" ); //DELETE ME

		core::pose::PoseOP mypose3( mypose2->clone() );
		FastRelaxOP frlx( new FastRelax( scorefxn, 2 ) );
		frlx->apply(*mypose3);
		//mypose3->dump_pdb( "vtemp_nmethyl_post_frlx.pdb" ); //DELETE ME
	}

};

} //namespace cyclic_peptide
} //namespace protocols
} //namespace test

#endif //INCLUDED_test_protocols_cyclic_peptide_NMethylation_functions_HH
