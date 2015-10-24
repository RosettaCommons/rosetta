// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/public/cyclic_peptide/simple_cycpep_predict.cc
/// @brief Predicts structure of N-to-C cyclic peptides from amino acid sequence.
/// @details
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifdef BOINC
#include <utility/boinc/boinc_util.hh>
#include <protocols/boinc/boinc.hh>
#include "boinc_zip.h"
#endif // BOINC -- keep this first.

//General includes
#include <basic/options/option.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <devel/init.hh>
#include <utility/exit.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/pose/PDBInfo.hh>
#include <numeric/conversions.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/constants.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/filters/Filter.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>
#include <protocols/filters/BasicFilters.hh>
#include <protocols/protein_interface_design/filters/HbondsToResidueFilter.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/scoring/Energies.hh>

//Application-specific includes
#include <protocols/cyclic_peptide_predict/SimpleCycpepPredictApplication.hh>

#include <utility/vector1.hh>
#include <stdio.h>

//Tracer:
static THREAD_LOCAL basic::Tracer TR( "apps.public.cyclic_peptide.simple_cycpep_predict" );

int
main( int argc, char * argv [] )
{
	try {
		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication::register_options(); //Has to be before init().
		devel::init(argc, argv);

		if ( TR.visible() ) {
			TR << "Starting simple_cycpep_predict.cc" << std::endl;
			TR << "Pilot app created 16 September 2015 by Vikram K. Mulligan, Ph.D., Baker laboratory." << std::endl;
			TR << "For questions, contact vmullig@uw.edu." << std::endl << std::endl;
			TR.flush();
		}

		protocols::cyclic_peptide_predict::SimpleCycpepPredictApplication this_app;
		this_app.run();

	} catch ( utility::excn::EXCN_Base& excn ) {
		TR.Error << "Exception caught: " << std::endl;
		excn.show( TR.Error );
		TR.Error.flush();
		return -1;
	}

	if ( TR.visible() ) {
		TR << "Finished simple_cycpep_predict.cc.  Exiting." << std::endl;
		TR.flush();
	}

	return 0;
}
