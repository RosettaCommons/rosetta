// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

///@brief Export all relevant torsions for each nt in an RNA

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/util.hh>
#include <numeric/xyz.functions.hh>

#include <core/chemical/rna/util.hh>
#include <core/id/TorsionID.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
//#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/options/option_macros.hh>

// C++ headers
#include <string>
#include <sstream>

using namespace core;
using namespace conformation;
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace core::id;
using namespace core::scoring;
using namespace pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using basic::Error;
using basic::Warning;
using utility::file::FileName;

// tracer - used to replace cout
static basic::Tracer TR("TorsionScoring");

OPT_KEY( String, resname )

int
main( int argc, char* argv[] )
{
	try {
		NEW_OPT( resname, "name of residue on which to focus", "ALL" );
		devel::init(argc, argv);

		scoring::ScoreFunctionOP score_fxn( new scoring::ScoreFunction );

		score_fxn->set_weight( core::scoring::rna_torsion, 1 );
		std::string const residue_name = option[ resname ].value();

		Pose pose;
		Size const n_files = option[in::file::s].value().size();
		auto const fns = option[in::file::s].value();
		for ( Size ii = 1; ii <= n_files; ++ii ) {
			std::string filename = fns[ii];

			import_pose::pose_from_file( pose, filename , core::import_pose::PDB_file);
			TR << "Importing pose from " << filename << std::endl;
			filename.erase( filename.size()-4 );

			Real score = ( *score_fxn )( pose );
			TR << "Total torsion score: " << score << std::endl;
			// TODO: chainbreaks
			for ( Size ii = 2; ii <= pose.size() - 1; ++ii ) {
				if ( !pose.residue_type( ii - 1 ).is_RNA() ) continue;
				if ( !pose.residue_type( ii ).is_RNA() ) continue;
				if ( !pose.residue_type( ii + 1 ).is_RNA() ) continue;

				// Skip non res
				//TR << "residue_name " << residue_name << " " << pose.residue_type( ii ).name3() << std::endl;
				if ( residue_name != "ALL" && pose.residue_type( ii ).name3() != residue_name ) continue;

				TR << ii << " " << pose.residue_type( ii ).name3() << " "
					<< pose.residue_type( ii-1 ).name3() << " " <<  pose.residue_type( ii+1 ).name3() // flanking
					<< " " << pose.torsion( TorsionID( ii, id::BB, BETA  ) ) // intra
					<< " " << pose.torsion( TorsionID( ii, id::BB, GAMMA ) ) // intra
					<< " " << pose.torsion( TorsionID( ii, id::BB, DELTA ) ) // intra
					<< " " << pose.torsion( TorsionID( ii, id::CHI, 1 ) ) // intra, chi1
					<< " " << pose.torsion( TorsionID( ii, id::CHI, 2 ) ) // intra, nu2
					<< " " << pose.torsion( TorsionID( ii, id::CHI, 3 ) ) // intra, nu1
					<< " " << pose.torsion( TorsionID( ii, id::CHI, 4 ) ) // intra, proton chi
					<< " " << pose.torsion( TorsionID( ii, id::BB, ALPHA ) ) // inter
					<< " " << pose.torsion( TorsionID( ii - 1, id::BB, ZETA    ) ) // inter
					<< " " << pose.torsion( TorsionID( ii - 1, id::BB, EPSILON ) ) // inter
					<< " " << pose.torsion( TorsionID( ii - 1, id::BB, DELTA   ) ) // inter
					// << " " << pose.torsion( TorsionID( ii, id::BB, DELTA ) ) // inter
					<< " " << pose.torsion( TorsionID( ii, id::BB, EPSILON ) ) // inter
					<< " " << pose.torsion( TorsionID( ii, id::BB, ZETA ) ) // inter
					<< " " << pose.torsion( TorsionID( ii + 1, id::BB, ALPHA ) ) // inter
					// LATER: make this just report one term.
					<< " " << pose.energies().residue_total_energy( ii ) << std::endl;

			}
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}//main

