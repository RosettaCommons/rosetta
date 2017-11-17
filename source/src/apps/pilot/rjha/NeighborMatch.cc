// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rjha/NeighborMatch.cc
/// @brief Finds residues neighboring a given residue across an interface.  Intended for creating posfiles suitable for RosettaMatch
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <protocols/toolbox/pose_metric_calculators/NeighborsByDistanceCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using basic::Error;
using basic::Warning;

//replaces cout
static basic::Tracer TR( "apps.pilot.rjha" );

//local options
namespace local {
basic::options::IntegerOptionKey const nbr_residue("local::nbr_residue");
basic::options::StringOptionKey const master("local::master");
}//local

int
main( int argc, char* argv[] )
{
	try {
		using basic::options::option;
		option.add( local::nbr_residue, "residue whose neighbors we care about" ).def(62);
		option.add( local::master, "master pose (other poses matched against this" ).def("ubc12.pdb");
		devel::init(argc, argv);

		core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

		//pose containing res
		core::pose::Pose master;
		core::import_pose::pose_from_file(master, basic::options::option[local::master].value(), core::import_pose::PDB_file);
		core::Size const mastersize(master.size());

		//containers for other poses
		core::pose::Pose pose; //fill it later...
		core::pose::Pose combined;
		utility::vector1< std::string > const pdbs( basic::options::start_files() );

		//create calculator
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( "nbr_calc", new protocols::toolbox::PoseMetricCalculators::NeighborsByDistanceCalculator(basic::options::option[local::nbr_residue].value()) );

		//iterate over all pdbs
		for ( utility::vector1< std::string >::const_iterator pdbname(pdbs.begin()), end(pdbs.end()); pdbname != end; ++pdbname ) {
			std::string const & pdb = *pdbname;
			TR << "we are on PDB: " << pdb << std::endl;
			core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);

			combined = master;
			combined.append_residue_by_jump(pose.residue(1), 1);
			for ( core::Size i=2; i<=pose.size(); ++i ) {
				combined.append_residue_by_bond(pose.residue(i));
			}
			//patching in disulfide compatibility
			combined.conformation().detect_disulfides();
			combined.conformation().detect_bonds();
			combined.conformation().detect_pseudobonds();
			for ( core::Size i=1; i<=combined.size(); ++i ) {
				combined.conformation().update_polymeric_connection(i);
			}

			(*score_fxn)(combined);

			//combined.dump_scored_pdb("test.pdb", *score_fxn);

			//query calculator
			typedef std::set< core::Size > SetSize;
			basic::MetricValue< SetSize > neighbors;
			combined.metric( "nbr_calc", "neighbors", neighbors );

			//print results
			TR << pdb << " ";
			for ( SetSize::const_iterator it(neighbors.value().begin()), end(neighbors.value().end()); it != end; ++it ) {
				if ( *it > mastersize ) {
					TR << (*it)-mastersize << " ";
				}
			}
			TR << std::endl;
		}

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
