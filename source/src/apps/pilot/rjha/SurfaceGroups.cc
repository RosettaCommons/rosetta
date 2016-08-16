// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/rjha/SurfaceGroups.cc
/// @brief defines groups of neighbors on the surface of a pose
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/PDBInfo.hh>
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

// C++ headers
//#include <iostream>
#include <string>
#include <sstream>

#include <core/import_pose/import_pose.hh>


using basic::T;
using basic::Error;
using basic::Warning;

//replaces cout
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.rjha" );

//local options
namespace local{
basic::options::IntegerOptionKey const surface_residue("surface_residue");
basic::options::BooleanOptionKey const use_pdb_numbering("use_pdb_numbering");
}//local

int
main( int argc, char* argv[] )
{
	try {

	using basic::options::option;
	option.add( local::surface_residue, "cutoff for surface residues ( <= # is surface)" ).def(18);
	option.add( local::use_pdb_numbering, "use pdb numbering instead of pose numbering").def(false);
	devel::init(argc, argv);

	//containers for other poses
	core::pose::Pose pose; //fill it later...
	utility::vector1< std::string > const pdbs( basic::options::start_files() );

	//iterate over all pdbs
	core::Size biggest_calc(0);
	std::string const calc_stem("nbr_dist_calc_");
	std::ostringstream calcname;
	for(utility::vector1< std::string >::const_iterator pdbname(pdbs.begin()), end(pdbs.end()); pdbname != end; ++pdbname){
		std::string const & pdb = *pdbname;
		TR << "we are on PDB: " << pdb << std::endl;
		core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);

		for(core::Size res(1); res <= pose.total_residue(); ++res){
			if(biggest_calc < res){ //create calculator
				calcname << calc_stem << res;
				core::pose::metrics::CalculatorFactory::Instance().register_calculator( calcname.str(), new protocols::toolbox::pose_metric_calculators::NeighborsByDistanceCalculator(res) );
				calcname.str("");
				++biggest_calc;
			}
		}

		basic::MetricValue< core::Size > num_n;

		typedef std::set< core::Size > SizeSet;
		SizeSet surface_res;

		//find surface residues
		TR << pdb << std::endl;
		for(core::Size i(1); i<=pose.total_residue(); ++i){
			calcname << calc_stem << i;
			pose.metric( calcname.str(), "num_neighbors", num_n);
			calcname.str("");

			TR << "residue  " << i << " num_neighbors " << num_n.value() << std::endl;
			TR << "pose2pdb " << pose.pdb_info()->pose2pdb(i) << "did this work" << std::endl;

			if( num_n.value() <= core::Size(option[local::surface_residue].value())) {
				TR << "adding " << i << " to surface set" << std::endl;
				surface_res.insert(i);
			}
		}

		//find groups
		TR << "sel resi ";
		for(SizeSet::const_iterator it(surface_res.begin()), end(surface_res.end()); it!=end; ++it){
			TR << pose.pdb_info()->pose2pdb(*it) << "+";
		}
		TR << std::endl;

		basic::MetricValue< SizeSet > neighbors;
		for(SizeSet::const_iterator it(surface_res.begin()), end(surface_res.end()); it!=end; ++it){
			calcname << calc_stem << *it;
			pose.metric( calcname.str(), "neighbors", neighbors );
			calcname.str("");

			if(option[local::use_pdb_numbering]) {
				TR << "residue " << pose.pdb_info()->pose2pdb(*it) << " surface neighbors ";
			}
			else {
				TR << "residue " << *it << " surface neighbors ";
			}

			for(SizeSet::const_iterator nit(neighbors.value().begin()), nend(neighbors.value().end()); nit!=nend; ++nit){

				if( surface_res.find(*nit) != end ){
					TR << " " << *nit;
				}
			}
			TR << std::endl;
		}

	}

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
