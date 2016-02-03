// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/smlewis/AlanineScanner.cc
/// @brief Alanine scanning hack.  Ron Jacak has a much better one.  This is really only intended for use with alanine scanning results from the AnchoredDesign app.
/// @author Steven Lewis

// Unit Headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

#include <protocols/loops/LoopClass.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/MetricValue.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>


// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <basic/prof.hh>

#include <ObjexxFCL/string.functions.hh>


// C++ headers
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>


using basic::T;
using basic::Error;
using basic::Warning;

//replaces cout
static THREAD_LOCAL basic::Tracer TR( "apps.pilot.smlewis.AlanineScanner" );

//local options
namespace local{
basic::options::IntegerOptionKey const chain1("local::chain1");
basic::options::IntegerOptionKey const chain2("local::chain2");
}//local

//helper function
void register_metrics() {
	core::Size const chain1(basic::options::option[ local::chain1 ].value());
	core::Size const chain2(basic::options::option[ local::chain2 ].value());
	TR << "pick what chains your interface is between with -chain1 and -chain2, currently using "
		 << chain1 << " and " << chain2 << std::endl;
	TR << "you are likely to get bizarre behaviors if these chains aren't in the poses!" << std::endl;

	core::pose::metrics::PoseMetricCalculatorOP int_calculator =
		new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator(chain1, chain2);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface", int_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int_sasa_calculator =
		new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator(chain1, chain2);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa_interface", int_sasa_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int_delta_energy_calculator =
		new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( "interface" );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface_delta_energies", int_delta_energy_calculator );

	return;
}

int
main( int argc, char* argv[] )
{
	clock_t starttime = clock();
	basic::prof_reset();
	using basic::options::option;
	option.add( local::chain1, "chain 1 for interface definition" ).def(1);
	option.add( local::chain2, "chain 2 for interface definition" ).def(2);
	devel::init(argc, argv);
	register_metrics();

	core::pose::Pose pose; //fill it later...
	core::pose::Pose backup;
	utility::vector1< std::string > const pdbs( basic::options::start_files() );
	core::Size const numpdbs(pdbs.size());

	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	protocols::loops::Loops loops;
	loops.read_loop_file( basic::options::option[ basic::options::OptionKeys::loops::loop_file ].value() );
	utility::vector1< core::Size > loopspos;
	for( protocols::loops::Loops::const_iterator it=loops.begin(), it_end=loops.end(); it != it_end; ++it ){
		for( core::Size here(it->start()), stop(it->stop()); here <= stop; ++here){
			loopspos.push_back(here);
		}
	}
	TR << "loops: ";
	core::Size const loopsposn(loopspos.size());
	for( core::Size i(1); i <= loopsposn; ++i){
		TR << loopspos[i] << " ";
	}
	TR << std::endl;


	utility::vector1< core::Energy > pdb_start_energies(pdbs.size(), 0);
	utility::vector1< core::Energy > interface_start_energies(pdbs.size(), 0);
	//utility::vector1< core::Energy > diff_start_energies(pdbs.size(), 0);
	utility::vector1< utility::vector1< core::Energy > > mut_energies(pdbs.size(), (utility::vector1< core::Energy>(loopsposn, 0)));
	utility::vector1< utility::vector1< core::Energy > > interface_mut_energies(pdbs.size(), (utility::vector1< core::Energy>(loopsposn, 0)));
	utility::vector1< utility::vector1< core::Energy > > diff_mut_energies(pdbs.size(), (utility::vector1< core::Energy>(loopsposn, 0)));
	utility::vector1< utility::vector1< core::Energy > > diff_interface_mut_energies(pdbs.size(), (utility::vector1< core::Energy>(loopsposn, 0)));

	/*	TR << "pdb start energies size " << pdb_start_energies.size() << std::endl;;
	TR << "mut energies size " << mut_energies.size() << std::endl;
	TR << "mut energies 1 size " << mut_energies[1].size() << std::endl;
	TR << "mut energies 10 size " << mut_energies[10].size() << std::endl;*/

	basic::MetricValue< core::Real > mv_delta_sasa;
	basic::MetricValue< core::Real > mv_delta_total;

	core::chemical::ResidueTypeSetCAP typeset(core::chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::FA_STANDARD));
	core::chemical::ResidueType const & ALA(typeset->name_map("ALA"));

	for( core::Size pdb_n(1); pdb_n <= numpdbs; ++pdb_n){
		std::string const & pdb = pdbs[pdb_n];
		T(pdb) << "start " << pdb << " number " << pdb_n << " of " << numpdbs << std::endl;
		core::import_pose::pose_from_file(pose, pdb, core::import_pose::PDB_file);
		backup = pose;
		core::Energy const start_score((*score_fxn)(pose));
		T(pdb) << "pose starting score: " << start_score << std::endl;
		pdb_start_energies[pdb_n] = start_score;

		T(pdb) << "Info for this starting pose: " << std::endl;

		pose.metric("sasa_interface","delta_sasa",mv_delta_sasa);
		T(pdb) << "delta_sasa is: " << mv_delta_sasa.value() << std::endl;

		pose.metric("interface_delta_energies","weighted_total",mv_delta_total);
		core::Energy const interface_start_score(mv_delta_total.value());
		T(pdb) << "delta_weighted_total (fixed bb binding energy) is: " << interface_start_score << std::endl;
		interface_start_energies[pdb_n] = interface_start_score;

		T(pdb) << "distance between residues 99 and 12 (proxy for two binding modes): "
					 << pose.residue(12).atom("CA").xyz().distance(pose.residue(99).atom("CA").xyz()) << std::endl;

		for( core::Size i(1); i <= loopsposn; ++i){
			std::string TR(pdb + " " + ObjexxFCL::right_string_of(loopspos[i],2,'0'));
			T(TR) << "mutating position " << loopspos[i] << std::endl;
			core::pose::replace_pose_residue_copying_existing_coordinates(pose, loopspos[i], ALA);
			//std::string name(pdb + "_" + ObjexxFCL::right_string_of(loopspos[i],2,'0') + ".pdb");
			//pose.dump_scored_pdb(name, *score_fxn);

			core::Energy const mut_energy((*score_fxn)(pose));
			mut_energies[pdb_n][i] = mut_energy;
			T(TR) << "mutant score: " << mut_energy << std::endl;
			pose.metric("interface_delta_energies","weighted_total",mv_delta_total);
			core::Energy const interface_mut_energy(mv_delta_total.value());
			interface_mut_energies[pdb_n][i] = interface_mut_energy;
			T(TR) << "mutant delta_weighted_total (fixed bb binding energy) is: " << interface_mut_energy << std::endl;

			diff_mut_energies[pdb_n][i] = mut_energy - start_score;
			T(TR) << "mutant-start complex energy " << diff_mut_energies[pdb_n][i] << std::endl;
			diff_interface_mut_energies[pdb_n][i] = interface_mut_energy - interface_start_score;
			T(TR) << "mutant-start interface energy " << diff_interface_mut_energies[pdb_n][i] << std::endl;

			pose = backup;
		}
	}

	basic::prof_show();
	clock_t stoptime = clock();
	TR << "Whole run took " << ((double) stoptime - starttime )/CLOCKS_PER_SEC << " seconds" << std::endl;
	TR << "************************d**o**n**e**************************************" << std::endl;

	return 0;
}
