// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
/// @brief
/// @author Ragul Gowthaman

// Project Headers
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>

#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <utility/excn/Exceptions.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>

//Metrics
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/SaltBridgeCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/CatPiCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PiPiCalculator.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <protocols/analysis/PackStatMover.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// C++ Headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <cmath>
#include <map>

//Auto Headers
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>

using namespace core;
using namespace core::pose::datacache;
using namespace core::optimization;
using namespace core::pose::metrics;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::id;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace conformation;
using namespace protocols::simple_moves;
using namespace protocols::rigid;
using namespace core::chemical;
using namespace core::conformation;

OPT_KEY( Boolean, print_unbound )
OPT_KEY( Integer, jump_num )

int main( int argc, char * argv [] ){
	try{
		NEW_OPT( print_unbound, "print the mimnimized protein for debugging", false );
		NEW_OPT( jump_num, "chain break number to separate bound and unbound complex", 1 );

		devel::init(argc, argv);

		//setup scorefxn
		scoring::ScoreFunctionOP scorefxn = get_score_function();
		scoring::ScoreFunctionOP repack_scorefxn = get_score_function();

		//setup the bound pose
		pose::Pose bound_pose;
		std::string const input_pdb_name( basic::options::start_file() );
		core::import_pose::pose_from_file( bound_pose, input_pdb_name , core::import_pose::PDB_file);
		(*scorefxn)(bound_pose);

		//setup the unbound pose
		core::pose::Pose unbound_pose = bound_pose;
		core::Real const unbound_dist = 100.;
		//Size const rb_jump = bound_pose.num_jump(); // use the LAST jump as the one between partners
		//rb_jump hardcoded to move last two chains 3rd and 4th (alpha and beta chain in TCR)
		protocols::rigid::RigidBodyTransMover trans_mover1( unbound_pose, 2 );
		trans_mover1.trans_axis( trans_mover1.trans_axis() );
		trans_mover1.step_size(unbound_dist);
		trans_mover1.apply( unbound_pose );
		protocols::rigid::RigidBodyTransMover trans_mover2( unbound_pose, 3 );
		trans_mover2.trans_axis( trans_mover1.trans_axis() );
		trans_mover2.step_size(unbound_dist);
		trans_mover2.apply( unbound_pose );

		(*scorefxn)(unbound_pose);

		//create tag for output filename
		int pfounddir = input_pdb_name.find_last_of("/\\");
		int pfounddot = input_pdb_name.find_last_of(".");
		std::string tag = input_pdb_name.substr((pfounddir+1),(pfounddot-(pfounddir+1)));
		std::string unbo_pdb = "unbo_" + tag + ".pdb";

		if ( option[ print_unbound ] ) {
			unbound_pose.dump_pdb( unbo_pdb );
		}

		//Register calculators
		std::string sasa_calc_name = "sasa";
		std::string hbond_calc_name = "hbond";
		std::string packstat_calc_name = "packstat";
		std::string unsat_calc_name = "unsat";
		std::string sb_calc_name = "sb";
		std::string catpi_calc_name = "catpi";
		std::string pipi_calc_name = "pipi";

		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );

		core::pose::metrics::PoseMetricCalculatorOP hbonds_calculator( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hbonds_calculator );

		core::pose::metrics::PoseMetricCalculatorOP packstat_calculator( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calculator );

		core::pose::metrics::PoseMetricCalculatorOP unsat_calculator( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( unsat_calc_name, unsat_calculator );

		core::pose::metrics::PoseMetricCalculatorOP sb_calculator( new protocols::toolbox::pose_metric_calculators::SaltBridgeCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( sb_calc_name, sb_calculator );

		core::pose::metrics::PoseMetricCalculatorOP catpi_calculator( new protocols::toolbox::pose_metric_calculators::CatPiCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( catpi_calc_name, catpi_calculator );

		core::pose::metrics::PoseMetricCalculatorOP pipi_calculator( new protocols::toolbox::pose_metric_calculators::PiPiCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( pipi_calc_name, pipi_calculator );


		// define containers for metrics for total complex
		basic::MetricValue<core::Real> tot_sasa_mval;
		basic::MetricValue<core::Size> tot_hb_mval;
		basic::MetricValue<core::Real> tot_packstat_mval;
		basic::MetricValue<core::Size> tot_unsat_mval;
		basic::MetricValue<core::Size> tot_pipi_mval;
		basic::MetricValue<core::Size> tot_catpi_mval;
		basic::MetricValue<core::Size> tot_sb_mval;

		// calculate and store total metrics for bound and unbound poses
		core::Real bound_energy = 0.0, unbound_energy = 0.0, Interface_Energy = 0.0;
		core::Real bound_sasa = 0.0, unbound_sasa = 0.0, Total_BSA = 0.0;
		core::Size bound_hb = 0, unbound_hb = 0, Interface_HB = 0;
		core::Real bound_packstat = 0.0, unbound_packstat = 0.0, Total_packstats = 0.0;
		core::Size bound_unsat = 0, unbound_unsat = 0, Interface_unsat = 0;
		core::Size bound_pipi = 0, unbound_pipi = 0, Interface_pipi = 0;
		core::Size bound_catpi = 0, unbound_catpi = 0, Interface_catpi = 0;
		core::Size bound_sb = 0, unbound_sb = 0, Interface_sb = 0;

		//calculate interface Energy
		bound_energy = bound_pose.energies().total_energy();
		unbound_energy = unbound_pose.energies().total_energy();
		Interface_Energy = bound_energy - unbound_energy;

		//delta sasa calculation
		bound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
		bound_sasa = tot_sasa_mval.value();
		unbound_pose.metric(sasa_calc_name,"total_sasa",tot_sasa_mval);
		unbound_sasa = tot_sasa_mval.value();
		Total_BSA = unbound_sasa - bound_sasa;

		//interface hb calculation
		bound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
		bound_hb = tot_hb_mval.value();
		unbound_pose.metric(hbond_calc_name,"all_Hbonds", tot_hb_mval);
		unbound_hb = tot_hb_mval.value();
		Interface_HB = bound_hb - unbound_hb;

		//packstat calculation
		bound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		bound_packstat = tot_packstat_mval.value();
		unbound_pose.metric(packstat_calc_name,"total_packstat", tot_packstat_mval);
		unbound_packstat = tot_packstat_mval.value();
		Total_packstats = bound_packstat - unbound_packstat;

		//unsat polar calculation
		bound_pose.metric(unsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
		bound_unsat = tot_unsat_mval.value();
		unbound_pose.metric(unsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
		unbound_unsat = tot_unsat_mval.value();
		Interface_unsat = bound_unsat - unbound_unsat;

		//pi_pi calculation
		bound_pose.metric( pipi_calc_name, "pi_pi", tot_pipi_mval );
		bound_pipi = tot_pipi_mval.value();
		unbound_pose.metric( pipi_calc_name, "pi_pi", tot_pipi_mval );
		unbound_pipi = tot_pipi_mval.value();
		Interface_pipi = bound_pipi - unbound_pipi;

		//cat_pi calculation
		bound_pose.metric( catpi_calc_name, "cat_pi", tot_catpi_mval );
		bound_catpi = tot_catpi_mval.value();
		unbound_pose.metric( catpi_calc_name, "cat_pi", tot_catpi_mval );
		unbound_catpi = tot_catpi_mval.value();
		Interface_catpi = bound_catpi - unbound_catpi;

		//cat_pi calculation
		bound_pose.metric( sb_calc_name, "salt_bridge", tot_sb_mval  );
		bound_sb = tot_sb_mval.value();
		unbound_pose.metric( sb_calc_name, "salt_bridge",tot_sb_mval  );
		unbound_sb = tot_sb_mval.value();
		Interface_sb = bound_sb - unbound_sb;

		std::cout << "Interface_Scores:\t"<< tag <<"\t"<< input_pdb_name <<"\t" << bound_energy <<"\t"<< Interface_Energy <<"\t"<< Total_BSA <<"\t"<< Interface_HB <<"\t"<< Total_packstats <<"\t"<< Interface_unsat <<"\t"<<Interface_pipi<<"\t"<<Interface_catpi<<"\t"<<Interface_sb<< std::endl;

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
