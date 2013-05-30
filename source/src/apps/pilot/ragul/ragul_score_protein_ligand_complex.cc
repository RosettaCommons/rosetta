// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @brief
/// @author jk
// Project Headers
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <basic/MetricValue.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/PackstatCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/pack_rotamers.hh>
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
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/rigid/RigidBodyMover.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
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
OPT_KEY( Real, cst_force_constant )
static basic::Tracer TR( "apps.pilot.ragul_darc_minimize.main" );
int main( int argc, char * argv [] ){
        NEW_OPT( cst_force_constant, "coordinate constraint force constant", 0.5 );
        devel::init(argc, argv);
        //setup scorefxn
        scoring::ScoreFunctionOP scorefxn(getScoreFunction());
        scoring::ScoreFunctionOP repack_scorefxn(getScoreFunction());
        //setup the bound pose
        pose::Pose bound_pose;
        std::string const input_pdb_name ( basic::options::start_file() );
        core::import_pose::pose_from_pdb( bound_pose, input_pdb_name );
  //create tag for output filename
  int dot_index1 = input_pdb_name.rfind(".", input_pdb_name.size());
  assert(dot_index1 != -1 && "No dot found in filename");
        std::string tag = input_pdb_name.substr(0,dot_index1);
        std::string init_pdb = "init_" + tag + ".pdb";
        std::string mini_pdb = "mini_" + tag + ".pdb";
        std::string unbo_pdb = "unbo_" + tag + ".pdb";
				(*scorefxn)(bound_pose);

        //setup the unbound pose
        core::pose::Pose unbound_pose = bound_pose;
        core::Real const unbound_dist = 80.;
				Size const rb_jump = 1; // use the first jump as the one between partners
        protocols::rigid::RigidBodyTransMover trans_mover( unbound_pose, rb_jump );
				trans_mover.trans_axis( trans_mover.trans_axis() );
				trans_mover.step_size(unbound_dist);
				trans_mover.apply( unbound_pose );
        (*scorefxn)(unbound_pose);
				unbound_pose.dump_pdb( unbo_pdb );
				std::string sasa_calc_name = "sasa";
        std::string hbond_calc_name = "hbond";
        std::string packstat_calc_name = "packstat";
        std::string burunsat_calc_name = "burunsat";
				//Register calculators
        core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculator;
        core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );
        core::pose::metrics::PoseMetricCalculatorOP hb_calc = new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
        core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );
        core::pose::metrics::PoseMetricCalculatorOP packstat_calc =     new protocols::toolbox::pose_metric_calculators::PackstatCalculator();
        core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calc );
        core::pose::metrics::PoseMetricCalculatorOP burunsat_calc = new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name);
        core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );
        // define containers for metrics for total complex
        basic::MetricValue<Real> tot_sasa_mval;
        basic::MetricValue<Size> tot_hb_mval;
        basic::MetricValue<Real> tot_packstat_mval;
        basic::MetricValue<Size> tot_unsat_mval;
        // calculate and store total metrics for bound and unbound poses
        core::Real bound_energy = 0.0, unbound_energy = 0.0, Interface_Energy = 0.0;
        core::Real bound_sasa = 0.0, unbound_sasa = 0.0, Total_BSA = 0.0;
        core::Size  bound_hb = 0,   unbound_hb = 0, Interface_HB = 0;
        core::Real bound_packstat = 0.0, unbound_packstat = 0.0, Total_packstats = 0.0;
        core::Size  bound_unsat = 0, unbound_unsat = 0, Interface_unsat = 0;
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
        bound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
        bound_unsat = tot_unsat_mval.value();
        unbound_pose.metric(burunsat_calc_name,"all_bur_unsat_polars", tot_unsat_mval);
        unbound_unsat = tot_unsat_mval.value();
        Interface_unsat = bound_unsat - unbound_unsat;
        std::cout << "Interface_Scores:" <<"    "<< input_pdb_name <<"  " << bound_energy <<"   " << Interface_Energy <<"       "<< Total_BSA <<"       "<< Interface_HB <<"    "<< Total_packstats <<" "<< Interface_unsat << std::endl;
        return 0;
}
