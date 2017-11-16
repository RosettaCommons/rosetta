// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief Computes a ligand's buried polar unsatisfied atoms
///
/// @param[in] -s <PDBFIL>, where <PDBFIL> is the path to the PDB file
///  containing the pose
/// @param[in] -extra_res_fa <PARFIL>, where <PARFIL> is the path to the params
///  file describing the ligand
///
/// @details It is assumed that the ligand is the last residue of the pose
///
/// @author Andrea Bazzoli (bazzoli@ku.edu)

#include <core/pose/metrics/simple_calculators/SasaCalculator.hh>
#include <protocols/toolbox/pose_metric_calculators/NumberHBondsCalculator.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <basic/MetricValue.hh>
#include <protocols/toolbox/pose_metric_calculators/BuriedUnsatisfiedPolarsCalculator.hh>
#include <basic/Tracer.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/after_opts.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <set>
#include <string>

static basic::Tracer TR( "apps.pilot.list_accpt_pos" );

using core::Size;
using core::Real;
using core::pose::Pose;
using std::string;


////////////////////////////////////////////////////////////////////////////////
//                                    MAIN                                    //
////////////////////////////////////////////////////////////////////////////////


int main( int argc, char * argv [] )
{
	devel::init(argc, argv);

	// create pose from pdb
	core::pose::Pose ps;
	std::string const input_pdb_name( basic::options::start_file() );
	core::import_pose::pose_from_file( ps, input_pdb_name , core::import_pose::PDB_file);

	// create score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	(*scorefxn)(ps);

	// register calculators needed to compute buried polar unsatisfied
	string sasa_calc_name = "sasa";
	core::pose::metrics::PoseMetricCalculatorOP sasa_calculator =
		new core::pose::metrics::simple_calculators::SasaCalculator;
	core::pose::metrics::CalculatorFactory::Instance().register_calculator(
		sasa_calc_name, sasa_calculator );

	string hbond_calc_name = "hbond";
	core::pose::metrics::PoseMetricCalculatorOP hb_calc =
		new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator();
	core::pose::metrics::CalculatorFactory::Instance().register_calculator(
		hbond_calc_name, hb_calc );

	std::set<Size> lig;
	lig.insert(ps.size());

	string burunsat_calc_name = "burunsat";
	core::pose::metrics::PoseMetricCalculatorOP burunsat_calc =
		new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(
		sasa_calc_name, hbond_calc_name, lig);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator(
		burunsat_calc_name, burunsat_calc );

	basic::MetricValue<Size> tot_unsat_mval;

	ps.metric(burunsat_calc_name,"special_region_bur_unsat_polars", tot_unsat_mval);
	Size nbpuns = tot_unsat_mval.value();

	std::cout << "number of unsatisfied polar atoms: " << nbpuns << std::endl;
}
