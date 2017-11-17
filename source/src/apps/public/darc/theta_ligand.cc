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
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/after_opts.hh>

#include <protocols/simple_moves/ScoreMover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
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
#include <core/scoring/sasa.hh>

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

OPT_KEY( String, input_protein_ligand_complex )
OPT_KEY( String, input_protein )
OPT_KEY( String, input_ligand )

int main( int argc, char * argv [] ){
	try{

		NEW_OPT( input_protein_ligand_complex, "bound protein-ligand file name", "bound.pdb" );
		NEW_OPT( input_protein, "unbound protein file name", "unbound.pdb" );
		NEW_OPT( input_ligand, "unbound ligand file name", "ligand.pdb" );

		devel::init(argc, argv);

		std::string const bound_protein = option[ input_protein_ligand_complex ];
		std::string const unbound_protein = option[ input_protein ];
		std::string const inp_ligand = option[ input_ligand ];

		//setup scorefxn
		scoring::ScoreFunctionOP scorefxn(get_score_function());
		scoring::ScoreFunctionOP repack_scorefxn(get_score_function());

		//Register calculators
		std::string sasa_calc_name = "sasa";
		std::string hbond_calc_name = "hbond";
		std::string packstat_calc_name = "packstat";
		std::string burunsat_calc_name = "burunsat";
		core::pose::metrics::PoseMetricCalculatorOP sasa_calculator( new core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( sasa_calc_name, sasa_calculator );

		core::pose::metrics::PoseMetricCalculatorOP hb_calc( new protocols::toolbox::pose_metric_calculators::NumberHBondsCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( hbond_calc_name, hb_calc );

		core::pose::metrics::PoseMetricCalculatorOP packstat_calc( new protocols::toolbox::pose_metric_calculators::PackstatCalculator() );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( packstat_calc_name, packstat_calc );

		core::pose::metrics::PoseMetricCalculatorOP burunsat_calc( new protocols::toolbox::pose_metric_calculators::BuriedUnsatisfiedPolarsCalculator(sasa_calc_name, hbond_calc_name) );
		core::pose::metrics::CalculatorFactory::Instance().register_calculator( burunsat_calc_name, burunsat_calc );

		pose::Pose bound_protein_pose, unbound_protein_pose, ligand_pose;
		core::import_pose::pose_from_file( bound_protein_pose, bound_protein , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( unbound_protein_pose, unbound_protein , core::import_pose::PDB_file);
		core::import_pose::pose_from_file( ligand_pose, inp_ligand , core::import_pose::PDB_file);
		basic::MetricValue<Real> total_sasa_mval;
		core::Real ligand_pose_sasa = 0.0, bound_pose_sasa = 0.0, unbound_pose_sasa = 0.0, Total_pose_exposed_SASA = 0.0;
		ligand_pose.metric(sasa_calc_name,"total_sasa",total_sasa_mval);
		ligand_pose_sasa = total_sasa_mval.value();
		bound_protein_pose.metric(sasa_calc_name,"total_sasa",total_sasa_mval);
		bound_pose_sasa = total_sasa_mval.value();
		unbound_protein_pose.metric(sasa_calc_name,"total_sasa",total_sasa_mval);
		unbound_pose_sasa = total_sasa_mval.value();
		Total_pose_exposed_SASA = 1 - ((((unbound_pose_sasa + ligand_pose_sasa) - bound_pose_sasa)/2)/ligand_pose_sasa);
		core::Real theta_ligand = Total_pose_exposed_SASA;

		//hydrophobic and polar SASA
		/*
		utility::vector1< core::Real > complex_rsd_sasa( bound_protein_pose.size(), 0.0 );
		utility::vector1< core::Real > separated_rsd_sasa( unbound_protein_pose.size(), 0.0 );
		utility::vector1< core::Real > complex_rsd_hsasa( bound_protein_pose.size(), 0.0 ); // hydrophobic SASA only
		utility::vector1< core::Real > separated_rsd_hsasa( unbound_protein_pose.size(), 0.0 ); // hydrophobic SASA only

		core::Real probe_radius = basic::options::option[basic::options::OptionKeys::pose_metrics::sasa_calculator_probe_radius];

		core::scoring::calc_per_res_hydrophobic_sasa( bound_protein_pose, complex_rsd_sasa, complex_rsd_hsasa, probe_radius, false );
		core::scoring::calc_per_res_hydrophobic_sasa( unbound_protein_pose, separated_rsd_sasa, separated_rsd_hsasa, probe_radius, false );

		core::Size lig_res_num = 0;
		for ( int j = 1, resnum = bound_protein_pose.size(); j <= resnum; ++j ) {
		if (!bound_protein_pose.residue(j).is_protein()){
		lig_res_num = j;
		break;
		}
		}
		if (lig_res_num == 0){
		std::cout<<"No ligand given in bound protein PDB structure.  Cannot identify interface."<<std::endl;
		exit (1);
		}

		core::Real complex_sasa = 0.0;
		core::Real complex_polar_sasa = 0.0;
		core::Real complex_hydrophobic_sasa = 0.0;
		core::Real separated_sasa = 0.0;
		core::Real separated_polar_sasa = 0.0;
		core::Real separated_hydrophobic_sasa = 0.0;

		for ( core::uint j = 1, resnum = bound_protein_pose.size(); j <= resnum; ++j ) {
		if ( j == lig_res_num ) continue;
		complex_sasa += complex_rsd_sasa[j];
		complex_hydrophobic_sasa += complex_rsd_hsasa[j];
		separated_sasa += separated_rsd_sasa[j];
		separated_hydrophobic_sasa += separated_rsd_hsasa[j];
		}
		complex_polar_sasa = complex_sasa - complex_hydrophobic_sasa;
		separated_polar_sasa = separated_sasa - separated_hydrophobic_sasa;

		core::Real interface_sasa = std::abs(complex_sasa - separated_sasa);
		core::Real interface_hydrophobic_sasa = (std::abs(complex_hydrophobic_sasa - separated_hydrophobic_sasa))/interface_sasa;
		core::Real interface_polar_sasa = (std::abs(complex_polar_sasa - separated_polar_sasa))/interface_sasa;
		std::cout << "Scores: interface_hydrophobic_sasa "<<bound_protein<<" "<<  interface_hydrophobic_sasa <<std::endl;
		std::cout << "Scores: interface_polar_sasa "<<bound_protein<<" "<<  interface_polar_sasa <<std::endl;
		*/
		std::cout << "Scores: theta_ligand "<<bound_protein<<"\t"<<  theta_ligand <<std::endl;


	} catch (utility::excn::Exception const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
