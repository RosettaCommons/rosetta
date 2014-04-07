// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
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

#include <protocols/simple_moves/ScoreMover.hh>

//#include <basic/options/util.hh>
//#include <basic/options/after_opts.hh>
//#include <core/pack/task/TaskFactory.hh>
//#include <core/pack/task/PackerTask.hh>
//#include <core/pack/pack_rotamers.hh>
//#include <core/scoring/TenANeighborGraph.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>
#include <core/pose/metrics/simple_calculators/InterfaceNeighborDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceSasaDefinitionCalculator.hh>
#include <core/pose/metrics/simple_calculators/InterfaceDeltaEnergeticsCalculator.hh>

#include <basic/options/option_macros.hh>

// Utility Headers
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/excn/Exceptions.hh>

// C++ Headers
#include <iostream>
#include <iomanip>
#include <map>

//Auto Headers
#include <core/import_pose/import_pose.hh>





using namespace core;
using namespace core::scoring;

OPT_KEY( Boolean, test_incorrect_type )

void register_metrics() {

	core::pose::metrics::PoseMetricCalculatorOP sasa_calculator = new core::pose::metrics::simple_calculators::SasaCalculatorLegacy;
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa", sasa_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int12_calculator =
		new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator((Size)1,(Size)2);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface_1_2", int12_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int23_calculator =
		new core::pose::metrics::simple_calculators::InterfaceNeighborDefinitionCalculator((Size)2,(Size)3);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface_2_3", int23_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int12_sasa_calculator =
		new core::pose::metrics::simple_calculators::InterfaceSasaDefinitionCalculator((Size)1,(Size)2);
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "sasa_interface_1_2", int12_sasa_calculator );

	core::pose::metrics::PoseMetricCalculatorOP int_delta_energy_calculator =
		new core::pose::metrics::simple_calculators::InterfaceDeltaEnergeticsCalculator( "interface_1_2" );
	core::pose::metrics::CalculatorFactory::Instance().register_calculator( "interface_delta_energies", int_delta_energy_calculator );

	return;
}

/// General testing code
int
main( int argc, char * argv [] )
{
    try {
	devel::init(argc, argv);

	register_metrics();

	NEW_OPT( test_incorrect_type, "test_incorrect_type", true );

	std::cout << "JK Starting hotspot test" << std::endl;

	pose::Pose pose;
	core::import_pose::pose_from_pdb( pose, "LJ111.pdb" );

	std::cout << "JK about to print total SASA directly as a string" << std::endl;
	std::string sasa_val = pose.print_metric("sasa","total_sasa");
	std::cout << "JK Total SASA is: " << sasa_val << std::endl;

	std::cout << "JK about to request total SASA in a MetricValue" << std::endl;
	basic::MetricValue<Real> mv_sasa;
	pose.metric("sasa","total_sasa",mv_sasa);
	std::cout << "JK Total SASA is: " << mv_sasa.value() << std::endl;

	std::cout << "JK about to print total SASA from a MetricValue" << std::endl;
	pose.metric("sasa","total_sasa",mv_sasa);
	std::cout << "JK Total SASA is: " << mv_sasa.print() << std::endl;

	std::cout << "JK about to request list_interface" << std::endl;
	basic::MetricValue< std::pair< utility::vector1<Size>, utility::vector1<Size> > > mv_li;
	pose.metric("interface_1_2","list_interface",mv_li);
	std::cout << "JK li first is: " << (mv_li.value()).first << std::endl;
	std::cout << "JK li second is: " << (mv_li.value()).second << std::endl;

	std::cout << "JK about to request delta SASA" << std::endl;
	basic::MetricValue< Real > mv_delta_sasa;
	pose.metric("sasa_interface_1_2","delta_sasa",mv_delta_sasa);
	std::cout << "JK delta_sasa is: " << mv_delta_sasa.value() << std::endl;

	// Setup for scoring/repacking
	std::cout << "JK About to score the pose" << std::endl;
	protocols::simple_moves::ScoreMover *score_mover = new protocols::simple_moves::ScoreMover;
	score_mover->apply( pose );
	Real score = pose.energies().total_energy();
	std::cout << "JK Total score is: " << score << std::endl;

	std::cout << "JK about to request interface_delta weighted_total" << std::endl;
	basic::MetricValue< Real > mv_delta_total;
	pose.metric("interface_delta_energies","weighted_total",mv_delta_total);
	std::cout << "JK delta_weighted_total is: " << mv_delta_total.value() << std::endl;

	std::cout << "JK about to request interface_delta fa_atr" << std::endl;
	basic::MetricValue< Real > mv_delta_atr;
	pose.metric("interface_delta_energies","fa_atr",mv_delta_atr);
	std::cout << "JK delta_atr is: " << mv_delta_atr.value() << std::endl;

	std::cout << "JK about to request interface_delta fa_rep" << std::endl;
	basic::MetricValue< Real > mv_delta_rep;
	pose.metric("interface_delta_energies","fa_rep",mv_delta_rep);
	std::cout << "JK delta_rep is: " << mv_delta_rep.value() << std::endl;

	std::cout << "JK about to request interface_delta fa_sol" << std::endl;
	basic::MetricValue< Real > mv_delta_sol;
	pose.metric("interface_delta_energies","fa_sol",mv_delta_sol);
	std::cout << "JK delta_sol is: " << mv_delta_sol.value() << std::endl;

	std::cout << "JK about to request interface_delta fa_pair" << std::endl;
	basic::MetricValue< Real > mv_delta_pair;
	pose.metric("interface_delta_energies","fa_pair",mv_delta_pair);
	std::cout << "JK delta_pair is: " << mv_delta_pair.value() << std::endl;

	//	std::cout << "JK about to request interface_delta fa_junk" << std::endl;
	//	basic::MetricValue< Real > mv_delta_junk;
	//	pose.metric("interface_delta_energies","fa_junk",mv_delta_junk);
	//	std::cout << "JK delta_junk is: " << mv_delta_junk.value() << std::endl;

	//	std::cout << "JK Moving phi angle and rescoring" << std::endl;
	//	pose.set_psi( 5, 4.5 );
	//	std::cout << "JK about to request total SASA after moving phi angle" << std::endl;
	//	pose.metric("sasa","total_sasa",mv_sasa);
	//	std::cout << "JK Total SASA is: " << mv_sasa.value() << std::endl;
	//	std::cout << "JK about to request total SASA again" << std::endl;
	//	pose.metric("sasa","total_sasa",mv_sasa);
	//	std::cout << "JK Total SASA is: " << mv_sasa.value() << std::endl;
	//	std::cout << "JK about to request total SASA again" << std::endl;
	//	pose.metric("sasa","total_sasa",mv_sasa);
	//	std::cout << "JK Total SASA is: " << mv_sasa.value() << std::endl;

	if ( basic::options::option[ basic::options::OptionKeys::test_incorrect_type ].user() ) {
		std::cout << "JK about to request total SASA as an int" << std::endl;
		basic::MetricValue<int> mv_int_sasa;
		pose.metric("sasa","total_sasa",mv_int_sasa);
		std::cout << "JK Total SASA is: " << mv_int_sasa.value() << std::endl;
	}

	std::cout << "JK about to request first_chain_first_resnum" << std::endl;
	basic::MetricValue<Size> mv_fcfr;
	pose.metric("interface_1_2","first_chain_first_resnum",mv_fcfr);
	std::cout << "JK mv_fcfr is: " << mv_fcfr.value() << std::endl;

	std::cout << "JK about to request first_chain_last_resnum" << std::endl;
	basic::MetricValue<Size> mv_fclr;
	pose.metric("interface_1_2","first_chain_last_resnum",mv_fclr);
	std::cout << "JK mv_fclr is: " << mv_fclr.value() << std::endl;

	std::cout << "JK about to request second_chain_first_resnum" << std::endl;
	basic::MetricValue<Size> mv_scfr;
	pose.metric("interface_1_2","second_chain_first_resnum",mv_scfr);
	std::cout << "JK mv_scfr is: " << mv_scfr.value() << std::endl;

	std::cout << "JK about to request second_chain_last_resnum" << std::endl;
	basic::MetricValue<Size> mv_sclr;
	pose.metric("interface_1_2","second_chain_last_resnum",mv_sclr);
	std::cout << "JK mv_sclr is: " << mv_sclr.value() << std::endl;

	std::cout << "JK Successfully finishing hotspot test" << std::endl;

    } catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
    }
	return 0;
}

