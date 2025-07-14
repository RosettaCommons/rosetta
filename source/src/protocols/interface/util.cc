// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody_design/util.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/interface/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/chains_util.hh>
#include <protocols/toolbox/CalcInterNeighborGroup.hh>

#include <utility/vector1.hh>
#include <utility/string_util.hh>
#include <basic/Tracer.hh>

#include <set>

static basic::Tracer TR("protocols.interface");

namespace protocols {
namespace interface {

utility::vector1<bool>
select_interface_residues(core::pose::Pose const & pose, std::string interface, core::Size interface_distance) {

	//Adapted from IAM
	using namespace core;
	using namespace utility;
	using namespace core::pose::metrics;
	using namespace protocols::toolbox;

	using one_group = std::set<core::Size>;
	typedef std::pair< one_group, one_group > group_pair;
	using group_set = utility::vector1<group_pair>;

	if ( interface.find('_') ==  std::string::npos ) {
		utility_exit_with_message("Unrecognized interface: "+interface+" must have side1 and side2, ex: LH_A or L_H to calculate interface residues");
	}



	std::set<core::Size> side1_chains;
	std::set<core::Size> side2_chains;
	vector1<std::string> sides = utility::string_split(interface, '_');
	//TR <<"Interface:"<< interface <<":"<<std::endl;
	//TR << "side1:" << sides[1] << ":" << std::endl;
	//TR << "side3:" << sides[2] << ":" << std::endl;

	// Because of the input format, we only support single chain letters.
	for ( core::Size i = 0; i <= sides[1].length() -1; ++i ) {
		//TR <<"C:"<<utility::to_string(sides[1][i]) << std::endl;
		side1_chains.insert(pose::get_chain_id_from_chain(std::string{sides[1][i]}, pose));
	}
	for ( core::Size i = 0; i <= sides[2].length() -1; ++i ) {
		//TR <<"C:"<<utility::to_string(sides[2][i]) << std::endl;
		side2_chains.insert(pose::get_chain_id_from_chain(std::string{sides[2][i]}, pose));
	}

	debug_assert (!side1_chains.empty()/*size() >= 1*/);
	debug_assert (!side2_chains.empty()/*size() >= 1*/);

	std::set<core::Size> side1_residues, side2_residues;

	for ( core::Size ii = 1; ii<= pose.size(); ++ii ) {

		if ( side1_chains.count( pose.chain( ii ) ) ) {
			side1_residues.insert( ii );
		} else if ( side2_chains.count(pose.chain(ii)) ) {
			side2_residues.insert( ii );
		}
	}

	//prep a vector of a pair of these residue sets for Steven's calculator
	std::pair< std::set<core::Size>, std::set<core::Size> > side_pairs;
	side_pairs.first = side1_residues;
	side_pairs.second = side2_residues;
	group_set chain_groups;
	chain_groups.push_back( side_pairs );

	/*
	std::string calc = "interface_res_calc" ;
	if ( ! CalculatorFactory::Instance().check_calculator_exists( calc ) ) {
	CalculatorFactory::Instance().register_calculator(  calc, utility::pointer::make_shared< InterGroupNeighborsCalculator >(chain_groups, interface_distance ) );
	}

	//std::set<core::Size> multichain_interface;
	basic::MetricValue< std::set<core::Size> > mv_interface_set;
	pose.metric(  calc, "neighbors", mv_interface_set);
	//set_interface_set( mv_interface_set.value() );

	*/

	CalcInterNeighborGroup calc = CalcInterNeighborGroup(chain_groups, interface_distance);
	calc.compute(pose);
	std::set<core::Size> interface_residues = calc.neighbors();

	vector1<bool> residues(pose.size(), false);
	std::set<core::Size>::const_iterator it, end ;
	for ( it = interface_residues.begin(), end = interface_residues.end(); it != end; ++it ) {
		residues[*it] = true;
	}

	return residues;

}


}
}
