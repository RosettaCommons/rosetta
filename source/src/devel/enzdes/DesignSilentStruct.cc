// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/enzdes/DesignSilentStruct.cc
///
/// @brief protein silent-file structures for designs, also contains functionality to query
/// @brief a given pose and extract some more data to print
/// @author Florian Richter


// C++ Headers
#include <iostream>
#include <vector>
#include <map>

#include <core/chemical/ResidueConnection.hh>
#include <core/pose/Pose.hh>

// mini headers

#include <basic/Tracer.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <devel/enzdes/DesignSilentStruct.hh>
#include <basic/MetricValue.hh>

//Auto Headers
#include <core/io/silent/ProteinSilentStruct.tmpl.hh>
#include <utility/vector1.hh>

//#include <devel/PoseMetricCalculators/InterfaceDeltaEnergeticsCalculator.hh> //needed for type checking

using namespace core;
using namespace core::io::silent;
using namespace ObjexxFCL::format;

namespace devel {
namespace enzdes {

static THREAD_LOCAL basic::Tracer tr( "devel.enzdes.DesignSilentStruct" );


DesignSilentStruct::DesignSilentStruct(
	core::io::silent::SilentFileOptions const & opts,
	core::pose::Pose const & pose,
	std::string tag,
	bool const add_in,
	bool const onlyadd_in
) :
	ProteinSilentStruct( opts )
{
	decoy_tag( tag );
	sequence( pose.sequence() );
	nres( pose.size() );

	print_additional_ = add_in;
	print_only_additional_ = onlyadd_in;
}


DesignSilentStruct::DesignSilentStruct(
	core::io::silent::SilentFileOptions const & opts,
	pose::Pose const & pose,
	std::string tag, // = "empty_tag",
	utility::vector1<Size> const & spec_res_in,
	utility::vector1< std::string > const & rel_score_in,
	bool const add_in,
	bool const onlyadd_in
) :
	DesignSilentStruct( opts, pose, tag, spec_res_in, rel_score_in,
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > >(), /* empty calculator map */
	add_in, onlyadd_in)
{
}


DesignSilentStruct::DesignSilentStruct(
	core::io::silent::SilentFileOptions const & opts,
	pose::Pose const & pose,
	std::string tag,
	utility::vector1<Size> const & spec_res_in,
	utility::vector1< std::string > const & rel_score_in,
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators,
	bool const add_in,
	bool const onlyadd_in
) :
	ProteinSilentStruct( opts )
{

	if ( (add_in == false) && (onlyadd_in == true) ) {
		std::cerr << "do you want a silent structure or not?" << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	decoy_tag( tag );
	sequence( pose.sequence() );
	nres( pose.size() );

	print_additional_ = add_in;
	print_only_additional_ = onlyadd_in;

	if ( (add_in == true) || (onlyadd_in == true) ) {
		calculate_additional_info( pose, spec_res_in, rel_score_in, calculators );
	}

}

void
DesignSilentStruct::print_header( std::ostream& out ) const {

	out << LJ( 50, "#Design_name");
	for ( auto const & it : additional_info_silent_energy_ ) {
		out << " " << A( it.width(), it.name() );
	}

	out << "\n";

}

void
DesignSilentStruct::print_scores( std::ostream& out ) const
{

	if ( !print_only_additional_ ) {
		parent::print_scores( out );
	}

	if ( print_additional_ ) {
		print_additional_info( out );
	}

}

void
DesignSilentStruct::add_to_additional_silent_energies(
	utility::vector1< core::io::silent::SilentEnergy > const & silent_Es)
{

	for ( auto const & silent_E : silent_Es ) {
		additional_info_silent_energy_.push_back( silent_E );
	}

}

void
DesignSilentStruct::print_additional_info( std::ostream& out) const
{
	int precision = 2; //number of digits after decimal place

	out << LJ( 50, decoy_tag() );

	for ( auto const & it : additional_info_silent_energy_ ) {
		out << " " << F( it.width(), precision, it.value() );
	}

	out << "\n";

}

void
DesignSilentStruct::calculate_additional_info(
	pose::Pose const & pose,
	utility::vector1<Size> const & special_res,
	utility::vector1< std::string > const & score_terms,
	std::map< Size, utility::vector1< std::pair< std::string, std::string > > > const & calculators )
{

	bool separate_out_constraints = false;
	auto cstfind = find( score_terms.begin(), score_terms.end(),"all_cst");
	if ( cstfind != score_terms.end() ) separate_out_constraints = true;

	//first write out the relevant score terms for the pose total
	for ( auto const & score_term : score_terms ) {
		std::string sco_name = score_term;
		int width = std::max( 10, (int) sco_name.length() + 3 );

		SilentEnergy new_se;
		if ( score_term == "all_cst" ) {
			new_se = SilentEnergy ( sco_name, sum_constraint_terms(pose, -1 ), 1, width);
		} else if ( separate_out_constraints && ( score_term == "total_score" ) ) {
			core::Real desired_value = pose.energies().total_energies()[ core::scoring::score_type_from_name( score_term ) ] - sum_constraint_terms(pose, -1 );
			new_se = SilentEnergy(sco_name, desired_value, 1, width);
		} else {
			new_se = SilentEnergy ( sco_name, pose.energies().total_energies()[ core::scoring::score_type_from_name( score_term ) ] * pose.energies().weights()[ core::scoring::score_type_from_name( score_term ) ], 1 ,width);
		}
		additional_info_silent_energy_.push_back( new_se );

	}

	//pose metric calculators for pose total
	auto totcalc_it = calculators.find( 0 );
	if ( totcalc_it != calculators.end() ) {

		utility::vector1< std::pair< std::string, std::string > > const & tot_calculators = totcalc_it->second;
		for ( auto const & tot_calculator : tot_calculators ) {

			std::string calc_name = "tot_" + tot_calculator.first;
			int width = std::max( 10, (int) calc_name.length() + 3 );

			core::Real calc_value;

			//following lines fairly hacky, but don't know a better solution at the moment
			if ( tot_calculator.first == "hbond_pm" || tot_calculator.first == "burunsat_pm" || tot_calculator.first == "NLconts_pm" ) {
				basic::MetricValue< core::Size > mval_size;
				pose.metric( tot_calculator.first, tot_calculator.second, mval_size );
				calc_value = mval_size.value();
			} else {
				basic::MetricValue< core::Real > mval_real;
				pose.metric( tot_calculator.first, tot_calculator.second, mval_real );
				calc_value = mval_real.value();
			}

			SilentEnergy new_se( calc_name, calc_value, 1, width);
			additional_info_silent_energy_.push_back( new_se );
		}

	}

	//done with pose totals

	//then write out the relevant scoreterms (and potentially pose metrics) for each of the special residues
	Size spec_res_counter(0);
	for ( core::Size special_re : special_res ) {

		spec_res_counter++;
		//for convenience, the sequence number of the residue will be written out
		std::stringstream temp;
		temp << spec_res_counter;
		std::string spec_res_name = "SR_" + temp.str();
		//std::cerr << "name for res " << *res_it << " is " << spec_res_name ;
		SilentEnergy res_name(spec_res_name, special_re, 1, 10);
		additional_info_silent_energy_.push_back( res_name );

		for ( auto const & score_term : score_terms ) {

			std::string sco_name = spec_res_name + "_" +  score_term ;
			int width = std::max( 10, (int) sco_name.length() + 3 );

			SilentEnergy new_se;
			if ( score_term == "all_cst" ) {
				new_se = SilentEnergy ( sco_name, sum_constraint_terms(pose, special_re), 1, width);
			} else if ( separate_out_constraints && ( score_term == "total_score" ) ) {
				core::Real desired_value = pose.energies().residue_total_energies( special_re )[ core::scoring::score_type_from_name( score_term ) ] - sum_constraint_terms(pose, special_re );
				new_se = SilentEnergy(sco_name, desired_value, 1, width);
			} else {
				new_se = SilentEnergy ( sco_name, pose.energies().residue_total_energies( special_re )[ core::scoring::score_type_from_name( score_term ) ] * pose.energies().weights()[ core::scoring::score_type_from_name( score_term ) ], 1 ,width);
			}

			additional_info_silent_energy_.push_back( new_se );
		}//loop over relevant scoreterms


		//if there are calculators that need to be evaluated for this residue, let's do that now
		//note: section still under development, right now only calculators that return reals are supported
		auto res_calc_it = calculators.find( special_re );
		if ( res_calc_it != calculators.end() ) {

			utility::vector1< std::pair< std::string, std::string > > calculators_this_res = res_calc_it->second;
			for ( auto & calculators_this_re : calculators_this_res ) {

				std::string res_calc_name = spec_res_name + "_" + calculators_this_re.first;
				int width = std::max( 10, (int) res_calc_name.length() + 3 );

				core::Real calc_value;

				basic::MetricValue< core::Real > mval_real;
				basic::MetricValue< utility::vector1< core::Size > >mval_sizevec;
				basic::MetricValue< utility::vector1< core::Real > >mval_realvec;

				//following lines fairly hacky, but don't know a better solution at the moment
				if ( ( calculators_this_re.first == "hbond_pm") || ( calculators_this_re.first == "burunsat_pm") ) {
					pose.metric( calculators_this_re.first, calculators_this_re.second, mval_sizevec );
					calc_value = mval_sizevec.value()[special_re];
				} else if ( (calculators_this_re.first == "pstat_pm") || (calculators_this_re.first == "nlpstat_pm" ) ) {
					pose.metric( calculators_this_re.first, calculators_this_re.second, mval_realvec );
					calc_value = mval_realvec.value()[special_re];
				} else {
					pose.metric( calculators_this_re.first, calculators_this_re.second, mval_real );
					calc_value = mval_real.value();
				}
				//std::cerr << " hehe, just executed pose metric for " << calc_it->first << " calculator   ";


				//special case: if this is an interface calculation, we do not want to include constraint terms
				//hacky at the moment, haven't figured out yet how to do this really clean
				if ( separate_out_constraints && ( calculators_this_re.first == "interf_E_1_2") ) {
					calc_value = calc_value - ( 2 * sum_constraint_terms(pose, special_re) );
				}
				SilentEnergy new_se( res_calc_name, calc_value, 1, width);
				additional_info_silent_energy_.push_back( new_se );

			} // for calculators this res

		}// if calculators for this res
		//else std::cerr << "resi " << *res_it << " has no calcs." << std::endl;

	}//loop over special residues

}

core::Real
DesignSilentStruct::sum_constraint_terms( pose::Pose const & pose, int which_res){

	using namespace core::scoring;

	EnergyMap all_weights = pose.energies().weights();
	EnergyMap scores;

	if ( which_res == -1 ) { //means we want to know stuff for the whole pose
		scores = pose.energies().total_energies();
	} else { scores = pose.energies().residue_total_energies( which_res ); }

	return scores[ coordinate_constraint ] * all_weights[ coordinate_constraint ] + scores[atom_pair_constraint] * all_weights[ atom_pair_constraint] +
		scores[ angle_constraint ] * all_weights[ angle_constraint ] + scores[ dihedral_constraint ] * all_weights[ dihedral_constraint ];

}


} //namespace enzdes
} //namespace devel
