// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/moves/EvolutionaryDynamicsMover.cc
/// @brief Overwrites the boltzmann function in GenericMonteCarloMover to sample according to
/// evolutionary dynamics instead.
/// @author Christoffer Norn ( chnorn@gmail.com )


// Unit Headers
#include <protocols/simple_moves/EvolutionaryDynamicsMover.hh>
#include <protocols/simple_moves/EvolutionaryDynamicsMoverCreator.hh>
#include <basic/datacache/DataMapObj.hh>
#include <basic/datacache/DataMap.hh>

// C/C++ headers
#include <iostream>
#include <iterator>
#include <string>

// External headers
#include <boost/bind.hpp>
#include <boost/foreach.hpp>
#include <boost/function.hpp>
#include <algorithm>

// Utility headers
#include <basic/Tracer.hh>
#include <ObjexxFCL/format.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/ParsedProtocol.hh>

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <protocols/filters/Filter.hh>


// Package Headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverStatus.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>
#include <cmath>
#include <core/pose/util.hh>
#include <protocols/simple_filters/OperatorFilter.hh>
#include <protocols/filters/BasicFilters.hh>
//////////////////////////////////////////////////
// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
///////////////////////////////////////////////////
static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.EvolutionaryDynamicsMover" );
static THREAD_LOCAL basic::Tracer TR_energies( "protocols.simple_moves.EvolutionaryDynamicsMover.individual_energies" );

using namespace core;

namespace protocols {
namespace simple_moves {

using namespace ObjexxFCL::format;

std::string
EvolutionaryDynamicsMoverCreator::keyname() const
{
	return EvolutionaryDynamicsMoverCreator::mover_name();
}

protocols::moves::MoverOP
EvolutionaryDynamicsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new EvolutionaryDynamicsMover );
}

std::string
EvolutionaryDynamicsMoverCreator::mover_name()
{
	return "EvolutionaryDynamics";
}

/// @brief default constructor
EvolutionaryDynamicsMover::EvolutionaryDynamicsMover():
	Super(),
	population_size_( 1000000 ),
	disable_fitness_evaluation_( false )
{
	Super::initialize();
}


/// @brief destructor
EvolutionaryDynamicsMover::~EvolutionaryDynamicsMover()= default;

/// @brief clone this object
EvolutionaryDynamicsMover::MoverOP
EvolutionaryDynamicsMover::clone() const
{
	return EvolutionaryDynamicsMover::MoverOP( new EvolutionaryDynamicsMover( *this ) );
}

/// @brief create this type of object
EvolutionaryDynamicsMover::MoverOP
EvolutionaryDynamicsMover::fresh_instance() const
{
	return EvolutionaryDynamicsMover::MoverOP( new EvolutionaryDynamicsMover() );
}

bool
EvolutionaryDynamicsMover::boltzmann( Pose & pose, utility::vector1< core::Real > const & random_nums )
{
	++trial_counter_;
	TR.Debug <<"filters.size() "<<filters().size()<<std::endl;

	core::Real filter_val(0.0);

	runtime_assert( filters().size() == temperatures().size() );
	runtime_assert( filters().size() == sample_types().size() );
	runtime_assert( filters().size() == num_rejections().size() );
	runtime_assert( filters().size() == random_nums.size() );
	bool accepted( false );
	utility::vector1< Real > provisional_scores;
	provisional_scores.clear();
	for ( core::Size index( 1 ); index <= filters().size(); ++index ) {
		runtime_assert( random_nums[ index ] >= 0.0 );
		runtime_assert( random_nums[ index ] <= 1.0 );
		TR.Debug <<"Filter #"<<index<<std::endl;
		protocols::filters::FilterCOP filter( filters()[ index ] );
		Real const temp( temperatures()[ index ] );
		Real const flip( sample_types()[ index ] == "high" ? -1 : 1 );
		filter_val = filter->report_sm( pose );
		TR<<"Filter "<<index<<" reports "<<filter_val<<" ( best="<<lowest_scores()[index]<<"; last="<<last_accepted_scores()[index]<<" )"<<std::endl;

		provisional_scores.push_back( flip * filter_val );

		Real const boltz_factor = ( last_accepted_scores()[ index ] - provisional_scores[ index ] ) / temp;

		TR_energies.Debug <<"energy index, last_accepted_score, current_score "<<index<<" "<< last_accepted_scores()[ index ]<<" "<<provisional_scores[ index ]<<std::endl;
		TR.Debug <<"Current, best, boltz "<<provisional_scores[ index ]<<" "<<last_accepted_scores()[ index ]<<" "<<boltz_factor<<std::endl;


		bool reject_filter;

		// Here I'm using Crow and Kimuras fixation probability for diploid organisms: f = [1-exp(-2s)] / [1 - exp(-4Ns)],
		// where s is the selection coefficient and N is the population size.
		// The equation requires some thought to compute correctly for small selection coefficients due to
		// underflow problems. All problems can be avoided in the relevant range by using expm1 to avoid underflow.
		// If the selection coefficient is 0 f is undefined. As it I'll never be truely 0, I set it to
		// f = 1/(2N) for s=0.

#if __cplusplus>=201103L
		Real const selection_coefficient = ( provisional_scores[ index ] -  last_accepted_scores()[ index ]) / last_accepted_scores()[ index ];
		Real const x1 = -2 * selection_coefficient ;
		Real const numerator = -expm1(x1);
		Real const x2 = -4 * selection_coefficient * population_size_;
		Real const denominator = -expm1(x2);
		Real const fixation_probability = ( selection_coefficient == 0.0 ) ? 1 / ( 2 * population_size_ ) : numerator / denominator;


		Real const max_selection_coefficient = ( 1 - last_accepted_scores()[ index ]) / last_accepted_scores()[ index ];
		Real const max_x1 = -2 * max_selection_coefficient;
		Real const max_numerator = -expm1( max_x1 );
		Real const max_x2 = -4 * max_selection_coefficient * population_size_;
		Real const max_denominator = -expm1( max_x2 );
		Real const fix_p_max = max_numerator / max_denominator;

		Real const fix_p_normalizer = 1/fix_p_max;
#else
        utility_exit_with_message( "this code relies on expm1, which is not implemented in C98 currently (Aug 2016) used for building PyRosetta in windows");


        Real const fix_p_normalizer = 0.0;
        Real const fixation_probability = 0.0;
        Real const max_selection_coefficient = 0.0;



#endif

		if ( max_selection_coefficient == 0.0 ) {
			utility_exit_with_message( "The current sequence is perfectly fit! Some weird stuff might be going on. Check if the implementation is fit to handle this situation!" );
		}

		std::ostringstream curr_score;
		curr_score.precision(16);
		curr_score << last_accepted_scores()[ index ];
		pose::add_comment(pose, "curr_score", curr_score.str());

		std::ostringstream proposed_score;
		proposed_score.precision(16);
		proposed_score << provisional_scores[ index ];
		pose::add_comment(pose, "proposed_score", proposed_score.str());

		std::ostringstream fixp;
		fixp.precision(16);
		fixp << fixation_probability;
		pose::add_comment(pose, "fixation_probability", fixp.str());

		std::ostringstream fix_p_norm;
		fix_p_norm.precision(16);
		fix_p_norm << fix_p_normalizer;
		pose::add_comment(pose, "fix_p_norm", fix_p_norm.str());

		if ( fixation_probability * fix_p_normalizer > 1.0 ) {
			utility_exit_with_message( "Fixation probability greater than 1. This should never happen. There is a capping problem" );
		}

		reject_filter = ( random_nums[ index ] >= fixation_probability * fix_p_normalizer ); // we need to normalize here, as the accept rate otherwise will be too low (mutations will very rarely be accepted).

		if ( disable_fitness_evaluation_ || !reject_filter ) {
			accepted = true;
		} else {
			accepted = false;
			++num_rejections_[index];
			break;
		}

	}//for index

	if ( progress_file() != "" ) { //write progress data to file
		/// write a table to a progress file that has the following structure
		/// Trial# accept? filter_val score pose_comments protein_sequence

		std::ofstream data;
		data.open( progress_file().c_str(), std::ios::app );
		if ( !data.good() ) {
			utility_exit_with_message( "Unable to open EvolutionaryDynamics progress file for writing: " + progress_file() + "\n" );
		}

		std::string pose_sequence( "" );
		for ( core::Size chaini = 1 ; chaini <= pose.conformation().num_chains(); ++chaini ) {
			pose_sequence += pose.chain_sequence( chaini );
		}

		using namespace std;
		string stringed_comments("");
		map< string, string > const comments = core::pose::get_all_comments(pose);
		for ( auto const & comment : comments ) {
			stringed_comments += comment.first + ":" + comment.second + " ";
		}

		core::Real energy = scoring(pose);
		data.precision(12);
		data << trial_counter() << " " << accepted << " " << filter_val << " " << energy << " " << stringed_comments << " " <<pose_sequence<<'\n';

		data.flush();
	}

	if ( accepted ) {
		TR<<"Accept"<<std::endl;
		accept( pose, provisional_scores, MCA_accepted_thermally );
		return true;
	} else { // fi accepted
		TR<<"Reject"<<std::endl;
		set_mc_accepted( MCA_rejected );
		return false;
	}
} // boltzmann


/// @Brief
///comment
void
EvolutionaryDynamicsMover::apply( Pose & pose )
{
	Super::apply( pose );
}// apply


std::string
EvolutionaryDynamicsMover::get_name() const {
	return EvolutionaryDynamicsMoverCreator::mover_name();
}

/// @brief parse xml file
void
EvolutionaryDynamicsMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &filters, Movers_map const &movers, Pose const & pose )
{
	set_keep_filters( true );
	Super::parse_my_tag( tag, data, filters, movers, pose );
	population_size( tag->getOption< core::Size >( "population_size", 1000000 ) );
	disable_fitness_evaluation( tag->getOption< bool >( "disable_fitness_evaluation", false ) );


}


} // ns simple_moves
} // ns protocols
