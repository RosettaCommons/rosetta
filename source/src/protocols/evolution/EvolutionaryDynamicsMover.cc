// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/evolution/EvolutionaryDynamicsMover.cc
/// @brief Overwrites the boltzmann function in GenericMonteCarloMover to sample according to
/// evolutionary dynamics instead.
/// @author Christoffer Norn ( chnorn@gmail.com )


// Unit Headers
#include <protocols/evolution/EvolutionaryDynamicsMover.hh>
#include <protocols/evolution/EvolutionaryDynamicsMoverCreator.hh>
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
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>
///////////////////////////////////////////////////
static THREAD_LOCAL basic::Tracer TR( "protocols.evolution.EvolutionaryDynamicsMover" );
static THREAD_LOCAL basic::Tracer TR_energies( "protocols.evolution.EvolutionaryDynamicsMover.individual_energies" );

using namespace core;
using namespace protocols::moves;

namespace protocols {
namespace evolution {

using namespace ObjexxFCL::format;

// XRW TEMP std::string
// XRW TEMP EvolutionaryDynamicsMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return EvolutionaryDynamicsMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP EvolutionaryDynamicsMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new EvolutionaryDynamicsMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP EvolutionaryDynamicsMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "EvolutionaryDynamics";
// XRW TEMP }

/// @brief default constructor
EvolutionaryDynamicsMover::EvolutionaryDynamicsMover():
	Super(),
	population_size_( 1000000 ),
	disable_fitness_evaluation_( false ),
	n_nucleotide_mut_trials_corrected_( 0.0 ),
	total_trials_( 0 )
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
		Real const flip( sample_types()[ index ] == "high" ? -1 : 1 );
		filter_val = filter->report_sm( pose );
		TR<<"Filter "<<index<<" reports "<<filter_val<<" ( best="<<lowest_scores()[index]<<"; last="<<last_accepted_scores()[index]<<" )"<<std::endl;

		std::map< std::string, std::string > comments = get_all_comments( pose );
		std::string key = "stop_codon";
		bool is_stop_codon = ( "1" == comments[ key ] );

		if ( is_stop_codon ) {
			provisional_scores.push_back( 0.0 );
		} else {
			provisional_scores.push_back( flip * filter_val );
		}
		TR_energies.Debug <<"energy index, last_accepted_score, current_score "<<index<<" "<< last_accepted_scores()[ index ]<<" "<<provisional_scores[ index ]<<std::endl;

		bool reject_filter;

		// Here I'm using Crow and Kimuras fixation probability for diploid organisms: f = [1-exp(-2s)] / [1 - exp(-4Ns)],
		// where s is the selection coefficient and N is the population size.
		// The equation requires some thought to compute correctly for small selection coefficients due to
		// underflow problems. All problems can be avoided in the relevant range by using expm1 to avoid underflow.
		// If the selection coefficient is 0 f is undefined. As it I'll never be truely 0, I set it to
		// f = 1/(2N) for s=0.

#if __cplusplus>=201103L
		Real const selection_coefficient = ( provisional_scores[ index ] -  last_accepted_scores()[ index ]) / last_accepted_scores()[ index ];
		Real const x1 = -2 * selection_coefficient;
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

		Real fix_p_normalizer;
		if ( total_trials_ != 0 ) { // branch length and mutation_rate is provided by the user and is respected in the below
			core::Real remaining_trials = total_trials_ - n_nucleotide_mut_trials_corrected_;
			if ( 1/fix_p_max < remaining_trials ) { // we cant do more trials than we have remaining.
				fix_p_normalizer = 1/fix_p_max;
			} else {
				fix_p_normalizer = remaining_trials;
				// Here we need to communicate to the MC mover that it is time to stop the trajectory
				set_stop_sampling( true ); // This will stop the MC trajectory
			}
		} else { // We are running a trajectory with a certain number of target accepts
			fix_p_normalizer = 1/fix_p_max;
		}
		n_nucleotide_mut_trials_corrected_ += fix_p_normalizer;

		std::ostringstream ntrials_out;
		ntrials_out.precision(16);
		ntrials_out << n_nucleotide_mut_trials_corrected_;
		pose::add_comment(pose, "elapsed_nt_trials_corrected", ntrials_out.str());

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
		pose::add_comment(pose, "fitness_current_seq", curr_score.str());

		std::ostringstream proposed_score;
		proposed_score.precision(16);
		proposed_score << provisional_scores[ index ];
		pose::add_comment(pose, "fitness_proposed_seq", proposed_score.str());

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

		data.precision(12);
		data << "trial:" << trial_counter() << " accept:" << accepted << " " << stringed_comments <<'\n';

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


// XRW TEMP std::string
// XRW TEMP EvolutionaryDynamicsMover::get_name() const {
// XRW TEMP  return EvolutionaryDynamicsMover::mover_name();
// XRW TEMP }

/// @brief parse xml file
void
EvolutionaryDynamicsMover::parse_my_tag( TagCOP const tag, basic::datacache::DataMap & data, Filters_map const &filters, Movers_map const &movers, Pose const & pose )
{
	set_keep_filters( true );

	Super::parse_my_tag( tag, data, filters, movers, pose );
	population_size( tag->getOption< core::Size >( "population_size", 1000000 ) );
	disable_fitness_evaluation( tag->getOption< bool >( "disable_fitness_evaluation", false ) );
	if ( tag->hasOption("branch_length") && tag->hasOption("mutation_rate") ) {
		branch_length( tag->getOption< core::Real >( "branch_length", 100000.0 ) );
		mutation_rate( tag->getOption< core::Real >( "mutation_rate", 0.001 ) );
		total_trials( branch_length() / mutation_rate() ); // This is to total nt substitution trials we need to try to simulate the expected branch length
		TR << "Branch length " << branch_length() << " mutation rate " << mutation_rate() << " Total nt substitution attempts: " << total_trials_  << std::endl;;
		set_maxtrials( total_trials_ ); // here we overwrite whatever was given for trials for the MC mover.
		set_max_accepted_trials( total_trials_ );
		TR << "If sequences are near perfectly fit you can expect " << total_trials_ / ( Real (2*population_size()) ) << " trials." << std::endl;
	}
}

std::string EvolutionaryDynamicsMover::get_name() const {
	return mover_name();
}

std::string EvolutionaryDynamicsMover::mover_name() {
	return "EvolutionaryDynamics";
}

void EvolutionaryDynamicsMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::attribute_w_default("population_size", xsct_positive_integer, "Set population size. Default = 1000000", "1000000")
		+ XMLSchemaAttribute::attribute_w_default("disable_fitness_evaluation", xsct_rosetta_bool, "Disable fitness evaluation. Default = false", "false")
		+ XMLSchemaAttribute::attribute_w_default("branch_length", xsct_real, "Set branch_length. ", "1000000")
		+ XMLSchemaAttribute::attribute_w_default("mutation_rate", xsct_real, "Set mutation_rate. ", "0.001");

	utility::tag::XMLSchemaComplexTypeGeneratorOP ct_gen = define_composition_schema( xsd );
	ct_gen->element_name( mover_name() )
		.description( "Overwrites the boltzmann function in GenericMonteCarloMover to sample according to evolutionary dynamics instead." )
		.add_attributes( attlist )
		.write_complex_type_to_schema( xsd );
}

std::string EvolutionaryDynamicsMoverCreator::keyname() const {
	return EvolutionaryDynamicsMover::mover_name();
}

protocols::moves::MoverOP
EvolutionaryDynamicsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new EvolutionaryDynamicsMover );
}

void EvolutionaryDynamicsMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	EvolutionaryDynamicsMover::provide_xml_schema( xsd );
	//GenericMonteCarloMover::define_composition_schema( xsd );
}



} // ns evolution
} // ns protocols
