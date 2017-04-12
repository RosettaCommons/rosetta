// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/sidechain_moves/SidechainMCMover.cc
/// @brief implementation of SidechainMCMover class and functions
/// @author Colin A. Smith (colin.smith@ucsf.edu)

#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMoverCreator.hh>

#include <basic/prof.hh>

// Core Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/tag/Tag.hh>

// Numeric Headers
#include <numeric/random/random.hh>


// C++ Headers
#include <ostream>
#include <sstream>

#include <utility/string_util.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <core/pack/task/operation/TaskOperation.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

#ifdef WIN_PYROSETTA
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


using namespace core;
using namespace core::pose;
using namespace basic;
using namespace protocols::moves;


static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.sidechain_moves.SidechainMCMover" );

namespace protocols {
namespace simple_moves {
namespace sidechain_moves {

// XRW TEMP std::string
// XRW TEMP SidechainMCMoverCreator::keyname() const {
// XRW TEMP  return SidechainMCMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP SidechainMCMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new SidechainMCMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP SidechainMCMover::mover_name() {
// XRW TEMP  return "SidechainMC";
// XRW TEMP }

SidechainMCMover::SidechainMCMover():
	protocols::simple_moves::sidechain_moves::SidechainMover(),
	current_(),
	previous_(),
	best_(),
	temperature_(0),
	ntrials_(0),
	best_energy_(0),
	sfxn_( core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
	inherit_scorefxn_temperature_(false),
	ig_(/* 0 */),
	accepts_(0),
	current_ntrial_(0),
	score_pre_apply_(0),
	score_post_apply_(0),
	metropolis_hastings_mover_(/* 0 */)
{}

SidechainMCMover::SidechainMCMover(
	core::pack::dunbrack::RotamerLibrary const & rotamer_library
):
	protocols::simple_moves::sidechain_moves::SidechainMover(rotamer_library),
	current_(),
	previous_(),
	best_(),
	temperature_(0),
	ntrials_(0),
	best_energy_(0),
	sfxn_( core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
	inherit_scorefxn_temperature_(false),
	ig_(/* 0 */),
	accepts_(0),
	current_ntrial_(0),
	score_pre_apply_(0),
	score_post_apply_(0),
	metropolis_hastings_mover_(/* 0 */)
{}

SidechainMCMover::~SidechainMCMover() {}


void
SidechainMCMover::setup( core::scoring::ScoreFunctionCOP sfxn ){
	ig_ = core::pack::interaction_graph::SimpleInteractionGraphOP( new core::pack::interaction_graph::SimpleInteractionGraph() ); //commented out debug
	//(*sfxn)(pose); //gets called in apply
	set_scorefunction( *sfxn );
	ig_->set_scorefunction( *sfxn );
	//ig_->initialize( pose ); //gets called in apply. we can't count on the graph being up-to-date between setup and apply
}

void
SidechainMCMover::show_counters( std::ostream & out){
	out << "SCMCMover: trials " << current_ntrial_ << " accepts= " << (accepts_/current_ntrial_) << std::endl;
}

protocols::moves::MoverOP
SidechainMCMover::clone() const {
	return( protocols::moves::MoverOP( new protocols::simple_moves::sidechain_moves::SidechainMCMover( *this ) ) );
}

protocols::moves::MoverOP
SidechainMCMover::fresh_instance() const {
	return (protocols::moves::MoverOP( new SidechainMCMover ));
}

/// @details
void
SidechainMCMover::apply(
	Pose & pose
)
{

	bool const DEBUG = false;

	using namespace core::scoring;
	using namespace protocols;
	using namespace protocols::moves;
	//using namespace protocols::fast_sc_mover;

	utility::vector1< core::Real > new_chi;
	core::Real current_energy = sfxn_->score(pose);
	score_pre_apply_ = current_energy;
	init_task(pose);

	current_.resize(pose.size());
	previous_.resize(pose.size());
	best_.resize(pose.size());
	runtime_assert(temperature_ != 0);
	runtime_assert(ntrials_ != 0);
	runtime_assert(packed_residues().size() > 0);

	if ( inherit_scorefxn_temperature_ ) {
		protocols::canonical_sampling::MetropolisHastingsMoverOP metropolis_hastings_mover( metropolis_hastings_mover_ );
		runtime_assert(metropolis_hastings_mover != 0);
		// update temperature every time in case temperature changes are implemented in the future
		set_temperature(metropolis_hastings_mover->monte_carlo()->temperature());
	}

	// for debugging
	pose::Pose temp(pose);
	pose::Pose dummy(pose);


	for ( core::Size itr = 1; itr <= pose.size(); itr++ ) {
		current_[ itr ] = core::conformation::ResidueOP( new core::conformation::Residue(pose.residue( itr )) );
	}

	// PROF_START( SIMPLEINTGRAPH );
	//SimpleInteractionGraphOP ig(new SimpleInteractionGraph()); //commented out debug
	//ig->set_scorefunction( sfxn_ ); //commented out debug
	ig_->initialize( pose ); //commented out debug
	SidechainMover::init_task( pose );
	// PROF_STOP( SIMPLEINTGRAPH  );
	//runtime_assert(ig_ != 0);

	for ( core::Size iter_i = 1; iter_i <= ntrials_; iter_i++ ) {
		//pick randomly
		core::Size rand_res;
		do {
			//pick residue, respecting underlying packer task
			rand_res = packed_residues()[ numeric::random::rg().random_range(1, packed_residues().size()) ];
		}while ( pose.residue( rand_res ).name1() == 'P'); //SidechainMover cannot sample proline rotamers? (is this true?)

		core::conformation::ResidueOP new_state( new core::conformation::Residue( pose.residue( rand_res ) ) );

		/// APL Note: you have to remove output if you're trying to optimize.
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "old-chi-angles: ";
			for ( unsigned int i = 1; i <= new_state->nchi(); i++ ) {
				TR.Debug << new_state->chi( i ) << " ";
			}
			TR.Debug << std::endl;
		}
		//  if( !task_initialized() ){
		//   init_task( pose );
		//  }
		new_state = make_move( new_state );
		//new_state->update_actcoord();
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "new-chi-angles: ";
			for ( unsigned int i = 1; i <= new_state->nchi(); i++ ) {
				TR.Debug << new_state->chi( i ) << " ";
			}
			TR.Debug << std::endl;
		}


		PROF_START( SIMPLEINTGRAPH );
		core::Real delta_energy = ig_->consider_substitution( rand_res, new_state, *new_state->nonconst_data_ptr() );
		PROF_STOP( SIMPLEINTGRAPH );

		if ( DEBUG ) {
			//debug
			core::Real s1= sfxn_->score(temp);
			dummy = temp;
			//for(unsigned int i = 1; i <= new_state->nchi(); i++ ){
			//dummy.set_chi( i, rand_res, new_state->chi( i ) );
			//}
			dummy.replace_residue( rand_res, *new_state, true );
			core::Real s2 = sfxn_->score(dummy);

			if ( (s1 - s2) - delta_energy > 0.05 ) {
				TR.Warning << "ENERGIES DON'T MATCH UP! " << s1 << " " << s2 << " " << (s1 - s2) << " " << delta_energy << std::endl;
				dummy.dump_pdb("dummy.pdb");
				temp.dump_pdb("temp.pdb");
				//exit(1);

			} else {
				TR.Debug << "energies match up " << std::endl;
			}
		}//debug


		if ( pass_metropolis( delta_energy, SidechainMover::last_proposal_density_ratio() ) ) { //ek
			if ( DEBUG ) {
				//debug//
				core::Real s1= sfxn_->score(temp);
				temp.replace_residue( rand_res, *new_state, true );
				core::Real s2 = sfxn_->score(dummy);
				TR.Debug << "current energy after accept: " << sfxn_->score(temp) << " delta " << (s2-s1) << std::endl;
				//for( core::Size itr_i = 1; itr_i <= new_state->nchi(); itr_i++ ){
				//temp.set_chi( itr_i, rand_res, new_state->chi(itr_i) );
				//}
			}

			if ( TR.visible( basic::t_debug ) ) {
				TR.Debug << "passed metropolis! assigning new move to current " << std::endl;
			}
			previous_[ rand_res ] = current_[ rand_res ] ;

			PROF_START( SIMPLEINTGRAPH );
			ig_->commit_change( rand_res );
			PROF_STOP( SIMPLEINTGRAPH );
			current_energy -= delta_energy;
			current_[ rand_res ] =   new_state;
			if ( ( current_energy ) <  best_energy_ ) {
				best_[ rand_res ] = new_state;
				best_energy_ = current_energy;
			}
		} else { //rejected metropolis criterion
			PROF_START( SIMPLEINTGRAPH );
			ig_->reject_change( rand_res, new_state, *new_state->nonconst_data_ptr() );
			PROF_STOP( SIMPLEINTGRAPH );
		}
	} // n_iterations


	for ( core::Size res_i = 1; res_i <= current_.size(); res_i++ ) {
		//for(core::Size chi_i = 1; chi_i <= current_[ res_i ]->nchi(); chi_i++){
		//pose.set_chi( chi_i, res_i, current_[ res_i ]->chi( chi_i ) );
		//}
		pose.replace_residue( res_i, (*current_[ res_i ]), true );
	}

	score_post_apply_ = current_energy;
	if ( !metropolis_hastings_mover_.expired() ) {
		score_post_apply_ = sfxn_->score(pose);
		//TR << "Score Actual: " << score_post_apply_ << " Accumulated: " << current_energy << " Delta: " << current_energy - score_post_apply_ << std::endl;
	}

	type("sc_mc");
}

// XRW TEMP std::string
// XRW TEMP SidechainMCMover::get_name() const {
// XRW TEMP  return "SidechainMCMover";
// XRW TEMP }

bool
SidechainMCMover::pass_metropolis(core::Real delta_energy , core::Real last_proposal_density_ratio ){

	core::Real boltz_factor = delta_energy / temperature_;
	core::Real probability = std::exp( std::min( 40.0, std::max( -40.0, boltz_factor ))) *  last_proposal_density_ratio ;

	if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {
		current_ntrial_++;
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "delta energy is " << delta_energy << " probability: " << probability << " accepted: FALSE " << std::endl;
		}
		return false;
	} else {
		accepts_++;
		current_ntrial_++;
		if ( TR.visible( basic::t_debug ) ) {
			TR.Debug << "delta energy is " << delta_energy << " probability: " << probability << " accepted: TRUE" << std::endl;
		}
		return true;
	}
}

void
SidechainMCMover::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose) {

	// code duplication: should really call SidechainMover::parse_my_tag() instead of having most of the code below
	if ( tag->hasOption("task_operations") ) {

		core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );

		std::string const t_o_val( tag->getOption<std::string>("task_operations") );
		typedef utility::vector1< std::string > StringVec;
		StringVec const t_o_keys( utility::string_split( t_o_val, ',' ) );
		for ( StringVec::const_iterator t_o_key( t_o_keys.begin() ), end( t_o_keys.end() );
				t_o_key != end; ++t_o_key ) {
			if ( data.has( "task_operations", *t_o_key ) ) {
				new_task_factory->push_back( data.get_ptr< core::pack::task::operation::TaskOperation >( "task_operations", *t_o_key ) );
			} else {
				throw utility::excn::EXCN_RosettaScriptsOption("TaskOperation " + *t_o_key + " not found in basic::datacache::DataMap.");
			}
		}

		set_task_factory(new_task_factory);

	} else {

		pack::task::PackerTaskOP pt = core::pack::task::TaskFactory::create_packer_task( pose );
		set_task( pt );
		pt->restrict_to_repacking();
		// probably should call init_task(), but it's cleaner if we wait untill it's called in apply()
	}

	ntrials_ = tag->getOption<core::Size>( "ntrials", 10000 );
	set_preserve_detailed_balance( tag->getOption<bool>( "preserve_detailed_balance", 1 ) );
	temperature_ = tag->getOption<core::Real>( "temperature", 1.0 );
	set_inherit_scorefxn_temperature( tag->getOption<bool>( "inherit_scorefxn_temperature", inherit_scorefxn_temperature() ) );

	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", 0.0 ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", 0.0 ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", 0.0 ) );
	core::Real between_rot = 1.0 - prob_uniform() - prob_withinrot () - prob_random_pert_current();

	setup( rosetta_scripts::parse_score_function( tag, data )->clone() );

	TR
		<< "Initialized SidechainMCMover from .xml file:"
		<< " ntrials=" << ntrials_
		<< " temperature= " << temperature_
		<< " detailed balance= " << (preserve_detailed_balance()?"true":"false")
		<< " scorefunction=" << rosetta_scripts::get_score_function_name(tag)
		<< " Probablities uniform/withinrot/random_pert/betweenrot: "
		<< prob_uniform() << '/' <<prob_withinrot() << '/' << prob_random_pert_current() << '/' << between_rot
		<< std::endl;
}

void
SidechainMCMover::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
)
{
	SidechainMover::initialize_simulation(pose, metropolis_hastings_mover,cycle);

	if ( inherit_scorefxn_temperature_ ) {
		protocols::canonical_sampling::MetropolisHastingsMoverOP metropolis_hastings_mover( metropolis_hastings_mover_ );
		runtime_assert(metropolis_hastings_mover != 0);
		set_scorefunction(metropolis_hastings_mover->monte_carlo()->score_function());
		set_temperature(metropolis_hastings_mover->monte_carlo()->temperature());
	}
}

std::string SidechainMCMover::get_name() const {
	return mover_name();
}

std::string SidechainMCMover::mover_name() {
	return "SidechainMC";
}

void SidechainMCMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	attlist + XMLSchemaAttribute::attribute_w_default("ntrials", xsct_non_negative_integer, "Number of trials.", "10000" )
		+ XMLSchemaAttribute::attribute_w_default("preserve_detailed_balance", xsct_rosetta_bool, "Should the simulation preserve detailed balance?", "1" )
		+ XMLSchemaAttribute::attribute_w_default("temperature", xsct_real, "Simulation temperature.", "1.0" )
		+ XMLSchemaAttribute( "inherit_scorefxn_temperature", xsct_rosetta_bool, "Should the simulation inherit the temperature provided by the scorefunction?" )
		+ XMLSchemaAttribute::attribute_w_default("prob_uniform", xsct_real, "The probability of uniform chi sampling.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("prob_withinrot", xsct_real, "The probability of sampling within-rotamer.", "0.0" )
		+ XMLSchemaAttribute::attribute_w_default("prob_random_pert_current", xsct_real, "The probability of making a random perturbation to the current chi value.", "0.0" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "In a Monte Carlo simulation, moves the side chain for a set of residues identified by a task operation in a manner that can be totally independent of rotamer assignments", attlist );
}

std::string SidechainMCMoverCreator::keyname() const {
	return SidechainMCMover::mover_name();
}

protocols::moves::MoverOP
SidechainMCMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new SidechainMCMover );
}

void SidechainMCMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SidechainMCMover::provide_xml_schema( xsd );
}



} // sidechain_moves
} // simple_moves
} // protocols
