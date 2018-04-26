// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

// Headers
#include <protocols/loop_modeling/LoopProtocol.hh>
#include <protocols/loop_modeling/LoopProtocolCreator.hh>

// Protocols headers
#include <protocols/filters/Filter.hh>
#include <protocols/loop_modeling/LoopMover.hh>
#include <protocols/loop_modeling/utilities/rosetta_scripts.hh>
#include <protocols/loop_modeling/utilities/LoopMoverGroup.hh>
#include <protocols/loop_modeling/utilities/AcceptanceCheck.hh>
#include <protocols/loop_modeling/utilities/TrajectoryLogger.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/mover_schemas.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/util.hh>

// Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <boost/lexical_cast.hpp>
#include <boost/xpressive/xpressive.hpp>
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tools/make_vector1.hh>

// C++ headers
#include <iostream>
#include <cmath>

// Namespaces
using namespace std;

using core::Real;
using core::Size;
using core::pose::Pose;
using core::scoring::ScoreFunctionOP;
using core::scoring::ScoreFunctionCOP;
using protocols::loops::Loop;
using protocols::loops::Loops;
using protocols::moves::MoverOP;
using protocols::moves::MonteCarloOP;
using protocols::loop_modeling::utilities::LoopMoverGroup;
using protocols::loop_modeling::utilities::LoopMoverGroupOP;
using protocols::loop_modeling::utilities::TrajectoryLogger;
using protocols::loop_modeling::utilities::TrajectoryLoggerOP;
using utility::tools::make_vector1;

using utility::tag::TagCOP;
using basic::datacache::DataMap;
using protocols::filters::Filters_map;
using protocols::moves::Movers_map;

namespace protocols {
namespace loop_modeling {

static basic::Tracer TR("protocols.loop_modeling.LoopProtocol");

LoopProtocol::LoopProtocol() {
	protocol_ = add_child( LoopMoverGroupOP( new LoopMoverGroup ) );
	movers_ = protocol_->add_mover_group();
	refiners_ = protocol_->add_mover_group();
	logger_ = TrajectoryLoggerOP( new TrajectoryLogger );
	monte_carlo_ = nullptr;

	sfxn_cycles_ = 1;
	temp_cycles_ = 1;
	mover_cycles_ = 1;
	ramp_sfxn_rep_ = false;
	ramp_sfxn_rama_ = false;
	ramp_sfxn_cst_ = false;
	ramp_temp_ = true;
	scale_temp_cycles_ = false;
	initial_temp_ = 1.5;
	final_temp_ = 0.5;
	original_rama_weight_ = 0.0;       // Initialized in start_protocol().
	original_rama2b_weight_ = 0.0;     // Ditto.
	original_repulsive_weight_ = 0.0;  // Ditto.
	test_run_ = false;
}

LoopProtocol::~LoopProtocol() = default;

void LoopProtocol::parse_my_tag(
	TagCOP tag,
	DataMap & data,
	Filters_map const & filters,
	Movers_map const & movers,
	Pose const & pose) {


	LoopMover::parse_my_tag(tag, data, filters, movers, pose);
	utilities::set_scorefxn_from_tag(*this, tag, data);

	sfxn_cycles_ = tag->getOption<Size>("sfxn_cycles", sfxn_cycles_);
	mover_cycles_ = tag->getOption<Size>("mover_cycles", mover_cycles_);
	ramp_sfxn_rama_ = tag->getOption<bool>("ramp_rama", ramp_sfxn_rama_);
	ramp_sfxn_rep_ = tag->getOption<bool>("ramp_rep", ramp_sfxn_rep_);
	ramp_sfxn_cst_ = tag->getOption<bool>("ramp_down_cst", ramp_sfxn_cst_);
	ramp_temp_ = tag->getOption<bool>("ramp_temp", ramp_temp_);
	initial_temp_ = tag->getOption<Real>("initial_temp", initial_temp_);
	final_temp_ = tag->getOption<Real>("final_temp", final_temp_);

	// Parse the 'temp_cycles' tag.

	// The temp_cycles tag is complicated because it's a number optionally
	// followed by an 'x'.  If the 'x' is present, it means the number of cycles
	// should be the specified number multiplied by the loop length.  Otherwise
	// the number of cycles will just be the specified number.

	if ( tag->hasOption("temp_cycles") ) {
		using namespace boost::xpressive;

		string cycles = tag->getOption<string>("temp_cycles");
		sregex regex = (s1 = +digit) >> (s2 = !as_xpr('x'));
		smatch match;

		if ( regex_match(cycles, match, regex) ) {
			temp_cycles_ = boost::lexical_cast<Size>(std::string(match[1]));
			scale_temp_cycles_ = (match[2] == "x");
		} else {
			stringstream message;
			message << "getOption: key = temp_cycles stream extraction failed! ";
			message << "Tried to parse '" << cycles << "'\n";
			throw CREATE_EXCEPTION(utility::excn::Exception, message.str());
		}
	}

	// Parse the 'auto_refine' tag.

	// By default, the standard set of refiners will be applied after any
	// manually specified movers.  This makes it easy to experiment with
	// different movers without messing up the core protocol.  However, if you
	// want to specify your own refinement algorithm, you need to disable this
	// behavior via the auto_refine flag.

	if ( ! tag->getOption<bool>("auto_refine", true) ) {
		clear_refiners();
	}

	// Parse the 'fast' tag.

	if ( tag->getOption<bool>("fast", false) ) {
		mark_as_test_run();
	}

	// Add loop movers to this protocol by parsing the subtags.

	for ( TagCOP subtag : tag->getTags() ) {

		// Ignore <Loop> subtags (parsed by parent class).

		if ( subtag->getName() == "Loop" ) { continue; }

		// Parse <AcceptanceCheck> subtags.

		else if ( subtag->getName() == "AcceptanceCheck" ) {
			string name = subtag->getOption<string>("name", "loop_modeling");
			add_acceptance_check(name);
		} else {
			// Parse LoopMover subtags.
			LoopMoverOP loop_mover = utilities::loop_mover_from_tag(subtag, data, filters, movers, pose);
			add_mover(loop_mover);
		}
	}
}

void LoopProtocol::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen = complex_type_generator_for_loop_protocol( xsd );
	ct_gen->write_complex_type_to_schema( xsd );
}

utility::tag::AttributeList
LoopProtocol::loop_protocol_attlist( ) {
	using namespace utility::tag;

	AttributeList attlist;

	// Create a complex type and  get the LoopMover attributes, as parse_my_tag calls LoopMover::parse_my_tag

	attlist
		+ XMLSchemaAttribute("sfxn_cycles", xsct_non_negative_integer, "Number of iterations to make in the sfxn loop.")
		+ XMLSchemaAttribute("mover_cycles", xsct_non_negative_integer, "The number of iterations to make in the mover loop")
		+ XMLSchemaAttribute("ramp_rama", xsct_rosetta_bool,
		"If enabled, the Ramachandran weight will start near zero and"
		" will finish at whatever it was in the original score function.")
		+ XMLSchemaAttribute("ramp_rep", xsct_rosetta_bool,
		"If enabled, the repulsive weight will start near zero and"
		"will finish at whatever it was in the original score function")
		+ XMLSchemaAttribute("ramp_down_cst", xsct_rosetta_bool, "If enabled, the constraint weight will start at whatever value in the scorefunction and end at 0.")
		//+ XMLSchemaAttribute("ramp_down_cst_fast", xsct_rosetta_bool, "If enabled, the constraint weight will fluctuate between its starting value in the scorefunction and close to 0 along with temperature.")
		+ XMLSchemaAttribute("ramp_temp", xsct_rosetta_bool, "Ramp the temperature during the temp loop.")
		+ XMLSchemaAttribute("initial_temp", xsct_real, "Initial temperature. Ignored if temperature ramping is disabled.")
		+ XMLSchemaAttribute("final_temp", xsct_real, "Final temperature. Ignored if temperature ramping is disabled.");
	return attlist;
}

utility::tag::XMLSchemaComplexTypeGeneratorOP
LoopProtocol::complex_type_generator_for_loop_protocol( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;


	AttributeList attlist = loop_protocol_attlist();

	// Restriction for temp_cycles
	XMLSchemaRestriction restriction_type;
	restriction_type.name( "begin_w_int" );
	restriction_type.base_type( xs_string );
	restriction_type.add_restriction( xsr_pattern , "^[0-9]" );
	xsd.add_top_level_element( restriction_type );
	attlist
		+ XMLSchemaAttribute("temp_cycles", "begin_w_int",
		"The number of iterations to make in the temp loop. This number may optionally be followed by an \"x\","
		"in which case the number of iterations will be the given number times the length of the loop being sampled.")
		+ XMLSchemaAttribute::attribute_w_default("auto_refine", xsct_rosetta_bool,
		"Invoke the built-in refiners after any user-specified movers.","true")
		+ XMLSchemaAttribute::attribute_w_default("fast", xsct_rosetta_bool, "Test mode, reduces number of cycles.","false");

	//subelements
	XMLSchemaSimpleSubelementList subelement_list;

	AttributeList subelement_attributes;
	subelement_attributes
		+ XMLSchemaAttribute::attribute_w_default("name", xs_string, "\"name\" attribute of mover after which to do AcceptanceCheck.", "loop_modeling");
	subelement_list.add_simple_subelement("AcceptanceCheck", subelement_attributes,
		"Add a Monte Carlo score function evaluation"
		"and acceptance check between any of your movers.");

	XMLSchemaComplexTypeGeneratorOP ct_gen( new XMLSchemaComplexTypeGenerator );
	LoopMover::define_composition_schema( xsd, *ct_gen, subelement_list );
	ct_gen->element_name( mover_name() )
		.description( "Optimizes a region of protein backbone using a simulated annealing MonteCarlo simulation." )
		.add_attributes( attlist  )
		.set_subelements_repeatable( subelement_list );

	return ct_gen;
}

bool LoopProtocol::do_apply(Pose & pose) {
	start_protocol(pose);

	for ( Size i = 1; i <= get_sfxn_cycles(); i++ ) {
		monte_carlo_->recover_low(pose);
		logger_->record_new_pose(TR, pose);

		ramp_score_function(i);
		logger_->record_new_score_function(TR, pose);

		for ( Size j = 1; j <= get_temp_cycles(); j++ ) {
			ramp_temperature(j);

			for ( Size k = 1; k <= get_mover_cycles(); k++ ) {
				attempt_loop_move(pose, i, j, k);
			}
		}
		// The legacy code repacks and minimizes here.  This doesn't really
		// translate to the new framework, because the "kic" loop mover isn't
		// distinct from the "repack" and "minimize" loop movers.
	}

	pose = monte_carlo_->lowest_score_pose();
	finish_protocol(pose);

	return true;
}

void LoopProtocol::start_protocol(Pose & pose) {
	using core::scoring::fa_rep;
	using core::scoring::rama;
	using core::scoring::rama2b;
	using core::scoring::constraints::add_fa_constraints_from_cmdline;
	using core::scoring::constraints::add_constraints_from_cmdline;
	using core::scoring::EnergyMap;

	if ( movers_->empty() && refiners_->empty() ) {
		utility_exit_with_message("Refusing to run LoopProtocol without any movers or refiners.");
	}

	// Prepare the score function.

	ScoreFunctionOP standard_sfxn = core::scoring::get_score_function();
	ScoreFunctionOP scorefxn = get_score_function();

	protocols::loops::add_cutpoint_variants(pose);
	protocols::loops::loop_mover::loops_set_chainbreak_weight(scorefxn, 1);

	if ( pose.is_fullatom() ) {
		add_fa_constraints_from_cmdline(pose, *scorefxn);
	} else {
		add_constraints_from_cmdline(pose, *scorefxn);
	}

	monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo(
		pose, *scorefxn, initial_temp_) );

	// Use the rama2b scoring term instead of the rama term

	original_repulsive_weight_ = scorefxn->get_weight(fa_rep);
	original_rama_weight_ = 0;
	original_rama2b_weight_ = standard_sfxn->get_weight(rama);

	original_sfxn_weights_ = scorefxn->weights();

	// Describe the protocol to the user.

	vector1<string> algorithm_names;
	movers_->get_children_names(algorithm_names);
	refiners_->get_children_names(algorithm_names);

	TR << "sfxn_cycles:    " << get_sfxn_cycles() << endl;
	TR << "temp_cycles:    " << get_temp_cycles() << endl;
	TR << "mover_cycles:   " << get_mover_cycles() << endl;
	TR << "ramp_sfxn_rep:  " << ramp_sfxn_rep_ << endl;
	TR << "ramp_sfxn_rama: " << ramp_sfxn_rama_ << endl;
	TR << "ramp_sfxn_cst:  " << ramp_sfxn_cst_ << endl;
	TR << "ramp_temp:      " << ramp_temp_ << endl;
	TR << "initial_temp:   " << initial_temp_ << endl;
	TR << "final_temp:     " << final_temp_ << endl;
	for ( string const & name : algorithm_names ) {
		TR << "loop_movers:    " << name << endl;
	}

	// Initialize the trajectory reporter.

	logger_->init(this);
}

void LoopProtocol::ramp_score_function(Size iteration) {
	using namespace core::scoring;

	ScoreFunctionOP score_function = get_score_function();
	Real ramp_factor = 1.0 / (get_sfxn_cycles() - iteration + 1);
	Real ramp_factor_down = 1.0 - (1.0 / (get_sfxn_cycles() - iteration + 1) ) + (1.0 / (1.0 + get_sfxn_cycles() ) );

	// Ramp the chainbreak score term.

	score_function->set_weight(chainbreak, 3.33 * iteration);

	// Ramp the repulsive score term.

	if ( ramp_sfxn_rep_ ) {
		score_function->set_weight(
			fa_rep, original_repulsive_weight_ * ramp_factor);
	}

	// Ramp constraint score terms.

	if ( ramp_sfxn_cst_ ) {
		utility::vector1 < core::scoring::ScoreType > constraint_list {
			atom_pair_constraint,
			base_pair_constraint,
			coarse_chainbreak_constraint,
			constant_constraint,
			coordinate_constraint,
			angle_constraint,
			dihedral_constraint,
			big_bin_constraint,
			dunbrack_constraint,
			site_constraint,
			metalhash_constraint,
			metalbinding_constraint
			};
		for ( utility::vector1< core::scoring::ScoreType >::iterator it = constraint_list.begin(); it != constraint_list.end(); ++ it ) {
			score_function->set_weight(
				*it, original_sfxn_weights_[*it] * ramp_factor_down);
		}
	}

	// Ramp the rama and/or rama2b score terms.

	if ( ramp_sfxn_rama_ ) {
		score_function->set_weight(
			rama, original_rama_weight_ * ramp_factor);
		score_function->set_weight(
			rama2b, original_rama2b_weight_ * ramp_factor);
	}

	// The MonteCarlo object has to be explicitly updated with the new score
	// function.  You might think that since everything is sharing the same
	// pointer that this step would be unnecessary, but it seems that the
	// MonteCarlo object makes a copy or something.

	monte_carlo_->score_function(*score_function);
}

void LoopProtocol::ramp_temperature(Size /*iteration*/) {
	if ( ramp_temp_ ) {
		Real temperature = monte_carlo_->temperature();
		Real ramp_factor = std::pow(
			final_temp_ / initial_temp_,
			1.0 / (get_sfxn_cycles() * get_temp_cycles()));

		monte_carlo_->set_temperature(temperature * ramp_factor);

		logger_->record_new_temperature(TR, monte_carlo_->temperature());
	}
}

void LoopProtocol::attempt_loop_move(Pose & pose, Size i, Size j, Size k) {
	protocol_->apply(pose);

	bool move_proposed = protocol_->was_successful();
	bool move_accepted = false;

	if ( move_proposed ) {
		move_accepted = monte_carlo_->boltzmann(pose);
	} else {
		pose = monte_carlo_->last_accepted_pose();
	}

	logger_->record_move(TR, pose, i, j, k, move_proposed, move_accepted);
}

void LoopProtocol::finish_protocol(Pose & pose) {
	logger_->record_endpoint(TR, pose);
}

void LoopProtocol::add_mover(LoopMoverOP mover) {
	movers_->add_mover(mover);
}

void LoopProtocol::add_refiner(LoopMoverOP refiner) {
	refiners_->add_mover(refiner);
}

void LoopProtocol::add_filter(protocols::filters::FilterOP filter) {
	movers_->add_filter(filter);
}

void LoopProtocol::add_acceptance_check(string name) {
	movers_->add_mover(LoopMoverOP( new utilities::AcceptanceCheck(monte_carlo_, name) ));
}

void LoopProtocol::clear_movers() {
	movers_->clear();
}

void LoopProtocol::clear_refiners() {
	refiners_->clear();
}

void LoopProtocol::clear_movers_and_refiners() {
	clear_movers();
	clear_refiners();
}

void LoopProtocol::mark_as_default() {
	movers_->mark_as_default();
	refiners_->mark_as_default();
}

ScoreFunctionOP LoopProtocol::get_score_function() const {
	return get_tool<ScoreFunctionOP>(ToolboxKeys::SCOREFXN);
}

void LoopProtocol::set_score_function(ScoreFunctionOP score_function) {
	set_tool(ToolboxKeys::SCOREFXN, score_function);
}

Size LoopProtocol::get_sfxn_cycles() const {
	return test_run_ ? 3 : sfxn_cycles_;
}

Size LoopProtocol::get_temp_cycles() const {
	if ( test_run_ ) { return 3; }
	return scale_temp_cycles_ ?
		temp_cycles_ * get_loops()->loop_size() :
		temp_cycles_;
}

Size LoopProtocol::get_mover_cycles() const {
	return test_run_ ? 2 : mover_cycles_;
}

void LoopProtocol::set_sfxn_cycles(Size x) {
	sfxn_cycles_ = x;
}

void LoopProtocol::set_temp_cycles(Size x, bool times_loop_length) {
	temp_cycles_ = x;
	scale_temp_cycles_ = times_loop_length;
}

void LoopProtocol::set_mover_cycles(Size x) {
	mover_cycles_ = x;
}

void LoopProtocol::mark_as_test_run() {
	test_run_ = true;
}

void LoopProtocol::set_repulsive_term_ramping(bool value) {
	ramp_sfxn_rep_ = value;
}

void LoopProtocol::set_rama_term_ramping(bool value) {
	ramp_sfxn_rama_ = value;
}

void LoopProtocol::set_cst_term_ramping(bool value) {
	ramp_sfxn_cst_ = value;
}

void LoopProtocol::set_temperature_ramping(bool value) {
	ramp_temp_ = value;
}

void LoopProtocol::set_temperature_schedule(Real start, Real stop) {
	initial_temp_ = start;
	final_temp_ = stop;
}

TrajectoryLoggerOP LoopProtocol::get_logger() const {
	return logger_;
}

std::string LoopProtocol::get_name() const {
	return mover_name();
}

std::string LoopProtocol::mover_name() {
	return "LoopProtocol";
}

std::string LoopProtocolCreator::keyname() const {
	return LoopProtocol::mover_name();
}

protocols::moves::MoverOP
LoopProtocolCreator::create_mover() const {
	return protocols::moves::MoverOP( new LoopProtocol );
}

void LoopProtocolCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LoopProtocol::provide_xml_schema( xsd );
}


}
}
