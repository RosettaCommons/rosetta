// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/monte_carlo/MonteCarloInterface.cc
/// @brief A MonteCarlo object for optimizing the interface dG as defined using InterfaceAnalyzer. The dG and Total energy can be weighted.  This is so that the interface energy itself can be optimized through a protocol.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/monte_carlo/MonteCarloInterface.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>


#include <basic/options/option.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh> // AUTO IWYU For ScoreFunction

static basic::Tracer TR( "protocols.monte_carlo.MonteCarloInterface" );


namespace protocols {
namespace monte_carlo {

using namespace protocols::moves;
using namespace core::scoring;


MonteCarloInterface::MonteCarloInterface(
	Pose const & init_pose, // PoseCOP init_pose,
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature,
	core::pose::DockingPartners const & interface,
	bool detect_disulfide_in_separated_pose
):

	MonteCarlo( init_pose, scorefxn, temperature),
	detect_disulfide_in_separated_pose_(detect_disulfide_in_separated_pose)

{
	init_interface_analyzer();
	set_interface( interface );
	read_cmd_line_options();
	MonteCarloInterface::reset(init_pose); // Note that the baseclass reset() has already been called in its constructor
}


MonteCarloInterface::MonteCarloInterface(
	ScoreFunction const & scorefxn, // ScoreFunctionCOP scorefxn,
	Real const temperature,
	core::pose::DockingPartners const & interface):

	protocols::moves::MonteCarlo( scorefxn, temperature )

{

	init_interface_analyzer();
	set_interface( interface );
	read_cmd_line_options();

}

MonteCarloInterface::~MonteCarloInterface(){}

MonteCarloInterface::MonteCarloInterface( MonteCarloInterface const & src):
	MonteCarlo( src ),
	interface_definition_( src.interface_definition_ ),
	dG_weight_( src.dG_weight_ ),
	total_weight_( src.total_weight_ ),
	last_dG_( src.last_dG_ ),
	last_score_( src.last_score_ )

{
	init_interface_analyzer();
}

void
MonteCarloInterface::init_interface_analyzer() {
	analyzer_ = utility::pointer::make_shared< analysis::InterfaceAnalyzerMover >();

	//Turn most computations off here.
	analyzer_->set_compute_interface_delta_hbond_unsat( false );
	analyzer_->set_compute_packstat( false );
	analyzer_->set_compute_interface_sc( false );
	analyzer_->set_calc_hbond_sasaE( false );
	analyzer_->set_use_centroid_dG( false );
	analyzer_->set_detect_disulfide_in_separated_pose( detect_disulfide_in_separated_pose_ );
	TR << "Detection of disulfides in the separated pose is set to " << detect_disulfide_in_separated_pose_ << std::endl;

	//Make sure pack separated is set to TRUE!!
	analyzer_->set_pack_separated( repack_separated_ );
	ScoreFunctionOP scoreop =  moves::MonteCarlo::score_function().clone();
	analyzer_->set_scorefunction( scoreop );

	//Add SKIP dSASA
	analyzer_->set_calc_dSASA( false );

	//analyzer_->set_pack_rounds( 2 );

}

void
MonteCarloInterface::read_cmd_line_options() {
	using namespace basic::options;

	dG_weight_ = option[ OptionKeys::score::mc_interface_weight].value();
	total_weight_ = option[ OptionKeys::score::mc_total_weight].value();

}


void
MonteCarloInterface::set_interface( core::pose::DockingPartners const & interface ){
	interface_definition_ = interface;
	analyzer_->set_interface( interface );
}

void
MonteCarloInterface::set_dG_weight(core::Real dG_weight){
	dG_weight_ = dG_weight;
}
void
MonteCarloInterface::set_total_weight(core::Real total_weight){
	total_weight_ = total_weight;
}

void
MonteCarloInterface::set_pack_interface_separated(bool pack_separated){
	analyzer_->set_pack_separated( pack_separated );
}

void
MonteCarloInterface::score_function(const ScoreFunction &scorefxn){
	ScoreFunctionOP scoreop = scorefxn.clone();
	analyzer_->set_scorefunction( scoreop );
	MonteCarlo::score_function( scorefxn );

	core::Real score = calculate_score( last_accepted_pose_non_const() );
	set_last_accepted( score );
	set_lowest( score );


}

core::Real
MonteCarloInterface::calculate_score(Pose &pose){

	analyzer_->apply( pose );
	core::Real interface_score = analyzer_->get_interface_dG();
	last_dG_ = interface_score;
	TR.Debug << "Interface Score " << interface_score << std::endl;

	core::scoring::ScoreFunction const & scorefxn = moves::MonteCarlo::score_function();
	core::Real total_score = scorefxn.score( pose );
	TR.Debug << "Total Score " << total_score << std::endl;

	core::Real final_score = (dG_weight_ * interface_score) + (total_weight_ * total_score );
	last_score_ = final_score;
	TR.Debug << "Final Score " << final_score << std::endl;

	return final_score;
}


bool
MonteCarloInterface::boltzmann(
	Pose & pose,
	std::string const & move_type,
	core::Real const proposal_density_ratio,
	core::Real const inner_score_delta_over_temperature
){

	//Calculate the score.

	core::Real score = calculate_score( pose );
	return MonteCarlo::boltzmann(score, pose, move_type, proposal_density_ratio, inner_score_delta_over_temperature);
}

void
MonteCarloInterface::reset( Pose const & pose )
{



	set_last_accepted_pose( pose );

	Real const score = calculate_score( last_accepted_pose_non_const() );

	set_last_accepted( score );
	set_lowest_score_pose( last_accepted_pose() , score);
	set_last_score( score );

}

void
MonteCarloInterface::set_repack_separated( bool repack_separated ){
	repack_separated_ = repack_separated;
	analyzer_->set_pack_separated( repack_separated_ );
}

MonteCarloOP
MonteCarloInterface::clone() {
	return utility::pointer::make_shared< MonteCarloInterface >( *this );
}


} //protocols
} //monte_carlo






