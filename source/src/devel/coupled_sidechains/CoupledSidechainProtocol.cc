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
/// @author Oliver Lange

#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <devel/coupled_sidechains/CoupledSidechainProtocol.hh>
//#include <devel/coupled_sidechains/CoupledSidechainProtocolCreator.hh>

#include <basic/prof.hh>

// Core Headers
#include <core/types.hh>

#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/SingleResidueDunbrackLibrary.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>

#include <core/pack/interaction_graph/SimpleInteractionGraph.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/pose/Pose.hh>

#include <basic/Tracer.hh>
#include <basic/basic.hh>
//Auto Headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/canonical_sampling/SilentTrajectoryRecorder.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>

#include <utility/tag/Tag.hh>

// Numeric Headers
#include <numeric/angle.functions.hh>
#include <numeric/constants.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

#include <ObjexxFCL/format.hh>
// Auto-header: duplicate removed #include <basic/Tracer.hh>

// C++ Headers
#include <ostream>
#include <sstream>
#include <fstream>
#include <utility/fixedsizearray1.hh>

#ifdef WIN_PYROSETTA
	#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#endif


using namespace core;
using namespace core::pose;
using namespace basic;
using namespace protocols::moves;
using namespace protocols;

static THREAD_LOCAL basic::Tracer tr( "devel.coupled_sidechains.CoupledSidechainProtocol" );

OPT_1GRP_KEY(Integer,rotamers,traj_interval)
OPT_1GRP_KEY(Integer,rotamers,score_interval)
OPT_1GRP_KEY(Integer,rotamers,rot_interval)
OPT_1GRP_KEY(Real,rotamers,unif)
OPT_1GRP_KEY(Real,rotamers,pert)
OPT_1GRP_KEY(Real,rotamers,within)
OPT_1GRP_KEY(Real,rotamers,pert_magnitude)


bool devel::coupled_sidechains::CoupledSidechainProtocol::options_registered_( false );

void devel::coupled_sidechains::CoupledSidechainProtocol::register_options() {
  using namespace basic::options;
  using namespace OptionKeys;
  if ( options_registered_ ) return;
  options_registered_ = true;
	protocols::canonical_sampling::SimulatedTemperingObserver::register_options();
	//protocols::canonical_sampling::SilentTrajectoryRecorder::register_options();
	OPT( score::weights );
	OPT( score::patch );
	OPT( run::n_cycles );
	NEW_OPT( rotamers::score_interval," how often should the entire score-line be written to file", 100);
	NEW_OPT( rotamers::traj_interval," how often should the entire structure be written to file relative to score_interval", 10);
	NEW_OPT( rotamers::rot_interval," how often should rotamers be written to file... ", 100 );
	NEW_OPT( rotamers::unif,"",0.1);
	NEW_OPT( rotamers::pert,"",0);
	NEW_OPT( rotamers::pert_magnitude, "deviation for perturb_chi moves in degree", 10 );
	NEW_OPT( rotamers::within,"",0);
}

namespace devel {
namespace coupled_sidechains {
/*
std::string
CoupledSidechainProtocolCreator::keyname() const {
	return CoupledSidechainProtocolCreator::mover_name();
}

protocols::moves::MoverOP
CoupledSidechainProtocolCreator::create_mover() const {
	return new CoupledSidechainProtocol;
}

std::string
CoupledSidechainProtocolCreator::mover_name() {
	return "SidechainMC";
}
*/
CoupledSidechainProtocol::CoupledSidechainProtocol():
	protocols::simple_moves::sidechain_moves::SidechainMover(),
	current_(),
	previous_(),
	best_(),
	temperature_(0),
	ntrials_(0),
	best_energy_(0),
	sfxn_(),
	ig_(0),
	accepts_(0),
	current_ntrial_(0),
	score_pre_apply_(0),
	score_post_apply_(0)
{
	//	set_defaults();
	init_from_options();
	setup_objects();
}

CoupledSidechainProtocol::CoupledSidechainProtocol(
	core::pack::dunbrack::RotamerLibrary const & rotamer_library
):
	protocols::simple_moves::sidechain_moves::SidechainMover(rotamer_library),
	current_(),
	previous_(),
	best_(),
	temperature_(0),
	ntrials_(0),
	best_energy_(0),
	sfxn_(),
	ig_(0),
	accepts_(0),
	current_ntrial_(0),
	score_pre_apply_(0),
	score_post_apply_(0)
{
	//	set_defaults();
	init_from_options();
	setup_objects();
}

CoupledSidechainProtocol::~CoupledSidechainProtocol() = default;


void
CoupledSidechainProtocol::set_scorefunction( core::scoring::ScoreFunctionOP sfxn ){
	ig_ = new core::pack::interaction_graph::SimpleInteractionGraph(); //commented out debug
	ig_->set_scorefunction( *sfxn );
	sfxn_ = sfxn;
}

void
CoupledSidechainProtocol::show_counters( std::ostream & out){
	out << "SCMCMover: trials " << current_ntrial_ << " accepts= " << (accepts_/current_ntrial_) << std::endl;
}

protocols::moves::MoverOP
CoupledSidechainProtocol::clone() const {
  return( protocols::moves::MoverOP( new CoupledSidechainProtocol( *this ) ) );
}

protocols::moves::MoverOP
CoupledSidechainProtocol::fresh_instance() const {
	return (protocols::moves::MoverOP( new CoupledSidechainProtocol ));
}

void
CoupledSidechainProtocol::update_rotamers( Size resid ) {
	chi_vectors_[ resid ]=current_[ resid ]->chi();
	core::pack::dunbrack::rotamer_from_chi( *current_[ resid ], rot_vectors_[ resid ] );
}

void
CoupledSidechainProtocol::init_from_options() {
	using namespace options;
	using namespace options::OptionKeys;
	core::Real unif = option[ rotamers::unif ];
	core::Real pert = option[ rotamers::pert ];
	core::Real within = option[ rotamers::within ];
	ntrials_=option[ run::n_cycles ]();
	set_prob_uniform( unif );
	set_prob_withinrot(within );
	set_prob_random_pert_current( pert );
	set_pert_magnitude( option[ rotamers::pert_magnitude ] );

	score_stride_ = option[ rotamers::score_interval ]();
	traj_stride_ = option[ rotamers::traj_interval ]();
	rotamer_stride_ = option[ rotamers::rot_interval ]();

	core::pack::task::TaskFactoryOP new_task_factory( new core::pack::task::TaskFactory );
	new_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
	set_task_factory( new_task_factory );

	set_scorefunction( core::scoring::get_score_function() );
}

void
CoupledSidechainProtocol::setup_objects() {
	tempering_ = new protocols::canonical_sampling::SimulatedTemperingObserver();
	moves::MonteCarloOP mc_object = new moves::MonteCarlo( *sfxn_, 0.6 );
	tempering_->set_monte_carlo( mc_object );
	counters_.set_temperature_observer( tempering_ );
}

/// @details
void
CoupledSidechainProtocol::apply(
	Pose & pose
)
{
	using namespace core::scoring;
	using namespace protocols;
	using namespace protocols::moves;
	//using namespace protocols::fast_sc_mover;

	utility::vector1< core::Real > new_chi;
	core::Real current_energy = (*sfxn_)(pose);

	SidechainMover::init_task( pose );

	std::string const traj_file_tag( jd2::current_output_name() );

	counters_.reset();

	jd2::JobOP job;
	if ( jd2::jd2_used() ) {
		job = jd2::get_current_job();
	}


	current_.resize(pose.size());
	previous_.resize(pose.size());
	best_.resize(pose.size());
	rot_vectors_.resize( pose.size() );
	chi_vectors_.resize( pose.size() );

	runtime_assert(ntrials_ != 0);
	runtime_assert(packed_residues().size() > 0);
	runtime_assert( ig_ );

	for ( core::Size i = 1; i <= pose.size(); i++ ){
		current_[ i ] = new core::conformation::Residue(pose.residue( i ));
		update_rotamers( i );
	}

	if ( tempering_ ) {
		tempering_->initialize_simulation();
		set_temperature( tempering_->temperature() );
		set_sampling_temperature( tempering_->temperature() );
	}

	ig_->initialize( pose );

	for( core::Size ct = 1; ct <= ntrials_; ct++){
		core::Size resid;
		do { //pick residue, respecting underlying packer task
			resid = packed_residues()[Rg.random_range(1, packed_residues().size())];
		} while ( pose.residue( resid ).name1() == 'P'); //SidechainMover cannot sample proline rotamers? (is this true?)

		core::conformation::ResidueOP new_state( new core::conformation::Residue( pose.residue( resid ) ) );
		new_state = make_move( new_state );
		counters_.count_trial( type() );
		core::Real delta_energy = ig_->consider_substitution( resid, new_state, *new_state->nonconst_data_ptr() );
		if ( pass_metropolis( delta_energy, SidechainMover::last_proposal_density_ratio() ) ) { //ek
			previous_[ resid ] = current_[ resid ] ;
			ig_->commit_change( resid );
			current_energy -= delta_energy;
			current_[ resid ] = new_state;
			update_rotamers( resid );
			if( ( current_energy ) <  best_energy_ ) {
				best_[ resid ] = new_state;
				best_energy_ = current_energy;
			}
			counters_.count_accepted( type() );
		} else { //rejected metropolis criterion
			ig_->reject_change( resid, new_state, *new_state->nonconst_data_ptr() );
		}

		if ( tempering_ ) {
			set_temperature( tempering_->temperature_move( current_energy ) );
			set_sampling_temperature( tempering_->temperature() );
		}

		observe_rotamers( ct, traj_file_tag );

		core::Size model_count( output_count( ct ) );
		if ( model_count ) {
			for( core::Size res_i = 1; res_i <= current_.size(); res_i++ ){
				pose.replace_residue( res_i, (*current_[ res_i ]), true );
			}
			(*sfxn_)( pose );
			if ( job ) {
				job->add_string_real_pair( "prop_density", last_proposal_density_ratio() );
				job->add_string_string_pair( "move_type",type() );
			}
			jd2::output_intermediate_pose( pose, traj_file_tag, model_count, ( model_count % traj_stride_ ) != 0 && model_count > 1 ); //write always first a structure
		}
	} //for ntrials

	if ( tempering_ ) {
		tempering_->finalize_simulation( jd2::current_output_name() );
	}
	counters_.show( tr.Info );
	counters_.write_to_file( "trial_counts.stats", jd2::current_output_name() );
	score_post_apply_ = current_energy;

}

void CoupledSidechainProtocol::observe_rotamers( core::Size ct, std::string const& traj_file_tag ) {
	using namespace ObjexxFCL;

	if ( ct % rotamer_stride_ != 0 ) return;
	std::string filename( traj_file_tag+".rotamers");
	if (rotamer_stream_.filename() != filename) rotamer_stream_.open( filename );

	//format::F output
	Size const w( 8 );
	Size const d( 4 );
	for ( core::Size resi=1; resi<=chi_vectors_.size(); ++resi ) {
		rotamer_stream_ << format::I( w, resi ) << " ";
		for ( core::Size i=1; i<=4; ++i ) {
			if ( i<=chi_vectors_[ resi ].size() ) rotamer_stream_ << format::F( w, d, chi_vectors_[ resi ][i] ) << " ";
			else rotamer_stream_<< format::F( w, d, 0.0 ) << " ";
		}
		for ( core::Size i=1; i<=4; ++i ) {
			if ( i<=rot_vectors_[ resi ].size() ) rotamer_stream_ << format::I( w, rot_vectors_[ resi ][i] ) << " ";
			else rotamer_stream_<< format::I( w, 0.0 ) << " ";
		}
		rotamer_stream_ << format::I( w, ct );
		if ( tempering_ ) rotamer_stream_ << " " << format::F( w, d, tempering_->temperature() );
		rotamer_stream_	<< std::endl;
	}
}

std::string
CoupledSidechainProtocol::get_name() const {
	return "CoupledSidechainProtocol";
}

bool
CoupledSidechainProtocol::pass_metropolis(core::Real delta_energy , core::Real last_proposal_density_ratio ){

	core::Real boltz_factor = delta_energy / temperature_;
	core::Real probability = std::exp( std::min( 40.0, std::max( -40.0, boltz_factor ))) *  last_proposal_density_ratio ;
	if ( probability < 1 && numeric::random::rg().uniform() >= probability ) {
		current_ntrial_++;
		return false;
	} else {
		accepts_++;
		current_ntrial_++;
		return true;
	}
}

void
CoupledSidechainProtocol::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & pose) {
	ntrials_ = tag->getOption<core::Size>( "ntrials", 10000 );
	set_prob_uniform( tag->getOption<core::Real>( "prob_uniform", 0.0 ) );
	set_prob_withinrot( tag->getOption<core::Real>( "prob_withinrot", 0.0 ) );
	set_prob_random_pert_current( tag->getOption<core::Real>( "prob_random_pert_current", 0.0 ) );
	core::Real between_rot = 1.0 - prob_uniform() - prob_withinrot () - prob_random_pert_current();
	set_scorefunction(protocols::rosetta_scripts::parse_score_function(tag, data)->clone());
}

void
CoupledSidechainProtocol::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & metropolis_hastings_mover,
	core::Size cycle //default=0; non-zero if trajectory is restarted
)
{
	SidechainMover::initialize_simulation(pose, metropolis_hastings_mover);
}


} // moves
} // protocols
