// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/MetropolisHastingsMover.cc
/// @brief MetropolisHastingsMover methods implemented
/// @author


// Unit Headers
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMoverCreator.hh>


// protocols headers
#include <protocols/backrub/BackrubMover.hh>
#include <protocols/jd2/SilentFileJobOutputter.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>
#include <protocols/checkpoint/Checkpoint.hh>
#include <protocols/kinematic_closure/BalancedKicMover.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMover.hh>
#include <protocols/simple_moves/sidechain_moves/SidechainMCMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/canonical_sampling/ThermodynamicMover.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/TrajectoryRecorder.hh>
#include <protocols/canonical_sampling/BiasEnergy.hh>
#include <protocols/canonical_sampling/WTEBiasEnergy.hh>
#include <protocols/canonical_sampling/BiasedMonteCarlo.hh>

#include <protocols/rosetta_scripts/util.hh>

// core headers
#include <core/kinematics/MoveMap.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>
#include <core/pose/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentStruct.hh>
#include <basic/Tracer.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

#include <ObjexxFCL/format.hh>

//Auto Headers
#include <utility/excn/Exceptions.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

// C++ Headers
#include <vector>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.MetropolisHastingsMover" );

namespace protocols {
namespace canonical_sampling {

std::string
MetropolisHastingsMoverCreator::keyname() const {
	return MetropolisHastingsMoverCreator::mover_name();
}

protocols::moves::MoverOP
MetropolisHastingsMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new MetropolisHastingsMover );
}

std::string
MetropolisHastingsMoverCreator::mover_name() {
	return "MetropolisHastings";
}

MetropolisHastingsMover::MetropolisHastingsMover() :
	monte_carlo_(/*0*/),
	ntrials_(1000),
	checkpoint_count_(0)
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ run::checkpoint_interval ].user() ) {
		protocols::checkpoint::checkpoint_with_interval( option[ run::checkpoint_interval ] );
		job_outputter_ = jd2::SilentFileJobOutputterOP( new jd2::SilentFileJobOutputter() ); // only required when writing out checkpoints
	}
}

MetropolisHastingsMover::MetropolisHastingsMover(
	MetropolisHastingsMover const & metropolis_hastings_mover
) :
	//utility::pointer::ReferenceCount(),
	Mover(metropolis_hastings_mover),
	monte_carlo_( metropolis_hastings_mover.monte_carlo_->clone() ),
	ntrials_(metropolis_hastings_mover.ntrials_),
	output_name_(metropolis_hastings_mover.output_name_),
	weighted_sampler_(metropolis_hastings_mover.weighted_sampler_),
	job_outputter_(metropolis_hastings_mover.job_outputter_),
	checkpoint_count_(metropolis_hastings_mover.checkpoint_count_)
{
	for ( core::Size i = 1; i <= metropolis_hastings_mover.movers_.size(); ++i ) {
		movers_.push_back(utility::pointer::dynamic_pointer_cast<ThermodynamicMover>(metropolis_hastings_mover.movers_[i]));
	}

	for ( core::Size i = 1; i <= metropolis_hastings_mover.observers_.size(); ++i ) {
		observers_.push_back(utility::pointer::dynamic_pointer_cast<ThermodynamicObserver>(metropolis_hastings_mover.observers_[i]));
	}

	if ( metropolis_hastings_mover.tempering_ ) {
		tempering_ = utility::pointer::dynamic_pointer_cast<TemperatureController>(metropolis_hastings_mover.tempering_);
		if ( monte_carlo_ ) tempering_->set_monte_carlo(monte_carlo_);
	}

}

MetropolisHastingsMover::~MetropolisHastingsMover()= default;

/// @details The return value indicates the number of cycles that have already
/// been run, if the simulation is not being started or restarted.  I'm not
/// totally sure what this means though, and I couldn't see any way for this
/// function to return anything other than 0.  The necessary logic might be
/// commented out right now.
core::Size
MetropolisHastingsMover::prepare_simulation( core::pose::Pose & pose ) {
	if ( output_name() == "" ) {
		set_output_name(protocols::jd2::JobDistributor::get_instance()->current_output_name());
		tr.Info  << " obtained output name from JobDistributor " << std::endl;
		output_name_from_job_distributor_ = true;
	} else {
		tr.Info << " running with preset output name: " << output_name() << std::endl;
	}

	if ( !tempering_ ) {
		//get this done before "initialize_simulation" is called no movers and observers
		tr.Info << "no temperature controller in MetropolisHastings defined... generating FixedTemperatureController" << std::endl;
		tempering_ = TemperatureControllerOP( new protocols::canonical_sampling::FixedTemperatureController( monte_carlo_->temperature() ) );
		tempering_->set_monte_carlo( monte_carlo_ ); // MonteCarlo required for tempering_->observe_after_metropolis();
	}
	using namespace core;
	bool restart = false;
	core::Size cycle_number = 0;
	Size temp_level = 0;
	Real temperature = -1.0;
	//  for (core::Size i = 1; i <= observers_.size() && !restart; ++i) {
	//   tr.Info<< "Attempting restart using " << observers_[i]->get_name() << std::endl;
	//   restart = observers_[i]->restart_simulation(pose, *this, cycle_number, temp_level, temperature );
	//   if ( restart ) tr.Info<< "Restarted using " << observers_[i]->get_name() << std::endl;
	//  }

	if ( get_checkpoints() ) {
		for ( core::Size i = 1; i <= observers_.size() && !restart; ++i ) {
			if ( utility::pointer::dynamic_pointer_cast< TrajectoryRecorder >( observers_[i] ) ) { // first get the silent struct and cycle info
				tr.Info<< "Attempting restart using " << observers_[i]->get_name() << std::endl;
				restart = observers_[i]->restart_simulation(pose, *this, cycle_number, temp_level, temperature );
				if ( restart ) tr.Info<< "Restarted using " << observers_[i]->get_name() << " cycle_number " << cycle_number << " temp_level " << temp_level << std::endl;
			}
		}

		for ( core::Size i = 1; i <= observers_.size() && restart; ++i ) { // if restart-able from silent trajectory file, collect BiasEnergy info
			BiasEnergyOP bias_energy( utility::pointer::dynamic_pointer_cast< BiasEnergy >( observers_[i] ));
			if ( bias_energy ) {
				tr.Debug << "bias energy used " << std::endl;
				restart = bias_energy->restart_simulation( pose, *this, cycle_number, temp_level, temperature );
				if ( restart ) {
					tr.Info << "grid info collected, restart now!" << std::endl;
				} else {
					tr.Info << "failed collecting grid info, restart failed!" << std::endl;
				}
			}
		}
	}

	if ( !restart ) {
		cycle_number = 0; //make sure this is zero if we don't have a restart.
		tempering_->initialize_simulation(pose, *this, cycle_number );
	} else {
		tempering_->initialize_simulation(pose, *this, temp_level, temperature, cycle_number );
	}

	for ( core::Size i = 1; i <= movers_.size(); ++i ) {
		tr.Info << "Initializing " << movers_[i]->get_name() << std::endl;
		movers_[i]->set_metropolis_hastings_mover(
			MetropolisHastingsMoverAP( utility::pointer::static_pointer_cast< MetropolisHastingsMover >( get_self_ptr() ) )
		);
		movers_[i]->initialize_simulation(pose, *this, cycle_number);
	}

	runtime_assert( monte_carlo_ != nullptr );
	monte_carlo_->reset(pose);
	monte_carlo_->reset_counters();

	for ( core::Size i = 1; i <= observers_.size(); ++i ) {
		tr.Info << "Initializing " << observers_[i]->get_name() << std::endl;
		observers_[i]->initialize_simulation(pose, *this, cycle_number);
	}

	tr.Info << "Initial Score:\n";
	monte_carlo_->score_function().show(tr, pose);

	tr.Info << "Running " << ntrials_-cycle_number << " trials..." << std::endl;
	return cycle_number;
}

void
MetropolisHastingsMover::apply( core::pose::Pose& pose ) {
	output_name_from_job_distributor_ = false;

	Size start_cycle = prepare_simulation( pose );

	for ( current_trial_ = start_cycle+1; !tempering_->finished_simulation( current_trial_, ntrials() ); ++current_trial_ ) {
		write_checkpoint( pose );
		ThermodynamicMoverOP mover(random_mover());
		mover->apply(pose);
		bool accepted = monte_carlo_->boltzmann(
			pose,
			mover->type(),
			mover->last_proposal_density_ratio(),
			mover->last_inner_score_delta_over_temperature()
		);
		set_last_move( mover );
		set_last_accepted( accepted );
		tempering_->temperature_move( pose );
		mover->observe_after_metropolis(*this);
		tr.Trace << "current move accepted " << accepted <<std::endl;
		for ( core::Size i = 1; i <= observers_.size(); ++i ) {
			observers_[i]->observe_after_metropolis(*this);
		}
		//currently used to write tempering stats every so often...
		tempering_->observe_after_metropolis(*this);
	}
	wind_down_simulation( pose );
}

void
MetropolisHastingsMover::write_checkpoint( core::pose::Pose const & pose ) {
	using namespace ObjexxFCL;
	using namespace protocols::checkpoint;

	if ( !Timer::is_on() ) return;
	if ( !Timer::time_to_checkpoint() ) return;

	checkpoint_count_++;

	std::ostringstream checkpoint_id;
	// checkpoint_id << jd2::current_output_filename(); // decoys.out
	utility::file::FileName jd2_filename( jd2::current_output_filename() );
	checkpoint_id << jd2_filename.base(); // decoys
	checkpoint_id << "_" << jd2::current_output_name(); // protAB_0001
	checkpoint_id << "_" << jd2::current_replica(); // rep
	checkpoint_id << "_" << checkpoint_count_;

	checkpoint_ids_.push_back( checkpoint_id.str() ); // decoys_protAB_0001_1_n
	utility::io::ozstream out( checkpoint_id.str()+".out");
	// std::string const & checkpoint_id( jd2::current_output_name() + "_" + string_of( jd2::current_replica() ) + "_" + string_of( checkpoint_count_ ) ); // protAB_0001_3_N
	// utility::file::FileName jd2_filename( jd2::current_output_filename() );
	// checkpoint_ids_.push_back( jd2_filename.base()+"_"+checkpoint_id );
	// tr.Debug << "checkpoint_id: " << checkpoint_id << std::endl;  // "decoys_protAB_0001_3_n"
	// core::pose::Pose tmp_pose( pose );

	core::io::silent::SilentStructOP pss( new core::io::silent::BinarySilentStruct() );
	std::ostringstream tag;
	tag << jd2::current_output_name();
	tag << "_" << std::setfill('0') << std::setw(3) << jd2::current_replica();
	tag << "_" << std::setfill('0') << std::setw(8) << current_trial_;

	pss->fill_struct( pose, tag.str() );

	BiasedMonteCarloOP biased_mc = utility::pointer::dynamic_pointer_cast< protocols::canonical_sampling::BiasedMonteCarlo > ( monte_carlo_ );
	if ( biased_mc ) { // if BiasedMonteCarlo, write out bias grid info
		std::string str="";
		biased_mc->bias_energy()->write_to_string( str );
		//  core::pose::Pose tmp_pose( pose );
		//  core::pose::add_comment( tmp_pose, "BIASENERGY", str );
		pss->add_comment("BIASENERGY", str);
	}
	// write out snapshots for the replica, with joboutputter, the comment in pose will be not outputted
	//  job_outputter_->other_pose( jd2::get_current_job(), tmp_pose, checkpoint_id, current_trial_, false );

	protocols::jd2::add_job_data_to_ss( pss, jd2::get_current_job() );
	pss->print_header( out );
	pss->print_scores( out );
	pss->print_conformation( out );

	Timer::reset();
	// tr.Debug << "have done writing checkpoint " << checkpoint_ids_.back() << std::endl;

	/// maintain 5 snapshots for each replica, delete the rest
	if ( checkpoint_ids_.size() > 5 ) { ///keep 5 snapshots for each replica
		//  tr.Debug << "deleting old checkpoints, only keep the last 5... " << std::endl;
		utility::file::file_delete( checkpoint_ids_.front()+".out" );
		checkpoint_ids_.erase( checkpoint_ids_.begin() ); ///each time delete one
	}
}

bool
MetropolisHastingsMover::get_checkpoints() {
	using namespace utility::file;

	utility::file::FileName jd2_filename( jd2::current_output_filename() );
	std::string filename_base( jd2_filename.base()+"_"+jd2::current_output_name() );
	std::string filename_pattern( filename_base+"_"+ObjexxFCL::string_of( jd2::current_replica()) +"_"); // current_replica()=0 for FixedTemp

	utility::vector1< std::string > names;
	std::vector< int > checkpoint_indics;
	utility::file::list_dir( ".", names );
	for ( core::Size i=1; i<=names.size(); ++i ) {
		if ( names[i].find( filename_pattern ) != std::string::npos ) {
			FileName found_name( names[i] );
			core::Size ind = found_name.base().find_last_of('_');
			checkpoint_indics.push_back( utility::string2int( found_name.base().substr(ind+1) ) );
		}
	}
	tr.Debug << "found " << checkpoint_indics.size() << " checkpoints" << std::endl;
	if ( checkpoint_indics.size()==0 ) return false; // no checkpoint found
	std::sort( checkpoint_indics.begin(), checkpoint_indics.end() ); // sort by checkpoint index

	// starting from the last checkpoint, if all temp_levels are collected, then return
	auto rit = checkpoint_indics.rbegin();
	for ( ; rit != checkpoint_indics.rend(); ++rit ) {
		tr.Debug << "checkpoint_indics: " << *rit << std::endl;

		if ( tempering_->n_temp_levels() == 1 ) { // FixedTemperatureController
			core::io::silent::SilentFileData sfd( filename_base+"_"+ObjexxFCL::string_of( 0 )+"_"+ObjexxFCL::string_of( *rit )+".out");
			if ( utility::file::file_exists( sfd.filename() ) && utility::file::file_size( sfd.filename() ) ) {
				checkpoint_ids_.push_back( filename_pattern + ObjexxFCL::string_of( *rit ) );
				checkpoint_count_ = *rit;

				// delete any checkpoint file before the successfully collected checkpoint
				auto d_rit = ++rit;
				for ( ; d_rit != checkpoint_indics.rend(); ++d_rit ) utility::file::file_delete( filename_pattern+ObjexxFCL::string_of( *d_rit )+".out" );
				tr.Debug << "checkpoint for restart: " << checkpoint_ids_.front() << std::endl;
				return true;
			} else continue;

		} else { // replica exchange

			utility::vector1< core::Size > found_levels;
			// collect all the temp_levels from the file of this checkpoint
			for ( core::Size replica=1; replica <= tempering_->n_temp_levels(); ++replica ) {
				core::io::silent::SilentFileData sfd( filename_base+"_"+ObjexxFCL::string_of( replica )+"_"+ObjexxFCL::string_of( *rit )+".out");
				if ( utility::file::file_exists( sfd.filename() ) ) {
					sfd.read_file( sfd.filename() );
					found_levels.push_back( sfd.begin()->get_energy( "temp_level" ));
				} else break;
			}
			// check if any temp_level is missing, if not, then return.
			for ( core::Size replica=1; replica <= tempering_->n_temp_levels(); ++replica ) {
				if ( !found_levels.has_value( replica ) ) { // certain level missing
					tr.Debug << "temp_level: " << replica << " is missing" << std::endl;
					break;
				}
				if ( replica == tempering_->n_temp_levels() ) {
					checkpoint_ids_.push_back( filename_pattern + ObjexxFCL::string_of( *rit ) );
					checkpoint_count_ = *rit;
					// delete any checkpoint file before the successfully collected checkpoint
					auto d_rit = ++rit;
					for ( ; d_rit != checkpoint_indics.rend(); ++d_rit ) utility::file::file_delete( filename_pattern+ObjexxFCL::string_of( *d_rit )+".out" );
					tr.Debug << "checkpoint for restart: " << checkpoint_ids_.front() << std::endl;
					return true;
				}
			}
		}
	}
	std::cout << "complete restart info not collected from checkpoint, start from the begining..." << std::endl;
	return false; // if runs to here, then not found for sure;
}

std::string
MetropolisHastingsMover::get_last_checkpoint() const{
	if ( checkpoint_ids_.size()==0 ) {
		return "";
	}
	return checkpoint_ids_.back();
}

void
MetropolisHastingsMover::wind_down_simulation( core::pose::Pose& pose) {
	for ( core::Size i = 1; i <= movers_.size(); ++i ) {
		tr.Info << "Finalizing " << movers_[i]->get_name() << std::endl;
		movers_[i]->finalize_simulation(pose, *this);
	}

	for ( core::Size i = 1; i <= observers_.size(); ++i ) {
		tr.Info << "Finalizing " << observers_[i]->get_name() << std::endl;
		observers_[i]->finalize_simulation(pose, *this);
	}

	tempering_->finalize_simulation(pose, *this);

	tr.Debug << "Final Score:" << std::endl;
	monte_carlo_->score_function().show(tr.Debug, pose);
	tr.Debug << std::endl;

	monte_carlo_->show_counters();

	if ( output_name_from_job_distributor_ ) set_output_name("");
}

std::string
MetropolisHastingsMover::get_name() const
{
	return "MetropolisHastingsMover";
}

protocols::moves::MoverOP
MetropolisHastingsMover::clone() const
{
	return protocols::moves::MoverOP( new protocols::canonical_sampling::MetropolisHastingsMover(*this) );
}

protocols::moves::MoverOP
MetropolisHastingsMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new MetropolisHastingsMover );
}

bool
MetropolisHastingsMover::reinitialize_for_each_job() const { return false; }

bool
MetropolisHastingsMover::reinitialize_for_new_input() const { return false; }

void
MetropolisHastingsMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {
	///ntrials
	ntrials_ = tag->getOption< core::Size >( "trials", ntrials_ );

	///checkpoint interval
	if ( tag->hasOption("checkpoint_interval") ) {
		core::Size const checkpoint_interval( tag->getOption< core::Size >( "checkpoint_interval")); // every 6 minutes
		protocols::checkpoint::checkpoint_with_interval( checkpoint_interval );
		job_outputter_ = jd2::SilentFileJobOutputterOP( new jd2::SilentFileJobOutputter() );
	}

	//monte-carlo
	//read tag scorefxn
	core::scoring::ScoreFunctionOP score_fxn(protocols::rosetta_scripts::parse_score_function(tag, data));
	core::Real const temperature(tag->getOption< core::Real >( "temperature", 0.6 ) );

	bool wte_sampling( tag->getOption< bool >( "wte", false ) );
	if ( wte_sampling ) {
		BiasEnergyOP bias_energy( new WTEBiasEnergy() ); // stride, omega, gamma );
		bias_energy->parse_my_tag( tag, data, filters, movers, pose );
		add_observer( bias_energy );
		monte_carlo_ = protocols::moves::MonteCarloOP( new BiasedMonteCarlo( *score_fxn, temperature, bias_energy ) );
	} else {
		monte_carlo_ = protocols::moves::MonteCarloOP( new protocols::moves::MonteCarlo( *score_fxn, temperature ) );
	}

	//add movers and observers
	utility::vector0< utility::tag::TagCOP > const subtags( tag->getTags() );
	for (auto subtag : subtags) {
			protocols::moves::MoverOP mover;
		if ( subtag->getName() == "Add" ) { //add existing mover
			std::string mover_name = subtag->getOption<std::string>( "mover_name", "null" );
			auto mover_iter( movers.find( mover_name ) );
			if ( mover_iter == movers.end() ) {
				tr.Error<< "Mover not found for XML tag:\n" << subtag << std::endl;
				throw utility::excn::EXCN_RosettaScriptsOption("");
			}
			mover = mover_iter->second;
		} else { //generate new mover
			protocols::moves::MoverFactory * mover_factory(protocols::moves::MoverFactory::get_instance());
			mover = mover_factory->newMover(subtag, data, filters, movers, pose);
		}

		ThermodynamicMoverOP th_mover( utility::pointer::dynamic_pointer_cast< protocols::canonical_sampling::ThermodynamicMover > ( mover ) );
		ThermodynamicObserverOP th_observer( utility::pointer::dynamic_pointer_cast< protocols::canonical_sampling::ThermodynamicObserver > ( mover ) );
		TemperatureControllerOP temp_controller( utility::pointer::dynamic_pointer_cast< protocols::canonical_sampling::TemperatureController > ( mover ) );
		//figure out if ThermodynamicMover or ThermodynamicObserver
		if ( th_mover ) { //its a mover
			core::Real const weight( subtag->getOption< core::Real >( "sampling_weight", 1 ) );
			add_mover( th_mover, weight, subtag );
			//    add_mover( th_mover, weight );
		} else if ( th_observer ) { //its an observer
			//it might also be a tempering module...
			if ( temp_controller ) { // it is a temperature controller
				if ( tempering_ ) {
					throw utility::excn::EXCN_RosettaScriptsOption( "cannot define two TemperatureControllers" );
				}
				set_tempering( temp_controller );
			} else {  //no just an plain old observer
				add_observer( th_observer );
			}
		} else { //its something different
			tr.Error<< "Mover is not a ThermodynamicMover or ThermodynamicObserver for XML tag:\n" << subtag << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption("");
		}
	}
}

protocols::moves::MonteCarloCOP
MetropolisHastingsMover::monte_carlo() const {
	return monte_carlo_;
}

protocols::moves::MonteCarlo& MetropolisHastingsMover::nonconst_monte_carlo() {
	return *monte_carlo_;
}

void
MetropolisHastingsMover::set_monte_carlo(
	protocols::moves::MonteCarloOP monte_carlo
)
{
	monte_carlo_ = monte_carlo;
	if ( tempering_ ) tempering_->set_monte_carlo( monte_carlo_ );
}

void
MetropolisHastingsMover::set_tempering(
	TemperatureControllerOP tempering
) {
	tempering_=tempering;
	if ( monte_carlo_ ) tempering_->set_monte_carlo( monte_carlo_ );
}

TemperatureControllerCOP
MetropolisHastingsMover::tempering() const {
	return tempering_;
}

void
MetropolisHastingsMover::set_ntrials(
	core::Size ntrials
) {
	ntrials_ = ntrials;
}

std::string const &
MetropolisHastingsMover::output_name() const {
	return output_name_;
}

void
MetropolisHastingsMover::set_output_name(
	std::string const & output_name
) {
	output_name_ = output_name;
}

std::string
MetropolisHastingsMover::output_file_name(
	std::string const & suffix,
	bool cumulate_jobs, //= false
	bool cumulate_replicas //= false
) const {

	if ( cumulate_jobs || cumulate_replicas ) runtime_assert(utility::io::ozstream::MPI_reroute_rank() >= 0);

	std::ostringstream file_name_stream;

	if ( !cumulate_jobs ) {
		file_name_stream << output_name_;
	}

	core::Size const replica( protocols::jd2::current_replica() );
	if ( !cumulate_replicas && replica ) {
		if ( file_name_stream.str().length() ) file_name_stream << "_";
		file_name_stream << std::setfill('0') << std::setw(3) << replica;
	}

	if ( suffix.length() ) {
		if ( file_name_stream.str().length() ) file_name_stream << "_";
		file_name_stream << suffix;
	}

	return file_name_stream.str();
}

void
MetropolisHastingsMover::set_last_move(
	ThermodynamicMoverOP setting
) {
	last_move_ = setting;
}

ThermodynamicMover const&
MetropolisHastingsMover::last_move() const {
	assert( last_move_ );
	return *last_move_;
}

/// @details You can control the probability of selecting a particular mover
/// via the weight argument to the @ref add_mover() method.
ThermodynamicMoverOP
MetropolisHastingsMover::random_mover() const {
	return movers_[weighted_sampler_.random_sample(numeric::random::rg())];
}

/// @details Specify a weight to control how often this mover should be
/// invoked.  The weight does not have to be in any particular range.  The
/// probability of choosing any particular move will be the weight of that
/// sampler divided by the sum of the weights of all the moves.
void
MetropolisHastingsMover::add_mover(
	ThermodynamicMoverOP mover,
	core::Real weight
) {
	mover->set_preserve_detailed_balance(true);
	movers_.push_back(mover);
	weighted_sampler_.add_weight(weight);
}

/// @details In principle, information about the mover could be extracted from
/// the given XML tag, but currently this function is a simple alias for @ref
/// add_mover().
void
MetropolisHastingsMover::add_mover(
	ThermodynamicMoverOP mover,
	core::Real weight,
	utility::tag::TagCOP const&
) {
	add_mover( mover, weight );
}

void
MetropolisHastingsMover::add_backrub_mover(
	core::Real weight
)
{
	if ( !weight ) return;

	protocols::backrub::BackrubMoverOP backrub_mover( new protocols::backrub::BackrubMover );
	backrub_mover->branchopt().read_database();

	add_mover(backrub_mover, weight);
}

void
MetropolisHastingsMover::add_kic_mover(
	core::Real weight,
	protocols::loops::Loop const & loop
)
{
	if ( !weight ) return;

	protocols::kinematic_closure::BalancedKicMoverOP kic_mover( new protocols::kinematic_closure::BalancedKicMover );
	kic_mover->set_loop(loop);

	add_mover(kic_mover, weight);
}

void
MetropolisHastingsMover::add_small_mover(
	core::Real weight
)
{
	if ( !weight ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::SmallMoverOP small_mover( new protocols::simple_moves::SmallMover );
	small_mover->nmoves(1);
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	if ( utility::file::file_exists(option[ in::file::movemap ]) ) {
		movemap->init_from_file(option[ in::file::movemap ]);
	} else {
		movemap->set_bb(true);
	}
	small_mover->movemap(movemap);

	add_mover(small_mover, weight);
}

void
MetropolisHastingsMover::add_shear_mover(
	core::Real weight
)
{
	if ( !weight ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	protocols::simple_moves::ShearMoverOP shear_mover( new protocols::simple_moves::ShearMover );
	shear_mover->nmoves(1);
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	if ( utility::file::file_exists(option[ in::file::movemap ]) ) {
		movemap->init_from_file(option[ in::file::movemap ]);
	} else {
		movemap->set_bb(true);
	}
	shear_mover->movemap(movemap);

	add_mover(shear_mover, weight);
}

void
MetropolisHastingsMover::add_sidechain_mover(
	core::Real weight,
	core::Real prob_uniform,
	core::Real prob_withinrot,
	bool preserve_cbeta
)
{
	if ( !weight ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::pack::task::operation::TaskOperationCOP;

	core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
	main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );
	} else {
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	}
	if ( preserve_cbeta ) main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::PreserveCBeta ) );

	protocols::simple_moves::sidechain_moves::SidechainMoverOP sidechain_mover( new protocols::simple_moves::sidechain_moves::SidechainMover );
	sidechain_mover->set_task_factory(main_task_factory);
	sidechain_mover->set_prob_uniform(prob_uniform);
	sidechain_mover->set_prob_withinrot(prob_withinrot);

	add_mover(sidechain_mover, weight);
}

void
MetropolisHastingsMover::add_sidechain_mc_mover(
	core::Real weight,
	core::Real prob_uniform,
	core::Real prob_withinrot,
	bool preserve_cbeta,
	core::Size ntrials
)
{
	if ( !weight ) return;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using core::pack::task::operation::TaskOperationCOP;

	core::pack::task::TaskFactoryOP main_task_factory( new core::pack::task::TaskFactory );
	main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::ReadResfile ) );
	} else {
		main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking ) );
	}
	if ( preserve_cbeta ) main_task_factory->push_back( TaskOperationCOP( new core::pack::task::operation::PreserveCBeta ) );

	protocols::simple_moves::sidechain_moves::SidechainMCMoverOP sidechain_mc_mover( new protocols::simple_moves::sidechain_moves::SidechainMCMover );
	sidechain_mc_mover->set_task_factory(main_task_factory);
	sidechain_mc_mover->set_prob_uniform(prob_uniform);
	sidechain_mc_mover->set_prob_withinrot(prob_withinrot);
	sidechain_mc_mover->setup(monte_carlo_->score_function().clone());
	sidechain_mc_mover->set_ntrials(ntrials);

	add_mover(sidechain_mc_mover, weight);
}

void
MetropolisHastingsMover::add_observer(
	ThermodynamicObserverOP observer
)
{
	observers_.push_back(observer);
}

} //moves
} //protocols

