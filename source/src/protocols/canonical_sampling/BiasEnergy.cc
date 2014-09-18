// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/canonical_sampling/BiasEnergyMover.cc
/// @brief BiasEnergy methods implemented
/// @author



// Unit Headers
#include <protocols/canonical_sampling/BiasEnergy.hh>
#include <protocols/canonical_sampling/BiasEnergy.tmpl.hh>

//for trajctory output of bias energy
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/util.hh>


#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/io/ozstream.hh>

//#include <math.h>
//#include <iomanip>

static thread_local basic::Tracer tr( "protocols.canonical_sampling.BiasEnergy" );

namespace protocols {
namespace canonical_sampling {

using namespace core;

BiasEnergy::BiasEnergy() :
	bias_grid_( NULL ),
	count_grid_( NULL ),
	step_counter_( 0 ),
	stride_ ( 1 ),
	omega_ ( 1.0 ),
	gamma_ ( 5.0 ),
	deltaT_( 1.0 ),
	default_bias_( 1e4 )
{ }

BiasEnergy::BiasEnergy( core::Size stride, core::Real omega, core::Real gamma ) :
	bias_grid_( NULL ),
	count_grid_( NULL ),
	step_counter_( 0 ),
	stride_ ( stride ),
	omega_ ( omega ),
	gamma_ ( gamma ),
	deltaT_( 1.0 ),
	default_bias_( 1e4 )
{ }

BiasEnergy::~BiasEnergy() {
	if ( protocols::jd2::jd2_used() ) {
		// I have to say, it makes me a bit nervous assuming that the BiasEnergy is going to be destroyed
		// before the current_job is.
		protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer( this );
	}
}

#ifdef USEMPI
core::Real BiasEnergy::update_and_evaluate_replica(
  int exchange_partner,
	MPI_Comm const& mpi_comm,
	core::pose::Pose const& pose )
{
	bias_grid_->mpi_exchange( exchange_partner, mpi_comm );
	count_grid_->mpi_exchange( exchange_partner, mpi_comm );
	core::Real const colvar = extract_collective_var( pose );
	if ( !bias_grid_->check_range( colvar ) ) return default_bias_;

	bias_grid_->toggle();
	core::Size index = bias_grid_->cell_index( colvar );
	core::Real score = bias_grid_->at( index );
	bias_grid_->toggle();

	return score;
}
#endif

void BiasEnergy::initialize_simulation(
	core::pose::Pose &,
	MetropolisHastingsMover const &,
	core::Size //non-zero if trajectory is restarted
) {
	runtime_assert( bias_grid_ );
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_output_observer( this );
	}
	bias_grid_->reset();
	count_grid_->reset();
}

void
BiasEnergy::finalize_simulation(
	core::pose::Pose &,
	MetropolisHastingsMover const &
) {
	tr.Trace << "writing grid information to file " << grid_output_file_ << std::endl;
	utility::io::ozstream os( grid_output_file_ );
	os << jd2::current_output_name() << " ";
	os << "WTE_Bias_Energy_Grid for temperature " << temperature_ << std::endl;
	bias_grid_->write_to_stream( os );
	os << "WTE_Count_Grid for temperature " << temperature_ << std::endl;
	count_grid_->write_to_stream( os );
}

void BiasEnergy::add_values_to_job( core::pose::Pose const & pose, protocols::jd2::Job & job ) const {
	job.add_string_real_pair( "bias", evaluate( pose ) );
}

core::Real BiasEnergy::evaluate( core::pose::Pose const& pose ) const {
	core::Real const colvar = extract_collective_var( pose );
	tr.Debug << "bias_energy evaluated at pos: " << colvar << std::endl;

	if ( !bias_grid_->check_range( colvar ) ) {
		return default_bias_;
	}

	core::Size const index = bias_grid_->cell_index( colvar );
	tr.Debug << " index: " << index << " bias: " << bias_grid_->at( index ) << std::endl;
	return bias_grid_->at( index );
}

void BiasEnergy::update( core::pose::Pose const& pose ) {
	step_counter_++;
	if ( step_counter_ % stride_ != 0 ) return;

	core::Real const colvar = extract_collective_var( pose );
	if ( !bias_grid_->check_range( colvar ) ) return;

	core::Size const index = bias_grid_->cell_index( colvar );
	core::Real const Vbias( bias_grid_->at( index ) );
	core::Real const height( omega_ * exp( -Vbias/deltaT_ ) * stride_ );
	tr.Debug << " update bias at index " << index << " pos: " << colvar << " add " << height << std::endl;
	bias_grid_->at( index ) += height;
	count_grid_->at( index ) += 1;
}

void BiasEnergy::set_temperature( core::Real val ) {
	deltaT_ = val * ( gamma_ - 1 );
	temperature_ = val;
}

void BiasEnergy::swap_replicas() {
	bias_grid_->toggle();
	count_grid_->toggle();
}

void
BiasEnergy::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	stride_ = tag->getOption< core::Size >( "wte_stride", 10 );
	omega_ = tag->getOption< core::Real >( "wte_omega", 1.0 );
	gamma_ = tag->getOption< core::Real >( "wte_gamma", 10.0 );
	core::Real const grid_min( tag->getOption< core::Real >( "wte_grid_min", -10 ) );
	core::Real const grid_max( tag->getOption< core::Real >( "wte_grid_max", 1000 ) );
	core::Size const grid_size( tag->getOption< core::Size >( "wte_grid_size", 100 ) );
	bias_grid_ = new Histogram<float>( grid_min, grid_max, grid_size );
	count_grid_ = new Histogram<int>( grid_min, grid_max, grid_size );
	grid_output_file_ = tag->getOption< std::string >( "wte_output", "wte_bias.grid" );
}
	///ntrials

} //moves
} //protocols

