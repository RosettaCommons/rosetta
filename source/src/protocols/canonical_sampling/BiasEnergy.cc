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
#include <utility/io/izstream.hh>

//#include <math.h>
//#include <iomanip>

static thread_local basic::Tracer tr( "protocols.canonical_sampling.BiasEnergy" );

namespace protocols {
namespace canonical_sampling {

using namespace core;

BiasEnergy::BiasEnergy() :
	bias_grid_( /* NULL */ ),
	count_grid_( /* NULL */ ),
	step_counter_( 0 ),
	stride_ ( 1 ),
	omega_ ( 1.0 ),
	gamma_ ( 5.0 ),
	deltaT_( 1.0 ),
	default_bias_( 1e4 )
{ }

BiasEnergy::BiasEnergy( core::Size stride, core::Real omega, core::Real gamma ) :
	bias_grid_( /* NULL */ ),
	count_grid_( /* NULL */ ),
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
		// Should be now OK with proper weak pointers: if BiasEnergy is not longer valid,
		// the AP lock will fail and the output observer will not be called... but not seg fault expected.
		protocols::jd2::JobDistributor::get_instance()->current_job()->remove_output_observer(
			utility::pointer::dynamic_pointer_cast< protocols::jd2::JobOutputterObserver >( get_self_ptr() )
		);
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
	runtime_assert( bias_grid_ != 0 );
	if ( protocols::jd2::jd2_used() ) {
		protocols::jd2::JobDistributor::get_instance()->current_job()->add_output_observer(
			utility::pointer::dynamic_pointer_cast< protocols::jd2::JobOutputterObserver >( get_self_ptr() )
		);
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

bool
BiasEnergy::restart_simulation(
	core::pose::Pose &,
	MetropolisHastingsMover & metropolis_hastings_mover,
	core::Size& cycle,
	core::Size&,
	core::Real& temperature
) {
	utility::io::izstream in( metropolis_hastings_mover.get_last_checkpoint()+".out" );
	tr.Debug << "restarting from checkpoint file " << in.filename() << std::endl;
	if ( !in.good() ) {
		tr.Error << "cannot open checkpoint file " << metropolis_hastings_mover.get_last_checkpoint() << ".out" << std::endl;
		return false;
	}

	std::string line;
	while ( getline( in, line ) ) {
		if ( line.substr(0, 17) != "REMARK BIASENERGY" ) continue;
		std::istringstream line_stream( line );
		std::string tag_remark, tag_key;
		line_stream >> tag_remark >> tag_key;
		if ( line_stream.good() ) {
			std::string tag, tag_hist, tag_grid_start, tag_grid_end;
			core::Real grid_min, grid_max;
			core::Size ngrid_cells;

			// read in WTE_Bias_Energy_Grid
			line_stream >> tag >> tag_hist >> ngrid_cells >> grid_min >> grid_max >> tag_grid_start;
			tr.Debug << tag <<" "<<tag_hist<<" "<<ngrid_cells<<" "<<grid_min<<" "<<grid_max<<" "<<tag_grid_start<< std::endl;
			if ( line_stream.fail() ) {
				tr.Debug << "bias energy grid is not complete in checkpoint file, restart from checkpoint abort!" << std::endl;
				return false;
			}
			for ( core::Size i=1; i<= bias_grid_->size(); ++i ) {
				line_stream >> ( bias_grid_->at(i) );
				if ( line_stream.fail() ) {
					tr.Debug << "bias energy grid is not complete in checkpoint file, restart from checkpoint abort!" << std::endl;
					return false;
				}
			}
			line_stream >> tag_grid_end;
			if ( line_stream.fail() ) {
				tr.Debug << "count grid is not complete in checkpoint file, restart from checkpoint abort!" << std::endl;
				return false;
			}

			// read in WTE_Count_Grid
			line_stream >> tag >> tag_hist >> ngrid_cells >> grid_min >> grid_max >> tag_grid_start;
			tr.Debug << tag <<" "<<tag_hist<<" "<<ngrid_cells<<" "<<grid_min<<" "<<grid_max<<" "<<tag_grid_start<< std::endl;
			if ( line_stream.fail() ) {
				tr.Debug << "count grid is not complete in checkpoint file, restart from checkpoint abort!" << std::endl;
				return false;
			}
			for ( core::Size i=1; i<= count_grid_->size(); ++i ) {
				line_stream >> ( count_grid_->at(i) );
				if ( line_stream.fail() ) {
					tr.Debug << "i= "<< i <<std::endl;
					tr.Debug << "count grid is not complete in checkpoint file, restart from checkpoint abort!" << std::endl;
					return false;
				}
			}
			std::string tag_end;
			line_stream >> tag_end;
			if ( line_stream.fail() || (tag_end != "GRID_END") ) {
				tr.Debug << "grid info format error, restart from checkpoint abort!" << std::endl;
				return false;
			}
		}
		step_counter_ = cycle;
		temperature_ = temperature;
		return true;
	}
	return false; // it should never reach here if complete info has been obtained from checkpoint file
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

void BiasEnergy::write_to_string( std::string& str ) const{
	/// current_output_name current_replica trial_count temp_level all in checkpoint_id
	str = "WTE_Bias_Energy_Grid ";
	bias_grid_->write_to_string( str );
	str = str + " WTE_Count_Grid ";
	count_grid_->write_to_string( str );
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
	bias_grid_ = HistogramFloatOP( new Histogram<float>( grid_min, grid_max, grid_size ) );
	count_grid_ = HistogramIntOP( new Histogram<int>( grid_min, grid_max, grid_size ) );
	grid_output_file_ = tag->getOption< std::string >( "wte_output", "wte_bias.grid" );
}
	///ntrials

} //moves
} //protocols

