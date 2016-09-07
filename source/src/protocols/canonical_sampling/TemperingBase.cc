// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/canonical_sampling/TemperingBaseMover.cc
/// @brief TemperingBase methods implemented
/// @author Oliver Lange ( oliver.lange@tum.de )


// Unit Headers
#include <protocols/canonical_sampling/TemperingBase.hh>

// protocols headers
#include <basic/datacache/DataMap.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/canonical_sampling/ThermodynamicObserver.hh>
#include <protocols/canonical_sampling/MetropolisHastingsMover.hh>

#include <protocols/rosetta_scripts/util.hh>

//#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/jd2/Job.hh>

// core headers
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>

// numeric headers
#include <numeric/random/random.hh>

// utility headers
#include <utility/file/file_sys_util.hh>
#include <utility/string_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/tag/Tag.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

// C++ Headers
#include <cmath>

// cmd-line options
OPT_2GRP_KEY( File, tempering, temp, file )
OPT_2GRP_KEY( Integer, tempering, temp, levels )
OPT_2GRP_KEY( Real, tempering, temp, low )
OPT_2GRP_KEY( Real, tempering, temp, high )
OPT_2GRP_KEY( RealVector, tempering, temp, range )
OPT_2GRP_KEY( File, tempering, stats, file )
OPT_2GRP_KEY( Boolean, tempering, stats, silent )
OPT_2GRP_KEY( Boolean, tempering, stats, line_output )
OPT_1GRP_KEY( Integer, tempering, stride )

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer tr( "protocols.canonical_sampling.TemperingBase" );


bool protocols::canonical_sampling::TemperingBase::options_registered_( false );

//Mike: when you want to remove these Macros... leave them at least here as comment - since they provide documentation
void protocols::canonical_sampling::TemperingBase::register_options() {
	if ( !options_registered_ ) {
		options_registered_ = true;
		NEW_OPT( tempering::temp::file, "file with temperature definitions for simulated tempering, including weights and counts","" );
		NEW_OPT( tempering::temp::levels, "how many temp-levels, superseeded by -tempering::temp::file",10 );
		NEW_OPT( tempering::temp::range, "min and max temperature, superseeded by -tempering::temp::file", 0 );
		NEW_OPT( tempering::temp::high, "min and max temperature, superseeded by -tempering::temperatures or -tempering::temp::range", 3.0 );
		NEW_OPT( tempering::temp::low, "min and max temperature, superseeded by -tempering::temperatures or -tempering::temp::range", 0.5 );

		NEW_OPT( tempering::stats::file, "filename (postfix) for output of tempering statistics (i.e, counts) <job>_.tempering.stats", "tempering.stats" );
		NEW_OPT( tempering::stats::silent, "write all tempering information as lines into a single file -- similar to silent format", false );
		NEW_OPT( tempering::stats::line_output, "choose line-output as in silent mode, even when using individual files", false );

		NEW_OPT( tempering::stride, "how often should a temperature switch be attempted", 1);
	}
}

namespace protocols {
namespace canonical_sampling {
using namespace core;

TemperingBase::TemperingBase() :
	instance_initialized_( false )
{
	set_defaults();
}

TemperingBase::TemperingBase( TemperingBase const & other ) :
	protocols::canonical_sampling::TemperatureController(other)
{
	temperatures_ = other.temperatures_;
	temperature_stride_ = other.temperature_stride_;
	trust_current_temp_ = other.trust_current_temp_ ;
	stats_line_output_ = other.stats_line_output_;
	stats_silent_output_ = other.stats_silent_output_;
	stats_file_ = other.stats_file_;

	job_ = other.job_;
	instance_initialized_ = other.instance_initialized_;
	current_temp_ = other.current_temp_;
	temp_trial_count_ = other.temp_trial_count_;
}

Size TemperingBase::n_temp_levels() const { return temperatures_.size(); }

core::Real
TemperingBase::temperature() const {
	return monte_carlo()->temperature();
}

core::Real
TemperingBase::temperature( Size level ) const {
	return temperatures_[ level ];
}

std::string
TemperingBase::get_name() const
{
	return "TemperingBase";
}

void
TemperingBase::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	pose::Pose const &
) {
	//figure out temperatures...
	std::string temp_file = tag->getOption< std::string >( "temp_file", "" );
	bool success( false );
	if ( temp_file.size() ) {
		success=initialize_from_file( temp_file );
		if ( !success ) tr.Info << "cannot read temperatures from file, will initialize from options... " << std::endl;
	}
	if ( !success ) {
		Real temp_low = tag->getOption< Real >( "temp_low", 0.6 );
		Real temp_high = tag->getOption< Real >( "temp_high", 3.0 );
		Size temp_levels = tag->getOption< Size >( "temp_levels", 10 );
		InterpolationType temp_interpolation = interpolation_type_string_to_enum(tag->getOption< std::string >( "temp_interpolation", "linear" ));
		generate_temp_range( temp_low, temp_high, temp_levels, temp_interpolation );
	}

	//simple options
	temperature_stride_ = tag->getOption< Size >( "temp_stride", 100 );
	io_stride_ = tag->getOption< Size >("io_stride", 10000 );
	trust_current_temp_ = tag->getOption< bool >( "trust_crurrent_temp", true );
	stats_line_output_ = tag->getOption< bool >( "stats_line_output", false );
	stats_silent_output_ = tag->getOption< bool >( "stats_silent_output", false );
	stats_file_ = tag->getOption< std::string >( "stats_file", "tempering.stats" );
	instance_initialized_ = true;
}


/// handling of options including command-line
void TemperingBase::set_defaults() {
	//are we the only object controlling temperature in the MC object ?!
	trust_current_temp_ = true;
	//how often will no temperature jump occur
}

void TemperingBase::initialize_simulation(
	pose::Pose&,
	protocols::canonical_sampling::MetropolisHastingsMover const &,
	core::Size //default=0; non-zero if trajectory is restarted
) {
	tr.Trace << "initialize Tempering Base... " << std::endl;
	if ( !instance_initialized_ ) init_from_options();
	current_temp_=temperatures_.size();
	monte_carlo()->set_temperature( temperatures_[ current_temp_ ] );
	if ( jd2::jd2_used() ) {
		job_ = jd2::get_current_job();
	}
	tr.Debug << std::setprecision(2);
	if ( job_ ) {
		job_->add_string_real_pair( "temperature", monte_carlo()->temperature() );
		job_->add_string_real_pair( "temp_level", current_temp_ );
	}
	temp_trial_count_ = 0;
	tr.Trace << "initialized Tempering Base!!!" << std::endl;
	trial_counter_.set_temperature_observer( this );
}

void
TemperingBase::initialize_simulation(
	core::pose::Pose & pose,
	protocols::canonical_sampling::MetropolisHastingsMover const & mh_mover,
	core::Size level,
	core::Real temp_in,
	core::Size cycle //default=0; non-zero if trajectory is restarted
) {
	initialize_simulation( pose, mh_mover,cycle );
	set_current_temp( level );
	if ( std::abs( temperature() - temp_in ) > 0.01 ) {
		using namespace ObjexxFCL;
		utility_exit_with_message( "inconsistent temperature when trying to restart Tempering temp_in: "
			+string_of( temp_in )
			+" temperature expected at level "+string_of( level )+ " is "+string_of( temperature() ) );
	}
}

void
TemperingBase::observe_after_metropolis( MetropolisHastingsMover const & mhm ) {
	if ( mhm.current_trial() && mhm.current_trial() % io_stride_ == 0 ) {
#ifdef WIN32
		core::Size output_ct( floor( (double)mhm.current_trial() / (double)io_stride_ ) );
#else
		core::Size output_ct( floor( mhm.current_trial() / io_stride_ ) );
#endif
		trial_counter_.write_to_file( stats_file_, mhm.output_name() + utility::to_string( output_ct ) );
	}
}

void
TemperingBase::finalize_simulation(
	pose::Pose&,
	protocols::canonical_sampling::MetropolisHastingsMover const & mhm
) {
	//std::string tag( "no_tag" );
	tr.Trace << "write statistics to " << stats_file_ << "..." << std::endl;
	trial_counter_.write_to_file( stats_file_, mhm.output_name() );
	tr.Trace << "done" << stats_file_ << std::endl;
	job_ = nullptr;
}

bool TemperingBase::check_temp_consistency() {
	//we can trust our current temp variable or get it from protocols::moves::MonteCarlo object
	if ( !trust_current_temp_ ) {
		current_temp_=1;
		Real const mc_temp( monte_carlo()->temperature() );
		for ( utility::vector1< Real >::const_iterator it = temperatures_.begin();
				it != temperatures_.end(); ++it ) {
			if ( *it > mc_temp ) break;
			++current_temp_;
		}
		if ( current_temp_ > temperatures_.size() ) current_temp_ = temperatures_.size();
	} //if trusted
	runtime_assert( monte_carlo()->temperature() == temperatures_[ current_temp_ ] );
	return true;
}


void TemperingBase::init_from_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core;

	if ( !options_registered_ ) {
		utility_exit_with_message( "cannot call TemperingBase::init_from_options() unless TemperingBase::register_options() "
			"on application level before devel::init()" );
	}

	bool success( false );
	if ( option[ tempering::temp::file ].user() ) {
		success=initialize_from_file( option[ tempering::temp::file ]() );
		if ( !success ) tr.Info << "cannot read temperatures from file, will initialize from options... " << std::endl;
	}
	if ( !success ) {
		Real temp_high, temp_low;
		if ( option[ tempering::temp::range ].user() ) {
			temp_low=option[ tempering::temp::range ]().front();
			temp_high=option[ tempering::temp::range ]().back();
		} else {
			temp_low=option[ tempering::temp::low ]();
			temp_high=option[ tempering::temp::high ]();
		}
		Size const n_levels( option[ tempering::temp::levels ]() );
		generate_temp_range( temp_low, temp_high, n_levels );
	}
	stats_file_ = option[ tempering::stats::file ]();
	stats_silent_output_ =  option[ tempering::stats::silent ]();
	stats_line_output_ =  option[ tempering::stats::line_output ]() || stats_silent_output_;
	temperature_stride_ = option[ tempering::stride ]();
	instance_initialized_ = true;
}

void TemperingBase::generate_temp_range( Real temp_low, Real temp_high, Size n_levels, InterpolationType interpolation /*= linear*/ ) {
	temperatures_.clear();
	runtime_assert( n_levels >= 2 );
	tr.Info << "initializing temperatures from " << temp_low << " to " << temp_high << " with " << n_levels << " levels using "
		<< interpolation_type_enum_to_string(interpolation) << " interpolation." << std::endl;
	if ( interpolation == linear ) {
		Real const temp_step ( (temp_high-temp_low)/(n_levels-1) );
		for ( Size ct=0; ct<n_levels; ++ct ) {
			temperatures_.push_back( temp_low+ct*temp_step );
		}
	} else if ( interpolation == exponential ) {
		Real next_temp( temp_low );
		Real temp_factor( pow(temp_high/temp_low, 1/Real(n_levels-1)) );
		for ( Size ct=1; ct<n_levels; ++ct ) {
			temperatures_.push_back( next_temp );
			next_temp *= temp_factor;
		}
		temperatures_.push_back( temp_high );
	}
}

bool TemperingBase::initialize_from_file( std::string const& filename ) {
	clear();

	utility::io::izstream in( filename );
	if ( !in.good() ) {
		tr.Error << "cannot open file " << filename << std::endl;
		return false;
	}
	std::string line;
	getline( in, line );
	std::istringstream line_stream( line );
	std::string tag;
	Size n_levels( 0 );

	line_stream >> tag >> n_levels;

	bool line_format( false );
	if ( !line_stream.good() ) {
		tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
		tr.Error << "excpected TEMPERING or TEMPERING_TABLE" << std::endl;
		return false;
	}

	if ( tag == "TEMPERING" ) {
		line_format=true;
	} else if ( tag == "TEMPERING_TABLE" ) {
		line_format=false;
	} else {
		tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
		tr.Error << "excpected TEMPERING or TEMPERING_TABLE" << std::endl;
		return false;
	}

	Real temp;
	if ( line_format ) {
		for ( Size ct=1; ct <= n_levels; ++ct ) {
			line_stream >> temp;
			temperatures_.push_back( temp );
		}
		if ( !line_stream.good() ) {
			tr.Error << "format not recognized in temperature file: " << filename << " at line " << line << std::endl;
			tr.Error << "excpected TEMPERING N t1 t2 ... tN" << std::endl;
			return false;
		}
		return true;
	}

	// table format
	while ( getline( in, line ) ) {
		std::istringstream line_stream( line );
		Real temp;
		line_stream >> temp;
		if ( !line_stream.good() ) tr.Error << "format error in temperature file: " << filename << " at line " << line << std::endl;
		temperatures_.push_back( temp );
	}
	return true; //succesfully initialized
}

void TemperingBase::clear() {
	temperatures_.clear();
	instance_initialized_ = false;
}

void TemperingBase::write_to_file( std::string const& file_in, std::string const& output_name, utility::vector1< Real > const& ) {

	//either write all these things to a single file, or write to <jobname>.file_in
	std::string file;
	if ( stats_silent_output_ ) {
		file=file_in;
	} else {
		file=output_name+"."+file_in;
	}

	//open file
	utility::io::ozstream out( file, stats_line_output_ ? std::ios::app : std::ios::out );

	//write
	if ( stats_line_output_ ) { //line format
		out << "TEMPERING " << n_temp_levels();
		for ( Size i=1; i <= n_temp_levels(); ++i ) {
			out << " " << temperatures_[ i ];
		}
		out << " " << output_name << std::endl;
	} else { //table format
		out << std::setw( 10 );
		out << "TEMPERING_TABLE " << n_temp_levels() << " " << output_name << std::endl;
		for ( Size i=1; i <= n_temp_levels(); ++i ) {
			out << temperatures_[ i ] << std::endl;
		}
	}
}

void TemperingBase::set_temperatures( utility::vector1< Real > const& temps ) {
	temperatures_ = temps;
	set_current_temp( temps.size() );
}


void TemperingBase::set_current_temp( Size new_temp ) {
	current_temp_ = new_temp;
	Real real_temp = temperatures_[ current_temp_ ];
	if ( monte_carlo() ) {
		monte_carlo()->set_temperature( real_temp );
	}
	if ( job_ ) {
		job_->add_string_real_pair( "temperature", real_temp );
		job_->add_string_real_pair( "temp_level", new_temp );
	}
}


} //moves
} //protocols

