// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#include <devel/replica_docking/TempInterpolatorFactory.hh>
#include <devel/replica_docking/TempInterpolator.hh>

#include <protocols/canonical_sampling/TemperatureController.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption
#include <utility/tag/Tag.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

static basic::Tracer tr("devel.replica_docking.TempInterpolatorFactory");

namespace devel {
namespace replica_docking {


static basic::Tracer TR( "devel.replica_docking.TempInterpolatorFactory" );

TempInterpolatorFactory::TempInterpolatorFactory(){}

TempInterpolatorFactory::~TempInterpolatorFactory(){}

TempInterpolatorFactory * TempInterpolatorFactory::instance_( 0 );

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex TempInterpolatorFactory::singleton_mutex_;

std::mutex & TempInterpolatorFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
TempInterpolatorFactory * TempInterpolatorFactory::get_instance()
{
	boost::function< TempInterpolatorFactory * () > creator = boost::bind( &TempInterpolatorFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

TempInterpolatorFactory *
TempInterpolatorFactory::create_singleton_instance()
{
	return new TempInterpolatorFactory;
}


///@brief return new TempInterpolator by key lookup in mover_prototype_map_ (new TempInterpolator parses Tag if provided)
TempInterpolatorBaseOP
TempInterpolatorFactory::new_tempInterpolator( utility::tag::TagCOP const tag, core::Size n_levels )
{
	if ( !tag->hasOption( "curve" ) ){
		throw utility::excn::EXCN_RosettaScriptsOption("Error: interpolation curve type required !" );
	}
	std::string curve = tag->getOption< std::string >( "curve" );
	if ( curve == "const" ) {
		return new TempFixValue( tag->getOption< core::Real >( "value" ) );
	}
	else {
		if ( tag->hasOption( "start" ) && tag->hasOption( "end" ) ) {
			core::Real start = tag->getOption< core::Real >( "start" );
			core::Real end = tag->getOption< core::Real >( "end" );
			return new TempInterpolator( n_levels, start, end, curve );
		}
		else {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: start and end value must be given for linear and expotential interpolation !" );
		}
	}
}


} //namespace devel
} //namespace replica_docking
