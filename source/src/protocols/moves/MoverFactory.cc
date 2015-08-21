// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/moves/MoverFactory.cc
/// @brief
/// @author ashworth

#include <protocols/moves/MoverFactory.hh>
#include <protocols/moves/Mover.hh>

// required for passing to Mover::parse_my_tag

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


namespace protocols {
namespace moves {


static thread_local basic::Tracer TR( "protocols.moves.MoverFactory" );

#if defined MULTI_THREADED && defined CXX11
std::atomic< MoverFactory * > MoverFactory::instance_( 0 );
#else
MoverFactory * MoverFactory::instance_( 0 );
#endif

#ifdef MULTI_THREADED
#ifdef CXX11

std::mutex MoverFactory::singleton_mutex_;

std::mutex & MoverFactory::singleton_mutex() { return singleton_mutex_; }

#endif
#endif

/// @brief static function to get the instance of ( pointer to) this singleton class
MoverFactory * MoverFactory::get_instance()
{
	boost::function< MoverFactory * () > creator = boost::bind( &MoverFactory::create_singleton_instance );
	utility::thread::safely_create_singleton( creator, instance_ );
	return instance_;
}

MoverFactory *
MoverFactory::create_singleton_instance()
{
	return new MoverFactory;
}

MoverFactory::MoverFactory()
{
	forbidden_names_.clear();
	//the following 5 names are used in the original
	//LoopMoverFactory
	forbidden_names_.insert("quick_ccd");
	forbidden_names_.insert("sdwindow");
	forbidden_names_.insert("quick_ccd_moves");
	forbidden_names_.insert("perturb_ccd");
	forbidden_names_.insert("perturb_kic");
	//original loop mover names over
}

MoverFactory::~MoverFactory(){}

/// @brief add a Mover prototype, using its default type name as the map key
void
MoverFactory::factory_register( MoverCreatorOP creator )
{
	runtime_assert( creator != 0 );
	std::string const mover_type( creator->keyname() );
	if ( mover_type == "UNDEFINED NAME" ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Can't map derived Mover with undefined type name.");
	}
	if ( forbidden_names_.find( mover_type ) != forbidden_names_.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption("Name "+mover_type+" is not an allowed mover name, probably because it has historical meaning.");
	}
	if ( mover_creator_map_.find( mover_type ) != mover_creator_map_.end() ) {
		throw utility::excn::EXCN_RosettaScriptsOption("MoverFactory::factory_register already has a mover creator with name \"" + mover_type + "\".  Conflicting Mover names" );
	}
	mover_creator_map_[ mover_type ] = creator;
}


/// @brief return new Mover by key lookup in mover_prototype_map_ (new Mover parses Tag if provided)
MoverOP
MoverFactory::newMover( std::string const & mover_type )
{
	MoverMap::const_iterator iter( mover_creator_map_.find( mover_type ) );
	if ( iter != mover_creator_map_.end() ) {
		if ( ! iter->second ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "Error: MoverCreatorOP prototype for " + mover_type + " is NULL!" );
		}
		// use of cloning method would be faithful to pre-initialized prototypes
		//return iter->second->clone();
		// fresh_instance prevents propagation of pre-initialized prototypes, which may be safer(?)
		return iter->second->create_mover();
	} else {
		TR<<"Available movers: ";
		for ( MoverMap::const_iterator mover_it = mover_creator_map_.begin(); mover_it != mover_creator_map_.end(); ++mover_it ) {
			TR<<mover_it->first<<", ";
		}
		TR<<std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( mover_type + " is not known to the MoverFactory. Was it registered via a MoverRegistrator in one of the init.cc files (devel/init.cc or protocols/init.cc)?" );
		return NULL;
	}
}

/// @brief return new Mover by Tag parsing
MoverOP
MoverFactory::newMover(
	TagCOP const tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	moves::Movers_map const & movers,
	Pose const & pose )
{
	MoverOP mover( newMover( tag->getName() ) );
	runtime_assert( mover != 0 );
	mover->parse_my_tag( tag, data, filters, movers, pose );
	return mover;
}

} //namespace moves
} //namespace protocols
