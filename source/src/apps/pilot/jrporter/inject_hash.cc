// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file topo_driver.cc
/// @brief  Inject proper hashes for cuts into silent structs.
/// @author Justin Porter

#include <basic/Tracer.hh>

#include <core/id/types.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>
#include <basic/datacache/WriteableCacheableData.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.hh>

#include <protocols/environment/AutoCutData.hh>

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>

#include <devel/init.hh>

#include <utility/string_util.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/excn/Exceptions.hh>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>


static THREAD_LOCAL basic::Tracer tr( "main" );


class InjectorMover : public protocols::moves::Mover {
  typedef protocols::moves::Mover Parent;
public:

  InjectorMover():
    Parent() {}

  void apply( core::pose::Pose& pose ){

    using namespace basic::datacache;
    using namespace core::pose::datacache;

    core::Size const hash = this->hash( pose.fold_tree() );
    basic::datacache::BasicDataCache& datacache = pose.data();

    if( !datacache.has( CacheableDataType::WRITEABLE_DATA ) ){
      datacache.set( CacheableDataType::WRITEABLE_DATA, new WriteableCacheableMap() );
    }
    WriteableCacheableMapOP datamap = dynamic_cast< WriteableCacheableMap* >( datacache.get_raw_ptr( CacheableDataType::WRITEABLE_DATA ) );

    std::set< core::Size > auto_cuts;
    for( int i = 1; i <= pose.fold_tree().num_cutpoint(); ++i ){
      auto_cuts.insert( (Size) pose.fold_tree().cutpoint(i) );
    }

    protocols::environment::AutoCutDataOP data = new protocols::environment::AutoCutData( hash, auto_cuts );

    datamap->insert( data );
  }

  virtual std::string get_name() const {
    return "InjectorMover";
  }

  virtual protocols::moves::MoverOP	fresh_instance() const{
    return new InjectorMover();
  }

  virtual protocols::moves::MoverOP clone() const {
    return new InjectorMover( *this );
  }


private:

  core::Size hash( core::kinematics::FoldTree const& ft ){
    utility::vector1< std::string > jump_points;

    for( Size i = 1; i <= ft.nres(); ++i ){
      for( Size j = i + 1; j <= ft.nres(); ++j ){
        if( ft.jump_exists( i, j ) ){
          jump_points.push_back( utility::to_string( i ) );
          jump_points.push_back( utility::to_string( j ) );
        }
      }
    }

    core::Size hashvalue = boost::hash_value( utility::join( jump_points, "," ) );

    tr.Info << "Hash calculated: " << hashvalue << std::endl;

    return hashvalue;
  }

};

typedef utility::pointer::owning_ptr< InjectorMover > InjectorMoverOP;

int main( int argc, char** argv ){

  devel::init( argc, argv );

  InjectorMoverOP mover = new InjectorMover();

  protocols::jd2::JobDistributor::get_instance()->go( mover );


}
