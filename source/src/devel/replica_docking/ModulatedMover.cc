// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief rb_mover container, for n-replica, each modulatedMover holds n rb_mover with different
/// step size, then apply rb_mover according to the current temp
/// @author Zhe Zhang

#include <devel/replica_docking/ModulatedMover.hh>
#include <devel/replica_docking/ModulatedMoverCreator.hh>
#include <devel/replica_docking/TempInterpolatorFactory.hh>
#include <protocols/docking/RigidBodyInfo.hh>
#include <protocols/moves/MoverFactory.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <utility/tag/Tag.hh>
#include <core/pose/Pose.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>


static basic::Tracer tr( "devel.replica_docking.ModulatedMover" );

namespace devel {
namespace replica_docking {

std::string
ModulatedMoverCreator::keyname() const {
  return ModulatedMoverCreator::mover_name();
}

protocols::moves::MoverOP
ModulatedMoverCreator::create_mover() const {
  return new ModulatedMover;
}

std::string
ModulatedMoverCreator::mover_name() {
  return "ModulatedMover";
}

ModulatedMover::ModulatedMover() :
  tempering_( NULL )
{}

ModulatedMover::ModulatedMover( ModulatedMover const& other ) : ThermodynamicMover( other ) {
  ///copy value of every private variables
  tempering_ = other.tempering_;
  movers_ = other.movers_;
  interpolators_ = other.interpolators_;
}

ModulatedMover::~ModulatedMover() {}

std::string
ModulatedMover::get_name() const
{
  return "ModulatedMover";
}

protocols::moves::MoverOP
ModulatedMover::clone() const
{
  return new ModulatedMover(*this);
}

protocols::moves::MoverOP
ModulatedMover::fresh_instance() const
{
  return new ModulatedMover;
}

void
ModulatedMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data_map,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & pose
) {
  protocols::moves::MoverOP mover = protocols::rosetta_scripts::parse_mover( tag->getOption< std::string >( "tempering", "null" ), movers );
  tempering_ = dynamic_cast< protocols::canonical_sampling::HamiltonianExchange* >( mover() );
  if ( !tempering_ ) {
    throw utility::excn::EXCN_RosettaScriptsOption( "ModulatedMover requires an tempering argument" );
  }
  n_temp_levels_ = tempering_->n_temp_levels();
  if ( tag->hasOption( "type" ) ) {
    // I don't know yet how to use this type option, I mean even i got the mover type, so what?
    mover_name_ = tag->getOption< std::string >( "type" );
		set_type( mover_name_ ); //for proper reporting in trial_counters
  }

  utility::vector1< std::string > interpolated_options;
  utility::vector0< utility::tag::TagCOP > const subtags( tag->getTags() );


  for ( utility::vector0< utility::tag::TagCOP >::const_iterator subtag_it = subtags.begin(); subtag_it != subtags.end(); ++subtag_it ) {
    utility::tag::TagCOP const subtag = *subtag_it;
    tr.Debug << "subtag->getName() " << subtag->getName() << std::endl;
    if ( subtag->getName() == "Interp" && subtag->getOption< std::string >("key")!="weight" ) {
//       //      interpolator_[ ]=InterpolatorFactory->new_interpolator( subtag );
//       if ( tag->hasOption( "const" ) ) {
// 	//	interpolators_[ subtag->getOption< std::string >( "key" ) ] = new TempFixValue( subtag->getOption<core::Real>("const") );
//       } else {
// 	interpolators_[ subtag->getOption< std::string >( "key" ) ] = new TempInterpolator( subtag, tempering_ );
//       }
      interpolators_[ subtag->getOption< std::string >( "key" ) ] = TempInterpolatorFactory::get_instance()->new_tempInterpolator( subtag, tempering_->n_temp_levels() );
    } else {
      tr.Warning << "Cannot interpret subtag "<<subtag->getName() << std::endl;
    }
  }
  std::string our_name = tag->getOption< std::string >("name");
  //now crete mover list (maybe put this in extra method? )

  for ( core::Size temp_level = 1; temp_level <= n_temp_levels_; ++temp_level ) {
    //have a method that creates a tag for a certain temperature level
    utility::tag::TagCOP mover_tag = generate_mover_tag( temp_level, our_name );
    using namespace protocols::moves;
    using namespace protocols::canonical_sampling;
    MoverOP new_mover( MoverFactory::get_instance()->newMover( mover_tag, data_map, filters, movers, pose ) );
    if ( !new_mover ) {
      throw utility::excn::EXCN_RosettaScriptsOption( "you suck" );
    }
    ThermodynamicMoverOP th_mover( dynamic_cast<ThermodynamicMover*>( new_mover.get() ) );
    if ( !th_mover) {
      throw utility::excn::EXCN_RosettaScriptsOption( "you still suck: Class "+mover_name_+"is not a ThermodynamicMover" );
    }
    movers_.push_back( th_mover ); //
  } // end of for ( temp_level )
}

utility::tag::TagCOP
ModulatedMover::generate_mover_tag( core::Size temp_level, std::string const& prefix ) const {
  using namespace utility::tag;
  TagOP tag = new Tag;
  tag->setName(mover_name_);
  tag->setOption< std::string >("name",prefix+"_"+ ObjexxFCL::string_of(temp_level));
  for ( Interpolators::const_iterator it=interpolators_.begin(); it!=interpolators_.end(); ++it ) {
    //	string key=it->first
    //      InterpolatorOP interpolator_ptr=it->second
    tag->setOption< core::Real >( it->first, it->second->get_value(temp_level) );
    tr.Debug << "modulated " << it->first << " is " << tag->getOption< core::Real >( it->first ) << std::endl;
  }
  return tag;
}

void ModulatedMover::apply( core::pose::Pose& pose ) {
  core::Size temp_level( tempering_->temperature_level() );
  tr.Debug << "in ModulatedMover::apply current temp_level " << temp_level <<" n_temp_leves " << n_temp_levels_<< std::endl;
  if ( temp_level > n_temp_levels_ ) {
    utility_exit_with_message( "ERROR: temp_level should never be larger than n_temp_levels, something is wrong ");
  }
  movers_[ temp_level ]->apply( pose );
}

void ModulatedMover::initialize_simulation(
  core::pose::Pose & pose,
  protocols::canonical_sampling::MetropolisHastingsMover const & mhm,
  core::Size cycle
) {
  tr.Debug << "initialize_simulation, call Modulated movers' initialize_simulation"<< std::endl;
  for ( MoverOPs::const_iterator mover_it=movers_.begin(); mover_it !=movers_.end(); ++mover_it ) {
    (*mover_it)->initialize_simulation( pose, mhm, cycle );
  }
  tr.Debug << "Modulated movers' initialize_simulation done" << std::endl;
  //  Parent::initialize_simulation( pose, mhm, cycle );
//   for ( core::Size n_mover=1; n_mover<=movers_.size(); ++n_mover ) {
//     for ( core::Size i=1; i<=pose.num_jump(); ++i ) {
//       movers_[n_mover]->add_jump( i );
//       tr.Debug << "jump added into movers_[ " << n_mover << " ]" << std::endl;
//     }
//   }
}


void ModulatedMover::finalize_simulation(
  core::pose::Pose & pose,
  protocols::canonical_sampling::MetropolisHastingsMover const & mhm
) {
  // core::Size temp_level( tempering_->temperature_level() );
  // tr.Debug << "current temp_level " << temp_level << std::endl;
  // tr.Debug << "finalize_simulation of Modulated mover " << movers_[ temp_level ]->get_name() << std::endl;
  // movers_[ temp_level ]->finalize_simulation( pose, mhm );
  tr.Debug <<"finalize_simulation, call Modulated movers' finalize_simulation"<< std::endl;
  for ( MoverOPs::const_iterator mover_it=movers_.begin(); mover_it !=movers_.end(); ++mover_it ) {
    (*mover_it)->finalize_simulation( pose, mhm );
  }
  tr.Debug << "Modulated movers' finalize_simulation done" << std::endl;
}

// void ModulatedMover::observe_after_metropolis(
//   protocols::canonical_sampling::MetropolisHastingsMover const & mhm
// ) {
//   core::Size temp_level( tempering_->temperature_level() );
//   tr.Debug << "current temp_level " << temp_level << " observe_after_metropolis... " ;
//   movers_[ temp_level ]->observe_after_metropolis( mhm );
//   tr.Debug << "done!" << std::endl;
// }

} //devel
} //replica_docking
