// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/abinitio/abscript/FragmentCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/abinitio/abscript/FragmentCM.hh>
#include <protocols/abinitio/abscript/FragmentCMCreator.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <protocols/environment/DofUnlock.hh>

#include <protocols/environment/claims/TorsionClaim.hh>
#include <protocols/environment/claims/CutBiasClaim.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/FragmentIO.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <core/pack/task/residue_selector/ResidueSelector.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/simple_moves/SmoothFragmentMover.hh>
#include <protocols/simple_moves/GunnCost.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

// utility/basic headers
#include <utility/tag/Tag.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/abinitio.OptionKeys.gen.hh>

#include <basic/datacache/DataMap.hh>

#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.environment.movers.FragmentCM", basic::t_info );

namespace protocols {
namespace abinitio {
namespace abscript {

using namespace core::environment;
using namespace environment;


// creator
std::string
FragmentCMCreator::keyname() const {
  return FragmentCMCreator::mover_name();
}

protocols::moves::MoverOP
FragmentCMCreator::create_mover() const {
  return protocols::moves::MoverOP( new FragmentCM );
}

std::string
FragmentCMCreator::mover_name() {
  return "FragmentCM";
}


FragmentCM::FragmentCM():
  ClaimingMover(),
  mover_( /* NULL */ ),
  selector_( /* NULL */ ),
  bInitialize_( true ),
  bYieldCutBias_( false )
{}

FragmentCM::FragmentCM( simple_moves::FragmentMoverOP mover,
                        core::pack::task::residue_selector::ResidueSelectorCOP selector ):
  ClaimingMover(),
  selector_( selector ),
  bInitialize_( true ),
  bYieldCutBias_( false )
{
  set_mover( mover );
}

FragmentCM::~FragmentCM() {}

void FragmentCM::set_selector( core::pack::task::residue_selector::ResidueSelectorCOP selector ) {
  if( Parent::state_check( __FUNCTION__, ( selector.get() == selector_.get() ) ) ){
    selector_ = selector;
  }
}

void FragmentCM::set_mover( simple_moves::FragmentMoverOP mover ){
  if( Parent::state_check( __FUNCTION__, ( mover.get() == mover_.get() ) ) ){
    mover_ = mover;
    type( get_name() );
  }
}

core::Size determine_frag_size( std::string const& file ){
  using namespace core::fragment;
  using namespace basic::options;

  FragmentIO frag_io( option[ OptionKeys::abinitio::number_3mer_frags ](), 1,
                      option[ OptionKeys::frags::annotate ]() );

  FragSetOP fragset = frag_io.read_data( file );

  // Assumes constant length fragments. I don't think anyone uses these (or even knows its possible)
  // so this seems safe.
  return fragset->max_frag_length();
}

void FragmentCM::parse_my_tag( utility::tag::TagCOP tag,
                                  basic::datacache::DataMap& datamap,
                                  protocols::filters::Filters_map const&,
                                  protocols::moves::Movers_map const&,
                                  core::pose::Pose const& ) {
  using namespace core::fragment;
  using namespace basic::options;

  std::string frag_file_tag = "fragments";

  core::Size const frag_length = determine_frag_size( tag->getOption< std::string >( frag_file_tag ) );

  // The number of frags to be randomly selected out of the total fragments. I think.
  core::Size n_frags;
  if( frag_length == 3 ){
    n_frags = tag->getOption< core::Size >( "nfrags", option[ OptionKeys::abinitio::number_3mer_frags ]() );
  } else if ( frag_length == 9 ){
    n_frags = tag->getOption< core::Size >( "nfrags", option[ OptionKeys::abinitio::number_9mer_frags ]() );
  } else {
    n_frags = tag->getOption< core::Size >( "nfrags" );
  }

  FragmentIO frag_io( n_frags, 1, option[ OptionKeys::frags::annotate ]() );

  std::string const& frag_type = tag->getOption< std::string >( "frag_type", "classic" );
  if( frag_type == "classic" ){
    set_mover( simple_moves::FragmentMoverOP( new simple_moves::ClassicFragmentMover( frag_io.read_data( tag->getOption< std::string >( frag_file_tag ) ) ) ) );
  } else if( frag_type == "smooth" ){
    set_mover( simple_moves::FragmentMoverOP( new simple_moves::SmoothFragmentMover( frag_io.read_data( tag->getOption< std::string >( frag_file_tag ) ), protocols::simple_moves::FragmentCostOP( new simple_moves::GunnCost() ) ) ) );
  } else {
    std::ostringstream ss;
    ss << "The fragment type " << frag_type << " is not valid. The options "
       << " are 'classic' and 'smooth'.";
    throw utility::excn::EXCN_RosettaScriptsOption( ss.str() );
  }

  initialize( tag->getOption< bool >( "initialize", true ) );

  set_selector( datamap.get_ptr< core::pack::task::residue_selector::ResidueSelector const >( "ResidueSelector", tag->getOption<std::string>( "selector" ) ) );
}

claims::EnvClaims FragmentCM::yield_claims( core::pose::Pose const& pose,
                                            basic::datacache::WriteableCacheableMapOP ){
  using namespace claims;
  claims::EnvClaims claim_list;

  if( yield_cut_bias() ){
    core::fragment::SecondaryStructureOP ss( new core::fragment::SecondaryStructure( *( mover()->fragments() ) ) );
    claim_list.push_back( protocols::environment::claims::EnvClaimOP( new environment::claims::CutBiasClaim(
		utility::pointer::static_pointer_cast< ClaimingMover > ( get_self_ptr() ),
		"BASE",
		*ss ) ) );
  }

  int shift = 0;
  if( selector() ){
    utility::vector1< bool > torsion_mask = selector()->apply( pose );
    shift = torsion_mask.index( true )-1;
    mover()->set_fragments( mover()->fragments()->clone_shifted( shift ) );
  }

  TorsionClaimOP new_claim( new TorsionClaim( utility::pointer::static_pointer_cast< ClaimingMover > ( get_self_ptr() ),
                                               "BASE",
                                               std::make_pair( mover()->fragments()->min_pos(),
                                                               mover()->fragments()->max_pos() ) ) );

  if( initialize() ){
    new_claim->strength( CAN_CONTROL, CAN_CONTROL );
  } else {
    new_claim->strength( CAN_CONTROL, DOES_NOT_CONTROL );
  }
  claim_list.push_back( new_claim );

  return claim_list;
}

void FragmentCM::initialize( Pose& pose ){
  DofUnlock activation( pose.conformation(), passport() );
  mover()->apply_at_all_positions( pose );
}

void FragmentCM::apply( Pose& pose ) {
  if( pose.conformation().is_protected() ) {
    DofUnlock activation( pose.conformation(), passport() );
    mover()->apply( pose );
  } else {
    mover()->apply( pose );
  }
}

void FragmentCM::initialize( bool setting ) {
  if( Parent::state_check( __FUNCTION__, bInitialize_ == setting ) ){
    bInitialize_ = setting;
  }
}

void FragmentCM::yield_cut_bias( bool setting ){
  if( Parent::state_check( __FUNCTION__, bYieldCutBias_ == setting ) ){
    bYieldCutBias_ = setting;
  }
}

void FragmentCM::passport_updated() {
  assert( mover() );
  if( has_passport() ){
    mover()->set_movemap( passport()->render_movemap() );
  } else {
    core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap() );
    mm->set_bb( true );
    mover()->set_movemap( mm );
  }
}

std::string FragmentCM::get_name() const {
  return "FragmentCM("+mover()->get_name()+")";
}

} // abscript
} // abinitio
} // protocols
