// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/FragmentJumpCM.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/movers/FragmentJumpCM.hh>

// Package headers
#include <core/environment/DofPassport.hh>
#include <core/environment/LocalPosition.hh>
#include <core/environment/LocalPosition.hh>

#include <protocols/environment/claims/JumpClaim.hh>
#include <protocols/environment/claims/TorsionClaim.hh>

#include <protocols/environment/movers/JumpSampleData.hh>
#include <protocols/environment/movers/FragmentJumpCMCreator.hh>

// Project headers
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>

#include <core/kinematics/MoveMap.hh>

#include <core/fragment/FragSet.hh>
#include <core/fragment/SecondaryStructure.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/abinitio/PairingStatistics.hh>
#include <protocols/abinitio/TemplateJumpSetup.hh>
#include <protocols/abinitio/Templates.hh>

#include <protocols/jumping/JumpSetup.hh>
#include <protocols/jumping/JumpSample.hh>
#include <protocols/jumping/RandomSheetBuilder.hh>

#include <core/scoring/dssp/PairingsList.hh>

//Utility Headers
#include <utility/io/izstream.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

#include <boost/functional/hash.hpp>
#include <boost/foreach.hpp>

// tracer
#include <basic/Tracer.hh>

// C++ Headers

// ObjexxFCL Headers

static basic::Tracer tr("protocols.environment.movers.FragmentJumpCM", basic::t_info);

namespace protocols {
namespace environment {

using namespace core::environment;

// creator
std::string
FragmentJumpCMCreator::keyname() const {
  return FragmentJumpCMCreator::mover_name();
}

protocols::moves::MoverOP
FragmentJumpCMCreator::create_mover() const {
  return new FragmentJumpCM;
}

std::string
FragmentJumpCMCreator::mover_name() {
  return "FragmentJumpCM";
}

FragmentJumpCM::FragmentJumpCM():
  Parent(),
  moverkey_( "" )
{}

FragmentJumpCM::FragmentJumpCM( std::string const& topol_filename,
                                std::string const& label,
                                std::string const& moverkey ) {
  set_topology( topol_filename );
  Parent( mover(), label );
}

void FragmentJumpCM::parse_my_tag( utility::tag::TagCOP const tag,
                                   basic::datacache::DataMap&,
                                   protocols::filters::Filters_map const&,
                                   protocols::moves::Movers_map const&,
                                   core::pose::Pose const& ) {
  if( tag->hasOption( "topol_file" ) ){
    std::string const& topol_filename = tag->getOption< std::string >( "topol_file", "" );
    set_topology( topol_filename );
  } else if ( tag->hasOption( "ss_info" ) &&
              tag->hasOption( "pairing_file" ) &&
              tag->hasOption( "n_sheets" ) ) {
    std::string const& ss_info_file = tag->getOption< std::string >( "ss_info", "" );
    std::string const& pairing_file = tag->getOption< std::string >( "pairing_file", "" );
    core::Size const& n_sheets = tag->getOption< core::Size >( "n_sheets", 0 );
    bool bRandomSheets = tag->getOption< bool >( "random_sheets", true );

    set_topology( ss_info_file, pairing_file, n_sheets, bRandomSheets );
  } else {
    tr.Error << "[ERROR] Parsing problem in FragmentJumpCM: "
             << " option set not compatible. Valid sets are 'topol_file' "
             << " or 'ss_info', 'pairing_file', and 'n_sheets'." << std::endl;
    throw utility::excn::EXCN_RosettaScriptsOption( "FragmentJumpCM incompatible options." );
  }

  set_label( tag->getOption< std::string >( "label", "BASE" ) );
  set_moverkey( tag->getOption< std::string >( "name" , "" ) );

  tr.Debug << "Initialzed " << get_name() << " with label " << label()
           << " and moverkey " << moverkey() << std::endl;

}

claims::EnvClaims FragmentJumpCM::yield_claims( core::pose::Pose& in_pose ){
  using namespace core::pose::datacache;
  using namespace basic::datacache;

  BasicDataCache& datacache = in_pose.data();
  if( datacache.has( CacheableDataType::WRITEABLE_DATA ) ){
    WriteableCacheableMap const& map = datacache.get< WriteableCacheableMap >( CacheableDataType::WRITEABLE_DATA );
    if( map.find( "JumpSampleData") != map.end() ){
      std::set< WriteableCacheableDataOP > const& data_set = map.find( "JumpSampleData" )->second;

      BOOST_FOREACH( WriteableCacheableDataOP data_ptr, data_set ){
        JumpSampleDataOP jumpdata_ptr = dynamic_cast< JumpSampleData* >( data_ptr.get() );
        assert( jumpdata_ptr );

        if( jumpdata_ptr->moverkey() == moverkey() ){
          tr.Debug << "Found JumpSampleOP in cache with hash " << moverkey() << std::endl;
          return build_claims( jumpdata_ptr->jump_sample() );
        } else if( tr.Debug.visible() ) {
          tr.Debug << "Rejected JumpSampleDataOP with hash " << jumpdata_ptr->moverkey()
                   << ", which conflicts with present hash " << moverkey() << std::endl;
        }
      }
    }
  }
  //If we don't hit the return statement in the above block of code, we have to calculate a new jump set.
  tr.Debug << get_name() << " calculating new jumps based on input topology file" << std::endl;

  if( !mover() ){
    throw utility::excn::EXCN_BadInput( "A topology file is required for FragmentJumpCM claiming of a new file." );
  }

  // called for each to fold tree get random jumps from top file.
  jumping::JumpSample jump_sample = setup_fragments();

  // cache the result as a WriteableCacheable in the pose DataCache for later retrieval.
  // this is important for the restarting feature in Abinitio.
  JumpSampleDataOP data = new JumpSampleData( moverkey(), jump_sample );
  if( !datacache.has( CacheableDataType::WRITEABLE_DATA ) ){
    datacache.set( CacheableDataType::WRITEABLE_DATA, new WriteableCacheableMap() );
  }
  datacache.get< WriteableCacheableMap >( CacheableDataType::WRITEABLE_DATA )[ data->datatype() ].insert( data );

  return build_claims( jump_sample );
}

claims::EnvClaims FragmentJumpCM::build_claims( jumping::JumpSample const& jump_sample ) {
  claims::EnvClaims claim_list;

  for( int i = 1; i <= (int) jump_sample.size(); i++ ){
    Size const up( jump_sample.jumps()( 1, i ) );
    Size const dn( jump_sample.jumps()( 2, i ) );

    std::string jump_name = get_name() + "Jump" + utility::to_string( i );

    claims::JumpClaimOP jclaim = new claims::JumpClaim( this,
                                                        jump_name,
                                                        LocalPosition( label(), up ),
                                                        LocalPosition( label(), dn ) );
    jclaim->ctrl_strength( claims::EXCLUSIVE );
    jclaim->init_strength( claims::MUST_INITIALIZE );

    jclaim->set_atoms( jump_sample.jump_atoms()(1, i ), jump_sample.jump_atoms()(1, i ) );
    jclaim->physical( false ); //jumps are not physical, and should be scored as chainbreaks

    claim_list.push_back( jclaim );

    // Jump Fragments make the mover into a bit of an access primadonna because jump fragments
    // include torsions from the residues at takeoff and landing
    claims::TorsionClaimOP tclaim_up = new claims::TorsionClaim( this, LocalPosition( label(), up ) );
    claims::TorsionClaimOP tclaim_dn = new claims::TorsionClaim( this, LocalPosition( label(), dn ) );

    tclaim_up->ctrl_strength( claims::MUST_CONTROL );
    tclaim_dn->ctrl_strength( claims::MUST_CONTROL );

    tclaim_up->init_strength( claims::MUST_INITIALIZE );
    tclaim_dn->init_strength( claims::MUST_INITIALIZE );

    claim_list.push_back( tclaim_up );
    claim_list.push_back( tclaim_dn );
  }

  return claim_list;
}

std::string FragmentJumpCM::get_name() const {
  return "FragmentJumpCM("+mover()->get_name()+")";
}

void FragmentJumpCM::set_topology( std::string const& ss_info_file,
                                   std::string const& pairing_file,
                                   core::Size const& n_sheets,
                                   bool bRandomSheets ) {
  if( mover() ){
    tr.Warning << "[WARNING] internal ClassicFragmentMover overwritten during FragmentJumpCM::process_topology_file call." << std::endl;
  }

  jumping::SheetBuilder::SheetTopology sheets;
  sheets.push_back( n_sheets );

  core::fragment::SecondaryStructureOP ss_def = new core::fragment::SecondaryStructure;
  ss_def->read_psipred_ss2( ss_info_file );

  core::scoring::dssp::PairingList pairlist;
  core::scoring::dssp::read_pairing_list( pairing_file, pairlist );

  if( bRandomSheets ){
    jump_def_ = new jumping::RandomSheetBuilder( ss_def, pairlist, sheets );
  } else {
    jump_def_ = new jumping::SheetBuilder( ss_def, pairlist, sheets );
  }

  // Fail faster with better information on bad topology files.
  try{
    calculate_jump_sample();
  } catch( utility::excn::EXCN_BadInput ){
    std::stringstream ss;
    ss << "Was not able to construct a valid jump sample in 10 attempts using ss_info file "
       << ss_info_file << ", pairing file " << pairing_file << " and " << n_sheets
       << "( sheets are random = " << bRandomSheets << ".";

    tr.Error << ss << std::endl;

    throw utility::excn::EXCN_BadInput( ss.str() );
  }

  setup_fragments();
}

void FragmentJumpCM::set_topology( std::string const& topol_filename ){

  if( mover() ){
    tr.Warning << "[WARNING] internal ClassicFragmentMover overwritten during FragmentJumpCM::process_topology_file call." << std::endl;
  }

  utility::io::izstream is( topol_filename );
  if ( !is.good() )
    throw utility::excn::EXCN_FileNotFound(" Topology file '" + topol_filename +"' not found." );

  abinitio::PairingStatisticsOP ps = new abinitio::PairingStatistics;
  is >> *ps;
  tr.Info << *ps << std::endl;
  core::fragment::SecondaryStructureOP ss_def = new core::fragment::SecondaryStructure;
  ss_def->extend( 10000 ); //Set number of residues to unreasonably large.
  core::scoring::dssp::PairingList helix_pairings; // helix pairings not used, required by BaseJumpSetup.

  jump_def_ = new abinitio::TemplateJumpSetup( NULL, ss_def, ps, helix_pairings );

  // Fail faster with better information on bad topology files.
  try{
    calculate_jump_sample();
  } catch( utility::excn::EXCN_BadInput ){
    throw utility::excn::EXCN_BadInput( "Was not able to construct a valid jump sample in 10 attempts using topology file " + topol_filename + "." );
  }

  setup_fragments();
}

jumping::JumpSample FragmentJumpCM::setup_fragments() {

  jumping::JumpSample jump_sample = calculate_jump_sample();

  core::kinematics::MoveMapOP dummy_mm = new core::kinematics::MoveMap;
  dummy_mm->set_bb( true );
  dummy_mm->set_jump( true );

  core::fragment::FragSetOP jump_frags = jump_def_->generate_jump_frags( jump_sample, *dummy_mm );

  if( !mover() ){
    simple_moves::ClassicFragmentMoverOP mover =
    new simple_moves::ClassicFragmentMover( jump_frags, dummy_mm );
    mover->set_check_ss( false ); //for some reason it's all 'L' notated, which causes a rejection of fragments...
    mover->enable_end_bias_check( false );
    set_mover( mover );
  } else {
    mover()->set_fragments( jump_frags );
  }

  return jump_sample;
}

jumping::JumpSample FragmentJumpCM::calculate_jump_sample() const {
  jumping::JumpSample jump_sample;

  core::Size attempts = 10;
  while( !jump_sample.is_valid() ){
    attempts -= 1;
    if( attempts <= 0 )
      throw utility::excn::EXCN_BadInput( "Was not able to construct a valid jump sample in 10 attempts." );

    jump_sample = jump_def_->create_jump_sample();
  }

  return jump_sample;
}

moves::MoverOP FragmentJumpCM::fresh_instance() const {
  return new FragmentJumpCM();
}

moves::MoverOP FragmentJumpCM::clone() const{
  return new FragmentJumpCM( *this );
}

} // environment
} // protocols
