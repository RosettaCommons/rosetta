// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvClaimBroker.cc
/// @author Justin Porter

// Unit Headers
#include <protocols/environment/EnvClaimBroker.hh>

// Package headers
#include <core/chemical/VariantType.hh>

#include <core/environment/FoldTreeSketch.hh>

#include <core/conformation/util.hh>

#include <protocols/environment/claims/EnvClaim.hh>
#include <protocols/environment/claims/VirtResClaim.hh>
#include <protocols/environment/claims/JumpClaim.hh>

#include <protocols/environment/ProtectedConformation.hh>
#include <protocols/environment/ClaimingMover.hh>
#include <protocols/environment/Environment.hh>
#include <protocols/environment/AutoCutData.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/types.hh>

#include <utility/excn/Exceptions.hh>

#include <core/chemical/ResidueTypeSet.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Exceptions.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/tree/BondedAtom.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Atom.hh>

#include <core/pose/Pose.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/WriteableCacheableMap.hh>
#include <basic/datacache/WriteableCacheableData.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/PDBInfo.hh>

// utility headers
#include <utility/string_util.hh>
#include <boost/foreach.hpp>

#include <boost/foreach.hpp>
#include <boost/functional/hash.hpp>

// tracer
#include <basic/Tracer.hh>

// C++ Headers
#include <functional>
#include <algorithm>
#include <numeric>
#include <set>

// ObjexxFCL Headers

static thread_local basic::Tracer tr( "protocols.environment.EnvClaimBroker", basic::t_info );

namespace protocols {
namespace environment {

using namespace core;
using namespace protocols::environment::claims;

using core::environment::LocalPosition;
using core::environment::SequenceAnnotation;
using core::environment::SequenceAnnotationOP;
using core::environment::SequenceAnnotationCOP;
using core::environment::FoldTreeSketch;

using core::conformation::Conformation;

void update_pdb_info( core::pose::PDBInfoCOP input_pdb_info, core::pose::Pose& pose ){
  //Rebuild an appropriate PDBInfo object.
  if( pose.pdb_info() ){
    tr.Error << "Environment does not expect a PDBInfo object to be created during broking. Something has gone wrong!" << std::endl;
    utility_exit_with_message( "Problem in broking!" );
  } else if( input_pdb_info ) {
    //ASSUMPTION: all new residues are virtual
    core::Size const new_vrts = pose.total_residue() - input_pdb_info->nres();

    core::pose::PDBInfoOP new_info( new core::pose::PDBInfo( *input_pdb_info ) );

    for( Size i = 1; i <= new_vrts; ++i ){
      new_info->append_res( new_info->nres(), 3 );
    }

    tr.Debug << "Updating PDBInfo object to account for " << new_vrts << " (temporary) virtual residues in new pose of size "
             << pose.total_residue() << ". Old Size: " << input_pdb_info->nres() << "; New Size: "
             << new_info->nres() << std::endl;

    pose.pdb_info( new_info );
  } else {
    tr.Debug << "  PDBInfo processing being ignored as it is null in the input pose." << std::endl;
  }
}


void safe_set_conf( core::pose::Pose& pose, core::conformation::ConformationOP conf ){
  if( pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ){
    basic::datacache::WriteableCacheableMapOP data_map = utility::pointer::dynamic_pointer_cast< basic::datacache::WriteableCacheableMap > ( pose.data().get_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) );
    assert( data_map );
    pose.set_new_conformation( conf );
    pose.data().set( core::pose::datacache::CacheableDataType::WRITEABLE_DATA, data_map );
  } else {
    pose.set_new_conformation( conf );
  }
}

EnvClaimBroker::EnvClaimBroker( EnvironmentCAP env,
                                MoverPassMap const& movers_and_passes,
                                core::pose::Pose const& in_pose,
                                SequenceAnnotationOP ann ):
  movers_and_passes_( movers_and_passes ),
  ann_( ann ),
  env_( env )
{
  core::conformation::Conformation const & in_conf = in_pose.conformation();
  core::pose::Pose pose( in_pose );
  claims_ = collect_claims( movers_and_passes, pose );
  BOOST_FOREACH( EnvClaimOP claim, claims_ ){ claim->annotate( pose, ann ); }

  // Build a temporary conformation for which we manipulate the fold tree. After the fold tree is set,
  // we "seal it" as a ProtectedConformation for initializing DoFs.
  ConformationOP tmp_conf_op( new Conformation( in_conf ) );
  Conformation & tmp_conf = *tmp_conf_op;

  broker_fold_tree( tmp_conf, pose.data() );
  result_.ann = ann_;

  ProtectedConformationOP conf( new ProtectedConformation( env, tmp_conf ) );
  conf->attach_annotation( ann_ );

  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes ){
    pair.second->reference_conformation( conf );
  }

  // Bookkeeping for pose-associated info caches (datacache and pdb_info).
  safe_set_conf( pose, conf );
  update_pdb_info( in_pose.pdb_info(), pose );

  broker_dofs( pose );

  result_.pose = pose;

  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes ){
    pair.first->broking_finished( result() );
  }

  assert( utility::pointer::dynamic_pointer_cast< basic::datacache::WriteableCacheableMap >( pose.data().get_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ) );
}

EnvClaimBroker::~EnvClaimBroker() {}

EnvClaimBroker::BrokerResult const& EnvClaimBroker::result() const {
  return result_;
}

void EnvClaimBroker::broker_fold_tree( Conformation& conf,
                                       basic::datacache::BasicDataCache& datacache ){

  core::environment::FoldTreeSketch fts( conf.size() );
  result_.closer_ft = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree( conf.fold_tree() ) );

  //FTElements: virtual residues ---------------------------------------------------------------
  ResElemVect r_elems = collect_elements< ResidueElement >( fts );
  SizeToStringMap new_vrts;
  process_elements( r_elems, fts, new_vrts );
  result_.new_vrts = new_vrts;

  assert( fts.nres() == conf.size() + result().new_vrts.size() );

  //set up bogus attachments for new VRTs
  result_.closer_ft = core::kinematics::FoldTreeOP( new core::kinematics::FoldTree( conf.fold_tree() ) );
  for( Size i = 1; i <= result_.new_vrts.size(); ++i ){
    result_.closer_ft->insert_residue_by_jump( conf.size()+i, conf.size()+i-1 );
  }

  //FTElements: jumps --------------------------------------------------------------------------
  JumpElemVect j_elems = collect_elements< JumpElement >( fts );
  JumpDataMap new_jumps;
  process_elements( j_elems, fts, new_jumps );

  //FTElements cut -----------------------------------------------------------------------------
  CutElemVect c_elems = collect_elements< CutElement >( fts );
  process_elements( c_elems, fts );

  //FTElements CutBias -------------------------------------------------------------------------
  CutBiasElemVect cb_elems = collect_elements< CutBiasElement >( fts );
  utility::vector1< core::Real > bias( fts.nres(), 1.0 );
  process_elements( cb_elems, bias );

  //Render FoldTree ----------------------------------------------------------------------------
  tr.Info << "Claim collection complete. Constructing consensus FoldTree." << std::endl;
  core::kinematics::FoldTreeOP ft = render_fold_tree( fts, bias, datacache, conf );

  annotate_fold_tree( ft, new_jumps, ann_ );

  add_virtual_residues( conf, new_vrts, ann_ );

  tr.Debug << "EnvClaimBroker setting consensus fold tree " << *ft << std::endl;
  conf.fold_tree( *ft );

  //Swap out unphysical cut residues for chainbreak variants
  BOOST_FOREACH( Size cut , result().auto_cuts ){
    add_chainbreak_variants( cut, conf );
  }

  tr.Info << "Broking finished. Consensus fold tree: " << *ft << std::endl;
}

core::Size find_implied_cut( utility::vector1< core::Size > const& cycle,
                             core::conformation::Conformation const& conf ){
  utility::vector1< core::Size > cuts;

  core::Real const CA_CA_CUTOFF = 4.0;

  for( Size i = 2; i <= conf.size(); ++i ){
    core::Real const dist = conf.residue( i ).xyz( "CA" ).distance( conf.residue( i - 1 ).xyz( "CA" ) );
    if( dist > CA_CA_CUTOFF ){
      cuts.push_back( i );
    }
  }

  tr.Debug << "  Looking for implied cuts from spatial information indicating cuts at: " << cuts << std::endl;;
  for( Size i = 1; i <= cuts.size(); ++i ){
    core::Size const& cut = cuts[i];
    if( std::find( cycle.begin(), cycle.end(), cut ) != cycle.end() ){
      tr.Debug << "  Inheriting implied cut at " << cut << " to close ft cycle." << std::endl;
      return cut;
    }
  }

  return 0;
}

utility::vector1< core::Size > introduce_datamap_cuts( FoldTreeSketch const& fts,
                                                       basic::datacache::BasicDataCache& datacache ){

  using namespace basic::datacache;
  using namespace core::pose::datacache;

  // Calculate jump_points in prep for hash.
  utility::vector1< std::string > jump_points;
  for( Size i = 1; i <= fts.nres(); ++i ){
    for( Size j = i + 1; j <= fts.nres(); ++j ){
      if( fts.has_jump( i, j ) ){
        jump_points.push_back( utility::to_string( i ) );
        jump_points.push_back( utility::to_string( j ) );
      }
    }
  }

  //utility::join can't handle empty lists...
  core::Size const hash = ( jump_points.empty() ) ? boost::hash_value( "" ) :
  boost::hash_value( utility::join( jump_points, "," ) );

  if( !datacache.has( CacheableDataType::WRITEABLE_DATA ) ){
    datacache.set( CacheableDataType::WRITEABLE_DATA, DataCache_CacheableData::DataOP( new WriteableCacheableMap() ) );
  }

  WriteableCacheableMapOP datamap = utility::pointer::dynamic_pointer_cast< WriteableCacheableMap >( datacache.get_ptr( CacheableDataType::WRITEABLE_DATA ) );
  utility::vector1< core::Size > cuts;

  if( datamap && datamap->find( "AutoCutData" ) != datamap->end() ){
    BOOST_FOREACH( basic::datacache::WriteableCacheableDataOP data,
                  datamap->find( "AutoCutData" )->second ) {
      AutoCutDataOP auto_cut_data = utility::pointer::dynamic_pointer_cast< protocols::environment::AutoCutData > ( data );

      assert( auto_cut_data );

      if( auto_cut_data->hash() == hash ){
        tr.Debug << "  Found appropriate hash-valued AutoCutData. Repeating cut choices." << std::endl;

        BOOST_FOREACH( core::Size cut, auto_cut_data->cuts() ){
          cuts.push_back( cut );
        }

      } else {
        tr.Debug << "  Rejecting AutoCutData because its hash " << auto_cut_data->hash()
        << " does not match EnvClaimBroker hash " << hash << std::endl;
      }
    }
  } else {
    tr.Debug << "  No AutoCutData found. Cuts will not be generated from pose-cached values." << std::endl;
  }

  return cuts;
}

core::Size inherit_cuts( utility::vector1< core::Size > const& cycle,
                         core::kinematics::FoldTree const& input_ft ) {
  for( int i = 1; i <= input_ft.num_cutpoint(); ++i ){
    if( std::find( cycle.begin(), cycle.end(), input_ft.cutpoint( i ) ) != cycle.end() ){
      return input_ft.cutpoint( i );
    } else {
      tr.Trace << "  Inherited cut " << input_ft.cutpoint( i ) << " doesn't fit in cycle: "
      << cycle << std::endl;
    }
  }
  return 0;
}


core::kinematics::FoldTreeOP
EnvClaimBroker::render_fold_tree( FoldTreeSketch& fts,
                                  utility::vector1< core::Real > const& bias,
                                  basic::datacache::BasicDataCache& datacache,
                                  core::conformation::Conformation const& input_conf ) {



  utility::vector1< core::Size > const datamap_autocuts = introduce_datamap_cuts( fts, datacache );
  for( Size i = 1; i <= datamap_autocuts.size(); ++i ){
    core::Size const autocut = datamap_autocuts[ i ];
    fts.insert_cut( autocut );
    result_.auto_cuts.insert( autocut );
  }

  utility::vector1< Size > cycle = fts.cycle();
    
  // this is ugly, but seems to allow me to access the contents of env_
  // must some property of the new pointer system?
  EnvironmentCOP env( env_ );
    
  while( !cycle.empty() ){
    tr.Debug << "Trying to close cycle: " << cycle << std::endl;

    bool made_cut = false;
      
    if( env->inherit_cuts() ){
      core::Size const inherited_cut = inherit_cuts( cycle, input_conf.fold_tree() );
      if( inherited_cut ){
        fts.insert_cut( inherited_cut );
        made_cut = true;
        tr.Debug << "  Inheriting cut at " << inherited_cut << " to close ft cycle." << std::endl;
      }
    }

    if( !made_cut ){
      tr.Debug << "  Cut inheritance failed." << std::endl;
    }

    if( !made_cut && !datamap_autocuts.empty() && env->auto_cut() ){
      core::Size const cut = find_implied_cut( cycle, input_conf );
      if( cut ){
        fts.insert_cut( cut );
        result_.auto_cuts.insert( cut );
        made_cut = true;
        tr.Debug << "  Inheriting implied cut at " << cut << "." << std::endl;
      }
    }

    if( !made_cut ){
      tr.Debug << "  Automatic spatial inferred cut inheritance failed." << std::endl;
    }

    if( !made_cut && datamap_autocuts.empty() && env->auto_cut() ){
      utility::vector1< core::Real > masked_bias( bias.size(), 0 );

      // k must be strictly less than cycle.size() so that the automatic
      // cut doesn't get placed at the end of the cycle.
      core::Real bias_sum = 0.0;
      for( Size k = 1; k < cycle.size(); ++k ){
        if( fts.has_jump( cycle[k], cycle[k+1] ) ){
          // if a jump originating at residue k is part of the cycle, do not cut here because
          // a cut at this residue is actually *outside* of the loop, since cuts are placed after k.
          masked_bias[cycle[k]] = 0.0;
        } else {
          masked_bias[cycle[k]] = bias[cycle[k]];
          bias_sum += bias[cycle[k]];
        }
      }

      try {
        Size const cut = fts.insert_cut( masked_bias );
        tr.Debug << "  Inserted automatic cut at " << cut << std::endl;
        result_.auto_cuts.insert( cut );
        made_cut = true;
      } catch( ... ){
        tr.Error << "[ERROR] Cut insertion failed with biases: [";
        for( Size k = 1; k <= cycle.size(); ++k ){
          tr.Error << cycle[k] << ":" << bias[cycle[k]] << ", ";
        }
        tr.Error << std::endl;
        throw;
      }
    }

    if( !made_cut ) {
      std::ostringstream ss;
      ss << "Brokering failed (" << __FILE__ << ":" << __LINE__ << ") because the broker couldn't construct a fold tree. "
         << "It could not find a way to resolve the fold tree cycle: "  << cycle <<". Inherited cuts were "
         << ( env->inherit_cuts() ? "on" : "off" ) << " and automatic cuts were " << ( env->auto_cut() ? "on" : "off" )
         << ". At the point of failure were " << fts.num_jumps() << " jumps and "<< fts.num_cuts() << " cuts in the FoldTreeSketch.";
      if( env->inherit_cuts() ){
        ss << " Availiable inherited cuts are: " << input_conf.fold_tree().cutpoints() << ".";
      } if( env->auto_cut() ){
        ss << " Cut bias was " << bias << "." << std::endl;
      }
      throw utility::excn::EXCN_BadInput( ss.str() );
    }
    cycle = fts.cycle();
  }

  try {
    return fts.render();
  } catch ( ... ) {
    tr.Error << "[ERROR] A problem was encountered rendering fold tree"
             << "choices in EnvClaimBroker::render_fold_tree." << std::endl;
    throw;
  }
}

void EnvClaimBroker::annotate_fold_tree( core::kinematics::FoldTreeOP ft,
                                         JumpDataMap const& new_jumps,
                                         SequenceAnnotationOP ann ) {

  //Bind jump numbers to labels in annotations
  BOOST_FOREACH( JumpDataMap::value_type pair, new_jumps ){
    std::string const& label = pair.first;
    BrokeredJumpDataCOP jump_data = pair.second;

    Size jump_id = ft->jump_nr( jump_data->pos.first, jump_data->pos.second );

    if( jump_id != 0 ){ // The physical fold tree won't have some jumps, so we don't want to fail.
      if( ann ){ // Physical fold tree doesn't use the annotations.
        ann->add_jump_label( label, jump_id );
      }

      std::string const& a1 = jump_data->atoms.first;
      std::string const& a2 = jump_data->atoms.second;

      ft->set_jump_atoms( (int) jump_id, a1, a2,
                          jump_data->put_jump_stub_intra_residue );
    }
  }

  // like FoldTree::put_jump_stubs_intra_residue, except that it doesn't blindly
  // set all jumps without atom info.
  ft->reassign_atoms_for_intra_residue_stubs();
}

void EnvClaimBroker::add_virtual_residues( Conformation& conf,
    SizeToStringMap const& new_vrts,
    SequenceAnnotationOP ann )
{
  // Add new virtual residues into conformation.
  BOOST_FOREACH( SizeToStringMap::value_type pair, new_vrts ) {
  // Steal the residue type set of the first residue. Will obviously break if the conformation
  // has no residues. Is this a case I need to worry about?
  core::chemical::ResidueTypeSet const & rsd_set( conf.residue(1).residue_type_set() );
  core::conformation::ResidueOP rsd(
    core::conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );

  // where the jump goes doesn't matter since the current fold tree is about to be replaced by 'ft'.
  conf.append_residue_by_jump( *rsd, conf.size() );

  // This residue label should resolve to the VRT just added.
  assert( ann->resolve_seq( LocalPosition( pair.second, 1 ) ) == conf.size() );
  }
}

ControlStrength const& init_str_selector( std::pair< claims::DOFElement, ClaimingMoverOP > const& d ){
  return d.first.i_str;
}

ControlStrength const& ctrl_str_selector( std::pair< claims::DOFElement, ClaimingMoverOP > const& d ){
  return d.first.c_str;
}

void EnvClaimBroker::broker_dofs( core::pose::Pose& pose ){

  std::set< ClaimingMoverOP > initializers;

  //Broker Arbitrary DOFs ---------------------------------------------------------------------------
  DOFElemVect d_elems = collect_elements< DOFElement >( pose );
  setup_passports( d_elems, init_str_selector );

  // Now that passports are properly configured, allow each mover to initialize.
  // As soon as each mover initializes, access is revoked, so it can be built again
  // for the sampling phase.

  tr.Debug << "Beginning initialization:" << std::endl;
  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes_ ){
    if( pair.second->begin() != pair.second->end() ){
      tr.Debug << "  Invoking initialize: " << pair.first->get_name() << std::endl;
      pair.first->passport_updated();
      pair.first->initialize( pose );
      pair.second->revoke_all_access();
    }
  }

  setup_passports( d_elems, ctrl_str_selector );

  //Notify movers that passes have been updated.
  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes_ ){
    tr.Debug << "Passport for " << pair.first->get_name() << " has " << pair.second->active_jumps().size()
             << " allowed jumps and " << pair.second->active_dofs().size() << " dofs." << std::endl;
    pair.first->passport_updated();
  }
}

/// @brief A brief comparator object initialized with the correct strength accessor for reuse of setup_passports
class Comparator {
public:
  Comparator( claims::ControlStrength const& (*str_access)( std::pair< claims::DOFElement, ClaimingMoverOP > const& ) ):
    str_access_( str_access ) {}
  bool operator() ( std::pair< claims::DOFElement, ClaimingMoverOP > const& a,
                    std::pair< claims::DOFElement, ClaimingMoverOP > const& b ){
    return str_access_( a ) > str_access_( b );
  }
private:
  claims::ControlStrength const& (*str_access_)( std::pair< claims::DOFElement, ClaimingMoverOP > const& );
};

void EnvClaimBroker::setup_passports( DOFElemVect& elems,
                                      claims::ControlStrength const& (*str_access)( std::pair< claims::DOFElement, ClaimingMoverOP > const& ) ) {

  // Figure out if we're doing initialization or not for output purposes.
  std::string const style = ( str_access == *init_str_selector ) ? "initialization" : "sampling";

  // Similarly to setup_initialization_passports, sort by control strength
  std::sort( elems.begin(), elems.end(), Comparator( str_access ) );
  std::map< core::id::DOF_ID, std::pair< ControlStrength, ClaimingMoverOP > > max_strength;

  //Iterate through each element. If it's EXCLUSIVE, we know not to admit any further claims.
  BOOST_FOREACH( DOFElemVect::Value pair, elems ){
    DOFElement const& element = pair.first;
    ClaimingMoverOP const& owner = pair.second;
    ControlStrength const strength = str_access( pair );

    tr.Trace << "    considering: " << element.id << " from " << owner->get_name();

    if( !element.id.valid() ){
      tr.Error << "[ERROR] During " << style << " the DOF element with id " << element.id
               << " was produced by the mover '" << owner->get_name()
               << ", which is not a valid DOF_ID object." << std::endl;
      assert( element.id.valid() );
    }

    std::pair< ControlStrength, ClaimingMoverOP > prev_str = std::make_pair( DOES_NOT_CONTROL, ClaimingMoverOP( NULL ) );

    if( max_strength.find( element.id ) != max_strength.end() ){
      prev_str = max_strength[ element.id ];
    }

    if( prev_str.first == EXCLUSIVE ){
        if( strength >= MUST_CONTROL &&
            ( owner != prev_str.second ) ){
          assert( prev_str.second );
          std::ostringstream ss;
          ss << "While broking for " << style << ", there were conflicting claim strengths on DoF "
             << element.type << " " << element.id << ": "
             << prev_str.first << " from '" << prev_str.second->get_name() << "' and "
             << strength << " from '" << owner->get_name() << "'." << std::endl;
          throw utility::excn::EXCN_BadInput( ss.str() );
        } else {
          tr.Trace << " : not assigned" << std::endl;
        }
    } else if( str_access( pair ) > DOES_NOT_CONTROL ){
      grant_access( element, owner );
      if( prev_str.first < strength )
        max_strength[ element.id ] = std::make_pair( strength, owner );
      tr.Trace << " : assigned " << std::endl;
    } else { tr.Trace << " : not assigned" << std::endl; }
  }
}

void EnvClaimBroker::grant_access( DOFElement const& e, ClaimingMoverOP owner ) const {
  using namespace core::id;
  using namespace core::kinematics::tree;

  core::environment::DofPassportOP pass = movers_and_passes_.find( owner )->second;

  pass->add_dof_access( e.id );
}

template < typename T, typename I >
utility::vector1< std::pair< T, ClaimingMoverOP > > EnvClaimBroker::collect_elements( I const& info ) const {
  typedef utility::vector1< std::pair< T, ClaimingMoverOP > > ElementList;

  utility::vector1< T > tmp;
  ElementList elements;
  std::ostringstream mover_breakdown;

  BOOST_FOREACH( EnvClaimOP claim, claims_ ){
    claim->yield_elements( info , tmp );
    if( tmp.size() > 0 ){
      mover_breakdown << "  collected " << tmp.size() << " " << T::type << " elements from "
                      << claim->type() << "Claim owned by '" << claim->owner()->get_name() << "'."
                      << std::endl;
    }
    BOOST_FOREACH( T e, tmp ){ elements.push_back( std::make_pair( e, claim->owner() ) ); }
    tmp.clear();
  }

  tr.Info << "collected " << elements.size() << " " << T::type << "." << std::endl;
  tr.Debug << "\n" << mover_breakdown.str();
  tr.Debug.flush();

  return elements;
}

EnvClaims EnvClaimBroker::collect_claims( MoverPassMap const& movers_and_passes,
                                          core::pose::Pose& pose ) {
  using namespace basic::datacache;
  using namespace core::pose::datacache;
  typedef std::map< std::string, std::set< WriteableCacheableDataOP > > DataMap;
  typedef std::set< WriteableCacheableDataOP > DataSet;

  EnvClaims claims;

  WriteableCacheableMapOP bk_map;
  if( !pose.data().has( CacheableDataType::WRITEABLE_DATA ) ){
    pose.data().set( CacheableDataType::WRITEABLE_DATA, DataCache_CacheableData::DataOP( new WriteableCacheableMap() ) );
  }

  // Create a sandboxed copy of the map
  WriteableCacheableMapOP orig_map = utility::pointer::dynamic_pointer_cast< basic::datacache::WriteableCacheableMap > ( pose.data().get_ptr( CacheableDataType::WRITEABLE_DATA ) );

  // Store the newly inserted cache objects here.
  WriteableCacheableMapOP new_cached_data( new WriteableCacheableMap() );

  //Claiming
  for( MoverPassMap::const_iterator mp_it = movers_and_passes.begin();
      mp_it != movers_and_passes.end(); ++mp_it ){

    // a modifiable sandbox_map must be passed in separately, as pose is a const&.
    WriteableCacheableMapOP sandbox_map( new WriteableCacheableMap( *orig_map ) );

    claims::EnvClaims in_claims = mp_it->first->yield_claims( pose, sandbox_map );
    BOOST_FOREACH( EnvClaimOP claim, in_claims ){
      if( claim ){
        claim->annotate( pose, ann_ );
        claims.push_back( claim );
      } else {
        std::ostringstream ss;
	EnvironmentCOP env = this->env_.lock();
        ss << "The mover '" << claim->owner() << "' yielded a null pointer as one of its "
           << in_claims.size() << " claims during broking of the environment '" << (env ? env->name() : "(Unknown)")
           << "'." << std::endl;
        throw utility::excn::EXCN_NullPointer( ss.str() );
      }
    }

    // Copy any new data from the sandbox map into the "new" map.
    for( DataMap::const_iterator subset_it = sandbox_map->begin(); subset_it != sandbox_map->end(); ++subset_it ){
      for( DataSet::const_iterator data_it = subset_it->second.begin();
           data_it != subset_it->second.end(); ++data_it  ){
        if( !orig_map->has( *data_it ) ){
          new_cached_data->insert( *data_it );
        }
      }
    }
  }

  //Write the items from the new_cache into the old cache.
  for( DataMap::const_iterator newmap_it = new_cached_data->begin();
      newmap_it != new_cached_data->end(); ++newmap_it ){
    for( DataSet::const_iterator data_it = newmap_it->second.begin();
         data_it != newmap_it->second.end(); ++data_it ){
      orig_map->insert( *data_it );
    }
  }

  this->result_.cached_data = new_cached_data;

  tr.Debug << "collected " << claims.size() << " DoFClaims from "
           << movers_and_passes.size() << " movers." << std::endl;


  return claims;
}

void EnvClaimBroker::process_elements( ResElemVect const& elems, FoldTreeSketch& fts, SizeToStringMap& new_vrts ){

  for( ResElemVect::const_iterator e_it = elems.begin(); e_it != elems.end(); ++e_it ){
    ResidueElement const& element = e_it->first;
    ClaimingMoverCOP owner = e_it->second;

    if( !element.allow_duplicates ) {
      if( ann_->has_seq_label( element.label ) ){
        std::ostringstream ss;
        ss << "[ERROR] Failed broking process due to duplicate residue claims for residue "
           << element.label << " that do not allow duplicate labels.";
        throw utility::excn::EXCN_BadInput( ss.str() );
      }
      ann_->append_seq( element.label );
      fts.append_residue();
      new_vrts[ fts.nres() ] = element.label;
      tr.Debug << "  processed duplicate-disallowed residue element request (named '"
               << element.label << "') at " << fts.nres() << std::endl;
    }
  }

  for( ResElemVect::const_iterator e_it = elems.begin(); e_it != elems.end(); ++e_it ){
    ResidueElement const& element = e_it->first;
    ClaimingMoverCOP owner = e_it->second;

    if( element.allow_duplicates ) {
      if( !ann_->has_seq_label( element.label ) ){
        ann_->append_seq( element.label );
        fts.append_residue();
        new_vrts[ fts.nres() ] = element.label;

        tr.Debug << "  processed duplicate-possible residue element request (named '"
                 << element.label << "') at " << fts.nres() << std::endl;
      } else {
        // If we're allowing duplicates and we've already added this element
        // it's cool, don't throw an exception and don't add it again.
        tr.Debug << "  ignored duplicate-possible residue element request at " << fts.nres()
                 << " ( label '" << element.label << "' already exists)." << std::endl;
      }
    }
  }
}

void EnvClaimBroker::process_elements( JumpElemVect const& elems,
                                       FoldTreeSketch& fts,
                                       JumpDataMap& new_jumps ){
  typedef std::map< std::pair< core::Size, core::Size >, BrokeredJumpDataCOP > PositionDataMap;
  PositionDataMap brokered_jumps;

  BOOST_FOREACH( JumpElemVect::Value pair, elems ){
    JumpElement const& element = pair.first;
    ClaimingMoverCOP owner = pair.second;

    Size abs_p1 = ann_->resolve_seq( element.p1 );
    Size abs_p2 = ann_->resolve_seq( element.p2 );
    std::pair< Size, Size > const pos_pair = std::make_pair( abs_p1, abs_p2 );

    // Build the JumpData -- we use this for jump overlap checks.
    BrokeredJumpDataOP new_jump( new BrokeredJumpData( pos_pair,
                                                        std::make_pair( element.atom1, element.atom2 ),
                                                        element.force_stub_intra_residue ) );

    if( !fts.has_jump( abs_p1, abs_p2 ) ){
      fts.insert_jump( abs_p1, abs_p2 );
      brokered_jumps[ pos_pair ] = new_jump;
      new_jumps[ element.label ] = new_jump;
    } else {
      tr.Debug << "Ignoring jump element for jump at " << element.p1 << "->"
               << element.p2 << " because it already exists " << std::endl;

      if( brokered_jumps.find( pos_pair ) != brokered_jumps.end() &&
          !brokered_jumps[ pos_pair ]->operator==( *new_jump ) ){
        std::ostringstream ss;

        new_jump->operator<<( ss );

        ss << "The broker element claiming a jump between " << element.p1 << " and " << element.p2
           << " overlaps with an existing jump between absolute positions " << pos_pair.first
           << " and " << pos_pair.second << ". Normally, this would be ok (the jumps would be "
           << "the same jump and access would be offered to both ClaimingMovers), but the jump data "
           << "differs. This jump was: " << std::endl;

        new_jump->operator<<( ss );

        ss << std::endl << " whereas the existing was :" << std::endl;
        brokered_jumps[ pos_pair ]->operator<<( ss );
        ss << std::endl;

        throw utility::excn::EXCN_BadInput( ss.str() );
      }

      // if we make it through the above check, it's safe to double-label this jump with both.
      BrokeredJumpDataCOP old_jump_data = brokered_jumps.find( pos_pair )->second;
      new_jumps[ element.label ] = old_jump_data;
    }

    tr.Debug << "  processed jump element for jump " << abs_p1 << "->" << abs_p2 << std::endl;
  }
}

void EnvClaimBroker::process_elements( CutElemVect const& elems,
                                       FoldTreeSketch& fts ){
  BOOST_FOREACH( CutElemVect::Value pair, elems ){
    CutElement element = pair.first;

    Size abs_p = ann_->resolve_seq( element.p );
    try {
      fts.insert_cut( abs_p );
    } catch ( core::environment::EXCN_FTSketchGraph& excn ){
      std::ostringstream ss;
      EnvironmentCOP env = env_.lock();
      ss << "The Environment '" << (env ? env->name() : "(Unknown)" ) << "' had more than one"
         << " cut placed at absolute position " << abs_p << "." << std::endl;
      excn.add_msg( ss.str() );
      throw excn;
    }
  }
}

void EnvClaimBroker::process_elements( CutBiasElemVect const& elems, BiasVector& bias ){

  BOOST_FOREACH( CutBiasElemVect::Value pair, elems ){
    CutBiasElement element = pair.first;

    Size abs_p = ann_->resolve_seq( element.p );
    if( element.bias > 1 || element.bias < 0 ){
      std::ostringstream ss;
      ss << "Cut biases must be between 0 and 1. Cut was (" << element.p << "," << bias << ")";
      throw utility::excn::EXCN_RangeError( ss.str() );
    }
    bias[ abs_p ] *= element.bias;
  }
}

void EnvClaimBroker::add_chainbreak_variants( core::Size rsd_num_lower,
                                              core::conformation::Conformation& conf ) const {
  using namespace core::chemical;
  using namespace core::conformation;

  assert( conf.fold_tree().is_cutpoint( (int) rsd_num_lower ) );

  Residue const& rsd_lower( conf.residue( rsd_num_lower ) );
  Residue const& rsd_upper( conf.residue( rsd_num_lower + 1 ) );
  ResidueTypeSet const& rsd_set( rsd_lower.residue_type_set() );

  ResidueType const& new_type_lower( rsd_set.get_residue_type_with_variant_added( rsd_lower.type(),
                                                                                  CUTPOINT_LOWER ) );
  ResidueType const& new_type_upper( rsd_set.get_residue_type_with_variant_added( rsd_upper.type(),
                                                                                  CUTPOINT_UPPER ) );

  ResidueOP new_lower( ResidueFactory::create_residue( new_type_lower, rsd_lower, conf ) );
  ResidueOP new_upper( ResidueFactory::create_residue( new_type_upper, rsd_upper, conf ) );

  copy_residue_coordinates_and_rebuild_missing_atoms( rsd_lower, *new_lower, conf );
  copy_residue_coordinates_and_rebuild_missing_atoms( rsd_upper, *new_upper, conf );

  conf.replace_residue( rsd_lower.seqpos(), *new_lower, true );
  conf.replace_residue( rsd_upper.seqpos(), *new_upper, true );

}

EnvClaimBroker::BrokeredJumpData::BrokeredJumpData( std::pair< core::Size, core::Size > const& positions,
                                                    std::pair< std::string, std::string > const& atoms,
                                                    bool put_jump_stub_intra_residue  ) :
  pos( positions ),
  atoms( atoms ),
  put_jump_stub_intra_residue( put_jump_stub_intra_residue )
{
  std::pair< std::string, std::string > null_atoms = std::make_pair( "", "" );
  if( put_jump_stub_intra_residue &&
      atoms != null_atoms ){
    tr.Warning << "Jump at aboslute position (" << pos.first << "," << pos.second
               << ") made explicit atom choices " << atoms.first << " and "
               << atoms.second << ". These choices may be overwritten by "
               << "FoldTree::put_jump_stubs_intra_residue." << std::endl;
  }
}

bool EnvClaimBroker::BrokeredJumpData::operator==( BrokeredJumpData const& other ) const {
  if( this == &other ) return true;

  return ( other.atoms == this->atoms ) &&
         ( other.pos == this->pos ) &&
         ( this->put_jump_stub_intra_residue == other.put_jump_stub_intra_residue );
}

std::ostream& EnvClaimBroker::BrokeredJumpData::operator<<( std::ostream& os ) const {
  os << "BrokeredJumpData: position = {" << pos.first << "," << pos.second << "}, atoms = {"
     << atoms.first << "," << atoms.second << "}, put_stubs_intra_residue = "
     << ( put_jump_stub_intra_residue ? "true" : "false" );
  return os;
}


} // environment
} // protocols
