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

static basic::Tracer tr("protocols.environment.EnvClaimBroker", basic::t_info);

namespace protocols {
namespace environment {

using namespace protocols::environment::claims;

using core::environment::LocalPosition;

using core::environment::SequenceAnnotation;
using core::environment::SequenceAnnotationOP;
using core::environment::SequenceAnnotationCOP;

using core::environment::FoldTreeSketch;

using core::conformation::Conformation;

EnvClaimBroker::EnvClaimBroker( Environment const& env,
                                MoverPassMap const& movers_and_passes,
                                core::pose::Pose const& in_pose ):
  movers_and_passes_( movers_and_passes )
{
  core::conformation::ConformationCOP in_conf = &in_pose.conformation();
  core::pose::Pose pose( in_pose );
  claims_ = collect_claims( movers_and_passes, pose );

  // Initialize Annotations object:
  ProtectedConformationCOP prot_conf( dynamic_cast< ProtectedConformation const* >( in_conf.get() ) );
  if( prot_conf ){
    // Keep old annotations by copy-constructing a new annotations object.
    ann_ = new SequenceAnnotation( *( prot_conf->annotations() ) );
  } else {
    // Unprotected Conformation objects don't have annotations (yet?). We then build the annotation object.
    ann_ = new SequenceAnnotation( in_conf->size() );
  }

  // Build a temporary conformation for which we manipulate the fold tree. After the fold tree is set,
  // we "seal it" as a ProtectedConformation for initializing DoFs.
  Conformation tmp_conf( *in_conf );

  broker_fold_tree( tmp_conf, pose.data() );

  ProtectedConformationOP conf = new ProtectedConformation( &env, tmp_conf );
  conf->attach_annotation( ann_ );

  broker_dofs( conf );

  if( pose.data().has( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ){
    basic::datacache::WriteableCacheableMapOP data_map = dynamic_cast< basic::datacache::WriteableCacheableMap* >( pose.data().get_raw_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) );
    assert( data_map );
    pose.set_new_conformation( conf );
    pose.data().set( core::pose::datacache::CacheableDataType::WRITEABLE_DATA, data_map );
  } else {
    pose.set_new_conformation( conf );
  }

  result_.pose = pose;

  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes ){
    pair.first->broking_finished( result() );
  }

  assert( dynamic_cast< basic::datacache::WriteableCacheableMap* >( pose.data().get_raw_ptr( core::pose::datacache::CacheableDataType::WRITEABLE_DATA ) ) );
}

EnvClaimBroker::~EnvClaimBroker() {}

EnvClaimBroker::BrokerResult const& EnvClaimBroker::result() const {
  return result_;
}

void EnvClaimBroker::broker_fold_tree( Conformation& conf,
                                       basic::datacache::BasicDataCache& datacache ){
  core::environment::FoldTreeSketch fts( conf.size() );
  core::environment::FoldTreeSketch physical_fts( conf.size() );

  //FTElements: virtual residues ---------------------------------------------------------------
  ResidueElements r_elems;
  SizeToStringMap new_vrts;
  collect_elements( r_elems, fts );
  process_elements( r_elems, fts, new_vrts );
  result_.new_vrts = new_vrts;

  //FTElements: jumps --------------------------------------------------------------------------
  JumpElements j_elems;
  collect_elements( j_elems, fts );
  StringToSizePairMap new_jumps;
  StringToStringPairMap jump_atoms;
  process_elements( j_elems, fts, physical_fts, new_jumps, jump_atoms );

  //FTElements cut -----------------------------------------------------------------------------
  CutElements c_elems;
  collect_elements( c_elems, fts );
  std::set< core::Size > unphysical_cuts; //stores residue number of unphysical cuts
  process_elements( c_elems, fts, physical_fts, unphysical_cuts );

  //FTElements CutBias -------------------------------------------------------------------------
  CutBiasElements cb_elems;
  collect_elements( cb_elems, fts );
  utility::vector1< core::Real > bias( fts.nres(), 1.0 );
  process_elements( cb_elems, bias );

  //Render FoldTree ----------------------------------------------------------------------------
  core::kinematics::FoldTreeOP ft = render_fold_tree( fts, unphysical_cuts, bias, datacache );
  core::kinematics::FoldTreeOP physical_ft; // used to tell the loop closer what to do later.

  try {
    physical_ft = physical_fts.render();
  } catch( core::environment::EXCN_FTSketchGraph e ){
    tr.Error << "[ERROR] The physical fold tree construction failed. This happens when the EnvClaimBroker "
             << "can't figure out which cuts are physical (scored and removed at Environment close) and "
             << "aren't. Check your jump claims." << std::endl;
    throw;
  }

  annotate_fold_tree( ft, new_jumps, jump_atoms, ann_ );
  annotate_fold_tree( physical_ft, new_jumps, jump_atoms );

  ft->put_jump_stubs_intra_residue();
  physical_ft->put_jump_stubs_intra_residue();

  result_.closer_ft = physical_ft;

  add_virtual_residues( conf, new_vrts, ann_ );

  conf.fold_tree( *ft );

  if( tr.Debug.visible() ){
    tr.Debug << "EnvClaimBroker finished fold tree construction:" << *ft << std::endl;
  }

  //Swap out unphysical cut residues for chainbreak variants
  BOOST_FOREACH( Size cut , unphysical_cuts ){
    add_chainbreak_variants( cut, conf );
    result_.inserted_cut_variants.insert( cut );
  }
}

core::kinematics::FoldTreeOP
EnvClaimBroker::render_fold_tree( FoldTreeSketch& fts,
                                  std::set< Size >& unphysical_cuts,
                                  utility::vector1< core::Real > const& bias,
                                  basic::datacache::BasicDataCache& datacache ) {
  using namespace basic::datacache;
  using namespace core::pose::datacache;

  std::set< Size > auto_cuts;

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
    datacache.set( CacheableDataType::WRITEABLE_DATA, new WriteableCacheableMap() );
  }

  WriteableCacheableMapOP datamap = dynamic_cast< WriteableCacheableMap* >( datacache.get_raw_ptr( CacheableDataType::WRITEABLE_DATA ) );
  bool found_data = false;

  if( datamap && datamap->find( "AutoCutData" ) != datamap->end() ){
    BOOST_FOREACH( basic::datacache::WriteableCacheableDataOP data,
                   datamap->find( "AutoCutData" )->second ) {
      AutoCutDataOP auto_cut_data = dynamic_cast< AutoCutData* >( data.get() );
      assert( auto_cut_data );

      if( auto_cut_data->hash() == hash ){
        tr.Debug << "Found appropriate hash-valued AutoCutData. Repeating cut choices." << std::endl;
        found_data = true;

        BOOST_FOREACH( core::Size cut, auto_cut_data->cuts() ){
          auto_cuts.insert( cut );
          fts.insert_cut( cut );
        }

      } else {
        tr.Debug << "Rejecting AutoCutData because its hash " << auto_cut_data->hash()
                 << " does not match EnvClaimBroker hash " << hash << std::endl;
      }
    }
  } else {
    tr.Debug << "No AutoCutData found, generating cuts automatically." << std::endl;
  }

  if( !found_data ){
    // If we get here, the DataCache didn't include an appropriate AutoCutData instance.
    // Automatically generate cut choices. (This is the usually the case.)

    auto_cuts = fts.remove_cycles( bias ); //auto_cuts stores residue numbers
                                           // automatic cuts are assumed to be unphysical => add to list of unphysical cuts
    tr.Debug << "automatically inserted " << auto_cuts.size() << " cuts into fold tree." << std::endl;

    AutoCutDataOP data = new AutoCutData( hash, auto_cuts );
    datamap->insert( data );
  } else {
    //If we repeated cut choices, there should be no cuts.
    assert( fts.cycle().empty() );
  }

  std::copy( auto_cuts.begin(), auto_cuts.end(), std::inserter( unphysical_cuts, unphysical_cuts.begin() ) );

  try {
    return fts.render();
  } catch ( ... ) {
    tr.Error << "[ERROR] A problem was encountered rendering fold tree"
             << "choices in EnvClaimBroker::render_fold_tree." << std::endl;
    throw;
  }
}

void EnvClaimBroker::annotate_fold_tree( core::kinematics::FoldTreeOP ft,
                                         StringToSizePairMap const& jump_labels,
                                         StringToStringPairMap const& jump_atoms,
                                         SequenceAnnotationOP ann ) {

  //Bind jump numbers to labels in annotations
  BOOST_FOREACH( StringToSizePairMap::value_type pair, jump_labels ){
    std::string const& label = pair.first;
    std::pair< Size, Size > const& jump_resids = pair.second;

    Size jump_id = ft->jump_nr( jump_resids.first, jump_resids.second );

    if( jump_id != 0 ){ // The physical fold tree won't have some jumps, so we don't want to fail.
      if( ann ){ // Physical fold tree doesn't use the annotations.
        ann->add_jump_label( label, jump_id );
      }

      std::string const& a1 = jump_atoms.find( label )->second.first;
      std::string const& a2 = jump_atoms.find( label )->second.second;

      ft->set_jump_atoms( (int) jump_id, a1, a2,
                          (a1 == "" && a2 == "") );
    }
  }
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

void EnvClaimBroker::broker_dofs( ProtectedConformationOP conf ){

  std::set< ClaimingMoverOP > initializers;

  //Broker Arbitrary DOFs ---------------------------------------------------------------------------
  DOFElements d_elems;
  collect_elements( d_elems, conf );
  setup_initialization_passports( d_elems, initializers );

  // Once all initialization passports for both RTs and Torsions are configured, allow
  // movers to initialize

  core::pose::Pose p;
  p.set_new_conformation( conf );

  // Now that passports are properly configured, allow each mover to initialize.
  // As soon as each mover initializes, access is revoked, so it can be built again
  // for the sampling phase.
  BOOST_FOREACH( ClaimingMoverOP mover, initializers ){
    mover->passport_updated();
    mover->initialize( p );
    movers_and_passes_.find( mover )->second->revoke_all_access();
  }

  //Copy the changes from the initialized conformation into conf
  *conf = *( static_cast< ProtectedConformation const* >( &p.conformation() ) );

  setup_control_passports( d_elems );

  //Notify movers that passes have been updated.
  BOOST_FOREACH( MoverPassMap::value_type pair, movers_and_passes_ ){
    pair.first->passport_updated();
  }
}

template < typename T >
bool compare_init_strength( T const& a, T const& b ){
  return a.i_str > b.i_str;
}

bool is_initialized( std::set< std::string > const& initalized,
                     DOFElement e  ){
  return initalized.find( utility::to_string( e.id ) ) != initalized.end();
}

///@brief iterate through all dof elements, determine who gets initialization right, and modify
///       their DofPassport appropriately
template < typename T >
void EnvClaimBroker::setup_initialization_passports( utility::vector1< T >& elems,
                                                     std::set< ClaimingMoverOP >& initializers ){

  // Sort elements in descending strength order to guarantee highest prio gets first dibs.
  std::sort( elems.begin(), elems.end(), compare_init_strength< T > );
  std::map< core::id::DOF_ID, ClaimingMoverOP > initialized;

  BOOST_FOREACH( T element, elems ){
    std::map< core::id::DOF_ID, ClaimingMoverOP >::iterator prev_init = initialized.find( element.id );
    if( prev_init != initialized.end() && //we already assigned an initializer (of equal or higher prio) to this position
        element.i_str == MUST_INITIALIZE && //this initializer demanded to initialize this position
        initialized[ element.id ] != element.owner ){ //initializers are different (i.e. true conflict)
      std::string msg = "Cannot initialize " + element.type + " with id " + utility::to_string( element.id )
                        + " due conflicting init strengths (e.g. >1 MUST_INITIALIZE claim)."
                        + "  Conflicting initializers were: " + initialized[ element.id ]->get_name() + " and "
                        + element.owner->get_name();
      throw utility::excn::EXCN_BadInput( msg );
    } else if( initialized.find( element.id ) == initialized.end() &&
               element.i_str > DOES_NOT_INITIALIZE ) {
      grant_access( element );
      initialized[ element.id ] = element.owner;
      initializers.insert( element.owner );
    } else { /* do not grant access */ }
  }

  tr.Info << "  setting " << initialized.size() << " " << T::type
          << "-type DOF initial states." << std::endl;
}

template < typename T >
bool compare_ctrl_strength( T const& a, T const& b ){
  return a.c_str > b.c_str;
}

template < typename T >
void EnvClaimBroker::setup_control_passports( utility::vector1< T >& elems ) {

  // Similarly to setup_initialization_passports, sort by control strength
  std::sort( elems.begin(), elems.end(), compare_ctrl_strength< T > );
  std::map< core::id::DOF_ID, ControlStrength > max_strength;

  //Iterate through each element. If it's EXCLUSIVE, we know not to admit any further claims.
  BOOST_FOREACH( T element, elems ){
    ControlStrength prev_str = DOES_NOT_CONTROL;
    if( max_strength.find( element.id ) != max_strength.end() ){
      prev_str = max_strength[ element.id ];
    }

    if( prev_str == EXCLUSIVE ){
        if( element.c_str >= MUST_CONTROL ){
          std::string msg = "Cannot broker " + element.type + " " + utility::to_string( element.id ) +
                            "due conflicting control strengths (e.g. >1 EXCLUSIVE claim).";
          throw utility::excn::EXCN_BadInput( msg );
        } else {
          // do not grant access
        }
    } else if( element.c_str > DOES_NOT_CONTROL ){
      grant_access( element );
      max_strength[ element.id ] = std::max( element.c_str, prev_str );
    } else { /* do not grant access */ }
  }
}

void EnvClaimBroker::grant_access( DOFElement const& e ) const {
  using namespace core::id;
  using namespace core::kinematics::tree;

  core::environment::DofPassportOP pass = movers_and_passes_.find( e.owner )->second;

  pass->add_dof_access( e.id );
}

template < typename T, typename I >
void EnvClaimBroker::collect_elements( utility::vector1< T >& elements,
                                       I const& info ) const {
  typedef utility::vector1< T > ElementList;

  ElementList tmp;

  BOOST_FOREACH( EnvClaimOP claim, claims_ ){
    claim->yield_elements( info , tmp );
    std::copy( tmp.begin(), tmp.end(), std::back_inserter( elements ) );
    tmp.clear();
  }

  tr.Info << "collected " << elements.size() << " " << T::type << "." << std::endl;
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
    pose.data().set( CacheableDataType::WRITEABLE_DATA, new WriteableCacheableMap() );
  }

  // Create a sandboxed copy of the map
  WriteableCacheableMapOP orig_map = dynamic_cast< WriteableCacheableMap* >( pose.data().get_ptr( CacheableDataType::WRITEABLE_DATA ).get() );

  // Store the newly inserted cache objects here.
  WriteableCacheableMapOP new_cached_data = new WriteableCacheableMap();

  //Claiming
  for( MoverPassMap::const_iterator mp_it = movers_and_passes.begin();
      mp_it != movers_and_passes.end(); ++mp_it ){

    // a modifiable sandbox_map must be passed in separately, as pose is a const&.
    WriteableCacheableMapOP sandbox_map = new WriteableCacheableMap( *orig_map );

    claims::EnvClaims in_claims = mp_it->first->yield_claims( pose, sandbox_map );
    std::copy( in_claims.begin(), in_claims.end(), std::back_inserter( claims ) );

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

  this->result_.cached_data = WriteableCacheableMapCOP( new_cached_data.get() );

  tr.Debug << "collected " << claims.size() << " DoFClaims from "
           << movers_and_passes.size() << " movers." << std::endl;


  return claims;
}



void EnvClaimBroker::process_elements( ResidueElements const& elems, FoldTreeSketch& fts, SizeToStringMap& new_vrts ){
  for( ResidueElements::const_iterator e_it = elems.begin(); e_it != elems.end(); ++e_it ){
    fts.append_residue();
    ann_->append_seq( e_it->label );
    new_vrts[ fts.nres() ] = e_it->label;
  }
}

void EnvClaimBroker::process_elements( JumpElements const& elems,
                                       FoldTreeSketch& fts,
                                       FoldTreeSketch& physical_fts,
                                       StringToSizePairMap& new_jumps,
                                       StringToStringPairMap& jump_atoms ){
  BOOST_FOREACH( JumpElement element, elems ){
    Size abs_p1 = ann_->resolve_seq( element.p1 );
    Size abs_p2 = ann_->resolve_seq( element.p2 );

    if( !fts.has_jump( abs_p1, abs_p2 ) ){
      fts.insert_jump( abs_p1, abs_p2 );
      if( element.has_physical_cut ){
        physical_fts.insert_jump( abs_p1, abs_p2 );
      }

      new_jumps[  element.label ] = std::make_pair( abs_p1, abs_p2 );
      jump_atoms[ element.label ] = std::make_pair( element.atom1, element.atom2 );
    } else {
      tr.Debug << "Ignoring jump element for jump at " << element.p1 << "->"
               << element.p2 << " because it already exists " << std::endl;
    }

    tr.Debug << "  processed jump element for jump " << abs_p1 << "->" << abs_p2 << std::endl;
  }
}

void EnvClaimBroker::process_elements( CutElements const& elems,
                                       FoldTreeSketch& fts,
                                       FoldTreeSketch& physical_fts,
                                       std::set< core::Size >& unphysical_cuts ){
  BOOST_FOREACH( CutElement element, elems ){
    Size abs_p = ann_->resolve_seq( element.p );
    fts.insert_cut( abs_p );
    if( element.physical ){
      physical_fts.insert_cut( abs_p );
    } else {
      unphysical_cuts.insert( abs_p );
    }
  }
}

void EnvClaimBroker::process_elements( CutBiasElements const& elems, BiasVector& bias ){
  for( CutBiasElements::const_iterator e_it = elems.begin(); e_it != elems.end(); ++e_it ){
    Size abs_p = ann_->resolve_seq( e_it->p );
    if( e_it-> bias > 1 || e_it->bias < 0 ){
      throw utility::excn::EXCN_RangeError( "Cut biases must be between 0 and 1. Cut was ("+
                                           utility::to_string( e_it->p )+","+
                                           utility::to_string( bias )+")" );
    }
    bias[ abs_p ] *= e_it->bias;
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

} // environment
} // protocols
