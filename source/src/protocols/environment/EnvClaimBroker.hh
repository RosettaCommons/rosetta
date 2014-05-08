// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/environment/EnvClaimBroker.hh
/// @brief A class responsible for managing claims from ClaimingMovers, and retaining important information about the result of the brokering process.
///
/// @author Justin Porter

#ifndef INCLUDED_protocols_environment_EnvClaimBroker_hh
#define INCLUDED_protocols_environment_EnvClaimBroker_hh

// Unit Headers
#include <protocols/environment/EnvClaimBroker.fwd.hh>

// Package headers
#include <core/environment/DofPassport.fwd.hh>
#include <core/environment/SequenceAnnotation.fwd.hh>
#include <core/environment/FoldTreeSketch.hh>

#include <protocols/environment/Environment.fwd.hh>
#include <protocols/environment/ClaimingMover.fwd.hh>
#include <protocols/environment/ProtectedConformation.fwd.hh>

#include <protocols/environment/claims/BrokerElements.hh>
#include <protocols/environment/claims/EnvClaim.fwd.hh>

// Project headers
#include <basic/datacache/WriteableCacheableMap.fwd.hh>

#include <core/conformation/Conformation.hh>

#include <core/pose/Pose.hh>

#include <core/kinematics/FoldTree.fwd.hh>
#include <core/types.hh>

#include <utility/pointer/ReferenceCount.hh>

// C++ Headers
#include <map>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

class EnvClaimBroker : public utility::pointer::ReferenceCount {
  typedef std::map< ClaimingMoverOP, core::environment::DofPassportOP > MoverPassMap;
  typedef core::environment::SequenceAnnotationCOP SequenceAnnotationCOP;
  typedef core::environment::SequenceAnnotationOP SequenceAnnotationOP;
  typedef core::environment::FoldTreeSketch FoldTreeSketch;

  typedef core::conformation::Conformation Conformation;
  typedef core::conformation::ConformationOP ConformationOP;
  typedef core::conformation::ConformationCOP ConformationCOP;

  typedef std::map< core::Size, std::string > SizeToStringMap;
  typedef std::map< std::string, std::pair< Size, Size > > StringToSizePairMap;
  typedef std::map< std::string, std::pair< std::string, std::string > > StringToStringPairMap;
  typedef utility::vector1< core::Real > BiasVector;

public:
  EnvClaimBroker( Environment const& env,
                  MoverPassMap const& movers_and_passes,
                  core::pose::Pose const& in_pose );

  virtual ~EnvClaimBroker();

  struct BrokerResult {
    core::pose::Pose pose;
    SizeToStringMap new_vrts;
    core::kinematics::FoldTreeOP closer_ft;
    std::set< core::Size > inserted_cut_variants;
    basic::datacache::WriteableCacheableMapCOP cached_data;
  };

  BrokerResult const& result() const;

private:

  void broker_fold_tree( Conformation&, basic::datacache::BasicDataCache& );

  void broker_dofs( ProtectedConformationOP );

  core::kinematics::FoldTreeOP render_fold_tree( FoldTreeSketch& fts,
                                                 std::set< Size >& unphysical_cuts,
                                                 utility::vector1< core::Real > const& bias,
                                                 basic::datacache::BasicDataCache& );

  void annotate_fold_tree( core::kinematics::FoldTreeOP,
                           StringToSizePairMap const& jump_labels,
                           StringToStringPairMap const& jump_atoms,
                           SequenceAnnotationOP = NULL );

  void add_virtual_residues( Conformation&,
                             SizeToStringMap const& new_vrts,
                             SequenceAnnotationOP );

  void process_elements( claims::ResidueElements const& elems,
                         FoldTreeSketch& fts,
                         SizeToStringMap& new_vrts );

  void process_elements( claims::JumpElements const& elems,
                         FoldTreeSketch& fts,
                         FoldTreeSketch& physical_fts,
                         StringToSizePairMap& new_jumps,
                         StringToStringPairMap& jump_atoms );

  void process_elements( claims::CutElements const& elems,
                         FoldTreeSketch& fts,
                         FoldTreeSketch& physical_fts,
                         std::set< core::Size >& unphysical_cuts );

  void process_elements( claims::CutBiasElements const& elems, BiasVector& bias );

  void grant_access( claims::DOFElement const& e ) const;

  template < typename T >
  void setup_control_passports( utility::vector1< T >& elements );

  template < typename T >
  void setup_initialization_passports( utility::vector1< T >& elems,
                                       std::set< ClaimingMoverOP >& initializers );

  template < typename T, typename I >
  void collect_elements( utility::vector1< T >& elements,
                         I const& info ) const;

  claims::EnvClaims collect_claims( MoverPassMap const& movers_and_passes,
                                    core::pose::Pose& pose );

  ///@brief broker new residues
  void build_new_residues( claims::EnvClaims const& claims, FoldTreeSketch& fts, SequenceAnnotationOP ann );

  ///@brief use accepted claims to build DofPassport objects for movers.
  void assign_passports( claims::EnvClaims const&, ProtectedConformation const& );

  void add_chainbreak_variants( core::Size cut, core::conformation::Conformation& conf ) const;

  MoverPassMap const& movers_and_passes_;
  SequenceAnnotationOP ann_;
  claims::EnvClaims claims_;
  BrokerResult result_;

}; // EnvClaimBroker

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_EnvClaimBroker_hh
