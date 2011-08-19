// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/ScoreFunction.cc
/// @brief  Atom pair energy functions
/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)
/// @author Kevin P. Hinshaw (KevinHinshaw@gmail.com)
/// @author Christopher Miles (cmiles@uw.edu)

// Unit headers
#include <core/scoring/methods/LinearChainbreakEnergy.hh>
#include <core/scoring/methods/LinearChainbreakEnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/NamedAtomID.hh>

// utility headers
#include <basic/Tracer.hh>

// C++ headers
#include <iostream>

using core::Size;

static basic::Tracer tr("core.scoring.LinearChainbreak",basic::t_info);

namespace core {
namespace scoring {
namespace methods {

  LinearChainbreakEnergy::LinearChainbreakEnergy()
    : parent(new LinearChainbreakEnergyCreator) {
    initialize(core::SZ_MAX);
  }

  LinearChainbreakEnergy::LinearChainbreakEnergy(Size allowable_sequence_sep)
    : parent(new LinearChainbreakEnergyCreator) {
    initialize(allowable_sequence_sep);
  }

  LinearChainbreakEnergy::LinearChainbreakEnergy(const LinearChainbreakEnergy& o)
    : parent(o) {
    // assignment to non-cache-related variables
    allowable_sequence_sep_ = o.allowable_sequence_sep_;

    // assign cache-related variables to their default values. this will trigger
    // a cache miss on their first use, but greatly simplifies the design and
    // improves reasoning about code paths in the event of copying/cloning
    shortest_paths_.reset();
    previous_hash_value_ = 0;
  }

  LinearChainbreakEnergy& LinearChainbreakEnergy::operator=(const LinearChainbreakEnergy& o) {
    if (this != &o) {
      // base class assignments
      parent::operator=(o);

      // assignment to non-cache-related variables
      allowable_sequence_sep_ = o.allowable_sequence_sep_;

      // assign cache-related variables to their default values. this will trigger
      // a cache miss on their first use, but greatly simplifies the design and
      // improves reasoning about code paths in the event of copying/cloning
      shortest_paths_.reset();
      previous_hash_value_ = 0;
    }
    return *this;
  }

  LinearChainbreakEnergy::~LinearChainbreakEnergy() {}

  void LinearChainbreakEnergy::initialize(Size allowable_sequence_sep) {
    allowable_sequence_sep_ = allowable_sequence_sep;
    previous_hash_value_ = 0;
    shortest_paths_.reset();
  }

  /// called at the end of energy evaluation
  /// In this case (LinearChainbreakEnergy), all the calculation is done here
  void LinearChainbreakEnergy::finalize_total_energy(pose::Pose& pose,
                                                     ScoreFunction const &,
                                                     EnergyMap& totals) const {
    using conformation::Residue;
    using core::kinematics::ShortestPathInFoldTree;

    Real total_dev(0.0);
    Real total_ovp(0.0);

    // The cached ShortestPathInFoldTree instance is invalid and must be recomputed
    size_t hash_value = pose.fold_tree().hash_value();
    if (!shortest_paths_.get() || previous_hash_value_ != hash_value) {
      shortest_paths_.reset(new ShortestPathInFoldTree(pose.fold_tree()));
      previous_hash_value_ = hash_value;
    }

    for ( int n=1; n<= pose.fold_tree().num_cutpoint(); ++n ) {
      int const cutpoint( pose.fold_tree().cutpoint( n ) );
      Residue const & lower_rsd( pose.residue( cutpoint ) );
      if ( !lower_rsd.has_variant_type( chemical::CUTPOINT_LOWER ) ) continue;

      Residue const & upper_rsd( pose.residue( cutpoint+1 ) );
      // Don't stop the calculation by assert, because there is a case that
      // chain needs to be circularized. Nobu
      // assert( upper_rsd.has_variant_type( chemical::CUTPOINT_UPPER ) );
      if( !upper_rsd.has_variant_type( chemical::CUTPOINT_UPPER ) ) continue;
      Size const nbb( lower_rsd.mainchain_atoms().size() );

      // Determine whether the separation between <lowed_rsd> and <upper_rsd>,
      // as computed by ShortestPathInFoldTree, exceeds the current allowable
      // sequence separation
      Size separation = shortest_paths_->dist(cutpoint, cutpoint + 1);
      if (separation > allowable_sequence_sep_) {
        tr.Debug << "Chainbreak skipped-- "
                 << separation << " > " << allowable_sequence_sep_
                 << std::endl;
        continue;
      }

      // virtual N and CA on lower_rsd (OVL1 and OVL2) and a virtual C on upper_rsd "OVU1"
      total_dev += (upper_rsd.atom( upper_rsd.mainchain_atoms()[  1] ).xyz().distance( lower_rsd.atom( "OVL1"
 ).xyz() ) +
                    upper_rsd.atom( upper_rsd.mainchain_atoms()[  2] ).xyz().distance( lower_rsd.atom( "OVL2"
 ).xyz() ) +
                    lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz().distance( upper_rsd.atom( "OVU1"
 ).xyz() ) );

      // now compute the overlap score: this is done by comparing the stub   (lower side )  C, | N, CA (upper side)
      // the stub for the lower_rsd is already in the atom-tree at atom "OVL2 == CA*"
      kinematics::Stub lower_stub(pose.conformation().atom_tree().atom(id::AtomID( pose.residue(cutpoint).atom_index("OVL2"), cutpoint) ).get_stub() );

      // the upper stub... let's just compute it for now...
      // could be gained from AtomID( NamedAtomID( "OVU1", cutpoint+1 ).get_stub()
      // and then the correct reversal...
      kinematics::Stub upper_stub(upper_rsd.atom( upper_rsd.mainchain_atoms()[  2]).xyz(), //CA
                                  upper_rsd.atom( upper_rsd.mainchain_atoms()[  1]).xyz(),  //N
                                  upper_rsd.atom( "OVU1" ).xyz()); //virtual C

      //for double-checking... ( debug )
      kinematics::Stub manual_lower_stub(lower_rsd.atom( "OVL2" ).xyz(),  //virtual CA
                                         lower_rsd.atom( "OVL1" ).xyz(),  //virtual N
                                         lower_rsd.atom( lower_rsd.mainchain_atoms()[ nbb] ).xyz()); // C

      total_ovp +=
        manual_lower_stub.M.col_x().distance( upper_stub.M.col_x() ) +
        manual_lower_stub.M.col_y().distance( upper_stub.M.col_y() ) +
        manual_lower_stub.M.col_z().distance( upper_stub.M.col_z() );

      if ( distance( lower_stub, manual_lower_stub ) > 0.01 ) {
        tr.Warning << "WARNING: mismatch between manual computed and atom-tree stub: "
                   << lower_stub << " " << manual_lower_stub << std::endl;
      }
    }
    assert( std::abs( totals[ linear_chainbreak ] ) < 1e-3 );

    totals[ linear_chainbreak ] = total_dev/3.0;
    //division by 3 because we have 3 distances. Fix your weights to perfection.
    totals[ overlap_chainbreak ] = total_ovp;
  }

  /// called during gradient-based minimization inside dfunc
  /**
     F1 and F2 are not zeroed -- contributions from this atom are
     just summed in
  **/
  void LinearChainbreakEnergy::eval_atom_derivative(id::AtomID const & id,
                                                    pose::Pose const & pose,
                                                    kinematics::DomainMap const &, // domain_map,
                                                    ScoreFunction const &, // sfxn,
                                                    EnergyMap const & weights,
                                                    Vector & F1,
                                                    Vector & F2) const {
    using conformation::Residue;
    using chemical::CUTPOINT_LOWER;
    using chemical::CUTPOINT_UPPER;

    if ( pose.fold_tree().is_cutpoint( id.rsd() ) && pose.residue(id.rsd() ).has_variant_type( CUTPOINT_LOWER ) ) {
      Residue const & lower_rsd( pose.residue( id.rsd()     ) );
      Residue const & upper_rsd( pose.residue( id.rsd() + 1 ) );
      Vector const & xyz_moving( pose.xyz( id ) );
      bool match( true );
      Vector xyz_fixed;
      Size const nbb( lower_rsd.mainchain_atoms().size() );
      if ( id.atomno() == lower_rsd.mainchain_atoms()[nbb] ) {
        xyz_fixed = upper_rsd.atom( "OVU1" ).xyz();
      } else if ( lower_rsd.atom_name( id.atomno() ) == "OVL1" ) {
        xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[1] ).xyz();
      } else if ( lower_rsd.atom_name( id.atomno() ) == "OVL2" ) {
        xyz_fixed = upper_rsd.atom( upper_rsd.mainchain_atoms()[2] ).xyz();
      } else {
        match = false;
      }

      if ( match ) {
        // deriv = 1
        // factor = deriv / r = 1/r //sounds ugly!
        Vector const f2 ( xyz_moving - xyz_fixed );
        Real const dist ( f2.length() );
        if ( dist >= 0.01 ) { //avoid getting too close to singularity...
          Real const invdist( 1.0 / dist );
          F1 += weights[ linear_chainbreak ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
          F2 += weights[ linear_chainbreak ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
        }
      }
    }

    if ( id.rsd() > 1 && pose.fold_tree().is_cutpoint( id.rsd()-1 ) &&
         pose.residue(id.rsd() ).has_variant_type( CUTPOINT_UPPER )  ) {
      Residue const & lower_rsd( pose.residue( id.rsd() - 1 ) );
      Residue const & upper_rsd( pose.residue( id.rsd()     ) );
      Vector const & xyz_moving( pose.xyz( id ) );
      bool match( true );
      Vector xyz_fixed;
      Size const nbb( lower_rsd.mainchain_atoms().size() );
      if ( id.atomno() == upper_rsd.mainchain_atoms()[1] ) {
        xyz_fixed = lower_rsd.atom( "OVL1" ).xyz();
      } else if ( id.atomno() == upper_rsd.mainchain_atoms()[2] ) {
        xyz_fixed = lower_rsd.atom( "OVL2" ).xyz();
      } else if ( upper_rsd.atom_name( id.atomno() ) == "OVU1" ) {
        xyz_fixed = lower_rsd.atom( lower_rsd.mainchain_atoms()[nbb] ).xyz();
      } else {
        match = false;
      }

      if ( match ) {
        // factor = deriv / r = 1/r //sounds ugly!
        Vector const f2 ( xyz_moving - xyz_fixed );
        Real const dist ( f2.length() );
        if ( dist >= 0.01 ) { //avoid getting too close to singularity...
          Real const invdist( 1.0 / dist );
          F1 += weights[ linear_chainbreak ] * invdist * cross( xyz_moving, xyz_fixed ) / 3;
          F2 += weights[ linear_chainbreak ] * invdist * ( xyz_moving - xyz_fixed ) / 3;
        }
      }
    }
  }

  /// @brief LinearChainbreak Energy is context independent and thus indicates that no context graphs need to
  /// be maintained by class Energies
  void
  LinearChainbreakEnergy::indicate_required_context_graphs(utility::vector1<bool>&) const {} /*context_graphs_required*/
core::Size
LinearChainbreakEnergy::version() const
{
  return 1; // Initial versioning
}

} // namespace methods
} // namespace scoring
} // namespace core
