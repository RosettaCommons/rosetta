// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @author Oliver Lange


// Unit Headers
#include <protocols/noesy_assign/CovalentCompliance.hh>

// Utility headers
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

// Singleton instance and mutex static data members
namespace utility {

using protocols::noesy_assign::CovalentCompliance;

#if defined MULTI_THREADED && defined CXX11
template <> std::mutex utility::SingletonBase< CovalentCompliance > ::singleton_mutex_;
template <> std::atomic< CovalentCompliance * > utility::SingletonBase< CovalentCompliance >::instance_( 0 );
#else
template <> CovalentCompliance * utility::SingletonBase< CovalentCompliance >::instance_( 0 );
#endif

}

namespace protocols {
namespace noesy_assign {

bool fall_back( core::id::NamedAtomID const& _atom1, core::id::NamedAtomID const& _atom2 ) {
  bool flip = _atom1.rsd() > _atom2.rsd();
  core::id::NamedAtomID const& atom1( flip ? _atom2 : _atom1 );
  core::id::NamedAtomID const& atom2( flip ? _atom1 : _atom2 );
  //, core::pose::Pose const& pose, Real dmax ) {

  //compute upper distances for atoms less than 2 torsions apart.
  //if residues are not sequential this we will not happen
  if ( atom2.rsd() - atom1.rsd() > 1 ) return false;

  if ( atom1.rsd() != atom2.rsd() ) { //sequential --- only Ha - H
    if ( atom1.atom() == "HA" && atom2.atom() == "H" ) {
      return true;
    }
    //no more sequential cases
    return false;
  }

	///change Oct 12 2010: intra-residue is always true
	return true;

  //intra-residue
  if ( atom1.atom() == "H" || atom2.atom() == "H" ) {
    if ( atom1.atom() == "HA" || atom2.atom() == "HA" ) { // H - HA or HA-H
      return true;
    }
    if ( (atom1.atom().find("B") != std::string::npos) || (atom2.atom().find("B") !=std::string::npos )) {
			return true;
    }
    return false; // H- HX (HX not HA,HB ) that will be more than two torsion angles
  }

  if ( atom1.atom() == "HA" || atom2.atom() == "HA" ) {
    return true; // for now assume any HA-HX connection is close enough...
  }

	//both protons are not on backbone (call this good enough... ).
	if ( atom1.atom() != "H" && atom1.atom().find("A") == std::string::npos && atom2.atom() != "H" && atom2.atom().find("A") == std::string::npos ) {
		return true;
	}

	//two torsions between HA and HG
	if ( atom1.atom().find("A") != std::string::npos && atom2.atom().find("G") != std::string::npos ) return true;
	if ( atom2.atom().find("A") != std::string::npos && atom1.atom().find("G") != std::string::npos ) return true;



  return false;
  //now what to do with QB, QG, etc. treat as one of the H with same name...
  // heavy atoms as main dimensions ( some of paolos data sets? )

  //translate QB, QG, QD into HA, HB, HD, etc
  //	AmbiguousNMRDistanceConstraint translator( atom1, atom2, pose, NULL );
  //	AtomID num_atom1( translator.atom( 1 ) );
  //	AtomID num_atom2( translator.atom( translator.natoms( 1 )+1 ) );

  //how to find torsion angles between two atoms?
  //okay we are definitely intra-residue... (if we just look at protons... )
  //i guess one can take highest atom number and work up the atom-tree until number is lower or equal than the target.
  // if it is to low go down some other branch...

  //this is highly complicated.. and will not yield to much benefit. the local structure of decoys
  // should be pretty good and we will see those NOEs recommended from there.
}

CovalentCompliance *
CovalentCompliance::create_singleton_instance()
{
	return new CovalentCompliance;
}

CovalentCompliance::CovalentCompliance() :
  covalent_distances_( /* NULL */ )
{}

/// @details WARNING WARNING WARNING! THREAD UNSAFE!
void CovalentCompliance::load_dist_table( std::string const& file ) {
  covalent_distances_ = FragsToAtomDistOP( new FragsToAtomDist( file ) );
}

bool CovalentCompliance::is_compliant( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2, core::Real cutoff ) const {
  if ( covalent_distances_ ) {
    return covalent_distances_->distance( atom1, atom2 ) < cutoff;
  } else {
    return fall_back( atom1, atom2 );
  }
}

core::Real CovalentCompliance::distance( core::id::NamedAtomID const& atom1, core::id::NamedAtomID const& atom2 ) const {
  if ( covalent_distances_ ) {
    return covalent_distances_->distance( atom1, atom2 );
  } else {
    if ( fall_back( atom1, atom2 ) ) {
      return 3.0;
    } else {
      return 100.0;
    }
  }
}


}
}
