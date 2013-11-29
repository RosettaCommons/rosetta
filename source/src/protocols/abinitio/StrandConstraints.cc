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
#include <protocols/abinitio/StrandConstraints.hh>

// Package Headers
#include <protocols/abinitio/PairingStatistics.hh>
// AUTO-REMOVED #include <protocols/abinitio/Template.hh>
// AUTO-REMOVED #include <protocols/abinitio/Templates.hh>
// AUTO-REMOVED #include <protocols/abinitio/TemplateJumpSetup.hh>

// Project Headers
#include <core/types.hh>

// AUTO-REMOVED #include <core/pose/Pose.hh>
#include <core/pose/util.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// AUTO-REMOVED
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/fragment/FragSet.fwd.hh>
// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/FragData.hh> //to get secondary structure
// AUTO-REMOVED #include <core/fragment/SecstructSRFD.hh> //to get secondary structure
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.hh>

#include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>

// AUTO-REMOVED #include <core/sequence/util.hh>
// AUTO-REMOVED #include <core/id/SequenceMapping.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <core/fragment/SecondaryStructure.hh>
// AUTO-REMOVED #include <protocols/jumping/JumpSample.hh>

#ifdef WIN32
#include <core/scoring/dssp/PairingsList.hh>
#endif



//numeric headers
// AUTO-REMOVED #include <numeric/random/random.hh>

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <utility/io/izstream.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/vector1.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <cstdlib>
#include <string>
#include <vector>

#include <utility/vector1.hh>



static basic::Tracer tr("protocols.abinitio.StrandConstraints");
using namespace core;
using namespace basic;


namespace protocols {
namespace abinitio {

/// @details Auto-generated virtual destructor
StrandConstraints::~StrandConstraints() {}

//using namespace jumping;

Size register_cutoff( 5 );
Size residue_cutoff( 5 );

//static numeric::random::RandomGenerator RG(1123123544);  // <- Magic number, do not change it!
typedef utility::vector1< Size > SizeList;
bool AlternativePairings::compatible( core::scoring::dssp::StrandPairing const& strand_pairing ) const {
  if ( pairings_.size() == 0 ) return true;
  if ( strand_pairing.antiparallel() != antiparallel() ) return false;
  Size const new_reg ( strand_pairing.get_register() );

  //get list of pairs for all pairings that roughly match the register of the new strand_pairing
  for ( Pairings::const_iterator it=pairings_.begin(), eit = pairings_.end();
	it!=eit; ++it ) {
    //check register
    Size const reg ( it->pairing().get_register() );
    if ( std::abs( (int)new_reg  - (int)reg ) > (int)register_cutoff ) return false;
    //good register... get residues
    if ( (int) strand_pairing.end1() <  ( (int) it->pairing().begin1() - (int) residue_cutoff ) ) return false; //really no overlap
    if ( strand_pairing.begin1() > it->pairing().end1() + residue_cutoff ) return false;
  }
  // should one also check the begin2, end2 residues?

  /*at this point:
    register match up to register_cutoff_
    if none of the staring residues overlap the gap between is at least smaller than residue_cutoff_
  */
  return true;
}

bool AlternativePairings::add_pairing( PairingStatEntry const& pairing_entry ) {
	for ( Pairings::const_iterator it = pairings_.begin(), eit = pairings_.end();
				it != eit; ++it ) {
		it->has_model( pairing_entry.models()[ 1 ] );
	}
	if ( compatible( pairing_entry.pairing() ) ) {
		anti_ = pairing_entry.pairing().antiparallel();
		pairings_.push_back( pairing_entry );
		//evtl have to make sure that I don't add "mergeable" pairings.
		return true;
	};
	return false;
}

void AlternativePairings::show( std::ostream& out ) const {
  out << "\n\n Set of alternative strand pairings: \n";
  for ( Pairings::const_iterator it = pairings_.begin(), eit = pairings_.end();
	it!=eit; ++it ) {
    out << it->weight() << " " << it->pairing() << " ";
    for ( PairingStatEntry::ModelList::const_iterator mit = it->models().begin(), emit = it->models().end();
	  mit != emit; ++mit ) {
      out << *mit << " ";
    }
    out << "\n";
  }
}


void AlternativePairings::build_constraints( pose::Pose const& pose, scoring::constraints::ConstraintCOPs& pairing_csts ) const {
  using namespace core::scoring::constraints;
  using namespace core::id;

  //we only constrain the most "visited" residues of a strand ( as seen from the start -- for now)
  // example     residue     54 55 56 57 58 59 60 61 62
  //                 freq     1  2  4  4  4  4  3 2  1
  //   from strands going: 54-60, 55-59, 56-62, 56-61
  //   might be a problem, that I don't consider where the beta-pairing goes to...

  //        constraints on residues 56 57 58 59

  //  thought about two ways to make the constraints:
  //  1) CA-CA  with lengths 4.2-4.7 for A O
  //                       5.0-5.7 for A I
  //                       4.5-5.1 for P
//   HarmonicFunc const CAfuncAO( 4.45, 0.5 ); //pleating I, 2
//   HarmonicFunc const CAfuncAI 5.35, 0.7 ); //pleating O, 1
//   HarmonicFunc const CAfuncP( 4.8, 0.6 );  //all pleatings

  BoundFunc const CAfuncAO( 4.2, 4.8, 0.2, "STRAND_AO" ); //anti - pleating O, 1
  BoundFunc const CAfuncAI( 5.0, 5.7, 0.3, "STRAND_AI" ); //anti - pleating I, 2
  BoundFunc const CAfuncP( 4.5, 5.1, 0.3, "STRAND_P" );  //para - all pleatings

  //  2) N-O(=C)  constraints of 2.8-2.95A
  //        for A I residue n constrain N-O of n+1 and O-N n-1
  //        for A O residue n constrain N-O of n and O-N of n
  //        for P  ...
  HarmonicFunc const NOfunc( 2.875, 0.2 ); //pretty much the same for all different pleatings

  //   Although 2) will give tighter constraints it might have the disadvantage that it creates high-energy barriers between
  //    different registers: going from reg N -> N+1 requires change of pleating for one residue.
  //   the N-O group has to come around fro//      m the other side of the strand... that might be bad.
  //find frequency of start residues
  std::map< Size, Real > freq; //it will count the score values
  Real max_freq( 0 );
  for ( Pairings::const_iterator it = pairings_.begin(), eit = pairings_.end();
	it != eit; ++it ) {
    Size const start( it->pairing().begin1() );
    Size const end( it->pairing().end1() );
    for ( Size pos = start; pos <= end; pos++ ) {
      if ( ! it->pairing().is_bulge( pos ) ) freq[ pos ] += it->weight();
      if ( max_freq < freq[ pos ] ) max_freq = freq[ pos ];
    }
  }

  //*** no constraints if we don't have any confidence ***
  if ( max_freq <= 10 ) return;


  //for each start position
  for ( std::map< Size, Real>::const_iterator it = freq.begin(), eit = freq.end();
	it != eit; ++it ) {

    // is frequency == max_freq ?
    if ( it->second >= max_freq * 0.8 ) {
      //yes: a pairing that we want to enforce

      // start-atom
      Size pos( it->first );
      NamedAtomID atom1( "CA", pos );

      //make ambigous constraint with one distance constraint for each pairing
      ConstraintCOPs constraints;
      std::map< Size, bool > done;

      //for all end-atoms:
      for ( Pairings::const_iterator pit = pairings_.begin(), epit = pairings_.end();
	    pit != epit; ++pit ) {
	Size pair( pit->pairing().get_pair( pos ) );
	if ( pair ) { //we have a pairing make a constraint for this guy
	  if ( !done[ pair ] ) {
	    done[ pair ] = true;
	    NamedAtomID atom2( "CA", pair );

	    //get distance-function parameters right
	    FuncOP myfunc;
	    if ( pit->pairing().antiparallel() ) {
	      if ( pit->pairing().get_pleating( pos ) == 2 ) { // Inward
					myfunc = CAfuncAI.clone();
	      } else { // Outward
					myfunc = CAfuncAO.clone();
	      } //pleating
	    } else { //not anti-parallel
	      myfunc = CAfuncP.clone();
	    }

	    //make the constraint
	    tr.Debug << " add constraint for " << pos << "->" << pair << " " << pit->pairing().get_pleating( pos ) << std::endl;
	    constraints.push_back( new AtomPairConstraint(
			core::pose::named_atom_id_to_atom_id( atom1, pose),
			core::pose::named_atom_id_to_atom_id( atom2, pose),
			myfunc ) );

	  }
	}// if ( pair )
      } //for Pairings
      //      AmbiguousConstraint test_cst( constraints );
      //      test_cst.show_def( tr.Info, pose );

      pairing_csts.push_back( new AmbiguousConstraint( constraints ) );
    } // if freq == max_freq
  } // for each start position
}

std::ostream& operator<< ( std::ostream& out, AlternativePairings const& alt_pairs ) {
  alt_pairs.show( out );
  return out;
}

void StrandConstraints::add_pairing( core::scoring::dssp::StrandPairing const& pairing, std::string model ) {
  add_pairing( PairingStatEntry( pairing, model ) );
}

void StrandConstraints::add_pairing( PairingStatEntry const& pairing_entry ) {
  bool success( false );
  tr.Debug << "add pairing to FuzzyTopology "<< pairing_entry << std::endl;
  for ( FuzzyTopology::iterator it = fuzzy_topology_.begin(), eit = fuzzy_topology_.end();
	it != eit && !success; ++it ) {
    success = it->add_pairing( pairing_entry );
  }
  if ( !success ) {
    AlternativePairings new_strand;
    new_strand.add_pairing( pairing_entry );
    fuzzy_topology_.push_back( new_strand );
  }
}


StrandConstraints::StrandConstraints( PairingStatistics const& strand_stats ) {
  for ( PairingStatistics::const_iterator it = strand_stats.begin(), eit = strand_stats.end();
	it != eit; ++it ) {
    add_pairing( it->second );
  }
  tr.Info << (*this) << std::endl;
}

void StrandConstraints::show( std::ostream& out ) const {
  for (FuzzyTopology::const_iterator it = fuzzy_topology_.begin(), eit = fuzzy_topology_.end();
       it != eit; ++it ) {
    it->show( out );
  }
}

scoring::constraints::ConstraintCOPs StrandConstraints::build_constraints( pose::Pose const& pose ) const {
  scoring::constraints::ConstraintCOPs all_constraints;
  for (FuzzyTopology::const_iterator it = fuzzy_topology_.begin(), eit = fuzzy_topology_.end();
       it != eit; ++it ) {
    it->build_constraints( pose, all_constraints );
  }

//   for ( scoring::constraints::ConstraintCOPs::iterator it = all_constraints.begin(), eit = all_constraints.end();
// 	it != eit; ++it ) {
//     (*it)->show_def( tr.Info, pose );
//   }

  return all_constraints;
}

std::ostream& operator<< ( std::ostream& out, StrandConstraints const& st ) {
  st.show( out );
  return out;
}


}
}
