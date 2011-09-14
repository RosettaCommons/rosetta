// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/moves/HelixAssemblyMover.cc
/// @brief HelixAssemblyMover methods implemented
/// @author Tim Jacobs

//Unit headers
#include <devel/helixAssembly/HelixAssemblyMover.hh>

//Package headers
#include <devel/helixAssembly/HelixAssemblyJob.hh>

//core library
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionInfo.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/methods/EnergyMethod.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/etable/Etable.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/conformation/Residue.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArrayInitializer.hh>

//utility & numeric
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/HelixAssembly.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility Headers
#include <basic/Tracer.hh>

// C++ headers
#include <sstream>
#include <iostream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <utility>
#include <algorithm>
#include <iterator>

using namespace std;
using namespace core;
using namespace core::pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("HelixAssemblyMover");

///@brief
HelixAssemblyMover::HelixAssemblyMover() :
    single_helix_rmsd_cutoff_(1),
    helix_pair_rmsd_cutoff_(1.3),
    helix_cap_distance_cutoff_(12.0),
    helix_contact_distance_cutoff_(10.0),
    minimum_helix_contacts_(6)
{
  init_from_options();
}

HelixAssemblyMover::~HelixAssemblyMover(){}

core::scoring::ScoreFunctionOP HelixAssemblyMover::get_scorefxn() const
{
    return scorefxn_;
}

core::Size HelixAssemblyMover::get_frag1_start() const
{
    return frag1_start_;
}

core::Size HelixAssemblyMover::get_frag2_start() const
{
    return frag2_start_;
}

core::Size HelixAssemblyMover::get_frag1_end() const
{
    return frag1_end_;
}

core::Size HelixAssemblyMover::get_frag2_end() const
{
    return frag2_end_;
}

Real HelixAssemblyMover::get_helix_cap_distance_cutoff() const
{
    return helix_cap_distance_cutoff_;
}

Real HelixAssemblyMover::get_helix_contact_distance_cutoff() const
{
    return helix_contact_distance_cutoff_;
}

Real HelixAssemblyMover::get_helix_pair_rmsd_cutoff() const
{
    return helix_pair_rmsd_cutoff_;
}

core::Size HelixAssemblyMover::get_minimum_helix_contacts() const
{
    return minimum_helix_contacts_;
}

std::string HelixAssemblyMover::get_query_structure_path() const
{
    return query_structure_path_;
}

std::string HelixAssemblyMover::get_query_structure_string() const
{
    return query_structure_string_;
}

Real HelixAssemblyMover::get_single_helix_rmsd_cutoff() const
{
    return single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::set_frag1_start(core::Size frag1_start_)
{
    this->frag1_start_ = frag1_start_;
}

void HelixAssemblyMover::set_frag2_start(core::Size frag2_start_)
{
    this->frag2_start_ = frag2_start_;
}

void HelixAssemblyMover::set_frag1_end(core::Size frag1_end_)
{
    this->frag1_end_ = frag1_end_;
}

void HelixAssemblyMover::set_frag2_end(core::Size frag2_end_)
{
    this->frag2_end_ = frag2_end_;
}

void HelixAssemblyMover::set_scorefxn(core::scoring::ScoreFunctionOP scorefxn_)
{
    this->scorefxn_ = scorefxn_;
}

void HelixAssemblyMover::set_helix_cap_distance_cutoff(Real helix_cap_distance_cutoff_)
{
    this->helix_cap_distance_cutoff_ = helix_cap_distance_cutoff_;
}

void HelixAssemblyMover::set_helix_contact_distance_cutoff(Real helix_contact_distance_cutoff_)
{
    this->helix_contact_distance_cutoff_ = helix_contact_distance_cutoff_;
}

void HelixAssemblyMover::set_helix_pair_rmsd_cutoff(Real helix_pair_rmsd_cutoff_)
{
    this->helix_pair_rmsd_cutoff_ = helix_pair_rmsd_cutoff_;
}

void HelixAssemblyMover::set_minimum_helix_contacts(core::Size minimum_helix_contacts_)
{
    this->minimum_helix_contacts_ = minimum_helix_contacts_;
}

void HelixAssemblyMover::set_query_structure_path(std::string query_structure_path_)
{
    this->query_structure_path_ = query_structure_path_;
}

void HelixAssemblyMover::set_query_structure_string(std::string query_structure_string_)
{
    this->query_structure_string_ = query_structure_string_;
}

void HelixAssemblyMover::set_single_helix_rmsd_cutoff(Real single_helix_rmsd_cutoff_)
{
    this->single_helix_rmsd_cutoff_ = single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::init_from_options(){
  set_scorefxn(core::scoring::getScoreFunction());
  set_query_structure_path(basic::options::option[HelixAssembly::query_structure_path]);
  set_frag1_start(basic::options::option[HelixAssembly::frag1_start]);
  set_frag2_start(basic::options::option[HelixAssembly::frag2_start]);
  set_frag1_end(basic::options::option[HelixAssembly::frag1_end]);
  set_frag2_end(basic::options::option[HelixAssembly::frag2_end]);
  set_single_helix_rmsd_cutoff(basic::options::option[HelixAssembly::single_helix_rmsd_cutoff]);
  set_query_structure_path(basic::options::option[HelixAssembly::query_structure_path]);
  set_minimum_helix_contacts(basic::options::option[HelixAssembly::minimum_helix_contacts]);
  set_helix_pair_rmsd_cutoff(basic::options::option[HelixAssembly::helix_pair_rmsd_cutoff]);
  set_helix_contact_distance_cutoff(basic::options::option[HelixAssembly::helix_contact_distance_cutoff]);
  set_helix_cap_distance_cutoff(basic::options::option[HelixAssembly::helix_cap_distance_cutoff]);
}

Pose HelixAssemblyMover::combinePoses(Pose const & pose1, Pose const & pose2){
  Pose newPose;

  if(pose1.total_residue()>=1){
      newPose.append_residue_by_jump(pose1.residue(1), newPose.total_residue() , "", "", false/*start new chain*/);
      for(Size i=2; i<=pose1.total_residue(); i++){
          if(pose1.residue(i).is_lower_terminus() || pose1.residue(i).is_upper_terminus()){
              newPose.append_residue_by_jump(pose1.residue(i), newPose.total_residue(), "","", false);
          }
          else{
              newPose.append_residue_by_bond(pose1.residue(i));
          }
      }
  }

  if(pose2.total_residue()>=1){
      newPose.append_residue_by_jump(pose2.residue(1), newPose.total_residue() , "", "", true/*start new chain*/);
      for(Size i=2; i<=pose2.total_residue(); i++){
          if(pose2.residue(i).is_lower_terminus() || pose2.residue(i).is_upper_terminus()){
              newPose.append_residue_by_jump(pose2.residue(i), newPose.total_residue(), "","", false);
          }
          else{
              newPose.append_residue_by_bond(pose2.residue(i));
          }
      }
  }

  return newPose;
}

utility::vector1< std::pair< Size,Size > > HelixAssemblyMover::findHelices(Pose const & pose){

  utility::vector1< std::pair< Size,Size > > helix_endpts;
  for(Size i=1; i<=pose.total_residue(); i++){

      //find all the strands in the structure
      Size helix_start(0);
      Size helix_end(0);
      if(pose.secstruct(i) == 'H'){

          helix_start=i;
          while(i<=pose.total_residue() && pose.secstruct(i)=='H'){
              i++;
          }
          helix_end=i;

          helix_endpts.push_back(make_pair(helix_start, helix_end));
      }
  }
  return helix_endpts;
}

///@details return all poses from the targetPose that contain an RMSD match to the queryFragment
utility::vector1<Size> HelixAssemblyMover::findFragments(Pose const & targetPose, Pose const & queryFragment,
  utility::vector1< std::pair< Size,Size > > helix_endpts){

  utility::vector1<Size> foundFragments;
  cout << "target pose size: " << targetPose.total_residue() << endl;
  cout << "query fragment size: " << queryFragment.total_residue() << endl;

  for(Size j=1; j<=helix_endpts.size(); j++){
      if(helix_endpts[j].first + queryFragment.total_residue() - 1 <= helix_endpts[j].second){
          for(Size i=helix_endpts[j].first; i<=helix_endpts[j].second-queryFragment.total_residue()+1; i++){

              //make sure we don't make a test fragment out of two separate chains
              if(targetPose.total_residue() > i+queryFragment.total_residue() &&
                  targetPose.residue(i).chain() == targetPose.residue(i+queryFragment.total_residue()-1).chain()){

                  Pose testFragment(targetPose, i, i+queryFragment.total_residue()-1);

                  Real rmsd = core::scoring::bb_rmsd_including_O(queryFragment, testFragment);

                  if(rmsd <= single_helix_rmsd_cutoff_){
                      //cout << "found matching fragment(" << i << "): " << rmsd << endl;
                      foundFragments.push_back(i);
                  }
              }
          }
      }
  }

  return foundFragments;
}//findFragment

bool HelixAssemblyMover::checkHelixContacts(Pose const & pose, Pose const & fragment1, Pose const & fragment2,
            std::pair< Size,Size > helix_endpts){

  Pose testHelix(pose, helix_endpts.first, helix_endpts.second);
  Size distCutoff = pow((double)helix_contact_distance_cutoff_, 2);

  int interactionCounter = 0;
  for(Size i=1; i<=testHelix.total_residue(); i++){
      bool frag1Pass=false;
      bool frag2Pass=false;
      for(Size j=1; j<=fragment1.total_residue(); j++){
          if(testHelix.residue(i).atom("CA").xyz().distance_squared(fragment1.residue(j).atom("CA").xyz()) < distCutoff){
              frag1Pass=true;
              break;
          }
      }
      for(Size j=1; j<=fragment2.total_residue(); j++){
          if(testHelix.residue(i).atom("CA").xyz().distance_squared(fragment2.residue(j).atom("CA").xyz()) < distCutoff){
              frag2Pass=true;
              break;
          }
      }
      if(frag1Pass && frag2Pass){
          interactionCounter++;
      }
  }

  cout << "testHelix residues: " << testHelix.total_residue() << endl;
  cout << "interactionCounter: " << interactionCounter << endl;


  if(interactionCounter >= minimum_helix_contacts_){
      return true;
  }
  return false;
}

///@details search the pose for a helical segment that makes interactions with both of the pose fragments
utility::vector1< std::pair< Size,Size > > HelixAssemblyMover::findPartnerHelices(Pose const & pose, Pose const & fragment1,
    Pose const & fragment2, Size frag1Start, Size frag2Start, utility::vector1< std::pair< Size,Size > > helix_endpts){

  utility::vector1< std::pair< Size,Size > > partnerHelices;

  //look at each helix in the search structure as a potential helix to add
  for(Size i=1; i<=helix_endpts.size(); i++){

      //don't look at this helix if it contains either of the helices that we used in the initial fragment hits
      if(!((frag1Start >= helix_endpts[i].first && frag1Start <= helix_endpts[i].second) ||
          (frag2Start >= helix_endpts[i].first && frag2Start <= helix_endpts[i].second))){

          cout << "frag1 start: " << frag1Start << endl;
          cout << "frag2 start: " << frag2Start << endl;


          bool foundStart = false;
          bool foundEnd = false;
          Size helixStart;
          Size helixEnd;
          Size distCutoff = pow((double)helix_cap_distance_cutoff_,2);
          //Search for helical fragments that have n-terms of query helix close to n-terms frag1 & c-terms of frag2.
          for(Size helixOffset=0; helixOffset<=(helix_endpts[i].second-helix_endpts[i].first+1)/2; helixOffset++){
              if(!foundStart){
                  core::DistanceSquared startDistance1 = fragment1.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  cout << "residue " << helix_endpts[i].first+helixOffset << " start distance 1: " << startDistance1 << endl;

                  core::DistanceSquared startDistance2 = fragment2.residue(fragment2.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  cout << "residue " << helix_endpts[i].first+helixOffset << " start distance 2: " << startDistance2 << endl;

                  if(startDistance1 < distCutoff && startDistance2 < distCutoff){
                      helixStart = helix_endpts[i].first+helixOffset;
                      foundStart = true;
                      cout << "Found third helix start at residue " << helixStart << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
                  }

              }

              if(!foundEnd){
                  core::DistanceSquared endDistance1 = fragment1.residue(fragment1.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  cout << "residue " << helix_endpts[i].second-helixOffset << " end distance 1: " << endDistance1 << endl;

                  core::DistanceSquared endDistance2 = fragment2.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  cout << "residue " << helix_endpts[i].second-helixOffset << " end distance 2: " << endDistance2 << endl;

                  if(endDistance1 < distCutoff && endDistance2 < distCutoff){
                      helixEnd = helix_endpts[i].second-helixOffset;
                      foundEnd = true;
                      cout << "Found third helix end at residue " << helixEnd << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
                  }

              }
          }
          if(foundStart && foundEnd){
              //by now we've found a helical fragment that has ends close to each of the fragments, now we want to
              //check for more extensive interactions with *both* of the fragments
              std::pair<Size, Size> closeHelix(helixStart, helixEnd);

              if(checkHelixContacts(pose, fragment1, fragment2, closeHelix)){
                  partnerHelices.push_back(closeHelix);
              }
          }
          foundStart=false;
          foundEnd=false;
          //Search for helical fragments that have n-terms of query helic close to c-terms frag1 & n-terms of frag2.
          for(Size helixOffset=0; helixOffset<=(helix_endpts[i].second-helix_endpts[i].first+1)/2; helixOffset++){
              if(!foundStart){
                  core::DistanceSquared startDistance1 = fragment1.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  cout << "start distance 1: " << startDistance1 << endl;

                  core::DistanceSquared startDistance2 = fragment2.residue(fragment2.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  cout << "start distance 2: " << startDistance2 << endl;

                  if(startDistance1 < distCutoff && startDistance2 < distCutoff){
                      helixStart = helix_endpts[i].first+helixOffset;
                      foundStart = true;
                      cout << "Found third helix start at residue " << helixStart << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
                  }

              }

              if(!foundEnd){
                  core::DistanceSquared endDistance1 = fragment1.residue(fragment1.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  cout << "end distance 1: " << endDistance1 << endl;

                  core::DistanceSquared endDistance2 = fragment2.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  cout << "end distance 2: " << endDistance2 << endl;

                  if(endDistance1 < distCutoff && endDistance2 < distCutoff){
                      helixEnd = helix_endpts[i].second-helixOffset;
                      foundEnd = true;
                      cout << "Found third helix end at residue " << helixEnd << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
                  }

              }
          }
          if(foundStart && foundEnd){
              cout << "Found potential third helix from: " << helixStart << " to " << helixEnd << endl;

              //by now we've found a helical fragment that has ends close to each of the fragments, now we want to
              //check for more extensive interactions with *both* of the fragments
              std::pair<Size, Size> closeHelix(helixStart, helixEnd);

              if(checkHelixContacts(pose, fragment1, fragment2, closeHelix)){
                  partnerHelices.push_back(closeHelix);
              }
          }
      }
  }
  return partnerHelices;
}

void HelixAssemblyMover::superimposeBundles(Pose & query_structure, Pose const & results_structure){
  id::AtomID_Map< core::id::AtomID > atom_map;
  // maps every atomid to bogus atom. Atoms in the query_structure that are left mapped to the bogus atoms will not be used in the superimposition
  atom_map.clear();
  core::pose::initialize_atomid_map( atom_map, query_structure, id::BOGUS_ATOM_ID );

  //superimpose the found two-helix pose onto the query structure
  for(Size j=0; j<=get_frag1_end()-get_frag1_start(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      cout << "superimposing " << j+get_frag1_start() << " of query structure to " << j+1 << " of results structure" << endl;


      core::id::AtomID const id1( query_structure.residue(j+get_frag1_start()).atom_index("CA"), j+get_frag1_start());
      core::id::AtomID const id2( results_structure.residue(j+1).atom_index("CA"), j+1 );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( query_structure.residue(j+get_frag1_start()).atom_index("C"), j+get_frag1_start() );
      core::id::AtomID const id4( results_structure.residue(j+1).atom_index("C"), j+1 );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( query_structure.residue(j+get_frag1_start()).atom_index("N"), j+get_frag1_start() );
      core::id::AtomID const id6( results_structure.residue(j+1).atom_index("N"), j+1 );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( query_structure.residue(j+get_frag1_start()).atom_index("O"), j+get_frag1_start() );
      core::id::AtomID const id8( results_structure.residue(j+1).atom_index("O"), j+1 );
      atom_map[ id7 ] = id8;
  }

  for(Size j=0; j<=get_frag2_end()-get_frag2_start(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      core::Size results_frag2_start = get_frag1_end()-get_frag1_start()+2;

      cout << "superimposing " << j+get_frag2_start() << " of query structure to " << j+results_frag2_start << " of results structure" << endl;


      core::id::AtomID const id1( query_structure.residue(j+get_frag2_start()).atom_index("CA"), j+get_frag2_start());
      core::id::AtomID const id2( results_structure.residue(j+results_frag2_start).atom_index("CA"),j+results_frag2_start );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( query_structure.residue(j+get_frag2_start()).atom_index("C"), j+get_frag2_start() );
      core::id::AtomID const id4( results_structure.residue(j+results_frag2_start).atom_index("C"), j+results_frag2_start );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( query_structure.residue(j+get_frag2_start()).atom_index("N"), j+get_frag2_start() );
      core::id::AtomID const id6( results_structure.residue(j+results_frag2_start).atom_index("N"), j+results_frag2_start );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( query_structure.residue(j+get_frag2_start()).atom_index("O"), j+get_frag2_start() );
      core::id::AtomID const id8( results_structure.residue(j+results_frag2_start).atom_index("O"), j+results_frag2_start );
      atom_map[ id7 ] = id8;
    }

  //do superimpose
  scoring::superimpose_pose(query_structure, results_structure, atom_map);
}

core::Real HelixAssemblyMover::bb_score(pose::Pose & pose, core::Size unique_chain_num, core::scoring::ScoreFunctionOP & scorefxn){

  // score the bb-bb energy between the given chain num and the rest of the pose.
  // Adopted by code written by P.Doug Renfrew

  utility::vector1<core::conformation::Atom> unique_chain_bb_atoms;
  utility::vector1<core::conformation::Atom> other_bb_atoms;

  for( Size j = 1; j <= pose.total_residue(); ++j ) {
      core::conformation::Residue const & res( pose.residue(j) );
      core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
      //assert( bb_ai.size() == 4 );
      core::Size chain_num( res.chain() );
      for( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
          if( chain_num == unique_chain_num )
            unique_chain_bb_atoms.push_back( res.atom(jj) );

          else
            other_bb_atoms.push_back( res.atom(jj) );
      }
  }

  //NOW SCORE!
  // get instance of etable energy method
  core::scoring::methods::EnergyMethodOptions const & emo(scorefxn->energy_method_options());
  core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo.etable_type())));
  core::scoring::etable::EtableEnergy ete( et, emo );

  // iterate over both sets of atom and add into one emapvector
  //core::scoring::TwoBodyEMapVector tbemv;
  core::scoring::EMapVector tbemv;
  core::Real atr_wt( (*scorefxn).get_weight(core::scoring::fa_atr) );
  core::Real rep_wt( (*scorefxn).get_weight(core::scoring::fa_rep) );
  for ( Size ii = 1; ii <= unique_chain_bb_atoms.size(); ++ii ){
      for ( Size jj = 1; jj <= other_bb_atoms.size(); ++jj ) {
          //calc distr squared
          Real d2(unique_chain_bb_atoms[ii].xyz().distance_squared(other_bb_atoms[jj].xyz()));
          ete.atom_pair_energy( unique_chain_bb_atoms[ii], other_bb_atoms[jj], 1, tbemv, d2 );
      }
  }
  core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );

  cout<< "Backbone-backbone score: " << bb_energy << std::endl;

  return bb_energy;
}//end bb_score

///@details
utility::vector1<HelixAssemblyJob> HelixAssemblyMover::apply( HelixAssemblyJob & job, int rank ){

  try{
    utility::vector1<HelixAssemblyJob> new_jobs;

    cout << "working on file: " << job.get_job_name() << endl;

    Pose query_structure;
    core::import_pose::pose_from_pdbstring(query_structure, job.get_query_structure());
    cout << "Full Query Structure is " << query_structure.total_residue() << " residues" << endl;

    Pose search_structure;
    core::import_pose::pose_from_pdbstring(search_structure, job.get_search_structure());
    cout << "Search Structure is " << search_structure.total_residue() << " residues" << endl;

    //don't do anything with this round if the search structure is empty. Not testing for this causes all sorts of bad behavior
    if(search_structure.total_residue() <= 0){return new_jobs;}

    this->set_frag1_start(job.get_frag1_start());
    this->set_frag1_end(job.get_frag1_end());
    this->set_frag2_start(job.get_frag2_start());
    this->set_frag2_end(job.get_frag2_end());

    cout << "Frag1 (" << get_frag1_start() << ":" << get_frag1_end() << endl;
    Pose helixFragment1(query_structure, get_frag1_start(), get_frag1_end());
    cout << "Frag2 (" << get_frag2_start() << ":" << get_frag2_end() << endl;
    Pose helixFragment2(query_structure, get_frag2_start(), get_frag2_end());
    cout << "Fragment 1 is " << helixFragment1.total_residue() << " residues" << endl;
    cout << "Fragment 2 is " << helixFragment2.total_residue() << " residues" << endl;



    Pose combined_query_fragments = combinePoses(helixFragment1, helixFragment2);

    //set up dssp info - necessary in order to find helices based on secondary structure
    core::scoring::dssp::Dssp dssp( search_structure );
    dssp.insert_ss_into_pose( search_structure );

    utility::vector1< std::pair< Size,Size > > helix_endpts;
    helix_endpts = findHelices(search_structure);
    cout << "Found " << helix_endpts.size() << " helices in search structure" << endl;


    //search helix poses for all close RMSD matches to each query fragment
    utility::vector1<Size> fragment1Results = findFragments(search_structure, helixFragment1, helix_endpts);
    utility::vector1<Size> fragment2Results = findFragments(search_structure, helixFragment2, helix_endpts);

    cout << "Found " << fragment1Results.size() << " fragments for frag1 & " << fragment2Results.size() << " for frag2." << endl;


    //keep track of number of hits in structure, for output filename purposes
    Size resultsCounter;
    resultsCounter = 1;

    Size threeHelixCounter;
    threeHelixCounter = 1;

    //loop through lists of fragments for each helix and check rmsd to the full query structure
    for(Size i=1; i<=fragment1Results.size(); i++){
        for(Size j=1; j<=fragment2Results.size(); j++){
            Pose fragment1(search_structure, fragment1Results[i], fragment1Results[i]+helixFragment1.total_residue()-1);
            Pose fragment2(search_structure, fragment2Results[j], fragment2Results[j]+helixFragment2.total_residue()-1);
            Pose combinedResultFragments(combinePoses(fragment1, fragment2));

            Real bbrmsd = core::scoring::bb_rmsd_including_O(combined_query_fragments, combinedResultFragments);

            if(bbrmsd <= helix_pair_rmsd_cutoff_){

                cout << "combined bb_atom RMSD(" << i << ", " << j << "): " << bbrmsd << endl;


                utility::vector1< std::pair< Size,Size > > helix_partners = findPartnerHelices(search_structure, fragment1, fragment2,
                    fragment1Results[i], fragment2Results[j], helix_endpts);

                cout << "found " << helix_partners.size() << " suitable helix partners" << endl;


                for(Size k=1; k<=helix_partners.size(); k++){

                    Pose thirdHelix(search_structure, helix_partners[k].first, helix_partners[k].second);

                    cout << "New helical pose created." << endl;


                    //superimpose the found query helix pair onto the found pair
                    superimposeBundles(query_structure, combinedResultFragments);

                    cout << "Poses have been superimposed for bundle creation" << endl;


                    //combine query structure with third helix
                    Size third_helix_start(query_structure.total_residue()+1);
                    Pose new_bundle = combinePoses(query_structure, thirdHelix);
                    Size new_chain = new_bundle.chain(new_bundle.total_residue());
                    Size third_helix_end(new_bundle.total_residue());

                    cout << "New bundle created" << endl;


                    //check for backbone clashes introduced by adding the new helix
                    Real clash_score = bb_score(new_bundle, new_chain, scorefxn_);

                    //TODO (tjacobs) Turn this into an option
                    if(clash_score <= 5){
                        stringstream tempStream;
                        new_bundle.dump_pdb(tempStream, "");
                        std::string new_bundle_string = tempStream.str();

                        HelixAssemblyJob new_job1;
                        new_job1.set_query_structure(new_bundle_string);
                        new_job1.set_round(job.get_round()+1);

                        new_job1.set_frag1_start(third_helix_start);
                        new_job1.set_frag1_end(third_helix_end);

                        new_job1.set_frag2_start(job.get_frag1_start());
                        new_job1.set_frag2_end(job.get_frag1_end());

                        new_job1.set_job_name(job.get_job_name());

                        new_jobs.push_back(new_job1);

                        HelixAssemblyJob new_job2;
                        new_job2.set_query_structure(new_bundle_string);
                        new_job2.set_round(job.get_round()+1);

                        new_job2.set_frag1_start(third_helix_start);
                        new_job2.set_frag1_end(third_helix_end);

                        new_job2.set_frag2_start(job.get_frag2_start());
                        new_job2.set_frag2_end(job.get_frag2_end());

                        new_job2.set_job_name(job.get_job_name());

                        new_jobs.push_back(new_job2);
                    }
                }
                resultsCounter++;
            }
        }
    }
    return new_jobs;
  }
  catch(...){
      cout << "ERROR: Exception thrown in HelixAssemblyMover" << endl;
      utility::vector1<HelixAssemblyJob> blankJobs;
      return blankJobs;
  }
}//apply
