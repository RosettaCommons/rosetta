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

HelicalFragment HelixAssemblyMover::get_query_frag_1() const
{
    return query_frag_1_;
}

HelicalFragment HelixAssemblyMover::get_query_frag_2() const
{
    return query_frag_2_;
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

void HelixAssemblyMover::set_query_frag_1(HelicalFragment frag_1_)
{
    this->query_frag_1_ = frag_1_;
}

void HelixAssemblyMover::set_query_frag_2(HelicalFragment frag_2_)
{
    this->query_frag_2_ = frag_2_;
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

  set_single_helix_rmsd_cutoff(basic::options::option[HelixAssembly::single_helix_rmsd_cutoff]);
  set_query_structure_path(basic::options::option[HelixAssembly::query_structure_path]);
  set_minimum_helix_contacts(basic::options::option[HelixAssembly::minimum_helix_contacts]);
  set_helix_pair_rmsd_cutoff(basic::options::option[HelixAssembly::helix_pair_rmsd_cutoff]);
  set_helix_contact_distance_cutoff(basic::options::option[HelixAssembly::helix_contact_distance_cutoff]);
  set_helix_cap_distance_cutoff(basic::options::option[HelixAssembly::helix_cap_distance_cutoff]);

  HelicalFragment frag1(basic::options::option[HelixAssembly::frag1_start],
      basic::options::option[HelixAssembly::frag2_end]);

  HelicalFragment frag2(basic::options::option[HelixAssembly::frag2_start],
        basic::options::option[HelixAssembly::frag2_end]);

  this->set_query_frag_1(frag1);
  this->set_query_frag_2(frag2);
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
utility::vector1<HelicalFragment> HelixAssemblyMover::findFragmentMatches(Pose const & search_structure,
    Pose const & query_structure, HelicalFragment query_fragment, utility::vector1< std::pair< Size,Size > > helix_endpts){

  utility::vector1<HelicalFragment> frag_matches;

  for(Size j=1; j<=helix_endpts.size(); j++){

      for(Size i=helix_endpts[j].first; i<=helix_endpts[j].second-query_fragment.get_size()+1; i++){

          //make sure we don't make a test fragment out of two separate chains
          if(search_structure.total_residue() > i+query_fragment.get_size() &&
              search_structure.residue(i).chain() == search_structure.residue(i+query_fragment.get_size()-1).chain()){

              HelicalFragment test_fragment(i, i+query_fragment.get_size()-1);

              std::map<core::id::AtomID, core::id::AtomID> frag_map = getFragmentMap(search_structure,
                  query_structure, test_fragment, query_fragment);

              Real rmsd = core::scoring::rms_at_all_corresponding_atoms(search_structure, query_structure, frag_map);

              if(rmsd <= single_helix_rmsd_cutoff_){
                  frag_matches.push_back(test_fragment);
              }
          }
      }
  }

  return frag_matches;
}//findFragment

bool HelixAssemblyMover::checkHelixContacts(Pose const & query_structure, Pose const & fragment1, Pose const & fragment2,
            HelicalFragment helix_to_check){

  Size distCutoff = pow((double)helix_contact_distance_cutoff_, 2);

  int interactionCounter = 0;
  for(Size i=helix_to_check.get_start(); i<=helix_to_check.get_end(); i++){
      bool frag1Pass=false;
      bool frag2Pass=false;
      for(Size j=1; j<=fragment1.total_residue(); j++){
          if(query_structure.residue(i).atom("CA").xyz().distance_squared(fragment1.residue(j).atom("CA").xyz()) < distCutoff){
              frag1Pass=true;
              break;
          }
      }
      for(Size j=1; j<=fragment2.total_residue(); j++){
          if(query_structure.residue(i).atom("CA").xyz().distance_squared(fragment2.residue(j).atom("CA").xyz()) < distCutoff){
              frag2Pass=true;
              break;
          }
      }
      if(frag1Pass && frag2Pass){
          interactionCounter++;
      }
  }

  if(interactionCounter >= minimum_helix_contacts_){
      return true;
  }
  return false;
}

///@details search the pose for a helical segment that makes interactions with both of the pose fragments
utility::vector1<HelicalFragment> HelixAssemblyMover::findPartnerHelices(Pose const & search_structure,
    Pose const & fragment1, Pose const & fragment2, std::pair<HelicalFragment, HelicalFragment> helix_pair,
    utility::vector1< std::pair< Size,Size > > helix_endpts, bool first_round, bool direction_needed){

  utility::vector1<HelicalFragment> partner_helices;

  //look at each helix in the search structure as a potential helix to add
  for(Size i=1; i<=helix_endpts.size(); i++){

      //don't look at this helix if it contains either of the helices that we used in the initial fragment hits
      if(!((helix_pair.first.get_start() >= helix_endpts[i].first &&
          helix_pair.first.get_start() <= helix_endpts[i].second) ||
          (helix_pair.second.get_start() >= helix_endpts[i].first &&
          helix_pair.second.get_start() <= helix_endpts[i].second))){

          TR << "Checking helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;

          bool foundStart = false;
          bool foundEnd = false;
          Size helixStart;
          Size helixEnd;
          Real distCutoff = pow(helix_cap_distance_cutoff_,2);
          Real minStartDistance=distCutoff*2;
          Real minEndDistance=distCutoff*2;

          for(Size helixOffset=0; helixOffset<(helix_endpts[i].second-helix_endpts[i].first+1); helixOffset++){

              //Check distance between n-term of query fragment 1 and given residue of search fragment
              core::DistanceSquared startDistance1 = fragment1.residue(1).atom("CA").xyz().distance_squared(
                  search_structure.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());

              //Check distance between c-term of query fragment 2 and given residue of search fragment
              core::DistanceSquared startDistance2 = fragment2.residue(fragment2.total_residue()).atom("CA").xyz().distance_squared(
                  search_structure.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());

              if(startDistance1 < distCutoff && startDistance2 < distCutoff && (startDistance1 + startDistance2) < minStartDistance){
                  helixStart = helix_endpts[i].first+helixOffset;
                  minStartDistance = startDistance1 + startDistance2;
                  foundStart = true;
                  TR << "Found third helix start at residue " << helixStart << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
              }

              core::DistanceSquared endDistance1 = fragment1.residue(fragment1.total_residue()).atom("CA").xyz().distance_squared(
                  search_structure.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());

              core::DistanceSquared endDistance2 = fragment2.residue(1).atom("CA").xyz().distance_squared(
                  search_structure.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());

              //Check to make sure that this point in the potential 3rd helix is close to the "end" of both query helices
              if(endDistance1 < distCutoff && endDistance2 < distCutoff && (endDistance1 + endDistance2) < minEndDistance){
                  helixEnd = helix_endpts[i].second-helixOffset;
                  minEndDistance = endDistance1 + endDistance2;
                  foundEnd = true;
                  TR << "Found third helix end at residue " << helixEnd << " in helix (" << helix_endpts[i].first << ":" << helix_endpts[i].second << ")" << endl;
              }
          }
          if(foundStart && foundEnd){
              //if the round is 1 we can have a helix in either direction (ie, we can build off either the n or c-terminus).
              //if the round is not 1 then we must alternate directions
              HelicalFragment closeHelix;
              if(helixStart < helixEnd){
                  closeHelix.set_start(helixStart);
                  closeHelix.set_end(helixEnd);
                  closeHelix.set_direction(get_query_frag_1().get_direction());
              }
              else if(helixEnd < helixStart){
                  closeHelix.set_start(helixEnd);
                  closeHelix.set_end(helixStart);
                  closeHelix.set_direction(get_query_frag_2().get_direction());
              }
              if((first_round || closeHelix.get_direction() == direction_needed) &&
                  checkHelixContacts(search_structure, fragment1, fragment2, closeHelix)){
                  partner_helices.push_back(closeHelix);
              }
          }
      }
  }
  return partner_helices;
}

void HelixAssemblyMover::superimposeBundles(Pose & query_structure, Pose const & results_structure){
  id::AtomID_Map< core::id::AtomID > atom_map;
  // maps every atomid to bogus atom. Atoms in the query_structure that are left mapped to the bogus atoms will not be used in the superimposition
  atom_map.clear();
  core::pose::initialize_atomid_map( atom_map, query_structure, id::BOGUS_ATOM_ID );

  //superimpose the found two-helix pose onto the query structure
  for(Size j=0; j< get_query_frag_1().get_size(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      Size frag1_start = get_query_frag_1().get_start();

      core::id::AtomID const id1( query_structure.residue(j+frag1_start).atom_index("CA"), j+frag1_start);
      core::id::AtomID const id2( results_structure.residue(j+1).atom_index("CA"), j+1 );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( query_structure.residue(j+frag1_start).atom_index("C"), j+frag1_start);
      core::id::AtomID const id4( results_structure.residue(j+1).atom_index("C"), j+1 );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( query_structure.residue(j+frag1_start).atom_index("N"), j+frag1_start);
      core::id::AtomID const id6( results_structure.residue(j+1).atom_index("N"), j+1 );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( query_structure.residue(j+frag1_start).atom_index("O"), j+frag1_start);
      core::id::AtomID const id8( results_structure.residue(j+1).atom_index("O"), j+1 );
      atom_map[ id7 ] = id8;
  }

  for(Size j=0; j< get_query_frag_2().get_size(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      core::Size frag2_start = get_query_frag_2().get_start();
      core::Size results_frag2_start = get_query_frag_1().get_size()+1;

      TR << "superimposing " << j+frag2_start << " of query structure to " << j+results_frag2_start << " of results structure" << endl;


      core::id::AtomID const id1( query_structure.residue(j+frag2_start).atom_index("CA"), j+frag2_start);
      core::id::AtomID const id2( results_structure.residue(j+results_frag2_start).atom_index("CA"),j+results_frag2_start );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( query_structure.residue(j+frag2_start).atom_index("C"), j+frag2_start );
      core::id::AtomID const id4( results_structure.residue(j+results_frag2_start).atom_index("C"), j+results_frag2_start );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( query_structure.residue(j+frag2_start).atom_index("N"), j+frag2_start );
      core::id::AtomID const id6( results_structure.residue(j+results_frag2_start).atom_index("N"), j+results_frag2_start );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( query_structure.residue(j+frag2_start).atom_index("O"), j+frag2_start );
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

  TR<< "Backbone-backbone score: " << bb_energy << std::endl;

  return bb_energy;
}//end bb_score

bool HelixAssemblyMover::closenessCheck(const core::Distance maxRange, const core::Distance end1Dist, const core::Distance end2Dist,
    const core::pose::Pose & search_structure, HelicalFragment search_frag_1, HelicalFragment search_frag_2){

  //Calculate the maximum that a point can be off for the rmsd to still meet the threshold, then check the ends of the helices. This
  //is a quick and dirty way to quickly prune the number of helical pairs to do a full RMSD check on.

//***Turns out this doesn't do anything for us***//
//  for(core::Size i=0; i<frag1Size; i++){
//      for(core::Size j=0; j<frag2Size; j++){
//
//          core::Distance resDistance = search_structure.residue(frag1Start+i).atom("CA").xyz().distance_squared(
//              search_structure.residue(frag2Start+j).atom("CA").xyz());
//
//          if(resDistance < 0.5){
//              TR << "FOUND OVERLAPPING!" << endl;
//              return false;
//          }
//      }
//  }

  //if the query fragments are parallel treat the search fragments the same way & vice-versa
  bool parallel = get_query_frag_1().get_direction()==get_query_frag_2().get_direction();
  if(parallel){ //parallel
      core::Distance newEnd1Dist = search_structure.residue(search_frag_1.get_start()).atom("CA").xyz().distance(
            search_structure.residue(search_frag_2.get_start()).atom("CA").xyz());

        if(abs(newEnd1Dist - end1Dist) > maxRange){return false;}

        core::Distance newEnd2Dist = search_structure.residue(search_frag_1.get_end()).atom("CA").xyz().distance(
            search_structure.residue(search_frag_2.get_end()).atom("CA").xyz());

        if(abs(newEnd2Dist - end2Dist)+abs(newEnd1Dist - end1Dist) > maxRange){return false;}
  }
  else{//anti-parallel
      core::Distance newEnd1Dist = search_structure.residue(search_frag_1.get_start()).atom("CA").xyz().distance(
          search_structure.residue(search_frag_2.get_end()).atom("CA").xyz());

      if(abs(newEnd1Dist - end1Dist) > maxRange){return false;}

      core::Distance newEnd2Dist = search_structure.residue(search_frag_1.get_end()).atom("CA").xyz().distance(
          search_structure.residue(search_frag_2.get_start()).atom("CA").xyz());

      if(abs(newEnd2Dist - end2Dist)+abs(newEnd1Dist - end1Dist) > maxRange){return false;}
  }


  return true;
}//end closenessCheck

std::map<core::id::AtomID, core::id::AtomID> HelixAssemblyMover::getFragmentMap(const core::pose::Pose & pose_1,
    const core::pose::Pose & pose_2, HelicalFragment pose_1_fragment, HelicalFragment pose_2_fragment){

  //fragments must be same size
  assert(pose_1_fragment.get_size() == pose_2_fragment.get_size());

  std::map< core::id::AtomID, core::id::AtomID > atom_map;
  atom_map.clear();
  for(Size offset=0; offset<pose_1_fragment.get_size(); ++offset){

      core::id::AtomID const id1( pose_1.residue(pose_1_fragment.get_start()+offset).atom_index("CA"), pose_1_fragment.get_start()+offset);
      core::id::AtomID const id2( pose_2.residue(pose_2_fragment.get_start()+offset).atom_index("CA"), pose_2_fragment.get_start()+offset);
      atom_map.insert(std::pair<id::AtomID, id::AtomID>(id1, id2));

      core::id::AtomID const id3( pose_1.residue(pose_1_fragment.get_start()+offset).atom_index("C"), pose_1_fragment.get_start()+offset);
      core::id::AtomID const id4( pose_2.residue(pose_2_fragment.get_start()+offset).atom_index("C"), pose_2_fragment.get_start()+offset);
      atom_map.insert(std::pair<id::AtomID, id::AtomID>(id3, id4));

      core::id::AtomID const id5( pose_1.residue(pose_1_fragment.get_start()+offset).atom_index("N"), pose_1_fragment.get_start()+offset);
      core::id::AtomID const id6( pose_2.residue(pose_2_fragment.get_start()+offset).atom_index("N"), pose_2_fragment.get_start()+offset);
      atom_map.insert(std::pair<id::AtomID, id::AtomID>(id5, id6));

      core::id::AtomID const id7( pose_1.residue(pose_1_fragment.get_start()+offset).atom_index("O"), pose_1_fragment.get_start()+offset);
      core::id::AtomID const id8( pose_2.residue(pose_2_fragment.get_start()+offset).atom_index("O"), pose_2_fragment.get_start()+offset);
      atom_map.insert(std::pair<id::AtomID, id::AtomID>(id7, id8));
  }

  return atom_map;

}

//@brief Generate a map of all backbone atom ids between two pairs of helical fragments.
std::map<core::id::AtomID, core::id::AtomID> HelixAssemblyMover::getFragmentPairMap(const core::pose::Pose & pose_1,
    const core::pose::Pose & pose_2, const std::pair<HelicalFragment, HelicalFragment> & pose_1_fragments,
    const std::pair<HelicalFragment, HelicalFragment> & pose_2_fragments){

  std::map< core::id::AtomID, core::id::AtomID > atom_map_1 = getFragmentMap(pose_1, pose_2, pose_1_fragments.first,
      pose_2_fragments.first);
  std::map< core::id::AtomID, core::id::AtomID > atom_map_2 = getFragmentMap(pose_1, pose_2, pose_1_fragments.second,
      pose_2_fragments.second);
  atom_map_1.insert(atom_map_2.begin(), atom_map_2.end());

//  std::map< core::id::AtomID, core::id::AtomID >::iterator it;
//  for ( it=atom_map_1.begin() ; it != atom_map_1.end(); it++ )
//      TR << (*it).first << " => " << (*it).second << endl;

  return atom_map_1;
}

///@details
std::vector<HelixAssemblyJob> HelixAssemblyMover::apply(HelixAssemblyJob & job){

  try{
    std::vector<HelixAssemblyJob> new_jobs;

    TR << "working on file: " << job.get_name() << endl;

    Pose query_structure;
    core::import_pose::pose_from_pdbstring(query_structure, job.get_query_structure());
    TR << "Full Query Structure is " << query_structure.total_residue() << " residues" << endl;

    Pose search_structure;
    core::import_pose::pose_from_pdbstring(search_structure, job.get_search_structure());
    TR << "Search Structure is " << search_structure.total_residue() << " residues" << endl;

    //don't do anything with this round if the search structure is empty. Not testing for this causes all sorts of bad behavior
    if(search_structure.total_residue() <= 0){return new_jobs;}

    this->set_query_frag_1(job.get_query_frag_1());
    this->set_query_frag_2(job.get_query_frag_2());

    std::pair<HelicalFragment, HelicalFragment> query_fragments(get_query_frag_1(), get_query_frag_2());

    TR << "Fragment 1 is " << get_query_frag_1().get_size() << " residues" << endl;
    TR << "Fragment 2 is " << get_query_frag_2().get_size() << " residues" << endl;

    //set up dssp info - necessary in order to find helices based on secondary structure
    core::scoring::dssp::Dssp dssp( search_structure );
    dssp.insert_ss_into_pose( search_structure );

    utility::vector1< std::pair< Size,Size > > helix_endpts;
    helix_endpts = findHelices(search_structure);
    TR << "Found " << helix_endpts.size() << " helices in search structure" << endl;

    //search helix poses for all close RMSD matches to each query fragment
    utility::vector1<HelicalFragment> frag1_matches = findFragmentMatches(search_structure, query_structure,
        get_query_frag_1(), helix_endpts);
    utility::vector1<HelicalFragment> frag2_matches = findFragmentMatches(search_structure, query_structure,
        get_query_frag_2(), helix_endpts);

    TR << "Found " << frag1_matches.size() << " fragments for frag1 & " << frag2_matches.size() << " for frag2." << endl;

    //The max distance-squared a single point can be off an still meet the allowable RMSD requirements
    core::Distance max_point_distance = sqrt(get_helix_pair_rmsd_cutoff() * get_helix_pair_rmsd_cutoff() * query_structure.total_residue());


    bool parallel = get_query_frag_1().get_direction() == get_query_frag_2().get_direction();
    core::Distance end_1_dist;
    core::Distance end_2_dist;
    if(parallel){
        end_1_dist = query_structure.residue(get_query_frag_1().get_start()).atom("CA").xyz().distance(
            query_structure.residue(get_query_frag_2().get_start()).atom("CA").xyz());

        end_2_dist = query_structure.residue(get_query_frag_1().get_end()).atom("CA").xyz().distance(
            query_structure.residue(get_query_frag_2().get_end()).atom("CA").xyz());
    }
    else{
        end_1_dist = query_structure.residue(get_query_frag_1().get_start()).atom("CA").xyz().distance(
            query_structure.residue(get_query_frag_2().get_end()).atom("CA").xyz());

        end_2_dist = query_structure.residue(get_query_frag_1().get_end()).atom("CA").xyz().distance(
            query_structure.residue(get_query_frag_2().get_start()).atom("CA").xyz());
    }

    //Create a vector of only helical pairs that could possibly have an RMSD to the query structure within the allowable range
    utility::vector1< std::pair< HelicalFragment,HelicalFragment > > close_helix_pairs;
    for(Size i=1; i<=frag1_matches.size(); i++){
        for(Size j=1; j<=frag2_matches.size(); j++){
            if(closenessCheck(max_point_distance, end_1_dist, end_2_dist, search_structure, frag1_matches[i], frag2_matches[j])){
                close_helix_pairs.push_back(std::pair<HelicalFragment, HelicalFragment>(frag1_matches[i],frag2_matches[j]));
            }
        }
    }

    TR << close_helix_pairs.size() << " out of " << frag1_matches.size()*frag2_matches.size() <<
        " helix pairs were close enough for further investigation" << endl;

    //keep track of number of hits in structure, for output filename purposes
    Size resultsCounter(1);

    //loop through lists of fragments for each helix and check rmsd to the full query structure
    for(Size i=1; i<=close_helix_pairs.size(); ++i){

        std::map< core::id::AtomID, core::id::AtomID > atom_map = getFragmentPairMap(query_structure, search_structure,
            query_fragments, close_helix_pairs[i]);

        Real bbrmsd = core::scoring::rms_at_all_corresponding_atoms(query_structure, search_structure, atom_map);

        if(bbrmsd <= helix_pair_rmsd_cutoff_){

            TR << "Fragment pair (" << close_helix_pairs[i].first.get_start() << ":" << close_helix_pairs[i].first.get_end() <<
                ") - (" << close_helix_pairs[i].second.get_start() << ":" << close_helix_pairs[i].second.get_end() <<
                ") matches query (RMSD: " << bbrmsd << ")" << endl;

            Pose fragment1(search_structure, close_helix_pairs[i].first.get_start(), close_helix_pairs[i].first.get_end());
            Pose fragment2(search_structure, close_helix_pairs[i].second.get_start(), close_helix_pairs[i].second.get_end());
            Pose combinedResultFragments(combinePoses(fragment1, fragment2));

            utility::vector1<HelicalFragment> helix_partners = findPartnerHelices(search_structure, fragment1, fragment2,
                close_helix_pairs[i], helix_endpts, job.get_first_round(), job.get_direction_needed());

            TR << "found " << helix_partners.size() << " suitable helix partners" << endl;

            for(Size k=1; k<=helix_partners.size(); k++){

                Pose third_helix(search_structure, helix_partners[k].get_start(), helix_partners[k].get_end());

                TR << "New helical pose created." << endl;

                //superimpose the found query helix pair onto the found pair
                superimposeBundles(query_structure, combinedResultFragments);

                TR << "Poses have been superimposed for bundle creation" << endl;

                //combine query structure with third helix
                Size third_helix_start(query_structure.total_residue()+1);
                Pose new_bundle = combinePoses(query_structure, third_helix);
                Size new_chain = new_bundle.chain(new_bundle.total_residue());
                Size third_helix_end(new_bundle.total_residue());

                TR << "New bundle created" << endl;

                //check for backbone clashes introduced by adding the new helix
                Real clash_score = bb_score(new_bundle, new_chain, scorefxn_);

                //TODO (tjacobs) Turn this into an option
                if(clash_score <= 5){

                    HelicalFragment newFragment;
                    newFragment.set_start(third_helix_start);
                    newFragment.set_end(third_helix_end);

                    //capture native residue info from the matching query_fragment and the newly added helix
//                    for(Size m=close_helix_pairs[i].first.get_start(); m<=close_helix_pairs[i].first.get_end(); ++m){
//
//                    }
//
//                    for(Size m=close_helix_pairs[i].first.get_start(); m<=close_helix_pairs[i].first.get_end(); ++m){
//
//                    }

                    job.get_fragments().push_back(newFragment);

                    stringstream tempStream;
                    new_bundle.dump_pdb(tempStream, "");
                    std::string new_bundle_string = tempStream.str();

                    //Create a new job with query fragments set to the new helix and query helix 1
                    HelixAssemblyJob new_job1;
                    new_job1.set_query_structure(new_bundle_string);
                    new_job1.set_remaining_rounds(job.get_remaining_rounds()-1);
                    new_job1.set_query_frag_1_index(job.get_fragments().size()-1);
                    new_job1.set_query_frag_2_index(job.get_query_frag_1_index());
                    new_job1.set_name(job.get_name());
                    new_job1.set_fragments(job.get_fragments());
                    new_job1.set_direction_needed(!helix_partners[k].get_direction());//change direction for next helix

                    new_jobs.push_back(new_job1);

                    //If this was the last round, only return one new job (which is used only for output by the head node)
                    if(new_job1.get_remaining_rounds() > 0){

                        //Create a new job with query fragments set to the new helix and query helix 1
                        HelixAssemblyJob new_job2;
                        new_job2.set_query_structure(new_bundle_string);
                        new_job2.set_remaining_rounds(job.get_remaining_rounds()-1);
                        new_job2.set_query_frag_1_index(job.get_fragments().size()-1);
                        new_job2.set_query_frag_2_index(job.get_query_frag_2_index());
                        new_job2.set_name(job.get_name());
                        new_job2.set_fragments(job.get_fragments());
                        new_job2.set_direction_needed(!helix_partners[k].get_direction());//change direction for next helix

                        new_jobs.push_back(new_job2);
                    }
                }
            }
            resultsCounter++;
        }
    }
    return new_jobs;
  }
  catch(...){
      TR << "ERROR: Exception thrown in HelixAssemblyMover during job: " << job.get_name()  << endl;
      std::vector<HelixAssemblyJob> blankJobs;
      return blankJobs;
  }
}//apply
