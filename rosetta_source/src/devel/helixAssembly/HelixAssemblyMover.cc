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
#include <core/pose/util.tmpl.hh>

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
#include <core/conformation/Atom.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

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

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/helixAssembly.OptionKeys.gen.hh>
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
#include <vector>

using namespace std;
using namespace core;
using namespace core::pose;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer TR("HelixAssemblyMover");

///@brief
HelixAssemblyMover::HelixAssemblyMover(protocols::features::helixAssembly::HelicalFragment query_frag_1, protocols::features::helixAssembly::HelicalFragment query_frag_2) :
    scorefxn_(),
    query_frag_1_(query_frag_1),
    query_frag_2_(query_frag_2),
//    query_structure_path_(),
//    query_structure_string_(),
    single_helix_rmsd_cutoff_(1),
    helix_pair_rmsd_cutoff_(1.3),
    helix_cap_dist_cutoff_(12.0),
    helix_contact_dist_cutoff_(10.0),
    minimum_helix_contacts_(6)
{
  init_from_options();
}

HelixAssemblyMover::~HelixAssemblyMover(){}

core::scoring::ScoreFunctionOP HelixAssemblyMover::get_scorefxn() const
{
    return scorefxn_;
}

protocols::features::helixAssembly::HelicalFragment HelixAssemblyMover::get_query_frag_1() const
{
    return query_frag_1_;
}

protocols::features::helixAssembly::HelicalFragment HelixAssemblyMover::get_query_frag_2() const
{
    return query_frag_2_;
}

Real HelixAssemblyMover::get_helix_cap_dist_cutoff() const
{
    return helix_cap_dist_cutoff_;
}

Real HelixAssemblyMover::get_helix_contact_dist_cutoff() const
{
    return helix_contact_dist_cutoff_;
}

Real HelixAssemblyMover::get_helix_pair_rmsd_cutoff() const
{
    return helix_pair_rmsd_cutoff_;
}

core::Size HelixAssemblyMover::get_minimum_helix_contacts() const
{
    return minimum_helix_contacts_;
}

Real HelixAssemblyMover::get_single_helix_rmsd_cutoff() const
{
    return single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::set_query_frag_1(const protocols::features::helixAssembly::HelicalFragment & frag_1_)
{
    this->query_frag_1_ = frag_1_;
}

void HelixAssemblyMover::set_query_frag_2(const protocols::features::helixAssembly::HelicalFragment & frag_2_)
{
    this->query_frag_2_ = frag_2_;
}

void HelixAssemblyMover::set_scorefxn(const core::scoring::ScoreFunctionOP & scorefxn_)
{
    this->scorefxn_ = scorefxn_;
}

void HelixAssemblyMover::set_helix_cap_dist_cutoff(Real helix_cap_dist_cutoff_)
{
    this->helix_cap_dist_cutoff_ = helix_cap_dist_cutoff_;
}

void HelixAssemblyMover::set_helix_contact_dist_cutoff(Real helix_contact_dist_cutoff_)
{
    this->helix_contact_dist_cutoff_ = helix_contact_dist_cutoff_;
}

void HelixAssemblyMover::set_helix_pair_rmsd_cutoff(Real helix_pair_rmsd_cutoff_)
{
    this->helix_pair_rmsd_cutoff_ = helix_pair_rmsd_cutoff_;
}

void HelixAssemblyMover::set_minimum_helix_contacts(core::Size minimum_helix_contacts_)
{
    this->minimum_helix_contacts_ = minimum_helix_contacts_;
}

void HelixAssemblyMover::set_single_helix_rmsd_cutoff(Real single_helix_rmsd_cutoff_)
{
    this->single_helix_rmsd_cutoff_ = single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::init_from_options(){
  set_scorefxn(core::scoring::getScoreFunction());
//  set_query_structure_path(basic::options::option[HelixAssembly::query_structure_path]);

  set_single_helix_rmsd_cutoff(basic::options::option[helixAssembly::single_helix_rmsd_cutoff]);
  set_minimum_helix_contacts(basic::options::option[helixAssembly::minimum_helix_contacts]);
  set_helix_pair_rmsd_cutoff(basic::options::option[helixAssembly::helix_pair_rmsd_cutoff]);
  set_helix_contact_dist_cutoff(basic::options::option[helixAssembly::helix_contact_dist_cutoff]);
  set_helix_cap_dist_cutoff(basic::options::option[helixAssembly::helix_cap_dist_cutoff]);
}

void HelixAssemblyMover::combinePoses(Pose & pose1, Pose const & pose2){
  if(pose2.total_residue()>=1){
      pose1.append_residue_by_jump(pose2.residue(1), pose1.total_residue() , "", "", true/*start new chain*/);
      for(Size i=2; i<=pose2.total_residue(); i++){
          if(pose2.residue(i).is_lower_terminus()){
        	  pose1.append_residue_by_jump(pose2.residue(i), pose1.total_residue(), "","", true);
          }
          else{
        	  pose1.append_residue_by_bond(pose2.residue(i));
          }
      }
  }
}

utility::vector1<protocols::features::helixAssembly::HelicalFragment> HelixAssemblyMover::findHelices(Pose const & pose){

  utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices;
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

          all_helices.push_back(protocols::features::helixAssembly::HelicalFragment(helix_start, helix_end));
      }
  }
  return all_helices;
}

///@details return all poses from the targetPose that contain an RMSD match to the queryFragment
utility::vector1<protocols::features::helixAssembly::HelicalFragment> HelixAssemblyMover::findFragmentMatches(Pose const & search_structure,
    Pose const & query_structure, protocols::features::helixAssembly::HelicalFragment query_fragment, utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices){

  utility::vector1<protocols::features::helixAssembly::HelicalFragment> frag_matches;

  //iterate through each full-length helix in the search structure
  for(Size j=1; j<=all_helices.size(); j++){

      if(all_helices[j].get_end() > query_fragment.get_size()+1){ //I hate unsigned ints

          for(Size i=all_helices[j].get_start(); i<all_helices[j].get_end()-query_fragment.get_size(); i++){

              //make sure we don't make a test fragment out of two separate chains
              if(search_structure.total_residue() > i+query_fragment.get_size() &&
                  search_structure.residue(i).chain() == search_structure.residue(i+query_fragment.get_size()-1).chain()){

                  protocols::features::helixAssembly::HelicalFragment test_fragment(i, i+query_fragment.get_size()-1);

                  std::map<core::id::AtomID, core::id::AtomID> frag_map = getFragmentMap(search_structure,
                      query_structure, test_fragment, query_fragment);

                  Real rmsd = core::scoring::rms_at_all_corresponding_atoms(search_structure, query_structure, frag_map);

                  if(rmsd <= single_helix_rmsd_cutoff_){
                      frag_matches.push_back(test_fragment);
                  }
              }
          }
      }
  }

  return frag_matches;
}//findFragment

bool HelixAssemblyMover::checkHelixContacts(Pose const & search_structure, std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> helix_pair,
            protocols::features::helixAssembly::HelicalFragment helix_to_check){

  Size distCutoff = pow((double)helix_contact_dist_cutoff_, 2);

  int interactionCounter = 0;
  for(Size i=helix_to_check.get_start(); i<=helix_to_check.get_end(); i++){
      bool frag1Pass=false;
      bool frag2Pass=false;
      for(Size j=helix_pair.first.get_start(); j<=helix_pair.first.get_end(); j++){

          if(search_structure.residue(i).atom("CA").xyz().distance_squared(search_structure.residue(j).atom("CA").xyz()) < distCutoff){
              frag1Pass=true;
              break;
          }
      }
      for(Size j=helix_pair.second.get_start(); j<helix_pair.second.get_end(); j++){

          if(search_structure.residue(i).atom("CA").xyz().distance_squared(search_structure.residue(j).atom("CA").xyz()) < distCutoff){
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
utility::vector1<protocols::features::helixAssembly::HelicalFragment> HelixAssemblyMover::findPartnerHelices(Pose const & search_structure,
    std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> query_match, utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices,
    bool first_round, bool direction_needed){

  utility::vector1<protocols::features::helixAssembly::HelicalFragment> partner_helices;

  bool parallel = query_frag_1_.get_direction() == query_frag_2_.get_direction();

  //look at each helix in the search structure as a potential helix to add
  for(Size i=1; i<=all_helices.size(); i++){

      //don't look at this helix if it contains either of the helices that we used in the initial fragment hits
      if(!((query_match.first.get_start() >= all_helices[i].get_start() &&
          query_match.first.get_start() <= all_helices[i].get_end()) ||
          (query_match.second.get_start() >= all_helices[i].get_start() &&
          query_match.second.get_start() <= all_helices[i].get_end()))){

//          cout << "Checking helix (" << all_helices[i].get_start() << ":" << all_helices[i].get_end() << ") for end point proximity" << endl;

          bool foundStart = false;
          bool foundEnd = false;
          Size helixStart;
          Size helixEnd;
          Real distCutoff = pow(helix_cap_dist_cutoff_,2);
          Real minStartDistance=distCutoff*2;
          Real minEndDistance=distCutoff*2;

          //look for a sub-helix of this full-size helical fragment that minimizes the distance between each end of the
          //query match fragment pair
          for(Size helixOffset=0; helixOffset<all_helices[i].get_size(); helixOffset++){

              core::DistanceSquared startDistance1;
              core::DistanceSquared startDistance2;
              core::DistanceSquared endDistance1;
              core::DistanceSquared endDistance2;
              if(!parallel){//anti-parallel query

                  //distance between n-term of query fragment 1 match and n-term+offset of potential third-helix
                  startDistance1 = search_structure.residue(query_match.first.get_start()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_start()+helixOffset).atom("CA").xyz());

                  //distance between c-term of query fragment 2 match and n-term+offset of potential third-helix
                  startDistance2 = search_structure.residue(query_match.second.get_end()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_start()+helixOffset).atom("CA").xyz());

                  //distance between c-term of query fragment 1 match and c-term-offset of potential third-helix
                  endDistance1 = search_structure.residue(query_match.first.get_end()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_end()-helixOffset).atom("CA").xyz());

                  //distance between n-term of query fragment 2 match and c-term-offset of potential third-helix
                  endDistance2 = search_structure.residue(query_match.second.get_start()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_end()-helixOffset).atom("CA").xyz());
              }
              else{//parallel
                  //distance between n-term of query fragment 1 match and n-term+offset of potential third-helix
                  startDistance1 = search_structure.residue(query_match.first.get_start()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_start()+helixOffset).atom("CA").xyz());

                  //distance between n-term of query fragment 2 match and n-term+offset of potential third-helix
                  startDistance2 = search_structure.residue(query_match.second.get_start()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_start()+helixOffset).atom("CA").xyz());

                  //distance between c-term of query fragment 1 match and c-term-offset of potential third-helix
                  endDistance1 = search_structure.residue(query_match.first.get_end()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_end()-helixOffset).atom("CA").xyz());

                  //distance between c-term of query fragment 1 match and c-term-offset of potential third-helix
                  endDistance2 = search_structure.residue(query_match.second.get_end()).atom("CA").xyz().distance_squared(
                      search_structure.residue(all_helices[i].get_end()-helixOffset).atom("CA").xyz());
              }

              if(startDistance1 < distCutoff && startDistance2 < distCutoff && (startDistance1 + startDistance2) < minStartDistance){
                  helixStart = all_helices[i].get_start()+helixOffset;
                  minStartDistance = startDistance1 + startDistance2;
                  foundStart = true;
//                  cout << "Found third helix start at residue " << helixStart << " in helix (" << all_helices[i].get_start()
//                      << ":" << all_helices[i].get_end() << ")" << endl;
              }

              //Check to make sure that this point in the potential 3rd helix is close to the "end" of both query helices
              if(endDistance1 < distCutoff && endDistance2 < distCutoff && (endDistance1 + endDistance2) < minEndDistance){
                  helixEnd = all_helices[i].get_end()-helixOffset;
                  minEndDistance = endDistance1 + endDistance2;
                  foundEnd = true;
//                  cout << "Found third helix end at residue " << helixEnd << " in helix (" << all_helices[i].get_start()
//                      << ":" << all_helices[i].get_end() << ")" << endl;
              }
          }
          if(foundStart && foundEnd){
              //if the round is 1 we can have a helix in either direction (ie, we can build off either the n or c-terminus).
              //if the round is not 1 then we must alternate directions

              bool direction = query_frag_1_.get_direction();
              if(helixEnd < helixStart){
                  Size temp = helixStart;
                  helixStart = helixEnd;
                  helixEnd = temp;
                  direction = query_frag_2_.get_direction();
              }
              protocols::features::helixAssembly::HelicalFragment closeHelix(helixStart, helixEnd, direction);

//              cout << "Checking helix conacts for match helices (" << query_match.first.get_start() << "," <<
//                  query_match.first.get_end() << ") (" << query_match.second.get_start() << "," << query_match.second.get_end()
//                  << ") and potential third helix (" << closeHelix.get_start() << "," << closeHelix.get_end() << ")." << endl;
              if((first_round || closeHelix.get_direction() == direction_needed) &&
                  checkHelixContacts(search_structure, query_match, closeHelix)){
                  partner_helices.push_back(closeHelix);
              }
          }
      }
  }
  return partner_helices;
}

void HelixAssemblyMover::superimposeBundles(Pose const & query_structure, Pose & search_structure,
		std::pair< protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> matching_pair){
  id::AtomID_Map< core::id::AtomID > atom_map;
  // maps every atomid to bogus atom. Atoms in the query_structure that are left mapped to the bogus atoms will not be used in the superimposition
  atom_map.clear();
  core::pose::initialize_atomid_map( atom_map, search_structure, id::BOGUS_ATOM_ID );

  //superimpose the found two-helix pose onto the query structure
  for(Size j=0; j< query_frag_1_.get_size(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      Size query_frag1_start = query_frag_1_.get_start();
      Size search_frag1_start = matching_pair.first.get_start();

//      cout << "superimposing " << j+query_frag1_start << " of query structure to " << j+search_frag1_start << " of search structure" << endl;

      core::id::AtomID const id1( query_structure.residue(j+query_frag1_start).atom_index("CA"), j+query_frag1_start);
      core::id::AtomID const id2( search_structure.residue(j+search_frag1_start).atom_index("CA"), j+search_frag1_start );
      atom_map[ id2 ] = id1;

      core::id::AtomID const id3( query_structure.residue(j+query_frag1_start).atom_index("C"), j+query_frag1_start);
      core::id::AtomID const id4( search_structure.residue(j+search_frag1_start).atom_index("C"), j+search_frag1_start );
      atom_map[ id4 ] = id3;

      core::id::AtomID const id5( query_structure.residue(j+query_frag1_start).atom_index("N"), j+query_frag1_start);
      core::id::AtomID const id6( search_structure.residue(j+search_frag1_start).atom_index("N"), j+search_frag1_start );
      atom_map[ id6 ] = id5;

      core::id::AtomID const id7( query_structure.residue(j+query_frag1_start).atom_index("O"), j+query_frag1_start);
      core::id::AtomID const id8( search_structure.residue(j+search_frag1_start).atom_index("O"), j+search_frag1_start );
      atom_map[ id8 ] = id7;
  }

  for(Size j=0; j< query_frag_2_.get_size(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...

      core::Size query_frag2_start = query_frag_2_.get_start();
      core::Size search_frag2_start = matching_pair.second.get_start();

//      cout << "superimposing " << j+query_frag2_start << " of query structure to " << j+search_frag2_start << " of search structure" << endl;

      core::id::AtomID const id1( query_structure.residue(j+query_frag2_start).atom_index("CA"), j+query_frag2_start);
      core::id::AtomID const id2( search_structure.residue(j+search_frag2_start).atom_index("CA"),j+search_frag2_start );
      atom_map[ id2 ] = id1;

      core::id::AtomID const id3( query_structure.residue(j+query_frag2_start).atom_index("C"), j+query_frag2_start );
      core::id::AtomID const id4( search_structure.residue(j+search_frag2_start).atom_index("C"), j+search_frag2_start );
      atom_map[ id4 ] = id3;

      core::id::AtomID const id5( query_structure.residue(j+query_frag2_start).atom_index("N"), j+query_frag2_start );
      core::id::AtomID const id6( search_structure.residue(j+search_frag2_start).atom_index("N"), j+search_frag2_start );
      atom_map[ id6 ] = id5;

      core::id::AtomID const id7( query_structure.residue(j+query_frag2_start).atom_index("O"), j+query_frag2_start );
      core::id::AtomID const id8( search_structure.residue(j+search_frag2_start).atom_index("O"), j+search_frag2_start );
      atom_map[ id8 ] = id7;
    }

  //do superimpose
  scoring::superimpose_pose(search_structure, query_structure, atom_map);
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

//  cout<< "Backbone-backbone score: " << bb_energy << std::endl;

  return bb_energy;
}//end bb_score

bool HelixAssemblyMover::closenessCheck(const core::Distance maxRange, const core::Distance end1Dist, const core::Distance end2Dist,
    const core::pose::Pose & search_structure, protocols::features::helixAssembly::HelicalFragment search_frag_1, protocols::features::helixAssembly::HelicalFragment search_frag_2){

  //Calculate the maximum that a point can be off for the rmsd to still meet the threshold, then check the ends of the helices. This
  //is a quick and dirty way to quickly prune the number of helical pairs to do a full RMSD check on.

  //Don't consider this pair of helices if they have overlapping residues
  for(core::Size i=search_frag_1.get_start(); i<search_frag_1.get_end(); i++){
      for(core::Size j=search_frag_2.get_start(); j<search_frag_2.get_end(); j++){

          core::Distance resDistance = search_structure.residue(i).atom("CA").xyz().distance_squared(
              search_structure.residue(j).atom("CA").xyz());

          if(resDistance < 0.5){
              return false;
          }
      }
  }

  //if the query fragments are parallel treat the search fragments the same way & vice-versa
  bool parallel = query_frag_1_.get_direction()==query_frag_2_.get_direction();
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
    const core::pose::Pose & pose_2, protocols::features::helixAssembly::HelicalFragment pose_1_fragment, protocols::features::helixAssembly::HelicalFragment pose_2_fragment){

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
    const core::pose::Pose & pose_2, const std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> & pose_1_fragments,
    const std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> & pose_2_fragments){

  std::map< core::id::AtomID, core::id::AtomID > atom_map_1 = getFragmentMap(pose_1, pose_2, pose_1_fragments.first,
      pose_2_fragments.first);
  std::map< core::id::AtomID, core::id::AtomID > atom_map_2 = getFragmentMap(pose_1, pose_2, pose_1_fragments.second,
      pose_2_fragments.second);
  atom_map_1.insert(atom_map_2.begin(), atom_map_2.end());

//  std::map< core::id::AtomID, core::id::AtomID >::iterator it;
//  for ( it=atom_map_1.begin() ; it != atom_map_1.end(); it++ )
//      cout << (*it).first << " => " << (*it).second << endl;

  return atom_map_1;
}

void HelixAssemblyMover::removeDuplicateFragmentPairs(const core::pose::Pose & pose,
		utility::vector1< std::pair<protocols::features::helixAssembly::HelicalFragment,protocols::features::helixAssembly::HelicalFragment> > & helix_pairs){
	for(core::Size i=1; i<=helix_pairs.size(); ++i){
		for(core::Size j=i+1; j<=helix_pairs.size(); ++j){


			std::map<core::id::AtomID, core::id::AtomID> atom_map(getFragmentPairMap(pose, pose, helix_pairs[i], helix_pairs[j]));
			core::Real rms(core::scoring::rms_at_all_corresponding_atoms(pose, pose, atom_map));


			std::string frag_pair_1(pose.sequence().substr(helix_pairs[i].first.get_start(), helix_pairs[i].first.get_size()) + " " +
					pose.sequence().substr(helix_pairs[i].second.get_start(), helix_pairs[i].second.get_size()));
			std::string frag_pair_2(pose.sequence().substr(helix_pairs[j].first.get_start(), helix_pairs[j].first.get_size()) + " " +
					pose.sequence().substr(helix_pairs[j].second.get_start(), helix_pairs[j].second.get_size()));

//			cout << frag_pair_1 << endl;
//			cout << frag_pair_2 << endl;
//			cout << "RMS: " << rms << endl;

			if(rms < 0.1){
				helix_pairs.erase(helix_pairs.begin()+j-1);
			}
		}
	}
}

void HelixAssemblyMover::removeDuplicateFragments(const core::pose::Pose & pose,
		const utility::vector1<protocols::features::helixAssembly::HelicalFragment> & all_helix_fragments,
		utility::vector1<protocols::features::helixAssembly::HelicalFragment> & helix_fragments){
	for(core::Size i=1; i<=all_helix_fragments.size(); ++i){
		for(core::Size j=1; j<=helix_fragments.size(); ++j){
			bool remove_j=true;

			std::string frag_pair_1(pose.sequence().substr(all_helix_fragments[i].get_start(), all_helix_fragments[i].get_size()));
			std::string frag_pair_2(pose.sequence().substr(helix_fragments[j].get_start(), helix_fragments[j].get_size()));

			if(all_helix_fragments[i].get_size() == helix_fragments[j].get_size()){
				for(core::Size offset=0; offset<all_helix_fragments[i].get_size(); offset++){
					core::Size i_res=all_helix_fragments[i].get_start()+offset;
					core::Size j_res=helix_fragments[j].get_start()+offset;
					if(pose.residue(i_res).aa() != pose.residue(j_res).aa()){
						remove_j=false;
					}
				}
				if(remove_j){
					helix_fragments.erase(helix_fragments.begin()+j-1);
				}
			}
		}
	}
}

///@details
std::vector<HelixAssemblyJob> HelixAssemblyMover::apply(HelixAssemblyJob & job){

	try{
		std::vector<HelixAssemblyJob> new_jobs;
        int totFrag1Hits=0;
        int totFrag2Hits=0;
        int totCombinedHits=0;

		TR << "working on file: " << job.get_name() << endl;

		if(job.get_search_structure().length()>=1000000){
			return new_jobs;
		}

		Pose search_structure;
		core::import_pose::pose_from_pdbstring(search_structure, job.get_search_structure());
		TR << "Search Structure is " << search_structure.total_residue() << " residues" << endl;

		//don't do anything with this round if the search structure is empty. Not testing for this causes all sorts of bad behavior
		if(search_structure.total_residue() <= 0){return new_jobs;}

		//set up dssp info - necessary in order to find helices based on secondary structure
		core::scoring::dssp::Dssp dssp( search_structure );
		dssp.insert_ss_into_pose( search_structure );

		utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helices;
		all_helices = findHelices(search_structure);
		TR << "Found " << all_helices.size() << " helices in search structure" << endl;

		//If there aren't 3 separate helices in the structure then there's not enough information to use
		if(all_helices.size() <= 2){return new_jobs;}

		Pose query_structure;
		core::import_pose::pose_from_pdbstring(query_structure, job.get_query_structure());
		TR << "Full Query Structure is " << query_structure.total_residue() << " residues" << endl;

		//If this is the first round we need the original query structure sequence data
		if(job.get_first_round()){
			    	query_frag_1_.insertResiduesFromPose(query_structure, query_frag_1_.get_start(), query_frag_1_.get_end(), query_structure);
			    	query_frag_2_.insertResiduesFromPose(query_structure, query_frag_2_.get_start(), query_frag_2_.get_end(), query_structure);
		}

		std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment> query_fragments(query_frag_1_, query_frag_2_);
		cout << "Fragment 1 is " << query_frag_1_.get_size() << " residues" << endl;
		cout << "Fragment 2 is " << query_frag_2_.get_size() << " residues" << endl;

		//search helix poses for all close RMSD matches to each query fragment
		utility::vector1<protocols::features::helixAssembly::HelicalFragment> frag1_matches = findFragmentMatches(search_structure, query_structure,
				query_frag_1_, all_helices);
		utility::vector1<protocols::features::helixAssembly::HelicalFragment> frag2_matches = findFragmentMatches(search_structure, query_structure,
				query_frag_2_, all_helices);
		cout << "Found " << frag1_matches.size() << " fragments for frag1 & " << frag2_matches.size() << " for frag2." << endl;
        totFrag1Hits+=frag1_matches.size();
        totFrag2Hits+=frag2_matches.size();

		//The max distance-squared a single point can be off an still meet the allowable RMSD requirements
		core::Distance max_point_distance = sqrt(get_helix_pair_rmsd_cutoff() * get_helix_pair_rmsd_cutoff() * query_structure.total_residue());

		bool parallel = query_frag_1_.get_direction() == query_frag_2_.get_direction();
		core::Distance end_1_dist;
		core::Distance end_2_dist;
		if(parallel){
			end_1_dist = query_structure.residue(query_frag_1_.get_start()).atom("CA").xyz().distance(
					query_structure.residue(query_frag_2_.get_start()).atom("CA").xyz());

			end_2_dist = query_structure.residue(query_frag_1_.get_end()).atom("CA").xyz().distance(
					query_structure.residue(query_frag_2_.get_end()).atom("CA").xyz());
		}
		else{
			end_1_dist = query_structure.residue(query_frag_1_.get_start()).atom("CA").xyz().distance(
					query_structure.residue(query_frag_2_.get_end()).atom("CA").xyz());

			end_2_dist = query_structure.residue(query_frag_1_.get_end()).atom("CA").xyz().distance(
					query_structure.residue(query_frag_2_.get_start()).atom("CA").xyz());
		}

		//Create a vector of only helical pairs that could possibly have an RMSD to the query structure within the allowable range
		utility::vector1< std::pair< protocols::features::helixAssembly::HelicalFragment,protocols::features::helixAssembly::HelicalFragment > > close_helix_pairs;
		for(Size i=1; i<=frag1_matches.size(); i++){
			for(Size j=1; j<=frag2_matches.size(); j++){
				if(closenessCheck(max_point_distance, end_1_dist, end_2_dist, search_structure, frag1_matches[i], frag2_matches[j])){
					close_helix_pairs.push_back(std::pair<protocols::features::helixAssembly::HelicalFragment, protocols::features::helixAssembly::HelicalFragment>(frag1_matches[i],frag2_matches[j]));
				}
			}
		}

		TR.Debug << close_helix_pairs.size() << " out of " << frag1_matches.size()*frag2_matches.size() <<
				" helix pairs were close enough for further investigation" << endl;

		removeDuplicateFragmentPairs(search_structure, close_helix_pairs);

		TR.Debug << close_helix_pairs.size() << " helix pairs left after duplication removal" << endl;

		//loop through lists of fragments for each helix and check rmsd to the full query structure
		Size new_helix_counter(0);
		utility::vector1<protocols::features::helixAssembly::HelicalFragment> all_helix_partners;
		for(Size i=1; i<=close_helix_pairs.size(); ++i){

//			cout << "examining close helix pair (" << close_helix_pairs[i].first.get_start() << "," << close_helix_pairs[i].first.get_end() <<
//					") (" << close_helix_pairs[i].second.get_start() << "," << close_helix_pairs[i].second.get_end() <<
//					")" << endl;

			std::map< core::id::AtomID, core::id::AtomID > atom_map = getFragmentPairMap(query_structure, search_structure,
					query_fragments, close_helix_pairs[i]);

			Real bbrmsd = core::scoring::rms_at_all_corresponding_atoms(query_structure, search_structure, atom_map);
//			cout << "bbrmsd: " << bbrmsd << endl;

			if(bbrmsd <= helix_pair_rmsd_cutoff_){
                totCombinedHits++;

//				cout << "Fragment pair (" << close_helix_pairs[i].first.get_start() << ":" << close_helix_pairs[i].first.get_end() <<
//						") - (" << close_helix_pairs[i].second.get_start() << ":" << close_helix_pairs[i].second.get_end() <<
//						") matches query (RMSD: " << bbrmsd << ")" << endl;

				//Find all helices that can be used as the third helix for the current close_helix_pair
				utility::vector1<protocols::features::helixAssembly::HelicalFragment> helix_partners = findPartnerHelices(search_structure,
						close_helix_pairs[i], all_helices, job.get_first_round(), job.get_direction_needed());

//				cout << "found " << helix_partners.size() << " suitable helix partners" << endl;

				removeDuplicateFragments(search_structure, all_helix_partners, helix_partners);
				all_helix_partners.insert(all_helix_partners.end(), helix_partners.begin(), helix_partners.end());

//				cout << helix_partners.size() << " partners left after duplication removal" << endl;

				for(Size k=1; k<=helix_partners.size(); k++){

					//superimpose the search structure onto the query structure so a helix can be stolen for bundle creation
					superimposeBundles(query_structure, search_structure, close_helix_pairs[i]);
//					cout << "Poses have been superimposed for bundle creation" << endl;

					Pose third_helix(search_structure, helix_partners[k].get_start(), helix_partners[k].get_end());
//					cout << "New helical pose created." << endl;

					if(job.get_first_round()){
						if(helix_partners[k].get_direction() == query_frag_1_.get_direction()){
							job.set_n_term_growth(false);
						}
						else{
							job.set_n_term_growth(true);
						}
					}

					Pose new_bundle;
					Size new_chain;
					Size third_helix_start;
					Size third_helix_end;

					//if we are growing from the n-term then each new helix is the new pose start.
					if(job.get_n_term_growth()){
						new_bundle=third_helix;
						combinePoses(new_bundle, query_structure);
						new_chain = new_bundle.chain(1);
						third_helix_start = 1;
						third_helix_end = helix_partners[k].get_size();
					}
					//if we are growing from the c-term then each new helix is the last chain.
					else{
						new_bundle=query_structure;
						combinePoses(new_bundle, third_helix);
						new_chain = new_bundle.chain(new_bundle.total_residue());
						third_helix_start = query_structure.total_residue()+1;
						third_helix_end = new_bundle.total_residue();
					}

					//check for backbone clashes introduced by adding the new helix
					Real clash_score = bb_score(new_bundle, new_chain, scorefxn_);

					//TODO (tjacobs) Turn this into an option
					if(clash_score <= 5){

						new_helix_counter++;

						protocols::features::helixAssembly::HelicalFragment newFragment(third_helix_start, third_helix_end, helix_partners[k].get_direction());

						std::vector<protocols::features::helixAssembly::HelicalFragment> new_fragment_list = job.get_fragments();

						//capture native residue info from helix pair matching query_fragment_1

						protocols::features::helixAssembly::HelicalFragment query_frag_1_copy(query_frag_1_);
						protocols::features::helixAssembly::HelicalFragment query_frag_2_copy(query_frag_2_);

						query_frag_1_copy.insertResiduesFromPose(search_structure, close_helix_pairs[i].first.get_start(),
							close_helix_pairs[i].first.get_end(), new_bundle);
						new_fragment_list[job.get_query_frag_1_index()] = query_frag_1_copy;

//						cout << "Just stole from query frag 1" << endl;

						//                    search_structure.dump_pdb("reference.pdb");
						//                    new_bundle.dump_pdb("this_pose.pdb");
						//                    exit(1);

						query_frag_2_copy.insertResiduesFromPose(search_structure, close_helix_pairs[i].second.get_start(),
							close_helix_pairs[i].second.get_end(), new_bundle);
						new_fragment_list[job.get_query_frag_2_index()] = query_frag_2_copy;

//						cout << "Just stole from query frag 2" << endl;

						newFragment.insertResiduesFromPose(search_structure, helix_partners[k].get_start(),
							helix_partners[k].get_end(), new_bundle);
						new_fragment_list.push_back(newFragment);

//						cout << "Just stole from new fragment" << endl;

						stringstream tempStream;
						new_bundle.dump_pdb(tempStream, "");
						std::string new_bundle_string = tempStream.str();

						//Create a new job with query fragments set to the new helix and query helix 1
						HelixAssemblyJob new_job1;
						new_job1.set_query_structure(new_bundle_string);
						new_job1.set_remaining_rounds(job.get_remaining_rounds()-1);
						new_job1.set_query_frag_1_index(new_fragment_list.size()-1);//last element in std::vector
						new_job1.set_query_frag_2_index(job.get_query_frag_1_index());
						new_job1.set_name(job.get_name()+"_"+utility::to_string(new_helix_counter));
						new_job1.set_fragments(new_fragment_list);
						new_job1.set_direction_needed(!newFragment.get_direction());//change direction for next helix
						new_job1.set_n_term_growth(job.get_n_term_growth());
						new_job1.set_first_round(false);
                        


						new_jobs.push_back(new_job1);

//						cout << "new job 1: " << new_job1.get_query_frag_1().get_start() << "," << new_job1.get_query_frag_1().get_end() <<
//								" - " << new_job1.get_query_frag_2().get_start() << "," << new_job1.get_query_frag_2().get_end() << endl;

						//If this was the last round, only return one new job (which is used only for output by the head node)
						if(new_job1.get_remaining_rounds() > 0){

							//Create a new job with query fragments set to the new helix and query helix 1
							HelixAssemblyJob new_job2;
							new_job2.set_query_structure(new_bundle_string);
							new_job2.set_remaining_rounds(job.get_remaining_rounds()-1);
							new_job2.set_query_frag_1_index(new_fragment_list.size()-1);//last element in std::vector
							new_job2.set_query_frag_2_index(job.get_query_frag_2_index());
							new_job2.set_name(job.get_name()+"_"+utility::to_string(new_helix_counter));
							new_job2.set_fragments(new_fragment_list);
							new_job2.set_direction_needed(!newFragment.get_direction());//change direction for next helix
							new_job2.set_n_term_growth(job.get_n_term_growth());
							new_job2.set_first_round(false);

							new_jobs.push_back(new_job2);

//							cout << "new job 2: " << new_job2.get_query_frag_1().get_start() << "," << new_job2.get_query_frag_1().get_end() <<
//									" - " << new_job2.get_query_frag_2().get_start() << "," << new_job2.get_query_frag_2().get_end() << endl;
						}
					}
				}
			}
		}
        
        TR << "Total hits for fragment 1: " << totFrag1Hits << endl;
        TR << "Total hits for fragment 2: " << totFrag2Hits << endl;
        TR << "Total combined hits: " << totCombinedHits << endl;
        
		return new_jobs;
	}
	catch(...){
		cout << "ERROR: Exception thrown in HelixAssemblyMover during job: " << job.get_name()  << endl;
		std::vector<HelixAssemblyJob> blankJobs;
		return blankJobs;
	}
}//apply
