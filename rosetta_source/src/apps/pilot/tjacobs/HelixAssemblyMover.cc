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
#include <devel/init.hh>
#include <apps/pilot/tjacobs/HelixAssemblyMover.hh>

//core library
#include <core/init.hh>
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

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <core/scoring/hbonds/hbonds.hh>
#include <core/scoring/hbonds/HBondSet.hh>
#include <core/scoring/dssp/Dssp.hh>

#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.hh>

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

//utility & numeric
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/vector1.functions.hh>
#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <numeric/xyzVector.hh>

//protocols
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/StructureRestrictor.hh>
#include <core/pose/metrics/CalculatorFactory.hh>

//JD2
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// option key includes
#include <basic/options/util.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>
//#include <core/options/keys/packing.OptionKeys.gen.hh>
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

//static basic::Tracer TR( "HelixAssemblyMover" );
using namespace std;
using namespace core;
using namespace core::pose;

static basic::Tracer TR("HelixAssemblyMover");

//application specific options
namespace HelixAssembly{
        basic::options::FileOptionKey const query_structure_path( "query_structure" ); //
        basic::options::FileOptionKey const frag1_path( "frag1" ); //
        basic::options::FileOptionKey const frag2_path( "frag2" ); //
        basic::options::RealOptionKey const single_helix_rmsd_cutoff( "single_helix_rmsd_cutoff" );
        basic::options::RealOptionKey const helix_pair_rmsd_cutoff( "helix_pair_rmsd_cutoff" );
        basic::options::RealOptionKey const helix_cap_distance_cutoff( "helix_cap_distance_cutoff" );
        basic::options::RealOptionKey const helix_contact_distance_cutoff( "helix_contact_distance_cutoff" );
        basic::options::IntegerOptionKey const minimum_helix_contacts( "minimum_helix_contacts" );
}

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

std::string HelixAssemblyMover::get_frag1_path() const
{
    return frag1_path_;
}

std::string HelixAssemblyMover::get_frag2_path() const
{
    return frag2_path_;
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

Real HelixAssemblyMover::get_single_helix_rmsd_cutoff() const
{
    return single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::set_frag1_path(std::string frag1_path_)
{
    this->frag1_path_ = frag1_path_;
}

void HelixAssemblyMover::set_frag2_path(std::string frag2_path_)
{
    this->frag2_path_ = frag2_path_;
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

void HelixAssemblyMover::set_single_helix_rmsd_cutoff(Real single_helix_rmsd_cutoff_)
{
    this->single_helix_rmsd_cutoff_ = single_helix_rmsd_cutoff_;
}

void HelixAssemblyMover::init_from_options(){
  set_query_structure_path(basic::options::option[HelixAssembly::query_structure_path]);
  set_frag1_path(basic::options::option[HelixAssembly::frag1_path]);
  set_frag2_path(basic::options::option[HelixAssembly::frag2_path]);
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
      newPose.append_residue_by_jump(pose2.residue(1), newPose.total_residue() , "", "", false/*start new chain*/);
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
  TR << "target pose size: " << targetPose.total_residue() << endl;
  TR << "query fragment size: " << queryFragment.total_residue() << endl;
  for(Size j=1; j<=helix_endpts.size(); j++){
      if(helix_endpts[j].first + queryFragment.total_residue() - 1 <= helix_endpts[j].second){
          for(Size i=helix_endpts[j].first; i<=helix_endpts[j].second-queryFragment.total_residue()+1; i++){

              //make sure we don't make a test fragment out of two separate chains
              if(targetPose.total_residue() > i+queryFragment.total_residue() &&
                  targetPose.residue(i).chain() == targetPose.residue(i+queryFragment.total_residue()-1).chain()){

                  Pose testFragment(targetPose, i, i+queryFragment.total_residue()-1);

                  Real rmsd = core::scoring::bb_rmsd_including_O(queryFragment, testFragment);

                  if(rmsd <= single_helix_rmsd_cutoff_){
                      //TR << "found matching fragment(" << i << "): " << rmsd << endl;
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

  TR << "testHelix residues: " << testHelix.total_residue() << endl;
  TR << "interactionCounter: " << interactionCounter << endl;

  if(interactionCounter >= minimum_helix_contacts_){
      return true;
  }
  return false;
}

///@details search the pose for a helical segment that makes interactions with both of the pose fragments
utility::vector1< std::pair< Size,Size > > HelixAssemblyMover::findPartnerHelices(Pose const & pose, Pose const & fragment1,
    Pose const & fragment2, Size frag1Start, Size frag2Start, utility::vector1< std::pair< Size,Size > > helix_endpts){

  utility::vector1< std::pair< Size,Size > > partnerHelices;

  for(Size i=1; i<=helix_endpts.size(); i++){

      //don't look at this helix if it contains either of the helices that we used in the initial fragment hits
      if(!((frag1Start >= helix_endpts[i].first && frag1Start <= helix_endpts[i].second) ||
          (frag2Start >= helix_endpts[i].first && frag2Start <= helix_endpts[i].second))){

          bool foundStart = false;
          bool foundEnd = false;
          Size helixStart;
          Size helixEnd;
          Size distCutoff = pow((double)helix_cap_distance_cutoff_,2);
          //Search for helical fragments that have n-terms of query helix close to n-terms frag1 & c-terms of frag2.
          for(Size helixOffset=0; helixOffset<=(helix_endpts[i].second-helix_endpts[i].first+1)/2; helixOffset++){
              //                          TR << "helix -" << helix_endpts[k].first+helixOffset << ":" << helix_endpts[k].second-helixOffset << std::endl;

              if(!foundStart){
                  core::DistanceSquared startDistance1 = fragment1.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  //                              TR << "start distance1: " << startDistance1 << endl;

                  core::DistanceSquared startDistance2 = fragment2.residue(fragment2.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  //                              TR << "start distance2: " << startDistance2 << endl;

                  if(startDistance1 < distCutoff && startDistance2 < distCutoff){
                      helixStart = helix_endpts[i].first+helixOffset;
                      foundStart = true;
                  }
              }

              if(!foundEnd){
                  core::DistanceSquared endDistance1 = fragment1.residue(fragment1.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  //                              TR << "end distance1: " << endDistance1 << endl;

                  core::DistanceSquared endDistance2 = fragment2.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  //                              TR << "end distance2: " << endDistance2 << endl;

                  if(endDistance1 < distCutoff && endDistance2 < distCutoff){
                      helixEnd = helix_endpts[i].second-helixOffset;
                      foundEnd = true;
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
              //                          TR << "helix -" << helix_endpts[k].first+helixOffset << ":" << helix_endpts[k].second-helixOffset << std::endl;

              if(!foundStart){
                  core::DistanceSquared startDistance1 = fragment1.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  //                              TR << "start distance1: " << startDistance1 << endl;

                  core::DistanceSquared startDistance2 = fragment2.residue(fragment2.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].second-helixOffset).atom("CA").xyz());
                  //                              TR << "start distance2: " << startDistance2 << endl;

                  if(startDistance1 < distCutoff && startDistance2 < distCutoff){
                      helixStart = helix_endpts[i].first+helixOffset;
                      foundStart = true;
                  }
              }

              if(!foundEnd){
                  core::DistanceSquared endDistance1 = fragment1.residue(fragment1.total_residue()).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  //                              TR << "end distance1: " << endDistance1 << endl;

                  core::DistanceSquared endDistance2 = fragment2.residue(1).atom("CA").xyz().distance_squared(
                      pose.residue(helix_endpts[i].first+helixOffset).atom("CA").xyz());
                  //                              TR << "end distance2: " << endDistance2 << endl;

                  if(endDistance1 < distCutoff && endDistance2 < distCutoff){
                      helixEnd = helix_endpts[i].second-helixOffset;
                      foundEnd = true;
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
      }
  }
  return partnerHelices;
}

void HelixAssemblyMover::superimposeBundles(Pose & pose1, Pose const & pose2){
  id::AtomID_Map< core::id::AtomID > atom_map;
  // maps every atomid to bogus atom
  atom_map.clear();
  core::pose::initialize_atomid_map( atom_map, pose1, id::BOGUS_ATOM_ID );

  //superimpose the found two-helix pose onto the query structure
  for(Size j=1; j<=pose1.total_residue(); j++){
      //Sloppy way of adding all bb atoms, must be a better way...
      core::id::AtomID const id1( pose1.residue(j).atom_index("CA"), j);
      core::id::AtomID const id2( pose2.residue(j).atom_index("CA"), j );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( pose1.residue(j).atom_index("C"), j );
      core::id::AtomID const id4( pose2.residue(j).atom_index("C"), j );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( pose1.residue(j).atom_index("N"), j );
      core::id::AtomID const id6( pose2.residue(j).atom_index("N"), j );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( pose1.residue(j).atom_index("O"), j );
      core::id::AtomID const id8( pose2.residue(j).atom_index("O"), j );
      atom_map[ id7 ] = id8;
  }

  //do superimpose
  scoring::superimpose_pose(pose1, pose2, atom_map);
}

///@details
void HelixAssemblyMover::apply( Pose & pose ){

  utility::file::FileName filename (pose.pdb_info()->name());
  TR << "working on file: " << filename << "\n";
  std::string baseOutputName = filename.base() + "_hit_";

  Pose fullQueryStructure;
  Pose helixFragment1;
  Pose helixFragment2;

  core::import_pose::pose_from_pdb(fullQueryStructure, get_query_structure_path(), false);
  core::import_pose::pose_from_pdb(helixFragment1, get_frag1_path(), false);
  core::import_pose::pose_from_pdb(helixFragment2, get_frag2_path(), false);

  //set up dssp info - necessary in order to find helices based on secondary structure
  core::scoring::dssp::Dssp dssp( pose );
  dssp.insert_ss_into_pose( pose );

  utility::vector1< std::pair< Size,Size > > helix_endpts;
  helix_endpts = findHelices(pose);

  //search helix poses for all close RMSD matches to each query fragment
  utility::vector1<Size> fragment1Results = findFragments(pose, helixFragment1, helix_endpts);
  utility::vector1<Size> fragment2Results = findFragments(pose, helixFragment2, helix_endpts);

  //keep track of number of hits in structure, for output filename purposes
  Size resultsCounter;
  resultsCounter = 1;

  Size threeHelixCounter;
  threeHelixCounter = 1;

  //loop through lists of fragments for each helix and check rmsd to the full query structure
  for(Size i=1; i<=fragment1Results.size(); i++){

      for(Size j=1; j<=fragment2Results.size(); j++){
          Pose fragment1(pose, fragment1Results[i], fragment1Results[i]+helixFragment1.total_residue()-1);
          Pose fragment2(pose, fragment2Results[j], fragment2Results[j]+helixFragment2.total_residue()-1);

          //combine the poses in both possible ways to maximize RMSD
          Pose combinedResultFragments(combinePoses(fragment1, fragment2));

          Real bbrmsd = core::scoring::bb_rmsd_including_O(fullQueryStructure, combinedResultFragments);

          Real carmsd;

//          TR << "fragment 1: " << fragment1Results[i] << endl;
//          TR << "fragment 2: " << fragment2Results[j] << endl;
//          fragment1.dump_pdb("fragment1.pdb");
//          fragment2.dump_pdb("fragment2.pdb");
//          combinedResultFragments.dump_pdb("combined.pdb");
//          exit(1);

//          TR << "fragment 1: " << fragment1Results[i] << endl;
//          TR << "fragment 2: " << fragment2Results[j] << endl;


          //Do the combined fragments match the full query structure
//          if(bbrmsd <= rmsd_cutoff_){
          if(bbrmsd <= helix_pair_rmsd_cutoff_){

              TR << "combined bb_atom RMSD(" << i << ", " << j << "): " << bbrmsd << endl;

              utility::vector1< std::pair< Size,Size > > helix_partners = findPartnerHelices(pose, fragment1, fragment2,
                  fragment1Results[i], fragment2Results[j], helix_endpts);

              for(Size k=1; k<=helix_partners.size(); k++){

                  Pose thirdHelix(pose, helix_partners[k].first, helix_partners[k].second);

                  //superimpose the found query helix pair onto the found pair
                  superimposeBundles(fullQueryStructure, combinedResultFragments);

                  //combine query structure with third helix
                  Pose threeHelixBundle = combinePoses(fullQueryStructure, thirdHelix);

                  //output
                  std::string threeHelixOutput = baseOutputName + utility::to_string(resultsCounter) +
                      "_threeHelix_" + utility::to_string(threeHelixCounter) +".pdb";
                  TR << "outputting file: " << threeHelixOutput << "\n";
                  threeHelixBundle.dump_pdb(threeHelixOutput);
                  threeHelixCounter++;
              }

              resultsCounter++;
          }
      }
  }

}//apply

// run protocol
int
main( int argc, char * argv [] )
{
        using namespace basic::options;
        using namespace basic::options::OptionKeys;

        option.add( HelixAssembly::query_structure_path, "" );
        option.add( HelixAssembly::frag1_path, "" );
        option.add( HelixAssembly::frag2_path, "" );
        option.add( HelixAssembly::single_helix_rmsd_cutoff, "");
        option.add( HelixAssembly::helix_pair_rmsd_cutoff, "");
        option.add( HelixAssembly::helix_cap_distance_cutoff, "");
        option.add( HelixAssembly::helix_contact_distance_cutoff, "");
        option.add( HelixAssembly::minimum_helix_contacts, "");

        // initialize core
        devel::init(argc, argv);

        protocols::jd2::JobDistributor::get_instance()->go( new HelixAssemblyMover() );

        std::cout << "Done! -------------------------------\n";
}
