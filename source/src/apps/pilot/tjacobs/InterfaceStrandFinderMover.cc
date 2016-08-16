// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/jacobs/InterfaceStrandFinder.cc
/// @brief uses a library of known gtpase effectors to locate pdb files with beta strands that may be good
/// candidates for beta-strand interface design with gtpases
/// @author tim jacobs

//Unit headers
#include <devel/init.hh>

//core library
#include <devel/init.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/pdb/pdb_writer.hh>
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

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
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


static THREAD_LOCAL basic::Tracer TR( "InterfaceStrandFinder" );

using namespace core;
using namespace utility;
using namespace protocols;
using namespace protocols::moves;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::jd2;

// application specific options
namespace InterfaceStrandFinder {
  RealOptionKey const is_max_RMSD( "is_max_RMSD" ); //max rmsd allowed
  FileVectorOptionKey const is_known_strands( "is_known_strands" ); //file containing pdb names of known strands
  FileVectorOptionKey const is_complexes( "is_complexes" ); //file containing pdb names of known strands
  IntegerOptionKey const is_beta_length( "is_beta_length" ); // min beta lenght to consider
  RealOptionKey const is_vector_length( "is_vector_length" ); // length of the vector to draw from the carboxyl oxygen for strand exposure
  RealOptionKey const is_exposure_radius( "is_exposure_radius" ); // radius of sphere to draw at the query point for strand exposure
  RealOptionKey const is_maxE( "is_maxE" ); // radius of sphere to draw at the query point for strand exposure
}

// mover definition
class InterfaceStrandFinderMover : public Mover {
public:

  InterfaceStrandFinderMover();

  virtual core::Real bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn);

  virtual std::pair<pose::Pose, Size> dock_strands(Pose known_complex, Size known_start, Pose pose, Size found_start);

  virtual Size find_known_strand_start(Pose known_complex, Pose knownStrand);

  virtual bool is_strand_exposed(core::pose::Pose pose, Size startRes, Size endRes);

  virtual void apply( core::pose::Pose& pose );

  virtual MoverOP clone() const {
    return new InterfaceStrandFinderMover( *this );
  }

  virtual
  std::string
  get_name() const {
    return "InterfaceStrandFinderMover";
  }

  virtual	MoverOP	fresh_instance() const {
    return new InterfaceStrandFinderMover();
  }

private:
  scoring::ScoreFunctionOP scorefxn_;
  core::Real max_RMSD_;
  core::Size beta_length_;
  core::Real vector_length_;
  core::Real exposure_radius_;
  core::Real maxE_;
};

//constructor
InterfaceStrandFinderMover::InterfaceStrandFinderMover() {
  scorefxn_ = core::scoring::get_score_function();
  beta_length_ = option[InterfaceStrandFinder::is_beta_length].def(4);
  max_RMSD_= option[InterfaceStrandFinder::is_max_RMSD].def(1);
  vector_length_ = option[InterfaceStrandFinder::is_vector_length].def(3);
  exposure_radius_ = option[InterfaceStrandFinder::is_exposure_radius].def(3);
  maxE_ = option[InterfaceStrandFinder::is_maxE].def(5);
}

Size InterfaceStrandFinderMover::find_known_strand_start(Pose known_complex, Pose knownStrand){
  for(Size i=1; i<=known_complex.total_residue(); i++){
      if(known_complex.pdb_info()->chain(i) == knownStrand.pdb_info()->chain(1) &&
          known_complex.pdb_info()->number(i) == knownStrand.pdb_info()->number(1)){
          return i;
      }
  }
  utility_exit_with_message("known strand not found in known complex (are you sure they match)\n");
  return 0;
}

///////////////////////////////////////
// bb score
///////////////////////////////////////
core::Real InterfaceStrandFinderMover::bb_score(pose::Pose & pose, core::Size aligned_chain_num, core::scoring::ScoreFunctionOP & scorefxn){

  // score the bb-bb energy between chains
  // This part written by P.Doug Renfrew
  // make vectors of bb atoms in each chain individually
  // the master pose will always be chain 1.
  // need to make a vector of all atoms in the chain you are comparing too

  utility::vector1<core::conformation::Atom> chain1_bb_atoms;
  utility::vector1<core::conformation::Atom> chain2_bb_atoms;
  utility::vector1<core::conformation::Atom> all_bb_atoms;

  for( Size j = 1; j <= pose.total_residue(); ++j ) {
      core::conformation::Residue const & res( pose.residue(j) );
      core::chemical::AtomIndices bb_ai( res.mainchain_atoms() );
      //assert( bb_ai.size() == 4 );
      core::Size chain_num( res.chain() );
      for( Size jj = 1; jj <= bb_ai.size(); ++jj ) {
          if( chain_num == 1 )
            chain1_bb_atoms.push_back( res.atom(jj) );

          else if( chain_num == aligned_chain_num )
            chain2_bb_atoms.push_back( res.atom(jj) );
          //optional get all the atoms not in allinged chain,
          //only need to do if more than two chains in pose

          if( pose.conformation().num_chains() >= 3 && chain_num != 1 )
            all_bb_atoms.push_back( res.atom(jj) );
          //end optional
      }
  }

  //NOW SCORE!
  // get instance of etable energy method
  core::scoring::methods::EnergyMethodOptions const & emo(scorefxn->energy_method_options());
  core::scoring::etable::Etable const & et(*(core::scoring::ScoringManager::get_instance()->etable(emo.etable_type())));
  core::scoring::etable::TableLookupEtableEnergy ete( et, emo );

  // iterate over both sets of atom and add into one emapvector
  //core::scoring::TwoBodyEMapVector tbemv;
  core::scoring::EMapVector tbemv;
  core::Real atr_wt( (*scorefxn).get_weight(core::scoring::fa_atr) );
  core::Real rep_wt( (*scorefxn).get_weight(core::scoring::fa_rep) );
  for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ){
      for ( Size jj = 1; jj <= chain2_bb_atoms.size(); ++jj ) {
          //calc distr squared
          Real d2(chain1_bb_atoms[ii].xyz().distance_squared(chain2_bb_atoms[jj].xyz()));
          ete.atom_pair_energy( chain1_bb_atoms[ii], chain2_bb_atoms[jj], 1, tbemv, d2 );
      }
  }
  core::Real bb_energy (rep_wt * tbemv[core::scoring::fa_rep] + atr_wt * tbemv[core::scoring::fa_atr] );

  // begin optional  ie skip if not needed
  core::Real all_energy;
  if( pose.conformation().num_chains() >= 3){
      core::scoring::EMapVector tbemv_all;
      core::scoring::etable::TableLookupEtableEnergy ete_all( et, emo );
      for ( Size ii = 1; ii <= chain1_bb_atoms.size(); ++ii ){
          for ( Size jj = 1; jj <= all_bb_atoms.size(); ++jj ) {
              //calc distr squared
              Real d2_all(chain1_bb_atoms[ii].xyz().distance_squared(all_bb_atoms[jj].xyz() ));
              ete.atom_pair_energy( chain1_bb_atoms[ii], all_bb_atoms[jj], 1, tbemv_all, d2_all );
          }
      }
      all_energy = (rep_wt * tbemv_all[core::scoring::fa_rep] + atr_wt * tbemv_all[core::scoring::fa_atr] );
  } //end optional for many chains
  else{
      all_energy = bb_energy ;

  }
  TR<< "Number of chains2: " <<pose.conformation().num_chains()
		    <<"    Backbone-backbone score: " << all_energy << std::endl;
  return all_energy;
  //return bb_energy;
}//end bb_score

bool InterfaceStrandFinderMover::is_strand_exposed(core::pose::Pose pose, Size startRes, Size endRes){
  Size exposed_residues(0);
  Size radius_squared = exposure_radius_ * exposure_radius_;

  //iterate the found strand by every other residue (so as to only check residues on the interface-facing
  //side of the strand). If at least half of these residues are "exposed" then the strand itself is exposed.
  for(Size i=startRes; i<=endRes; i+=2){
      numeric::xyzVector<core::Real> C_xyz(pose.residue(i).atom("C").xyz());
      numeric::xyzVector<core::Real> O_xyz(pose.residue(i).atom("O").xyz());
      numeric::xyzVector<core::Real> CO_vector = O_xyz-C_xyz;
      CO_vector = CO_vector.normalize(vector_length_);
      numeric::xyzVector<core::Real> exposure_vector = O_xyz+CO_vector;

      bool residueExposed(true);

      //check to see if any atoms are within distance of this exposure point
      for(Size j=1; j<=pose.total_residue(); j++){
          if(pose.residue(j).atom("C").xyz().distance_squared(exposure_vector)<=radius_squared ||
              pose.residue(j).atom("N").xyz().distance_squared(exposure_vector)<=radius_squared ||
              pose.residue(j).atom("CA").xyz().distance_squared(exposure_vector)<=radius_squared ||
              pose.residue(j).atom("O").xyz().distance_squared(exposure_vector)<=radius_squared){
              residueExposed=false;
              break;
          }
      }
      if(residueExposed){
          exposed_residues++;
      }
  }
  return exposed_residues>=beta_length_/2;
}

std::pair<Pose, Size> InterfaceStrandFinderMover::dock_strands(Pose known_complex, Size known_start, Pose pose, Size found_start){
  id::AtomID_Map< core::id::AtomID > atom_map;
  //some initialization for each new run
  // maps every atomid to bogus atom
  atom_map.clear();
  core::pose::initialize_atomid_map( atom_map, pose, id::BOGUS_ATOM_ID );

  cout << "docking strands: " << known_start << " " << found_start << endl;

  //superimpose the known strand onto the complex
  for(Size j=0; j<beta_length_; j++){
      //Sloppy way of adding all bb atoms, must be a better way...
      core::id::AtomID const id1( pose.residue(found_start + j).atom_index("CA"), found_start + j);
      core::id::AtomID const id2( known_complex.residue(known_start + j).atom_index("CA"), known_start+j );
      atom_map[ id1 ] = id2;

      core::id::AtomID const id3( pose.residue(found_start + j).atom_index("C"),  found_start+j );
      core::id::AtomID const id4( known_complex.residue(known_start + j).atom_index("C"),  known_start+j );
      atom_map[ id3 ] = id4;

      core::id::AtomID const id5( pose.residue(found_start + j).atom_index("N"),  found_start+j );
      core::id::AtomID const id6( known_complex.residue(known_start + j).atom_index("N"),  known_start+j );
      atom_map[ id5 ] = id6;

      core::id::AtomID const id7( pose.residue(found_start + j).atom_index("O"),  found_start+j );
      core::id::AtomID const id8( known_complex.residue(known_start + j).atom_index("O"),  known_start+j );
      atom_map[ id7 ] = id8;
  }
  //do superimpose
  scoring::superimpose_pose(pose, known_complex, atom_map);

  //combine pdbs
  Size removal_chain = known_complex.chain(known_start);
  Size removal_chain_start = known_complex.conformation().chain_begin(removal_chain);
  Size removal_chain_end = known_complex.conformation().chain_end(removal_chain);
  known_complex.conformation().delete_residue_range_slow(removal_chain_start, removal_chain_end);

  known_complex.append_residue_by_jump(pose.residue(1), known_complex.total_residue(), "", "", true /*start new chain*/);
  Size new_strand(known_complex.chain(known_complex.total_residue()));
  for(Size i=2; i<=pose.total_residue(); i++){
      known_complex.append_residue_by_bond(pose.residue(i));
  }

  return make_pair(known_complex, new_strand);
}

//begin mover apply
void InterfaceStrandFinderMover::apply (pose::Pose& pose ) {

  //////////////////////////////
  // handle inputs
  //////////////////////////////

  // get the job info
  //  protocols::jd2::JobOP const job_me( JobDistributor::get_instance()->current_job() );
  //	std::string job_name (JobDistributor::get_instance()->job_outputter()->output_name( job_me ) );

  utility::file::FileName filename (pose.pdb_info()->name());
  cout << "working on file: " << filename << "\n";

  //import the list of known poses
  utility::vector1<std::string> strand_filenames;
  utility::vector1<std::string> complex_filenames;

  //utility::vector1<core::pose::Pose> temp_strands;
  utility::vector1<core::pose::Pose> known_complexes;
  utility::vector1< std::pair<Pose, Size> > known_strands;

  //ensure both known lists are given, and that each given file exists
  if ( option[InterfaceStrandFinder::is_known_strands].user() && option[InterfaceStrandFinder::is_complexes].user()) {
      utility::vector1< std::string > strand_files ( option[ InterfaceStrandFinder::is_known_strands ].value() );
      utility::vector1< std::string > complex_files ( option[ InterfaceStrandFinder::is_complexes ].value() );

      //two lists must be the same size
      if(strand_files.size() != complex_files.size()){
          utility_exit_with_message("is_known_strands and is_complexes must contain the same number of files\n");
      }

      for ( size_t i=1; i<= strand_files.size(); ++i ) {
          utility::io::izstream data( strand_files[i].c_str() );
          if ( !data.good() ) {
              utility_exit_with_message("unable to open known strands file: "+strand_files[i]+"\n");
          }
          while ( data.good() ) {
              std::string line;
              data.getline(line);
              istringstream iss(line);

              if ( data.good() ){
                  //split strand line into tokens (first is the name, second is the index of the first space-pointing carboxyl)
                  vector1<string> tokens;
                  copy(istream_iterator<string>(iss),
                      istream_iterator<string>(),
                      back_inserter<vector1<string> >(tokens));

                  std::string name = tokens[1];
                  core::Size outwardCarboxylIndex = atoi(tokens[2].c_str());
                  strand_filenames.push_back( name );

                  Pose strandPose;
                  core::import_pose::pose_from_file(strandPose, name, false, core::import_pose::PDB_file);
                  std::pair<Pose,Size> strand_pair = make_pair(strandPose, outwardCarboxylIndex);

                  known_strands.push_back(strand_pair);
              }
          }
          data.close();

          utility::io::izstream data2( complex_files[i].c_str() );
          if ( !data2.good() ) {
              utility_exit_with_message("unable to open complexes file: "+complex_files[i]+"\n");
          }

          while ( data2.good() ) {
              std::string name;
              data2.getline(name);
              if ( data2.good() ) complex_filenames.push_back( name );
          }
          data2.close();
      }
  }
  else{
      utility_exit_with_message_status( "need to define -is_known_strands & -is_complexes\n", 1 );
  }
  //known_strands = core::import_pose::poses_from_files(strand_filenames, false/*read_fold_tree?*/, core::import_pose::PDB_file);
  known_complexes = core::import_pose::poses_from_files(complex_filenames, false/*read_fold_tree?*/, core::import_pose::PDB_file);
  //	for (core::Size p = 1; p <= known_complexes.size(); p++){
  //	    std::stringstream temp;
  //	    temp << "TEST_" << p << ".pdb";
  //	    known_complexes[p].dump_pdb(temp.str());
  //	}

  //////////////////////////////////////////////////////
  // find all beta strands in pose
  //////////////////////////////////////////////////////

  //set up dssp info
  core::scoring::dssp::Dssp dssp( pose );
  dssp.insert_ss_into_pose( pose );

  vector1< std::pair< core::Size,core::Size > > strand_endpts;
  for(core::Size i=1; i<=pose.total_residue(); i++){

      //find all the strands in the structure
      core::Size strand_start(0);
      core::Size strand_end(0);
      if(pose.secstruct(i) == 'E'){

          strand_start=i;
          while(i<=pose.total_residue() && pose.secstruct(i)=='E'){
              i++;
          }
          strand_end=i;

          //ensure strand has at least beta_length residues
          if(strand_end-strand_start >= beta_length_+1){
              strand_endpts.push_back(make_pair(strand_start, strand_end));

              //DEBUG
              //cout << "found strand at: " << strand_start << ":" << strand_end << "\n";
          }
      }
  }

  //////////////////////////////////////////////////////
  // check rmsd to known effector-interface strands, check strand exposure, and dock possible matches
  ////////////////////////////////////////////////////////

  Pose testFragment;
  Pose knownFragment;

  //iterate through all of the known effector strands
  for(core::Size j=1; j<=known_strands.size(); j++){
      Size complexNum(0);

      //iterate through each of the found exposed strands
      for(core::Size i=1; i<=strand_endpts.size(); i++){
          core::Size strand_start = strand_endpts[i].first;
          core::Size strand_end = strand_endpts[i].second;

          //compare each beta_length_ fragment from the current found/exposed strand and current known strand
          for(core::Size k=strand_start; k<=(strand_end+1)-beta_length_; k++){
              testFragment.clear();
              for (core::Size resInc = 0 ; resInc < beta_length_; resInc++){
                  testFragment.append_residue_by_bond(pose.residue(k + resInc));
              }
              //We iterate by two here so that the starting residue is always a residue with an interface-facing carboxyl group. This
              //is only true because each of the known_strands must start with an interface-facing carboxyl residue
              for(core::Size l=known_strands[j].second; l<=(known_strands[j].first.total_residue()+1)-beta_length_; l+=2){
                  knownFragment.clear();
                  for (core::Size resInc = 0 ; resInc < beta_length_; resInc++){
                      knownFragment.append_residue_by_bond(known_strands[j].first.residue(l + resInc));
                  }
                  core::Real rms = core::scoring::bb_rmsd(testFragment, knownFragment);

                  bool exposed = is_strand_exposed(pose, k, k+(beta_length_-1));

                  //DEBUG
                  //cout << "strand(" << strand_start << ":" << strand_end << ") fragment(" << k << ":" << k+(beta_length_-1) << ") known( " << l << ":" << l+(beta_length_-1) << ") - RMSD: " << rms << "- Exposed: " << exposed << "\n";

                  if(rms < max_RMSD_ && exposed){

                      //iterate comlexNum for naming purposes
                      complexNum++;

                      //position the potential interaction by bringing together the gtpase and new possible effector.
                      //This is done by superimposing the atoms in the found strant with the atoms in the known strand of the
                      //gtpase-effector pdb for which the rmsd requirement was met.

                      Size complex_strand_start = find_known_strand_start(known_complexes[j], known_strands[j].first);

                      std::pair<Pose, Size> pose_pair = dock_strands(known_complexes[j], (complex_strand_start+l)-1, pose, k);
                      Pose combined_pose = pose_pair.first;
                      Size aligned_chain_num = pose_pair.second;

                      Real bb_clash_score(bb_score(combined_pose, aligned_chain_num, scorefxn_));

                      if(bb_clash_score < maxE_){
                          utility::file::FileName filename2 (known_strands[j].first.pdb_info()->name());
                          cout << filename.base() << "(" << k << "-" << k+beta_length_ << ")" << " --- " <<
                              filename2.base() << "(" << l << "-" << l+beta_length_ << "):\n\trms: " << rms << "\n\tclash score: " << bb_clash_score << "\n";

                          string outputName = filename.base() + "_" + filename2.base() + "_";
                          stringstream out;
                          out << complexNum;
                          outputName = outputName + out.str() + ".pdb";
                          cout << "outputting file: " << outputName << "\n";
                          combined_pose.dump_pdb(outputName);
                      }

                  }
              }
          }
      }
  }

}//end apply


// run protocol
int
main( int argc, char * argv [] )
{

	try {

  using namespace protocols::jd2;
  using namespace protocols::moves;
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  option.add( InterfaceStrandFinder::is_known_strands, "File containing list of pdbs containing known GTPase effector interface strands");
  option.add( InterfaceStrandFinder::is_complexes, "File containing list of pdbs containing GTPase-effector complexes (same effectors from list of known strands)");
  option.add( InterfaceStrandFinder::is_beta_length, "The min length of a beta sheet to count as exposed." );
  option.add( InterfaceStrandFinder::is_max_RMSD, "Max RMSD to allow between known strands and searched strands" );
  option.add( InterfaceStrandFinder::is_vector_length, "length of vector to draw from carboxyl oxyen for strand exposure tests" );
  option.add( InterfaceStrandFinder::is_exposure_radius, "radius of sphere at the end of the exposure vector to search for atoms" );
  option.add( InterfaceStrandFinder::is_maxE, "radius of sphere at the end of the exposure vector to search for atoms" );

  // initialize core
  devel::init(argc, argv);

  protocols::jd2::JobDistributor::get_instance()->go( new InterfaceStrandFinderMover() );

  std::cout << "Done! -------------------------------\n";

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}


