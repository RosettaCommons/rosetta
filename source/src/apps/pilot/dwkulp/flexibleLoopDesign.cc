
#include <basic/options/util.hh>
#include <devel/init.hh>
#include <basic/options/option.hh>

#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/dwkulp.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/conformation/Residue.hh>

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/operation/TaskOperations.fwd.hh>
#include <core/pack/task/operation/NoRepackDisulfides.hh>



#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/MoveMap.fwd.hh>
//#include <protocols/relax/SequenceRelax.hh>
//#include <protocols/relax/ClassicRelax.hh>
//#include <protocols/relax/FastMultiRelax.hh>
//#include <protocols/relax_protocols.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.fwd.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/CCDLoopClosureMover.hh>
#include <protocols/simple_moves/BackboneMover.hh> //23
#include <protocols/moves/MoverContainer.hh> //Sequence Mover
#include <protocols/loops/kinematic_closure/KinematicMover.hh>
#include <protocols/moves/kinematic_closure/KinematicPerturber.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <numeric/random/random.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/rms_util.hh>


#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/FragmentPicker.fwd.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallProvider.fwd.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentCandidate.fwd.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>
#include <protocols/frag_picker/DiversifyCrmsdByClustering.hh>

#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>

#include <core/fragment/ConstantLengthFragSet.hh>

#include <protocols/simple_moves/FragmentMover.hh>
#include <protocols/loops/ccd_closure.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/excn/Exceptions.hh>


#include <numeric/random/random.hh>
#include <core/import_pose/import_pose.hh>
#include <iostream>
#include <fstream>
#include <numeric>
#include <functional>


using namespace std;
using namespace core;
using namespace core::sequence;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace protocols::loops;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace core::chemical;
using core::import_pose::pose_from_pdb;

double annealTemperature(double initialTemp, double finalTemp, int step, int totalsteps, int numAnnealCycles);
void setup_fragment_mover(core::pose::Pose &_pose,    protocols::loops::Loop &_loop,   vector<vector<protocols::moves::MoverOP> > &_bbMovers, string _polyAminoAcidOverride);
unsigned int hd(const std::string& s1, const std::string& s2); // hamming distance
static basic::Tracer TR("flexibleLoopDesign");
static numeric::random::RandomGenerator RG(222578262);



int main( int argc, char * argv[] ) {
    try {
  // init option system
  devel::init(argc,argv);

  // Read in a pdb (assume pre-relaxed structure)
  core::pose::Pose pose;
  core::import_pose::pose_from_pdb(pose,basic::options::start_file());

  cout << "Get loops from file.."<<endl;
  // Get loops
  protocols::loops::Loops my_loops = protocols::loops::get_loops_from_file();

  cout << "Done getting loops"<<endl;
  // Get CCD closure movers (one per loop)
  vector<protocols::loops::CCDLoopClosureMoverOP> ccds;
  //vector<protocols::simple_moves::BackboneMoverOP> bbMovers;
  vector<vector<protocols::moves::MoverOP> > bbMovers;

  // poly-Val or poly-Ala sequence for testing fragment picker while designing.
  string forcePolyAAsequence = "";
  if (basic::options::option[basic::options::OptionKeys::dwkulp::forcePolyAAfragments].user()){
    forcePolyAAsequence = basic::options::option[basic::options::OptionKeys::dwkulp::forcePolyAAfragments]();
  }

  map<int,bool> neighboringResidues;
  map<int,bool> loopResidues;
  for (core::uint i = 1 ; i <= my_loops.num_loop();i++){
      cout << "Loop: "<<i<<endl;

    protocols::loops::Loop &l = my_loops[i];

    core::kinematics::MoveMapOP moveMap = new core::kinematics::MoveMap();
    moveMap->set_bb(false);
    moveMap->set_chi(false);

    for (core::uint j = l.start(); j <=l.stop();j++){
      cout << "Residue "<<j<<" is being marked as mobile\n";
      moveMap->set_bb(j,true);
      moveMap->set_chi(j,true);
      loopResidues[j] = true;


      // Only search for loop-neighbors if no resfile has been specified
      if (!basic::options::option[basic::options::OptionKeys::packing::resfile].user()){

	// Find neighboring residues
	for (core::Size r = 1;r <= pose.n_residue();r++){
	  if (r >= l.start() && r <= l.stop()) continue; //skip loop residues

	  // Distance from this residue (r) to the current loop residue (j)
	  double dist = pose.residue(r).xyz("CA").distance_squared(pose.residue(j).xyz("CA"));

	  // This residue is a neighbor!
	  if (dist < 64){
	    neighboringResidues[r] = true;
	  }

	}
      }



    }

    protocols::loops::CCDLoopClosureMoverOP closer = new protocols::loops::CCDLoopClosureMover(l , moveMap );
    closer->set_bRama_check(true);

    // Keep list of CCD Closure Movers
    ccds.push_back(closer);
    /*
    protocols::simple_moves::BackboneMoverOP small = new protocols::simple_moves::SmallMover(moveMap, 10, 200); //huge moves for sampling
    small->angle_max( 'H', 180.0 );
    small->angle_max( 'E', 180.0 );
    small->angle_max( 'L', 180.0 );

    bbMovers.push_back(small);
    */

    // Get sequence and residues from loop definition
    setup_fragment_mover(pose,l,bbMovers,forcePolyAAsequence);
  }

  cout << "Fold tree now.."<<endl;
  // Define the fold tree from the loops..
  core::kinematics::FoldTree ft;
  protocols::loops::fold_tree_from_loops(pose,my_loops,ft,false);
  pose.fold_tree(ft);
  cout <<  pose.fold_tree() <<endl;

  // Define moves

  // Get a scoring function
  core::scoring::ScoreFunctionOP scorefxn;
  scorefxn = core::scoring::get_score_function();

  // A Repacking Mover
  core::pack::task::TaskFactoryOP task = new core::pack::task::TaskFactory;

  task->push_back( new core::pack::task::operation::NoRepackDisulfides);
  //task->push_back( new core::pack::task::operation::RestrictToRepacking );
  task->push_back( new core::pack::task::operation::InitializeFromCommandline );

  basic::options::option[basic::options::OptionKeys::packing::use_input_sc](true);
  basic::options::option[basic::options::OptionKeys::packing::ex1::ex1](true);
  basic::options::option[basic::options::OptionKeys::packing::ex2::ex2](true);
  //basic::options::option[basic::options::OptionKeys::packing::ex3::ex3](true);
  //basic::options::option[basic::options::OptionKeys::packing::ex4::ex4](true);
  basic::options::option[basic::options::OptionKeys::packing::extrachi_cutoff](1);  // want to relax previously buried side-chains, that are now exposed in unbound state

  if (basic::options::option[basic::options::OptionKeys::packing::resfile].user()){
    cout << "READ RESFILE PLEASE!!!"<<endl;
    task->push_back( new core::pack::task::operation::ReadResfile );
  }

  core::pack::task::PackerTaskOP ptask = task->create_task_and_apply_taskoperations( pose );
  //core::pack::task::PackerTaskOP ptask = core::pack::task::TaskFactory::create_packer_task( pose );

  cout << "Num to be repacked11: "<<ptask->num_to_be_packed()<<endl;


  // If no resfile, make loops repack, loop-neighbors design and fix everything else.
  if (!basic::options::option[basic::options::OptionKeys::packing::resfile].user()){

    core::pack::task::operation::RestrictResidueToRepackingOP restrict_to_repack_taskop( new core::pack::task::operation::RestrictResidueToRepacking() );
    core::pack::task::operation::PreventRepackingOP prevent_repack_taskop( new core::pack::task::operation::PreventRepacking() );
    core::pack::task::operation::RestrictAbsentCanonicalAASOP design_taskop(new core::pack::task::operation::RestrictAbsentCanonicalAAS() );

    utility::vector1< bool > allowed_aas( chemical::num_canonical_aas, true );
    design_taskop->keep_aas(allowed_aas);
    for ( core::Size i=1; i <= pose.total_residue(); ++i ) {

      // A loop-neighbor residue?
      if (neighboringResidues[i]){
	cout << "Residue "<<i<<" is a neighbor."<<endl;
	ptask->nonconst_residue_task( i ).restrict_absent_canonical_aas( allowed_aas );
	design_taskop->include_residue(i);
	continue;
      }

      // A loop residue ?
      if (loopResidues[i]){
	cout << "Residue "<<i<<" is a loop."<<endl;
	ptask->nonconst_residue_task( i ).restrict_to_repacking();
	restrict_to_repack_taskop->include_residue(i);
	continue;
      }

      // Default is no repacking, no nothin..
      ptask->nonconst_residue_task( i ).prevent_repacking();
      prevent_repack_taskop->include_residue(i);
    }

    // Add TaskOperators to task
    task->push_back(design_taskop);
    task->push_back(restrict_to_repack_taskop);
    task->push_back(prevent_repack_taskop);
  } //else {
//	ptask = task->create_task_and_apply_taskoperations( pose );
//  }

  cout << "Num to be repacked22: "<<ptask->num_to_be_packed()<<endl;


  protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover();
  pack_mover->task_factory( task );
  pack_mover->score_function( scorefxn );

  // Define a design-movemap
  core::kinematics::MoveMapOP designMoveMap = new core::kinematics::MoveMap();
  designMoveMap->set_bb(false);
  designMoveMap->set_chi(false);


  int numTotalDesignablePositions = 0;
 for ( core::Size i=1; i <= pose.total_residue(); ++i ) {
   if (ptask->being_designed( i ) ||
       ptask->being_packed( i )){
     designMoveMap->set_bb(i,true);
     designMoveMap->set_chi(i,true);
   }

   if (ptask->being_designed( i )) {
       cout << "Position "<<i <<" is being designed."<<endl;
     numTotalDesignablePositions++;
   }


 }
  cout << "Num to be designed22: "<<numTotalDesignablePositions<<endl;

  // Define a sequence of moves
  //protocols::moves::SequenceMoverOP seqmover = new protocols::moves::SequenceMover();
  //for (core::uint i = 0; i < ccds.size();i++){
      //seqmover->add_mover(bbMovers[i]);
      //seqmover->add_mover(ccds[i]);

  //cout << "Sequence mover add bbMover and CCD "<<i<<"\n";
  //}

  // Add 5 pack movers. This allows 5 mutations per structure predicted loop
  //seqmover->add_mover(pack_mover);
  //seqmover->add_mover(pack_mover);
  //seqmover->add_mover(pack_mover);
  //seqmover->add_mover(pack_mover);
  //seqmover->add_mover(pack_mover);


  core::optimization::AtomTreeMinimizerOP minimizer;
  core::optimization::MinimizerOptions options( "dfpmin", 0.001, true /*use_nblist*/, false /*deriv_check*/ );
  minimizer = new core::optimization::AtomTreeMinimizer;


  // Monte Carlo loop
  core::Real temperature = 100;
  protocols::moves::MonteCarlo mc( pose, *scorefxn, temperature);
  int numMCSteps = 100000;

  core::Real wtScorePreMin = scorefxn->score(pose);
  pose.dump_pdb("wtPre.pdb");

  // Repack and minimize the wild type protein
  pack_mover->apply(pose);
  minimizer->run( pose, *designMoveMap, *scorefxn, options );

  core::Real wtScorePostMin = scorefxn->score(pose);
  TR << "WT Score: "<<wtScorePreMin<<" "<<wtScorePostMin<<std::endl;
  pose.dump_pdb("wtPost.pdb");

  ofstream fout("data.txt");

  core::pose::Pose startingPose = pose;
  string startingSequence = startingPose.sequence();
  for (unsigned int i = 0; i < (unsigned int)numMCSteps;i++){

	  // Anneal temperature
	  mc.set_temperature(annealTemperature(temperature,0.0,i,numMCSteps,10));

	  // Modify system (loop conformation)
	  for (core::uint c = 0; c < ccds.size();c++){
	      int randomLoopSize = RG.random_range(0,bbMovers[c].size()-1);

	      cout << "Fragment insertion on loop "<<c<<" index of fragment size(0=3,1=9): "<<randomLoopSize<<endl;
	      bbMovers[c][randomLoopSize]->apply(pose);
	      ccds[c]->apply(pose);

	      // test for proper closure ... chain energy.
	  }



	  // Design/Repack a large number of positions (default 20 outer_iterations, inner_iterations 6443 with 22 designable positions)
	  pack_mover->apply(pose);

	  // Minimize system
	  minimizer->run( pose, *designMoveMap, *scorefxn, options );

   	  bool accepted = mc.boltzmann( pose );
 	  if (accepted) {
	       core::Real score =  scorefxn->score(pose);

	      (*scorefxn)(pose);
	      scorefxn->show(std::cout,pose);

	      stringstream RMSDs;
	      double averageRMSD;
	      for (core::uint loopIndex = 1 ; loopIndex <= my_loops.num_loop();loopIndex++){
		protocols::loops::Loop &l = my_loops[loopIndex];
		core::Real rmsd = core::scoring::CA_rmsd(startingPose,pose,l.start(),l.stop());

	        cout << "RMSD["<<loopIndex<<"]: "<<rmsd<<endl;
		char tmpRMSD[80];
		sprintf(tmpRMSD, "%8.3f ",rmsd);
		RMSDs << tmpRMSD;

		averageRMSD += rmsd;
	      }

	      // Re-pick fragments based on % change from WT.
	      cout << "compute % change"<<endl;
	      int numDifferentAA = hd(startingPose.sequence(),pose.sequence());

	      double percentSequenceChange  = 100 - ( ((double) (numTotalDesignablePositions - numDifferentAA) / (double) numTotalDesignablePositions) * 100 );
	      cout << "Number of positions designed is "<<numDifferentAA<<" out of "<<numTotalDesignablePositions<<" percentSequenceChange: "<<percentSequenceChange<<endl;
	      char tmp[80];
	      sprintf(tmp, "DATA: accepted.%06d.pdb %06d %8.3f %8.3f %5.0f %8.3f %8.3f %s",
		      i, i,
		      score,
		      score - wtScorePostMin,
		      mc.temperature(),
		      averageRMSD / (float)my_loops.num_loop(),
		      percentSequenceChange,
		      RMSDs.str().c_str());

	      TR << tmp <<endl;

	      fout << tmp <<endl;

	      char tmp2[80];
	      sprintf(tmp2,"accepted.%06d.pdb",i);

	      pose.dump_pdb(tmp2);

	      cout << "shall we change fragments? "<<endl;
	      // After 50 sequence change, we pick more fragments
	      if (forcePolyAAsequence == "" && percentSequenceChange > 50) {

		cout << "New fragments please.."<<endl;
		bbMovers.clear();
		for (core::uint i = 1 ; i <= my_loops.num_loop();i++){
		  cout << "Loop: "<<i<<endl;
		  protocols::loops::Loop &l = my_loops[i];

		  // Get sequence and residues from loop definition
		  setup_fragment_mover(pose,l,bbMovers,forcePolyAAsequence);
		}

/*
		seqmover->clear();
		for (core::uint i = 0; i < ccds.size();i++){
		  seqmover->add_mover(bbMovers[i]);
		  seqmover->add_mover(ccds[i]);

		  cout << "Sequence mover add bbMover and CCD "<<i<<"\n";
		}
*/

	      }
	  }

	  // Recover WT each iteration?
	  cout << "Recover WT"<<endl;
	  pose = startingPose;

  }

  // Recover low and spit out pdb file
  mc.recover_low(pose);
  pose.dump_pdb("winnerMC.pdb");

  fout.close();

    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}



double annealTemperature(double initialTemp, double finalTemp, int step, int totalsteps, int numAnnealCycles){
		/*
		  FREQ = 5; // a parameter to set from outside the object..
		  cycleNumber = floor(double(step) / double(FREQ)) ;
		  expFactor   = double(FREQ) / 5;
		  expVal  =  (step - (FREQ * cycleNumber)) / expFactor;
		  setCurrentTemp(startTemp* exp(-expVal));

		*/
		//
		int numStepsInCycle = int(totalsteps / numAnnealCycles);
		int lastStart       = (int(step / numStepsInCycle) * numStepsInCycle)+1;
		if (step % numStepsInCycle == 0){
			lastStart = step;
		}
		double zeroingFactor  = double(numStepsInCycle)/10.0 + double(numStepsInCycle)/20.0; // seems to work ok for deciding when we are close to zero

		return (initialTemp * exp( - (step-lastStart) / zeroingFactor));
}


void setup_fragment_mover(core::pose::Pose &_pose,    protocols::loops::Loop &_loop,   vector< vector<protocols::moves::MoverOP> > &_bbMovers, string _polyAminoAcidSequence){

  // Fragment move map says all positions in fragSets can move (or be inserted on)?
  core::kinematics::MoveMapOP moveMapFrags3 = new core::kinematics::MoveMap();
  moveMapFrags3->set_bb(false);
  moveMapFrags3->set_chi(false);

  core::kinematics::MoveMapOP moveMapFrags9 = new core::kinematics::MoveMap();
  moveMapFrags9->set_bb(false);
  moveMapFrags9->set_chi(false);

  string sequence;
  string secstruct;
  for (unsigned int j = _loop.start(); j <=_loop.stop();j++){
    if (j <= _loop.stop()){// - 3){
      cout << "Loop residue "<<j<<" is being set to moveable for FRG3"<<endl;
      moveMapFrags3->set_bb(j,true);
      moveMapFrags3->set_chi(j,true);
    }
    if (j <= _loop.stop()){// - 9){
      cout << "Loop residue "<<j<<" is being set to moveable for FRG9"<<endl;
      moveMapFrags9->set_bb(j,true);
      moveMapFrags9->set_chi(j,true);
    }


    if (_polyAminoAcidSequence == ""){
      sequence += oneletter_code_from_aa(_pose.aa(j));
    } else {
      sequence += _polyAminoAcidSequence;
    }
      secstruct += "L";
  }

  cout << "Sequence is "<<sequence<<" sectstruct is "<<secstruct<<endl;

    // Setup a new FragmentPicker
  FragmentPickerOP pickIt = new FragmentPicker();

  // Create a VallProvider, which reads in a database and creates VallResidues and VallChunks
  VallProviderOP chunks = new VallProvider();
  chunks->vallChunksFromLibrary(basic::options::option[basic::options::OptionKeys::in::file::vall]()[1]);

  // Pass the VallProvider to FragmentPicker
  pickIt->set_vall(chunks);

  // Setup the size fragments
  pickIt->frag_sizes_.push_back(3);
  if (sequence.size() > 9){
    pickIt->frag_sizes_.push_back(9);
  }

  // Setup the number of pre-filtered fragments (candidates)
  pickIt->n_candidates_ = 100;

  // Setup the number of saved fragments
  pickIt->n_frags_ = 50;

  // Set sequence string
  pickIt->set_query_seq(sequence);

  // Set the SS string
  pickIt->add_query_ss(secstruct,"loop");

  // Set the picked position range
  //pickIt->set_picked_positions(1,10); // range within query sequence

  // Prefix for output to a file..., in the end I will remove this..
  pickIt->prefix_ = "fragments_design";

  // Read in the set of FragmentScoreMethods.  Look at FragmentScoreManager::create_scores() if you want to avoid reading in this file and do it by hand.
  //  A reason for doing it by hand is that not all FragmentScoreMethods will work in this "manual" mode, because many rely on a SequenceProfile object
  //  which normally is read in from a file .. in design mode this won't work.  In place of that, we use a ProfileScoreSubMatrix (BLOSUM62)
  pickIt->get_score_manager()->create_scores(basic::options::option[basic::options::OptionKeys::frags::scoring::config](),pickIt);

  // Setup a comparator, which sorts by the total score.
  CompareTotalScore comparator( pickIt->get_score_manager() );


  // 3 residue and 9 residue collector
  CandidatesCollectorOP collector3 = new BoundedCollector<CompareTotalScore> (
      //10, 	  // collector must know the size of query
	         sequence.size(), 	  // collector must know the size of query
		pickIt->n_candidates_, // how many candidates to collect
		comparator,		  // yes, here the comparator comes to sort fragments within the collector
		pickIt->get_score_manager()->count_components());
  pickIt->set_candidates_collector(3, collector3);

  if (sequence.size() > 9){
    CandidatesCollectorOP collector9 = new BoundedCollector<CompareTotalScore> (
	    			                sequence.size(), 	  // collector must know the size of query
						pickIt->n_candidates_, // how many candidates to collect
						comparator,		  // yes, here the comparator comes to sort fragments within the collector
						pickIt->get_score_manager()->count_components());

    pickIt->set_candidates_collector(9, collector9);
  }


   // Setup a fragment selector
   pickIt->selector_ = new BestTotalScoreSelector(pickIt->n_frags_, pickIt->get_score_manager());


  // Simplest selecting fragments and writing a fragment database file out.... we will have to modify this function later to avoid "save_fragments" call.
  pickIt->bounded_protocol();


  utility::vector1<core::fragment::ConstantLengthFragSetOP> fragSets = pickIt->getFragSet(_loop.start());

  // Setup a fragment mover
  vector<protocols::moves::MoverOP> multipleSizeMovers;

  using protocols::simple_moves::ClassicFragmentMover;

  protocols::simple_moves::ClassicFragmentMoverOP fragMover3  = new ClassicFragmentMover(fragSets[1],moveMapFrags3);
  fragMover3->enable_end_bias_check(false);
  multipleSizeMovers.push_back(fragMover3);

  if (sequence.size() > 9){
    protocols::simple_moves::ClassicFragmentMoverOP fragMover9  = new ClassicFragmentMover(fragSets[2],moveMapFrags9);
    fragMover9->enable_end_bias_check(false);
    multipleSizeMovers.push_back(fragMover9);
  }

  _bbMovers.push_back(multipleSizeMovers);
}


unsigned int hd(const std::string& s1, const std::string& s2){

    if (s1.size() != s2.size()){
	throw std::invalid_argument(
          "Strings passed to hd() must have the same length"
	    );
    }

    return std::inner_product(
        s1.begin(), s1.end(), s2.begin(),
        0, std::plus<unsigned int>(),
        std::not2(std::equal_to<std::string::value_type>())
	);
}
