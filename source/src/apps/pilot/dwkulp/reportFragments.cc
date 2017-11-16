
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
#include <core/io/pdb/pdb_writer.hh>
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

#include <core/fragment/FragData.hh>
#include <core/fragment/ConstantLengthFragSet.hh>

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
using core::import_pose::pose_from_file;

typedef utility::pointer::owning_ptr< AnnotatedFragData > AnnotatedFragDataOP;
typedef utility::pointer::owning_ptr< AnnotatedFragData const > AnnotatedFragDataCOP;

void get_fragment_sets(core::pose::Pose &_pose, string _polyAminoAcidSequence,   utility::vector1<core::fragment::ConstantLengthFragSetOP> &fragSets);

static basic::Tracer TR( "reportFragments" );


int main( int argc, char * argv[] ) {
	try {
		// init option system
		devel::init(argc,argv);

		// Read in a pdb (assume pre-relaxed structure)
		core::pose::Pose pose;
		core::import_pose::pose_from_file(pose,basic::options::start_file(), core::import_pose::PDB_file);

		// poly-Val or poly-Ala sequence for testing fragment picker while designing.
		string forcePolyAAsequence = "";
		if ( basic::options::option[basic::options::OptionKeys::dwkulp::forcePolyAAfragments].user() ) {
			forcePolyAAsequence = basic::options::option[basic::options::OptionKeys::dwkulp::forcePolyAAfragments]();
		}

		// Get sequence and residues from loop definition
		utility::vector1<core::fragment::ConstantLengthFragSetOP> fragSets;
		get_fragment_sets(pose,forcePolyAAsequence,fragSets);

		ofstream fout("data.txt");
		core::fragment::FrameIterator fIt;

		cout << "Positions affected by fragments: "<<fragSets[1]->min_pos()<<" "<<fragSets[1]->max_pos()<<" max frag length: "<<fragSets[1]->max_frag_length()<<endl;
		//for (fIt = fragSets[1]->begin(); fIt != fragSets[1]->end();++fIt){
		FrameList frames;
		fragSets[1]->region_simple(fragSets[1]->min_pos(),fragSets[1]->max_pos(),frames);
		//fragSets[1]->region_simple(1,fragSets[1]->max_pos() - fragSets[1]->min_pos() ,frames);
		cout << "Frames Size: "<<frames.size()<<endl;
		for ( unsigned a = 1; a<= frames.size(); a++ ) {
			FrameOP f = frames[a];
			f->show(std::cout);
			for ( unsigned int i = 1; i <= f->nr_frags(); i++ ) {
				FragDataOP fd = f->fragment_ptr(i);
				//      //fd->apply(pose,*(*fIt));
				//
				//      // Apply whole fragdata using apply(pose,frame)
				//      for (unsigned int j = 1; j <= fd->size();j++){
				// SingleResidueFragDataCOP res = fd->get_residue(j);
				//
				// res->apply(pose,j);
				// char ss  = res->secstruct();
				// char seq = res->sequence();
				//
				// Real phi = pose.torsion(id::TorsionID( j, id::BB, 1 ));
				// Real psi = pose.torsion(id::TorsionID( j, id::BB, 2 ));
				//
				// std::cout << a<<" "<<i<<" "<<j<<" "<<ss<<" "<<seq<<" "<<phi<<" "<<psi<<std::endl;
				//      }
			}
			a++;
		}

		fout.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}


void get_fragment_sets(core::pose::Pose &_pose, string _polyAminoAcidSequence,   utility::vector1<core::fragment::ConstantLengthFragSetOP> &fragSets){

	string sequence;
	string secstruct;
	for ( unsigned int j = 1; j <=_pose.size(); j++ ) {

		if ( _polyAminoAcidSequence == "" ) {
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
	if ( sequence.size() > 9 ) {
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
		//10,    // collector must know the size of query
		sequence.size(),    // collector must know the size of query
		pickIt->n_candidates_, // how many candidates to collect
		comparator,    // yes, here the comparator comes to sort fragments within the collector
		pickIt->get_score_manager()->count_components());
	pickIt->set_candidates_collector(3, collector3);

	if ( sequence.size() > 9 ) {
		CandidatesCollectorOP collector9 = new BoundedCollector<CompareTotalScore> (
			sequence.size(),    // collector must know the size of query
			pickIt->n_candidates_, // how many candidates to collect
			comparator,    // yes, here the comparator comes to sort fragments within the collector
			pickIt->get_score_manager()->count_components());

		pickIt->set_candidates_collector(9, collector9);
	}


	// Setup a fragment selector
	pickIt->selector_ = new BestTotalScoreSelector(pickIt->n_frags_, pickIt->get_score_manager());


	// Simplest selecting fragments and writing a fragment database file out.... we will have to modify this function later to avoid "save_fragments" call.
	pickIt->bounded_protocol();

	cout << "Get Fragments for "<<_pose.pdb_info()->number(1)<<endl;
	fragSets = pickIt->getFragSet(_pose.pdb_info()->number(1));
}


