/// @file
/// @brief

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <devel/init.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <protocols/moves/Mover.hh>

//options
#include <basic/basic.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cluster.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

/// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>


#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/conformation/Residue.hh>

#include <protocols/scoring/Interface.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//docking flags
#include <protocols/docking/util.hh>
#include <protocols/docking/types.hh>
#include <protocols/docking/DockMCMProtocol.hh>
#include <protocols/docking/DockingHighResLegacy.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/ChemicalManager.hh>
#include <protocols/docking/DockFilters.hh>

//relax
#include <protocols/relax/FastRelax.hh>
#include <protocols/relax/cst_util.hh>

//Constraints
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

//boost
#include <boost/lexical_cast.hpp>

//docking metrics
#include <protocols/docking/metrics.hh>

//rms
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/rms_util.hh>

//sequence alignment
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/id/SequenceMapping.hh>

//secondary structures
#include <core/scoring/dssp/Dssp.hh>

//#include <protocols/medal/MedalMain.hh>

using basic::T;
using basic::Error;
using basic::Warning;
        
using namespace ObjexxFCL;
using namespace core;
using namespace conformation;
//using namespace scoring;
//using namespace protocols;

///////////////////////////////////////////////////////////////////////////////
class complex_interface_optimize : public protocols::moves::Mover {
public:
	complex_interface_optimize(){}
	void apply( pose::Pose & pose) {

	//check if pose has two chains only
	if (pose.conformation().num_chains()!=2) {
		utility_exit_with_message( "expect pose to have TWO chains only, this one has " + pose.conformation().num_chains());
	}

        //set up foldtree, which one?
        //The one from hybridize
        //Chris Miles' star fold-tree, medal mover
        //detect interface residues/secondary structures
        //figure out interface (does it work on centroid models? )


	//core::util::switch_to_residue_type_set( pose, core::chemical::FA_STANDARD);

        //sequence alignment
        core::sequence::SequenceOP seq1 ( new core::sequence::Sequence(pose) );
        core::sequence::SequenceOP seq2 ( new core::sequence::Sequence(pose) );
        core::sequence::NWAligner nw_aligner;
	utility::file::FileName blosum62( basic::database::full_name("sequence/substitution_matrix/BLOSUM62"));
        core::sequence::ScoringSchemeOP blosum_score( new core::sequence::MatrixScoringScheme( -11, -1, blosum62 ) );
        core::sequence::SequenceAlignment global_align = nw_aligner.align( seq1, seq2, blosum_score ) ;        
	std::cout << global_align;
        std::cout << std::endl;

        core::scoring::ScoreFunctionOP tmpscorefxn = core::scoring::ScoreFunctionFactory::create_score_function(core::scoring::TALARIS_2013);
        (*tmpscorefxn)(pose);
        protocols::scoring::Interface interface( 1 );
        interface.distance( 8 );
        interface.calculate( pose );

        ObjexxFCL::FArray1D_bool is_interface ( (pose).total_residue(), false );
        std::string strseq2=(*global_align.sequence(2)).sequence();
        for ( Size i=1; i<= (pose).total_residue(); ++i ) {
                if (interface.is_interface(i))  {
                        is_interface(i) = true;
                }
        }

        //print out interface
        std::cout << "interface" << std::endl;  
        for ( Size i=1; i<= (pose).total_residue(); ++i ) {
                std::cout << is_interface(i);
        }        
	std::cout<< std::endl;


        //identify secondary structures
        ObjexxFCL::FArray1D_char ssPose ( (pose).total_residue() );
        core::scoring::dssp::Dssp dssp_obj( pose );
        dssp_obj.dssp_reduced(ssPose);
        std::cout << "ss structure" << std::endl;
        for ( Size i=1; i<= (pose).total_residue(); ++i ) {                
		std::cout << ssPose(i);
        }
        std::cout<< std::endl;

	}//end of apply

	virtual std::string get_name() const {
		return "complex_interface_optimize";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	//seq->add_mover(protocols::medal::Medal_main); 
	seq->add_mover( new complex_interface_optimize() );

	try{
		protocols::jd2::JobDistributor::get_instance()->go( seq );
	} catch ( utility::excn::EXCN_Base& excn ) {
		std::cerr << "Exception: " << std::endl;
		excn.show( std::cerr );
	}

	return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {

	// initialize option and random number system
	devel::init( argc, argv );

	protocols::viewer::viewer_main( my_main );
}
