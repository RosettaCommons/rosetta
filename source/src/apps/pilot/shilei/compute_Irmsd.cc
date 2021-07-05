/// @file
/// @brief

#include <protocols/viewer/viewers.hh>
#include <devel/init.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/moves/Mover.hh>

//options
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

/// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/stream_util.hh>
#include <utility/tools/make_vector1.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/scoring/Interface.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

//docking flags
#include <protocols/docking/types.hh>


//relax

//Constraints

//boost

//docking metrics
#include <protocols/docking/metrics.hh>

//rms
#include <core/scoring/rms_util.hh>

//sequence alignment
#include <core/sequence/NWAligner.hh>
#include <core/sequence/SequenceAlignment.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/ScoringScheme.fwd.hh>
#include <core/sequence/MatrixScoringScheme.hh>
#include <core/id/SequenceMapping.hh>


using basic::Error;
using basic::Warning;

using namespace ObjexxFCL;
using namespace core;
using namespace conformation;
//using namespace scoring;
//using namespace protocols;
using utility::operator <<;

///////////////////////////////////////////////////////////////////////////////

class compute_Irmsd : public protocols::moves::Mover {
public:
	compute_Irmsd()= default;
	void apply( pose::Pose & pose) override {

		core::pose::PoseOP native_pose( new core::pose::Pose() );

		//set native pose if asked for
		if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
			core::import_pose::pose_from_file( *native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
			set_native_pose( native_pose );
		} else {
			throw ( CREATE_EXCEPTION(utility::excn::BadInput, "native expected for this app") );
			set_native_pose(nullptr);
		}

		//std::string partners_="_";
		protocols::docking::DockJumps movable_jumps_;
		movable_jumps_.push_back( 1 );

		//perform sequence alignment
		core::sequence::SequenceOP seq1( new core::sequence::Sequence(*native_pose) );
		core::sequence::SequenceOP seq2( new core::sequence::Sequence(pose) );
		//std::cout<<"seq1: "<<*seq1<<std::endl;
		//std::cout<<"seq2: "<<*seq2<<std::endl;

		//utility::io::izstream stream;
		utility::file::FileName blosum62( basic::database::full_name("sequence/substitution_matrix/BLOSUM62"));
		//utility::file::FileName blosum62( "/work/shilei/fresh/rosetta/rosetta_database/sequence/substitution_matrix/BLOSUM62" );
		//basic::database::open( stream, "sequence/substitution_matrix/BLOSUM62" );
		//core::sequence::MatrixScoringScheme matrix_blosum_score;
		//matrix_blosum_score.read_data(stream);
		//core::sequence::ScoringSchemeOP blosum_score( matrix_blosum_score );
		core::sequence::ScoringSchemeOP blosum_score( new core::sequence::MatrixScoringScheme( -11, -1, blosum62 ) );

		core::sequence::NWAligner nw_aligner;
		core::sequence::SequenceAlignment global_align = nw_aligner.align( seq1, seq2, blosum_score ) ;
		std::cout << global_align;
		std::cout << std::endl;

		core::id::SequenceMapping mapping = global_align.sequence_mapping(1, 2);
		//std::cout << mapping;
		//std::cout << std::endl;

		//correspondence map for rms analysis
		std::map<core::Size, core::Size> correspondence;
		std::map<core::Size, core::Size> correspondence_interface;

		//Figure out the interface at the native structure
		core::scoring::ScoreFunctionCOP scorefxn (core::scoring::get_score_function());
		(*scorefxn)(*native_pose);
		//Assume binary complex with jump at 1
		protocols::scoring::Interface interface( 1 );
		interface.distance( 8.0 );
		interface.calculate( *native_pose );
		ObjexxFCL::FArray1D_bool is_interface ( (*native_pose).size(), false );
		for ( Size i=1; i<= (*native_pose).size(); ++i ) {
			if ( interface.is_interface(i) ) {
				is_interface(i) = true;
			}
		}

		//set up correspondence map for alignment
		for ( core::Size i = 1; i <= (*seq1).length(); ++i ) {
			if ( mapping[i] ) {
				correspondence[i] = mapping[i];
				if ( is_interface(i) ) {
					correspondence_interface[i]=mapping[i];
				}
			}
		}

		std::cout << correspondence;
		std::cout << std::endl;

		std::cout << correspondence_interface;
		std::cout << std::endl;


		core::Real Isc=protocols::docking::calc_interaction_energy(pose,scorefxn,utility::tools::make_vector1<core::Size>(1));

		//std::cout << "CA_rms: " << core::scoring::CA_rmsd( *native_pose, pose, correspondence) << " Irms_CA: " << core::scoring::CA_rmsd( *native_pose, pose, correspondence_interface) << std::endl;
		protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );
		job->add_string_real_pair("rms_CA", core::scoring::CA_rmsd( *native_pose, pose, correspondence));
		job->add_string_real_pair("Irms_CA", core::scoring::CA_rmsd( *native_pose, pose, correspondence_interface));
		job->add_string_real_pair("Isc", Isc);

	}//end of apply

	std::string get_name() const override {
		return "compute_Irmsd";
	}
};

///////////////////////////////////////////////////////////////////////////////

void*
my_main( void* ) {
	using namespace protocols::moves;

	SequenceMoverOP seq( new SequenceMover() );
	seq->add_mover( utility::pointer::make_shared< compute_Irmsd >() );

	protocols::jd2::JobDistributor::get_instance()->go( seq );

	return nullptr;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] ) {

	try{
		// initialize option and random number system
		devel::init( argc, argv );

		protocols::viewer::viewer_main( my_main );
	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}
	return 0;

}
