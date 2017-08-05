#include <devel/init.hh>

#include <core/sequence/Sequence.hh>
#include <basic/Tracer.hh>

// option key includes
#include <core/options/option.hh>
#include <core/options/option_macros.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/options/keys/frags.OptionKeys.gen.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/constraints.OptionKeys.gen.hh>

// utility headers
#include <utility/vector1.hh>

#include <core/fragment/picking/FragmentPicker.hh>
#include <core/fragment/picking/VallProvider.hh>
#include <core/fragment/picking/VallChunk.hh>
#include <core/fragment/picking/VallResidue.hh>
#include <core/fragment/picking/FragmentCandidate.hh>
#include <core/fragment/picking/FragmentSelectingRule.hh>
#include <core/fragment/picking/VallChunkFilter.hh>
#include <core/fragment/picking/BoundedCollector.hh>
#include <core/fragment/picking/BestTotalScoreSelector.hh>

#include <core/fragment/picking/scores/FragmentScoringMethod.hh>
#include <utility>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/io/mpistream.hh>


static THREAD_LOCAL basic::Tracer trace("PickSSE");

using namespace core;
using namespace core::fragment;
using namespace core::fragment::picking;
using namespace core::fragment::picking::scores;
using namespace options;
using namespace options::OptionKeys;

class PickSSE {

public:
	/// \brief constructor parses command line and sets up picking
	PickSSE() {

		using namespace core;
		using namespace core::sequence;
		using namespace core::options;
		using namespace core::options::OptionKeys;

		VallProviderOP chunks = new VallProvider();
		chunks->vallChunksFromLibrary(option[in::file::vall]()[1]);
		my_picker.set_vall(chunks);

		SequenceProfileOP q_prof(new SequenceProfile);
		q_prof->read_from_checkpoint(option[in::file::checkpoint]());
		my_picker.set_query_seq(q_prof);

		n_res_ = q_prof->size();

		my_picker.n_candidates_ = option[frags::n_candidates]();
		my_picker.n_frags_ = option[frags::n_frags]();

		FragmentScoreManagerOP scoring = my_picker.get_score_manager();
		CompareTotalScore comparator( my_picker.get_score_manager() );
		CandidatesCollectorOP collector = new BoundedCollector<CompareTotalScore> (
			n_res_,my_picker.n_candidates_,comparator);
	}

private:

	Size n_res_;
	utility::vector1<Size> sse_start;
	utility::vector1<Size> sse_end;
	utility::vector1<char> sse_type;
	utility::vector1<Real> sse_dist;
	FragmentPicker my_picker;
	CompareTotalScore *comparator;

	/// Example format: ALPHA_HELIX    8   21 -16.000   4.000   7.000
	void readSSE(std::strig filename) {

		std::string line;
		while ( getline( input, line ) ) {
			if ( line.substr(0,1) == "#" ) continue;
			istringstream line_stream( line );
			std::string type;
			Size from, to;
			Real x,y,z;
			line_stream >> type >> from >> to >> x >> y >> z;
			if ( type.find("HELIX") !=string::npos ) {
				sse_type.push_back( 'H' );
				sse_start.push_back( from );
				sse_end.push_back( to );
				sse_dist.push_back( sqrt( x*x + y*y + z*z ) );
			}
			if ( type.find("SHEET") !=string::npos ) {
				sse_type.push_back( 'E' );
				sse_start.push_back( from );
				sse_end.push_back( to );
				sse_dist.push_back( sqrt( x*x + y*y + z*z ) );
			}
		}
	}

	void pick_sse() {

		CompareTotalScore comparator( my_picker->get_score_manager() );
		for ( Size fi=1; fi<=sse_type.size(); fi++ ) {
			CandidatesCollectorOP collector = new BoundedCollector<CompareTotalScore> (
				n_res_,my_picker.n_candidates_,comparator);
			my_picker.set_candidates_collector(3, collector);
			my_picker.selector_ = new BestTotalScoreSelector(my_picker.n_frags_, my_picker.get_score_manager() );
		}
	}
};

void register_options() {

	OPT(in::file::native);
	OPT(in::file::checkpoint);
	OPT(frags::scoring::config);
	OPT(constraints::cst_file);
	OPT(in::path::database);
	OPT(frags::allowed_pdb);
	OPT(frags::denied_pdb);
}


int main(int argc, char * argv[]) {
	try {
		register_options();
		devel::init(argc, argv);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;


	// my_picker->set_candidates_collector(3, collector);
	// my_picker->set_candidates_collector(9, collector);
	// my_picker->selector_ = new BestTotalScoreSelector(my_picker->n_frags_, scoring);

	// my_picker->prefix_ = "fragments";

	// my_picker->bounded_protocol();
}
