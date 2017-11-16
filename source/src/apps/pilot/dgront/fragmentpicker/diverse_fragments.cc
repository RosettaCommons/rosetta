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


static basic::Tracer trace("fragmentpicker_integration_demo");

using namespace core;
using namespace core::fragment;
using namespace core::fragment::picking;
using namespace core::fragment::picking::scores;
using namespace options;
using namespace options::OptionKeys;

void register_options() {

	OPT(in::file::native);
	OPT(in::file::psipred_ss2);
	OPT(in::file::checkpoint);
	OPT(frags::scoring::config);
	OPT(constraints::cst_file);
	OPT(in::path::database);
	OPT(frags::denied_pdb);
}

void setup_nnmake_style_quota() {
}

int main(int argc, char * argv[]) {
	try {
		using namespace core;
		using namespace core::sequence;
		using namespace core::options;
		using namespace core::options::OptionKeys;

		register_options();
		devel::init(argc, argv);

		//------------ CREATE A PICKER OBJECT
		// FragmentPickerOP my_picker = new FragmentPicker("PValuedFragmentScoreManager");
		FragmentPickerOP my_picker = new FragmentPicker();

		//------------ PLUG-IN SEQUENCE PROFILE
		SequenceProfileOP q_prof(new SequenceProfile);
		q_prof->read_from_checkpoint(option[in::file::checkpoint]());
		my_picker->set_query_profile(q_prof);

		//------------ LOAD SECONDARY STRUCTURE(s) - it is possible to load more than one
		//------------ TWO PARAMETERS: SS_FILE AND A TAG
		//----------- THE TAG MUST BE CONSISTENT WITH THE ONE GIVEN IN THE SCORE CONFIG FILE!!!!!
		my_picker->read_ss_file(option[in::file::psipred_ss2](),"johny");

		//------------ READ VALL  AND PLUG IT INTO THE PICKER
		VallProviderOP chunks = new VallProvider();
		chunks->vallChunksFromLibrary(option[in::file::vall]()[1]);
		my_picker->set_vall(chunks);

		//------------ FRAGMENT SIZE: WE NEED 9-MERS AND 3-MERS
		my_picker->frag_sizes_.push_back(3);
		my_picker->frag_sizes_.push_back(9);

		//------------ HOW MANY CANDIDATES, HOW MANY FRAGMENTS
		my_picker->n_candidates_ = option[frags::n_candidates]();
		my_picker->n_frags_ = option[frags::n_frags]();
		trace.Info << "Picking " << my_picker->n_frags_ << " fragments based on "<<
			my_picker->n_candidates_ << " candidates" << std::endl;

		//----------- SETUP SCORING SYSTEM
		FragmentScoreManagerOP scoring = my_picker->get_score_manager();
		scoring->create_scores(option[frags::scoring::config](), my_picker);

		//----------- SETUP COLLECTOR (for candidates) AND SELECTOR (for final fragments)
		//- this comparator is used for collecting
		CompareTotalScore comparator( my_picker->get_score_manager() );
		CandidatesCollectorOP collector = new BoundedCollector<CompareTotalScore> (
			q_prof->length(),    // collector must know the size of query
			my_picker->n_candidates_, // how many candidates to collect
			comparator);    // yes, here the comparator comes to sort fragments within the collector
		my_picker->set_candidates_collector(3, collector);
		my_picker->set_candidates_collector(9, collector);
		my_picker->selector_ = new BestTotalScoreSelector(my_picker->n_frags_, scoring);

		my_picker->prefix_ = "div";

		//----------- WE SET UP BOUNDED COLLECTOR, WE RUN BOUNDED PROTOCOL
		//----------- TO RUN QUOTA PROTOCOL, ONE HAS TO SET UP QUOTA SELECTOR AND QUOTA COLLECTOR
		my_picker->bounded_protocol();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
