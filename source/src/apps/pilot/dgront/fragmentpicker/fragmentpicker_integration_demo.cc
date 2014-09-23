#include <devel/init.hh>

#include <core/sequence/Sequence.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>

// utility headers
#include <utility/vector1.hh>

#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/FragmentSelectingRule.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>
#include <protocols/frag_picker/DiversifyCrmsdByClustering.hh>

#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <utility>

#include <utility/excn/Exceptions.hh>

static thread_local basic::Tracer trace( "fragmentpicker_integration_demo" );

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( Boolean, frags, cluster_by_rms )

void register_options() {

	OPT(out::file::frag_prefix);
	OPT(frags::describe_fragments);
	OPT(in::file::native);
	OPT(in::file::psipred_ss2);
	OPT(in::file::checkpoint);
	OPT(frags::scoring::config);
	OPT(constraints::cst_file);
	OPT(in::path::database);
	OPT(frags::denied_pdb);
	OPT(frags::picking::query_pos);
	NEW_OPT(frags::cluster_by_rms,"",true);
}

int main(int argc, char * argv[]) {
    try {
	using namespace core;
	using namespace core::sequence;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	register_options();
	devel::init(argc, argv);

//------------ CREATE A PICKER OBJECT
//	FragmentPickerOP my_picker = new FragmentPicker("PValuedFragmentScoreManager");
	FragmentPickerOP my_picker( new FragmentPicker() );

//------------ PLUG-IN SEQUENCE PROFILE
	SequenceProfileOP q_prof( new SequenceProfile );
	q_prof->read_from_checkpoint(option[in::file::checkpoint]());
	my_picker->set_query_profile(q_prof);

//------------ LOAD SECONDARY STRUCTURE(s) - it is possible to load more than one
//------------ TWO PARAMETERS: SS_FILE AND A TAG
//----------- THE TAG MUST BE CONSISTENT WITH THE ONE GIVEN IN THE SCORE CONFIG FILE!!!!!
	my_picker->read_ss_file(option[in::file::psipred_ss2](),"johny");

//------------ READ VALL  AND PLUG IT INTO THE PICKER
	VallProviderOP chunks( new VallProvider() );
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
	Size n_scores = scoring->count_components();
	CompareTotalScore comparator( my_picker->get_score_manager() );
	CandidatesCollectorOP collector3( new BoundedCollector<CompareTotalScore> (
		q_prof->length(), 	  // collector must know the size of query
		my_picker->n_candidates_, // how many candidates to collect
		comparator,		  // yes, here the comparator comes to sort fragments within the collector
		n_scores) );
	CandidatesCollectorOP collector9( new BoundedCollector<CompareTotalScore> (
		q_prof->length(), 	  // collector must know the size of query
		my_picker->n_candidates_, // how many candidates to collect
		comparator,		  // yes, here the comparator comes to sort fragments within the collector
		n_scores) );
	my_picker->set_candidates_collector(3, collector3);
	my_picker->set_candidates_collector(9, collector9);

	if(option[frags::cluster_by_rms].user()) {
	    my_picker->selector_ = FragmentSelectingRuleOP( new DiversifyCrmsdByClustering(my_picker->n_frags_) );
	}
	else
    	    my_picker->selector_ = FragmentSelectingRuleOP( new BestTotalScoreSelector(my_picker->n_frags_, scoring) );

	my_picker->prefix_ = "fragments";
	if (option[out::file::frag_prefix].user())
	    my_picker->prefix_ = option[out::file::frag_prefix]();

//----------- SETUP QUERY POSITIONS
	if (option[frags::picking::query_pos].user())
	    my_picker->set_picked_positions( option[frags::picking::query_pos]() );

//----------- WE SET UP BOUNDED COLLECTOR, WE RUN BOUNDED PROTOCOL
//----------- TO RUN QUOTA PROTOCOL, ONE HAS TO SET UP QUOTA SELECTOR AND QUOTA COLLECTOR
	my_picker->bounded_protocol();
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}
