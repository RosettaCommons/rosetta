#include <devel/init.hh>

#include <core/sequence/Sequence.hh>
#include <basic/Tracer.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
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

#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <utility>

#include <utility/excn/Exceptions.hh>

//Auto Headers
#include <utility/io/mpistream.hh>


static basic::Tracer trace( "rescore_fragments" );

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

void register_options() {

	OPT(in::file::native);
	OPT(in::file::s);
	OPT(in::file::psipred_ss2);
	OPT(in::file::fasta);
	OPT(in::file::pssm);
	OPT(in::file::checkpoint);
	OPT(frags::scoring::config);
	OPT(constraints::cst_file);
	OPT(in::path::database);
	OPT(in::file::fragA);
	OPT(in::file::fragB);
	OPT(in::file::frag3);
	OPT(in::file::frag9);
	OPT(frags::allowed_pdb);
	OPT(frags::denied_pdb);
	OPT(frags::scoring::profile_score);
}

int main(int argc, char * argv[]) {
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init(argc, argv);

		if ( option[in::file::native].user() ) {
			trace.Debug << option[in::file::native]() << std::endl;
		}

		//---------- Set up a picker.
		FragmentPickerOP pickIt = new FragmentPicker();
		pickIt->parse_command_line();
		trace << "After setup; size of a query is: " << pickIt->size_of_query()
			<< std::endl;
		scores::FragmentScoreManagerOP manager = pickIt->get_score_manager();
		VallProviderOP vall = pickIt->get_vall();

		//-------- read fragments
		utility::vector1<std::string> frag_files;
		if ( option[in::file::fragA].user() ) {
			frag_files.push_back( option[in::file::fragA]() );
		}
		if ( option[in::file::fragB].user() ) {
			frag_files.push_back( option[in::file::fragB]() );
		}
		if ( option[in::file::frag9].user() ) {
			frag_files.push_back( option[in::file::frag9]() );
		}
		if ( option[in::file::frag3].user() ) {
			frag_files.push_back( option[in::file::frag3]() );
		}
		for ( Size i_file=1; i_file<=frag_files.size(); ++i_file ) {
			std::cerr << "# ------------ Fragments from "<<frag_files[i_file]<<" ---------------"<<std::endl;
			utility::vector1<FragmentCandidateOP> frags = read_fragment_candidates(
				frag_files[i_file], vall);

			utility::vector1<std::pair<FragmentCandidateOP,
				scores::FragmentScoreMapOP> > frags_rescored;
			for ( Size i = 1; i <= frags.size(); i++ ) {
				scores::FragmentScoreMapOP map = manager->create_empty_map();
				manager->do_caching(frags[i]->get_chunk());
				trace.Info << "rescoring: " << *(frags[i]) << std::endl;
				manager->score_fragment(frags[i], map);
				manager->clean_up();
				std::pair<FragmentCandidateOP, scores::FragmentScoreMapOP> p(
					frags[i], map);
				frags_rescored.push_back(p);
			}
			pickIt->get_score_manager()->describe_fragments(frags_rescored,
				std::cerr);
		}
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}

