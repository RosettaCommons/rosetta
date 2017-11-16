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

#include <core/fragment/picking/FragmentPicker.hh>
#include <core/fragment/picking/VallProvider.hh>
#include <core/fragment/picking/VallChunk.hh>
#include <core/fragment/picking/VallResidue.hh>
#include <core/fragment/picking/FragmentCandidate.hh>
#include <core/fragment/picking/FragmentSelectingRule.hh>
#include <core/fragment/picking/VallChunkFilter.hh>
#include <core/fragment/picking/BoundedCollector.hh>
#include <core/fragment/picking/BestTotalScoreSelector.hh>
#include <core/fragment/picking/PdbIdChunkFilter.hh>

#include <core/fragment/picking/scores/FragmentScoringMethod.hh>
#include <core/fragment/picking/scores/ProfileScoreL1.hh>
#include <core/fragment/picking/scores/SecondarySimilarity.hh>
#include <core/fragment/picking/scores/FragmentScoreManager.hh>
#include <core/fragment/picking/scores/RamaScore.hh>
#include <utility>

#include <core/fragment/picking/quota/QuotaCollector.hh>
#include <core/fragment/picking/quota/QuotaPool.hh>
#include <core/fragment/picking/quota/ABEGO_SS_Pool.hh>
#include <core/fragment/picking/quota/ABEGO_SS_Config.hh>

//Auto Headers
#include <utility/io/mpistream.hh>

// boost
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <utility/excn/Exceptions.hh>

static basic::Tracer trace( "quotapicker_demo" );

using namespace core;
using namespace core::fragment;
using namespace core::fragment::picking;
using namespace core::fragment::picking::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

OPT_1GRP_KEY( String, frags, abego_logratios )
OPT_1GRP_KEY( String, frags, wildparam )

void register_options() {

	OPT(in::file::native);
	OPT(in::file::psipred_ss2);
	OPT(in::file::checkpoint);
	OPT(frags::scoring::config);
	OPT(frags::picking::query_pos);
	OPT(constraints::cst_file);
	OPT(in::path::database);
	OPT(frags::denied_pdb);
	OPT(frags::ss_pred);
	OPT(frags::denied_pdb);
	OPT(frags::allowed_pdb);
	OPT(in::file::torsion_bin_probs);
	NEW_OPT(frags::abego_logratios,"","");
	NEW_OPT(frags::wildparam,"","");
}

void attach_abego_pools(FragmentPickerOP picker,Size n_candidates_) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string quota_config_file("UNKNOWN-QUOTA-PROBABILITY_FILE");
	if ( option[in::file::torsion_bin_probs].user() ) {
		quota_config_file = option[in::file::torsion_bin_probs]();
	}
	quota::ABEGO_SS_Config q_config(quota_config_file);

	std::string quota_logratios_file("UNKNOWN-QUOTA-LOGRATIO_FILE");
	if ( option[frags::abego_logratios].user() ) {
		quota_logratios_file = option[frags::abego_logratios]();
	}
	quota::ABEGO_SS_Config l_config(quota_logratios_file);

	scores::FragmentScoreManagerOP scores_ = picker->get_score_manager();
	utility::vector1<Size> components;
	utility::vector1<Real> weights;
	utility::vector1<Real> scoring_weights = scores_->get_weights();
	bool has_ss = false;
	for ( Size i = 1; i <= scores_->count_components(); ++i ) {
		ProfileScoreL1 *s1 =
			dynamic_cast<ProfileScoreL1*> (scores_->get_component(i).get());
		if ( s1 != 0 ) {
			components.push_back( i );
			weights.push_back( scoring_weights[s1->get_id()] );
		}

		RamaScore *s2 =
			dynamic_cast<RamaScore*> (scores_->get_component(i).get());
		if ( s2 != 0 ) {
			components.push_back( i );
			weights.push_back( scoring_weights[s2->get_id()] );
		}

		SecondarySimilarity *s4 =
			dynamic_cast<SecondarySimilarity*> (scores_->get_component(i).get());
		if ( (s4 != 0)&&(has_ss==false) ) {
			components.push_back( i );
			weights.push_back( scoring_weights[s4->get_id()] );
			has_ss = true;
		}
	}

	Real global_fraction = 0.25;
	Size buffer_factor = 5;
	Real score_cutoff = 1.5;
	Real cutoffs[8] = {0.0,10.0,10.0,score_cutoff,score_cutoff,score_cutoff,score_cutoff,score_cutoff};
	for ( Size f=1; f<=picker->frag_sizes_.size(); f++ ) {
		Size f_size = picker->frag_sizes_[f];
		quota::QuotaCollectorOP collector = dynamic_cast<quota::QuotaCollector*> (picker->get_candidates_collector(f_size).get());
		Size middle = f_size / 2 + 1;
		for ( Size j=1; j<=picker->size_of_query() - f_size +1; j++ ) {
			//  Size i = l_config.most_probable_bin(j+middle-1);
			//  Real l = l_config.probability(j+middle-1,i);  // get log-ratio for the center
			for ( Size i=1; i<=7; i++ ) {
				Real l =  l_config.probability(j+middle-1,i);            // get log-ratio for the center
				if ( l>cutoffs[i] ) {
					//      std::cerr<<l<<" "<<cutoffs[i]<<" "<<i<<" "<<q_config.get_pool_name(i)<<"\n";
					Real ff = q_config.probability(j+middle-1,i); // get true posterior probability
					ff = ff * global_fraction;
					quota::QuotaPoolOP p = (quota::QuotaPool*) new quota::ABEGO_SS_Pool(n_candidates_,q_config.get_pool_name(i),
						q_config.get_pool_bins((i)),components,weights,ff,scores_->count_components(),buffer_factor);
					collector->add_pool(j,p);
					trace.Info << "Attached ABEGO pool for label "<<q_config.get_pool_name(i)
						<<" with fraction: "<<q_config.probability(j+middle-1,i)<<" x "<<global_fraction
						<<" = "<<ff<<" Its logratio is: "<<l <<std::endl;
				}
			}
		}
	}
}


void attach_simple_abego_pools(FragmentPickerOP picker,Size n_candidates_,
	std::string ss1, std::string abego1, Real f1,
	std::string ss2, std::string abego2, Real f2,
	std::string ss3, std::string abego3, Real f3) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	Size buffer_factor = 5;
	utility::vector1< quota::ABEGO_SS_Config > q_cfgs;
	quota::ABEGO_SS_Config q_config1(abego1);
	quota::ABEGO_SS_Config q_config2(abego2);
	quota::ABEGO_SS_Config q_config3(abego3);
	q_cfgs.push_back( q_config1 );
	q_cfgs.push_back( q_config2 );
	q_cfgs.push_back( q_config3 );
	utility::vector1< Real > fractions;
	fractions.push_back( f1 );
	fractions.push_back( f2 );
	fractions.push_back( f3 );
	utility::vector1< std::string > prediction_names;
	prediction_names.push_back( ss1 );
	prediction_names.push_back( ss2 );
	prediction_names.push_back( ss3 );
	utility::vector1< std::string > abego_names;
	abego_names.push_back( abego1 );
	abego_names.push_back( abego2 );
	abego_names.push_back( abego3 );

	scores::FragmentScoreManagerOP scores_ = picker->get_score_manager();
	for ( Size i_cfg = 1; i_cfg <= 3; i_cfg++ ) {

		utility::vector1<Size> components;
		utility::vector1<Real> weights;
		utility::vector1<Real> scoring_weights = scores_->get_weights();
		bool has_ss = false;
		for ( Size i = 1; i <= scores_->count_components(); ++i ) {
			ProfileScoreL1 *s1 =
				dynamic_cast<ProfileScoreL1*> (scores_->get_component(i).get());
			if ( s1 != 0 ) {
				components.push_back( i );
				weights.push_back( scoring_weights[s1->get_id()] );
			}

			RamaScore *s2 =
				dynamic_cast<RamaScore*> (scores_->get_component(i).get());
			if ( s2 != 0 ) {
				components.push_back( i );
				weights.push_back( scoring_weights[s2->get_id()] );
			}

			SecondarySimilarity *s4 =
				dynamic_cast<SecondarySimilarity*> (scores_->get_component(i).get());
			if ( (s4 != 0)&&(has_ss==false) ) {
				if ( s4->get_prediction_name().compare( prediction_names[i_cfg] ) == 0 ) {
					components.push_back( i );
					weights.push_back( scoring_weights[s4->get_id()] );
					has_ss = true;
				}
			}
		}

		for ( Size f=1; f<=picker->frag_sizes_.size(); f++ ) {
			trace.Info << "Creating a pool from the data: "<<prediction_names[i_cfg]<<" "<<
				abego_names[i_cfg]<<" "<<fractions[i_cfg]<<std::endl;
			Size f_size = picker->frag_sizes_[f];

			quota::QuotaCollectorOP collector = new quota::QuotaCollector( picker->size_of_query(), f_size );
			picker->set_candidates_collector(f_size,collector);

			Size middle = f_size / 2 + 1;
			quota::ABEGO_SS_Config & q_config = q_cfgs[i_cfg];
			for ( Size j=1; j<=picker->size_of_query() - f_size +1; j++ ) {

				trace.Trace <<"Probabilities for middle pos. " << j+middle-1 << " in config no. " << i_cfg<< ":\n";
				for ( Size i=1; i<=5; i++ ) trace.Trace << q_cfgs[i_cfg].probability(j+middle-1,i) << " ";
				trace.Trace <<std::endl;

				for ( Size i=1; i<=5; i++ ) {
					Real ff = q_config.probability(j+middle-1,i); // get true posterior probability
					ff = ff * fractions[i_cfg];
					if ( n_candidates_ * ff < 1.0 ) {
						trace.Info << "Skipping a pool " << q_config.get_pool_name(i) << " at pos " << j
							<< "because quota fraction is: " << ff << std::endl;
					} else {
						quota::QuotaPoolOP p = (quota::QuotaPool*) new quota::ABEGO_SS_Pool(n_candidates_,q_config.get_pool_name(i),
							q_config.get_pool_bins((i)),components,weights,ff,scores_->count_components(),buffer_factor);
						collector->add_pool(j,p);
						trace.Info << "Attached ABEGO pool for label "<<q_config.get_pool_name(i)
							<<" with fraction: "<<q_config.probability(j+middle-1,i)<<" x "<<fractions[i_cfg]
							<<" = "<<ff<<std::endl;
					}
				}
			}
		}
	}
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
		FragmentPickerOP my_picker = new FragmentPicker();

		//------------ PLUG-IN SEQUENCE PROFILE
		SequenceProfileOP q_prof(new SequenceProfile);
		q_prof->read_from_checkpoint(option[in::file::checkpoint]());
		my_picker->set_query_profile(q_prof);

		//------------ LOAD SECONDARY STRUCTURE(s) - it is possible to load more than one
		//------------ TWO PARAMETERS: SS_FILE AND A TAG (possibly repeated)
		my_picker->read_ss_files(option[frags::ss_pred]());

		//------------ READ VALL  AND PLUG IT INTO THE PICKER
		VallProviderOP chunks = new VallProvider();
		chunks->vallChunksFromLibrary(option[in::file::vall]()[1]);
		my_picker->set_vall(chunks);

		//---__------- setup chunk filters
		if ( option[frags::allowed_pdb].user() ) {
			AllowPdbIdFilterOP allow = new AllowPdbIdFilter();
			allow->load_pdb_id_from_file(option[frags::allowed_pdb]());
			my_picker->add_chunk_filter(allow);
			trace.Info << "Allowed PDB chains:\n";
			allow->show_pdb_ids(trace.Info);
		}

		if ( option[frags::denied_pdb].user() ) {
			DenyPdbIdFilterOP deny = new DenyPdbIdFilter();
			deny->load_pdb_id_from_file(option[frags::denied_pdb]());
			my_picker->add_chunk_filter(deny);
			trace.Info << "Excluded PDB chains:\n";
			deny->show_pdb_ids(trace.Info);
		}

		//------------ FRAGMENT SIZE: WE NEED 9-MERS AND 3-MERS
		my_picker->frag_sizes_.push_back(3);
		my_picker->frag_sizes_.push_back(9);

		//------------ HOW MANY CANDIDATES, HOW MANY FRAGMENTS
		my_picker->n_candidates_ = option[frags::n_candidates]();
		my_picker->n_frags_ = option[frags::n_frags]();
		trace.Info << "Picking " << my_picker->n_frags_ << " fragments based on "<<
			my_picker->n_candidates_ << " candidates" << std::endl;

		my_picker->prefix_ = "fragments";

		//----------- SETUP SCORING SYSTEM
		FragmentScoreManagerOP scoring = my_picker->get_score_manager();
		scoring->create_scores(option[frags::scoring::config](), my_picker);

		//----------- SETUP QUOTA COLLECTOR (for candidates).
		//----------- QUOTA ALWAYS USES QUOTA SELECTOR (we don't need to specify a selector)
		//----------- This is done by set_up_quota_nnmake_style() method
		if ( option[frags::wildparam].user() ) {
			typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
			boost::char_separator<char> sep(", :");
			tokenizer tokens(option[frags::wildparam](), sep);
			utility::vector1<std::string> t;
			for ( tokenizer::iterator tok_iter = tokens.begin(); tok_iter != tokens.end(); ++tok_iter ) {
				t.push_back(*tok_iter);
			}

			attach_simple_abego_pools(my_picker,my_picker->n_candidates_,
				t[1],t[2],boost::lexical_cast<double>(t[3]),
				t[4],t[5],boost::lexical_cast<double>(t[6]),
				t[7],t[8],boost::lexical_cast<double>(t[9]));

		} else {
			my_picker->set_up_quota_nnmake_style();
			if ( option[frags::abego_logratios].user() ) {
				attach_abego_pools(my_picker,my_picker->n_candidates_);
			}
		}
		// my_picker->set_up_ss_abego_quota();
		//----------- SETUP QUERY POSITIONS
		if ( option[frags::picking::query_pos].user() ) {
			my_picker->set_picked_positions( option[frags::picking::query_pos]() );
		}

		//----------- TO RUN QUOTA PROTOCOL, ONE HAS TO SET UP QUOTA SELECTOR AND QUOTA COLLECTOR
		my_picker->quota_protocol();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;

}
