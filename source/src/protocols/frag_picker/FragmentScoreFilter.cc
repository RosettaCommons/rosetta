// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file --path--/--class--.cc
/// @brief --brief--
/// @author --name-- (--email--)

#include <protocols/frag_picker/FragmentScoreFilter.hh>
#include <protocols/frag_picker/FragmentScoreFilterCreator.hh>

//XSD includes
#include <numeric>      // std::accumulate
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/sequence/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/ResidueIndexDescription.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/TaskFactory.hh>
#include <utility/string_util.hh>

#include <core/scoring/dssp/Dssp.hh>

#include <protocols/filters/Filter.hh>

#include <basic/options/option.hh>

// Utility includes
#include <utility/vector1.hh>
#include <sys/stat.h>
#include <fstream>
#include <sys/param.h>
#include <unistd.h>
#include <boost/lexical_cast.hpp>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/prof.hh>
// XSD XRW Includes

#include <protocols/frag_picker/FragmentPicker.hh>
#include <protocols/frag_picker/scores/FragmentCrmsd.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>

static basic::Tracer TR( "protocols.frag_picker.FragmentScoreFilter" );

namespace protocols {
namespace frag_picker {

FragmentScoreFilter::FragmentScoreFilter():
	protocols::filters::Filter( "FragmentScoreFilter" )
{

}

FragmentScoreFilter::~FragmentScoreFilter() = default;

void
FragmentScoreFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & ,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	threshold_ = tag->getOption< core::Real >("threshold",1);
	direction_ = tag->getOption< std::string >("direction","-");
	score_type_ = tag->getOption< std::string >("scoretype","FragmentCrmsd");
	if ( tag->hasOption( "start_res" ) ) {
		start_res_ = core::pose::parse_resnum( tag->getOption< std::string >("start_res") );
	}
	if ( tag->hasOption( "end_res" ) ) {
		end_res_ = core::pose::parse_resnum( tag->getOption< std::string >("end_res") );
	}
	compute_ = tag->getOption< std::string >("compute","maximum");
	fragment_size_ = tag->getOption< core::Size >("fragment_size",9);
	sort_by_ = tag->getOption< std::string >("sort_by", "FragmentCrmsd" );
	outputs_folder_ = tag->getOption< std::string >("outputs_folder", "");
	outputs_name_ = tag->getOption<std::string >("outputs_name", "pose");
	csblast_ = tag->getOption< std::string >("csblast","");
	blast_pgp_ = tag->getOption< std::string >("blast_pgp","");
	placeholder_seqs_ = tag->getOption< std::string >("placeholder_seqs","");
	psipred_ = tag->getOption< std::string >("psipred","");
	sparks_x_ = tag->getOption< std::string >("sparks-x","");
	sparks_x_query_ = tag->getOption< std::string >("sparks-x_query","");
	vall_path_ = tag->getOption< std::string >("vall_path");
	frags_scoring_config_ = tag->getOption< std::string >("frags_scoring_config");
	n_frags_ = tag->getOption< core::Size >("n_frags",200);
	n_candidates_ = tag->getOption< core::Size >("n_candidates",1000);
	print_out_info_to_pdb_ = tag->getOption< bool >("print_to_pdb",false);

}

void convert_binary_checkpoint( std::string const & check_filename ){
	utility::vector1< core::Size > column_map{1, 5, 4, 7, 14, 8, 9, 10, 12, 11, 13, 3, 15, 6, 2, 16, 17, 20, 18, 19};

	double x = 0;
	char bb4[4];
	std::ifstream binary_check ( check_filename.c_str(), std::ios::in|std::ios::binary);
	if ( !binary_check ) {
		utility_exit_with_message("ERROR: Unable to open file " + check_filename);
	}

	binary_check.read(bb4, 4);

	int seqLength = bb4[0];
	int b2 = ((int) bb4[1]) << 8;
	int b3 = ((int) bb4[2]) << 16;
	int b4 = ((int) bb4[3]) << 24;
	seqLength = (seqLength < 0) ? 256 + seqLength + b2 + b3 + b4 : seqLength + b2 + b3 + b4;
	std::cout << seqLength << " " <<(int)bb4[0]<<" " <<(int) bb4[1] << "\n";

	utility::vector1<char> seq;
	for ( int i = 0; i<seqLength; ++i ) {
		seq.push_back( binary_check.get() );
	}
	std::string strSeq(seq.begin(), seq.end());
	TR.Debug << "Read sequence " << strSeq << " from " << check_filename << std::endl;


	std::ofstream checkpoint_file;
	checkpoint_file.open(check_filename + "point", std::ios_base::trunc);
	checkpoint_file << seqLength;
	checkpoint_file << "\n";
	char bb8[8];
	for ( int i = 0; i < seqLength; i++ ) {
		utility::vector1< core::Real > prof_row;
		utility::vector1< core::Real > arranged_row;
		checkpoint_file << strSeq.at(i);
		checkpoint_file << " ";
		for ( int j = 0; j < 20; j++ ) {
			binary_check.read(bb8, 8);
			std::copy(bb8, bb8 + sizeof(double), reinterpret_cast<char*>(&x));
			//     x = swap(x);
			std::cout << std::fixed<<std::setw(7)<<std::setprecision(4)<<x<<" ";
			prof_row.push_back(x);
		}
		for ( utility::vector1< core::Size >::const_iterator it = column_map.begin(); it != column_map.end(); ++it ) {
			arranged_row.push_back(prof_row[*it]);
		}
		for ( utility::vector1< core::Real >::const_iterator it = arranged_row.begin(); it!= arranged_row.end(); ++it ) {
			checkpoint_file << std::to_string(*it) + " ";
		}
		checkpoint_file << "\n";
	}
	checkpoint_file << "END";
	checkpoint_file.close();
}

/// @brief Function is required to get FragmentCrmsd scores from the
/// RosettaScripts pose; if FragmentCrmsd is calculated using the
/// ScoreMap, it tries to get a pose from the command line. Similar
/// functions will likely need to be implemented for pose-dependent
/// fragment scores.
/// TO DO: Make this function take weight and priority based on weights
/// file, then have it call a constructor with the RosettaScripts pose
/// for ANY fragment scoretype. This will rescore everything with the new
/// pose no matter what scoretype is used.
void FragmentScoreFilter::rescore_fragment_crmsd( core::pose::PoseOP const pose_chain_OP, utility::vector1<frag_picker::Candidate> candidates,  frag_picker::scores::FragmentCrmsd* fc )const{
	// Loop over best-scoring candidates and get the score component
	// specified by the user
	for ( frag_picker::Candidate ct : candidates ) {
		frag_picker::FragmentCandidateOP candidate_op = ct.first;
		frag_picker::scores::FragmentScoreMapOP scoremap_op = ct.second;

		// Change scoringmethod's pose
		fc->set_pose(pose_chain_OP);

		// Refill ScoreMap with FramgnetCrmsd score using RosettaScripts
		// pose
		TR.Debug << "FragmentScoreMap->get_score_components() before re-scoring: " << std::endl;
		TR.Debug << scoremap_op->get_score_components() << std::endl;
		fc->score(candidate_op, scoremap_op);
		TR.Debug << "FragmentScoreMap->get_score_components() after re-scoring: " << std::endl;
		TR.Debug << scoremap_op->get_score_components() << std::endl;

	}
}


/// @brief Run a command in a terminal
void run_command( std::string const & cmd ){
	// for( utility::vector1< std::string >::iterator cmd = commands.begin(); cmd != commands.end(); ++cmd ){
	TR.Debug << "Running command " << cmd << std::endl;
#ifdef __native_client__
		core::Size retval = 1;
#else
	core::Size retval = system( cmd.c_str() );
#endif
	if ( retval != 0 ) {
		utility_exit_with_message( "Failed to run command \"" + cmd + "\". Make sure you specify the full path to the relevant command in your XML file, and make sure the command works outside of Rosetta. Return code was " + boost::lexical_cast<std::string>( retval ) );
	}

}

/// @brief Setup the FragmentPicker. Replicates much of the
/// functionality from the FragmentPicker's parse_my_tag() function.
void FragmentScoreFilter::setup_fragment_picker( core::pose::PoseOP const pose_chain_OP, frag_picker::FragmentPickerOP picker ) const{
	// TODO: Make outputs_folder_directory


	// Record cwd and switch to outputs_folder_
	char buffer[MAXPATHLEN];
	char *path = getcwd(buffer, MAXPATHLEN);
	std::string original_path = path;
	if ( !outputs_folder_.empty() ) {
		int make_dir = mkdir( outputs_folder_.c_str(), 0700 );
		if ( make_dir == -1 ) {
			if ( errno != EEXIST ) {
				utility_exit_with_message( "Failed to create directory for sequence profile outputs." );
			}
		}
		int rc = chdir( outputs_folder_.c_str() );
		if ( rc < 0 ) {
			utility_exit_with_message("Unable to change directories to " + outputs_folder_ + " for FragmentPicker setup.");
		}
	}
	// check pose split to see if there are protein residues present
	std::string pose_seq = pose_chain_OP->sequence();
	std::string seq = "";
	for ( char j : pose_seq ) {
		if ( ( j == 'A' ) || ( j == 'C' ) || ( j == 'D' ) || ( j == 'E' ) ||
				( j == 'F' ) || ( j == 'G' ) || ( j == 'H' ) || ( j == 'I' ) ||
				( j == 'K' ) || ( j == 'L' ) || ( j == 'M' ) || ( j == 'N' ) ||
				( j == 'P' ) || ( j == 'Q' ) || ( j == 'R' ) || ( j == 'S' ) ||
				( j == 'T' ) || ( j == 'V' ) || ( j == 'W' ) || ( j == 'Y' ) ) {
			seq += j;
		} else {
			utility_exit_with_message("Chain contained non-protein residue, possibly a ligand. Make sure your start and end residues span only amino acid residues according to Rosetta numbering.");
		}
	}

	// Initialize filename variables
	std::string blast_filename;
	std::string pssm_filename;
	std::string check_filename;
	std::string checkpoint_filename;
	std::string fasta_pssm_filename;
	std::string ss2_filename;

	core::sequence::SequenceProfileOP q_prof( new core::sequence::SequenceProfile );


	// Generate fasta file
	// core::sequence::output_fasta_file(fasta_filename, pose);
	std::string fasta_filename = core::sequence::create_fasta_file( outputs_name_, seq );

	std::string fasta_sequence;
	std::ifstream in_fasta(fasta_filename);
	core::Size fasta_size = 0;
	// Count sequence length from fasta file to ensure the lengths
	// match (this is probably not necessary and may be removed later)
	std::string line;
	while ( getline( in_fasta, line ) ) {
		if ( line[0] != '>' ) {
			std::istringstream iss( line );
			char character = '\0';
			while ( iss.get(character) ) {
				fasta_size++;
				fasta_sequence += character;
			}
		}
	}

	if ( !csblast_.empty() ) {
		check_filename = outputs_name_  + ".check";
		std::string cmd = csblast_ + "/bin/csbuild -i " + fasta_filename + " -I fas -D " + csblast_ + "/data/K4000.crf -o " + check_filename + " -O chk";
		run_command( cmd );
		//q_prof->read_from_binary_chk(check_filename);
		//picker->set_query_seq(q_prof);
		//TR.Debug << "Added query from binary checkpoint with size " << picker->size_of_query() << " and sequence " << picker->get_query_seq_string() << std::endl;
		TR.Debug << "Sequence length (from fasta file): " << fasta_size << std::endl;

		checkpoint_filename = outputs_name_ + ".checkpoint";
		// std::string cmd = "python " + convert_checkpoint_ + " " + check_filename + " " + checkpoint_filename + " " + std::to_string(fasta_size);
		// run_command( cmd );
		convert_binary_checkpoint( check_filename );
		q_prof->read_from_checkpoint(checkpoint_filename);
		picker->set_query_seq(q_prof);
		TR.Debug << "Added query from checkpoint file" << std::endl;

	}

	std::string f_seq = core::sequence::read_fasta_file(fasta_filename)[1]->sequence();
	picker->set_query_seq(f_seq);
	TR.Debug << "Read in fasta sequence " << f_seq << "; query size is now " << picker->size_of_query() << " and FragmentPicker has a sequence of " << picker->get_query_seq_string() << std::endl;

	// Generate psipred secondary structure prediction file
	if ( !psipred_.empty() ) {
		std::string cmd = psipred_ + " " + fasta_filename;
		run_command( cmd );
		ss2_filename = outputs_name_ + ".ss2";
		TR.Debug << "Reading secondary structure prediction file " << ss2_filename << std::endl;
		picker->read_ss_file(ss2_filename, "predA");
	} else {
		// If unable to create psipred file, use dssp to get secondary
		// structure
		core::scoring::dssp::Dssp dssp( *pose_chain_OP );
		std::string secstruct = dssp.get_dssp_secstruct();
		picker->add_query_ss( secstruct, "predA" );

	}

	// Generate blast and pssm files (used in later SparkX step)
	if ( !blast_pgp_.empty() and !placeholder_seqs_.empty() ) {
		blast_filename = outputs_name_ + ".blast";
		std::ofstream blastfile( blast_filename );
		if ( !blastfile.is_open() ) {
			utility_exit_with_message( "Error: could not open " + blast_filename + " for writing." );
		}
		blastfile << outputs_name_ << " " << fasta_sequence << std::endl;
		if ( blastfile.bad() ) {
			utility_exit_with_message( "Encountered an error while writing to " + blast_filename );
		}
		blastfile.close();
		pssm_filename = outputs_name_ + ".pssm";
		std::string cmd = blast_pgp_ + " -i " + fasta_filename + " -B " + blast_filename + " -Q " + pssm_filename + " -t 1 -j 1 -h 0.001 -e 0.001 -b 0 -k 0 -d " + placeholder_seqs_;
		run_command( cmd );
	}

	// Copy pssm file into fasta.pssm and run Sparks-X
	if ( !sparks_x_.empty() and !sparks_x_query_.empty() ) {
		fasta_pssm_filename = outputs_name_ + ".fasta.pssm";
		std::ifstream  pssm(pssm_filename, std::ios::binary);
		std::ofstream  fasta_pssm(fasta_pssm_filename,   std::ios::binary);
		// Copy *.pssm to *.fasta.pssm so that sparks-x won't run a BLAST
		// query
		fasta_pssm << pssm.rdbuf();

		std::string cmd = "export SPARKSXDIR=" + sparks_x_ + " && " + sparks_x_query_ + " " + fasta_filename;
		run_command( cmd );
		picker->read_spine_x(outputs_name_ + ".fasta.phipsi");
	}

	if ( !outputs_folder_.empty() ) {
		int rc = chdir( original_path.c_str() );
		if ( rc < 0 ) {
			utility_exit_with_message("Unable to change directories back to " + original_path + " after FragmentPicker setup.");

		}
	}
	// Setup VALL
	PROF_START( basic::FRAGMENTPICKING_READ_VALL );
	picker->read_vall( vall_path_ );
	PROF_STOP( basic::FRAGMENTPICKING_READ_VALL );

	// Setup frag sizes
	utility::vector1< core::Size > size_vector;
	size_vector.push_back(fragment_size_);
	picker->set_frag_sizes(size_vector, true);

	// Setup scoring methods
	picker->create_scores( frags_scoring_config_ );

	// Setup number of fragments and candidates
	picker->set_n_frags( n_frags_ );
	core::Size n_candidates_scored = n_candidates_;
	if ( n_frags_ > n_candidates_scored ) {
		n_candidates_scored = n_frags_;
	}
	picker->set_n_candidates( n_candidates_scored );
	TR.Info << "Picking " << n_frags_ << " fragments based on " << n_candidates_ << " candidates" << std::endl;

	// Comparator used for collecting & selecting fragments
	protocols::frag_picker::CompareTotalScore comparator( picker->get_score_manager() );

	// Setup scoring scheme for selection step
	TR.Info << "Creating fragment scoring scheme for the selection step" << std::endl;
	protocols::frag_picker::scores::FragmentScoreManagerOP selection_scoring;
	protocols::frag_picker::FragmentSelectingRuleOP selector( new protocols::frag_picker::BestTotalScoreSelector( n_frags_, selection_scoring ) );
	picker->set_selector( selector );

	// Collectors & selector setup for bounded protocol
	protocols::frag_picker::CandidatesCollectorOP collector( new protocols::frag_picker::BoundedCollector<protocols::frag_picker::CompareTotalScore> (picker->size_of_query(), picker->get_n_candidates(), comparator, picker->get_score_manager()->count_components() ) );
	picker->set_candidates_collector( fragment_size_, collector, 1 );


}

/// @brief Compute filter value
core::Real FragmentScoreFilter::compute( core::pose::Pose const pose ) const
{
	core::Size start_res_orig = 1;
	if ( start_res_ != nullptr ) {
		start_res_orig = start_res_->resolve_index( pose );
	}
	core::Size end_res_orig = pose.size();
	if ( end_res_ != nullptr ) {
		end_res_orig = end_res_->resolve_index( pose );
	}

	// Set up pose for chain of interest only
	if ( pose.chain(start_res_orig) != pose.chain(end_res_orig) ) {
		utility_exit_with_message( "Start and end residues must be on the same chain. To score a whole pose, use a separate filter for each chain." );
	}
	core::pose::PoseOP pose_chain_OP = pose.split_by_chain( pose.chain( start_res_orig ) );
	core::Size residue_adjust = 0;
	for ( core::Size i = 1; i < pose.chain( start_res_orig ); i++ ) {
		TR.Debug << "The size of chain " << i << " is " << pose.chain_sequence( i ).size() << std::endl;
		residue_adjust += pose.chain_sequence( i ).size();
	}
	TR.Debug << "Adjusting start and end residues by " << residue_adjust << ", with initial start and end residues of " << start_res_orig << " and " << end_res_orig << std::endl;
	core::Size start_res = start_res_orig - residue_adjust;
	core::Size end_res = end_res_orig - residue_adjust;
	TR.Debug << "Start residue is now " << start_res << " and end residue is " << end_res << std::endl;

	// Get secondary structure from pose
	TR.Debug << "Full input sequence is: " << std::endl;
	TR.Debug << pose_chain_OP->sequence() << std::endl;

	// Setting up FragmentPicker
	TR.Info << "Starting FragmentScorefilter..." << std::endl;
	frag_picker::FragmentPickerOP picker(new frag_picker::FragmentPicker);
	// std::string query = pose_chain_OP->sequence();
	// picker->set_query_seq(query);
	setup_fragment_picker( pose_chain_OP, picker );
	picker->set_picked_positions(start_res, end_res);


	// Running FragmentPicker
	picker->pick_candidates();
	TR.Info << "Computing average fragment score..." << std::endl;
	frag_picker::scores::FragmentScoringMethodOP score_component = picker->get_score_manager()->get_component_by_name(score_type_);
	TR.Info << "Filtering based on the following metric:" << score_component << std::endl;
	utility::vector1<frag_picker::Candidates> final_fragments;

	// Loop over residues in pose segment
	core::Size maxqpos = picker->size_of_query() - fragment_size_ + 1;
	TR.Debug << "Query sequence is: " << picker->get_picked_positions() << std::endl;
	TR.Debug << "Query size is: " << picker->size_of_query() << std::endl;

	utility::vector1< core::Size > picked_positions = picker->get_picked_positions();
	for ( core::Size iqpos = 1; iqpos <= picked_positions.size(); ++iqpos ) {
		core::Size qPos = picked_positions[iqpos];

		if ( qPos > maxqpos ) {
			TR.Info << "Attempted to get fragments for residue " << qPos << " but this is too close to the end of the query sequence. Skipping this position." << std::endl;
			continue;
		}

		TR.Info << "Using user-defined fragment size" << fragment_size_ << " at position " << qPos << std::endl;
		frag_picker::Candidates out;
		picker->get_final_fragments(qPos,fragment_size_,out);
		final_fragments.push_back(out);
	}

	//Rescore Crmsd with RosettaScripts pose
	if ( score_type_ == "FragmentCrmsd" || sort_by_ == "FragmentCrmsd" ) {
		TR.Debug << "Recalculating FragmentCrmsd" << std::endl;
		frag_picker::scores::FragmentCrmsd * fragment_crmsd = static_cast<frag_picker::scores::FragmentCrmsd*>( &(*score_component) );
		for ( frag_picker::Candidates cs : final_fragments ) {
			FragmentScoreFilter::rescore_fragment_crmsd(pose_chain_OP, cs, fragment_crmsd);
		}
	}


	utility::vector1<frag_picker::Candidate> total_scores;
	frag_picker::scores::FragmentScoreManagerOP ms = picker->get_score_manager();



	// For each residue position, sort final fragments by either total
	// score or some score chosen in the XML file.
	if ( sort_by_ == "TotalScore" ) {
		TR.Debug << "Sorting by total score" <<std::endl;
		// Loop over final candidates and find the best-scoring ones
		for ( frag_picker::Candidates cs : final_fragments ) {
			frag_picker::Candidate best_candidate = cs.front();
			core::Real best_total_score = ms->total_score(best_candidate.second);
			for ( frag_picker::Candidate ct : cs ) {
				core::Real ct_total_score = ms->total_score(ct.second);
				if ( ct_total_score < best_total_score ) {
					best_candidate = ct;
					best_total_score = ct_total_score;
				}
			}
			total_scores.push_back(best_candidate);
		}
	} else {
		core::Size sortby_id = picker->get_score_manager()->get_component_by_name(sort_by_)->get_id();
		TR.Debug << "Sorting by " << sort_by_ << " whose component ID is " << sortby_id << std::endl;
		core::Size res = 1;
		for ( frag_picker::Candidates cs : final_fragments ) {
			frag_picker::Candidate best_candidate = cs.front();
			frag_picker::scores::FragmentScoreMap initial_scoremap = *best_candidate.second;
			utility::vector1<core::Real> initial_score_components = initial_scoremap.get_score_components();
			core::Real best_score = initial_score_components[sortby_id];
			for ( frag_picker::Candidate ct : cs ) {
				frag_picker::scores::FragmentScoreMap candidate_scoremap = *ct.second;
				utility::vector1<core::Real> score_components = candidate_scoremap.get_score_components();
				core::Real score = score_components[sortby_id];
				TR.Debug << "Evaluating residue " << res << std::endl;
				TR.Debug << "Fragment PDBID: "<< ct.first->get_pdb_id() << std::endl;
				TR.Debug << "FRAG_ID: " << std::to_string(ct.first->key()) << std::endl;
				TR.Debug << "Scoremap evaluated is " << std::endl;
				TR.Debug  << score_components;
				TR.Debug << "; relevant score is " << std::to_string(score) << std::endl;
				if ( score < best_score ) {
					best_candidate = ct;
					best_score = score;
					TR.Debug << "====Found new best score: " << score << "====" << std::endl;
				}
				TR.Debug << "-----------------------------------------" << std::endl;
			}
			total_scores.push_back(best_candidate);
			res += 1;
		}
	}

	utility::vector1<core::Real> scores;
	int score_id = score_component->get_id();

	// Loop over best-scoring candidates and get the score component
	// specified by the user
	for ( frag_picker::Candidate ct : total_scores ) {
		frag_picker::scores::FragmentScoreMap scoremap = *ct.second;
		utility::vector1<core::Real> score_components = scoremap.get_score_components();
		scores.push_back(score_components[score_id]);
	}

	core::Real result = get_result( scores );
	TR.Debug << "Scores:" << scores << std::endl;


	std::ostringstream oss;
	if ( print_out_info_to_pdb_ ) {
		std::string filter_name = this->name();
		std::string user_name = this->get_user_defined_name();
		oss << std::endl << filter_name << " " << user_name + ": " << std::endl;

		core::Size index = 1 + start_res;
		for ( core::Real score : scores ) {
			core::conformation::Residue residue = pose_chain_OP->residue(index);
			std::string restype = residue.name3();
			std::string resnum = utility::to_string(index - start_res + start_res_orig);
			std::string temp_str = "FragmentScoreFilter_metric " + restype + " " + resnum + " " + score_type_ + ": " + utility::to_string(score);
			TR.Info << temp_str << std::endl;
			oss << temp_str << std::endl;
			++index;
		}
		runtime_assert( index == end_res );

		core::Real min_score = *std::min_element( scores.begin(), scores.end() );
		core::Size min_residue = std::distance( scores.begin(), std::min_element( scores.begin(), scores.end() ) ) + start_res_orig;

		core::Real max_score = *std::max_element( scores.begin(), scores.end() );
		core::Size max_residue = std::distance( scores.begin(), std::max_element( scores.begin(), scores.end() ) ) + start_res_orig;

		core::Real average_score = std::accumulate( scores.begin(), scores.end(), 0.0 ) / scores.size();

		oss << "FragmentScoreFilter_metric Max " << score_type_ << " res: " << utility::to_string(max_residue) << std::endl;
		oss << "FragmentScoreFilter_metric Max " << score_type_ << " score: " << utility::to_string(max_score) << std::endl;
		oss << "FragmentScoreFilter_metric Min " << score_type_ << " res: " << utility::to_string(min_residue) << std::endl;
		oss << "FragmentScoreFilter_metric Min " << score_type_ << " score: " << utility::to_string(min_score) << std::endl;
		oss << "FragmentScoreFilter_metric Avg " << score_type_ << " score: " << utility::to_string(average_score) << std::endl << std::endl;

		protocols::jd2::JobOP job(protocols::jd2::JobDistributor::get_instance()->current_job());
		job->add_string( oss.str() );
	}

	return result;
}

protocols::filters::FilterOP
FragmentScoreFilter::clone() const
{
	return protocols::filters::FilterOP( new FragmentScoreFilter( *this ) );
}


protocols::filters::FilterOP
FragmentScoreFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new FragmentScoreFilter );
}

/// @brief Given a list of scores, find the one the user is interested in.
core::Real FragmentScoreFilter::get_result( utility::vector1<core::Real> scores ) const {
	core::Real result;
	if ( compute_ == "minimum" || compute_ == "min" ) {
		result = *std::min_element(scores.begin(),scores.end());
	} else if ( compute_ == "maximum" || compute_ == "max" ) {
		result = *std::max_element(scores.begin(),scores.end());
	} else {
		result = std::accumulate(scores.begin(),scores.end(),0.0)/scores.size();
	}
	return result;
}

bool
FragmentScoreFilter::apply( core::pose::Pose const & pose ) const
{

	if ( direction_ != "+" and direction_ != "-" ) {
		utility_exit_with_message( "\"direction\" must be set to either \"+\" or \"-\".");
	}
	core::Real average = compute(pose);

	TR.Debug << "direction_ is " << direction_ << ", score is " << average <<", and threshold is " << threshold_ << std::endl;
	if ( average > threshold_ ) {
		return direction_ == "+";
	} else if ( average < threshold_ ) {
		return direction_ == "-";
	}
	return true;
}

core::Real
FragmentScoreFilter::report_sm( core::pose::Pose const & pose ) const
{
	core::Real average = compute( pose );
	return average;
}

void
FragmentScoreFilter::report( std::ostream &, core::pose::Pose const & pose ) const
{
	compute( pose );
}

std::string FragmentScoreFilter::name() const {
	return class_name();
}

std::string FragmentScoreFilter::class_name() {
	return "FragmentScoreFilter";
}

void FragmentScoreFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	std::string dbpath =basic::options::option[basic::options::OptionKeys::in::path::database][1];
	using namespace utility::tag;
	AttributeList attlist;

	//here you should write code to describe the XML Schema for the class.  If it has only attributes, simply fill the probided AttributeList.
	attlist + XMLSchemaAttribute::attribute_w_default( "scoretype", xs_string, "Which attribute to filter on. See FragmentScoreManager::register_score_maker for options.", "FragmentCrmsd" )
		+ XMLSchemaAttribute::attribute_w_default( "threshold", xsct_real, "Filter threshold.", "6" )
		+ XMLSchemaAttribute::attribute_w_default( "direction", xs_string, "Choose whether you want outputs to be greater or less than the threshold. Right now the only options are greater than ('+') or less than ('-'); if you put anything else, all will pass.", "-" )
		+ XMLSchemaAttribute::attribute_w_default("compute",xs_string,"How to calculate filter value. Right now the options are average, maximum, or minimum of the scores collected.","average")
		+ XMLSchemaAttribute::attribute_w_default("fragment_size",xsct_non_negative_integer, "Size of fragments to be computed (default 3)", "3")
		+ XMLSchemaAttribute::attribute_w_default("sort_by", xs_string, "Choose how to pick the best fragment from the final list of candidates at each position. Default is by FragmentCrmsd, and RamaScore, SecondarySimilarity, and TotalScore are enabled by default. You can use any other fragment score if you include it via a weights file by using the command line option \"-frags::scoring::config\"", "FragmentCrmsd")
		+ XMLSchemaAttribute::attribute_w_default("outputs_folder",xs_string,"Folder where you would like external files to go, such as secondary structure prediction and blast checkpoint files. If not set, files will be placed in current working directory.","")
		+ XMLSchemaAttribute::attribute_w_default("outputs_name",xs_string,"Base name of intermediate files, i.e. the name of your design.","pose")
		+ XMLSchemaAttribute::attribute_w_default("csblast",xs_string,"Directory of csblast program, ex. software/csblast-2.2.3_linux64","")
		+ XMLSchemaAttribute::attribute_w_default("blast_pgp",xs_string,"Path to blastpgp for creating psipred file.","")
		+ XMLSchemaAttribute::attribute_w_default("placeholder_seqs",xs_string,"Path to placeholder sequence database for creating psipred file.","")
		+ XMLSchemaAttribute::attribute_w_default("sparks-x",xs_string,"Path to sparks-x directory","")
		+ XMLSchemaAttribute::attribute_w_default("sparks-x_query",xs_string,"Path to sparks-x script, \"buildinp_query.sh\"","")
		+ XMLSchemaAttribute::attribute_w_default("psipred",xs_string,"Path to runpsipred_single script.","")
		+ XMLSchemaAttribute::attribute_w_default("vall_path",xs_string,"Path to vall database.",dbpath + "/sampling/vall.jul19.2011.gz")
		+ XMLSchemaAttribute("frags_scoring_config",xs_string,"Path to scoring config file (required).")
		+ XMLSchemaAttribute::attribute_w_default("n_frags",xsct_non_negative_integer, "Number of fragments to be picked (default 200)", "200")
		+ XMLSchemaAttribute::attribute_w_default("n_candidates", xsct_non_negative_integer, "Number of candidates per position (default 1000)","1000")
		+ XMLSchemaAttribute::attribute_w_default("print_to_pdb", xs_boolean, "Prints scores for all residues analyzed to the pdb.", "false");

	core::pose::attributes_for_parse_resnum( attlist, "start_res", "The N-terminal residue of the piece of backbone to be analyzed.");
	core::pose::attributes_for_parse_resnum( attlist, "end_res", "The C-terminal residue of the piece of backbone to be analyzed. ");

	rosetta_scripts::attributes_for_parse_task_operations( attlist );

	protocols::filters::xsd_type_definition_w_attributes(
		xsd,
		class_name(),
		"Filter based on any score that can be calculated in fragment_picker.",
		attlist );
}

/////////////// Creator ///////////////

protocols::filters::FilterOP
FragmentScoreFilterCreator::create_filter() const
{
	return protocols::filters::FilterOP( new FragmentScoreFilter );
}

std::string
FragmentScoreFilterCreator::keyname() const
{
	return FragmentScoreFilter::class_name();
}

void FragmentScoreFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	FragmentScoreFilter::provide_xml_schema( xsd );
}

}}
