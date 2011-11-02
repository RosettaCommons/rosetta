// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/frag_picker/FragmentPicker.cc
/// @brief  Fragment picker - the core part of picking machinery
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// unit headers
#include <protocols/frag_picker/FragmentPicker.hh>

#include <protocols/frag_picker/VallProvider.hh>
#include <protocols/frag_picker/VallChunk.hh>
#include <protocols/frag_picker/VallResidue.hh>
#include <protocols/frag_picker/VallChunkFilter.hh>
#include <protocols/frag_picker/FragmentCandidate.hh>
#include <protocols/frag_picker/CandidatesCollector.hh>
#include <protocols/frag_picker/GrabAllCollector.hh>
#include <protocols/frag_picker/quota/QuotaSelector.hh>
#include <protocols/frag_picker/quota/QuotaConfig.hh>
#include <protocols/frag_picker/quota/QuotaCollector.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Config.hh>
#include <protocols/frag_picker/quota/ABEGO_SS_Pool.hh>
#include <protocols/frag_picker/quota/SecondaryStructurePool.hh>
#include <protocols/frag_picker/BoundedCollector.hh>
#include <protocols/frag_picker/PdbIdChunkFilter.hh>
#include <protocols/frag_picker/BestTotalScoreSelector.hh>
#include <protocols/frag_picker/CustomScoreSelector.hh>
#include <protocols/frag_picker/scores/FragmentScoringMethod.hh>
#include <protocols/frag_picker/scores/FragmentScoreManager.hh>
#include <protocols/frag_picker/scores/SecondarySimilarity.hh>
#include <protocols/frag_picker/scores/PartialSecondarySimilarity.hh>
#include <protocols/frag_picker/scores/ProfileScoreL1.hh>
#include <protocols/frag_picker/scores/RamaScore.hh>
#include <protocols/frag_picker/scores/CSScore.hh>
#include <protocols/frag_picker/scores/ABEGO_SS_Score.hh>
#include <protocols/frag_picker/scores/TorsionBinSimilarity.hh>
#include <protocols/frag_picker/scores/FragmentScoreMap.hh>
#include <protocols/frag_picker/TorsionBinIO.hh>

#include <core/fragment/FragID.hh> // required for windows build
#include <core/fragment/SecondaryStructure.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <core/fragment/util.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/frags.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/sequence/util.hh>
#include <basic/prof.hh>
#include <basic/Tracer.hh>

#include <utility/exit.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <utility>
#include <sstream>
#include <iostream>
#include <fstream>

#include <core/fragment/ConstantLengthFragSet.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace frag_picker {

using namespace core;
using namespace core::fragment;
using namespace protocols::frag_picker;
using namespace protocols::frag_picker::scores;
using namespace basic::options;
using namespace basic::options::OptionKeys;

static basic::Tracer trPicker("protocols.frag_picker.FragmentPicker");

FragmentPicker::~FragmentPicker() {}

void FragmentPicker::bounded_protocol() {
	pick_candidates();
	save_fragments();
}

void FragmentPicker::quota_protocol() {
	using namespace ObjexxFCL;

	pick_candidates();

	FragmentScoreManagerOP ms = get_score_manager();
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		std::string out_file_name = prefix_ + "." + string_of(fragment_size)
			+ "mers";
		utility::io::ozstream output(out_file_name);
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);
		quota::QuotaCollector *c =
			dynamic_cast<quota::QuotaCollector*> (storage());
		if (c == 0) {
			utility_exit_with_message("Cant' cast candidates' collector to QuotaCollector. Is quota set up correctly?");
		}

		log_25_.setup_summary(c);
		log_200_.setup_summary(c);

		std::ofstream out_file;
		if (option[frags::describe_fragments].user()) {
			std::string describe_name = option[frags::describe_fragments]()+"."+string_of(fragment_size)+"mers";
			out_file.open(describe_name.c_str());
		}

		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {

			Size qPos = query_positions_[iqpos];
			if ( qPos > size_of_query() - fragment_size + 1 ) continue;

			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > dummy_input;	// we don't need them, but some parameter has to be passed to the method
			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > out;

			quota::QuotaSelector selector(n_frags_,qPos, c );

			selector.select_fragments(dummy_input,out);
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,
				out.size()) << std::endl << std::endl;
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
			if( ms->if_late_scoring_for_zeros() )  {
				for (Size fi = 1; fi <= out.size(); ++fi)
					ms->score_zero_scores(out[fi].first,out[fi].second);
			}
			if (option[frags::describe_fragments].user()) {
				ms->describe_fragments(out, out_file);
			}
		}
		log_25_.write_summary();
		log_200_.write_summary();
		//		storage->print_report(trPicker.Debug, ms);
		output.close();
	}
}

void FragmentPicker::keep_all_protocol() {
	using namespace ObjexxFCL;

	CompareTotalScore comparator(get_score_manager());
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];

		std::ofstream out_file;
		if (option[frags::describe_fragments].user()) {
			std::string describe_name = option[frags::describe_fragments]()+"."+string_of(fragment_size)+"mers";
			out_file.open(describe_name.c_str());
		}

		std::string out_file_name = prefix_ + "." + string_of(fragment_size)
			+ "mers";
		utility::io::ozstream output(out_file_name);
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);

		for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) {

			Size qPos = query_positions_[iqpos];
			if ( qPos > size_of_query() - fragment_size + 1 ) continue;

			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> >	out;
			pick_candidates(qPos,fragment_size);
			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> >  candidates = storage->get_candidates(qPos);
			std::sort(candidates.begin(),candidates.end(),comparator);
			selector_->select_fragments(candidates, out);
			if(out.size() == 0) continue;
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,
				out.size()) << std::endl << std::endl;
			FragmentScoreManagerOP ms = get_score_manager();
			if( ms->if_late_scoring_for_zeros() )  {
				for (Size fi = 1; fi <= out.size(); ++fi)
					ms->score_zero_scores(out[fi].first,out[fi].second);
			}
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
			if (option[frags::describe_fragments].user()) {
				get_score_manager()->describe_fragments(out, out_file);
			}
			trPicker.Info << "Collected candidates of size "<<fragment_size<<" at pos"<<qPos<<std::endl;
			storage->clear();
			trPicker.Debug<< storage->count_candidates()<<" candidates left in a sink after flushing"<<std::endl;
		}
		output.close();
		out_file.close();
	}
	trPicker.Info<<std::endl;
}

void FragmentPicker::nonlocal_pairs_protocol() {
	using namespace ObjexxFCL;

	typedef std::pair<Size,Size> pos_pair;

	// always print ca coords
	bool orig_opt = option[frags::write_ca_coordinates]();
	option[frags::write_ca_coordinates].value(true);

	// defaults
	Size min_seq_sep = 12;
	Size ca_dist_squared = 100;
	Size min_contacts_per_res =  1;
	// non-local contact definition options
	if (option[ frags::nonlocal::min_seq_sep ].user()) {
		min_seq_sep = option[ frags::nonlocal::min_seq_sep ]();
	}
	if (option[ frags::nonlocal::ca_dist ].user()) {
		Size min_dist = option[ frags::nonlocal::ca_dist ]();
		ca_dist_squared = min_dist*min_dist;
	}
	if (option[ frags::nonlocal::min_contacts_per_res ].user()) {
		min_contacts_per_res = option[ frags::nonlocal::min_contacts_per_res ]();
	}
	// for valid interacting pair
	/*
		Real max_rmsd_after_relax =  1.5;
		Real max_ddg_score = -4.0;
		if (option[ frags::nonlocal::max_ddg_score ].user()) {
		max_ddg_score = option[ frags::nonlocal::max_ddg_score ]();
		}
		if (option[ frags::nonlocal::max_rmsd_after_relax ].user()) {
		max_rmsd_after_relax = option[ frags::nonlocal::max_rmsd_after_relax ]();
		}
	*/

	std::string query_seq = get_query_seq_string();
	utility::vector1<char> query_sequence;
  for (Size i=0;i<query_seq.length();++i) query_sequence.push_back(query_seq[i]);

	// use quota protocol to select initial single fragment candidates
	// pick candidates as normal
	pick_candidates();

	FragmentScoreManagerOP ms = get_score_manager();
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		Size min_pdb_seq_sep = fragment_size+min_seq_sep;
		Size min_contacts = min_contacts_per_res*fragment_size;
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);
		quota::QuotaCollector *c =
			dynamic_cast<quota::QuotaCollector*> (storage());
		CandidatesCollectorOP storagej = get_candidates_collector(fragment_size);
		if (c == 0)
			utility_exit_with_message("Cant' cast candidates' collector to QuotaCollector. Is quota set up correctly?");
		Size maxiqpos = query_positions_.size()-min_seq_sep-fragment_size-fragment_size+1;
		Size maxjqpos = query_positions_.size()-fragment_size+1;

		std::string out_file_name = prefix_ + "." + string_of(fragment_size) + "mers.nonlocal_pairs";
    std::string contacts_out_file_name = prefix_ + "." + string_of(fragment_size) + "mers.nonlocal_pairs.contacts";
    const std::string silent_out_file_name = prefix_ + "." + string_of(fragment_size) + "mers.nonlocal_pairs.out";
		core::io::silent::SilentFileData sfd;

		utility::io::ozstream output(out_file_name);
		// double(...) is added for MSVC compatability
		output << "# min_seq_sep: " << min_seq_sep << " ca_dist: " << sqrt( double(ca_dist_squared) ) << " min_contacts: " << min_contacts << std::endl << std::endl;
		log_25_.setup_summary(c);
		log_200_.setup_summary(c);

		std::map<pos_pair,Size> contacts_map;
		std::map<pos_pair,Size>::iterator contacts_map_it;
		// loop through query positions, qPosi
		utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > dummy_input;	// we don't need them, but some parameter has to be passed to the method
		for (Size iqpos = 1; iqpos <= maxiqpos; ++iqpos) {
			Size qPosi = query_positions_[iqpos];
			trPicker.Info << "Query position: " << qPosi << std::endl;
			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > outi;
			quota::QuotaSelector selectori(n_frags_,qPosi, c );
			selectori.select_fragments(dummy_input,outi);
			// loop through nonlocal query positions, qPosj
			for (Size jqpos = iqpos+fragment_size+min_seq_sep+1; jqpos <= maxjqpos; ++jqpos) {
				Size qPosj = query_positions_[jqpos];
				utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > outj;
				quota::QuotaSelector selectorj(n_frags_,qPosj, c );
				selectorj.select_fragments(dummy_input,outj);
				for (Size fi = 1; fi <= outi.size(); ++fi) { // loop through selected fragments at qPosi
					for (Size fj = 1; fj <= outj.size(); ++fj) { // loop through selected fragments at qPosj
						if (!outi[fi].first->same_chain( outj[fj].first )) continue; // skip if not from same pdb chain
						//if (outi[fi].first->get_residue(1)->resi() >= outj[fj].first->get_residue(1)->resi()) continue; // skip inverse pairs
						if (std::abs(int(outi[fi].first->get_residue(1)->resi()-outj[fj].first->get_residue(1)->resi())) < min_pdb_seq_sep) continue; // skip if too local in PDB
						Size contacts = 0;
						Size qpi = qPosi;
						for (Size i=1; i<=fragment_size;++i) {
							VallResidueOP ri = outi[fi].first->get_residue(i);
							Size qpj = qPosj;
							for (Size j=1; j<=fragment_size;++j) {
								if (ri->distance_squared(outj[fj].first->get_residue(j)) < ca_dist_squared) {
									pos_pair contact = std::make_pair(qpi,qpj);
									contacts_map[contact]++;
									contacts++;
								}
								qpj++;
							}
							qpi++;
						}
						if (contacts >= min_contacts) {

							trPicker.Info <<	"	Pair " << qPosi << " " << qPosj << " " << outi[fi].first->get_residue(1)->resi() <<
								" " << outj[fj].first->get_residue(1)->resi() << " " << outi[fi].first->get_pdb_id() << outi[fi].first->get_chain_id() << " " << contacts << std::endl;

							output << "pair: " << qPosi << " " << qPosj << " " << contacts << std::endl << std::endl;
							outi[fi].first->print_fragment(output);
							output << std::endl;
							outj[fj].first->print_fragment(output);
							output << std::endl;

							// make pose from frag
							utility::vector1<fragment::FragDataCOP> fragdatapair;
							fragdatapair.push_back(outi[fi].first->get_frag_data());
							fragdatapair.push_back(outj[fj].first->get_frag_data());
							pose::Pose pose;
							std::string sequence;    // = outi[fi].first->sequence() + outj[fj].first->sequence();
							for (Size i=0; i<fragment_size;++i) sequence += query_sequence[qPosi+i];
							for (Size i=0; i<fragment_size;++i) sequence += query_sequence[qPosj+i];

							fragment::make_pose_from_frags( pose, sequence, fragdatapair, false );
							/*
								std::stringstream outputpdb;
								outputpdb << "frags_pose_dump_" << fragment_size << "_" << qPosi << "_" << qPosj << "_" << outi[fi].first->get_residue(1)->resi() <<
								"_" << outj[fj].first->get_residue(1)->resi() << "_" << outi[fi].first->get_pdb_id() << outi[fi].first->get_chain_id()	<< ".pdb";
								pose.dump_pdb(outputpdb.str());
							*/
							core::io::silent::SilentStructOP ss = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct_out( pose );
							std::stringstream tag;
							tag << outi[fi].first->get_pdb_id() << outi[fi].first->get_chain_id() << "_" << qPosi << "_" << qPosj <<
								"_" << outi[fi].first->get_residue(1)->resi() << "_" << outj[fj].first->get_residue(1)->resi();
							ss->fill_struct( pose, tag.str() );
							sfd.write_silent_struct( *ss, silent_out_file_name );

						}

					}
				}
			}
		}
		output.close();

		// save pair contact counts
		utility::io::ozstream contacts_output(contacts_out_file_name);
		// double(...) is added for MSVC compatability
		contacts_output << "# min_seq_sep: " << min_seq_sep << " ca_dist: " << sqrt( double(ca_dist_squared) ) << " min_contacts: " << min_contacts << std::endl << std::endl;
		for ( contacts_map_it=contacts_map.begin() ; contacts_map_it != contacts_map.end(); contacts_map_it++ ) {
			pos_pair p = (*contacts_map_it).first;
			contacts_output << p.first << " " << p.second << " " << (*contacts_map_it).second << std::endl;
		}
		contacts_output.close();

	}
	option[frags::write_ca_coordinates].value(orig_opt);
}

void FragmentPicker::pick_candidates() {

	PROF_START( basic::FRAGMENTPICKING );
	// Size n_total = chunks_->size()/100;

	scores::FragmentScoreMapOP empty_map = scores_->create_empty_map();
	for (Size i = 1; i <= chunks_->size(); i++) { // loop over provided chunks
		VallChunkOP chunk = chunks_->at(i); // For each chunk from a provider...
		if (chunk->size() < max_frag_size_) // This fragment is too short
			continue;
		bool flag = true;
		for (Size iFilter = 1; iFilter <= filters_.size(); iFilter++) {
			if ((flag = filters_[iFilter]->test_chunk(chunk)) == false) {
				trPicker.Debug << "Chunk: " << chunk->get_pdb_id()
											 << " didn't pass a filter" << std::endl;
				break;
			}
		}
		if (!flag)
			continue;
		trPicker.Debug << "Processing sequence from vall: "
									 << chunk->get_sequence() << std::endl;

		// cache the new chunk
		scores_->do_caching(chunk);

		for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
			Size fragment_size = frag_sizes_[iFragSize];
			if (chunk->size() < fragment_size) // This fragment is too short
				continue;

			CandidatesCollectorOP sink = candidates_sink_[fragment_size];
			trPicker.Debug << "Picking fragments of size "<<fragment_size<<
				" at "<<query_positions_.size()<<" query positions"<<std::endl;
			for (Size iqpos = 1; iqpos <= query_positions_.size(); ++iqpos) { // loop over positions in a query

				Size iPos = query_positions_[iqpos];
				if ( iPos > size_of_query() - fragment_size + 1 ) continue;

				// split chunk into fragment candidates and score them
				for (Size j = 1; j <= chunk->size() - fragment_size + 1; j++) {
					FragmentCandidateOP f = new FragmentCandidate(iPos, j,
						chunk, fragment_size);
					if (scores_->score_fragment_from_cache(f, empty_map)) {
						std::pair<FragmentCandidateOP,scores::FragmentScoreMapOP> p(f,empty_map);
						if(sink->add(p))
							empty_map  = scores_->create_empty_map();
						/*
							scores::FragmentScoreMapOP new_map = empty_map->clone();
							std::pair<FragmentCandidateOP,
							scores::FragmentScoreMapOP> p(f, new_map);
							sink->add(p);
						*/
					}
				}
			} // all query positions done
		} // all fragment sizes done
		scores_->clean_up();
		trPicker.Debug << chunk->get_pdb_id() << " done" << std::endl;
		if( (i*100) % chunks_->size() == 0 ) trPicker.Info << (i*100) / chunks_->size()
																											 << "% done at "<<chunk->get_pdb_id()<< std::endl;
		trPicker.flush();
	} // all chunks done
	PROF_STOP( basic::FRAGMENTPICKING );
}

double FragmentPicker::total_score(scores::FragmentScoreMapOP f) {

	utility::vector1<Real> components = f->get_score_components();
	utility::vector1<Real> weights = scores_->get_weights();
	Real total = 0.0;
	for (Size i = 1; i <= components.size(); i++)
		total += components.at(i) * weights.at(i);

	return total;
}


void FragmentPicker::read_ss_files(utility::vector1<std::string> sec_str_input) {

	trPicker.Debug << sec_str_input.size() / 2
								 << " secondary structure assignment(s):\n";
	for (Size i = 1; i <= sec_str_input.size(); i += 2) {
		trPicker.Debug << i / 2 << " " << sec_str_input[i]
									 << " file will be loaded under \"" << sec_str_input[i + 1]
									 << "\" name\n";
		read_ss_file(sec_str_input[i], sec_str_input[i + 1]);
	}
	trPicker.Debug << std::endl;
}

void FragmentPicker::read_ss_file(std::string const & file_name,
	std::string prediction_name) {

	utility::io::izstream data( file_name.c_str() );
	if ( !data ) {
		data.close();
		utility_exit_with_message( "Can't read secondary structure file: "+file_name );
	}

	std::string line, l1, l2, l3, l4, l5;
	getline( data, line );
	data.close();

	std::istringstream line_stream( line );
	line_stream >> l1 >> l2 >> l3 >> l4 >> l5;

	if ( (l1 == "#") && (l2 == "PSIPRED") && (l3 == "VFORMAT")
		&& (l4 == "(PSIPRED") ) {
		read_psipred_ss2( file_name, prediction_name);
	} else {
		if ( (l1 == "REMARK") && (l2 == "Neural") && (l3 == "network")
			&& (l4 == "secondary") && (l5 == "structure") ) {
			read_talos_ss( file_name, prediction_name);
		} else {
			utility_exit_with_message( "Can't identify secondary structure file type (needs vertical psipred_ss2 or talos+ pred.ss): "+file_name );
		}
	}
}

void FragmentPicker::read_psipred_ss2(std::string const & file_name,
	std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile =
		new core::fragment::SecondaryStructure();
	ss_profile->read_psipred_ss2(file_name);

	std::string query_ss_as_string;
	for (Size i = 1; i <= ss_profile->total_residue(); i++)
		query_ss_as_string += ss_profile->secstruct(i);

	query_ss_as_string_[prediction_name] = query_ss_as_string;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::read_talos_ss(std::string const & file_name,
	std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile =
		new core::fragment::SecondaryStructure();
	ss_profile->read_talos_ss(file_name);

	std::string query_ss_as_string;
	for (Size i = 1; i <= ss_profile->total_residue(); i++)
		query_ss_as_string += ss_profile->secstruct(i);

	query_ss_as_string_[prediction_name] = query_ss_as_string;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::add_query_ss(std::string query_secondary,
	std::string prediction_name) {

	core::fragment::SecondaryStructureOP ss_profile =
		new core::fragment::SecondaryStructure();
	ss_profile->extend(query_secondary.length());

	for (Size i = 1; i <= query_secondary.length(); ++i) {
		char ss = query_secondary[i - 1];
		if (ss == 'E')
			ss_profile->set_fractions(i, 0.0, 1.0, 0.0);
		else if (ss == 'L')
			ss_profile->set_fractions(i, 0.0, 0.0, 1.0);
		else
			ss_profile->set_fractions(i, 1.0, 0.0, 0.0);
	}
	query_ss_as_string_[prediction_name] = query_secondary;
	query_ss_profile_[prediction_name] = ss_profile;
}

void FragmentPicker::save_fragments() {
	using namespace ObjexxFCL;

	FragmentScoreManagerOP ms = get_score_manager();
	CompareTotalScore comparator(ms);
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments

		Size fragment_size = frag_sizes_[iFragSize];

		std::ofstream out_file;
		if (option[frags::describe_fragments].user()) {
			std::string describe_name = option[frags::describe_fragments]()+"."+string_of(fragment_size)+"mers";
			out_file.open(describe_name.c_str());
		}

		std::string out_file_name = prefix_ + "." + string_of(fragment_size)
			+ "mers";
		utility::io::ozstream output(out_file_name);
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);
		for (Size qPos = 1; qPos <= size_of_query(); ++qPos) {
			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> >
				out;
			if(storage->get_candidates(qPos).size() == 0) continue;
			selector_->select_fragments(storage->get_candidates(qPos), out);
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,
				out.size()) << std::endl << std::endl;
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
			if( ms->if_late_scoring_for_zeros() )  {
				for (Size fi = 1; fi <= out.size(); ++fi)
					ms->score_zero_scores(out[fi].first,out[fi].second);
			}
			if (option[frags::describe_fragments].user()) {
				ms->describe_fragments(out, out_file);
			}
		}
		storage->print_report(trPicker.Info, get_score_manager());
		output.close();
		out_file.close();
	}
}

void FragmentPicker::save_candidates() {
	using namespace ObjexxFCL;

	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) { // Loop over various sizes of fragments
		Size fragment_size = frag_sizes_[iFragSize];
		std::string out_file_name = prefix_ + "." + string_of(fragment_size)
			+ "mers";
		utility::io::ozstream output(out_file_name);
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);

		for (Size qPos = 1; qPos <= size_of_query(); ++qPos) {
			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> >
				out(storage->get_candidates(qPos));
			output << "position: " << I(12, qPos) << " neighbors:   " << I(10,
				out.size()) << std::endl << std::endl;
			for (Size fi = 1; fi <= out.size(); ++fi) {
				out[fi].first->print_fragment(output);
				output << std::endl;
			}
		}
		storage->print_report(trPicker.Debug, get_score_manager());
		output.close();
	}
}

void FragmentPicker::pick_candidates(Size i_pos,Size frag_len) {

	scores::FragmentScoreMapOP empty_map = scores_->create_empty_map();
	for (Size i = 1; i <= chunks_->size(); i++) { // loop over provided chunks
		VallChunkOP chunk = chunks_->at(i); // For each chunk from a provider...
		if (chunk->size() < frag_len) // This fragment is too short
			continue;
		bool flag = true;
		for (Size iFilter = 1; iFilter <= filters_.size(); iFilter++) {
			if ((flag = filters_[iFilter]->test_chunk(chunk)) == false) {
				trPicker.Debug << "Chunk: " << chunk->get_pdb_id()
											 << " didn't pass a filter" << std::endl;
				break;
			}
		}
		if (!flag)
			continue;
		trPicker.Debug << "Processing sequence from vall: "
									 << chunk->get_sequence() << std::endl;

		CandidatesCollectorOP sink = candidates_sink_[frag_len];

		// split chunk into fragment candidates and score them
		for (Size j = 1; j <= chunk->size() - frag_len + 1; j++) {
			FragmentCandidateOP f = new FragmentCandidate(i_pos, j,
				chunk, frag_len);
			if (scores_->score_fragment(f, empty_map)) {
				scores::FragmentScoreMapOP new_map = empty_map->clone();
				std::pair<FragmentCandidateOP,
					scores::FragmentScoreMapOP> p(f, new_map);
				sink->add(p);
			}
		} // All chunk locations done
		trPicker.Debug << chunk->get_pdb_id() << " done" << std::endl;
		trPicker.Debug << sink->count_candidates()<<" candidates stored at pos. "
									 <<i_pos<<", "<<sink->count_candidates()<<" in total"<< std::endl;
		trPicker.flush();
	} // all chunks done
}


void FragmentPicker::parse_command_line() {

	//## -------- setup query profile
	if (option[in::file::checkpoint].user()) {
		core::sequence::SequenceProfileOP q_prof(
																						 new core::sequence::SequenceProfile);
		trPicker.Info << "reading a query profile from: "
									<< option[in::file::checkpoint]() << std::endl;
		q_prof->read_from_checkpoint(option[in::file::checkpoint]());
		set_query_seq(q_prof);
		trPicker.Info << "picking fragments for query profile: "
									<< get_query_seq_string() << std::endl;
	}
	if (option[in::file::pssm].user()) {
		core::sequence::SequenceProfileOP q_prof(
																						 new core::sequence::SequenceProfile);
		trPicker.Info << "reading a query profile from: "
									<< option[in::file::pssm]()[1] << std::endl;
		q_prof->read_from_file(option[in::file::pssm]()[1], 1.0);
		set_query_seq(q_prof);
		trPicker.Info << "picking fragments for query profile: "
									<< get_query_seq_string() << std::endl;
	}

	//Fasta file trumps sequence profile as far as query_seq_string_ is concerned
	if (option[in::file::fasta].user()) {
		std::string q_seq = core::sequence::read_fasta_file(
																												option[in::file::fasta]()[1])[1]->sequence();
		trPicker.Info << "reading a query sequence from: "
									<< option[in::file::fasta]()[1] << std::endl;

		set_query_seq(q_seq);
		trPicker.Info << "picking fragments for query sequence: "
									<< get_query_seq_string() << std::endl;
	}

	// --------- setup query secondary structure
	if (option[frags::ss_pred].user()) {
		utility::vector1<std::string> sec_str_input(option[frags::ss_pred]());
		read_ss_files(sec_str_input);
	}
	//---------- setup chunk filters
	if (option[frags::allowed_pdb].user()) {
		AllowPdbIdFilterOP allow = new AllowPdbIdFilter();
		allow->load_pdb_id_from_file(option[frags::allowed_pdb]());
		add_chunk_filter(allow);
		trPicker.Info << "Allowed PDB chains:\n";
		allow->show_pdb_ids(trPicker.Info);
	}

	if (option[frags::denied_pdb].user()) {
		DenyPdbIdFilterOP deny = new DenyPdbIdFilter();
		deny->load_pdb_id_from_file(option[frags::denied_pdb]());
		add_chunk_filter(deny);
		trPicker.Info << "Excluded PDB chains:\n";
		deny->show_pdb_ids(trPicker.Info);
	}

	// ##--------- setup VALL
	PROF_START( basic::FRAGMENTPICKING_READ_VALL );
	if (option[in::file::vall].user()) {
		read_vall(option[in::file::vall]());
	}
	PROF_STOP( basic::FRAGMENTPICKING_READ_VALL );

	// -------- fragment sizes
	if (option[frags::frag_sizes].user()) {
		utility::vector1<Size> frag_sizes_tmp = option[frags::frag_sizes]();
		for (Size i = 1; i <= frag_sizes_tmp.size(); ++i) {
			if(frag_sizes_tmp[i] > max_frag_size_)
				max_frag_size_ = frag_sizes_tmp[i];
			frag_sizes_.push_back(frag_sizes_tmp[i]);
		}
	} else {
		max_frag_size_ = 9;
		frag_sizes_.push_back(3);
		frag_sizes_.push_back(9);
	}
	trPicker.Info << "Will pick fragments of size:";
	for (Size i = 1; i <= frag_sizes_.size(); ++i)
		trPicker.Info << frag_sizes_[i] << " ";
	trPicker.Info << std::endl;

	//---------- setup scoring scheme
	trPicker.Info << "Creating fragment scoring scheme" << std::endl;
	FragmentScoreManagerOP scoring = get_score_manager();
	if (option[frags::scoring::config].user()) {
		scoring->create_scores(option[frags::scoring::config](), this);
	}

	// -------- how many fragments and candidates
	n_frags_ = option[frags::n_frags]();
	n_candidates_ = option[frags::n_candidates]();

	if (n_frags_ > n_candidates_) n_candidates_ = n_frags_;

	trPicker.Info << "Picking " << n_frags_ << " fragments based on "
								<< n_candidates_ << " candidates" << std::endl;

	//-------- this comparator is used both for collecting and selecting fragments
	CompareTotalScore comparator(get_score_manager());


	//---------- setup scoring scheme for the selection step
	trPicker.Info << "Creating fragment scoring scheme for the selection step" << std::endl;
	FragmentScoreManagerOP selection_scoring;
	if (option[frags::picking::selecting_scorefxn].user()) {
		selection_scoring  = new FragmentScoreManager();
		selection_scoring->create_scores(option[frags::picking::selecting_scorefxn](), this);
		selector_ = new CustomScoreSelector(n_frags_, selection_scoring);
	} else {
		selector_ = new BestTotalScoreSelector(n_frags_, scoring);
	}

	//-------- collector & selector set up
	if (option[frags::quota_protocol].user() || option[frags::picking::quota_config_file].user()) {
		// This setup is a bit more complicated when user needs quota.
		// The quota version of this code was moved into a separate method
		parse_quota_command_line();
		// This setup is a bit more complicated, when user needs quota. The quota version of this code was moved into a separate method
	} else {
		if (option[frags::keep_all_protocol].user()) {
			for (Size i = 1; i <= frag_sizes_.size(); ++i) {
		    CandidatesCollectorOP collector = new GrabAllCollector(size_of_query());
		    set_candidates_collector(frag_sizes_[i], collector);
		    trPicker.Info << "Collector for fragment size: " << frag_sizes_[i]
											<< " set to: GrabAllCollector" << std::endl;
			}
		} else {
			for (Size i = 1; i <= frag_sizes_.size(); ++i) {
				CandidatesCollectorOP collector = new BoundedCollector<
				CompareTotalScore> (size_of_query(), n_candidates_,
					comparator,get_score_manager()->count_components());
		    set_candidates_collector(frag_sizes_[i], collector);
		    trPicker.Info << "Collector for fragment size: " << frag_sizes_[i]
											<< " set to: BoundedCollector" << std::endl;
			}
		}
		//-------- Selecting fragments from candidates
		/*	    if (option[frags::picking::selecting_rule].user()) {
						std::string type = option[frags::picking::selecting_rule]();
						if (type.compare("BestTotalScoreSelector")==0) {
						selector_ = new BestTotalScoreSelector(n_frags_, selection_scoring);
						trPicker.Info << "Fragment selector: BestTotalScoreSelector"
						<< std::endl;
						} else {
						utility_exit_with_message("[ERROR]: unknown fragment selecting rule: " + type + "!");
						}
						} else {
						selector_ = new BestTotalScoreSelector(n_frags_, selection_scoring);
						trPicker.Info << "Fragment selector: BestTotalScoreSelector" << std::endl;
						}*/
	}
	// # ---------- output file prefix:
	if (option[out::file::frag_prefix].user()) {
		prefix_ = option[out::file::frag_prefix]();
	}

	if (option[frags::picking::query_pos].user()) {
		set_picked_positions( option[frags::picking::query_pos]() );
	}

	show_scoring_methods(trPicker);
	trPicker << std::endl;
}

void FragmentPicker::set_up_ss_abego_quota() {

	std::string quota_config_file("UNKNOWN-QUOTA-CONFIG_FILE");
	if (option[frags::picking::quota_config_file].user())
		quota_config_file = option[frags::picking::quota_config_file]();
	quota::ABEGO_SS_Config q_config(quota_config_file);

	utility::vector1<Size> components;
	utility::vector1<Real> weights;
	utility::vector1<Real> scoring_weights = scores_->get_weights();
	for (Size i = 1; i <= scores_->count_components(); ++i) {
		ABEGO_SS_Score *s0 =
			dynamic_cast<ABEGO_SS_Score*> (scores_->get_component(i).get());
		if (s0 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s0->get_id()] );
		}
		ProfileScoreL1 *s1 =
			dynamic_cast<ProfileScoreL1*> (scores_->get_component(i).get());
		if (s1 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s1->get_id()] );
		}

		RamaScore *s2 =
			dynamic_cast<RamaScore*> (scores_->get_component(i).get());
		if (s2 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s2->get_id()] );
		}

		CSScore *s3 =
			dynamic_cast<CSScore*> (scores_->get_component(i).get());
		if (s3 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s3->get_id()] );
		}

		SecondarySimilarity *s4 =
			dynamic_cast<SecondarySimilarity*> (scores_->get_component(i).get());
		if (s4 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s4->get_id()] );
		}
	}

	trPicker.Debug<<"Scoring scheme for ABEGO_SS quota pool sorting is:";
	for(Size l=1;l<=weights.size();l++)
		trPicker.Debug<<"\n\t"<<components[l]<<"\t"<<weights[l];
	trPicker.Debug<<std::endl;
	Size buffer_factor = 5;
	for(Size f=1;f<=frag_sizes_.size();f++) {
		quota::QuotaCollectorOP collector = new quota::QuotaCollector( size_of_query(), frag_sizes_[f] );
		set_candidates_collector(frag_sizes_[f],collector);
		Size middle = frag_sizes_[f] / 2 + 1;
		assert( size_of_query() == q_config.size() ); // Test if the abego-ss table has the same size as the query sequence
		for(Size j=1;j<=size_of_query()-frag_sizes_[f]+1;j++) {
			trPicker.Trace<<"Creating "<<q_config.n_columns()<<" quota pools at pos "<<j<<std::endl;
			for(Size i=1;i<=q_config.n_columns();i++) {
		    Real f = q_config.probability(j+middle-1,i);
		    quota::QuotaPoolOP p = new quota::ABEGO_SS_Pool(n_candidates_,q_config.get_pool_name(i),
			    q_config.get_pool_bins((i)),components,weights,f,scores_->count_components(),buffer_factor);
		    collector->add_pool(j,p);
			}
		}
	}
}

void FragmentPicker::set_up_quota_nnmake_style() {

	std::string quota_config_file("UNKNOWN-QUOTA-CONFIG_FILE");
	if (option[frags::picking::quota_config_file].user())
		quota_config_file = option[frags::picking::quota_config_file]();
	quota::QuotaConfig q_config(quota_config_file);

	utility::vector1<Size> components;
	utility::vector1<Real> weights;
	components.push_back( 0 );		// the free entry in the vector is for secondary structure score (only one for each pool)
	weights.push_back( 0.0 );		// score weight for SecondarySimilarity; will be changed later
	components.push_back( 0 );		// this free entry in the vector is for RamaScore
	weights.push_back( 0.0 );		// score weight for RamaScore; will be changed later
	utility::vector1<Real> scoring_weights = scores_->get_weights();
	for (Size i = 1; i <= scores_->count_components(); ++i) {
		ProfileScoreL1 *s1 =
			dynamic_cast<ProfileScoreL1*> (scores_->get_component(i).get());
		if (s1 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s1->get_id()] );
		}

		/****** RamaScore is a special case, dispatched below
						RamaScore *s2 =
						dynamic_cast<RamaScore*> (scores_->get_component(i).get());
						if (s2 != 0) {
						components.push_back( i );
						weights.push_back( scoring_weights[s2->get_id()] );
						}
		*********/
		CSScore *s3 =
			dynamic_cast<CSScore*> (scores_->get_component(i).get());
		if (s3 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s3->get_id()] );
		}
		ABEGO_SS_Score *s4 =
			dynamic_cast<ABEGO_SS_Score*> (scores_->get_component(i).get());
		if (s4 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s4->get_id()] );
		}
		TorsionBinSimilarity *s5 =
			dynamic_cast<TorsionBinSimilarity*> (scores_->get_component(i).get());
		if (s5 != 0) {
			components.push_back( i );
			weights.push_back( scoring_weights[s5->get_id()] );
		}
	}

	utility::vector1<core::fragment::SecondaryStructureOP> predictions;
	utility::vector1<Real> ss_weights;
	std::map<std::string, core::fragment::SecondaryStructureOP>::iterator it;
	Real weight = 1.0 / ((Real) query_ss_profile_.size());
	for ( it=query_ss_profile_.begin() ; it != query_ss_profile_.end(); it++ ) {
		predictions.push_back((*it).second);
		ss_weights.push_back(weight);
	}
	//	core::fragment::SecondaryStructureOP avg_ss = new core::fragment::SecondaryStructure(predictions,ss_weights);

	for(Size f=1;f<=frag_sizes_.size();f++) {
		quota::QuotaCollectorOP collector = new quota::QuotaCollector( size_of_query(), frag_sizes_[f] );
		set_candidates_collector(frag_sizes_[f],collector);

		// --------- This part puts RamaScore into quota scoring; each Rama is based on a certain SS prediction and this part of the code
		// --------- dispatches each Rama into a proper pool
		for (Size i = 1; i <= scores_->count_components(); ++i) {
			RamaScore *sr =
				dynamic_cast<RamaScore*> (scores_->get_component(i).get());
			if (sr != 0) {
				std::string & name = sr->get_prediction_name();
				if( ! q_config.is_valid_quota_pool_name( name ) )
					continue;
				components[2] = i;
				weights[2] = scoring_weights[sr->get_id()];
				trPicker.Warning<<"RamaScore with ID "<<sr->get_id()<<" named "<<name<<
					" has been attached to its quota pool with weight "<<weights[2]<<std::endl;
			}
		}
		// ---------- end of RamaScore dispatch

		// Create secondary structure pools (if any)
		for (Size i = 1; i <= scores_->count_components(); ++i) {
			SecondarySimilarity *s =
				dynamic_cast<SecondarySimilarity*> (scores_->get_component(i).get());

			//PartialSecondarySimilarity is a variant of SecondarySimilarity, this means they're not compatible
			//So what is it compatible with???
			//		if (s == 0) {
			//			PartialSecondarySimilarity *s =
			//				dynamic_cast<PartialSecondarySimilarity*> (scores_->get_component(i).get());
			//		}

			if (s != 0) {
				std::string & name = s->get_prediction_name();
				if( ! q_config.is_valid_quota_pool_name( name ) )
					continue;
				components[1] = i;
				weights[1] = scoring_weights[s->get_id()];
				Size size = (Size)(q_config.get_fraction( name ) * n_candidates_);
				if( size == 0 ) {
					trPicker.Warning<<"Config file couldn't provide quota fraction for the pool named "
													<<name<<". Skipping the pool"<<std::endl;
					continue;
				}
				collector->attach_secondary_structure_pools(q_config.get_fraction( name ) ,
			    get_query_ss( name ),name,n_candidates_,components,weights,scores_->count_components());
				//			    avg_ss,name,n_candidates_,components,weights,scores_->count_components());
			}
		}
	}
}

void FragmentPicker::parse_quota_command_line() {

	set_up_quota_nnmake_style();

	if (option[frags::picking::query_pos].user()) {
		set_picked_positions( option[frags::picking::query_pos]() );
	}
}

void FragmentPicker::read_vall( utility::vector1< std::string > const & fns ) {
	chunks_ = new VallProvider();
	chunks_->vallChunksFromLibraries(fns);
}

void FragmentPicker::read_vall( std::string const & fn ) {
	chunks_ = new VallProvider();
	chunks_->vallChunksFromLibrary(fn);
}

void FragmentPicker::set_picked_positions(Size from,Size to) {

	query_positions_.clear();
	for(Size i=from;i<=to;i++)
		query_positions_.push_back( i );
}

void FragmentPicker::set_picked_positions(utility::vector1<Size> q_positions) {

	query_positions_.clear();
	for(Size i=1;i<=q_positions.size();i++)
		query_positions_.push_back( q_positions[i] );

}


Size QuotaDebug::max_pools() {

	return tags_.size();
}

void QuotaDebug::write_summary() {

	trPicker<< "Quota report: difference between the total expected and total picked foreach pool"<<std::endl;
	trPicker<< "This table is for first "<<nFrags_<<" fragments"<<std::endl;
	trPicker<< "Negative value says that was picked more that expected."<<std::endl;
	trPicker<< this->str()<<std::endl;
	this->str("");
}

void QuotaDebug::log(Size frag_len,Size q_pos,utility::vector1<Real> data) {

	*this  << std::setw(4)<<q_pos<< std::setw(4)<<frag_len;
	for(Size i=1;i<=data.size();i++)
		if(data[i]<1000000)
	    *this  << std::setw(10)<<std::setprecision(3)<<data[i];
		else
	    *this  << std::setw(10)<<"  --- ";
	*this  << std::endl;
}

void QuotaDebug::setup_summary(quota::QuotaCollector* collector_) {

	Size last_tag = 0;
	for(Size i=1;i<=collector_->query_length();++i) {
		for(Size j=1;j<=collector_->count_pools(i);++j) {
	    if( tag_map_.find(collector_->get_pool(i,j)->get_pool_name())==tag_map_.end() ) {
				tags_.push_back(collector_->get_pool(i,j)->get_pool_name());
				last_tag++;
				tag_map_[collector_->get_pool(i,j)->get_pool_name()] = last_tag;
	    }
		}
	}

	*this <<"\n#len pos ";
	for(Size i=1;i<=tags_.size();i++) {
		*this << std::setw(10)<<tags_[i];
	}
	*this<<std::endl;
}

utility::vector1<ConstantLengthFragSetOP> FragmentPicker::getFragSet(int residueInPose_){



	// Storage for the resulting FragSets
	utility::vector1<ConstantLengthFragSetOP> result;

  // Loop over various sizes of fragments
	for (Size iFragSize = 1; iFragSize <= frag_sizes_.size(); ++iFragSize) {

		ConstantLengthFragSetOP myFragSet = new ConstantLengthFragSet();

		Size fragment_size = frag_sizes_[iFragSize];
		CandidatesCollectorOP storage = get_candidates_collector(fragment_size);
		for (Size qPos = 1; qPos <= size_of_query(); ++qPos) {
			if(storage->get_candidates(qPos).size() == 0) continue;

			utility::vector1<std::pair<FragmentCandidateOP, FragmentScoreMapOP> > out;
			selector_->select_fragments(storage->get_candidates(qPos), out);

			//FrameOP frame = new Frame(out[1].first->get_residue(1)->resi()); // start pos = residue id from the fragment residue (sequence from original pdb file?) or is this some internal index?

			/*
				start pos = residue id from the fragment residue (sequence from original pdb file?) or is this some internal index?

				In  ConstantLengthFragSet, this is set from insertion_pos of fragment file.

				Therefore,  I think this is the connection between Pose and Fragment File. In design mode,
				we do not have this, therefore we just add an incremented index for each qPos
			*/
			FrameOP frame = new Frame(residueInPose_++);

			for (Size fi = 1; fi <= out.size(); ++fi) {

				FragDataOP current_fragment( NULL );

				for (Size i = 1; i <= out[1].first->get_length(); ++i) {
					VallResidueOP r   =  out[fi].first->get_residue(i);
					string pdbid      = out[fi].first->get_pdb_id();
					char chainid      = out[fi].first->get_chain_id();
					Size index        = r->resi();
					char aa           = toupper(r->aa());
					char ss           = r->ss();
					Real phi          = r->phi();
					Real psi          = r->psi();
					Real omega        = r->omega();

					if (i == 1){
						current_fragment = new AnnotatedFragData( pdbid, index );
					}
					utility::pointer::owning_ptr< BBTorsionSRFD > res_torsions( new BBTorsionSRFD(3,ss,aa) ); // 3 protein torsions
					res_torsions->set_torsion   ( 1, phi   ); // ugly numbers 1-3, but pose.set_phi also uses explicit numbers
					res_torsions->set_torsion   ( 2, psi   );
					res_torsions->set_torsion   ( 3, omega );
					res_torsions->set_secstruct ( ss );

					// Add residue to fragment
					current_fragment->add_residue( res_torsions );

				} // End VallResidue loop

				if (current_fragment) { // != NULL) {
					current_fragment->set_valid(); //it actually containts data

					// Add fragment to frame
					if (!frame->add_fragment(current_fragment)){
						cerr << "ERROR Bad fragment : "<<endl;
						current_fragment->show(cout);
						exit(1111);
					}
				}
			} // End FragmentCandidate loop

			// Add frame to myFragSet.
			myFragSet->add(frame);

		} // End size of query

		// For each size fragment add to vector of ConstantLengthFragSet
		result.push_back(myFragSet);

	} // End size of frags


	return (result);

}


QuotaDebug FragmentPicker::log_25_(25);
QuotaDebug FragmentPicker::log_200_(200);

} // frag_picker
} // protocols
