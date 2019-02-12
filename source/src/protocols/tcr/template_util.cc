// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/template_util.cc
/// @brief template identification utility functions for protocols/tcr/TCRmodel
/// @author Ragul (ragul@umd.edu)

#include <protocols/tcr/template_util.hh>
#include <basic/Tracer.hh>
#include <protocols/tcr/grafting_util.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/tcrmodel.OptionKeys.gen.hh>
#include <set>
#include <utility/file/file_sys_util.hh>

static basic::Tracer TR( "protocols.tcr.template_util" );
/// @brief MINSCORE const value with large negative number
/// @details used for checking and ignoring templates
/// @details also used for initializing template score value
const core::Real MINSCORE = -99999.0;
/// @brief MAXRESOLUTION const value with large possitive number
/// @details used for checking/ignoring templates
/// @details also used for initializing template resolution value
const core::Real MAXRESOLUTION = 10.0;

namespace protocols {
namespace tcr {

bool ignore_template(std::string const &tmpltid, std::list< std::set< std::string > > const &ignore_lists) {
	bool ignore = false;
	for ( std::list< std::set< std::string > > ::const_iterator li = ignore_lists.begin(); li != ignore_lists.end(); ++li ) {
		if ( (li->find(tmpltid) != li->end()) == true ) {
			ignore = true;
			break;
		}
	}
	return ignore;
}

std::pair<core::Real, std::string> find_template( std::string const &query, std::string const &db, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	std::ifstream inp(db);
	if ( !inp.good() ) {
		utility_exit_with_message("Error opening template database : "+db);
	}
	std::string line, name, content;
	std::string best_name, best_content;
	core::Real best_score = MINSCORE;
	core::Real best_resolution = MAXRESOLUTION;
	while ( std::getline(inp, line) ) {
		if ( line.empty() || line[0] == '>' ) {
			if ( !name.empty() ) {
				std::string first_four = name.substr(0, 4);
				if ( !ignore_template(first_four, ignore_lists) ) {
					if ( calculate_identity_score(query,content) <= cutoff ) {
						core::Real curr_resolution = MAXRESOLUTION;
						if ( !(name.substr(9) == "None") ) {
							curr_resolution = std::stof(name.substr(9));
						}
						core::Real curr_score = score_alignment( query, content, tcr_ss );
						if ( curr_score > best_score ) {
							best_score = curr_score;
							best_content = content;
							best_name = name;
							best_resolution = curr_resolution;
						} else if ( curr_score == best_score ) {
							best_score = curr_score;
							best_content = content;
							best_name = name;
							best_resolution = curr_resolution;
						} else if ( curr_score == best_score ) {
							if ( curr_resolution < best_resolution ) {
								best_score = curr_score;
								best_content = content;
								best_name = name;
								best_resolution = curr_resolution;
							}
						}
					}
				}
				name.clear();
			}
			if ( !line.empty() ) {
				name = line.substr(1);
			}
			content.clear();
		} else {
			content += line;
		}
	}
	std::pair<core::Real, std::string> best_pair(best_score, best_name.substr(0,8));
	return best_pair;
}


std::pair<core::Real, std::string> find_template_from_multiple_input_db( std::string const &query, utility::vector1< std::string > const &multidb, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	std::pair<core::Real, std::string> best_pair(MINSCORE, "");
	for ( core::Size i=1; i<=multidb.size(); ++i ) {
		std::pair<core::Real,std::string> curr_best_template = find_template(query, multidb[i], ignore_lists, cutoff, tcr_ss);
		if ( curr_best_template.first > best_pair.first ) {
			best_pair = curr_best_template;
		}
	}
	return best_pair;
}

utility::vector1<std::string> update_template_files(utility::vector1<std::string> const &file_names, std::string const &tcr_seq_db) {
	utility::vector1<std::string> template_db_list;
	std::string template_db_file;
	template_db_list.clear();
	for ( core::Size i=1; i<=file_names.size(); ++i ) {
		template_db_file = tcr_seq_db + file_names[i];
		if ( utility::file::file_exists(template_db_file) ) {
			template_db_list.push_back(template_db_file);
		}
	}
	return template_db_list;
}

void initialize_template_db_files(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, TCRmodel::tcrtmplts &atmplt, TCRmodel::tcrtmplts &btmplt, std::string const &tcr_seq_db) {
	struct {
		utility::vector1< std::string > name;
		utility::vector1< std::string > &value;
	} J[] {
	{{"A_TCR_GM.fasta"}, atmplt.gm.tdb},
	{{"B_TCR_GM.fasta"}, btmplt.gm.tdb},
	{{"A_TCR_FW.fasta"}, atmplt.fr.tdb},
	{{"B_TCR_FW.fasta"}, btmplt.fr.tdb},
	{{"TCR_FW_ORIENTATION.seq"}, atmplt.ori.tdb},
	{{"TCR_FW_ORIENTATION.seq"}, btmplt.ori.tdb},
	{{"A_TCR_CDR1_" + std::to_string(atcr.cdr1.size()) + ".fasta", "B_TCR_CDR1_" + std::to_string(atcr.cdr1.size()) + ".fasta"}, atmplt.cdr1.tdb},
	{{"A_TCR_CDR1_" + std::to_string(btcr.cdr1.size()) + ".fasta", "B_TCR_CDR1_" + std::to_string(btcr.cdr1.size()) + ".fasta"}, btmplt.cdr1.tdb},
	{{"A_TCR_CDR3_" + std::to_string(atcr.cdr3.size()) + ".fasta", "B_TCR_CDR3_" + std::to_string(atcr.cdr3.size()) + ".fasta"}, atmplt.cdr3.tdb},
	{{"A_TCR_CDR3_" + std::to_string(btcr.cdr3.size()) + ".fasta", "B_TCR_CDR3_" + std::to_string(btcr.cdr3.size()) + ".fasta"}, btmplt.cdr3.tdb},
	{{"A_TCR_CDR2HV4_" + std::to_string(atcr.cdr2.size()) + ".fasta", "B_TCR_CDR2HV4_" + std::to_string(atcr.cdr2.size()) + ".fasta"}, atmplt.cdr2hv4.tdb},
	{{"A_TCR_CDR2HV4_" + std::to_string(btcr.cdr2.size()) + ".fasta", "B_TCR_CDR2HV4_" + std::to_string(btcr.cdr2.size()) + ".fasta"}, btmplt.cdr2hv4.tdb},
		};
	for ( auto &j : J ) {
		j.value = update_template_files(j.name, tcr_seq_db);
	}
	return;
}

std::set< std::string > create_ignore_list(std::string const query, std::string const blast_db, core::Real const &blast_cutoff) {
	run_blast(query, blast_db);
	std::string out_file = query+".blastp.out";
	std::ifstream bpofile (out_file);
	std::string line;
	std::string query_id, subject_id;
	core::Real percent_identity, alignment_length, mismatches, gap_opens, qstart, qend, sstart, send, evalue, bit_score;
	std::set< std::string > blastp_ignore_list;
	blastp_ignore_list.clear();
	while ( getline(bpofile, line) ) {
		if ( line[0] == '#' ) continue;
		std::istringstream ss(line);
		ss >> query_id >> subject_id >> percent_identity >> alignment_length >> mismatches >> gap_opens >> qstart >> qend >> sstart >> send >> evalue >> bit_score;
		if ( percent_identity > blast_cutoff ) {
			blastp_ignore_list.insert(subject_id.substr(0,4));
		}
	}
	return blastp_ignore_list;
}

void initialize_template_ignore_list(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, std::list< std::set<std::string> > &ignore_lists, core::Real const &blast_cutoff, std::string const &tcr_seq_db) {
	ignore_lists.clear();
	if ( blast_cutoff < 100 ) {
		std::string atcrdomain_file = "atcr.fasta";
		std::ofstream atcrdomain(atcrdomain_file);
		atcrdomain << "> " << atcrdomain_file << '\n' << atcr.truncdomain << '\n';
		atcrdomain.close();
		std::string btcrdomain_file = "btcr.fasta";
		std::ofstream btcrdomain(btcrdomain_file);
		btcrdomain << "> " << btcrdomain_file << '\n' << btcr.truncdomain << '\n';
		btcrdomain.close();
		std::string a_blast_db = tcr_seq_db + "A_TCR_VD";
		std::set< std::string > a_ignore_list;
		a_ignore_list = create_ignore_list(atcrdomain_file, a_blast_db, blast_cutoff);
		std::string b_blast_db = tcr_seq_db + "B_TCR_VD";
		std::set< std::string > b_ignore_list;
		b_ignore_list = create_ignore_list(btcrdomain_file, b_blast_db, blast_cutoff);
		//combine both sets, keep unique elements
		std::set< std::string > template_ignore_list;
		std::merge(a_ignore_list.begin(), a_ignore_list.end(),
			b_ignore_list.begin(), b_ignore_list.end(),
			std::inserter(template_ignore_list, template_ignore_list.begin()));
		ignore_lists.push_back(template_ignore_list);
	}
	//read user ignore_list file
	using namespace basic::options;
	std::string user_ignore_list_filename = option[OptionKeys::tcrmodel::ignore_list ];
	if ( !user_ignore_list_filename.empty() ) {
		std::string igpdbid;
		std::set< std::string > user_ignore_list;
		user_ignore_list.clear();
		std::ifstream igfile (user_ignore_list_filename);
		std::string igline;
		while ( getline(igfile, igline) ) {
			std::istringstream igss(igline);
			igss >> igpdbid;
			user_ignore_list.insert(utility::lower(igpdbid));
		}
		ignore_lists.push_back(user_ignore_list);
	}
	return;
}

void setup_template_from_pdbfile(TCRmodel::tmpltinfo &currtmplt) {
	if ( !utility::file::file_exists(currtmplt.tpdb) ) {
		utility_exit_with_message(currtmplt.tpdb + " : File not exists!");
	}
	core::pose::Pose tmplt_pose;
	core::import_pose::pose_from_file( tmplt_pose, currtmplt.tpdb, core::import_pose::PDB_file);
	currtmplt.tid = currtmplt.tpdb;
	currtmplt.tpiece = core::pose::PoseOP( new core::pose::Pose(tmplt_pose) );
	return;
}

void setup_cdr_template_from_pdbid(TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path) {
	std::string template_pdb_path = tpdb_path + currtmplt.tid.substr(0,4) + "_aho.pdb";
	graft_cdr(currtmplt, template_pdb_path, cdrstart, cdrend);
	return;
}

void setup_cdr_template_from_db(std::string const &seq, TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	std::pair<core::Real,std::string> tmplt_info = find_template_from_multiple_input_db(seq, currtmplt.tdb, ignore_lists, cutoff, tcr_ss);
	currtmplt.tscore = tmplt_info.first;
	currtmplt.tid = tmplt_info.second;
	if ( currtmplt.tscore > MINSCORE ) {
		std::string template_pdb_path = tpdb_path + tmplt_info.second.substr(0,4) + "_aho.pdb";
		graft_cdr(currtmplt, template_pdb_path, cdrstart, cdrend);
	} else {
		TR <<"Sequence: "<<seq<<std::endl;
		utility_exit_with_message("No Template found for TCR sequence");
	}
	return;
}

void setup_cdr_template(std::string const &seq, TCRmodel::tmpltinfo &currtmplt, core::Size const cdrstart, core::Size const cdrend, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	if ( currtmplt.tpdb != "" ) {
		setup_template_from_pdbfile(currtmplt);
	} else if ( currtmplt.tid != "" ) {
		setup_cdr_template_from_pdbid(currtmplt, cdrstart, cdrend, tpdb_path);
	} else {
		setup_cdr_template_from_db(seq, currtmplt, cdrstart, cdrend, tpdb_path, ignore_lists, cutoff, tcr_ss);
	}
	return;
}

void setup_fw_template_from_pdbid(TCRmodel::tmpltinfo &currtmplt, TCRseqInfo::tcrposi const &aho_pos, std::string const &tpdb_path) {
	std::string template_pdb_path = tpdb_path + currtmplt.tid.substr(0,4) + "_aho.pdb";
	graft_framework(currtmplt, template_pdb_path, aho_pos);
	return;
}

void setup_fw_template_from_db(std::string const &seq, TCRseqInfo::tcrposi const &aho_pos, TCRmodel::tmpltinfo &currtmplt, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	std::pair<core::Real, std::string> tmplt_info = find_template_from_multiple_input_db(seq, currtmplt.tdb, ignore_lists, cutoff, tcr_ss);
	currtmplt.tscore = tmplt_info.first;
	currtmplt.tid = tmplt_info.second;
	if ( currtmplt.tscore > MINSCORE ) {
		std::string template_pdb_path = tpdb_path + tmplt_info.second.substr(0,4) + "_aho.pdb";
		graft_framework(currtmplt, template_pdb_path, aho_pos);
	} else {
		TR <<"Sequence: "<<seq<<std::endl;
		utility_exit_with_message("No Template found for TCR sequence");
	}
	return;
}

void setup_fw_template(std::string const &seq, TCRseqInfo::tcrposi const &aho_pos, TCRmodel::tmpltinfo &currtmplt, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	if ( currtmplt.tpdb != "" ) {
		setup_template_from_pdbfile(currtmplt);
	} else if ( currtmplt.tid != "" ) {
		setup_fw_template_from_pdbid(currtmplt, aho_pos, tpdb_path);
	} else {
		setup_fw_template_from_db(seq, aho_pos, currtmplt, tpdb_path, ignore_lists, cutoff, tcr_ss);
	}
	return;
}

std::pair<std::string, std::string> get_orientation_template_from_db(std::string const &fra, std::string const &frb, TCRmodel::tmpltinfo const &oritmplt, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &tcr_ss) {
	std::string best_ori_a;
	std::string best_ori_b;
	std::string fwline;
	std::string curr_name_a;
	std::string curr_name_b;
	core::Real best_fw_score = MINSCORE;
	core::Real curr_fw_score = MINSCORE;
	core::Real best_ori_resolution = MAXRESOLUTION;
	char delim = ' ';
	//at present, ori template db is same for both atmplt and btmpt
	utility::vector1< std::string > multidb = oritmplt.tdb;
	for ( core::Size i=1; i<=multidb.size(); ++i ) {
		std::ifstream fwinp(multidb[i]);
		while ( std::getline(fwinp, fwline) ) {
			std::stringstream ss(fwline);
			std::string item;
			std::vector<std::string> tokens;
			while ( std::getline(ss, item, delim) ) { tokens.push_back(item); }
			if ( (fra.size() != tokens[2].size()) || (frb.size() != tokens[3].size()) ) continue;
			std::string fw_seq_querry = fra + frb;
			std::string fw_seq_content = tokens[2] + tokens[3];
			std::string first_four = tokens[0].substr(0, 4);
			if ( !ignore_template(first_four, ignore_lists) ) {
				if ( calculate_identity_score(fw_seq_querry,fw_seq_content) <= cutoff ) {
					core::Real fwa_score = score_alignment( fra, tokens[2], tcr_ss );
					core::Real fwb_score = score_alignment( frb, tokens[3], tcr_ss );
					curr_fw_score = fwa_score + fwb_score;
					curr_name_a = tokens[0];
					curr_name_b = tokens[1];
				}
			}
			core::Real curr_ori_resolution = MAXRESOLUTION;
			if ( !(tokens[4] == "None") ) {
				curr_ori_resolution = std::stof(tokens[4]);
			}
			if ( curr_fw_score > best_fw_score ) {
				best_fw_score = curr_fw_score;
				best_ori_a = curr_name_a;
				best_ori_b = curr_name_b;
				best_ori_resolution = curr_ori_resolution;
			} else if ( curr_fw_score == best_fw_score ) {
				if ( curr_ori_resolution < best_ori_resolution ) {
					best_fw_score = curr_fw_score;
					best_ori_a = curr_name_a;
					best_ori_b = curr_name_b;
					best_ori_resolution = curr_ori_resolution;
				}
			}
		}
	}
	if ( !(best_fw_score > MINSCORE) ) {
		utility_exit_with_message("No Template found for Orientation");
	}
	std::pair<std::string, std::string> best_ori_tmplt(best_ori_a, best_ori_b);
	return best_ori_tmplt;
}

void setup_orientation_template(TCRseqInfo::tcrsegs const &atcr, TCRseqInfo::tcrsegs const &btcr, TCRmodel::tmpltinfo &oria_tmplt, TCRmodel::tmpltinfo &orib_tmplt, TCRseqInfo::tcrposi const &aaho_posi, TCRseqInfo::tcrposi const &baho_posi, std::string const &tpdb_path, std::list< std::set< std::string > > const &ignore_lists, core::Real const &cutoff, core::sequence::ScoringSchemeOP const &ss) {
	std::string template_pdb_path;
	if ( ( oria_tmplt.tpdb != "" ) && ( orib_tmplt.tpdb != "" ) ) {
		setup_template_from_pdbfile(oria_tmplt);
		setup_template_from_pdbfile(orib_tmplt);
	} else  if ( ( oria_tmplt.tid != "" ) && ( oria_tmplt.tid != "" ) ) {
		template_pdb_path = tpdb_path + oria_tmplt.tid.substr(0,4) + "_aho.pdb";
		graft_framework(oria_tmplt, template_pdb_path, aaho_posi);
		template_pdb_path = tpdb_path + orib_tmplt.tid.substr(0,4) + "_aho.pdb";
		graft_framework(orib_tmplt, template_pdb_path, baho_posi);
	} else {
		std::pair<std::string, std::string> ori_tmplt = get_orientation_template_from_db(atcr.fr, btcr.fr, oria_tmplt, ignore_lists, cutoff, ss);
		oria_tmplt.tid = ori_tmplt.first;
		orib_tmplt.tid = ori_tmplt.second;
		template_pdb_path = tpdb_path + oria_tmplt.tid.substr(0,4) + "_aho.pdb";
		graft_framework(oria_tmplt, template_pdb_path, aaho_posi);
		template_pdb_path = tpdb_path + orib_tmplt.tid.substr(0,4) + "_aho.pdb";
		graft_framework(orib_tmplt, template_pdb_path, baho_posi);
	}
	return;
}

void dump_templates(TCRmodel::tcrtmplts const &atmplt, TCRmodel::tcrtmplts const &btmplt, bool const use_gma, bool const use_gmb){
	if ( !atmplt.ori.tpiece->empty() ) atmplt.ori.tpiece->dump_pdb("oria_tmplt.pdb");
	if ( !btmplt.ori.tpiece->empty() ) btmplt.ori.tpiece->dump_pdb("orib_tmplt.pdb");
	if ( use_gma ) {
		if ( !atmplt.gm.tpiece->empty() ) atmplt.gm.tpiece->dump_pdb("gma_tmplt.pdb");
	} else {
		if ( !atmplt.fr.tpiece->empty() ) atmplt.fr.tpiece->dump_pdb("fra_tmplt.pdb");
		if ( !atmplt.cdr1.tpiece->empty() ) atmplt.cdr1.tpiece->dump_pdb("cdr1a_tmplt.pdb");
		if ( !atmplt.cdr2hv4.tpiece->empty() ) atmplt.cdr2hv4.tpiece->dump_pdb("cdr2a_tmplt.pdb");
	}
	if ( !atmplt.cdr3.tpiece->empty() ) atmplt.cdr3.tpiece->dump_pdb("cdr3a_tmplt.pdb");
	if ( use_gmb ) {
		if ( !btmplt.gm.tpiece->empty() ) btmplt.gm.tpiece->dump_pdb("gmb_tmplt.pdb");
	} else {
		if ( !btmplt.fr.tpiece->empty() ) btmplt.fr.tpiece->dump_pdb("frb_tmplt.pdb");
		if ( !btmplt.cdr1.tpiece->empty() ) btmplt.cdr1.tpiece->dump_pdb("cdr1b_tmplt.pdb");
		if ( !btmplt.cdr2hv4.tpiece->empty() ) btmplt.cdr2hv4.tpiece->dump_pdb("cdr2b_tmplt.pdb");
	}
	if ( !btmplt.cdr3.tpiece->empty() ) btmplt.cdr3.tpiece->dump_pdb("cdr3b_tmplt.pdb");
	return;
}

} //protocols
} //tcr


