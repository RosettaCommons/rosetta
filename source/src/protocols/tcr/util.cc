// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/tcr/util.cc
/// @brief utility functions for protocols/tcr
/// @author Ragul (ragul@umd.edu)

#include <protocols/antibody/grafting/util.hh>

#ifdef __ANTIBODY_GRAFTING__

#include <protocols/tcr/util.hh>
#include <protocols/tcr/TCRseqInfo.hh>
#include <protocols/tcr/TCRseqInfo.fwd.hh>
#include <protocols/tcr/template_util.hh>
#include <core/pose/util.hh>
#include <core/pose/PDBInfo.hh>
#include <basic/Tracer.hh>
#include <basic/execute.hh>
#include <core/types.hh>
#include <regex>
#include <basic/options/option.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <protocols/relax/FastRelax.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/sequence/Sequence.hh>

static basic::Tracer TR( "protocols.tcr.util" );
/// @brief MINSCORE const value with large negative number
/// @details used for checking and ignoring templates
/// @details also used for initializing template score value
const core::Real MINSCORE = -99999.0;


namespace protocols {
namespace tcr {

void initialize_aho_numbers(TCRseqInfo::tcrposi &a_aho_pos, TCRseqInfo::tcrposi &b_aho_pos) {
	//adopted from AntibodyNumberingParser
	utility::vector1< utility::vector1< std::string > > tcr_num_scheme;
	utility::file::FileName tcr_cdr_definition_file(basic::database::full_name("sampling/antibodies/numbering_schemes/tcr_cdr_definitions.txt"));
	std::ifstream tcr_cdr_numbering_file(tcr_cdr_definition_file);
	std::string line;
	while ( getline(tcr_cdr_numbering_file, line) ) {
		if ( utility::startswith(line, "#") ) continue;
		if ( utility::startswith(line, "\n") ) continue;
		if ( line == "" ) continue;
		utility::trim(line, "\n");
		utility::vector1< std::string > lineSP = utility::string_split_multi_delim(line);
		if ( (lineSP.size()) == 0 ) continue;
		if ( lineSP[1] == "DEFINES" ) {
			for ( core::Size i = 2; i <= lineSP.size(); ++i ) {
				TR <<"Reading numbering scheme: "<< lineSP[i]<<std::endl;
			}
		} else {
			for ( core::Size i = 3; i <= lineSP.size(); ++i ) {
				utility::vector1< std::string > raw_landmark = utility::string_split(lineSP[i], ':');
				if ( raw_landmark.size() != 3 ) {
					utility_exit_with_message("Error while reading scheme file. landmark must define: chain:resnum:insertion code. If no insertion code, please use a period (.)");
				}
				tcr_num_scheme.push_back(raw_landmark);
			}
		}
	}
	a_aho_pos.cdr1.begin = std::stoi(tcr_num_scheme[1][2]);
	a_aho_pos.cdr1.end = std::stoi(tcr_num_scheme[2][2]);
	a_aho_pos.cdr2hv4.begin = std::stoi(tcr_num_scheme[3][2]);
	a_aho_pos.cdr2hv4.end = std::stoi(tcr_num_scheme[4][2]);
	a_aho_pos.cdr3.begin = std::stoul(tcr_num_scheme[5][2]);
	a_aho_pos.cdr3.end = std::stoi(tcr_num_scheme[6][2]);
	b_aho_pos.cdr1.begin = std::stoi(tcr_num_scheme[7][2]);
	b_aho_pos.cdr1.end = std::stoi(tcr_num_scheme[8][2]);
	b_aho_pos.cdr2hv4.begin = std::stoi(tcr_num_scheme[9][2]);
	b_aho_pos.cdr2hv4.end = std::stoi(tcr_num_scheme[10][2]);
	b_aho_pos.cdr3.begin = std::stoi(tcr_num_scheme[11][2]);
	b_aho_pos.cdr3.end = std::stoi(tcr_num_scheme[12][2]);
	//These are the numbers to cap the TCR sequence and restrict to variable domain
	//Users may enter whole sequence for variable and constant domain, but we model only the variable domain
	a_aho_pos.cap.begin = 20;
	a_aho_pos.cap.end = 10;
	b_aho_pos.cap.begin = 20;
	b_aho_pos.cap.end = 10;
	return;
}

core::Real calculate_identity_score( std::string const &query, std::string const &content) {
	if ( query.size() != content.size() ) return MINSCORE;
	core::Real iden_sc(0.);
	core::Size num_iden_res(0);
	for ( core::Size x = 0; x < query.size(); ++x ) {
		if ( query[x] == content[x] ) {
			num_iden_res += 1;
		}
	}
	iden_sc = ((core::Real)num_iden_res/(core::Real)query.size())*100;
	return iden_sc;
}

core::Real score_alignment( std::string const &query, std::string const &content, core::sequence::ScoringSchemeOP const &tcr_ss) {
	if ( query.size() != content.size() ) return MINSCORE;
	core::Real alignment_score(0.);
	core::sequence::SequenceOP query_seq( new core::sequence::Sequence( query, "query", 1 ) );
	core::sequence::SequenceOP content_seq( new core::sequence::Sequence( content, "content", 1 ) );
	for ( core::Size x = 1; x <= query.size(); ++x ) {
		alignment_score += tcr_ss->score(query_seq, content_seq, x, x);
	}
	return alignment_score;
}

void relax_tcr_model(core::pose::Pose &model_pose) {
	//Fast relax
	protocols::relax::FastRelaxOP frelax( new protocols::relax::FastRelax(core::scoring::get_score_function()) );
	//movemap
	core::kinematics::MoveMapOP mm( new core::kinematics::MoveMap );
	mm->set_chi(true);
	mm->set_bb( true );
	mm->set_jump( true );
	frelax->set_movemap( mm );
	//task factory
	core::pack::task::TaskFactoryOP task_factory( new core::pack::task::TaskFactory() );
	task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::IncludeCurrent() ) );
	task_factory->push_back( core::pack::task::operation::TaskOperationCOP( new core::pack::task::operation::RestrictToRepacking() ) );
	frelax->set_task_factory( task_factory );
	//fast relax options
	frelax->constrain_relax_to_start_coords( true );
	frelax->coord_constrain_sidechains( true );
	frelax->ramp_down_constraints( false );
	frelax->apply(model_pose);
}

void repack_tcr_model(core::pose::Pose &model_pose) {
	//repack
	int nres = model_pose.total_residue();
	core::pack::task::PackerTaskOP base_packer_task( core::pack::task::TaskFactory::create_packer_task( model_pose ));
	base_packer_task->initialize_from_command_line().or_include_current( true );
	for ( int ii = 1; ii <= nres; ++ii ) {
		base_packer_task->nonconst_residue_task(ii).restrict_to_repacking();
	}
	using namespace core::scoring;
	core::scoring::ScoreFunctionOP repack_scorefxn = core::scoring::get_score_function();
	core::pack::pack_rotamers( model_pose, *repack_scorefxn, base_packer_task );
}

void minimize_tcr_model(core::pose::Pose &model_pose) {
	using namespace core::id;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::scoring::constraints;
	core::scoring::ScoreFunctionOP minimize_scorefxn = core::scoring::get_score_function();
	//repack
	repack_tcr_model(model_pose);
	//Constrained minimization of model pose
	//AtomTreeMinimizer minimizer;
	int nres = model_pose.total_residue();
	core::scoring::constraints::ConstraintSetOP cst_set( new core::scoring::constraints::ConstraintSet() );
	core::scoring::func::HarmonicFuncOP spring(
		new core::scoring::func::HarmonicFunc( 0 , 1 ));
	core::conformation::Conformation const &conformation( model_pose.conformation() );
	core::Size const my_anchor = 1;
	core::kinematics::FoldTree original_ftree( model_pose.fold_tree() );
	core::kinematics::FoldTree rerooted_fold_tree = original_ftree;
	rerooted_fold_tree.reorder( my_anchor );
	model_pose.fold_tree( rerooted_fold_tree);
	for ( int i=1; i <= nres; i++ ) {
		core::conformation::Residue const &nat_i_rsd( model_pose.residue(i) );
		for ( core::Size ii = 1; ii<= nat_i_rsd.nheavyatoms(); ++ii ) {
			AtomID CAi ( ii, i );
			ConstraintOP cst(new CoordinateConstraint(
				CAi, AtomID(1,my_anchor), conformation.xyz( CAi ), spring ));
			cst_set->add_constraint(cst);
		}
	}
	core::Real cst_weight = 1;
	model_pose.constraint_set( cst_set );
	minimize_scorefxn->set_weight( coordinate_constraint, cst_weight );
	core::optimization::AtomTreeMinimizer minimizer;
	//core::optimization::MinimizerOptions min_options( "dfpmin", 0.00001, true, false );
	core::optimization::MinimizerOptions min_options( "lbfgs_armijo_nonmonotone", 0.001, true, false );
	//min_options.nblist_auto_update(true);
	core::kinematics::MoveMap mm_all;
	mm_all.set_chi( true );
	mm_all.set_bb( true );
	mm_all.set_jump( true );
	minimizer.run( model_pose, mm_all, *minimize_scorefxn, min_options );
	model_pose.remove_constraints();
	model_pose.fold_tree( original_ftree );
}

void assign_achain_CDRs_using_REGEX( std::string const &seq, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi ) {
	if( ! protocols::antibody::grafting::antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex.");
	}
	utility::vector1< std::regex > tcra_regex_list;
	std::regex tcra_regex;
	tcra_regex_list.clear();
	tcra_regex = "^[A-Z]{0,10}(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*";
	tcra_regex_list.push_back(tcra_regex);
	tcra_regex = "^[A-Z]{0,10}(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{10}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*";
	tcra_regex_list.push_back(tcra_regex);
	tcra_regex = "^[A-Z]*(([A-Z]{19}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|V])([A-Z]{1,36})([L|I|F][A-Z]{12}Y[A-Z][C|W])([A-Z]{1,32})([F|W]G[A-Z]G[A-Z]{6}))[A-Z]*";
	tcra_regex_list.push_back(tcra_regex);

	std::smatch tcr_res;
	std::regex tcr_regex;
	bool found_regex_match = false;
	for ( core::Size i = 1; i <= tcra_regex_list.size(); ++i ) {
		tcr_regex = tcra_regex_list[i];
		if ( std::regex_match(seq, tcr_res, tcr_regex) ) {
			parsedseqs.truncdomain = tcr_res[1].str();
			parsedseqs.cdr1 = tcr_res[3].str();
			parsedseqs.cdr2 = tcr_res[5].str();
			parsedseqs.cdr3 = tcr_res[7].str();
			parsedseqs.fr = tcr_res[2].str() + tcr_res[4].str() + tcr_res[6].str() + tcr_res[8].str();
			parsedseqs.gm = tcr_res[2].str() + tcr_res[3].str() + tcr_res[4].str() + tcr_res[5].str() + tcr_res[6].str() + tcr_res[8].str();
			found_regex_match = true;
			break;
		}
	}
	if ( found_regex_match ) {
		pose_posi.cdr1.begin = tcr_res.position(3)-tcr_res.position(1)+1;
		pose_posi.cdr1.end = tcr_res.position(4)-tcr_res.position(1);
		pose_posi.cdr2hv4.begin = tcr_res.position(5)-tcr_res.position(1)+1;
		pose_posi.cdr2hv4.end = tcr_res.position(6)-tcr_res.position(1);
		pose_posi.cdr3.begin = tcr_res.position(7)-tcr_res.position(1)+1;
		pose_posi.cdr3.end = tcr_res.position(8)-tcr_res.position(1);
	}
	if ( !found_regex_match ) {
		utility_exit_with_message("No RegEx match for TCR chain sequence : " + seq);
	}
	return;
}

void assign_bchain_CDRs_using_REGEX( std::string const &seq, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi ) {
	if( ! protocols::antibody::grafting::antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex.");
	}
	utility::vector1< std::regex > tcrb_regex_list;
	std::regex tcrb_regex;
	tcrb_regex_list.clear();
	tcrb_regex = "^[A-Z]{0,5}(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C)([A-Z]C[A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*";
	tcrb_regex_list.push_back(tcrb_regex);
	tcrb_regex = "^[A-Z]{0,5}(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*";
	tcrb_regex_list.push_back(tcrb_regex);
	tcrb_regex = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{12}C)([A-Z]C[A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*";
	tcrb_regex_list.push_back(tcrb_regex);
	tcrb_regex = "^[A-Z]*(([A-Z]{2}Q[A-Z]{16}C)([A-Z]{1,19})(W[A-Z]{11}[L|I|M])([A-Z]{1,36})([L|M][A-Z]{14}C)([A-Z]{1,32})(FG[A-Z]G[A-Z]{2}L[A-Z]{3}))[A-Z]*";
	tcrb_regex_list.push_back(tcrb_regex);

	std::smatch tcr_res;
	std::regex tcr_regex;
	bool found_regex_match = false;
	for ( core::Size i = 1; i <= tcrb_regex_list.size(); ++i ) {
		tcr_regex = tcrb_regex_list[i];
		if ( std::regex_match(seq, tcr_res, tcr_regex) ) {
			parsedseqs.truncdomain = tcr_res[1].str();
			parsedseqs.cdr1 = tcr_res[3].str();
			parsedseqs.cdr2 = tcr_res[5].str();
			parsedseqs.cdr3 = tcr_res[7].str();
			parsedseqs.fr = tcr_res[2].str() + tcr_res[4].str() + tcr_res[6].str() + tcr_res[8].str();
			parsedseqs.gm = tcr_res[2].str() + tcr_res[3].str() + tcr_res[4].str() + tcr_res[5].str() + tcr_res[6].str() + tcr_res[8].str();
			found_regex_match = true;
			break;
		}
	}
	if ( found_regex_match ) {
		pose_posi.cdr1.begin = tcr_res.position(3)-tcr_res.position(1)+1;
		pose_posi.cdr1.end = tcr_res.position(4)-tcr_res.position(1);
		pose_posi.cdr2hv4.begin = tcr_res.position(5)-tcr_res.position(1)+1;
		pose_posi.cdr2hv4.end = tcr_res.position(6)-tcr_res.position(1);
		pose_posi.cdr3.begin = tcr_res.position(7)-tcr_res.position(1)+1;
		pose_posi.cdr3.end = tcr_res.position(8)-tcr_res.position(1);
	}
	if ( !found_regex_match ) {
		utility_exit_with_message("No RegEx match for TCR chain sequence : " + seq);
	}
	return;
}

void adjust_position_for_chain(TCRseqInfo::tcrposi &pose_posi, core::Size const &chainlen ) {
	pose_posi.cdr1.begin += chainlen;
	pose_posi.cdr1.end += chainlen;
	pose_posi.cdr2hv4.begin += chainlen;
	pose_posi.cdr2hv4.end += chainlen;
	pose_posi.cdr3.begin += chainlen;
	pose_posi.cdr3.end += chainlen;
	return;
}

void assign_CDRs_using_numbers( std::string const &seq, antibody::grafting::CDR_Bounds const &cdr1seqpos, antibody::grafting::CDR_Bounds const &cdr2seqpos, antibody::grafting::CDR_Bounds const &cdr3seqpos, antibody::grafting::CDR_Bounds const &cap, TCRseqInfo::tcrsegs &parsedseqs, TCRseqInfo::tcrposi &pose_posi ) {

	//second parameter of substr() should be the length of the substring
	//The first character is denoted by a value of 0
	core::Size startcap = cap.begin;
	core::Size endcap = cap.end;
	core::Size seq_begin = cdr1seqpos.begin - startcap - 1;
	core::Size seq_end = ( cdr3seqpos.end + endcap ) - seq_begin;
	core::Size cdr1_start = cdr1seqpos.begin - seq_begin;
	core::Size cdr1_end = cdr1seqpos.end - seq_begin;
	core::Size cdr2_start = cdr2seqpos.begin - seq_begin;
	core::Size cdr2_end = cdr2seqpos.end - seq_begin;
	core::Size cdr3_start = cdr3seqpos.begin - seq_begin;
	core::Size cdr3_end = cdr3seqpos.end - seq_begin;
	parsedseqs.truncdomain = seq.substr (seq_begin, seq_end);
	parsedseqs.cdr1 = seq.substr ( cdr1seqpos.begin - 1, cdr1seqpos.end - cdr1seqpos.begin + 1 );
	parsedseqs.cdr2 = seq.substr ( cdr2seqpos.begin - 1, cdr2seqpos.end - cdr2seqpos.begin + 1 ) ;
	parsedseqs.cdr3 = seq.substr ( cdr3seqpos.begin - 1, cdr3seqpos.end - cdr3seqpos.begin + 1 );
	parsedseqs.fr = parsedseqs.truncdomain;
	parsedseqs.fr.erase(cdr3_start - 1, parsedseqs.cdr3.length());
	parsedseqs.fr.erase(cdr2_start - 1, parsedseqs.cdr2.length());
	parsedseqs.fr.erase(cdr1_start - 1, parsedseqs.cdr1.length());
	parsedseqs.gm = parsedseqs.truncdomain;
	parsedseqs.gm.erase(cdr3_start - 1, parsedseqs.cdr3.length());
	//core::Size seq_len = parsedseqs.truncdomain.length();
	pose_posi.cdr1.begin = cdr1_start;
	pose_posi.cdr1.end = cdr1_end;
	pose_posi.cdr2hv4.begin = cdr2_start;
	pose_posi.cdr2hv4.end = cdr2_end;
	pose_posi.cdr3.begin = cdr3_start;
	pose_posi.cdr3.end = cdr3_end;
	return;
}

void assign_CDRs_using_anarci( std::string const &seq, std::string const &anarci_path, TCRseqInfo::tcrposi &aho_pos, TCRseqInfo::tcrposi &pose_posi, TCRseqInfo::tcrsegs &parsedseqs ) {

	antibody::grafting::CDR_Bounds cdr1seqpos;
	antibody::grafting::CDR_Bounds cdr2seqpos;
	antibody::grafting::CDR_Bounds cdr3seqpos;
	std::string myline;
	std::smatch reres;
	std::smatch res_cdr1_start, res_cdr1_end, res_cdr2_start, res_cdr2_end, res_cdr3_start, res_cdr3_end;
	int domain_start = 0, domain_end = 0;
	std::string domain_tag;
	int cdr1_seq_start_pos = 0, cdr1_seq_end_pos = 0, cdr2_seq_start_pos = 0, cdr2_seq_end_pos = 0, cdr3_seq_start_pos = 0, cdr3_seq_end_pos = 0;
	if( ! protocols::antibody::grafting::antibody_grafting_usable() ) {
		utility_exit_with_message("ERROR: Your compiler does not have full support for C++11 regex.");
	}
	std::regex base_regex("#[/|].*[/|](.*)[/|].*[/|].*[/|]([[:d:]]+)[/|]([[:d:]]+)[/|]");
	std::string anarci_command (anarci_path+" -r tr -s a -i "+seq+" -o anarci.out");
	basic::execute("Running ANARCI... ", anarci_command);
	std::ifstream myfile ("anarci.out");
	int count = 0;
	if ( myfile.is_open() ) {
		while ( std::getline (myfile,myline) ) {
			if ( (myline[0] == '#') && (std::regex_match(myline, reres, base_regex)) ) {
				domain_tag = reres[1].str();
				domain_start =  std::stoi(reres[2]);
				domain_end =  std::stoi(reres[3]);
			}
			std::regex cdr1_start_regex(domain_tag+" "+std::to_string(aho_pos.cdr1.begin)+" .*");
			std::regex cdr1_end_regex(domain_tag+" "+std::to_string(aho_pos.cdr1.end)+" .*");
			std::regex cdr2_start_regex(domain_tag+" "+std::to_string(aho_pos.cdr2hv4.begin)+" .*");
			std::regex cdr2_end_regex(domain_tag+" "+std::to_string(aho_pos.cdr2hv4.end)+" .*");
			std::regex cdr3_start_regex(domain_tag+" "+std::to_string(aho_pos.cdr3.begin)+" .*");
			std::regex cdr3_end_regex(domain_tag+" "+std::to_string(aho_pos.cdr3.end)+" .*");
			if ( myline[0] == domain_tag[0] ) {
				if ( myline[10] == '-' ) continue;
				count +=1;
				if ( std::regex_match(myline, res_cdr1_start, cdr1_start_regex) ) cdr1_seq_start_pos = count;
				if ( std::regex_match(myline, res_cdr1_end, cdr1_end_regex) ) cdr1_seq_end_pos = count;
				if ( std::regex_match(myline, res_cdr2_start, cdr2_start_regex) ) cdr2_seq_start_pos = count;
				if ( std::regex_match(myline, res_cdr2_end, cdr2_end_regex) ) cdr2_seq_end_pos = count;
				if ( std::regex_match(myline, res_cdr3_start, cdr3_start_regex) ) cdr3_seq_start_pos = count;
				if ( std::regex_match(myline, res_cdr3_end, cdr3_end_regex) ) cdr3_seq_end_pos = count;
			}
		}
		myfile.close();
	} else {
		utility_exit_with_message("Unable to open ANARCI output file");
	}
	cdr1seqpos.begin = cdr1_seq_start_pos;
	cdr1seqpos.end = cdr1_seq_end_pos;
	cdr2seqpos.begin = cdr2_seq_start_pos;
	cdr2seqpos.end = cdr2_seq_end_pos;
	cdr3seqpos.begin = cdr3_seq_start_pos;
	cdr3seqpos.end = cdr3_seq_end_pos;
	//core::Size startcap = aho_pos.cap.begin;
	//core::Size endcap = aho_pos.cap.end;
	std::string domainseq = seq.substr(domain_start, domain_end - domain_start + 1);
	assign_CDRs_using_numbers( domainseq, cdr1seqpos, cdr2seqpos, cdr3seqpos, aho_pos.cap, parsedseqs, pose_posi );
	return;
}

void run_blast( std::string const &query, std::string const &blast_db) {
	std::string out_file = query+".blastp.out";
	std::string query_file = query;
	std::string blastp_exec = "blastp";
	std::string blast_command (blastp_exec+" -db "+blast_db+" -query "+query_file+" -out "+out_file+" -word_size 7 -outfmt 7 -max_target_seqs 1024 -evalue 0.00001 -matrix BLOSUM62");
	basic::execute("Running BLASTP... ", blast_command);
	return;
}

antibody::grafting::CDR_Bounds string_to_CDRbounds(std::string const &position ) {
	utility::vector1< std::string > splitposi = utility::string_split(position, ':');
	antibody::grafting::CDR_Bounds myseg;
	myseg.begin = std::stoi( splitposi[1] );
	myseg.end = std::stoi( splitposi[2] );
	return myseg;
}


bool check_seq_match( std::string const &query, std::string const &db, std::list< std::set< std::string > > const &ignore_lists) {
	bool exact_match = false;
	std::ifstream inp(db);
	if ( !inp.good() ) {
		utility_exit_with_message("Error opening template database : "+db);
	}
	std::string line, name, content;
	while ( std::getline(inp, line) ) {
		if ( line.empty() || line[0] == '>' ) {
			if ( !name.empty() ) {
				std::string first_four = name.substr(0, 4);
				if ( !ignore_template(first_four, ignore_lists) ) {
					if ( query == content ) {
						exact_match = true;
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
	return exact_match;
}

bool check_seq_match_from_multiple_input_db( std::string const &query, utility::vector1< std::string > const &multidb, std::list< std::set< std::string > > const &ignore_lists ) {
	bool seq_match = false;
	for ( core::Size i=1; i<=multidb.size(); ++i ) {
		seq_match = check_seq_match(query, multidb[i], ignore_lists);
		if ( seq_match ) {
			break;
		}
	}
	return seq_match;
}

utility::vector1< std::pair<core::Size, core::Size> > get_cdr_pdb_pos_from_input_pose( core::pose::Pose const &inpose, TCRseqInfo::tcrposi const &selepos) {
	char chainid = inpose.pdb_info()->chain( 1 );
	core::Size val_cdr1_begin = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr1.begin );
	core::Size val_cdr1_end = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr1.end );
	core::Size val_cdr2hv4_begin = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr2hv4.begin );
	core::Size val_cdr2hv4_end = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr2hv4.end );
	core::Size val_cdr3_begin = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr3.begin );
	core::Size val_cdr3_end = inpose.pdb_info()->pdb2pose( chainid, selepos.cdr3.end );
	core::Size chain_length = (inpose.chain_sequence(1)).length();
	utility::vector1< std::pair<core::Size, core::Size> > segments;
	segments.resize(4);
	segments[1] = std::make_pair(1,val_cdr1_begin);
	segments[2] = std::make_pair(val_cdr1_end,val_cdr2hv4_begin);
	segments[3] = std::make_pair(val_cdr2hv4_end,val_cdr3_begin);
	segments[4] = std::make_pair(val_cdr3_end,chain_length);
	return segments;
}

} //protocols
} //tcr

#endif // __ANTIBODY_GRAFTING__

