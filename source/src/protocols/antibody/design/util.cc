// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/util.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/database/CDRSetOptions.hh>
#include <protocols/antibody/database/CDRSetOptionsParser.hh>
#include <protocols/antibody/design/CDRSeqDesignOptions.hh>
#include <protocols/antibody/design/CDRGraftDesignOptions.hh>
#include <protocols/antibody/database/AntibodyDatabaseManager.hh>
#include <protocols/antibody/design/NativeAntibodySeq.hh>

#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/grafting/util.hh>
#include <protocols/loops/util.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/datacache/CacheableDataType.hh>

#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>
#include <basic/datacache/CacheableData.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>

#include <boost/algorithm/string.hpp>
#include <string>
#include <utility/exit.hh>
#include <iostream>
#include <basic/Tracer.hh>

#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR("protocols.antibody.design.util");

namespace protocols {
namespace antibody {
namespace design {
using namespace core::pack::task::operation;
using namespace core::scoring;
using namespace protocols::antibody;
using namespace protocols::antibody::clusters;
using namespace utility;
using namespace basic::options;
using namespace core::pose::datacache;

ScoreFunctionOP
get_ab_design_global_scorefxn() {
	ScoreFunctionOP global_scorefxn = get_score_function();

	//Add any atom pair constraints
	if ( option[OptionKeys::antibody::design::global_atom_pair_cst_scoring]() ) {
		global_scorefxn->set_weight(atom_pair_constraint, option[ OptionKeys::antibody::design::atom_pair_cst_weight]());
	}
	if ( option[ OptionKeys::antibody::design::global_dihedral_cst_scoring]() ) {
		global_scorefxn->set_weight(dihedral_constraint, option[ OptionKeys::antibody::design::dihedral_cst_weight]());
	}
	return global_scorefxn;
}



ScoreFunctionOP
get_ab_design_global_scorefxn(utility::tag::TagCOP tag, basic::datacache::DataMap & data) {

	ScoreFunctionOP global_scorefxn;
	if ( tag->hasOption("scorefxn") ) {
		global_scorefxn = ( protocols::rosetta_scripts::parse_score_function( tag, data ) );
	} else {
		global_scorefxn = get_score_function();
	}

	core::Real atom_pair_cst_weight = tag->getOption<core::Real>("atom_pair_cst_weight",
		option[ OptionKeys::antibody::design::atom_pair_cst_weight]());

	core::Real dihedral_cst_weight = tag->getOption<core::Real>("dihedral_cst_weight",
		option[ OptionKeys::antibody::design::dihedral_cst_weight]());

	//Add any atom pair constraints
	bool global_atom_pair_cst_scoring = tag->getOption<bool>("global_atom_pair_cst_scoring",
		option[OptionKeys::antibody::design::global_atom_pair_cst_scoring]());

	bool global_dihedral_cst_scoring = tag->getOption<bool>("global_dihedral_cst_scoring",
		option[ OptionKeys::antibody::design::global_dihedral_cst_scoring]());


	if ( global_atom_pair_cst_scoring ) {
		global_scorefxn->set_weight(atom_pair_constraint, atom_pair_cst_weight);
	}

	if ( global_dihedral_cst_scoring ) {
		global_scorefxn->set_weight(dihedral_constraint, dihedral_cst_weight);
	}

	return global_scorefxn;
}

void
attributes_for_get_ab_design_global_scorefxn(utility::tag::AttributeList& attlist) {
	using namespace utility::tag;

	rosetta_scripts::attributes_for_parse_score_function(attlist);

	if ( !attribute_w_name_in_attribute_list("atom_pair_cst_weight", attlist) ) {
		attlist + XMLSchemaAttribute(
			"atom_pair_cst_weight", xsct_real,
			"Weight for atom pair constraints to use in the global antibody design score function" );
	}

	if ( !attribute_w_name_in_attribute_list("dihedral_cst_weight", attlist) ) {
		attlist + XMLSchemaAttribute(
			"dihedral_cst_weight", xsct_real,
			"Weight for atom pair constraints to use in the global antibody design score function" );
	}

	if ( !attribute_w_name_in_attribute_list("global_atom_pair_cst_scoring", attlist) ) {
		attlist + XMLSchemaAttribute(
			"global_atom_pair_cst_scoring", xsct_rosetta_bool,
			"Score atom pair constraints to use in the global antibody design score function?");
	}

	if ( !attribute_w_name_in_attribute_list("global_dihedral_cst_scoring", attlist) ) {
		attlist + XMLSchemaAttribute(
			"global_dihedral_cst_scoring", xsct_rosetta_bool,
			"Score dihedral constraints to use in the global antibody design score function?");
	}
}

ScoreFunctionOP
get_ab_design_dock_high_scorefxn() {
	ScoreFunctionOP docking_scorefxn_high = ScoreFunctionFactory::create_score_function( "docking", "docking_min" );
	docking_scorefxn_high->set_weight(atom_pair_constraint, option[ OptionKeys::antibody::design::atom_pair_cst_weight]());

	//High res dock will not affect dihedral angles - so no use adding csts for them here.
	return docking_scorefxn_high;

}

ScoreFunctionOP
get_ab_design_dock_low_scorefxn(){
	ScoreFunctionOP docking_scorefxn_low = ScoreFunctionFactory::create_score_function("interchain_cen");
	docking_scorefxn_low->set_weight(atom_pair_constraint, option[ OptionKeys::antibody::design::atom_pair_cst_weight]());

	//Low res dock will not effect dihedral angles - so no use adding  csts for them here.
	return docking_scorefxn_low;
}

ScoreFunctionOP
get_ab_design_min_scorefxn(){
	ScoreFunctionOP min_scorefxn = get_score_function();
	min_scorefxn->set_weight(dihedral_constraint, option[ OptionKeys::antibody::design::dihedral_cst_weight]());
	if ( option[ OptionKeys::antibody::design::global_atom_pair_cst_scoring]() ) {
		min_scorefxn->set_weight(atom_pair_constraint, option[ OptionKeys::antibody::design::atom_pair_cst_weight]());
	}
	return min_scorefxn;
}

ScoreFunctionOP
get_ab_design_min_scorefxn(utility::tag::TagCOP tag, basic::datacache::DataMap & data){
	ScoreFunctionOP min_scorefxn;
	if ( tag->hasOption("min_scorefxn") ) {
		std::string min_scorefxn_name = tag->getOption<std::string>("min_scorefxn");
		min_scorefxn = data.get_ptr<core::scoring::ScoreFunction>("scorefxns", min_scorefxn_name);
	} else {
		min_scorefxn = get_score_function();
	}

	core::Real atom_pair_cst_weight = tag->getOption<core::Real>("atom_pair_cst_weight",
		option[ OptionKeys::antibody::design::atom_pair_cst_weight]());

	core::Real dihedral_cst_weight = tag->getOption<core::Real>("dihedral_cst_weight",
		option[ OptionKeys::antibody::design::dihedral_cst_weight]());

	//Add any atom pair constraints
	bool global_atom_pair_cst_scoring = tag->getOption<bool>("global_atom_pair_cst_scoring",
		option[OptionKeys::antibody::design::global_atom_pair_cst_scoring]());

	min_scorefxn->set_weight(dihedral_constraint, dihedral_cst_weight);
	if ( global_atom_pair_cst_scoring ) {
		min_scorefxn->set_weight(atom_pair_constraint, atom_pair_cst_weight);
	}
	return min_scorefxn;

}

void
attributes_for_get_ab_design_min_scorefxn(utility::tag::AttributeList& attlist) {
	using namespace utility::tag;
	attlist + XMLSchemaAttribute(
		"min_scorefxn", xs_string,
		"Name of score function to use for minimization");

	if ( !attribute_w_name_in_attribute_list("atom_pair_cst_weight", attlist) ) {
		attlist + XMLSchemaAttribute(
			"atom_pair_cst_weight", xsct_real,
			"Weight for atom pair constraints to use in the antibody design minimization score function" );
	}

	if ( !attribute_w_name_in_attribute_list("dihedral_cst_weight", attlist) ) {
		attlist + XMLSchemaAttribute(
			"dihedral_cst_weight", xsct_real,
			"Weight for dihedral constraints to use in the antibody design minimization score function" );
	}

	if ( !attribute_w_name_in_attribute_list("global_atom_pair_cst_scoring", attlist) ) {
		attlist + XMLSchemaAttribute(
			"global_atom_pair_cst_scoring", xsct_rosetta_bool,
			"Use atom pair constraints in the antibody design minimization score function?" );
	}
}

void
insert_cdr_into_antibody(AntibodyInfoCOP ab_info, CDRNameEnum const cdr, core::pose::Pose & pose, core::pose::Pose & cdr_piece, core::Size overhang) {
	core::Size cdr_start = ab_info->get_CDR_start(cdr, pose);
	core::Size cdr_end = ab_info->get_CDR_end(cdr, pose);

	protocols::grafting::delete_region(cdr_piece, 1, overhang);
	protocols::grafting::delete_region(cdr_piece, cdr_piece.size() - overhang + 1, cdr_piece.size());

	//core::Size insert_length = cdr_piece.size();
	protocols::grafting::delete_region(pose, cdr_start+1, cdr_end - 1);
	pose = protocols::grafting::insert_pose_into_pose(pose, cdr_piece, cdr_start, cdr_start+1);

	pose.pdb_info()->copy(*(cdr_piece.pdb_info()), 1, cdr_piece.size(), cdr_start+1);
	pose.pdb_info()->obsolete(false);

}


vector1< vector1< core::Size > >
get_all_graft_permutations(
	vector1< vector1< core::Size> > permutations,
	vector1<core::Size > totals,
	core::Size const n
) {
	//Current index is what is being worked on.

	vector1<vector1<core::Size> > all_permutations;


	if ( n > totals.size() ) return permutations;

	else if ( n == 1 ) {

		if ( totals[1] == 0 ) {
			vector1<core::Size> zero_vector(1, 0);
			all_permutations.push_back(zero_vector);
		} else {
			for ( core::Size i = 1; i <= totals[1]; ++i ) {
				vector1<core::Size> start_vector(1, i);
				all_permutations.push_back(start_vector);
			}
		}
	} else {
		//No CDR in CDR set.  Set index 0, move on to next CDR
		if ( totals[n] == 0 ) {
			for ( core::Size i = 1; i <= permutations.size(); ++i ) {
				utility::vector1<core::Size> c = permutations[i];
				utility::vector1<core::Size> c2(c);
				c2.push_back(0);
				all_permutations.push_back(c2);
			}
		} else {
			for ( core::Size ind = 1; ind <= totals[n]; ++ind ) {

				for ( core::Size i = 1; i <= permutations.size(); ++i ) {
					utility::vector1<core::Size> c = permutations[i];
					utility::vector1<core::Size> c2(c);
					c2.push_back(ind);
					all_permutations.push_back(c2);
				}
			}
		}
	}

	vector1< vector1< core::Size > > combos = get_all_graft_permutations(all_permutations, totals, n+1);
	return combos;

}



AntibodyDesignProtocolEnum
design_protocol_to_enum(std::string const & design_type){
	std::string type = design_type;
	boost::to_upper(type);

	if ( type == "GENERALIZED_MONTE_CARLO" || type == "GEN_MC" ) {
		return generalized_monte_carlo;
	} else if ( type == "DETERMINISTIC_GRAFT" ) {
		return deterministic_graft;
	} else if ( type == "EVEN_CLUSTER_MC" ) {
		return even_cluster_monte_carlo;
	} else {
		utility_exit_with_message("DesignProtocol unrecognized.  Please check AntibodyDesign settings.");
	}
}

std::string
design_protocol_to_string(AntibodyDesignProtocolEnum const design_type){
	if ( design_type == generalized_monte_carlo ) {
		return "generalized_monte_carlo";
	} else if ( design_type==deterministic_graft ) {
		return "deterministic_graft";
	} else if ( design_type == even_cluster_monte_carlo ) {
		return "even_cluster_monte_carlo";
	} else {
		utility_exit_with_message("DesignProtocolunrecognized.  Please check AntibodyDesign settings.");
	}
}

SeqDesignStrategyEnum
seq_design_strategy_to_enum(std::string strategy){
	std::string option = strategy;
	boost::to_upper(option);

	if ( option == "CONSERVATIVE" ||  option=="CONSERVED" || option == "SEQ_DESIGN_CONSERVATIVE" ) {
		return seq_design_conservative;
	} else if ( option == "PROFILE" || option == "PROFILES" || option == "SEQ_DESIGN_PROFILES" ) {
		return seq_design_profiles;
	} else if ( option == "PROFILE_SETS" || option == "PROFILESETS" || option == "SEQ_DESIGN_PROFILE_SETS" ) {
		return seq_design_profile_sets;
	} else if ( option == "PROFILES_COMBINED" || option == "PROFILES_AND_SETS" || option == "PROFILES_AND_PROFILE_SETS" ) {
		return seq_design_profile_sets_combined;
	} else if ( option == "BASIC_DESIGN" || option == "BASIC"|| option == "SEQ_DESIGN_BASIC" || option == "NORMAL" ) {
		return seq_design_basic;
	} else if ( option == "SEQ_DESIGN_NONE" || option == "NO_DESIGN" || option == "OFF" ) {
		return seq_design_none;
	} else {
		utility_exit_with_message("Primary design strategy not understood: "+option);
	}
}

std::string
seq_design_strategy_to_string(SeqDesignStrategyEnum strategy){
	if ( strategy == seq_design_conservative ) {
		return "seq_design_conservative";
	} else if ( strategy == seq_design_profiles ) {
		return "seq_design_profiles";
	} else if ( strategy == seq_design_basic ) {
		return "seq_design_basic";
	} else if ( strategy == seq_design_none ) {
		return "seq_design_none";
	} else {
		utility_exit_with_message("Design strategy not understood");
	}
}

std::string
get_dock_chains_from_ab_dock_chains(AntibodyInfoCOP ab_info, std::string ab_dock_chains){
	vector1<char> chains;
	for ( core::Size i = 0; i <= ab_dock_chains.length()-1; ++i ) {
		char ab_chain = ab_dock_chains[i];
		if ( ab_chain == 'A' ) {
			vector1<char>antigen = ab_info->get_antigen_chains();
			if ( antigen.size() == 0 ) {
				TR <<" Antigen not present to dock. skipping addition of Antigen. - setting L_H dock instead" << std::endl;
				return "L_H";

			}

			for ( core::Size x = 1; x <= antigen.size(); ++x ) {
				chains.push_back(antigen[x]);
			}
		} else if ( ab_chain == 'L' || ab_chain == 'H' || ab_chain == '_' ) {
			chains.push_back(ab_chain);
		} else {
			utility_exit_with_message("ab_dock_chains must be L H or A: " + ab_dock_chains);
		}
	}
	std::string dock_chains(chains.begin(), chains.end());
	return dock_chains;
}

PDBNumbering
get_pdb_numbering_from_single_string( std::string const & pdb_residue){
	PDBNumbering number;

	vector1<std::string> res_str = utility::string_split(pdb_residue, ':');
	if ( res_str.size()== 1 ) {
		assert(res_str[1].length() >= 2);
		number.icode = ' ';
		number.chain = res_str[1][res_str[1].length() - 1];
		std::stringstream(res_str[1].substr(0, res_str[1].length() - 1)) >> number.resnum;
	} else if ( res_str.size() == 2 ) {
		assert(res_str[1].length() >= 2);
		assert(res_str[2].length() == 1);

		if ( res_str[2][0] == '~' ) {
			TR << "Found blank Icode as we should "<< std::endl;
			number.icode = ' ';
		} else {
			number.icode = res_str[2][0];
		}
		number.chain = res_str[1][res_str[1].length() - 1];
		std::stringstream(res_str[1].substr(0, res_str[1].length() - 1)) >> number.resnum;
	} else if ( res_str.size() == 3 ) {
		number.chain = res_str[1][0];
		number.resnum = utility::from_string(res_str[2], core::Size(0));
		if ( res_str[3][0] == '~' ) {
			number.icode = ' ';
		} else {
			number.icode = res_str[3][0];
		}
	} else {
		utility_exit_with_message("Cannot convert string to pdb string: "+pdb_residue);
	}




	return number;
}

core::Size
get_resnum_from_single_string(core::pose::Pose const & pose, std::string const & pdb_residue){
	PDBNumbering pdb_numbering = get_pdb_numbering_from_single_string(pdb_residue);
	return pose.pdb_info()->pdb2pose(pdb_numbering.chain, pdb_numbering.resnum, pdb_numbering.icode);
}

core::Size
get_resnum_from_single_string_w_landmark(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	std::string const & pdb_residue,
	AntibodyNumberingSchemeEnum const & scheme)
{
	PDBNumbering pdb_numbering = get_pdb_numbering_from_single_string(pdb_residue);
	return ab_info->get_landmark_resnum(pose, scheme, pdb_numbering.chain, pdb_numbering.resnum, pdb_numbering.icode, false);
}


vector1<PDBNumbering>
get_pdb_numbering_from_strings(vector1<std::string> const & pdb_residues) {
	vector1<PDBNumbering> numbering;


	for ( core::Size i = 1; i <= pdb_residues.size(); ++i ) {

		PDBNumbering number = get_pdb_numbering_from_single_string(pdb_residues[ i ]);
		numbering.push_back(number);
	}
	return numbering;
}

vector1<bool>
get_resnums_from_strings_with_ranges(core::pose::Pose const & pose, vector1<std::string> const & pdb_residues) {
	vector1<bool> numbers(pose.size(), false);


	for ( core::Size i = 1; i <= pdb_residues.size(); ++i ) {
		std::size_t found = pdb_residues[ i ].find("-");
		TR << pdb_residues[ i ] << std::endl;
		if ( found != std::string::npos ) {

			vector1<std::string> start_stop = utility::string_split(pdb_residues[ i ], '-');
			if ( start_stop.size() != 2 ) {
				utility_exit_with_message("Cannot parse numbering. Too many '-' "+ pdb_residues[ i ]);
			}

			PDBNumbering start_pdb_res = get_pdb_numbering_from_single_string(start_stop[ 1 ]);
			PDBNumbering end_pdb_res = get_pdb_numbering_from_single_string(start_stop[ 2 ]);

			core::Size start_resnum = pose.pdb_info()->pdb2pose(start_pdb_res.chain, start_pdb_res.resnum, start_pdb_res.icode);
			core::Size end_resnum = pose.pdb_info()->pdb2pose(end_pdb_res.chain, end_pdb_res.resnum, end_pdb_res.icode);
			for ( core::Size xi = start_resnum; xi <= end_resnum; ++xi ) {
				numbers[ xi ] = true;
				TR << "Adding: " << xi << std::endl;
			}
		} else {
			PDBNumbering numbering = get_pdb_numbering_from_single_string(pdb_residues[ i ]);
			core::Size resnum = pose.pdb_info()->pdb2pose(numbering.chain, numbering.resnum, numbering.icode);
			numbers[ resnum ] = true;
			TR << "Adding: " << resnum << std::endl;
		}
	}
	return numbers;
}

vector1<bool>
get_resnum_from_pdb_numbering(core::pose::Pose const & pose, vector1<PDBNumbering> const & pdb_residues){

	vector1<bool> residues(pose.size(), false);
	for ( core::Size i = 1; i <= pdb_residues.size(); ++i ) {
		PDBNumbering numbering = pdb_residues[i];
		core::Size resnum = pose.pdb_info()->pdb2pose(numbering.chain, numbering.resnum, numbering.icode);
		residues[resnum] = true;
	}
	return residues;
}

vector1<bool>
get_resnum_from_strings(core::pose::Pose const & pose, utility::vector1<std::string> const & pdb_residues) {
	vector1<bool> resnums = get_resnum_from_pdb_numbering(pose, get_pdb_numbering_from_strings(pdb_residues));
	return resnums;
}

void
add_loops_from_bool_vector(loops::Loops & loops, utility::vector1< bool > residues, bool add_cutpoints) {

	core::Size i = 1;
	while ( i <= residues.size() ) {
		//TR << i << std::endl;
		if ( residues[ i ] ) {
			//Search
			core::Size start = i; core::Size end = i;


			for ( core::Size x = start; x <= residues.size(); ++x ) {
				if ( residues[ x ] ) { end = x; }
				else { break; }
			}
			//TR << "Start: " << start<<"End: " << end << std::endl;

			//Add the loop
			if ( add_cutpoints ) {
				loops.add_loop(loops::Loop(start, end, start+ end-start/2));
			} else {
				loops.add_loop(loops::Loop(start, end));
			}

			//Increment I to end of segment:
			i = end;
		}
		i += 1;
	}
}

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, utility::vector1<bool> const & residues){
	assert(residues.size() == pose.size());

	std::pair<bool, core::Size> cb = std::make_pair(false, 0);
	for ( core::Size i = 1; i <= residues.size(); ++i ) {
		if ( ! residues[i] ) continue;
		cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, i, true, true, 1.5, 15, 15);
		if ( cb.first ) {
			return cb;
		}
	}

	return cb;
}

std::pair<bool, core::Size>
check_cb(core::pose::Pose const & pose, protocols::loops::Loops const & loops){
	std::pair<bool, core::Size> cb = std::make_pair(false, 0);
	for ( protocols::loops::Loops::const_iterator it = loops.begin(); it != loops.end(); ++it ) {
		cb = protocols::loops::has_severe_pep_bond_geom_issues(pose, it->cut(), true, true, 1.5, 15, 15);
		if ( cb.first == true ) {
			return cb;
		}
	}
	return cb;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_region(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	AntibodyRegionEnum region,
	bool cdr4_as_framework /* true */)
{
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ab_info->get_region_of_residue(pose, i, cdr4_as_framework) == region ) {
			restrict->include_residue(i);
		}
	}
	return restrict;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_antigen(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	return disable_design_region(ab_info, pose, antigen_region);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_framework(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	bool cdr4_as_framework /* true */)
{
	return disable_design_region(ab_info, pose, framework_region, cdr4_as_framework);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdrs(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	return disable_design_region(ab_info, pose, cdr_region);
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_design_cdr(
	AntibodyInfoCOP ab_info,
	CDRNameEnum cdr,
	const core::pose::Pose & pose) {

	if ( ! ab_info->has_CDR( cdr ) ) {
		utility_exit_with_message("disable_design_cdr antibody::design utility function failed \n"+ab_info->get_CDR_name( cdr ) +" is not present in antibody");
	}
	//One restrict op per CDR.  That way we can pop them off the TF  individually if we need to.
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	core::Size start = ab_info->get_CDR_start(cdr, pose);
	core::Size end = ab_info->get_CDR_end(cdr, pose);
	for ( core::Size i = start; i <= end; ++i ) {
		restrict->include_residue(i);
	}
	return restrict;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_conserved_framework_positions(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose)
{
	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	utility::vector1< core::Size > conserved_positions;

	if ( core::pose::has_chain( "L", pose ) ) {
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 23, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 43, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 106, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'L', 139, ' ', false ) );
	}

	if ( core::pose::has_chain( "H", pose ) ) {
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 23, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 43, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 106, ' ', false ) );
		conserved_positions.push_back( ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 139, ' ', false ) );
	}

	for ( core::Size i = 1; i <= conserved_positions.size(); ++i ) {
		if ( conserved_positions[ i ] != 0 ) {
			restrict->include_residue( conserved_positions[ i ] );
		}
	}
	return restrict;
}

core::pack::task::operation::RestrictResidueToRepackingOP
disable_h3_stem_positions(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose & pose,
	core::Size nter_stem,
	core::Size cter_stem)
{

	RestrictResidueToRepackingOP restrict( new RestrictResidueToRepacking() );
	utility::vector1< core::Size > stem_positions;

	core::Size start_cdr = ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 107, ' ', false);
	core::Size stop_cdr = ab_info->get_landmark_resnum( pose, AHO_Scheme, 'H', 138, ' ', false);

	for ( core::Size i = start_cdr; i <= start_cdr+nter_stem-1; ++i ) {
		stem_positions.push_back(i);
	}

	for ( core::Size i = stop_cdr-cter_stem+1; i <= stop_cdr; ++i ) {
		stem_positions.push_back(i);
	}

	for ( core::Size i = 1; i <= stem_positions.size(); ++i ) {
		if ( stem_positions[ i ] != 0 ) {
			restrict->include_residue( stem_positions[ i ] );
		}
	}
	return restrict;

}

AntibodyCDRSetOptions
get_cdr_set_options(){

	CDRSetOptionsParser parser = CDRSetOptionsParser();
	AntibodyCDRSetOptions options_settings;

	std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_cdr_instructions]();
	TR << "Reading CDRSetOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_options(cdr, filename));
	}
	return options_settings;
}

AntibodyCDRSetOptions
get_cdr_set_options(std::string filename){

	if ( filename.empty() || filename == "NA" ) {
		return get_cdr_set_options();
	}
	
	TR <<" Should be reading from: " << filename << std::endl;
	CDRSetOptionsParser parser = CDRSetOptionsParser();
	AntibodyCDRSetOptions options_settings;

	TR << "Reading CDRSetOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));

	}
	return options_settings;
}

AntibodyCDRGraftDesignOptions
get_graft_design_options(){

	CDRGraftDesignOptionsParser parser = CDRGraftDesignOptionsParser();
	AntibodyCDRGraftDesignOptions options_settings;

	std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_cdr_instructions]();
	TR << "Reading CDRGraftDesignOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_options(cdr, filename));

	}
	return options_settings;
}

AntibodyCDRGraftDesignOptions
get_graft_design_options(std::string filename){

	if ( filename.empty() || filename == "NA") {
		return get_graft_design_options();
	}
	CDRGraftDesignOptionsParser parser = CDRGraftDesignOptionsParser();
	AntibodyCDRGraftDesignOptions options_settings;

	TR << "Reading CDRGraftDesignOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <= core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));
	}
	return options_settings;
}

AntibodyCDRSeqDesignOptions
get_seq_design_options(){

	CDRSeqDesignOptionsParser parser = CDRSeqDesignOptionsParser();
	AntibodyCDRSeqDesignOptions options_settings;

	std::string filename = basic::options::option [basic::options::OptionKeys::antibody::design::base_cdr_instructions]();
	TR << "Reading CDRSeqDesignOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <=core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_options(cdr, filename));
	}
	return options_settings;
}

AntibodyCDRSeqDesignOptions
get_seq_design_options(std::string filename){

	if ( filename.empty() || filename == "NA") {
		return get_seq_design_options();
	}
	CDRSeqDesignOptionsParser parser = CDRSeqDesignOptionsParser();
	AntibodyCDRSeqDesignOptions options_settings;

	TR << "Reading CDRSeqDesignOptions from: " << filename << std::endl;
	for ( core::Size i = 1; i <=core::Size(CDRNameEnum_proto_total); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		options_settings.push_back(parser.parse_default_and_user_options(cdr, filename));
	}

	return options_settings;
}

std::map< core::Size, std::map< core::chemical::AA, core::Real > >
get_cluster_profile_probability_data(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose& pose,
	utility::vector1< bool > const & cdrs,
	utility::vector1< bool > & no_data_cdrs,
	const core::Size prob_cutoff,
	const bool use_outliers,
	const bool force_north_db,
	const bool ignore_light_chain){

	AntibodyDatabaseManager manager = AntibodyDatabaseManager( ab_info, force_north_db );
	manager.set_outlier_use(use_outliers);
	manager.ignore_light_chain( ignore_light_chain );

	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set;
	no_data_cdrs = manager.load_cdr_design_data_for_cdrs( cdrs, pose, prob_set, prob_cutoff );

	return prob_set;
}

std::map< core::Size, std::map< core::chemical::AA, core::Real > >
get_cluster_profile_probability_data(
	AntibodyInfoCOP ab_info,
	const core::pose::Pose& pose,
	AntibodyCDRSeqDesignOptions const & seq_design_options,
	utility::vector1< bool > & no_data_cdrs,
	const core::Size prob_cutoff,
	const bool use_outliers,
	const bool force_north_db,
	const bool ignore_light_chain){


	AntibodyDatabaseManager manager = AntibodyDatabaseManager( ab_info, force_north_db );
	manager.set_outlier_use(use_outliers);
	manager.ignore_light_chain( ignore_light_chain);
	
	std::map< core::Size, std::map< core::chemical::AA, core::Real > > prob_set;
	no_data_cdrs = manager.load_cdr_design_data( seq_design_options, pose, prob_set, prob_cutoff );

	return prob_set;
}

std::map<core::Size, core::chemical::AA>
transform_sequence_to_mutation_set(
	AntibodyInfoCOP ab_info,
	core::pose::Pose const & pose,
	CDRNameEnum const cdr,
	std::string const & sequence)
{

	if ( sequence.length()!= ab_info->get_CDR_length(cdr, pose, North) ) {
		utility_exit_with_message("Cannot transform sequence to mutation set.  Length mismatch " +
			utility::to_string(sequence.length()) +" vs "+ utility::to_string(ab_info->get_CDR_length(cdr, pose, North)));
	}

	std::map<core::Size, core::chemical::AA> mutation_set;
	core::Size start = ab_info->get_CDR_start(cdr, pose, North);

	for ( core::Size i = 0; i < sequence.length() -1; ++i ) {
		core::chemical::AA amino = core::chemical::aa_from_oneletter_code( sequence[ i ]);
		mutation_set[ start + i] = amino;
	}

	return mutation_set;
}

void
set_native_cdr_sequence( AntibodyInfoCOP ab_info, CDRNameEnum cdr, core::pose::Pose & pose){
	using basic::datacache::DataCache_CacheableData;
	if ( pose.data().has(CacheableDataType::NATIVE_ANTIBODY_SEQ) ) {

		//This is always a pretty line:
		//NativeAntibodySeq & seq = static_cast< NativeAntibodySeq & >(pose.data().get(CacheableDataType::NATIVE_ANTIBODY_SEQ));

		NativeAntibodySeqOP data
			=  utility::pointer::dynamic_pointer_cast< NativeAntibodySeq >
			( pose.data().get_ptr(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ) );

		if ( ! data ) {
			utility_exit_with_message("Pose said it had a NativeAntibodySeq, but it's bad!");
		}

		TR << "Setting CDR sequence in NativeAntibodySeq" << std::endl;
		data->set_from_cdr(pose, cdr );

	} else {
		pose.data().set(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ, DataCache_CacheableData::DataOP( new NativeAntibodySeq( pose, ab_info) ));
	}
}

std::string
get_native_sequence( core::pose::Pose const & pose){

	if ( ! pose.data().has(CacheableDataType::NATIVE_ANTIBODY_SEQ) ) {
		utility_exit_with_message("Pose does not have a valid NativeAntibodySeq set!");
	}

	NativeAntibodySeqCOP data
		=  utility::pointer::dynamic_pointer_cast< NativeAntibodySeq const >
		( pose.data().get_const_ptr(core::pose::datacache::CacheableDataType::NATIVE_ANTIBODY_SEQ) );

	if ( ! data ) {
		utility_exit_with_message("Pose does not have a valid NativeAntibodySeq!");
	}

	std::string seq = data->get_sequence(pose);
	//TR << "Getting sequence: " << seq << std::endl;
	return seq;
}

bool
has_native_sequence( core::pose::Pose const & pose){
	return pose.data().has(CacheableDataType::NATIVE_ANTIBODY_SEQ);

}


} //design
} //antibody
} //protocols
