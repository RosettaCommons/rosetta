// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/antibody/design/MutateFrameworkForCluster.cc
/// @brief
/// @author Jared Adolf-Bryfogle

#include <protocols/antibody/design/MutateFrameworkForCluster.hh>
#include <protocols/antibody/design/MutateFrameworkForClusterCreator.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/design/util.hh>
#include <protocols/antibody/clusters/util.hh>
#include <protocols/rosetta_scripts/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <utility/tag/Tag.hh>

#include <boost/algorithm/string.hpp>
#include <utility/string_util.hh>

#include <fstream>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.designMutateFrameworkForCluster");

namespace protocols {
namespace antibody {
namespace design {

using namespace core::scoring;
using namespace protocols::antibody::clusters;

MutateFrameworkForCluster::MutateFrameworkForCluster() :
	protocols::moves::Mover("MutateFrameworkForCluster"),
	ab_info_(/* NULL */)
{
	set_defaults();
	load_data();
}

MutateFrameworkForCluster::MutateFrameworkForCluster(AntibodyInfoCOP ab_info):
	protocols::moves::Mover("MutateFrameworkForCluster"),
	ab_info_(ab_info)
{
	set_defaults();
	load_data();
}

MutateFrameworkForCluster::~MutateFrameworkForCluster() {}



MutateFrameworkForCluster::MutateFrameworkForCluster(MutateFrameworkForCluster const & src) :
	protocols::moves::Mover(src),
	mutant_info_(src.mutant_info_),
	cdrs_(src.cdrs_),
	pack_shell_(src.pack_shell_),
	keep_current_(src.keep_current_)
{
	if ( src.ab_info_ ) ab_info_ = AntibodyInfoOP( new AntibodyInfo( *src.ab_info_ ));
	if ( src.scorefxn_ ) scorefxn_ = scorefxn_->clone();
	
}

void
MutateFrameworkForCluster::set_defaults() {
	cdrs_.clear();
	cdrs_.resize(6, true);
	pack_shell_ = 6.0;
	keep_current_ = true;
	scorefxn_ = get_score_function();

}

//moves::MoverOP
//MutateFrameworkForCluster::clone() const {
// return moves::MoverOP(new MutateFrameworkForCluster(*this));
//}

protocols::moves::MoverOP
MutateFrameworkForCluster::fresh_instance() const {
	return moves::MoverOP(new MutateFrameworkForCluster());
}

std::string
MutateFrameworkForCluster::get_name() const {
	return "MutateFrameworkForCluster";
}

void
MutateFrameworkForCluster::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap& data,
	const Filters_map&,
	const moves::Movers_map&,
	const Pose& pose)
{
	ab_info_ = AntibodyInfoOP( new AntibodyInfo( pose ));
	cdrs_ = get_cdr_bool_from_tag(tag, "cdrs");
	keep_current_ = tag->getOption("keep_current", keep_current_);
	pack_shell_ = tag->getOption("pack_shell", pack_shell_);

	scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, "scorefxn", data ) ;

}

void
MutateFrameworkForCluster::set_cdr_only(CDRNameEnum const & cdr){
	cdrs_.clear();
	cdrs_.resize(6, false);
	cdrs_[ cdr ] = true;
}

void
MutateFrameworkForCluster::set_cdrs(utility::vector1<bool> const & cdrs){
	cdrs_ = cdrs;
}

void
MutateFrameworkForCluster::set_pack_shell(core::Real const pack_shell){
	pack_shell_ = pack_shell;
}

void
MutateFrameworkForCluster::set_scorefxn(core::scoring::ScoreFunctionCOP scorefxn) {
	scorefxn_ = scorefxn;
}

void
MutateFrameworkForCluster::set_custom_data(const std::map<CDRClusterEnum,utility::vector1<MutantPosition> >& mutant_info){
	mutant_info_ = mutant_info;
}

void
MutateFrameworkForCluster::keep_current(bool keep_current){
	keep_current_ = keep_current;
}

std::map<clusters::CDRClusterEnum, utility::vector1<MutantPosition> >
MutateFrameworkForCluster::get_data(){
	return mutant_info_;
}

bool
MutateFrameworkForCluster::has_framework_dependant_cluster(const core::pose::Pose& pose, const CDRNameEnum cdr){
	CDRClusterEnum cluster = get_cluster_from_cache_or_ab_info(ab_info_, pose, cdr);
	if ( mutant_info_.count(cluster) > 0 ) {
		return true;
	} else {
		return false;
	}
}

bool
MutateFrameworkForCluster::has_framework_dependant_clusters(const core::pose::Pose& pose){
	for ( core::Size i = 1; i <= static_cast<core::Size>(ab_info_->get_total_num_CDRs()); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( has_framework_dependant_cluster(pose, cdr) ) {
			return true;
		} else {
			continue;
		}
	}
	return false;
}

utility::vector1<CDRClusterEnum>
MutateFrameworkForCluster::framework_dependant_clusters(){
	utility::vector1<CDRClusterEnum> clusters;
	for ( std::map<CDRClusterEnum, utility::vector1<MutantPosition> >::iterator iter = mutant_info_.begin(); iter != mutant_info_.end(); ++iter ) {
		clusters.push_back(iter->first);
	}
	return clusters;
}

utility::vector1<bool>
MutateFrameworkForCluster::framework_dependant_positions(const core::pose::Pose& pose){

	utility::vector1<bool> positions(pose.size(), false);

	for ( std::map<CDRClusterEnum, utility::vector1<MutantPosition> >::iterator iter = mutant_info_.begin(); iter != mutant_info_.end(); ++iter ) {

		for ( core::Size i = 1; i <= iter->second.size(); ++i ) {
			core::Size resnum = get_resnum_from_single_string_w_landmark(ab_info_, pose, iter->second[ i ].pdb_position_, iter->second[ i ].numbering_scheme_);
			positions[resnum] = true;
		}
	}
	return positions;
}

utility::vector1<bool>
MutateFrameworkForCluster::framework_dependant_positions(const core::pose::Pose& pose, const clusters::CDRClusterEnum cluster){

	utility::vector1<bool> positions(pose.size(), false);


	for ( core::Size i = 1; i <= mutant_info_[ cluster ].size(); ++i ) {
		core::Size resnum = get_resnum_from_single_string_w_landmark(ab_info_, pose, mutant_info_[ cluster ][ i ].pdb_position_, mutant_info_[ cluster ][ i ].numbering_scheme_);
		positions[resnum] = true;
	}
	return positions;
}

utility::vector1<bool>
MutateFrameworkForCluster::framework_dependant_mutations(core::pose::Pose const & pose, const clusters::CDRClusterEnum cluster, core::Size const resnum){

	utility::vector1<bool> mutants(20, false);
	for ( core::Size i = 1; i <= mutant_info_[ cluster ].size(); ++i ) {
		core::Size present_resnum = get_resnum_from_single_string_w_landmark(ab_info_, pose, mutant_info_[ cluster ][ i ].pdb_position_, mutant_info_[ cluster ][ i ].numbering_scheme_);
		if ( resnum == present_resnum ) {
			mutants = mutant_info_[ cluster ][ i ].mutants_allowed_;
			break;
		}
	}
	return mutants;
}

void
MutateFrameworkForCluster::load_data() {

	AntibodyEnumManager manager = AntibodyEnumManager();

	std::string line;
	std::string filename = basic::database::full_name("sampling/antibodies/design/cluster_framework_mutations.txt");
	std::fstream INFILE(filename.c_str(), std::ios::in); //Why a copy constructor with = will not work here is beyond me.  The std library is pathetic.

	if ( ! INFILE ) {
		utility_exit_with_message("Cannot open mutate cluster framework file.");
	}

	while ( getline(INFILE, line) ) {
		//TR << line << std::endl;
		utility::trim(line, "\n"); //Remove trailing line break
		boost::algorithm::trim(line); //Remove any whitespace on either side of the string

		//Continue to next line on empty string, comment
		if ( utility::startswith(line, "#") || utility::startswith(line, "\n") || line.empty()  ||  (line.find_first_not_of(' ') == std::string::npos) ) {
			continue;
		}

		utility::vector1< std::string > lineSP = utility::string_split_multi_delim(line); //Split on space or tab
		//TR << utility::to_string(lineSP) << std::endl;

		std::string cluster_str = lineSP[1];
		std::string resnum_str = lineSP[2];
		std::string aas_str = lineSP[3];
		std::string scheme_str = lineSP[4];

		//TR << cluster_str << " "<< resnum_str << std::endl;
		MutantPosition mut_pos;
		mut_pos.pdb_position_ = resnum_str;
		mut_pos.numbering_scheme_ = manager.numbering_scheme_string_to_enum(scheme_str);
		mut_pos.mutants_allowed_.clear();
		mut_pos.mutants_allowed_.resize(20, false); //Standard cannonical aas.

		for ( std::string::iterator iter_aa = aas_str.begin(); iter_aa != aas_str.end(); ++iter_aa ) {
			//TR << *iter_aa << std::endl;
			mut_pos.mutants_allowed_[ core::chemical::aa_from_oneletter_code( *iter_aa ) ] = true;

		}


		CDRClusterEnum cluster = ab_info_->get_cluster_enum(cluster_str);
		if ( mutant_info_.count(cluster) ) {
			mutant_info_[ cluster ].push_back(mut_pos);
		} else {
			mutant_info_[ cluster ];
			mutant_info_[ cluster ].push_back(mut_pos);
		}
	}

	INFILE.close();

	/*
	for (std::map<CDRClusterEnum, utility::vector1<MutantPosition> >::iterator iter = mutant_info_.begin(); iter != mutant_info_.end(); ++iter){
	TR << ab_info_->get_cluster_name(iter->first) << std::endl;
	for (core::Size i = 1; i <= mutant_info_[ iter->first ].size(); ++i){
	TR << ab_info_->get_cluster_name(iter->first) << " " << iter->second[ i ].pdb_position_ << " " << utility::to_string(iter->second[ i ].mutants_allowed_) << std::endl;
	}
	}
	*/
}


void
MutateFrameworkForCluster::apply(core::pose::Pose& pose) {

	using namespace core::pack::task;
	using namespace core::pack::task::operation;

	/*
	//Check to make sure everything is loaded correctly.
	for (std::map<CDRClusterEnum, utility::vector1<MutantPosition> >::iterator iter = mutant_info_.begin(); iter != mutant_info_.end(); ++iter){
	for (core::Size i = 1; i <= mutant_info_[ iter->first ].size(); ++i){
	TR << ab_info_->get_cluster_name(iter->first) << " " << iter->second[ i ].pdb_position_ << " " << utility::to_string(iter->second[ i ].mutants_allowed_);
	}
	}
	*/

	TaskFactoryOP tf = TaskFactoryOP( new TaskFactory());
	tf->push_back(TaskOperationCOP(new InitializeFromCommandline()) );

	PackerTaskOP task = tf->create_task_and_apply_taskoperations(pose);

	utility::vector1< bool > design_positions(pose.size(), false);

	bool framework_dependant_clusters = false;

	for ( core::Size i = 1; i <= core::Size( ab_info_->get_total_num_CDRs() ); ++i ) {
		CDRNameEnum cdr = static_cast<CDRNameEnum>( i );
		if ( ! cdrs_[ cdr ] ) continue;

		CDRClusterEnum current_cluster;

		if ( pose.data().has(core::pose::datacache::CacheableDataType::CDR_CLUSTER_INFO) ) {
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(
				pose.data().get(core::pose::datacache::CacheableDataType::CDR_CLUSTER_INFO));
			current_cluster = cluster_cache.get_cluster(cdr)->cluster();
		} else {
			current_cluster = ab_info_->get_CDR_cluster(cdr)->cluster();
		}
		if ( current_cluster  == NA ) {
			continue;
		}

		/// Cluster is found.  Now check that it actually has data we need.
		if ( ! mutant_info_.count( current_cluster ) ) continue;

		//TR << "Found cluster mutant. " << std::endl;
		framework_dependant_clusters = true;

		for ( core::Size i = 1; i <= mutant_info_[ current_cluster ].size(); ++i ) {
			MutantPosition mut_pos = mutant_info_[ current_cluster ][ i ];
			core::Size resnum = get_resnum_from_single_string_w_landmark(ab_info_, pose, mut_pos.pdb_position_, mut_pos.numbering_scheme_);
			if ( resnum == 0 ) {

				TR << "Resnum missing from pose "<< mut_pos.pdb_position_ <<
					" skipping mutation for " << ab_info_->get_cluster_name( current_cluster ) <<
					" skipping..." << std::endl;
				continue;

			}
			design_positions[ resnum ] = true;


			if ( keep_current_ ) {
				for ( core::Size aa_num = 1; aa_num <= 20; ++aa_num ) {

					if ( mut_pos.mutants_allowed_[ aa_num ] ) {

						core::chemical::AA amino = static_cast<core::chemical::AA>(aa_num);
						task->nonconst_residue_task(i).allow_aa(amino);
					}
				}
			} else {
				//Replace the current aminos
				task->nonconst_residue_task(i).restrict_absent_canonical_aas( mut_pos.mutants_allowed_);
			}

		}
	}

	//Now we get all the residues we will pack and design.
	// Currently at most we are talking about like 2 or 3 clusters that need a single mutation

	if ( ! framework_dependant_clusters ) {
		return;
	}

	core::pack::task::operation::RestrictResidueToRepacking turn_off_design;
	core::pack::task::operation::PreventRepacking turn_off_packing;

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! design_positions[ i ] ) {
			turn_off_design.include_residue( i );
		}
	}

	core::select::residue_selector::NeighborhoodResidueSelector neighbor_sel = core::select::residue_selector::NeighborhoodResidueSelector(design_positions, pack_shell_);
	utility::vector1<bool> pack_positions = neighbor_sel.apply(pose);

	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( ! pack_positions[ i ] ) {
			turn_off_packing.include_residue( i );
		}
	}

	turn_off_design.apply(pose, *task);
	turn_off_packing.apply(pose, *task);

	simple_moves::PackRotamersMover packer = simple_moves::PackRotamersMover(scorefxn_, task);
	packer.apply(pose);

}


protocols::moves::MoverOP
MutateFrameworkForClusterCreator::create_mover() const {
	return protocols::moves::MoverOP(new MutateFrameworkForCluster());
}

std::string
MutateFrameworkForClusterCreator::keyname() const {
	return MutateFrameworkForClusterCreator::mover_name();
}

std::string
MutateFrameworkForClusterCreator::mover_name() {
	return "MutateFrameworkForCluster";
}


} //design
} //antibody
} //protocols
