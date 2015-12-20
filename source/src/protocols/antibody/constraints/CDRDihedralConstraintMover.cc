// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/constraints/CDRDihedralConstraintMover.cc
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/constraints/CDRDihedralConstraintMover.hh>
#include <protocols/antibody/constraints/CDRDihedralConstraintMoverCreator.hh>
#include <protocols/antibody/clusters/CDRClusterEnumManager.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/constraints/util.hh>
#include <protocols/antibody/util.hh>

#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

#include <numeric/conversions.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/DataCache.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/tag/Tag.hh>

#include <iostream>
#include <fstream>
#include <cctype>

#include <boost/algorithm/string.hpp>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR("protocols.antibody.constraints.CDRDihedralConstraintMover");

namespace protocols {
namespace antibody {
namespace constraints {
using utility::vector1;
using namespace protocols::moves;
using namespace protocols::antibody::clusters;

CDRDihedralConstraintMover::CDRDihedralConstraintMover() :
	protocols::moves::Mover("CDRDihedralConstraintMover"),
	ab_info_(/*NULL*/)
{
	set_defaults();
	read_command_line_options();
}

CDRDihedralConstraintMover::CDRDihedralConstraintMover(AntibodyInfoCOP ab_info):
	protocols::moves::Mover("CDRDihedralConstraintMover"),
	ab_info_(ab_info)
{
	set_defaults();
	read_command_line_options();
}

CDRDihedralConstraintMover::CDRDihedralConstraintMover(AntibodyInfoCOP ab_info, CDRNameEnum cdr):
	protocols::moves::Mover("CDRDihedralConstraintMover"),
	ab_info_(ab_info)
{
	set_defaults();
	read_command_line_options();
	set_cdr(cdr);

}

CDRDihedralConstraintMover::~CDRDihedralConstraintMover() {}

CDRDihedralConstraintMover::CDRDihedralConstraintMover(CDRDihedralConstraintMover const & src) :
	protocols::moves::Mover(src),
	ab_info_(src.ab_info_),
	cdr_(src.cdr_),
	db_base_path_(src.db_base_path_),
	cdr_is_set_(src.cdr_is_set_),
	forced_cluster_(src.forced_cluster_),
	force_cluster_(src.force_cluster_),
	use_cluster_csts_(src.use_cluster_csts_),
	use_outliers_(src.use_outliers_),
	use_mean_cst_data_(src.use_mean_cst_data_),
	use_general_csts_on_failure_(src.use_general_csts_on_failure_),
	use_cluster_for_H3_(src.use_cluster_for_H3_),
	ignore_pose_datacache_(src.ignore_pose_datacache_),
	cluster_data_cutoff_(src.cluster_data_cutoff_),
	general_phi_sd_(src.general_phi_sd_),
	general_psi_sd_(src.general_psi_sd_)
{

}

void
CDRDihedralConstraintMover::set_defaults() {

	db_base_path_ = "sampling/antibodies/cluster_based_constraints/CircularHarmonic/";
	force_cluster_ = false;
	cdr_is_set_ = false;
	use_cluster_csts_ = true;
	use_general_csts_on_failure_ = true;
	use_cluster_for_H3_ = false;
	use_outliers_ = false;
	ignore_pose_datacache_ = false;

}

void
CDRDihedralConstraintMover::read_command_line_options(){
	using namespace  basic::options;
	use_mean_cst_data_ = option[OptionKeys::antibody::use_mean_cluster_cst_data]();
	cluster_data_cutoff_ = option[OptionKeys::antibody::cluster_csts_stats_cutoff]();
	general_phi_sd_ = option[OptionKeys::antibody::general_dihedral_cst_phi_sd]();
	general_psi_sd_ = option[OptionKeys::antibody::general_dihedral_cst_psi_sd]();

	bool use_outliers = option[OptionKeys::antibody::design::use_outliers]();
	bool force_outliers = option[OptionKeys::antibody::force_use_of_cluster_csts_with_outliers]();
	if ( use_outliers || force_outliers ) use_outliers_ = true;
}

void
CDRDihedralConstraintMover::parse_my_tag(
	TagCOP tag,
	basic::datacache::DataMap & ,
	Filters_map const & ,
	moves::Movers_map const & ,
	Pose const &
){
	AntibodyEnumManager manager = AntibodyEnumManager();
	clusters::CDRClusterEnumManager cluster_manager = clusters::CDRClusterEnumManager();

	if ( ! tag->hasOption("cdr") ) {
		utility_exit_with_message("Must have cdr option to add constraints");
	}

	set_cdr( manager.cdr_name_string_to_enum( tag->getOption< std::string >("cdr") ) );

	use_cluster_csts_ = tag->getOption< bool >("use_cluster_csts", use_cluster_csts_);
	use_general_csts_on_failure_ = tag->getOption< bool >("use_general_csts_on_failure", use_general_csts_on_failure_);

	if ( tag->hasOption("force_cluster") ) {
		set_force_cluster(cluster_manager.cdr_cluster_string_to_enum(tag->getOption< std::string >("force_cluster")));
	}

	cluster_data_cutoff_ = tag->getOption< core::Real >("cluster_data_required", cluster_data_cutoff_);

	use_outliers_ = tag->getOption< core::Real >("use_outliers", use_outliers_);

	general_phi_sd_ = tag->getOption< core::Real >("general_phi_sd", general_phi_sd_);
	general_psi_sd_ = tag->getOption< core::Real >("general_psi_sd", general_psi_sd_);

}

void
CDRDihedralConstraintMover::set_cdr(CDRNameEnum cdr){
	cdr_ = cdr;
	cdr_is_set_ = true;
}

void
CDRDihedralConstraintMover::set_use_cluster_csts(bool use_cluster_csts) {
	use_cluster_csts_ = use_cluster_csts;
}

void
CDRDihedralConstraintMover::set_use_general_csts_on_cluster_failure(bool use_general_csts_on_failure){
	use_general_csts_on_failure_ = use_general_csts_on_failure;
}

void
CDRDihedralConstraintMover::set_cluster_csts_data_cutoff(core::Size cutoff){
	cluster_data_cutoff_ = cutoff;
}

void
CDRDihedralConstraintMover::set_cluster_csts_use_mean_cst_data(bool use_mean_cst_data){
	use_mean_cst_data_ = use_mean_cst_data;
}

void
CDRDihedralConstraintMover::set_cluster_csts_use_outlier_data(bool use_outlier_data) {
	use_outliers_ = use_outlier_data;
}

void
CDRDihedralConstraintMover::set_force_cluster(clusters::CDRClusterEnum cluster) {
	forced_cluster_ = cluster;
	force_cluster_ = true;
}

void
CDRDihedralConstraintMover::set_remove_any_set_forced_cluster(){
	force_cluster_ = false;
}

void
CDRDihedralConstraintMover::set_use_cluster_for_H3(bool use_cluster_for_H3) {
	use_cluster_for_H3_ = use_cluster_for_H3;
}

void
CDRDihedralConstraintMover::set_general_phi_sd(core::Real phi_sd){
	general_phi_sd_ = phi_sd;
}

void
CDRDihedralConstraintMover::set_general_psi_sd(core::Real psi_sd) {
	general_psi_sd_ = psi_sd;
}

void
CDRDihedralConstraintMover::set_ignore_pose_datacache(bool ignore_pose_datacache){
	ignore_pose_datacache_ = ignore_pose_datacache;
}

//void
//CDRDihedralConstraintMover::parse_my_tag(
//  TagCOP tag,
//  basic::datacache::DataMap&,
//  const Filters_map&,
//  const Movers_map&,
//  const Pose&)
//{
//
//}

void
CDRDihedralConstraintMover::apply(core::pose::Pose& pose) {
	using namespace core::pose::datacache;

	if ( ! ab_info_ ) {
		ab_info_ = AntibodyInfoOP(new AntibodyInfo(pose));
	}

	if ( ! cdr_is_set_ ) throw utility::excn::EXCN_Msg_Exception("CDR not set for CDRDihedralConstraintMover!");

	if ( use_cluster_csts_ && (cdr_ != h3 || use_cluster_for_H3_) && (cdr_ != l4 || cdr_ != h4) )  {

		/////////////  Use Cluster from set cluster, Datacache, or AntibodyInfo
		CDRClusterEnum local_cluster;

		if ( force_cluster_ ) {
			local_cluster = forced_cluster_;
		} else if ( pose.data().has(CacheableDataType::CDR_CLUSTER_INFO) && (! ignore_pose_datacache_) ) {
			BasicCDRClusterSet const & cluster_cache = static_cast< BasicCDRClusterSet const & >(pose.data().get(CacheableDataType::CDR_CLUSTER_INFO));
			local_cluster  = cluster_cache.get_cluster(cdr_)->cluster();
		} else {
			local_cluster = ab_info_->get_CDR_cluster(cdr_)->cluster();
		}



		/////////////  Add the Constraint to the pose
		core::Size cluster_stats = get_number_of_struct_used_for_csts(local_cluster);
		TR << ab_info_->get_cluster_name(local_cluster) << " data: "<< cluster_stats <<" cutoff: " << cluster_data_cutoff_ << std::endl;
		bool constraint_add_successful;
		if ( cluster_stats >= cluster_data_cutoff_ ) {
			constraint_add_successful = add_harmonic_cluster_constraint( pose, local_cluster);
		} else {
			constraint_add_successful = false;
		}


		if ( ! constraint_add_successful && use_general_csts_on_failure_ ) {
			TR << "Adding general dihedral constraints for "<< ab_info_->get_CDR_name( cdr_ ) << std::endl;
			add_harmonic_dihedral_cst_general(ab_info_, pose, cdr_, general_phi_sd_, general_psi_sd_);
		}

	} else {
		TR << "Adding general dihedral constraints for "<< ab_info_->get_CDR_name( cdr_ ) << std::endl;
		add_harmonic_dihedral_cst_general(ab_info_, pose, cdr_, general_phi_sd_, general_psi_sd_);
	}

}

bool
CDRDihedralConstraintMover::add_harmonic_cluster_constraint(core::pose::Pose & pose, CDRClusterEnum const cluster){

	using namespace core::scoring::constraints;

	if ( ab_info_->get_current_AntibodyNumberingScheme() != "AHO_Scheme" ) {
		throw utility::excn::EXCN_Msg_Exception("CDRDihedralConstraintMover with cluster-based constraints "
			"only works with antibodies using the AHO_Scheme for numbering\n"
			"Please use a properly renumbered antibody or set the option set_use_cluster_csts(false) in the class");
	}


	std::string fname = get_harmonic_cluster_constraint_filename(cluster);
	if ( fname=="NA" ) { return false;}
	try {
		ConstraintSetOP cst = ConstraintIO::get_instance()->read_constraints(fname, ConstraintSetOP( new ConstraintSet ), pose);

		pose.add_constraints(cst->get_all_constraints());
		return true;
	}
catch(utility::excn::EXCN_Exception &excn){
	TR<< "Problem adding dihedral constraints for CDR cluster." <<std::endl;
	std::cerr << "Exception : " << std::endl;
	excn.show( std::cerr );
	excn.show( TR );
	return false;
}

}

core::Size
CDRDihedralConstraintMover::get_number_of_struct_used_for_csts(CDRClusterEnum const cluster){

	using namespace basic::options;

	//TODO: Refactor this to be read once at construction for each type of constraint.
	std::string path = get_harmonic_cluster_constraint_db_directory();
	std::string extension = ".txt";
	std::string specific_path = path + "/"+"MEAN_SD" + extension;
	std::string fname = basic::database::full_name( specific_path );

	if ( !utility::file::file_exists(fname) ) {
		throw utility::excn::EXCN_Msg_Exception(" "+fname+" does not exist.  Cannot load load dihedral cst mean_sd data");
	}

	std::string line;
	utility::io::izstream mean_sds(fname);
	if ( mean_sds.bad() ) {
		utility_exit_with_message("Unable to open "+fname);
	}

	core::Size nstruct = 0;

	while ( getline(mean_sds, line) ) {

		//Skip any comments + empty lines
		utility::trim(line, "\n"); //Remove trailing line break
		boost::algorithm::trim(line); //Remove any whitespace on either side of the string

		//Continue to next line on empty string, comment
		if ( utility::startswith(line, "#") || utility::startswith(line, "\n") || line.empty()  ||  (line.find_first_not_of(' ') == std::string::npos) ) {
			continue;
		}

		utility::vector1< std::string > lineSP = utility::string_split_multi_delim(line); //Split on space or tab
		if ( ab_info_->get_cluster_name(cluster) == lineSP[1] ) {
			nstruct = utility::string2Size(lineSP[4]);
			break;
		}
	}
	mean_sds.close();
	return nstruct;
}

std::string
CDRDihedralConstraintMover::get_harmonic_cluster_constraint_filename(CDRClusterEnum const cluster){

	using namespace basic::options;

	std::string cluster_type = ab_info_->get_cluster_name(cluster);
	if ( cluster_type=="NA" ) {
		TR<< "Cannot add cluster dihedral constraint to cdr cluster of type NA.  Skipping."<<std::endl;
		return "NA";
	}
	std::string path = get_harmonic_cluster_constraint_db_directory();
	std::string extension = ".txt";
	std::string specific_path = path + "/"+cluster_type + extension;
	std::string fname = basic::database::full_name( specific_path );
	if ( !utility::file::file_exists(fname) ) {
		TR<< "Fname "<<fname<<" Does not exist.  No constraint will be added."<<std::endl;
		return "NA";
	}
	return fname;
}

std::string
CDRDihedralConstraintMover::get_harmonic_cluster_constraint_db_directory() {
	std::string extension  = use_outliers_ ? "outliers_true" : "outliers_false_liberal";

	if ( ! use_mean_cst_data_ ) {
		return db_base_path_ + extension;
	} else {
		return db_base_path_ + extension + "_use_means";
	}
}
protocols::moves::MoverOP
CDRDihedralConstraintMover::clone() const {
	CDRDihedralConstraintMoverOP ptr(new CDRDihedralConstraintMover(*this));
	return ptr;
}

//CDRDihedralConstraintMover & operator=(CDRDihedralConstraintMover const & src) {
// return CDRDihedralConstraintMover(src);
//}

protocols::moves::MoverOP
CDRDihedralConstraintMover::fresh_instance() const {
	CDRDihedralConstraintMoverOP ptr(new CDRDihedralConstraintMover);
	return ptr;
}

std::string
CDRDihedralConstraintMover::get_name() const {
	return "CDRDihedralConstraintMover";
}

protocols::moves::MoverOP
CDRDihedralConstraintMoverCreator::create_mover() const {
	CDRDihedralConstraintMoverOP ptr(new CDRDihedralConstraintMover);
	return ptr;
}

std::string
CDRDihedralConstraintMoverCreator::keyname() const {
	return CDRDihedralConstraintMoverCreator::mover_name();
}

std::string
CDRDihedralConstraintMoverCreator::mover_name() {
	return "CDRDihedralConstraintMover";
}

} //constraints
} //antibody
} //protocols
