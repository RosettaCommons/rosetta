// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/jadolfbr/cluster_utilities/identify_cdr_clusters.cc
///
/// @brief This identifies CDR clusters of antibody.  North_AHO numbering required.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

#include <devel/init.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnum.hh>
#include <protocols/antibody/util.hh>
#include <protocols/antibody/clusters/util.hh>

#include <protocols/moves/Mover.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

//Options
//#include <basic/options/option.hh>
//#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>
#include <utility/string_util.hh>

using namespace protocols::antibody;
using namespace protocols::antibody::clusters;

//Documentation:  This application identifies the CDR cluster in an antibody, renumbered with North_AHO (Used by North clusters).  Works with one PDB.  Prints to screen
//  Use: Renumber antibody using [http://dunbrack.fccc.edu/IgClassify/].  Outputs info, and appends it to a new PDB that it will write.
//  Reference:
//    PyIgClassify: a database of antibody CDR structural classifications Jared Adolf-Bryfogle; Qifang Xu; Benjamin North; Andreas Lehmann; Roland L. Dunbrack Jr Nucleic Acids Research 2014; doi: 10.1093/nar/gku1106
//    North, B., A. Lehmann, et al. (2011). JMB 406(2): 228-256.
//  Docs:
//    https://www.rosettacommons.org/docs/wiki/application_documentation/antibody/CDR-Cluster-Identification
//
class IdentifyCDRClusters : public protocols::moves::Mover{
public:
	IdentifyCDRClusters(){};

	virtual ~IdentifyCDRClusters(){};

	virtual
	std::string
	get_name() const {
		return "IdentifyCDRClusters";
	}

	void
	apply(core::pose::Pose & pose){


		AntibodyInfoOP ab_info( new AntibodyInfo(pose, North) );
		ab_info->show(std::cout);
		ab_info->setup_CDR_clusters(pose);

		std::cout << std::endl;
		for ( core::Size i = 1; i<= core::Size(ab_info->get_total_num_CDRs()); ++i ) {
			CDRNameEnum cdr_name = static_cast<CDRNameEnum>(i);
			CDRClusterCOP result = ab_info->get_CDR_cluster(cdr_name);
			std::string output = "REMARK CLUSTER "+ ab_info->get_cluster_name(result->cluster()) +" "+utility::to_string(result->normalized_distance_in_degrees());
			std::cout << output << std::endl;
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string(output);

			std::string cdr_str = ab_info->get_cluster_name(result->cluster()).substr(0,2);
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair(cdr_str+"_cluster",ab_info->get_cluster_name(result->cluster()));
			protocols::jd2::JobDistributor::get_instance()->current_job()->add_string_string_pair(cdr_str+"_distance",utility::to_string(result->normalized_distance_in_degrees()));

			check_fix_aho_cdr_numbering(ab_info, cdr_name, pose);
		}
		std::cout << std::endl <<"Info added to any echo PDB" << std::endl << std::endl;
	}
};

int main(int argc, char* argv[])
{
	try{
		devel::init(argc, argv);

		protocols::jd2::JobDistributor::get_instance()->go(protocols::moves::MoverOP( new IdentifyCDRClusters ));
		std::cout<< std::endl << std::endl << "App Author: Jared Adolf-Bryfogle; PI: Roland Dunbrack " << std::endl<< std::endl;
		std::cout<< " CDR Classification done using the same methodology as PyIgClassify." << std::endl<<std::endl;
		std::cout<< " Please cite RosettaAntibody and the following:" << std::endl;
		std::cout<< "     Adolf-Bryfogle J, Xu Q, North B, Lehmann A, Dunbrack RL. PyIgClassify: a database of antibody CDR structural classifications. Nucleic Acids Res 2015; 43:D432-D438. " << std::endl;
		std::cout<< "     North B, Lehmann A, Dunbrack RL. A new clustering of antibody CDR loop conformations. J Mol Biol 2011; 406:228-256." << std::endl <<std::endl;

	} catch(utility::excn::EXCN_Base & excn){
		std::cout << "Exception"<<std::endl;
		excn.show(std::cerr);
		return -1;
	}

	return(0);
}

//To be included:
// 1) Multiple PDB's (done through JD2)
// 2) Output to an optional text file (done through additional info at end of PDB)
// 3) Optionally Output to a file instead of PDB
