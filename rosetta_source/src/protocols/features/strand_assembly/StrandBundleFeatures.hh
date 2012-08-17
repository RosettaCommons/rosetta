// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :notabs=false:tabSize=4:indentsize=4:
//
// (c) copyright rosetta commons member institutions.
// (c) this file is part of the rosetta software suite and is made available under license.
// (c) the rosetta software is developed by the contributing members of the rosetta commons.
// (c) for more information, see http://www.rosettacommons.org. questions about this can be
// (c) addressed to university of washington uw techtransfer, email: license@u.washington.edu.

/// @file strandBundleFeatures.hh
/// @brief
/// @author Doo Nam Kim (started with Tim Jacobs' code)

#ifndef INCLUDED_protocols_features_strand_assembly_StrandBundleFeatures_hh
#define INCLUDED_protocols_features_strand_assembly_StrandBundleFeatures_hh

//Unit
#include <protocols/features/strand_assembly/StrandBundleFeatures.fwd.hh>

//External
#include <boost/uuid/uuid.hpp>

//Protocols
#include <protocols/features/FeaturesReporter.hh>

//Devel
#include <protocols/features/strand_assembly/StrandFragment.hh>

//Utility
#include <utility/vector1.hh>

// for string return
#include <string>

namespace protocols {
namespace features {
namespace strand_assembly {

class StrandBundleFeatures : public protocols::features::FeaturesReporter {

public:

	StrandBundleFeatures();

	void init_from_options();

	virtual
	std::string
	type_name() const  {
		return "StrandBundleFeatures";
	}

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	///@brief collect all the feature data for the pose
	virtual
	core::Size
	report_features(
					core::pose::Pose const & pose,
					utility::vector1<bool> const & relevant_residues,
					boost::uuids::uuid struct_id,
					utility::sql_database::sessionOP db_session
					);

	utility::vector1<StrandFragment> get_full_strands(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session);
	utility::vector1<StrandFragment> get_i_strand_from_full_strand_pairs(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session);
	utility::vector1<StrandFragment> get_j_strand_from_full_strand_pairs(boost::uuids::uuid struct_id, utility::sql_database::sessionOP db_session);
	
	bool find_antiparallel                   (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);
	bool find_parallel                       (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);

	core::Real check_sheet_dis_antiparallel        (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);
	core::Real check_sheet_dis_parallel            (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);

	core::Real sheet_torsion                 (core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);	
	core::Real shortest_dis_sidechain	(core::pose::Pose const & pose, StrandFragment strand_i, StrandFragment strand_j);	
	core::Real shortest_dis_sidechain	(core::Real val_shortest_dis_sidechain_1,
										 core::Real val_shortest_dis_sidechain_2,
										 core::Real val_shortest_dis_sidechain_3,
										 core::Real val_shortest_dis_sidechain_4);	

	core::Real shortest_dis_pairs	(core::Real return_of_check_sheet_dis_antiparallel_1,
									 core::Real return_of_check_sheet_dis_antiparallel_2,
									 core::Real return_of_check_sheet_dis_antiparallel_3,
									 core::Real return_of_check_sheet_dis_antiparallel_4);	
private:

	core::Size min_strand_size_;
	core::Size max_strand_size_;
	core::Real min_O_N_dis_;
	core::Real max_O_N_dis_;
	core::Real min_sheet_dis_;
	core::Real max_sheet_dis_;
	core::Real min_sheet_torsion_;
	core::Real max_sheet_torsion_;	
};

} //namespace strand_assembly
} //namespace features
} //namespace protocols

#endif
