// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRClusterFeatures.cc
/// @brief 
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)



#include <protocols/antibody/clusters/CDRClusterFeatures.hh>
#include <protocols/antibody/clusters/CDRClusterEnum.hh>
#include <protocols/antibody/clusters/CDRCluster.hh>

#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/AntibodyEnumManager.hh>


#include <basic/database/sql_utils.hh>
#include <basic/database/schema_generator/PrimaryKey.hh>
#include <basic/database/schema_generator/ForeignKey.hh>
#include <basic/database/schema_generator/Column.hh>
#include <basic/database/schema_generator/Schema.hh>
#include <utility/string_util.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/datacache/DataMap.hh>

namespace protocols {
namespace antibody {
namespace clusters {
	using namespace protocols::features;
	using namespace protocols::antibody;

CDRClusterFeatures::CDRClusterFeatures():
	FeaturesReporter(),
	numbering_scheme_(AHO_Scheme)
{
	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast< CDRNameEnum >(i);
		cdrs_.push_back(cdr);
	}
}

//CDRClusterFeatures::CDRClusterFeatures(CDRClusterFeatures const & src):
	//FeaturesReporter(),
	//numbering_scheme_(src.numbering_scheme_),
	//cdrs_(src.cdrs_)
//{}

CDRClusterFeatures::~CDRClusterFeatures(){}

std::string
CDRClusterFeatures::type_name() const {
	return "CDRClusterFeatures";
}

void
CDRClusterFeatures::write_schema_to_db(utility::sql_database::sessionOP db_session) const {
	using namespace basic::database::schema_generator;
	
	Column struct_id("struct_id", new DbBigInt());
	Column chain("chain", new DbText());
	Column CDR("CDR", new DbText());
	Column length("length", new DbInteger());
	Column fullcluster("fullcluster", new DbText());
	Column dis("dis", new DbReal());
	Column normDis("normDis", new DbReal());
	Column normDis_deg("normDis_deg", new DbReal());
	Column resnum_begin("resnum_begin", new DbInteger());
	Column resnum_end("resnum_end", new DbInteger());
	Column sequence("sequence", new DbText());
	
	
	Columns primary_keys;
	
	primary_keys.push_back(struct_id);
	primary_keys.push_back(resnum_begin);
	primary_keys.push_back(resnum_end);
	
	PrimaryKey primary_key(primary_keys);
	ForeignKey foreign_key(struct_id, "residues", "struct_id", true);
	
	Schema table("CDR_clusters", primary_key);
	table.add_foreign_key(foreign_key);
	
	table.add_column(chain);
	table.add_column(CDR);
	table.add_column(length);
	table.add_column(fullcluster);
	table.add_column(dis);
	table.add_column(normDis);
	table.add_column(normDis_deg);
	table.add_column(sequence);
	table.write(db_session);
}

utility::vector1<std::string>
CDRClusterFeatures::features_reporter_dependencies() const {
	utility::vector1<std::string> dependencies;
	dependencies.push_back("StructureFeatures");
	return dependencies;
}

void
CDRClusterFeatures::parse_my_tag(const utility::tag::TagCOP tag, basic::datacache::DataMap&, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const &) {
	cdrs_.clear();
	AntibodyEnumManagerOP enum_manager = new AntibodyEnumManager();
	std::string cdrs = tag->getOption< std::string >("cdrs", "L1,L2,L3,H1,H2,H3");
	std::string scheme = tag->getOption< std::string >("numbering_scheme", "AHO_Scheme");
	
	if (! enum_manager->numbering_scheme_is_present(scheme)){
		throw utility::excn::EXCN_RosettaScriptsOption( "Numbering scheme not recognized: "+ scheme);
	}
	
	set_numbering_scheme( enum_manager->numbering_scheme_string_to_enum(scheme) );
	
	vector1< std::string > cdrsSP = utility::string_split(cdrs, ',');
	for (core::Size i = 1; i <= cdrsSP.size(); ++i){
		if (! enum_manager->cdr_name_is_present(cdrsSP[i])){
			throw utility::excn::EXCN_RosettaScriptsOption("CDR not recognized: " + cdrsSP[i]);
		}
		else {
			cdrs_.push_back(enum_manager->cdr_name_string_to_enum(cdrsSP[i]));
		}
	}
	
}

void
CDRClusterFeatures::set_numbering_scheme(AntibodyNumberingSchemeEnum const & numbering_scheme) {
	numbering_scheme_ = numbering_scheme;
}

void
CDRClusterFeatures::set_cdrs_to_use(vector1<CDRNameEnum> cdrs) {
	cdrs_ = cdrs;
}

core::Size
CDRClusterFeatures::report_features(core::pose::Pose const & pose, utility::vector1< bool > const & , StructureID struct_id, utility::sql_database::sessionOP db_session) {
	using cppdb::statement;
	
	AntibodyInfoOP ab_info = new AntibodyInfo(pose, numbering_scheme_, North);
	
	//utility::vector1<bool > relavent_residues = residues;
	
	std::string stmt_string = "INSERT INTO CDR_clusters( struct_id, resnum_begin, resnum_end, chain, CDR, length, fullcluster, dis, normDis, normDis_deg, sequence) VALUES (?,?,?,?,?,?,?,?,?,?,?);";
	statement stmt(basic::database::safely_prepare_statement(stmt_string, db_session));
	
	for (core::Size i = 1; i <= 6; ++i){
		CDRNameEnum cdr = static_cast<CDRNameEnum>(i);
		CDRClusterOP cluster = ab_info->get_CDR_cluster(cdr);
		
		//Short-circuit evaluation here:
		if (ab_info->has_cluster_for_cdr(cdr) && (std::find(cdrs_.begin(), cdrs_.end(), cluster->cdr()) != cdrs_.end())){
			std::string sequence = ab_info->get_CDR_sequence_with_stem(cdr, pose, North, 0, 0);
			std::stringstream str_chain;
			str_chain << ab_info->get_CDR_chain(cdr);
			stmt.bind(1, struct_id);
			stmt.bind(2, ab_info->get_CDR_start(cdr, pose, North));
			stmt.bind(3, ab_info->get_CDR_end(cdr, pose, North));
			stmt.bind(4, str_chain);
			stmt.bind(5, ab_info->get_CDR_Name(cdr));
			stmt.bind(6, ab_info->get_CDR_length(cdr, pose, North));
			stmt.bind(7, ab_info->get_cluster_name(cluster->cluster()) );
			stmt.bind(8, cluster->distance());
			stmt.bind(9, cluster->length_normalized_distance());
			stmt.bind(10, cluster->normalized_distance_in_degrees());
			stmt.bind(11, sequence);
			basic::database::safely_write_to_database(stmt);
		}
	}
	return 0;
}

	
	

	
} //clusters
} //antibody
} //protocols
