// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody_design/AntibodyFeatures.hh
/// @brief
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)


#ifndef INCLUDED_protocols_antibody_AntibodyFeatures_hh
#define INCLUDED_protocols_antibody_AntibodyFeatures_hh

#include <protocols/features/InterfaceFeatures.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <utility/sql_database/DatabaseSessionManager.hh>
#include <protocols/antibody/metrics.hh>

namespace protocols {
namespace antibody {
// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace protocols::features;
using namespace core::scoring;

/// @brief Collects data on an antibody including CDRs, interfaces of L_H, L_A, H_A, and LH_A (this can be set), and other metrics.
/// @details Adds .... tables to the database.  See also protocols/antibody/clusters/CDRClusterFeatures
class AntibodyFeatures : public InterfaceFeatures {

public:

	AntibodyFeatures();

	AntibodyFeatures(AntibodyInfoOP ab_info);

	AntibodyFeatures(AntibodyInfoOP ab_info, ScoreFunctionCOP scorefxn);


	/// @brief return string with class name
	virtual std::string
	type_name() const;

	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	void
	write_ab_metrics_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	/// @brief Write kink metrics schema.  Please add or modify as needed.
	void
	write_ab_H3_kink_metrics_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	void
	write_cdr_metrics_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	//void
	//write_cdr_definitions_schema_to_db(
	// utility::sql_database::sessionOP db_session) const;

	void
	write_cdr_residue_schema_to_db(
		utility::sql_database::sessionOP db_session) const;


	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	//virtual utility::vector1<std::string>
	//features_reporter_dependencies() const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose);

	/// @brief collect all the feature data for the pose
	virtual core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);



	void
	report_ab_metrics_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	report_ab_H3_kink_metrics_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	report_cdr_metrics_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		CDRNameEnum const & cdr);

	void
	report_cdr_residue_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		CDRNameEnum const & cdr,
		utility::vector1< bool > const & relevant_residues);

	void
	report_cdr_residue_features_row(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		CDRNameEnum const & cdr,
		core::Size resnum,
		core::Size position);

	//void
	//report_cdr_definitions_features();

public:

	void
	set_ab_info(AntibodyInfoOP ab_info);

	/// @brief Set intermediate interface chains: Example: L_A, H_A, L_H, LH_A (A stands for antigen)
	/// @brief Any other chains, use InterfaceFeatures.
	void
	set_interface_chains(utility::vector1< std::string > const & intermediate_interfaces);

private:

	//These can be moved to metrics, but they require IAM data already.

	core::Real
	calculate_cdr_totals(CDRNameEnum const cdr, const core::pose::Pose & pose, const utility::vector1<Real> & data) const ;


	/// @brief Calculate polar dSASA using IAM results.
	core::Real
	calculate_cdr_dpSASA(CDRNameEnum const cdr);



	core::Size
	calculate_cdr_aromatic_nres(const core::pose::Pose & pose, CDRNameEnum const cdr);

	/// @brief Calculate residue atomic contacts to antigen according to this metric:
	/// @details An atomic contact is defined as at least 5 atoms of the antigen that are within 5 angstroms of an atom
	void
	calculate_residue_atomic_contacts(const core::pose::Pose & pose, const utility::vector1<bool> & residues_to_match, const utility::vector1<bool> & antigen_residues);

	/// @brief Total number of cdr atomic contacts as defined above.
	core::Size
	calculate_cdr_contacts_total(const core::pose::Pose & pose, CDRNameEnum const cdr);

	/// @brief Total number of residues making at least one atomic contact with antigen.
	/// @details
	core::Size
	calculate_cdr_contacts_nres(const core::pose::Pose & pose, CDRNameEnum const cdr);

private:

	AntibodyInfoOP ab_info_;
	utility::vector1<std::string> intermediate_interfaces_; //Intermediate interfaces: L_A, H_A, LH_A, L_H
	bool skip_antigen_reports_;

	/////Collected Data - Maybe pass this around instead?

	//InterfaceData of H_A or LH_A depending on what kind of antibody we have.
	protocols::analysis::PerResidueInterfaceData interface_data_res_;
	protocols::analysis::InterfaceData interface_data_;
	std::pair<ParatopeMetric< core::Real >, ParatopeMetric<core::Real> > paratope_sasa_;
	ParatopeMetric< core::SSize>paratope_charge_;
	vector1<core::Size> ag_ab_atomic_contacts_;

	bool include_proto_cdr4_;

};

} //antibody
} //protocols


#endif //INCLUDED_protocols_antibody_AntibodyFeatures.hh

