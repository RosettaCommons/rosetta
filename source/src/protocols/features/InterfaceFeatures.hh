// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/features/InterfaceFeatures.hh
/// @brief Analyzes interfaces and interface residues of a pose mainly using InterfaceAnalayzerMover.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#ifndef INCLUDED_protocols_features_INTERFACEFEATURES_HH
#define INCLUDED_protocols_features_INTERFACEFEATURES_HH

#include <protocols/features/InterfaceFeatures.fwd.hh>
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.hh>
#include <core/scoring/ScoreFunction.hh>

namespace protocols {
namespace features {

	using namespace protocols::analysis;
	using utility::vector1;


///@brief Analyzes interfaces and interface residues of a pose mainly using InterfaceAnalayzerMover.
/// By default, will analyze every interface with and report any that have dSASA > cutoff.  Interfaces to report can be set via code or RS.
///
///@details  Should work (but untested) with most ligands if loaded, rna chains, and dna chains.
/// Note that interfacial waters and ions are currently ignored, but separate features reporters may soon be in the works to accomplish this.
///
/// Most main functions are virtual so you can derive from this and have more-specific interface analysis, such as an AntibodyInterfaceFeature class.
///
class InterfaceFeatures : public FeaturesReporter {



public:
	InterfaceFeatures();

	InterfaceFeatures(core::scoring::ScoreFunctionCOP scorefxn);


	// Undefined, commenting out to fix PyRosetta build  InterfaceFeatures( InterfaceFeatures const & );

	virtual ~InterfaceFeatures();

	///@breif return string with class name
	virtual std::string
	type_name() const;

	///@brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	virtual utility::vector1<std::string>
	features_reporter_dependencies() const;

	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	///@breif collect all the feature data for the pose
	virtual core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);


	////////////////////////////////////////////////////////////////////////////
	//Options
	//

	///@brief Set the fixed chain combinations that will be analyzed.  Default is all of them.
	///@details Example:  AB_C would analyze the interface between AB and C, keeping AB fixed while separating C.
	/// Note that here, you can also give A_C   and a new pose would be created with only A and C so that this interface can be tested.
	/// Note also that multiple pose interfaces can be set.
	virtual void
	set_interface_chains(vector1< std::string > interfaces);

	///@brief Pack the interface before separation? Default is false.
	void
	set_pack_separated(bool const pack_separated );

	///@brief Pack the interface after separation?  Default is true.
	void
	set_pack_together(bool const pack_together );

	///@brief Compute the packstat score?  Default true
	void
	set_compute_packstat(bool const compute_packstat);

	void
	set_defaults();

	///@brief Set the reporter to only include interfaces >dSASA_cutoff.
	void
	set_dSASA_cutoff(core::Real dSASA_cutoff);
	
	void
	set_scorefxn(core::scoring::ScoreFunctionOP scorefxn);
	
	////////////////////////////////////////////////////////////////////////////
	virtual void
	write_interface_schema_to_db(utility::sql_database::sessionOP db_session) const;

	virtual void
	write_interface_residues_schema_to_db(utility::sql_database::sessionOP db_session) const;

	virtual void
	write_interface_side_schema_to_db(utility::sql_database::sessionOP db_session) const;
	
	
	///@brief Report all features.  Called by report_features.  Easy interface for subclassing specific interfaces.
	///@details interface is the interface analyzed, db_interface is the name that is actually inserted into the database
	/// Usually this is the same, but useful when dealing with different chain ids but same interface type. db_interface should have sides as well (L_H))
	virtual void
	report_all_interface_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string const interface,
		std::string const db_interface);
	
	///@brief Add interfaces table data to table
	virtual void
	report_interface_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string const chains_side1,
		std::string const chains_side2) const;
	
	///@brief Add interface_sides table data to table
	virtual void
	report_interface_side_features(
		core::pose::Pose const & pose,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string const chains_side1,
		std::string const chains_side2,
		protocols::analysis::InterfaceRegion const region,
		std::string const region_string) const;
	
	///@brief Add interface_residues data to table
	virtual void
	report_interface_residue_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string const chains_side1,
		std::string const chains_side2) const;

	///@brief Gets all possible interface combinations of a given pose.
	void
	make_interface_combos(core::pose::Pose const & pose, vector1<std::string> & interfaces);

private:
	
	void
	write_interface_residue_data_row_to_db(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session,
		std::string const chains_side1,
		std::string const chains_side2,
		std::string const side,
		core::Size const resnum,
		protocols::analysis::PerResidueInterfaceData const & interface_data) const;


	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	///@brief Recursive function. Get all orders of ex ABCD.
	void
	get_all_string_combos(std::string & interface, std::string current, vector1<std::string> & chains) const;

	void
	get_length_combos(std::string current, vector1<std::string> & sizes) const;
	
	std::string
	get_all_pose_chains(core::pose::Pose const & pose);
	
protected:
	
	bool
	interface_exists(vector1<std::string> & interfaces, std::string & dock_chains) const;

	bool
	chains_exist_in_pose(core::pose::Pose const & pose, std::string const interface) const;
	
	protocols::analysis::InterfaceAnalyzerMoverOP interface_analyzer_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	vector1< std::string > interfaces_;

	bool pack_together_;
	bool pack_separated_;
	bool compute_packstat_;
	core::Real dSASA_cutoff_;

};
}
}

#endif	//#ifndef INCLUDED_protocols/antibody_design_INTERFACEFEATURES_HH
