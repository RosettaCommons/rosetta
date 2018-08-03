// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/chemically_conjugated_docking/UBQ_GTPaseMover.hh
/// @brief A mover typically used for binding ubiquitin to a substrate protein.
/// @author Steven Lewis and Hope Anderson

#ifndef INCLUDED_protocols_chemically_conjugated_docking_UBQ_GTPaseMover_HH
#define INCLUDED_protocols_chemically_conjugated_docking_UBQ_GTPaseMover_HH

// Unit headers
#include <protocols/chemically_conjugated_docking/Gp_quantification_metrics.hh>
#include <protocols/chemically_conjugated_docking/Gp_extra_bodies.hh>
#include <protocols/chemically_conjugated_docking/UBQ_GTPaseMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/analysis/InterfaceAnalyzerMover.fwd.hh>

// Core headers
#include <core/id/AtomID.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
//#include <utility/tag/XMLSchemaGeneration.fwd.hh> //transcluded from Mover

namespace protocols {
namespace chemically_conjugated_docking {

///@brief A mover typically used for binding ubiquitin to a substrate protein.
class UBQ_GTPaseMover : public protocols::moves::Mover {

public:

	/////////////////////
	/// Constructors  ///
	/////////////////////

	/// @brief Default constructor
	UBQ_GTPaseMover();

	/// @brief Copy constructor (not needed unless you need deep copies)
	UBQ_GTPaseMover( UBQ_GTPaseMover const & ) = default;

	/// @brief Destructor (important for properly forward-declaring smart-pointer members)
	~UBQ_GTPaseMover() override;

	/////////////////////
	/// Mover Methods ///
	/////////////////////

public:
	/// @brief Apply the mover
	void
	apply( core::pose::Pose & pose ) override;

	void
	show( std::ostream & output = std::cout ) const override;

	///////////////////////////////
	/// Rosetta Scripts Support ///
	///////////////////////////////

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & ) override;

	//UBQ_GTPaseMover & operator=( UBQ_GTPaseMover const & src );

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );


	// extra functions

	/// @brief sets up the pose with UBQ and GTPase
	void initialize(core::pose::Pose & GTPase);

	/// @brief applies filters and stores torsion angles near LYX bond
	void analyze_and_filter(core::pose::Pose & pose);

	/// @brief makes an internal ResidueIndexSelector given a residue index on the substrate.
	void make_index_selector(core::Size index);

	// getters / setters

	// original instance variables, primarily for internal use
	core::scoring::ScoreFunctionOP get_fullatom_scorefunction() const;

	void set_fullatom_scorefunction(core::scoring::ScoreFunctionOP const & sfxn);

	core::pack::task::TaskFactoryOP get_task_factory() const;

	void set_task_factory(core::pack::task::TaskFactoryOP const & tasks);

	core::kinematics::MoveMapOP get_amide_mm() const;

	void set_amide_mm(core::kinematics::MoveMapOP const & mm);

	protocols::loops::Loop const & get_loop() const;

	void set_loop(protocols::loops::Loop const & loop);

	utility::vector1<core::id::AtomID> const & get_atom_IDs() const;

	void set_atom_IDs(utility::vector1<core::id::AtomID> const & ids);

	core::pose::Pose const & get_starting_pose() const;

	std::string const & get_InterfaceSasaDefinition() const;

	protocols::analysis::InterfaceAnalyzerMoverOP get_IAM() const;

	void set_IAM(protocols::analysis::InterfaceAnalyzerMoverOP const & iam);

	utility::vector1<core::Size> const & get_extra_bodies_chains() const;

	void set_extra_bodies_chains(utility::vector1<core::Size> const & chains);

	bool get_extra_bodies() const;

	void set_extra_bodies(bool const hasExtra);

	core::Size get_n_tail_res() const;

	void set_n_tail_res(core::Size const numTails);

	// other getters and setters
	core::Size get_GTPase_lys() const;

	void set_GTPase_lys(core::Size const lys);

	core::Real get_scorefilter() const;

	void set_scorefilter(core::Real const maxScore);

	core::Real get_SASAfilter() const;

	void set_SASAfilter(core::Real const minSASA);

	std::string const & get_UBQpdb() const;

	void set_UBQpdb(std::string const & ubqFile);

	core::select::residue_selector::ResidueSelectorCOP get_selector() const;

	void set_selector(core::select::residue_selector::ResidueSelectorCOP const & select);

private: // methods

private: // data

	core::scoring::ScoreFunctionOP fullatom_scorefunction_;

	core::pack::task::TaskFactoryOP task_factory_;

	core::kinematics::MoveMapOP amide_mm_;
	//  core::kinematics::MoveMapOP loop_mm_;
	//  core::kinematics::MoveMapOP all_mm_;

	protocols::loops::Loop loop_;

	/// @brief vector contains atomIDs for isopeptide bond and atoms before/after bond to determine various torsions
	utility::vector1< core::id::AtomID > atomIDs;

	core::pose::Pose starting_pose_; //maintained from run to run

	std::string const InterfaceSasaDefinition_; //calculator name

	protocols::analysis::InterfaceAnalyzerMoverOP IAM_;

	core::Size GTPase_lys_ = 0; //converted to member data for sharing between setup and apply

	/// @brief used to track which chains are "extra" nonmoving bodies in extra bodies mode
	utility::vector1< core::Size > extra_bodies_chains_;

	// instance variables that were previously options

	bool extra_bodies_;

	core::Size n_tail_res_;

	core::Real scorefilter_;

	core::Real SASAfilter_;

	std::string UBQpdb_;

	core::select::residue_selector::ResidueSelectorCOP selector_;

	protocols::jd2::JobOP job_me_;
	//core::Size GTPase_residue_;
};

std::ostream &
operator<<( std::ostream & os, UBQ_GTPaseMover const & mover );

} //protocols
} //chemically_conjugated_docking

#endif //protocols_chemically_conjugated_docking_UBQ_GTPaseMover_HH
