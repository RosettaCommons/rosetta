// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh
/// @brief Mover that builds and folds a structure via fragment insertion
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_BuildDeNovoBackboneMover_hh
#define INCLUDED_protocols_denovo_design_movers_BuildDeNovoBackboneMover_hh

// Unit headers
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.fwd.hh>
#include <protocols/moves/Mover.hh>

// Protocol headers
#include <protocols/constraint_generator/ConstraintGenerator.fwd.hh>
#include <protocols/denovo_design/architects/DeNovoArchitect.fwd.hh>
#include <protocols/denovo_design/components/DivideAndConqueror.fwd.hh>
#include <protocols/denovo_design/components/SegmentPairing.fwd.hh>
#include <protocols/denovo_design/components/StructureData.fwd.hh>
#include <protocols/denovo_design/components/StructureDataPerturber.fwd.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.fwd.hh>
#include <protocols/denovo_design/components/PoseFolder.fwd.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.fwd.hh>
#include <protocols/denovo_design/types.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/loops/Loops.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueVector.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/excn/Exceptions.hh>

// C++ headers
#include <set>

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Mover that builds and folds a structure via fragment insertion
class BuildDeNovoBackboneMover : public protocols::moves::Mover {
public:
	typedef core::Real FoldScore;
	typedef utility::vector1< connection::ConnectionArchitectCOP > ConnectionArchitectCOPs;
	typedef utility::vector1< protocols::moves::MoverOP > MoverOPs;
	typedef utility::vector1< protocols::filters::FilterCOP > FilterCOPs;
	typedef core::select::residue_selector::ResidueVector ResidueVector;
	typedef protocols::constraint_generator::ConstraintGeneratorOP ConstraintGeneratorOP;
	typedef protocols::constraint_generator::ConstraintGeneratorCOP ConstraintGeneratorCOP;
	typedef protocols::constraint_generator::ConstraintGeneratorCOPs ConstraintGeneratorCOPs;
	typedef components::RegisterShift RegisterShift;
	typedef components::StrandOrientation StrandOrientation;

public:
	BuildDeNovoBackboneMover();

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~BuildDeNovoBackboneMover();

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	fresh_instance() const override;

	/// @brief required in the context of the parser/scripting scheme
	protocols::moves::MoverOP
	clone() const override;


	void
	apply( core::pose::Pose & pose ) override;

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose ) override;


public:
	std::string const &
	id() const;

	void
	set_id( std::string const & id_value );

	void
	set_architect( architects::DeNovoArchitect const & architect );

	void
	set_folder( components::PoseFolder const & folder );

	/// @brief sets names of segments to be included in the starting build phase
	/// @param[in] segments_csv Comma-separated string containing segment names
	void
	set_start_segments( std::string const & segments_csv );

	/// @brief sets names of segments to be included in the starting build phase
	/// @param[in] segments Set of segment names
	void
	set_start_segments( SegmentNameSet const & segments );

	/// @brief sets names of segments to be included in the final build phase
	/// @param[in] segments_csv Comma-separated string containing segment names
	void
	set_stop_segments( std::string const & segments_csv );

	/// @brief sets names of segments to be included in the final build phase
	/// @param[in] segments Set of segment names
	void
	set_stop_segments( SegmentNameSet const & segments );

	void
	set_build_overlap( core::Size const overlap_val );

	void
	set_score_filter( protocols::filters::Filter const & filter );

	void
	clear_prefold_movers();

	void
	add_prefold_mover( protocols::moves::Mover const & mover );

	void
	clear_postfold_movers();

	void
	add_postfold_mover( protocols::moves::Mover const & mover );

	void
	set_iterations_per_phase( core::Size const niter );

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

private:
	FoldScore
	folding_score( core::pose::Pose const & pose ) const;

	void
	fold( core::pose::Pose & pose ) const;

	/// @brief    runs user-provided filters on the pose
	/// @details  Filters are called for each build phase, and for every folding
	///           attempt, in order. If a filter fails, an exception is
	///           thrown and no later filters are run.
	/// @throws   EXCN_FilterFailed if any filter fails
	void
	check_pose( core::pose::Pose const & pose ) const;

	void
	check_and_accept(
		core::pose::Pose const & orig,
		core::pose::Pose & pose,
		FoldScore & bestscore ) const;

	protocols::loops::Loops
	create_loops(
		core::pose::Pose const & pose,
		ResidueVector & overlap_residues,
		SegmentNames & movable_segments ) const;

	core::select::residue_selector::ResidueSubset
	select_movable_residues(
		core::pose::Pose const & pose,
		SegmentNames const & movable_segments,
		ResidueVector const & overlap_residues ) const;

	SegmentNames
	segments_to_fold( components::StructureData const & sd ) const;

	MoverOPs
	prefold_movers(
		core::pose::Pose const & pose,
		protocols::loops::Loops const & loops,
		ConstraintGeneratorCOPs const & generators ) const;

	MoverOPs
	postfold_movers(
		core::pose::Pose const & pose,
		protocols::loops::Loops const & loops,
		ConstraintGeneratorCOPs const & generators ) const;

	void
	apply_movers( MoverOPs const & movers, core::pose::Pose & pose ) const;

	SegmentNames
	find_roots( core::pose::Pose const & pose, protocols::loops::Loops const & loops ) const;

	ConstraintGeneratorCOPs
	create_constraint_generators( ResidueVector const & residues ) const;

	ConstraintGeneratorOP
	create_overlap_constraint_generator( ResidueVector const & residues ) const;

	void
	remove_cutpoints( components::StructureData & sd, protocols::loops::Loops const & loops ) const;

private:
	// defaults
	static protocols::filters::FilterOP
	default_score_filter();

	void
	parse_architect( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	void
	parse_folder( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	void
	parse_perturber( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	static std::string
	build_denovo_backbone_ct_naming_func( std::string );

	static std::string
	prefold_ct_naming_func( std::string );

	static std::string
	postfold_ct_naming_func( std::string );

	static std::string
	filters_ct_naming_func( std::string );

	static std::string
	connections_ct_naming_func( std::string );

	static std::string
	folder_ct_namer( std::string );

	static std::string
	perturber_ct_namer( std::string );

	static void
	define_folder_group( utility::tag::XMLSchemaDefinition & );

	static void
	define_perturber_group( utility::tag::XMLSchemaDefinition & );

	static void
	define_filters_ct( utility::tag::XMLSchemaDefinition & );

	static void
	define_prefold_movers_ct( utility::tag::XMLSchemaDefinition & );

	static void
	define_postfold_movers_ct( utility::tag::XMLSchemaDefinition & );

	static std::string
	perturber_group_name();

	static std::string
	folder_group_name();

	void
	parse_filters( utility::tag::TagCOP tag, protocols::filters::Filters_map const & filter_map );

	void
	parse_prefold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers );

	void
	parse_postfold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers );

	MoverOPs
	parse_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers ) const;

	/// @brief builds/folds pose in phases using recursive algorithm
	/// @throws EXCN_Fold if we couldn't fold the pose
	core::pose::PoseOP
	build_in_phases(
		components::StructureData const & full_sd,
		components::BuildPhases const & phases,
		SegmentNameSet const & finished,
		core::Size const phase_num ) const;

	void
	fold_attempt( core::pose::Pose & pose ) const;




private:
	// objects
	architects::DeNovoArchitectCOP architect_;
	components::ExtendedPoseBuilderCOP builder_;
	components::PoseFolderCOP folder_;
	components::StructureDataPerturberCOP perturber_;
	MoverOPs prefold_movers_;
	MoverOPs postfold_movers_;
	FilterCOPs filters_;
	protocols::filters::FilterCOP score_filter_;

private:
	// options
	std::string id_;
	bool dry_run_;
	bool dump_pdbs_;
	core::Size build_overlap_;
	core::Size iterations_per_phase_;
	SegmentNameSet start_segments_;
	SegmentNameSet stop_segments_;
};

// helper mover
class SetPoseSecstructFromStructureDataMover : public protocols::moves::Mover {
public:
	static std::string const
	class_name() { return "SetPoseSecstructFromStructureDataMover"; }

	protocols::moves::MoverOP
	clone() const override;

	void
	apply( core::pose::Pose & pose ) override;

	std::string
	get_name() const override;
};

class EXCN_FilterFailed : public utility::excn::Exception {
public:
	EXCN_FilterFailed(char const *file, int line, std::string const & filter, core::Size const filter_num )
	: utility::excn::Exception(file, line, ""), filter_( filter ), filter_num_( filter_num ) {};

	std::string const &
	filter_name() const { return filter_; }

	core::Size
	filter_num() const { return filter_num_; }

private:
	std::string filter_;
	core::Size filter_num_;
};

class EXCN_NothingToFold : public utility::excn::Exception {
public:
	using utility::excn::Exception::Exception;
};

/// @brief goes through loops and adds overlapping residues
BuildDeNovoBackboneMover::ResidueVector
add_overlap_to_loops(
	protocols::loops::Loops & loops,
	core::Size const overlap,
	core::pose::Pose const & pose );

//void
//modify_for_check( components::StructureData & sd );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_BuildDeNovoBackboneMover_hh
