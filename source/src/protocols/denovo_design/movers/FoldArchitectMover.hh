// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/FoldArchitectMover.hh
/// @brief Mover that builds and folds a structure via fragment insertion
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_protocols_denovo_design_movers_FoldArchitectMover_hh
#define INCLUDED_protocols_denovo_design_movers_FoldArchitectMover_hh

// Unit headers
#include <protocols/denovo_design/movers/FoldArchitectMover.fwd.hh>
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

namespace protocols {
namespace denovo_design {
namespace movers {

///@brief Mover that builds and folds a structure via fragment insertion
class FoldArchitectMover : public protocols::moves::Mover {
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
	FoldArchitectMover();

	// destructor (important for properly forward-declaring smart-pointer members)
	~FoldArchitectMover() override;

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
		basic::datacache::DataMap & data
	) override;


public:
	std::string const &
	id() const;

	/// @brief Sets the ID of this mover
	void
	set_id( std::string const & id_value );

	/// @brief Sets the architect that will be used to build blueprints dynamically
	void
	set_architect( architects::DeNovoArchitect const & architect );

	/// @brief Sets the pose folder that will be used to sample conformations.
	/// @details Clones the PoseFolder
	void
	set_folder( components::PoseFolder const & folder );

	/// @brief Sets the pose folder that will be used to sample conformations
	/// @details Uses the provided PoseFolder pointer directly
	void
	set_folder( components::PoseFolderOP folder );

	/// @brief Sets the pose builder that will be used to create the extended-chain pose
	/// @details Clones the provided PoseBuilder
	void
	set_builder( components::PoseBuilder const & builder );

	/// @brief Sets the pose builder that will be used to create the extended-chain pose
	/// @details Uses/stores the provided PoseBuilder pointer directly
	void
	set_builder( components::PoseBuilderOP builder );

	/// @brief Sets the perturber that will be used to change the plans of the architects
	/// @details Clones the provided Perturber
	void
	set_perturber( components::StructureDataPerturber const & perturber );

	/// @brief Sets the perturber that will be used to change the plans of the architects
	/// @details Uses/stores the provided Perturber pointer directly
	void
	set_perturber( components::StructureDataPerturberOP perturber );

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

	/// @brief Sets the number of residues of overlap to use between build phases
	void
	set_build_overlap( core::Size const overlap_val );

	/// @brief Sets the filter that will be used to evaluate the generated conformations. Used to return a score, not to evaluate true/false
	void
	set_score_filter( protocols::filters::Filter const & filter );

	/// @brief Reads/parses an XML file and uses it to set up the mover
	void
	setup_from_xml_file( std::string const & xml_file, basic::datacache::DataMap & data );

	/// @brief Reads/parses an XML string and uses it to set up the mover
	void
	setup_from_xml_string( std::string const & xml_string, basic::datacache::DataMap & data );

	/// @brief Reads/parses a Tag object and uses it to set up the mover
	void
	setup_from_xml_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	/// @brief Clears the list of prefold movers
	void
	clear_prefold_movers();

	/// @brief Adds a mover to the list of prefold movers
	void
	add_prefold_mover( protocols::moves::Mover const & mover );

	/// @brief Clears the list of postfold movers
	void
	clear_postfold_movers();

	/// @brief Adds a mover to the list of postfold movers
	void
	add_postfold_mover( protocols::moves::Mover const & mover );

	/// @brief Adds a filter to the list of filters that indicate whether a given folding attempt is was successful or not. All filters must pass for a successful folding attempt
	void
	add_filter( protocols::filters::Filter const & filter );

	/// @brief Clears the list of filters used to indicate whether a given folding attempt is successful
	void
	clear_filters();

	/// @brief If true, extra PDBs will be outputted for debugging purposes. Could be useful for figuring out why folding attempts are failing. If false, no extra files are outputted
	void
	set_dump_pdbs( bool const dump ) { dump_pdbs_ = dump; }

	/// @brief The debug value of the poseBuilder will be set based on this
	void
	set_debug( bool const debug );

	/// @brief Sets the number of folding attempts for each build phase
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

	static
	void
	provide_xml_schema_with_name(
		utility::tag::XMLSchemaDefinition & xsd,
		std::string const & mover_name );

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
	build_denovo_backbone_ct_naming_func( std::string const & );

	static std::string
	prefold_ct_naming_func( std::string const & );

	static std::string
	postfold_ct_naming_func( std::string const & );

	static std::string
	filters_ct_naming_func( std::string const & );

	static std::string
	connections_ct_naming_func( std::string const & );

	static std::string
	folder_ct_namer( std::string const & );

	static std::string
	perturber_ct_namer( std::string const & );

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
	parse_filters( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	void
	parse_prefold_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	void
	parse_postfold_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	MoverOPs
	parse_movers( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) const;

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
	components::PoseBuilderOP builder_;
	components::PoseFolderOP folder_;
	components::StructureDataPerturberOP perturber_;
	MoverOPs prefold_movers_;
	MoverOPs postfold_movers_;
	FilterCOPs filters_;
	protocols::filters::FilterCOP score_filter_;

private:
	// options
	std::string id_;
	bool dry_run_;
	bool dump_pdbs_;
	bool debug_;
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

// backwards compatibility
class BuildDeNovoBackboneMover : public FoldArchitectMover {
public:
	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );
};

/// @brief goes through loops and adds overlapping residues
FoldArchitectMover::ResidueVector
add_overlap_to_loops(
	protocols::loops::Loops & loops,
	core::Size const overlap,
	core::pose::Pose const & pose );

//void
//modify_for_check( components::StructureData & sd );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_FoldArchitectMover_hh
