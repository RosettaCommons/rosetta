// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

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
#include <protocols/denovo_design/components/StructureData.fwd.hh>
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
#include <utility/excn/EXCN_Base.hh>

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
	typedef core::select::residue_selector::ResidueVector ResidueVector;
	typedef protocols::constraint_generator::ConstraintGeneratorOP ConstraintGeneratorOP;
	typedef protocols::constraint_generator::ConstraintGeneratorCOP ConstraintGeneratorCOP;
	typedef protocols::constraint_generator::ConstraintGeneratorCOPs ConstraintGeneratorCOPs;

public:
	BuildDeNovoBackboneMover();

	// destructor (important for properly forward-declaring smart-pointer members)
	virtual ~BuildDeNovoBackboneMover();

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	fresh_instance() const;

	/// @brief required in the context of the parser/scripting scheme
	virtual protocols::moves::MoverOP
	clone() const;

	virtual std::string
	get_name() const;

	virtual void
	apply( core::pose::Pose & pose );

	/// @brief parse XML tag (to use this Mover in Rosetta Scripts)
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		protocols::filters::Filters_map const & filters,
		protocols::moves::Movers_map const & movers,
		core::pose::Pose const & pose );

	static std::string
	class_name();

public:
	std::string const &
	id() const;

	void
	set_id( std::string const & id_value );

	void
	set_architect( architects::DeNovoArchitect const & architect );

	void
	set_folder( components::PoseFolder const & folder );

	void
	set_score_filter( protocols::filters::Filter const & filter );

	void
	set_connection_overlap( core::Size const overlap_val );

private:
	FoldScore
	folding_score( core::pose::Pose const & pose ) const;

	void
	fold( core::pose::Pose & pose ) const;

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
	clear_prefold_movers();

	void
	add_prefold_mover( protocols::moves::Mover const & mover );

	void
	clear_postfold_movers();

	void
	add_postfold_mover( protocols::moves::Mover const & mover );

	void
	apply_movers( MoverOPs const & movers, core::pose::Pose & pose ) const;

	SegmentNames
	find_roots( core::pose::Pose const & pose, protocols::loops::Loops const & loops ) const;

	ConstraintGeneratorCOPs
	create_constraint_generators( ResidueVector const & residues ) const;

	ConstraintGeneratorOP
	create_coordinate_constraint_generator( ResidueVector const & residues ) const;

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
	parse_prefold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers );

	void
	parse_postfold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers );

	MoverOPs
	parse_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & movers ) const;

private:
	// objects
	architects::DeNovoArchitectCOP architect_;
	components::ExtendedPoseBuilderCOP builder_;
	components::PoseFolderCOP folder_;
	MoverOPs prefold_movers_;
	MoverOPs postfold_movers_;
	protocols::filters::FilterCOP score_filter_;

private:
	// options
	std::string id_;
	bool dry_run_;
	core::Size connection_overlap_;
};

class EXCN_NothingToFold : public utility::excn::EXCN_Base {
public:
	EXCN_NothingToFold( std::string const & msg ):
		utility::excn::EXCN_Base(), msg_( msg ) {}

	virtual void
	show( std::ostream & os ) const { os << msg_ << std::endl; }

private:
	std::string msg_;
	EXCN_NothingToFold() {};
};

/// @brief goes through loops and adds overlapping residues
BuildDeNovoBackboneMover::ResidueVector
add_overlap_to_loops(
	protocols::loops::Loops & loops,
	core::Size const overlap,
	core::pose::Pose const & pose );

} //protocols
} //denovo_design
} //movers

#endif //protocols/denovo_design/movers_BuildDeNovoBackboneMover_hh

