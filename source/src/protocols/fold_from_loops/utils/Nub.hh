// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   Nub.hh
/// @brief
/// @author jaumebonet (jaume.bonet@gmail.com), Correia's LPDI/EPFL

#ifndef INCLUDED_protocols_fold_from_loops_utils_Nub_hh
#define INCLUDED_protocols_fold_from_loops_utils_Nub_hh

// Unit headers
#include <protocols/fold_from_loops/utils/Nub.fwd.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/select/residue_selector/ResidueRanges.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/movemap/MoveMapFactory.hh>
#include <core/kinematics/MoveMap.fwd.hh>
#include <core/scoring/constraints/ConstraintSet.hh>

// Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace fold_from_loops {
namespace utils {

struct NubSegment {
	core::Size order;
	utility::vector1< core::Size > coldspots = utility::vector1< core::Size >();
	core::Size c_flex = 0;
	core::Size n_flex = 0;
};

class Nub : public utility::pointer::ReferenceCount {

public:

	/// @brief Empty Constructor
	Nub();
	/// @brief Destructor
	~Nub();

	/// @brief Produces the unfolded pose
	void apply(
		core::pose::Pose const & pose,
		core::select::residue_selector::ResidueSubset const & template_selection,
		utility::vector1< core::pose::PoseOP > const & template_pieces,
		utility::vector1< std::pair< core::Size, core::Size > > const & disulfides
	);
	// @brief Transfers the unfolded conformation to the provided pose
	void transfer_unfolded_conformation( core::pose::Pose & pose );
	// @brief Retrieves the sidechains from the motif and refits them into the
	// after-abinitio pose
	void refit_motif_sidechains( core::pose::Pose & pose );

private:
	/// @brief Extracts the regions of interest from the nub
	void get_nub_pieces( core::pose::Pose & pose );
	/// @brief Attaches the template unfolded regions to the Nub to create the unfolded pose.
	void join_pieces( utility::vector1< core::pose::PoseOP > const & template_pieces );
	/// @brief Add selecter binders if provided.
	void add_binders();
	/// @brief Add labels to the unfolded pose
	void add_labels();
	/// @brief Creates the movemap to submit to ab initio protocol.
	void make_movemap();
	/// @brief Fixes and assigns disulfides to the unfolded pose
	void assign_disulfides( utility::vector1< std::pair< core::Size, core::Size > > const & disulfides );
	/// @brief Fixes fragments in case the insertion has changed the count of the residues of the template.
	void fix_fragments();
	/// @brief Fixes a FragSet
	core::fragment::FragSetOP fix_fragment( core::fragment::FragSetOP fragset );
	/// @brief Get closest binder residue to the FoldTree root
	core::Size closest_binder( core::pose::Pose const & pose );

public:
	// -- SETTERS/GETTERS -- //
	core::pose::PoseOP reference_pose() const;
	void reference_pose( core::pose::PoseOP const & pose );
	core::pose::PoseOP unfolded_pose() const;
	core::fragment::FragSetOP small_fragments() const;
	void small_fragments( core::fragment::FragSetOP fragset );
	core::fragment::FragSetOP large_fragments() const;
	void large_fragments( core::fragment::FragSetOP fragset );
	core::id::SequenceMapping template_to_unfolded_mapping() const;
	core::select::residue_selector::ResidueSelectorCOP selector() const;
	void selector( core::select::residue_selector::ResidueSelectorCOP const & selector );
	core::select::residue_selector::ResidueSelectorCOP binder_selector() const;
	void binder_selector( core::select::residue_selector::ResidueSelectorCOP const & selector );
	core::Size default_flexibility() const;
	void default_flexibility( core::Size pick );
	bool has_binder() const;
	bool is_multisegment() const;
	core::Size unfolded_jump() const;
	core::Size design_size() const;
	utility::vector1< core::pose::PoseOP > pieces() const;
	core::select::residue_selector::ResidueSubset apply_selector( core::pose::Pose const & pose ) const;
	core::select::residue_selector::ResidueRanges apply_selector_as_ranges( core::pose::Pose const & pose ) const;
	core::select::residue_selector::ResidueSubset apply_binder_selector( core::pose::Pose const & pose ) const;
	core::select::residue_selector::ResidueRanges apply_binder_selector_as_ranges( core::pose::Pose const & pose ) const;
	core::kinematics::MoveMapOP movemap() const;
	core::select::movemap::MoveMapFactoryOP movemapfactory() const;
	std::string design_chain() const;
	utility::vector1< std::pair< core::Size, core::Size > > disulfides() const;

public:
	// -- ROSETTASCRIPTS -- //
	static void provide_xml_definition( utility::tag::XMLSchemaDefinition & xsd, utility::tag::XMLSchemaSimpleSubelementList & elements );
	void parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &data, core::pose::Pose const & reference_pose );
	static std::string object_name();
	static std::string complex_object_name( std::string tag_name );
	static std::string segment_object_name();

public:
	// -- DEFAULTS -- //
	static core::select::residue_selector::ResidueSelectorCOP default_selector();

	static std::string motif_label();
	static std::string flexible_label();
	static std::string template_label();
	static std::string hotspot_label();
	static std::string coldspot_label();
	static std::string context_label();
	static std::string disulfide_label();
	static std::string contact_label();

private:
	core::pose::PoseOP reference_pose_;
	core::pose::PoseOP unfolded_pose_;
	core::fragment::FragSetOP small_;
	core::fragment::FragSetOP large_;
	core::id::SequenceMapping seqmap_;                   // Automatic, no outside access, map between the template and the unfolded pose
	core::select::residue_selector::ResidueSelectorCOP selector_;      // User provided, marks the region to keep from the Pose
	core::select::residue_selector::ResidueSelectorCOP binder_selector_;  // User provided, marks which chains of the Pose are binders
	core::select::residue_selector::ResidueSelectorCOP insertion_selector_; // Automatic, no outside access, marks the kept segments in the Unfolded Pose
	core::Size unfolded_jump_;                       // Automatic, no outside access, marks the last jump before inserting binders
	bool has_binder_;                            // Automatic, tracks if a binder was added.
	core::select::movemap::MoveMapFactoryOP movemap_;            // Automatic, generates the movemap that needs to be provided to the ab initio protocol
	utility::vector1< core::pose::PoseOP > pieces_;             // Automatic, stores the pieces of the motif.
	utility::vector1< NubSegment > piece_properties_;            // User provided, stores properties of the pieces.
	std::string design_chain_;                       // Chain ID of the design protein (chain of the first -after user order- piece of the motif).
	core::Size design_size_;                        // Automatic, size of the designed chain
	utility::vector1< std::pair< core::Size, core::Size > > disulfides_;  // Automatic, stores expected disulfides for the design chain.
	core::scoring::constraints::ConstraintSetOP working_cst_;         // Automatic, not expected to be accessed by the user
	core::Size default_flexibility_;                    // If not said otherwise, how many flexible residues per side for each segment?
	core::Size nub_disulfides_;                       // Automatic, counts internal disulfides between nub residues.
};

}
}
}

#endif
