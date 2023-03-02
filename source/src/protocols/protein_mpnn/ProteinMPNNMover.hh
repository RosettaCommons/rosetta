// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/protein_mpnn/ProteinMPNNMover.hh
/// @brief  Mover for the ProteinMPNN PyTorch model developed by Dauparas et al.
/// @author Frederick Chan (fredchan@uw.edu)

#ifndef INCLUDED_protocols_protein_mpnn_ProteinMPNNMover_hh
#define INCLUDED_protocols_protein_mpnn_ProteinMPNNMover_hh

// Unit headers
#include <protocols/protein_mpnn/ProteinMPNNMover.fwd.hh>

// Project headers
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>
#include <protocols/moves/Mover.hh>

// Utility headers
#include <utility/vector1.hh>

namespace protocols {
namespace protein_mpnn {

/// @brief
/// This mover uses the given pose and threads the predicted sequence from the
/// ProteinMPNN model onto it.
class ProteinMPNNMover : public::protocols::moves::Mover {

public:

	/// @brief Default constructor.
	ProteinMPNNMover() = default;

	/// @brief Destructor.
	~ProteinMPNNMover() override;

	void
	apply(core::pose::Pose & pose) override;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data
	) override;

	moves::MoverOP
	clone() const override;

	moves::MoverOP
	fresh_instance() const override;

	std::string
	get_name() const override;

	static
	std::string
	mover_name();

	static
	void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Provides the citation for this module.
	void
	provide_citation_info( basic::citation_manager::CitationCollectionList & ) const override;

	///// Option setters

	/// @brief Set to true to run in deterministic mode. Will override temperature.
	void
	set_deterministic_flag(bool deterministic_flag);

	/// @brief Set the sampling temperature of the model.
	void
	set_temperature( core::Real temperature );

	/// @brief Set globally disallowed AA types. Each character is a single-letter AA code.
	/// @note "X" for unknown residue type should be included.
	void
	set_omit_AAs( utility::vector1< char > omit_AAs );

	/// @brief Set designable residues
	void
	set_design_selector_rs( core::select::residue_selector::ResidueSelectorCOP design_selector_rs );

	/// @brief Set residues whose coordinates should be passed to ProteinMPNN
	void
	set_coord_selector_rs( core::select::residue_selector::ResidueSelectorCOP coord_selector_rs );

	/// @brief Set tied positions. This is a collection of lists, each list containing residue selectors that should be tied.
	/// @details The first residues of each selector will be tied together, then the second, etc.
	///          Each residue selector must have the same number of residues.
	void
	set_tied_pos_rs( utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs );

private:

	bool deterministic_flag_ = false;

	core::Real temperature_ = 0.1;

	utility::vector1< char > omit_AAs_;

	core::select::residue_selector::ResidueSelectorCOP design_selector_rs_;

	core::select::residue_selector::ResidueSelectorCOP coord_selector_rs_;

	utility::vector1< utility::vector1< core::select::residue_selector::ResidueSelectorCOP > > tied_pos_rs_;

	core::pack::task::TaskFactoryOP task_factory_;

	void
	validate_residue_selectors(
		utility::vector1< utility::vector1< core::Size > > const & selection_positions,
		core::Size const expected_size
	);

	/// @brief Turn a bunch of residue selectors into a bunch of lists of resnums selected by those residue selectors
	utility::vector1< utility::vector1< core::Size > >
	residue_selectors_to_indices(
		core::pose::Pose const & pose,
		utility::vector1< core::select::residue_selector::ResidueSelectorCOP > const & residue_selectors
	);

};

} // protein_mpnn
} // protocols

#endif //INCLUDED_protocols_protein_mpnn_ProteinMPNNMover_hh
