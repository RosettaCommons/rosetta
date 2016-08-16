// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/denovo_design/components/PoseFolder.hh
/// @brief Given a pose with all residues, assign a 3D conformation to that pose
/// @author Tom Linsky (tlinsky@uw.edu)
/// @note   This is interface: it has no fields, and only
///         pure virtual methods.  No further constructors should
///         be defined.

#ifndef INCLUDED_protocols_denovo_design_components_PoseFolder_hh
#define INCLUDED_protocols_denovo_design_components_PoseFolder_hh

// Unit headers
#include <protocols/denovo_design/components/PoseFolder.fwd.hh>

// Protocol headers
#include <protocols/loops/Loops.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/select/residue_selector/ResidueSelector.fwd.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.fwd.hh>
#include <utility/excn/EXCN_Base.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/Tag.fwd.hh>

namespace protocols {
namespace denovo_design {
namespace components {

/// @brief Exception thrown by PoseFolder::apply to indicate folding was not successful
class EXCN_Fold : public utility::excn::EXCN_Base {
public:
	EXCN_Fold( std::string const & msg ):
		EXCN_Base(), msg_( msg ) {}

	virtual void
	show( std::ostream & os ) const { os << msg_; }

	std::string msg_;
private:
	EXCN_Fold();
};

/// @brief Given a pose with all residues, and a StructureData object,
///        assign a 3D conformation to the pose
class PoseFolder : public utility::pointer::ReferenceCount {
public: // Creation
	/// @brief Default constructor
	PoseFolder( std::string const & type );

	/// @brief Destructor
	virtual
	~PoseFolder();

	// Pure virtuals
public:
	virtual PoseFolderOP
	clone() const = 0;

	/// @brief performs folding
	/// @param pose    - The pose to be folded, with all residues added.  The pose should be prepared with
	///                  any necessary cutpoints added before giving to the PoseFolder. Torsions in the pose
	///                  should be adjusted, and no residues should be added or removed.
	/// @param movable - Subset of residues for which new backbone conformations will be sampled. Residues
	///                  specified as 'True' in movable must also be present in one or more Loops in order
	///                  to be folded. Movable's size must match pose.total_residue()
	/// @param loops   - Loops to be folded.  Cutpoints specified here must be match the cutpoints found in
	///                  the pose. Residues not within any loop should not be folded. Residues contained
	///                  in a loop but not in the movable set should not be folded.
	/// @throws EXCN_Fold if anything goes wrong in folding. Derived classes should throw this.
	virtual void
	apply(
		core::pose::Pose & pose,
		core::select::residue_selector::ResidueSubset const & movable,
		protocols::loops::Loops const & loops ) const = 0;

protected:
	virtual void
	parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data ) = 0;

	// Public functions
public:
	void
	parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data );

	std::string const &
	type() const;

private:
	std::string type_;

private:
	/// @brief Prevent direct instantiation: No other constructors allowed.
	PoseFolder() {};

}; // PoseFolder

} //protocols
} //denovo_design
} //components

#endif //INCLUDED_protocols_denovo_design_components_PoseFolder_hh
