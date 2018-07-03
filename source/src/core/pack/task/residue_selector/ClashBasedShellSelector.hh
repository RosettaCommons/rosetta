// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/task/residue_selector/ClashBasedShellSelector.hh
/// @brief  The ClashBasedShellSelector identifies all residues that clash with at least one rotamer of a design position
/// @author Noah Ollikainen (nollikai@gmail.com)
/// @author Roland A. Pache, PhD
/// @author Kale Kundert (kale@thekunderts.net)

#ifndef INCLUDED_core_pack_task_residue_selector_ClashBasedShellSelector_HH
#define INCLUDED_core_pack_task_residue_selector_ClashBasedShellSelector_HH

// Unit headers
#include <core/pack/task/residue_selector/ClashBasedShellSelector.fwd.hh>

// Package headers
#include <core/types.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pack/task/TaskFactory.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility Headers
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <utility/vector1.hh>

// C++ headers
#include <set>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

/// @brief Define a shell based on potential sidechain clashes with a set of
/// residues being repacked or designed.
///
/// The most common reason to use a clash-based shell is to pick a minimal set
/// of residues to repack around some positions you want to design.  For most
/// simple cases, the ClashBasedRepackShell task operation is the best way to
/// do this.  Use ClashBasedShellSelector if you want to combine the repack
/// shell with other selections (e.g. repack everything that's in the first two
/// repack shells but not more than 6Å from a certain residue), or do something
/// else entirely.
///
/// The best way to specify which residues to build the shell around (called
/// the "focus" residues, see set_focus()) is to use a TaskFactory.  This is
/// because task factories know which rotamers go in each position, so each one
/// can be checked for clashes.  You can also specify focus residues using a
/// ResidueSelector.  In this case rotamers will be generated for only the
/// native amino acid (i.e. assuming the positions are being repacked) and
/// checked for clashes as usual.
///
/// If you ask for more than one shell (see set_num_shells()), the shell
/// calculation will be iterated, with the residues found at each step included
/// as inputs for the next.  The residues added in this manner are treated as
/// repackable for the purpose of making shells (i.e. only rotamers of the
/// native amino acid are considered).
///
/// This class mainly exists to provide a RosettaScripts interface for
/// selecting clash-based shells.  If you want to select such a shell from C++,
/// the find_clashing_shells() function will be simpler and easier to use.
///
/// In order for a residue the be included in the clash-based shell, its
/// sidechain must have a steric clash with one of the rotamers allowed for one
/// of the focus residues.  Furthermore, that rotamer must not have any
/// backbone clashes with any other residue.  Backbone clashes can't be
/// resolved by the packer, so there's no point considering rotamers that clash
/// with the backbone.  Clashes are defined entirely by inter-atomic distance
/// and do not depend on any particular score terms.
///
/// Think twice before using ClashBasedShellSelector or ClashBasedRepackShell
/// if you're running many simulations with different backbones, because you
/// may end up picking different shells in different simulations.  Repacking
/// different sets of residues can make scores hard to compare, especially if
/// your input structure isn't relaxed.
///
/// If you want to understand exactly why ClashBasedShellSelector either
/// included or excluded a particular residue, run rosetta with the
/// `-out:levels core.pack.task.residue_selector.ClashBasedShellSelector:trace`
/// flag.
class ClashBasedShellSelector
	: public core::select::residue_selector::ResidueSelector {

public:
	static core::Size const DEFAULT_NUM_SHELLS;
	static core::Real const DEFAULT_BUMP_FACTOR;
	static bool const DEFAULT_INVERT;
	static bool const DEFAULT_INCLUDE_FOCUS;
	static bool const DEFAULT_FOCUS_ON_DESIGNABLE;

public:
	/// @brief Default constructor.
	ClashBasedShellSelector() = default;

	/// @brief Default copy constructor.
	ClashBasedShellSelector(ClashBasedShellSelector const &) = default;

	/// @brief Default destructor.
	~ClashBasedShellSelector() = default;

	/// Initialize from a string containing a list of indices, e.g.
	/// "34-46,199-202".
	ClashBasedShellSelector(
		std::string resnums);

	/// @brief Initialize from a residue selector.
	ClashBasedShellSelector(
		core::select::residue_selector::ResidueSelectorCOP selector);

	/// @brief Initialize from a task factory.
	ClashBasedShellSelector(
		core::pack::task::TaskFactoryCOP task_factory,
		bool focus_on_designable=true);

	/// @brief Initialize from a packer task.
	/// @details It's usually a bad idea to use a packer task instead of a task
	/// factory, because you'll get bugs if the pose used to create the packer
	/// task is different than the one passed to this selector's apply() method.
	ClashBasedShellSelector(
		core::pack::task::PackerTaskCOP task,
		bool focus_on_designable=true);

	/// @brief Initialize from a set of residue indices.
	ClashBasedShellSelector(
		std::set<core::Size> resnums);

	/// @brief Initialize from a boolean vector.
	/// @details The vector should have the same number of indices as the pose
	/// ultimately passed to apply().
	ClashBasedShellSelector(
		utility::vector1<bool> mask);

	/// @brief Copy this object and return an owning pointer to the result.
	virtual
	core::select::residue_selector::ResidueSelectorOP
	clone() const;

	/// @brief Select the residues that are clashing with any rotamer of the
	/// previously specified focus residues.
	/// @details Rotamers that are clashing with the backbone are ignored, since
	/// the idea is to identify residues that might need to be repacked, and
	/// backbone clashes can't be solved by repacking.
	virtual
	core::select::residue_selector::ResidueSubset
	apply( core::pose::Pose const & pose ) const;

	/// @brief Initialize from a RosettaScripts tag.
	virtual void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	);

	/// @brief Define the expected RosettaScripts options.
	static void
	provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd );

	/// @brief Define the attributes understood by this RosettaScripts tag.
	/// @details This is provided separately from provide_xml_schema() so other
	/// entities (namely ClashBasedRepackShell) can easily support all the same
	/// attributes.
	static void
	provide_xml_schema_attributes( utility::tag::AttributeList & attributes );

	/// @brief Return "ClashBasedShell".
	virtual std::string get_name() const;

	/// @brief Return "ClashBasedShell".
	static std::string class_name();

	/// @brief Select the residues that @b aren't part of the clash-based shell.
	/// @details This is useful for creating a repack shell, because you
	/// typically start with everything being repackable and your goal is to
	/// freeze the residues that aren't in a clash-based shell.
	void invert( bool inverted=true );

	/// @brief Reset the focus residues.
	void clear_focus();

	/// @brief Specify the residues to build the shell around from a string
	/// containing a list of indices, e.g. "34-46,199-202".
	/// @details The shell will be built considering all rotamers of the native
	/// amino acids of the specified residues.  Specify a TaskFactory if you want
	/// to account for design.
	void set_focus( std::string resnums );

	/// @brief Specify the residues to build the shell around from a residue
	/// selector.
	/// @details The shell will be built considering all rotamers of the native
	/// amino acids of the specified residues.  Specify a TaskFactory if you want
	/// to account for design.
	void set_focus( core::select::residue_selector::ResidueSelectorCOP selector );

	/// @brief Specify the residues to build the shell around from a task
	/// factory.
	///
	/// @details The shell will be built considering all rotamers allowed by the
	/// given task factory.  By default, only positions that can design are
	/// considered (because usually you want to build a repack around a set of
	/// designable positions), but you can also include repack positions by
	/// setting @a focus_on_designable to false.
	void set_focus(
		core::pack::task::TaskFactoryCOP task_factory,
		bool focus_on_designable=DEFAULT_FOCUS_ON_DESIGNABLE );

	/// @brief Specify the residues to build the shell around from a packer task.
	///
	/// @details It's usually a bad idea to use a packer task instead of a task
	/// factory, because you'll get bugs if the pose used to create the packer
	/// task is different than the one passed to this selector's apply() method.
	///
	/// The shell will be built considering all rotamers allowed by the given
	/// packer task.  By default, only positions that can design are considered
	/// (because usually you want to build a repack around a set of designable
	/// positions), but you can also include repack positions by setting @a
	/// focus_on_designable to false.
	void set_focus(
		core::pack::task::PackerTaskCOP task,
		bool focus_on_designable=DEFAULT_FOCUS_ON_DESIGNABLE );

	/// @brief Specify the residues to build the shell around from a set of
	/// residue indices.
	/// @details The shell will be built considering all rotamers of the native
	/// amino acids of the specified residues.  Specify a TaskFactory if you want
	/// to account for design.
	void set_focus( std::set<core::Size> resnums );

	/// @brief Specify the residues to build the shell around from a boolean
	/// vector.
	/// @details The vector should have the same number of indices as the pose
	/// ultimately passed to apply().  The shell will be built considering all
	/// rotamers of the native amino acids of the specified residues.  Specify a
	/// TaskFactory if you want to account for design.
	void set_focus( utility::vector1<bool> bool_mask );

	/// @brief Set the score function to use when generating rotamers.
	/// @details I'm a little skeptical that we couldn't just use a
	/// default-constructed score function for this. -KBK
	void set_scorefxn( core::scoring::ScoreFunctionOP scorefxn );

	/// @brief If true, include the focus residues in the ultimate selection.
	void set_include_focus( bool include_focus );

	/// @brief Specify the number of shells to calculate.
	/// @details For example, if you asked for two shells, the first shell would
	/// include all the residues that could clash with the focus, and the second
	/// shell would contain all the residues that could clash with anything in
	/// the first shell.  The ultimate selection would be the union of the two
	/// shells.
	void set_num_shells( core::Size num_shells );

	/// @brief Control how close two atoms can be before they are considered to
	/// be clashing.
	/// @details Two atoms are considered to clash when the squared distance
	/// between them is less than: bump_overlap_factor * {sum of Lennard-Jones
	/// radii}^2.
	void set_bump_overlap_factor( core::Real set_bump_overlap_factor );

	/// @brief Get the score function to use when generating rotamers.
	core::scoring::ScoreFunctionOP get_scorefxn() const;

	/// @brief Get the number of shells to calculate.
	core::Size get_num_shells() const;

	/// @brief Get how close two atoms can be before they are considered to be
	/// clashing.
	core::Real get_bump_overlap_factor() const;

	/// @brief Return whether the focus residues will be included in the ultimate
	/// selection.
	bool get_include_focus() const;

	/// @brief Return whether to selection will be inverted.
	bool is_inverted() const;

private:

	std::string focus_str_;
	std::set<core::Size> focus_set_;
	utility::vector1<bool> focus_mask_;
	core::select::residue_selector::ResidueSelectorCOP focus_selector_ = nullptr;
	core::pack::task::PackerTaskCOP focus_task_ = nullptr;
	core::pack::task::TaskFactoryCOP focus_task_factory_ = nullptr;
	core::scoring::ScoreFunctionOP scorefxn_ = nullptr;
	core::Size num_shells_ = DEFAULT_NUM_SHELLS;
	core::Real bump_overlap_factor_ = DEFAULT_BUMP_FACTOR;
	bool is_inverted_ = DEFAULT_INVERT;
	bool include_focus_ = DEFAULT_INCLUDE_FOCUS;
	bool focus_on_designable_ = DEFAULT_FOCUS_ON_DESIGNABLE;

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};


/// @brief Return the indices of any residues that clash with any of the
/// rotamers allowed by the given PackerTask.
///
/// See ClashBasedShellSelector for more details about how clash-based shell
/// selection works, although most of the logic is implemented in this
/// function.
std::set<core::Size>
find_clashing_shells(
	core::pack::task::PackerTaskCOP const focus,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::Size const num_shells=ClashBasedShellSelector::DEFAULT_NUM_SHELLS,
	bool const include_focus=ClashBasedShellSelector::DEFAULT_INCLUDE_FOCUS,
	core::Real const bump_factor=ClashBasedShellSelector::DEFAULT_BUMP_FACTOR);

/// @brief Find a single shell of clashing residues.
std::set<core::Size>
find_clashing_shell(
	core::pack::task::PackerTaskCOP const focus,
	core::pose::Pose const & pose,
	core::scoring::ScoreFunction const & scorefxn,
	core::Real const bump_factor=ClashBasedShellSelector::DEFAULT_BUMP_FACTOR);

/// @brief Add any residue that clashes with the given rotamer to the given
/// shell, unless the given rotamer also clashes with the backbone somewhere.
void
add_clashes_to_shell(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rot1,
	std::set<core::Size> const & focus,
	std::set<core::Size> & shell);

void
add_clashes_to_shell(
	core::pose::Pose const & pose,
	core::conformation::Residue const & rot1,
	core::Real const bump_factor,
	std::set<core::Size> const & focus,
	std::set<core::Size> & shell);

/// @brief Return true if there is a sidechain clash between the two given
/// residues.
bool
is_sc_sc_clash(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	core::Real bump_factor=ClashBasedShellSelector::DEFAULT_BUMP_FACTOR);

/// @brief Return true if a sidechain atom in the first residue clashes with a
/// backbone atom in the second residue.
bool
is_sc_bb_clash(
	core::conformation::Residue const & sc_rsd,
	core::conformation::Residue const & bb_rsd,
	core::Real bump_factor=ClashBasedShellSelector::DEFAULT_BUMP_FACTOR);


/// @brief Return a set containing the indices with true values in the given
/// boolean mask.
std::set<core::Size>
resnums_from_bool_mask(
	utility::vector1<bool> mask);

/// @brief Return a set containing the positions that are not fixed (i.e. that
/// are packable or designable) in the given PackerTask.
std::set<core::Size>
resnums_from_task(
	core::pack::task::PackerTaskCOP task);

/// @brief Return a boolean mask where the indices in the given set are true
/// and all others are false.
utility::vector1<bool>
bool_mask_from_resnums(
	core::pose::Pose const & pose,
	std::set<core::Size> const & resis);

/// @brief Return a boolean mask where the indices that are not fixed (i.e.
/// that are packable or designable) in the given PackerTask are true and all
/// others are false.
utility::vector1<bool>
bool_mask_from_task(
	core::pack::task::PackerTaskCOP task);

/// @brief Return a PackerTask where the indices that are true in the given
/// boolean mask are repackable and all others are frozen.
core::pack::task::PackerTaskOP
task_from_bool_mask(
	core::pose::Pose const & pose,
	utility::vector1<bool> mask);

/// @brief Return a PackerTask where the indices that in the given set are
/// repackable and all others are frozen.
core::pack::task::PackerTaskOP
task_from_resnums(
	core::pose::Pose const & pose,
	std::set<core::Size> resnums);


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( core_pack_task_residue_selector_ClashBasedShellSelector )
#endif // SERIALIZATION


#endif
