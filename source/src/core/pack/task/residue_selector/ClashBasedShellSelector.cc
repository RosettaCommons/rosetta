// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file     core/pack/task/residue_selector/ClashBasedShellSelector.cc
/// @brief    The ClashBasedShellSelector identifies all residues that clash with at least one rotamer of a design position
/// @detailed The ClashBasedShellSelector identifies all residues that clash with at least one rotamer of a design position. There are two recommended ways to use this selector.
/// @detailed (Usage 1) identify all residues that could clash with residues designated as design positions by the packer. The constructor for this usage takes (a) task, (b) scorefxn, and (c) pose.
/// @detailed (Usage 2) identify all residues that could clash with an arbitrary set of residues. The constructor for this usage takes one argument of type std::set<core::Size>.
/// @detailed In either usage, once you create the residue selector object, you then get result of which positions could clash using the apply function e.g.  utility::vector1< bool > clashing_residues = ClashBasedShellSelector_object->apply( *pose_copy );
/// @author   Noah Ollikainen (nollikai@gmail.com)
/// @author   Roland A. Pache, PhD
/// @author   Amanda Loshbaugh (aloshbau@gmail.com)
/// @author   Kale Kundert (kale@thekunderts.net)

// Unit headers
#include <core/pack/task/residue_selector/ClashBasedShellSelector.hh>
#include <core/pack/task/residue_selector/ClashBasedShellSelectorCreator.hh>

// Core headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/find_neighbors.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/PointGraph.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/task/xml_util.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pose/selection.hh>
#include <core/pose/util.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/xml_util.hh>
#include <core/scoring/EnergyGraph.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/util.hh>
#include <core/select/residue_selector/ResidueSelectorCreators.hh>
#include <core/types.hh>

// Basic headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/exit.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector1.hh>

// C++ headers
#include <utility/assert.hh>
#include <boost/format.hpp>

// Serialization headers
#ifdef SERIALIZATION
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>
#include <cereal/types/set.hpp>
#include <cereal/types/polymorphic.hpp>
#endif

namespace core {
namespace pack {
namespace task {
namespace residue_selector {

using namespace std;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::ScoreFunctionOP;
using core::select::residue_selector::ResidueSubset;
using core::select::residue_selector::ResidueSelectorOP;
using core::select::residue_selector::ResidueSelectorCOP;
using core::select::residue_selector::xsd_type_definition_w_attributes_and_optional_subselector;
using utility::vector1;

static basic::Tracer TR( "core.pack.task.residue_selector.ClashBasedShellSelector" );
Size const ClashBasedShellSelector::DEFAULT_NUM_SHELLS = 1;
Real const ClashBasedShellSelector::DEFAULT_BUMP_FACTOR = 0.5;
bool const ClashBasedShellSelector::DEFAULT_INVERT = false;
bool const ClashBasedShellSelector::DEFAULT_INCLUDE_FOCUS = true;
bool const ClashBasedShellSelector::DEFAULT_FOCUS_ON_DESIGNABLE = true;

set<Size> find_clashing_shells(
	PackerTaskCOP const focus_task,
	Pose const & pose,
	ScoreFunction const & scorefxn,
	Size const num_shells,
	bool const include_focus,
	Real const bump_factor) {

	set<Size> combined_shells;
	PackerTaskOP task = focus_task->clone();

	if ( include_focus ) {
		combined_shells = resnums_from_task(focus_task);
	}

	for ( Size i = 1; i <= num_shells; i++ ) {
		set<Size> shell = find_clashing_shell(task, pose, scorefxn, bump_factor);

		// Keep track of all the clashing residues we found.
		for ( Size j: shell ) combined_shells.insert(j);

		// Update the packer task so that on the next iteration, we'll look for the
		// residues that are clashing with the shell we just found.
		task = task_from_resnums(pose, shell);
	}

	return combined_shells;
}

set<Size> find_clashing_shell(
	PackerTaskCOP task,
	Pose const & const_pose,
	ScoreFunction const & scorefxn,
	Real const bump_factor) {

	using core::pack::rotamer_set::RotamerSetFactory;
	using core::pack::rotamer_set::RotamerSetOP;
	using utility::graph::GraphOP;

	Pose pose(const_pose);
	set<Size> shell, focus = resnums_from_task(task);

	for ( Size i : focus ) {
		// Create a new packer task that will freeze everything except residue #i.
		PackerTaskOP bump_check_task = task->clone();
		bump_check_task->set_bump_check(false);                  // Don't throw out rotamers during packer graph generation
		bump_check_task->temporarily_fix_everything();           // NATRO for all residues...
		bump_check_task->temporarily_set_pack_residue(i, true);  // ...except residue #i.
		scorefxn.setup_for_packing(
			pose,
			bump_check_task->repacking_residues(),
			bump_check_task->designing_residues());

		// Create the rotamer set for the given design residue.
		// Residue const & res = pose.residue(i);
		RotamerSetOP rotamer_set = RotamerSetFactory::create_rotamer_set( pose );
		rotamer_set->set_resid( i );
		GraphOP packer_graph = core::pack::create_packer_graph(
			pose, scorefxn, bump_check_task);

		rotamer_set->set_resid(i);
		rotamer_set->build_rotamers(pose, scorefxn, *bump_check_task, packer_graph);

		// Find positions where any rotamer of the residue in question is clashing
		// with the sidechain, but not the backbone.
		for ( Size j = 1; j <= rotamer_set->num_rotamers(); j++ ) {
			add_clashes_to_shell(
				pose, *rotamer_set->rotamer(j), bump_factor, focus, shell);
		}
	}

	// Let the user know which residues were selected.
	if ( focus.empty() ) {
		TR.Warning << "No residues were given, the shell will be empty!" << endl;
	}
	TR << "Building a shell around the following residues: " << focus << endl;
	TR << "Adding the following residues to the shell: " << shell << endl;

	return shell;
}

void add_clashes_to_shell(
	Pose const & pose,
	Residue const & rot1,
	set<Size> const & focus,
	set<Size> & shell) {

	add_clashes_to_shell(pose, rot1, ClashBasedShellSelector::DEFAULT_BUMP_FACTOR, focus, shell);
}

void add_clashes_to_shell(
	Pose const & pose,
	Residue const & rot1,
	Real const bump_factor,
	set<Size> const & focus,
	set<Size> & shell) {

	bool checked_for_bb_clashes = false;

	for ( Size i = 1; i <= pose.size(); ++i ) {
		// Skip the bump check if we already know about the residue in question.
		if ( i == rot1.seqpos() ) continue;
		if ( focus.find( i ) != focus.end() ) continue;
		if ( shell.find( i ) != shell.end() ) continue;

		// Check for clashes between all pairs of side-chain atoms between the
		// query rotamer and the current neighbor.
		if ( is_sc_sc_clash(rot1, pose.residue(i), bump_factor) ) {

			TR.Trace.visible() && TR.Trace <<
				boost::format("Found sidechain clash: %s%s with %s%s rotamer")
				% pose.residue(i).name1() % i
				% rot1.name1() % rot1.seqpos() << endl;

			// If we haven't checked for backbone clashes, do that.  We don't want
			// any rotamers that are clashing with the backbone, because repacking
			// can't fix that.  Making this check after the sidechain check saves a
			// lot of time in the case that a significant fraction of the protein is
			// either in the focus or the shell already.
			if ( ! checked_for_bb_clashes ) {
				for ( Size j = 1; j <= pose.size(); ++j ) {
					if ( j == rot1.seqpos() ) continue;
					if ( is_sc_bb_clash(rot1, pose.residue(j), bump_factor) ) {
						TR.Trace.visible() && TR.Trace <<
							boost::format("  Skipping due to backbone clash: %s%s with above rotamer")
							% pose.residue(j).name1() % j << endl;
						return;
					}
				}

				checked_for_bb_clashes = true;
			}

			TR.Trace.visible() && TR.Trace <<
				boost::format("  Adding %s%s to the shell.")
				% pose.residue(i).name1() % i << endl;

			shell.insert(i);
		}
	}
}

bool is_sc_sc_clash(
	Residue const & rsd1,
	Residue const & rsd2,
	Real bump_factor) {

	for ( Size m = rsd1.first_sidechain_atom(); m <= rsd1.nheavyatoms(); ++m ) {
		core::Vector const & atom1_coords( rsd1.xyz( m ) );

		for ( Size n = rsd2.first_sidechain_atom(); n <= rsd2.nheavyatoms(); ++n ) {
			core::Vector const & atom2_coords( rsd2.xyz( n ) );

			// Use the sum of the LJ radii as the distance cutoff with a bump overlap factor
			Real lj_sum = (rsd1.atom_type( m ).lj_radius() + rsd2.atom_type( n ).lj_radius());
			Real distance_cutoff_squared = (lj_sum * lj_sum) * bump_factor;
			Real distance_squared = ( atom1_coords - atom2_coords ).length_squared();

			if ( distance_squared < distance_cutoff_squared ) return true;
		}
	}
	return false;
}

bool is_sc_bb_clash(
	Residue const & sc_rsd,
	Residue const & bb_rsd,
	Real bump_factor) {

	const core::Size first_backbone_atom = 1;

	for ( core::Size m = sc_rsd.first_sidechain_atom(); m <= sc_rsd.nheavyatoms(); ++m ) {
		core::Vector const & atom1_coords( sc_rsd.xyz( m ) );

		for ( core::Size n = first_backbone_atom; n <= bb_rsd.last_backbone_atom(); ++n ) {
			core::Vector const & atom2_coords( bb_rsd.xyz( n ) );

			// Use the sum of the LJ radii as the distance cutoff with a bump overlap factor
			core::Real lj_sum = (sc_rsd.atom_type( m ).lj_radius() + bb_rsd.atom_type( n ).lj_radius());
			core::Real distance_cutoff_squared = (lj_sum * lj_sum) * bump_factor;
			core::Real distance_squared = ( atom1_coords - atom2_coords ).length_squared();

			if ( distance_squared < distance_cutoff_squared ) return true;
		}
	}
	return false;
}

set<Size> resnums_from_bool_mask(
	vector1<bool> mask) {

	set<Size> resnums;

	for ( Size i = 1; i <= mask.size(); i++ ) {
		if ( mask[i] ) {
			resnums.insert(i);
		}
	}

	return resnums;
}

set<Size> resnums_from_task(
	PackerTaskCOP task) {

	return resnums_from_bool_mask(bool_mask_from_task(task));
}

vector1<bool> bool_mask_from_resnums(
	Pose const & pose,
	set<Size> const & resnums) {

	vector1<bool> mask( pose.size(), false );

	for ( Size i : resnums ) {
		runtime_assert( i > 0 && i <= pose.size() );
		mask[ i ] = true;
	}

	return mask;
}

vector1<bool> bool_mask_from_task(
	PackerTaskCOP task) {

	return task->repacking_residues();
}

PackerTaskOP task_from_resnums(
	Pose const & pose,
	set<Size> resnums) {

	return task_from_bool_mask(pose, bool_mask_from_resnums(pose, resnums));
}

PackerTaskOP task_from_bool_mask(
	Pose const & pose,
	vector1<bool> mask) {

	using namespace core::pack::task::operation;
	runtime_assert( pose.size() == mask.size() );

	// Restrict every position to repacking (no design).
	TaskFactory factory;
	factory.push_back(TaskOperationOP(new RestrictToRepacking));
	PackerTaskOP task = factory.create_task_and_apply_taskoperations(pose);

	// Prevent the positions that are outside the mask from moving at all.
	task->restrict_to_residues(mask);

	// Return a task where the positions in the mask can repack and everything
	// else is frozen.
	return task;
}

void freeze_repack_positions(PackerTaskOP task) {
	task->restrict_to_residues(task->designing_residues());
}


ClashBasedShellSelector::ClashBasedShellSelector(
	string resnums) {

	set_focus(resnums);
}

ClashBasedShellSelector::ClashBasedShellSelector(
	ResidueSelectorCOP selector) {

	set_focus(selector);
}

ClashBasedShellSelector::ClashBasedShellSelector(
	TaskFactoryCOP task_factory,
	bool focus_on_designable) {

	set_focus(task_factory, focus_on_designable);
}

ClashBasedShellSelector::ClashBasedShellSelector(
	PackerTaskCOP task,
	bool focus_on_designable) {

	set_focus(task, focus_on_designable);
}

ClashBasedShellSelector::ClashBasedShellSelector(
	set<Size> resnums) {

	set_focus(resnums);
}

ClashBasedShellSelector::ClashBasedShellSelector(
	vector1<bool> mask) {

	set_focus(mask);
}

/// @brief Clone operator.
/// @details Copy this object and return an owning pointer to the new object.
ResidueSelectorOP ClashBasedShellSelector::clone() const {
	return ResidueSelectorOP( new ClashBasedShellSelector(*this) );
}

ResidueSubset
ClashBasedShellSelector::apply( Pose const & pose ) const {

	// Work out which residues to build the repack shell around.
	//
	// We want to support lots of different ways of specifying the residues to
	// build the shell(s) around, to make this class nice and user-friendly.  The
	// purpose of this code is to get the information for any of these different
	// sources into a consistent format, specifically a PackerTask.  The setter
	// methods will ensure that only one of the "focus" variables can have a
	// value at a time (the ones that don't will be null pointers or empty
	// containers).

	PackerTaskOP task;

	if ( ! focus_str_.empty() ) {
		set<Size> resnums = core::pose::get_resnum_list(focus_str_, pose);
		task = task_from_resnums(pose, resnums);
	} else if ( focus_selector_ ) {
		task = task_from_bool_mask(pose, focus_selector_->apply(pose));
	} else if ( focus_task_factory_ ) {
		task = focus_task_factory_->create_task_and_apply_taskoperations(pose);
		if ( focus_on_designable_ ) freeze_repack_positions(task);
	} else if ( focus_task_ ) {
		if ( focus_task_->total_residue() != pose.size() ) {
			auto message = boost::format("PackerTask (%d res) out of sync with pose (%d res).") % focus_task_->total_residue() % pose.size();
			utility_exit_with_message(message.str());
		}
		task = focus_task_->clone();
		if ( focus_on_designable_ ) freeze_repack_positions(task);
	} else if ( ! focus_mask_.empty() ) {
		if ( focus_mask_.size() != pose.size() ) {
			auto message = boost::format("Boolean mask (%d res) out of sync with pose (%d res).") % focus_mask_.size() % pose.size();
			utility_exit_with_message(message.str());
		}
		task = task_from_bool_mask(pose, focus_mask_);
	} else {
		// Have this be the default case, because an empty set is a sane default.
		task = task_from_resnums(pose, focus_set_);
	}

	// Make a score function if we weren't given one.

	ScoreFunctionOP scorefxn = scorefxn_;

	if ( ! scorefxn ) {
		scorefxn = core::scoring::get_score_function();
	}

	// Create and return the clash-based repack shell.

	set<Size> shells = find_clashing_shells(
		task, pose, *scorefxn,
		num_shells_, include_focus_, bump_overlap_factor_);

	ResidueSubset selection = bool_mask_from_resnums(pose, shells);
	if ( is_inverted_ ) selection.flip();

	return selection;
}

void ClashBasedShellSelector::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap) {

	using core::pack::task::parse_task_operations;
	using core::select::residue_selector::parse_residue_selector;
	using core::select::residue_selector::get_embedded_residue_selector;
	using core::scoring::parse_score_function;

	int set_focus_counter = 0;

	if ( tag->hasOption("resnums") ) {
		set_focus( tag->getOption<string>("resnums") );
		set_focus_counter += 1;
	}
	if ( tag->hasOption("residue_selector") ) {
		set_focus( parse_residue_selector(tag, datamap) );
		set_focus_counter += 1;
	}
	if ( tag->hasOption("task_operations") ) {
		set_focus(
			parse_task_operations(tag, datamap),
			tag->getOption<bool>("focus_on_designable", true) );
		set_focus_counter += 1;
	}
	if ( tag->size() > 1 ) {
		set_focus( get_embedded_residue_selector(tag, datamap) );
		set_focus_counter += 1;
	}

	// Make sure nothing ends up getting silently ignored.
	if ( set_focus_counter > 1 ) {
		TR.Error << "May not specify more than one of 'resnums', 'residue_selector', 'task_operations', or a residue selector subtag." << endl;
		utility_exit();
	}

	set_scorefxn( parse_score_function(tag, datamap) );
	set_num_shells( tag->getOption<Size>("num_shells", DEFAULT_NUM_SHELLS) );
	set_bump_overlap_factor( tag->getOption<Real>("bump_overlap_factor", DEFAULT_BUMP_FACTOR) );
	set_include_focus( tag->getOption<bool>("include_focus", DEFAULT_INCLUDE_FOCUS) );
	invert( tag->getOption<bool>("invert", DEFAULT_INVERT) );
}

void
ClashBasedShellSelector::provide_xml_schema(
	utility::tag::XMLSchemaDefinition & xsd ) {

	utility::tag::AttributeList attributes;
	provide_xml_schema_attributes(attributes);

	xsd_type_definition_w_attributes_and_optional_subselector(
		xsd, class_name(),
		"The ClashBasedShellSelector identifies all residues that clash with at least one rotamer of a given residue selection.",
		attributes );
}

void
ClashBasedShellSelector::provide_xml_schema_attributes(
	utility::tag::AttributeList & attributes ) {

	using namespace utility::tag;
	using core::pack::task::attributes_for_parse_task_operations;
	using core::select::residue_selector::attributes_for_parse_residue_selector;
	using core::scoring::attributes_for_parse_score_function;

	attributes_for_parse_score_function(attributes);
	attributes_for_parse_task_operations(attributes);
	attributes_for_parse_residue_selector(attributes);

	attributes
		+ XMLSchemaAttribute::attribute_w_default(
		"resnums", xsct_int_cslist,
		"The residues to build the shell around.  These can also be specified via the 'selector' or 'task_operations' options, or via a ResidueSelector subtag.",
		"")
		+ XMLSchemaAttribute::attribute_w_default(
		"invert", xsct_rosetta_bool,
		"Select the residues that @b aren't part of the clash-based shell.",
		"false")
		+ XMLSchemaAttribute::attribute_w_default(
		"num_shells", xsct_non_negative_integer,
		"The number of shells to calculate.",
		"1")
		+ XMLSchemaAttribute::attribute_w_default(
		"include_focus", xsct_rosetta_bool,
		"If true, include the focus residues in the ultimate selection.",
		"true")
		+ XMLSchemaAttribute::attribute_w_default(
		"focus_on_designable", xsct_rosetta_bool,
		"Only applies if 'task_operations' is given.  By default, only residues that can design are considered when building the shell.  Set to 'false' to also consider residues that can repack.",
		"true")
		+ XMLSchemaAttribute::attribute_w_default(
		"bump_overlap_factor", xsct_real,
		"Two atoms are considered to clash when the squared distance between them is less than: bump_overlap_factor * {sum of Lennard-Jones radii}^2",
		"0.5");
}

std::string ClashBasedShellSelector::get_name() const {
	return ClashBasedShellSelector::class_name();
}

std::string ClashBasedShellSelector::class_name() {
	return "ClashBasedShell";
}

void ClashBasedShellSelector::invert( bool invert ) {
	is_inverted_ = invert;
}

void ClashBasedShellSelector::clear_focus() {
	focus_str_.clear();
	focus_set_.clear();
	focus_mask_.clear();
	focus_selector_ = nullptr;
	focus_task_ = nullptr;
	focus_task_factory_ = nullptr;
	focus_on_designable_ = true;
}

void ClashBasedShellSelector::set_focus( string resnum_list ) {
	clear_focus();
	focus_str_ = resnum_list;
}

void ClashBasedShellSelector::set_focus( ResidueSelectorCOP selector ) {
	clear_focus();
	focus_selector_ = selector;
}

void ClashBasedShellSelector::set_focus(
	TaskFactoryCOP task_factory,
	bool focus_on_designable) {

	clear_focus();
	focus_task_factory_ = task_factory;
	focus_on_designable_ = focus_on_designable;
}

void ClashBasedShellSelector::set_focus(
	PackerTaskCOP task,
	bool focus_on_designable) {

	clear_focus();
	focus_task_ = task;
	focus_on_designable_ = focus_on_designable;
}

void ClashBasedShellSelector::set_focus( set<Size> resnums ) {
	clear_focus();
	focus_set_ = resnums;
}

void ClashBasedShellSelector::set_focus( vector1<bool> bool_mask ) {
	clear_focus();
	focus_mask_ = bool_mask;
}

void ClashBasedShellSelector::set_scorefxn( ScoreFunctionOP scorefxn ) {
	scorefxn_ = scorefxn;
}

void ClashBasedShellSelector::set_num_shells( Size num_shells ) {
	num_shells_ = num_shells;
}

void ClashBasedShellSelector::set_bump_overlap_factor( Real bump_factor ) {
	bump_overlap_factor_ = bump_factor;
}

void ClashBasedShellSelector::set_include_focus( bool include_focus ) {
	include_focus_ = include_focus;
}

ScoreFunctionOP ClashBasedShellSelector::get_scorefxn() const {
	return scorefxn_;
}

Size ClashBasedShellSelector::get_num_shells() const {
	return num_shells_;
}

Real ClashBasedShellSelector::get_bump_overlap_factor() const {
	return bump_overlap_factor_;
}

bool ClashBasedShellSelector::get_include_focus() const {
	return include_focus_;
}

bool ClashBasedShellSelector::is_inverted() const {
	return is_inverted_;
}


core::select::residue_selector::ResidueSelectorOP
ClashBasedShellSelectorCreator::create_residue_selector() const {
	return core::select::residue_selector::ResidueSelectorOP( new ClashBasedShellSelector );
}

std::string
ClashBasedShellSelectorCreator::keyname() const {
	return ClashBasedShellSelector::class_name();
}

void
ClashBasedShellSelectorCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	ClashBasedShellSelector::provide_xml_schema( xsd );
}


} //namespace residue_selector
} //namespace task
} //namespace pack
} //namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pack::task::residue_selector::ClashBasedShellSelector::save( Archive & arc ) const {
	// TaskFactory can't be serialized, because doing so would require that all
	// TaskOperations be serializable, and that's too much of a burden for
	// supporting just this one case.  So instead we'll error out if we have a
	// TaskFactory, and simply ignore it otherwise.
	if ( focus_task_factory_ ) {
		utility_exit_with_message( "ClashBasedShellSelector: have a task factory, but task factories cannot be serialized." );
	}

	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( CEREAL_NVP( focus_str_ ) ); // std::string
	arc( CEREAL_NVP( focus_set_ ) ); // std::set<core::Size>
	arc( CEREAL_NVP( focus_mask_ ) ); // utility::vector1<bool>
	arc( CEREAL_NVP( focus_selector_ ) ); // core::select::residue_selector::ResidueSelectorCOP
	arc( CEREAL_NVP( focus_task_ ) ); // core::pack::task::PackerTaskCOP
	// EXEMPT focus_task_factory_
	arc( CEREAL_NVP( scorefxn_ ) ); // core::scoring::ScoreFunctionOP
	arc( CEREAL_NVP( num_shells_ ) ); // core::Size
	arc( CEREAL_NVP( bump_overlap_factor_ ) ); // core::Real
	arc( CEREAL_NVP( is_inverted_ ) ); // bool
	arc( CEREAL_NVP( include_focus_ ) ); // bool
	arc( CEREAL_NVP( focus_on_designable_ ) ); // bool
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pack::task::residue_selector::ClashBasedShellSelector::load( Archive & arc ) {
	arc( cereal::base_class< core::select::residue_selector::ResidueSelector >( this ) );
	arc( focus_str_ ); // std::string
	arc( focus_set_ ); // std::set<core::Size>
	arc( focus_mask_ ); // utility::vector1<bool>
	arc( focus_selector_ ); // core::select::residue_selector::ResidueSelectorCOP
	arc( focus_task_ ); // core::pack::task::PackerTaskCOP
	// EXEMPT focus_task_factory_
	arc( scorefxn_ ); // core::scoring::ScoreFunctionOP
	arc( num_shells_ ); // core::Size
	arc( bump_overlap_factor_ ); // core::Real
	arc( is_inverted_ ); // bool
	arc( include_focus_ ); // bool
	arc( focus_on_designable_ ); // bool
}

SAVE_AND_LOAD_SERIALIZABLE( core::pack::task::residue_selector::ClashBasedShellSelector );
CEREAL_REGISTER_TYPE( core::pack::task::residue_selector::ClashBasedShellSelector )
CEREAL_REGISTER_DYNAMIC_INIT( core_pack_task_residue_selector_ClashBasedShellSelector )

#endif // SERIALIZATION
