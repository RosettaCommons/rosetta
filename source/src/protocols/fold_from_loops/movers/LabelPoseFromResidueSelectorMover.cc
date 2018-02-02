// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/fold_from_loops/LabelPoseFromResidueSelectorMover.cc
/// @brief Adds labels to residues selected through a ResidueSelector
/// @author Jaume Bonet (jaume.bonet@gmail.com)

#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMover.hh>
#include <protocols/fold_from_loops/movers/LabelPoseFromResidueSelectorMoverCreator.hh>
#include <protocols/fold_from_loops/utils/utils.hh>

// Protocol headers
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/PDBInfo.hh>
#include <core/select/residue_selector/ResidueRanges.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/NotResidueSelector.hh>
#include <core/select/residue_selector/util.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/pose/selection.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/import_pose/import_pose.hh>

// Basic/Utility headers
#include <basic/datacache/DataMap.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <utility/string_util.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.fold_from_loops.LabelPoseFromResidueSelectorMover", basic::t_trace );

namespace protocols {
namespace fold_from_loops {
namespace movers {

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover():
	protocols::moves::Mover( mover_name() ),
	selector_( default_selector() ),
	label_( default_label() ),
	unlabel_( default_unlabel() ),
	from_remark_( default_from_remark() )
{}

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover( core::select::residue_selector::ResidueSelectorCOP selector ):
	protocols::moves::Mover( mover_name() ),
	selector_(std::move( selector )),
	label_( default_label() ),
	unlabel_( default_unlabel() ),
	from_remark_( default_from_remark() )
{}

LabelPoseFromResidueSelectorMover::LabelPoseFromResidueSelectorMover( core::select::residue_selector::ResidueSelectorCOP selector, std::string const & label ):
	protocols::moves::Mover( mover_name() ),
	selector_(std::move( selector )),
	label_( label ),
	unlabel_( default_unlabel() ),
	from_remark_( default_from_remark() )
{}

LabelPoseFromResidueSelectorMover::~LabelPoseFromResidueSelectorMover()= default;

void
LabelPoseFromResidueSelectorMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{

	if ( !tag->hasOption("from_remark") ) {
		if ( !tag->hasOption("property") and !tag->hasOption("remove_property") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "You need to provide a label to add (property) or remove (remove_property) or both.");
		}
		if ( !tag->hasOption("residue_selector") ) {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "A ResidueSelector must be provided to decide the residues to label");
		}
		core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::parse_residue_selector( tag, data );
		if ( selector ) residue_selector( selector );

		label( tag->getOption< std::string >( "property", default_label() ) );
		unlabel( tag->getOption< std::string >( "remove_property", default_label() ) );
	} else {
		from_remark( tag->getOption< bool >( "from_remark", default_from_remark() ) );
	}
}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMover::clone() const
{
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover( *this ) );
}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover );
}

void
LabelPoseFromResidueSelectorMover::apply_no_remarks( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;
	ResidueSubset const subset = selector_->apply( pose );
	if ( all_false_selection( subset ) ) {
		TR.Debug << "No residues meet the ResidueSelector criteria" << std::endl;
		return;
	}
	runtime_assert_msg( not label_.empty() or not unlabel_.empty(),
		"At least one label needs to be provided to insert/remove!" );

	if ( label_.size() > 0 ) {
		TR << "Add label " << label_ << std::endl;
	}
	if ( unlabel_.size() > 0 ) {
		TR << "Remove label " << unlabel_ << std::endl;
	}
	for ( core::Size i=1; i<=subset.size(); ++i ) {
		if ( subset[i] ) {
			if ( label_.size() > 0 ) {
				TR << "Adding label " << label_ << " to residue " << i << std::endl;
				pose.pdb_info()->add_reslabel(i, label_ );
			}
			if ( unlabel_.size() > 0 ) {
				if ( pose.pdb_info()->res_haslabel(i, unlabel_ ) ) {
					TR << "Removing label " << unlabel_ << " from residue " << i << std::endl;
					utility::vector1< std::string> labels = pose.pdb_info()->get_reslabels( i );
					pose.pdb_info()->clear_reslabel( i );
					for ( auto const & label : labels ) {
						if ( label != unlabel_ ) {
							pose.pdb_info()->add_reslabel(i, label );
						}
					}
				}
			}
		}
	}
}

void
LabelPoseFromResidueSelectorMover::apply_remarks( core::pose::Pose & pose )
{
	using namespace core::select::residue_selector;
	std::string alllabs;
	if ( core::pose::get_comment( pose, "LABELS", alllabs) ) {
		TR << alllabs << std::endl;
		utility::vector1< std::string > thelabs = utility::string_split( alllabs, ';');
		for ( auto lab : thelabs ) {
			TR << "  " << lab << std::endl;
			utility::vector1< std::string > thislab = utility::string_split( lab, ':');
			TR << "    " << thislab[1] << std::endl;
			TR << "    " << thislab[2] << std::endl;
			ResidueSelectorCOP selector( new ResidueIndexSelector( thislab[2] ) );
			LabelPoseFromResidueSelectorMoverOP labeler( new LabelPoseFromResidueSelectorMover );
			labeler->residue_selector( selector );
			labeler->label( thislab[1] );
			labeler->apply( pose );
		}
	} else {
		TR << "No remarks found" << std::endl;
	}
}

void
LabelPoseFromResidueSelectorMover::apply( core::pose::Pose & pose )
{
	if ( from_remark_ ) {
		TR << "Apply labels from SilentFile Remarks" << std::endl;
		apply_remarks( pose );
	} else {
		apply_no_remarks( pose );
	}
}

void
LabelPoseFromResidueSelectorMover::label( std::string const & label )
{
	label_ = label;
}

void
LabelPoseFromResidueSelectorMover::unlabel( std::string const & label )
{
	unlabel_ = label;
}

void
LabelPoseFromResidueSelectorMover::residue_selector( core::select::residue_selector::ResidueSelectorCOP const & selector )
{
	selector_ = selector;
}

void
LabelPoseFromResidueSelectorMover::residue_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	selector_ = selector.clone();
}
void
LabelPoseFromResidueSelectorMover::from_remark( bool const & pick )
{
	from_remark_ = pick;
}

std::string LabelPoseFromResidueSelectorMover::get_name() const {
	return mover_name();
}

std::string LabelPoseFromResidueSelectorMover::mover_name() {
	return "LabelPoseFromResidueSelectorMover";
}

void LabelPoseFromResidueSelectorMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "property", xs_string, "New label to add to the residues." )
		+ XMLSchemaAttribute( "remove_property", xs_string, "Label to remove from the selected residues (if present)." )
		+ XMLSchemaAttribute::attribute_w_default( "from_remark", xsct_rosetta_bool,
		"Read the LABELS from REMARKs in the input structure or silent file (when called, ignores all other attributes, "
		"so actively setting it to false makes the mover do nothing).", std::to_string( default_from_remark() ) );
	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "Selector specifying residues to be labeled." );
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Add/Remove labels to/from residues defined by ResidueSelector.", attlist );
}

core::select::residue_selector::ResidueSelectorCOP
LabelPoseFromResidueSelectorMover::default_selector()
{
	using namespace core::select::residue_selector;
	ResidueSelectorCOP true_selector( new TrueResidueSelector );
	ResidueSelectorCOP false_selector( new NotResidueSelector( true_selector ) );
	return false_selector;
}

std::string
LabelPoseFromResidueSelectorMover::default_label()
{
	return std::string();
}

std::string
LabelPoseFromResidueSelectorMover::default_unlabel()
{
	return std::string();
}

bool
LabelPoseFromResidueSelectorMover::default_from_remark()
{
	return false;
}


std::string LabelPoseFromResidueSelectorMoverCreator::keyname() const {
	return LabelPoseFromResidueSelectorMover::mover_name();
}

protocols::moves::MoverOP
LabelPoseFromResidueSelectorMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new LabelPoseFromResidueSelectorMover );
}

void LabelPoseFromResidueSelectorMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	LabelPoseFromResidueSelectorMover::provide_xml_schema( xsd );
}

}
} //protocols
} //fold_from_loops
