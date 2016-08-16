// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/ExtendedPoseBuilder.cc
/// @brief Builds a pose using a blueprint from a structure architect
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit header
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>

// Protocol headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers
#include <iterator>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.ExtendedPoseBuilder" );

namespace protocols {
namespace denovo_design {
namespace components {

ExtendedPoseBuilder::ExtendedPoseBuilder():
	utility::pointer::ReferenceCount()
{
}

ExtendedPoseBuilder::~ExtendedPoseBuilder()
{
}

core::pose::PoseOP
ExtendedPoseBuilder::apply( StructureData const & sd ) const
{
	SegmentNames const template_segments = find_template_segments( sd );
	TR.Debug << "Found template segments = " << template_segments << std::endl;

	core::pose::PoseOP newpose = create_template_pose( sd, template_segments );

	// rebuild fold tree
	if ( !template_segments.empty() ) {
		StructureDataFactory::get_instance()->clear_from_pose( *newpose );
		StructureData cur_sd = *StructureDataFactory::get_instance()->create_from_pose( *newpose );
		FoldGraph fg( cur_sd );
		newpose->fold_tree( fg.fold_tree( boost::assign::list_of (*cur_sd.segments_begin()) ) );
	}

	newpose->dump_pdb( "frag_templates.pdb" );
	extend_pose( *newpose, sd, template_segments );
	newpose->dump_pdb( "frag_built.pdb" );

	// rebuild fold tree
	FoldGraph fg( sd );
	newpose->fold_tree( fg.fold_tree( boost::assign::list_of(*sd.segments_begin()) ) );

	StructureDataFactory::get_instance()->save_into_pose( *newpose, sd );
	return newpose;
}

SegmentNames
ExtendedPoseBuilder::find_template_segments( StructureData const & sd ) const
{
	SegmentNames with_template;
	for ( SegmentNameList::const_iterator c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		if ( sd.segment( *c ).template_pose() ) with_template.push_back( *c );
	}
	return with_template;
}

core::pose::PoseOP
ExtendedPoseBuilder::create_template_pose( StructureData const & sd, SegmentNames const & template_segments ) const
{
	core::pose::PoseOP newpose( new core::pose::Pose );

	// look for segments with template pose
	// add these segments to the pose
	SegmentNames::const_iterator prev = template_segments.end();
	for ( SegmentNames::const_iterator s=template_segments.begin(); s!=template_segments.end(); ++s ) {
		Segment const & seg = sd.segment( *s );

		bool const new_chain = ( prev == template_segments.end() ) || ( *prev != seg.lower_segment() );
		if ( new_chain ) {
			TR.Debug << "Adding new chain from " << *s << std::endl;
			append_new_chain_from_template_segment( *newpose, seg );
		} else {
			TR.Debug << "Appending residues from " << *s << std::endl;
			append_residues_from_template_segment( *newpose, sd.segment( *prev ), seg );
		}
		prev = s;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Newpose: " << std::endl;
		for ( core::Size r=1; r<=newpose->total_residue(); ++r ) {
			TR.Debug << newpose->residue(r).name() << " " << r << std::endl;
		}
	}
	return newpose;
}

void
add_cutpoints( core::pose::Pose & pose, StructureData const & sd )
{
	for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		core::Size const cut = sd.segment( *s ).cutpoint();
		if ( cut ) {
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, cut );
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, cut+1 );
		}
	}
}

void
add_terminus_variants( core::pose::Pose & pose, ExtendedPoseBuilder::Resids const & endings )
{
	core::pose::add_lower_terminus_type_to_pose_residue( pose, 1 );
	for ( ExtendedPoseBuilder::Resids::const_iterator r=endings.begin(); r!=endings.end(); ++r ) {
		core::pose::add_upper_terminus_type_to_pose_residue( pose, *r );
		if ( *r + 1 <= pose.total_residue() ) {
			core::pose::add_lower_terminus_type_to_pose_residue( pose, *r+1 );
		}
	}
}

void
ExtendedPoseBuilder::extend_pose(
	core::pose::Pose & pose,
	StructureData const & sd,
	SegmentNames const & template_segments ) const
{
	Resids const endings = calc_chain_endings( sd );
	TR.Debug << "Chain endings = " << endings << std::endl;

	SegmentNames::const_iterator prev_template = template_segments.end();
	SegmentNames::const_iterator next_template = template_segments.begin();
	core::Size cur_resid = 1;
	core::Size extend_size = 0;
	for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		// if this is a template segment, prepend residues and skip ahead
		if ( ( next_template != template_segments.end() ) &&
				( *next_template == *s ) ) {
			TR.Debug << "Template reached: " << *next_template << " extending by " << extend_size << std::endl;
			core::Size const anchor = cur_resid + sd.segment( *s ).elem_length() - 1;
			prepend_new_residues( pose, extend_size, cur_resid, anchor, "VAL", sd.segment( *s ).lower_dihedrals() );
			if ( ( prev_template != template_segments.end() ) && ( !sd.segment( *prev_template ).has_free_upper_terminus() ) ) {
				TR.Debug << "Removing terminal variants starting at " << cur_resid - 1 << std::endl;
				core::pose::remove_upper_terminus_type_from_pose_residue( pose, cur_resid - 1 );
				core::pose::remove_lower_terminus_type_from_pose_residue( pose, cur_resid );
			}
			cur_resid += sd.segment( *s ).length();
			extend_size = 0;
			prev_template = next_template;
			++next_template;
			continue;
		}

		// if this segment contains a cutpoint, append residues up to the cut
		// and then start a new extension segment
		core::Size const cut = sd.segment( *s ).cutpoint();
		if ( cut != 0 ) {
			TR.Debug << "Cutpoint reached at residue " << cut << " in segment " << *s << std::endl;
			extend_size += sd.segment( *s ).n_residues_before_cutpoint();
			if ( prev_template != template_segments.end() ) {
				append_new_residues( pose,
					extend_size,
					cur_resid - 1,
					cur_resid - sd.segment( *prev_template ).elem_length(),
					"VAL",
					sd.segment( *prev_template ).upper_dihedrals() );
			} else {
				append_new_residues( pose, extend_size, cur_resid - 1, cur_resid - 2, "VAL", ResidueDihedrals() );
			}
			cur_resid += sd.segment( *s ).n_residues_before_cutpoint();
			extend_size = sd.segment( *s ).n_residues_after_cutpoint();
			continue;
		}

		TR.Debug << "Extending by segment " << *s << std::endl;
		extend_size += sd.segment( *s ).length();
	}
	if ( extend_size ) {
		debug_assert( next_template == template_segments.end() );
		if ( prev_template != template_segments.end() ) {
			append_new_residues( pose,
				extend_size,
				pose.total_residue(),
				pose.total_residue()-sd.segment( *prev_template ).elem_length()+1,
				"VAL",
				sd.segment( *prev_template ).upper_dihedrals() );
		} else {
			append_new_residues( pose,
				extend_size,
				pose.total_residue(),
				pose.total_residue() - 1,
				"VAL",
				ResidueDihedrals() );
		}
	}

	add_cutpoints( pose, sd );
	add_terminus_variants( pose, endings );
	pose.conformation().chains_from_termini();
}

ExtendedPoseBuilder::Resids
ExtendedPoseBuilder::calc_chain_endings( StructureData const & sd ) const
{
	Resids endings;
	for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		if ( ! sd.segment( *s ).cterm_included() ) {
			endings.push_back( sd.segment( *s ).upper() );
		}
	}
	return endings;
}

/// @brief helper function to create a ValineOP in the same typeset as the pose
core::conformation::ResidueOP
rsd_op( core::pose::Pose const & pose, std::string const & res_type )
{
	if ( pose.total_residue() == 0 ) {
		return core::conformation::ResidueFactory::create_residue(
			core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" )->name_map( res_type )
		);
	}
	core::chemical::ResidueType const & rtype = pose.residue(1).residue_type_set()->name_map( res_type );
	return core::conformation::ResidueFactory::create_residue( rtype );
}

void
append_new_chain_from_template_segment(
	core::pose::Pose & pose,
	Segment const & segment )
{
	debug_assert( segment.template_pose() );
	core::Size const chain_start = pose.total_residue() + 1;

	TR.Debug << "Appending " << segment.template_pose()->total_residue() << " template residues as new chain" << std::endl;
	if ( pose.total_residue() == 0 ) {
		pose = *segment.template_pose();
	} else {
		pose.conformation().insert_conformation_by_jump( segment.template_pose()->conformation(), pose.total_residue()+1, pose.num_jump()+1, 1 );
	}

	// add padding residues
	for ( core::Size resid=1; resid<=segment.lower_padding(); ++resid ) {
		//modify_ft_for_residue_insertion( pose, pose.total_residue() );
		TR.Debug << "Prepending residue before " << chain_start << std::endl;
		pose.prepend_polymer_residue_before_seqpos( segment.lower_residue(), chain_start, false );
		segment.lower_dihedrals().set_in_pose( pose, chain_start );
		//prepend_new_residues( pose, 1, chain_start, pose.total_residue(), segment.lower_residue().name3(), segment.lower_dihedrals() );
	}

	for ( core::Size resid=1; resid<=segment.upper_padding(); ++resid ) {
		//modify_ft_for_residue_insertion( pose, chain_start );
		TR.Debug << "Appending residue after " << pose.total_residue() << std::endl;
		pose.append_polymer_residue_after_seqpos( segment.upper_residue(), pose.total_residue(), false );
		segment.upper_dihedrals().set_in_pose( pose, pose.total_residue() - 1 );
		//append_new_residues( pose, 1, pose.total_residue(), 1, segment.upper_residue().name3(), segment.upper_dihedrals() );
	}

	core::pose::add_lower_terminus_type_to_pose_residue( pose, chain_start );
	core::pose::add_upper_terminus_type_to_pose_residue( pose, pose.total_residue() );
	if ( segment.upper_padding() ) {
		TR.Debug << "Setting upper_psi " << pose.total_residue() << " to " << segment.upper_dihedrals().upper_psi() << std::endl;
		pose.set_psi( pose.total_residue(), segment.upper_dihedrals().upper_psi() );
	}
}

void
append_residues_from_template_segment(
	core::pose::Pose & pose,
	Segment const & prev_segment,
	Segment const & segment )
{
	debug_assert( segment.template_pose() );
	TR.Debug << "Appending " << segment.template_pose()->total_residue() << " template residues into current chain" << std::endl;;
	core::Size const insert_pos = pose.total_residue();

	if ( pose.total_residue() == 0 ) {
		pose = *segment.template_pose();
	} else {
		if ( segment.lower_segment().empty() ) {
			std::stringstream msg;
			msg << "Segment does not have a lower segment, but we are appending residues.  This should never happen. Segment = "
				<< segment << std::endl;
			utility_exit_with_message( msg.str() );
		}
		bool const terminus = pose.residue( insert_pos ).is_upper_terminus();
		if ( terminus ) {
			core::pose::remove_upper_terminus_type_from_pose_residue( pose, insert_pos );
		}

		/*core::pose::Pose nonconst_template_pose( template_pose );
		for ( core::Size template_resid=1; template_resid<=nonconst_template_pose.total_residue(); ++template_resid ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( nonconst_template_pose, template_resid );
		core::pose::remove_upper_terminus_type_from_pose_residue( nonconst_template_pose, template_resid );
		}

		// removing terminus types leaves bogus phi/omega for residue 1
		// these angles are critical for correct placement of HN
		nonconst_template_pose.set_phi( 1, template_lower_phi );
		nonconst_template_pose.set_omega( 1, template_lower_omega );
		*/

		//modify_ft_for_residue_insertion( pose, insert_pos );

		core::Size pose_resid = insert_pos;
		for ( core::Size template_resid=1; template_resid<=segment.template_pose()->total_residue(); ++template_resid ) {
			pose.append_polymer_residue_after_seqpos( segment.template_pose()->residue( template_resid ), pose_resid, false );
			++pose_resid;
		}

		for ( core::Size resid=1; resid<=segment.upper_padding(); ++resid ) {
			pose.append_polymer_residue_after_seqpos( segment.upper_residue(), pose_resid, false );
			++pose_resid;
		}

		if ( terminus ) {
			core::pose::add_upper_terminus_type_to_pose_residue( pose, pose_resid );
		}
	}

	// set interface psi/psi
	pose.set_psi( insert_pos, prev_segment.upper_dihedrals().psi() );
	pose.set_omega( insert_pos, prev_segment.upper_dihedrals().omega() );
	pose.set_phi( insert_pos + 1, segment.lower_dihedrals().phi() );
	//pose.set_omega( insert_pos + 1, segment.lower_dihedrals().upper_omega() );
	TR.Debug << "Set psi of " << insert_pos << " to " << prev_segment.upper_dihedrals().psi() << std::endl;
	TR.Debug << "Set omega of " << insert_pos << " to " << prev_segment.upper_dihedrals().omega() << std::endl;
	TR.Debug << "Set phi of " << insert_pos + 1 << " to " << segment.lower_dihedrals().phi() << std::endl;
	//TR.Debug << "Set omega of " << insert_pos + 1 << " to " << segment.lower_dihedrals().omega() << std::endl;
}

/// @brief modifies teh ft in the pose, returns the original
core::kinematics::FoldTree
modify_ft_for_residue_insertion( core::pose::Pose & pose, core::Size const safe_res )
{
	TR.Debug << "re-rooting foldtree at residue " << safe_res << std::endl;
	core::kinematics::FoldTree const orig = pose.fold_tree();
	core::kinematics::FoldTree ft = pose.fold_tree();
	ft.reorder( safe_res );
	debug_assert( ft.check_fold_tree() );
	pose.fold_tree( ft );
	return orig;
}

void
set_dihedrals( core::pose::Pose & pose, core::Size const resid,  ResidueDihedrals const & dihedrals )
{
	pose.set_phi( resid, dihedrals.phi() );
	pose.set_psi( resid, dihedrals.psi() );
	pose.set_omega( resid, dihedrals.omega() );
}

void
prepend_new_residues(
	core::pose::Pose & pose,
	core::Size num_residues,
	core::Size insert_pos,
	core::Size const ,
	std::string const & type,
	ResidueDihedrals const & lower_dihedrals )
{
	TR.Debug << "Prepending " << num_residues << " at insert position " << insert_pos << std::endl;
	if ( num_residues < 1 ) return;

	//modify_ft_for_residue_insertion( pose, anchor );

	core::conformation::ResidueOP rsd = rsd_op( pose, type );

	// account for case of empty pose
	if ( pose.total_residue() == 0 ) {
		pose.append_residue_by_jump( *rsd, 1 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, 1 );
		++insert_pos;
		--num_residues;
	}

	// store lower terminus state so it can be restored below
	bool const terminus = pose.residue( insert_pos ).is_lower_terminus();
	if ( terminus ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, insert_pos );
	}

	for ( core::Size res=1; res<=num_residues; ++res ) {
		pose.prepend_polymer_residue_before_seqpos( *rsd, insert_pos, true );
	}

	if ( terminus ) {
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, insert_pos );
	}

	// set all new phi/psi/omegas
	// insert_pos is pointing to the first thing in the new pose, so
	// new_start is the start resid of the new segment.
	core::Size const stop = insert_pos + num_residues;
	for ( core::Size resid=insert_pos; resid<stop; ++resid ) {
		pose.set_phi( resid, 180.0 );
		pose.set_psi( resid, 180.0 );
		pose.set_omega( resid, 180.0 );
	}
	TR.Debug << "setting psi for " << stop - 1 << " to " << lower_dihedrals.psi() << " phi for " << lower_dihedrals.phi() << std::endl;
	lower_dihedrals.set_in_pose( pose, stop - 1 );
}

void
append_new_residues(
	core::pose::Pose & pose,
	core::Size num_residues,
	core::Size insert_pos,
	core::Size const anchor,
	std::string const & type,
	ResidueDihedrals const & upper_dihedrals )
{
	TR.Debug << "Appending " << num_residues << " at insert position " << insert_pos << std::endl;
	if ( num_residues < 1 ) return;

	core::conformation::ResidueOP rsd = rsd_op( pose, type );

	// account for case of empty pose
	if ( pose.total_residue() == 0 ) {
		pose.append_residue_by_jump( *rsd, 1 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::LOWER_TERMINUS_VARIANT, 1 );
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, 1 );
		++insert_pos;
		--num_residues;
	}

	bool const terminus = pose.residue( insert_pos ).is_upper_terminus();
	if ( terminus ) {
		core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, insert_pos );
	}

	if ( ( pose.total_residue() > 1 ) && ( anchor > 0 ) ) {
		//modify_ft_for_residue_insertion( pose, anchor );
	}

	core::Size const orig_pos = insert_pos;
	for ( core::Size res=1; res<=num_residues; ++res ) {
		pose.append_polymer_residue_after_seqpos( *rsd, insert_pos, true );
		++insert_pos;
	}

	if ( terminus ) {
		core::pose::add_variant_type_to_pose_residue( pose, core::chemical::UPPER_TERMINUS_VARIANT, insert_pos );
	}

	// set all new phi/psi/omegas
	// insert_pos is pointing to the first thing in the orig pose, so
	// we exclude it.
	// new_start is the start resid of the new segment.
	for ( core::Size resid=orig_pos+1; resid<=insert_pos; ++resid ) {
		pose.set_phi( resid, 180.0 );
		pose.set_psi( resid, 180.0 );
		pose.set_omega( resid, 180.0 );
	}
	TR.Debug << "setting psi/phi dihedrals for " << orig_pos << "/" << orig_pos+1 << " to " << upper_dihedrals.psi() << " " << upper_dihedrals.phi() << std::endl;
	upper_dihedrals.set_in_pose( pose, orig_pos );
}

} //protocols
} //denovo_design
} //components

