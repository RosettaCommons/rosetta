// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/movers/AlignResiduesMover.cc
/// @brief Aligns one residue onto another
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/AlignResiduesMover.hh>
#include <protocols/denovo_design/movers/AlignResiduesMoverCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/grafting/simple_movers/DeleteRegionMover.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/RT.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>
#include <core/select/residue_selector/ResidueVector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.AlignResiduesMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

AlignResiduesMover::AlignResiduesMover():
	protocols::moves::Mover( AlignResiduesMover::class_name() ),
	id_( "Default_Align_Residues_Mover_Name" ),
	template_selectors_(),
	target_selectors_()
{
}

AlignResiduesMover::~AlignResiduesMover(){}

AlignResiduesMover::ResidueSelectorCOPs
parse_residue_selectors( std::string const & selector_name_str, basic::datacache::DataMap const & data )
{
	AlignResiduesMover::ResidueSelectorCOPs selectors;
	if ( selector_name_str.empty() ) return selectors;

	utility::vector1< std::string > const selector_names = utility::string_split( selector_name_str, ',' );
	for ( utility::vector1< std::string >::const_iterator sname=selector_names.begin(); sname!=selector_names.end(); ++sname ) {
		core::select::residue_selector::ResidueSelectorCOP s = protocols::rosetta_scripts::get_residue_selector( *sname, data );
		if ( !s ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "target_selector named " + *sname + " could not be found in the RosettaScripts XML" );
		}
		selectors.push_back( s->clone() );
	}
	return selectors;
}

void
AlignResiduesMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & ,
	core::pose::Pose const & )
{
	id_ = tag->getOption< std::string >( "name" );

	std::string const template_selector_names = tag->getOption< std::string >( "template_selectors" );
	template_selectors_ = parse_residue_selectors( template_selector_names, data );


	std::string const target_selector_names = tag->getOption< std::string >( "target_selectors" );
	target_selectors_ = parse_residue_selectors( target_selector_names, data );

	if ( template_selectors_.size() != target_selectors_.size() ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "Number of template ids must match number of target ids." );
	}
}

protocols::moves::MoverOP
AlignResiduesMover::clone() const
{
	return protocols::moves::MoverOP( new AlignResiduesMover( *this ) );
}

protocols::moves::MoverOP
AlignResiduesMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new AlignResiduesMover );
}

std::string
AlignResiduesMover::get_name() const
{
	return AlignResiduesMover::class_name();
}

std::string
AlignResiduesMover::class_name()
{
	return "AlignResiduesMover";
}

void
AlignResiduesMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, AlignResiduesMover const & mover )
{
	mover.show(os);
	return os;
}

void
AlignResiduesMover::apply( core::pose::Pose & pose )
{
	if ( template_selectors_.size() != target_selectors_.size() ) {
		std::stringstream msg;
		msg << "AlignResiduesMover::apply(): Number of template selectors (" << template_selectors_.size()
			<< " must match number of target selectors (" << target_selectors_.size() << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( template_selectors_.empty() ) {
		utility_exit_with_message( "AlignResiduesMover::apply(): Template selectors must not be empty.\n" );
	}
	if ( target_selectors_.empty() ) {
		utility_exit_with_message( "AlignResiduesMover::apply(): Target selectors must not be empty.\n" );
	}

	ResidueVectors const template_subsets = compute_template_residues( pose );
	ResidueVectors const target_subsets = compute_target_residues( pose );

	core::select::residue_selector::ResidueVector template_resids;
	for ( ResidueVectors::const_iterator temp=template_subsets.begin(), targ=target_subsets.begin();
			(temp!=template_subsets.end()) && (targ!=target_subsets.end()); ++temp, ++targ ) {
		template_resids.push_back( align_residues( pose, template_resids.size() + 1, *temp, *targ ) );
	}

	// delete template residues
	delete_residues( pose, template_resids );
}

AlignResiduesMover::ResidueVectors
compute_residues_from_selectors( core::pose::Pose const & pose, AlignResiduesMover::ResidueSelectorCOPs const & selectors )
{
	AlignResiduesMover::ResidueVectors retval;
	for ( AlignResiduesMover::ResidueSelectorCOPs::const_iterator s=selectors.begin(); s!=selectors.end(); ++s ) {
		debug_assert( *s );
		retval.push_back( core::select::residue_selector::ResidueVector( (*s)->apply( pose ) ) );
	}
	return retval;
}

void
AlignResiduesMover::copy_residue(
	core::pose::Pose & pose,
	core::Size const resid_has_info,
	core::Size const resid_wants_info ) const
{
	//core::conformation::Residue const orig_rsd = pose.residue( resid_wants_info );
	bool lower_term = ( pose.residue( resid_has_info ).is_lower_terminus() );
	bool upper_term = ( pose.residue( resid_has_info ).is_upper_terminus() );

	if ( lower_term ) {
		core::pose::remove_lower_terminus_type_from_pose_residue( pose, resid_has_info );
	}
	if ( upper_term ) {
		core::pose::remove_upper_terminus_type_from_pose_residue( pose, resid_has_info );
	}
	//rsd_with_info.clear_residue_connections();
	//rsd_with_info.copy_residue_connections( pose.residue( resid_wants_info ) );
	pose.replace_residue( resid_wants_info, pose.residue( resid_has_info ), true );

	if ( lower_term ) {
		core::pose::add_lower_terminus_type_to_pose_residue( pose, resid_has_info );
	}
	if ( upper_term ) {
		core::pose::add_upper_terminus_type_to_pose_residue( pose, resid_has_info );
	}
	TR.Debug << "Copied residue " << resid_has_info << " into " << resid_wants_info << std::endl;
	/*
	for ( core::Size conn_id=1; conn_id<=orig_rsd.n_current_residue_connections(); ++conn_id ) {
	TR.Debug << "Connections=" << pose.residue(resid_wants_info).n_current_residue_connections()
	<< " conn_id=" << conn_id << std::endl;
	core::Size const partner_resid = orig_rsd.connected_residue_at_resconn( conn_id );
	TR.Debug << "Updating partner=" << partner_resid << std::endl;
	debug_assert( partner_resid );
	core::conformation::Residue partner = pose.residue( partner_resid );
	partner.update_connections_to_other_residue( pose.residue( resid_wants_info ) ); //Update the connections
	pose.replace_residue( partner_resid, partner, false );
	}
	*/
}

void
AlignResiduesMover::delete_residues(
	core::pose::Pose & pose,
	core::select::residue_selector::ResidueVector const & resids ) const
{
	if ( resids.empty() ) return;

	std::stringstream resid_str;
	for ( core::select::residue_selector::ResidueVector::const_iterator r=resids.begin(); r!=resids.end(); ++r ) {
		if ( r != resids.begin() ) resid_str << ',';
		resid_str << *r;
	}

	core::select::residue_selector::ResidueSelectorCOP index_sel(
		new core::select::residue_selector::ResidueIndexSelector( resid_str.str() ) );

	protocols::grafting::simple_movers::DeleteRegionMover delete_mover;
	delete_mover.set_residue_selector( index_sel );

	debug_assert( components::StructureDataFactory::get_instance()->observer_attached( pose ) );
	delete_mover.apply( pose );
}

AlignResiduesMover::ResidueVectors
AlignResiduesMover::compute_target_residues( core::pose::Pose const & pose ) const
{
	return compute_residues_from_selectors( pose, target_selectors_ );
}

AlignResiduesMover::ResidueVectors
AlignResiduesMover::compute_template_residues( core::pose::Pose const & pose ) const
{
	return compute_residues_from_selectors( pose, template_selectors_ );
}

core::Size
AlignResiduesMover::align_residues(
	core::pose::Pose & pose,
	core::Size const align_count,
	core::select::residue_selector::ResidueVector const & template_resids,
	core::select::residue_selector::ResidueVector const & target_resids ) const
{
	if ( template_resids.size() != 1 ) {
		std::stringstream msg;
		msg << class_name() << "::align_residues(): Exactly one residue must be selected by the template selector. "
			<< "You have selected the following residues: " << template_resids << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size const idx = numeric::random::rg().random_range( 1, target_resids.size() );

	align_residues( pose, align_count, *template_resids.begin(), target_resids[idx] );

	return *template_resids.begin();
}

void
AlignResiduesMover::align_residues(
	core::pose::Pose & pose,
	core::Size const align_count,
	core::Size const template_resid,
	core::Size const target_resid ) const
{
	TR << "Aligning " << target_resid << " to " << template_resid << std::endl;

	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	components::StructureData sd = factory.get_from_pose( pose );
	std::string const template_segment = sd.segment_name( template_resid );
	std::string const target_segment = sd.segment_name( target_resid );
	SegmentNames const roots = boost::assign::list_of
		(template_segment)(target_segment);

	FoldTreeFromFoldGraphMover create_ft( roots, protocols::loops::Loops() );
	create_ft.apply( pose );

	int const jump_idx = find_jump_rec( pose.fold_tree(), sd.segment( target_segment ).safe() );
	debug_assert( jump_idx );

	// change target_resid if doing reverse aligning
	align_residues( pose, jump_idx, template_resid, target_resid, target_resid );

	// add metadata to StructureData
	add_metadata( sd, align_count, template_segment, target_segment, target_resid );

	// copy aligned residue to aligned position
	copy_residue( pose, template_resid, target_resid );

	// save SD back into pose
	factory.save_into_pose( pose, sd );
}

/// @brief aligns two residues so their backbones are completely superimposed
void
AlignResiduesMover::align_residues(
	core::pose::Pose & pose,
	core::Size const jump_idx,
	core::Size const template_resid,
	core::Size const target_resid,
	core::Size const res_with_torsions ) const
{
	// move the jump temporarily
	core::kinematics::FoldTree const ft_backup = pose.fold_tree();
	pose.fold_tree( slide_jump( ft_backup, jump_idx, template_resid, target_resid ) );

	//setup torsions
	if ( res_with_torsions ) {
		pose.set_phi( template_resid, pose.phi(res_with_torsions) );
		pose.set_psi( template_resid, pose.psi(res_with_torsions) );
		pose.set_phi( target_resid, pose.phi(res_with_torsions) );
		pose.set_psi( target_resid, pose.psi(res_with_torsions) );
	}

	core::id::StubID stub1(
		core::id::AtomID(pose.residue(template_resid).type().atom_index("CA"), template_resid),
		core::id::AtomID(pose.residue(template_resid).type().atom_index("N"), template_resid),
		core::id::AtomID(pose.residue(template_resid).type().atom_index("C"), template_resid) );

	core::id::StubID stub2(
		core::id::AtomID(pose.residue(target_resid).type().atom_index("CA"), target_resid),
		core::id::AtomID(pose.residue(target_resid).type().atom_index("N"), target_resid),
		core::id::AtomID(pose.residue(target_resid).type().atom_index("C"), target_resid) );

	pose.conformation().set_stub_transform( stub1, stub2, core::kinematics::RT() );
	debug_assert( ft_backup.check_fold_tree() );
	pose.fold_tree( ft_backup );
}

void
AlignResiduesMover::add_metadata(
	components::StructureData & sd,
	core::Size const align_count,
	std::string const & template_segment,
	std::string const & target_segment,
	core::Size const target_resid ) const
{
	sd.renumber_movable_group(
		sd.segment( target_segment ).movable_group(),
		sd.segment( template_segment ).movable_group() );

	std::stringstream name_ss;
	name_ss << id_ << PARENT_DELIMETER << align_count;

	sd.set_alias( name_ss.str(), target_resid );
}

void
AlignResiduesMover::add_template_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	template_selectors_.push_back( selector.clone() );
}

void
AlignResiduesMover::add_target_selector( core::select::residue_selector::ResidueSelector const & selector )
{
	target_selectors_.push_back( selector.clone() );
}

void
AlignResiduesMover::set_id( std::string const & idval )
{
	id_ = idval;
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
AlignResiduesMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new AlignResiduesMover );
}

std::string
AlignResiduesMoverCreator::keyname() const
{
	return AlignResiduesMover::class_name();
}

} //protocols
} //denovo_design
} //movers

