// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/StructureDataWithPose.cc
/// @brief StructureData with pose attached
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/components/StructureDataWithPose.hh>

// Protocol headers
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/forge/methods/pose_mod.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/moves/Mover.hh>

// Core headers
#include <core/conformation/Conformation.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/OptCysHG.hh>
#include <core/pose/util.hh>
#include <core/scoring/Energies.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Basic/Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.StructureDataWithPose" );

namespace protocols {
namespace denovo_design {
namespace components {

StructureDataWithPose::StructureDataWithPose( StructureData const & sd, core::pose::Pose const & pose ):
	StructureData( sd ),
	pose_( pose.clone() )
{
}

/// @brief copy constructor -- poseOP is cloned
StructureDataWithPose::StructureDataWithPose( StructureDataWithPose const & perm ):
	StructureData( perm ),
	pose_( pose().clone() )
{
}

StructureDataWithPose::~StructureDataWithPose()
{}

StructureDataOP
StructureDataWithPose::fresh_instance() const
{
	return StructureDataOP( new StructureDataWithPose( id(), core::pose::Pose() ) );
}

StructureDataOP
StructureDataWithPose::clone() const
{
	return StructureDataOP( new StructureDataWithPose( *this ) );
}

/// @brief just return the pose
core::pose::Pose const &
StructureDataWithPose::pose() const
{
	debug_assert( pose_ );
	return *pose_;
}

/// @brief return nonconst pose
core::pose::Pose &
StructureDataWithPose::pose_nonconst()
{
	debug_assert( pose_ );
	return *pose_;
}

/// @brief sets the pose and does checks to ensure data is consistent
/// @throws EXCN_PoseInconsistent if the pose doesn't match up with SD
void
StructureDataWithPose::set_pose( core::pose::Pose const & new_pose )
{
	set_pose( new_pose.clone() );
}

/// @brief sets the pose and does checks to ensure data is consistent
/// @throws EXCN_PoseInconsistent if the pose doesn't match up with SD
void
StructureDataWithPose::set_pose( core::pose::PoseOP new_pose )
{
	if ( new_pose && new_pose->size() ) {
		TR.Debug << id() << ": new_pose len=" << new_pose->size() << " pose residue =" << pose_length() << " perm len=" << length() << std::endl;

		check_pose_consistency( *new_pose );

		// put terminal variants onto chain endings
		utility::vector1< core::Size > const conformation_chainend = new_pose->conformation().chain_endings();
		utility::vector1< core::Size > chain_endings = conformation_chainend;
		chain_endings.push_back( new_pose->size() );
		utility::vector1< core::Size > chain_starts;
		if ( new_pose->size() ) {
			chain_starts.push_back( 1 );
		}
		for ( utility::vector1< core::Size >::const_iterator r = conformation_chainend.begin(); r != conformation_chainend.end(); ++r ) {
			if ( *r == new_pose->size() ) continue;
			chain_starts.push_back( *r + 1 );
		}
		TR.Debug << "Chain endings: " << chain_starts << chain_endings << std::endl;

		for ( utility::vector1< core::Size >::const_iterator r = chain_starts.begin(); r != chain_starts.end(); ++r ) {
			debug_assert( *r > 0 );
			debug_assert( *r <= new_pose->size() );
			if ( new_pose->residue( *r ).is_polymer() &&
					( ! new_pose->residue( *r ).is_lower_terminus() ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) ) {
				core::pose::add_lower_terminus_type_to_pose_residue( *new_pose, *r );
				TR.Debug << "Added lower terminus type to residue " << *r << std::endl;
			} else {
				TR.Debug << "Residue " << new_pose->residue( *r ).name() << *r << " is already a lower terminus or non-polymer." << std::endl;
			}
		}
		for ( utility::vector1< core::Size >::const_iterator r = chain_endings.begin(); r != chain_endings.end(); ++r ) {
			debug_assert( *r > 0 );
			debug_assert( *r <= new_pose->size() );
			if ( new_pose->residue( *r ).is_polymer() &&
					( ! new_pose->residue( *r ).is_upper_terminus() ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) &&
					( ! new_pose->residue( *r ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) ) {
				core::pose::add_upper_terminus_type_to_pose_residue( *new_pose, *r );
				TR.Debug << "Added upper terminus type to residue " << *r << std::endl;
			} else {
				TR.Debug << "Residue " << new_pose->residue( *r ).name() << *r << " is already a upper terminus or non-polymer." << std::endl;
			}
		}
		new_pose->conformation().chains_from_termini();
	}
	pose_ = new_pose;
	save_into_pose();
}

void
StructureDataWithPose::on_change()
{
	save_into_pose();
}

/// @brief checks consistency of the data
/// @throws EXCN_PoseInconsistent if there is a problem
void
StructureDataWithPose::check_consistency() const
{
	check_pose_consistency( pose() );
}

/// @brief merge given data and segments from "other" into this StructureData
void
StructureDataWithPose::merge( StructureDataWithPose const & other, StringList const & segments )
{
	core::Size const movable_group_new = choose_new_movable_group();
	for ( StringList::const_iterator s=segments.begin(); s!=segments.end(); ++s ) {
		Segment newseg = other.segment( *s );
		newseg.set_movable_group( movable_group_new );
		core::pose::Pose newpose( other.pose(), other.segment( *s ).lower(), other.segment( *s ).upper() );
		add_segment( *s, newseg, segments_end_nonconst(), newpose );
		TR.Debug << "Added " << NamedSegment( *s, newseg ) << std::endl;
	}
	copy_data( other );
}

/// @brief merge given data and segments from "other" into this StructureData
void
StructureDataWithPose::merge_before( StructureDataWithPose const & other, std::string const & position, StringList const & segments )
{
	core::Size const movable_group_new = choose_new_movable_group();
	for ( StringList::const_iterator s=segments.begin(); s!=segments.end(); ++s ) {
		Segment newseg = other.segment( *s );
		newseg.set_movable_group( movable_group_new );
		core::pose::Pose newpose( other.pose(), other.segment( *s ).lower(), other.segment( *s ).upper() );
		add_segment( *s, newseg, find_segment_name( position ), newpose );
		TR.Debug << "Added " << NamedSegment( *s, newseg ) << std::endl;
	}
	copy_data( other );
}

void
StructureDataWithPose::save_into_pose()
{
	StructureDataFactory::get_instance()->save_into_pose( pose_nonconst(), *this );
}

/// @brief stores a string in the pose's datacache
bool
StructureDataWithPose::has_cached_string() const
{
	return StructureDataFactory::get_instance()->has_cached_string( pose() );
}

/// @brief stores a string in the pose's datacache
std::string
StructureDataWithPose::cached_string() const
{
	return StructureDataFactory::get_instance()->retrieve_cached_string( pose() );
}

/// @brief stores a string in the pose's datacache
void
StructureDataWithPose::set_cached_string( std::string const & ss )
{
	StructureDataFactory::get_instance()->set_cached_string( pose_nonconst(), ss );
}

/// @brief retrieves cached remarks from pose datacache
core::io::Remarks
StructureDataWithPose::retrieve_remarks() const
{
	return StructureData::retrieve_remarks( pose() );
}

/// @brief applies an arbitrary mover to the contained pose
void
StructureDataWithPose::apply_mover( protocols::moves::Mover & mover )
{
	core::Size const orig_len = pose().size();
	mover.apply( pose_nonconst() );
	debug_assert( orig_len == pose().size() );
}

/// @brief applies an arbitrary mover to the contained pose
void
StructureDataWithPose::apply_mover( protocols::moves::MoverOP mover )
{
	debug_assert( mover );
	apply_mover( *mover );
}

/// @brief switches residue type set of the contained pose
void
StructureDataWithPose::switch_residue_type_set( std::string const & typeset )
{
	set_fold_tree( remove_missing_jump_atoms( pose(), pose().fold_tree() ) );
	core::util::switch_to_residue_type_set( pose_nonconst(), typeset );
}

void
StructureDataWithPose::set_fold_tree( core::kinematics::FoldTree const & ft )
{
	debug_assert( ft.check_fold_tree() );
	pose_nonconst().fold_tree( ft );
}

/// @brief returns the chain number of the given residue
core::Size
StructureDataWithPose::chain( core::Size const resid ) const
{
	return pose().chain( resid );
}

/// @brief based on the pose, determines the jump number for the jump pointing to the segment specified
int
StructureDataWithPose::find_jump( std::string const & seg ) const
{
	int const j = find_jump_rec( pose().fold_tree(), segment( seg ).safe() );
	debug_assert( j >= 0 );
	return j;
}

/// @brief creates a new jump and cutpoint to build the given loop object
int
StructureDataWithPose::new_jump_and_cutpoint( protocols::loops::Loop const & loop, core::Size const loop_overlap )
{
	core::Size const startres = loop_start_without_overlap( pose(), loop.start(), loop_overlap );
	core::Size const stopres = loop_stop_without_overlap( pose(), loop.stop(), loop_overlap );
	core::Size const cutres = loop.cut();
	core::Size saferes1_dist = 0;
	core::Size saferes1 = 0;

	// find the segment containing the cutpoint and set it
	for ( StringList::const_iterator s=segments_begin(); s!=segments_end(); ++s ) {
		Segment & res = segment_nonconst( *s );
		if ( res.contains( cutres ) ) {
			res.set_cutpoint( cutres - res.lower() + 1 );
			break;
		}
	}
	// find the closest two safe residues to the loop and add a jump between them.
	for ( SegmentMap::const_iterator res=segments().begin(); res!=segments().end(); ++res ) {
		Segment const & resis = res->second;
		if ( resis.safe() >= startres ) {
			continue;
		}
		if ( !saferes1 || ( startres-resis.safe() < saferes1_dist ) ) {
			saferes1 = resis.safe();
			saferes1_dist = startres-saferes1;
		}
	}
	core::Size saferes2_dist = 0;
	core::Size saferes2 = 0;
	for ( SegmentMap::const_iterator res=segments().begin(); res!=segments().end(); ++res ) {
		Segment const & resis = res->second;
		if ( resis.safe() <= stopres ) {
			continue;
		}
		if ( !saferes2 || ( resis.safe()-stopres < saferes2_dist ) ) {
			saferes2 = resis.safe();
			saferes2_dist = saferes2-stopres;
		}
	}

	debug_assert( saferes1 < cutres );
	debug_assert( saferes2 > cutres );
	TR << "Inserting jump " << saferes1 << "__" << saferes2 << " with cut at " << cutres << std::endl;

	// modify fold tree
	core::kinematics::FoldTree ft = pose().fold_tree();
	int const newjump = ft.new_jump( saferes1, saferes2, cutres );
	if ( !ft.check_fold_tree() ) {
		TR.Error << "FOLDTREE=" << ft << std::endl;
		TR.Error << "StructureData=" << *this << std::endl;
	}
	debug_assert( ft.check_fold_tree() );
	pose_nonconst().fold_tree( ft );

	// remove terminal variants (if applicable)
	if ( pose_->residue(cutres).is_upper_terminus() ) {
		remove_upper_terminus_variant_type( cutres );
	}
	if ( pose_->residue(cutres+1).is_lower_terminus() ) {
		remove_lower_terminus_variant_type( cutres+1 );
	}
	// add cutpoint variant types
	add_cutpoint_variants( cutres );
	rebuild_missing_atoms( *pose_, cutres );
	rebuild_missing_atoms( *pose_, cutres+1 );

	chains_from_termini();

	// they should ALWAYS be on the same chain at this point
	TR.Debug << pose_->fold_tree() << std::endl;
	TR.Debug << "chain of " << saferes1 << " = " << pose_->chain( saferes1 ) << " chain of " << saferes2 << " = " << pose_->chain( saferes2 ) << std::endl;
	for ( core::Size i=1, endi=pose_->size(); i<=endi; ++i ) {
		TR.Debug << i << " " << pose_->residue(i).name() << std::endl;
	}
	debug_assert( pose_->chain( saferes1 ) == pose_->chain( saferes2 ) );

	return newjump;
}

/// @brief deletes the residues between the segment N terminus and the N anchor point
void
StructureDataWithPose::delete_leading_residues( std::string const & seg )
{
	core::Size const start = segment( seg ).lower();
	core::Size const stop = segment( seg ).start() - 1;
	TR.Debug << "Lead delete start=" << start << " stop=" << stop << std::endl;
	if ( start > stop ) return;
	StructureData::delete_leading_residues( seg );
	delete_residues_in_pose( start, stop );
}

/// @brief deletes the residues between the segment C terminus and the C anchor point
void
StructureDataWithPose::delete_trailing_residues( std::string const & seg )
{
	core::Size const start = segment( seg ).stop() + 1;
	core::Size const stop = segment( seg ).upper();
	TR.Debug << "Trail delete start=" << start << " stop=" << stop << std::endl;
	if ( start > stop ) return;
	StructureData::delete_trailing_residues( seg );
	delete_residues_in_pose( start, stop );
}

/// @brief determine pose chains based on termini
void
StructureDataWithPose::chains_from_termini()
{
	pose_nonconst().conformation().chains_from_termini();
}

/// @brief adds upper terminal variant to a residue
void
StructureDataWithPose::add_upper_terminus_variant_type( core::Size const resi )
{
	if ( !pose().residue( resi ).is_upper_terminus() ) {
		core::pose::add_variant_type_to_pose_residue( pose_nonconst(), core::chemical::UPPER_TERMINUS_VARIANT, resi );
	}
}

/// @brief adds lower terminal variant to a residue
void
StructureDataWithPose::add_lower_terminus_variant_type( core::Size const resi )
{
	if ( !pose().residue( resi ).is_lower_terminus() ) {
		core::pose::add_variant_type_to_pose_residue( pose_nonconst(), core::chemical::LOWER_TERMINUS_VARIANT, resi );
	}
}

/// @brief removes upper terminal variant from a residue
void
StructureDataWithPose::remove_upper_terminus_variant_type( core::Size const resi )
{
	if ( pose().residue( resi ).is_upper_terminus() ) {
		core::pose::remove_variant_type_from_pose_residue( pose_nonconst(), core::chemical::UPPER_TERMINUS_VARIANT, resi );
	}
}

/// @brief removes lower terminal variant from a residue
void
StructureDataWithPose::remove_lower_terminus_variant_type( core::Size const resi )
{
	if ( pose().residue( resi ).is_lower_terminus() ) {
		core::pose::remove_variant_type_from_pose_residue( pose_nonconst(), core::chemical::LOWER_TERMINUS_VARIANT, resi );
	}
}

/// @brief add lower cutpoint to residue cut and upper cutpoint to residue cut+1
void
StructureDataWithPose::add_cutpoint_variants( core::Size const cut_res )
{
	protocols::forge::methods::add_cutpoint_variants( pose_nonconst(), cut_res );
}

/// @brief removes cutpoint variants from residues cut and cut+1
void
StructureDataWithPose::remove_cutpoint_variants( core::Size const cut_res )
{
	debug_assert( cut_res + 1 <= pose().size() );
	TR.Debug << "Removing cutpoints from " << cut_res << " and " << cut_res + 1 << std::endl;

	if ( pose().residue( cut_res ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ) {
		core::pose::remove_variant_type_from_pose_residue( pose_nonconst(), core::chemical::CUTPOINT_LOWER, cut_res );
	}

	if ( pose().residue( cut_res + 1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
		core::pose::remove_variant_type_from_pose_residue( pose_nonconst(), core::chemical::CUTPOINT_UPPER, cut_res + 1 );
	}
}

/// @brief returns true of the last residue of segment1 contains a covalent bond to the first residue of segment2
bool
StructureDataWithPose::polymer_bond_exists( std::string const & segment1, std::string const & segment2 ) const
{
	return pose().residue( segment(segment1).upper() ).is_polymer_bonded( segment(segment2).lower() );
}

/// @brief (re-)detect disulfides, will convert CYD to CYS if disulfide bond is lost
void
StructureDataWithPose::detect_disulfides( core::scoring::ScoreFunctionOP sfx )
{
	debug_assert( sfx );
	if ( basic::options::option[ basic::options::OptionKeys::in::detect_disulf ].user() ?
			basic::options::option[ basic::options::OptionKeys::in::detect_disulf ]() : // detect_disulf true
			pose().is_fullatom() // detect_disulf default but fa pose
			) {
		pose_nonconst().conformation().detect_disulfides();
	}

	// fix HG of CYS to relieve clashes of any newly converted CYS
	core::pack::task::TaskFactoryOP tf( new core::pack::task::TaskFactory() );
	tf->push_back( core::pack::task::operation::TaskOperationOP( new core::pack::task::operation::OptCysHG() ) );
	core::pack::pack_rotamers( pose_nonconst(), *sfx, tf->create_task_and_apply_taskoperations( pose() ) );

	// safety
	pose_nonconst().energies().clear();
}

void
StructureDataWithPose::update_covalent_bonds_in_pose()
{
	for ( BondInfos::const_iterator p=covalent_bonds_begin(); p!=covalent_bonds_end(); ++p ) {
		TR.Debug << "Updating covalent bond " << p->seg1 << "-->" << p->seg2 << std::endl;
		declare_covalent_bond_in_pose( pose_residue( p->seg1, p->res1 ), p->atom1, pose_residue( p->seg2, p->res2 ), p->atom2 );
	}
}

/// @brief declares a covalent bond using pose residues
void
StructureDataWithPose::declare_covalent_bond(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	StructureData::declare_covalent_bond( res1, atom1, res2, atom2 );
	declare_covalent_bond_in_pose( res1, atom1, res2, atom2 );
}

/// @brief simply declares a covalent bond in the pose
void
StructureDataWithPose::declare_covalent_bond_in_pose(
	core::Size const res1, std::string const & atom1,
	core::Size const res2, std::string const & atom2 )
{
	if ( pose().residue( res1 ).has( atom1 ) && pose().residue( res2 ).has( atom2 ) ) {
		pose_nonconst().conformation().declare_chemical_bond( res1, atom1, res2, atom2 );
		//Rebuild the connection atoms:
		for ( core::Size i=1; i<=pose().conformation().residue(res1).n_possible_residue_connections(); ++i ) {
			pose_nonconst().conformation().rebuild_residue_connection_dependent_atoms( res1, i );
		}
		for ( core::Size i=1; i<=pose().conformation().residue(res2).n_possible_residue_connections(); ++i ) {
			pose_nonconst().conformation().rebuild_residue_connection_dependent_atoms( res2, i );
		}
	} else {
		TR.Debug << "Skipping covalent bond declaration between " << res1 << ":" << atom1
			<< " and " << res2 << ":" << atom2 << " due to missing atoms. The residue identity may have changed." << std::endl;
	}
}

/// @brief declares a covalent bond between the c-terminal atom of segment 1 and the n-terminal atom of segment 2
void
StructureDataWithPose::declare_polymer_bond( std::string const & segment1, std::string const & segment2 )
{
	core::Size const pos1 = segment(segment1).upper();
	core::Size const pos2 = segment(segment2).lower();
	debug_assert( pose().residue(pos1).is_polymer() );
	debug_assert( pose().residue(pos2).is_polymer() );
	debug_assert( ! pose().residue(pos1).is_polymer_bonded( pos2 ) );

	declare_covalent_bond( pos1, pose().residue(pos1).atom_name( pose().residue(pos1).upper_connect_atom() ),
		pos2, pose().residue(pos2).atom_name( pose().residue(pos2).lower_connect_atom() ) );

	debug_assert( pose().residue(pos1).is_polymer_bonded( pos2 ) );
}

/// @brief connects the given chains together -- don't call this, doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureDataWithPose::connect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	StructureData::connect_segments( segment1_c, segment2_n );
	remove_upper_terminus_variant_type( segment(segment1_c).upper() );
	remove_lower_terminus_variant_type( segment(segment2_n).lower() );
	if ( ! polymer_bond_exists( segment1_c, segment2_n ) ) {
		declare_polymer_bond( segment1_c, segment2_n );
	}
	chains_from_termini();
	TR.Debug << "Marked " << segment1_c << " and " << segment2_n << " as connected." << std::endl;
	save_into_pose();
}

/// @brief disconnects the given chains doesn't update anything -- don't call this on its own unless you know what you're doing.
void
StructureDataWithPose::disconnect_segments(
	std::string const & segment1_c,
	std::string const & segment2_n )
{
	StructureData::disconnect_segments( segment1_c, segment2_n );
	add_upper_terminus_variant_type( segment( segment1_c ).upper() );
	add_lower_terminus_variant_type( segment( segment2_n ).lower() );
	chains_from_termini();
}

/// @brief merges two segments into one that has the name new_name. They must be next to each other in sequence.
void
StructureDataWithPose::merge_segments(
	std::string const & segment1,
	std::string const & segment2,
	std::string const & new_name )
{
	StructureData::merge_segments( segment1, segment2, new_name );

	remove_upper_terminus_variant_type( segment(segment1).upper() );
	remove_lower_terminus_variant_type( segment(segment2).lower() );
	chains_from_termini();
	save_into_pose();
}

/// @brief sets jump with index jumpidx to datastructure j
void
StructureDataWithPose::set_jump( int const jumpidx, core::kinematics::Jump const & j )
{
	pose_nonconst().set_jump( jumpidx, j );
}

/// @brief sets phi for a specific residue
void
StructureDataWithPose::set_phi( core::Size const seqpos, core::Real const phi_val )
{
	pose_nonconst().set_phi( seqpos, phi_val );
}

/// @brief sets psi for a specific residue
void
StructureDataWithPose::set_psi( core::Size const seqpos, core::Real const psi_val )
{
	pose_nonconst().set_psi( seqpos, psi_val );
}

/// @brief sets omega for a specific residue
void
StructureDataWithPose::set_omega( core::Size const seqpos, core::Real const omega_val )
{
	pose_nonconst().set_omega( seqpos, omega_val );
}

/// @brief sets bond length in pose
void
StructureDataWithPose::set_bond_length(
	std::string const & segmentname,
	core::Size const resid,
	std::string const & atom1,
	std::string const & atom2,
	core::Real const newlength )
{
	debug_assert( newlength >= 0 );
	debug_assert( has_segment( segmentname ) );
	Segment const & res = segment(segmentname);
	debug_assert( resid );
	debug_assert( resid <= res.length() );
	core::Size const poseres = res.segment_to_pose(resid);
	debug_assert( poseres );
	debug_assert( poseres <= pose().size() );
	debug_assert( pose().residue(poseres).type().has( atom1 ) );
	debug_assert( pose().residue(poseres).type().has( atom2 ) );
	core::id::AtomID const atom1id( pose().residue(poseres).type().atom_index(atom1), poseres );
	core::id::AtomID const atom2id( pose().residue(poseres).type().atom_index(atom2), poseres );
	pose_nonconst().conformation().set_bond_length( atom1id, atom2id, newlength );

}

/// @brief moves around jumps so they are located in the center of rigid residue groups
void
StructureDataWithPose::move_jumps_to_safety()
{
	core::kinematics::FoldTree ft = pose().fold_tree();

	for ( core::Size i=1; i<=ft.num_jump(); ++i ) {
		core::kinematics::Edge e = ft.jump_edge(i);
		TR << "Jump " << e << " connects " << segment_name( e.start() ) << " to " << segment_name( e.stop() ) << std::endl;
		std::string const upname = segment_name( e.start() );
		e.start() = segment( upname ).safe();
		std::string const downname = segment_name( e.stop() );
		e.stop() = segment( downname ).safe();
		TR.Debug << "FT: " << ft << std::endl;
		TR.Debug << "Sliding jump " << e << std::endl;
		ft.slide_jump( e.label(), e.start(), e.stop() );
	}
	if ( !ft.check_fold_tree() ) {
		TR << "FAILING fold tree check:" << *this << std::endl;
		debug_assert( ft.check_fold_tree() );
	}
	pose_nonconst().fold_tree( ft );
}

/// @brief moves around jumps so that movable groups will all move together during folding
void
StructureDataWithPose::consolidate_movable_groups( utility::vector1< std::string > const & root_segments )
{
	consolidate_movable_groups( pose_nonconst(), root_segments );
}

/// @brief moves around jumps so that movable groups will all move together during folding
void
StructureDataWithPose::consolidate_movable_groups(
	core::pose::Pose & pose,
	utility::vector1< std::string > const & root_segments )
{
	// make fold graph
	FoldGraph fg( *this );
	core::kinematics::FoldTree ft = fg.fold_tree( root_segments );
	if ( ft.nres() != pose.size() ) {
		TR.Error << "FT nres " << ft.nres() << " != pose nres " << pose.size() << " FT: " << ft << std::endl;
		throw utility::excn::EXCN_Msg_Exception( "Bad nres" );
	}
	if ( pose.is_centroid() ) {
		TR << "Pose is centroid" << std::endl;
		ft = remove_missing_jump_atoms( pose, ft );
	}
	TR.Debug << "Created fold tree: " << ft << std::endl;
	TR.Debug << "FT Perm = " << *this << std::endl;
	pose.fold_tree( ft );
}


/// @brief removes jump and cutpoint between the two segments to create a single polymer chain
void
StructureDataWithPose::delete_jump_and_intervening_cutpoint( std::string const & segment1, std::string const & segment2 )
{
	StructureData::delete_jump_and_intervening_cutpoint( segment1, segment2 );

	core::Size const res1 = segment(segment1).upper();
	core::Size const res2 = segment(segment2).lower();
	int const jnum = find_jump(segment2);
	// check for cyclic case -- same jump for both segments.  Exit if cyclic
	if ( ( jnum == find_jump(segment1) ) || ( segment( segment1 ).upper() + 1 != segment( segment2 ).lower() ) ) {
		TR << "Can't delete jump/cutpoint between " << segment1 << ", res " << segment( segment1 ).safe() << " and " << segment2 << ", res " << segment( segment2 ).safe() << " because the polymer is apparently cyclic. Jump1=" << find_jump( segment1 ) << " Jump2=" << jnum << " FT=" << pose().fold_tree() << " perm=" << *this << std::endl;
		return;
	}
	delete_jump_and_intervening_cutpoint( jnum, res1, res2 );
}

/// @brief given a jump and a cutpoint residue, move the jump around the cut residue and delete to form a single edge
void
StructureDataWithPose::delete_jump_and_intervening_cutpoint(
	int const jnum,
	core::Size const cut_resi1,
	core::Size const cut_resi2 )
{
	// remove cutpoint variants from pose
	remove_cutpoint_variants( cut_resi1 );

	//move jump so it surrounds cutpoint
	core::kinematics::FoldTree ft = slide_jump( pose().fold_tree(), jnum, cut_resi1, cut_resi2 );

	// delete jump+cutpoint
	ft.delete_jump_and_intervening_cutpoint( jnum );
	debug_assert( ft.check_fold_tree() );
	pose_nonconst().fold_tree( ft );
}

/// @brief replaces one residue with another
void
StructureDataWithPose::replace_residue(
	std::string const & target_segment,
	core::Size const target_res,
	core::conformation::Residue const & res_in )
{
	core::Size const resnum = segment(target_segment).segment_to_pose(target_res);
	replace_residue( resnum, res_in, true );
}

/// @brief replaces one residue with another
void
StructureDataWithPose::replace_residue(
	core::Size const resnum,
	core::conformation::Residue const & res_in,
	bool const orient_bb )
{
	bool is_upper_term = pose().residue(resnum).is_upper_terminus();
	bool is_lower_term = pose().residue(resnum).is_lower_terminus();
	pose_nonconst().replace_residue( resnum, res_in, orient_bb );
	pose_nonconst().energies().clear();
	pose_nonconst().conformation().update_polymeric_connection( resnum-1 );
	pose_nonconst().conformation().update_polymeric_connection( resnum );

	if ( !is_lower_term ) remove_lower_terminus_variant_type( resnum );
	if ( !is_upper_term ) remove_upper_terminus_variant_type( resnum );
}

/// @brief aligns the upper-terminal residue of segment1 to the start "anchor" residue of segment2
void
StructureDataWithPose::align_segments(
	std::string const & segment1,
	std::string const & segment2 )
{
	core::Size const align_res_target( segment(segment1).upper() );
	core::Size const align_res_movable( segment(segment2).start() );
	core::Size const align_torsions = align_res_movable;
	int const jump_idx = find_jump( segment2 );
	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	TR.Debug << "perm=" << *this << std::endl;
	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief aligns the lower-terminal residue of segment1 to the end "anchor" residue of segment2
void
StructureDataWithPose::align_segments_rev(
	std::string const & segment1,
	std::string const & segment2 )
{
	core::Size const align_res_movable( segment(segment1).upper() );
	core::Size const align_res_target( segment(segment2).start() );
	core::Size const align_torsions = align_res_target;
	int const jump_idx = find_jump( segment2 );
	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	TR.Debug << "perm=" << *this << std::endl;
	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief aligns the template_res'th residue of template_seg to the movable_res'th residue of movable_seg
/// re-arranges fold tree if necessary
void
StructureDataWithPose::align_residues(
	std::string const & template_seg,
	core::Size const template_res,
	std::string const & movable_seg,
	core::Size const movable_res )
{
	debug_assert( template_res );
	debug_assert( template_res <= segment(template_seg).length() );
	debug_assert( movable_res );
	debug_assert( movable_res <= segment(movable_seg).length() );
	core::Size const align_res_target = segment(template_seg).segment_to_pose(template_res);
	core::Size const align_res_movable = segment(movable_seg).segment_to_pose(movable_res);
	core::Size const align_torsions = 0;

	TR.Debug << "aligning " << align_res_movable << " to " << align_res_target << " with torsions from " << align_torsions << std::endl;
	utility::vector1< std::string > roots;
	roots.push_back( template_seg );
	roots.push_back( movable_seg );
	consolidate_movable_groups( roots );
	int jump_idx = find_jump( movable_seg );
	debug_assert( jump_idx );

	align_residues( jump_idx, align_res_target, align_res_movable, align_torsions );
}

/// @brief aligns two residues so their backbones are completely superimposed
void
StructureDataWithPose::align_residues(
	core::Size const jump_idx,
	core::Size const align_res_target,
	core::Size const align_res_movable,
	core::Size const res_with_torsions )
{
	// move the jump temporarily
	core::kinematics::FoldTree const ft_backup = pose().fold_tree();
	pose_nonconst().fold_tree( slide_jump( ft_backup, jump_idx, align_res_target, align_res_movable ) );

	//setup torsions
	if ( res_with_torsions ) {
		set_phi( align_res_target, pose().phi(res_with_torsions) );
		set_psi( align_res_target, pose().psi(res_with_torsions) );
		set_phi( align_res_movable, pose().phi(res_with_torsions) );
		set_psi( align_res_movable, pose().psi(res_with_torsions) );
	}

	core::id::StubID stub1(
		core::id::AtomID(pose().residue(align_res_target).type().atom_index("CA"), align_res_target),
		core::id::AtomID(pose().residue(align_res_target).type().atom_index("N"), align_res_target),
		core::id::AtomID(pose().residue(align_res_target).type().atom_index("C"), align_res_target) );

	core::id::StubID stub2(
		core::id::AtomID(pose().residue(align_res_movable).type().atom_index("CA"), align_res_movable),
		core::id::AtomID(pose().residue(align_res_movable).type().atom_index("N"), align_res_movable),
		core::id::AtomID(pose().residue(align_res_movable).type().atom_index("C"), align_res_movable) );

	pose_nonconst().conformation().set_stub_transform( stub1, stub2, core::kinematics::RT() );
	debug_assert( ft_backup.check_fold_tree() );
	pose_nonconst().fold_tree( ft_backup );
}

/// @brief moves jump pointing to the first segment so that it will move with the second segment
void
StructureDataWithPose::slide_jump(
	std::string const & child_segment,
	std::string const & parent_segment )
{
	int const childjump = find_jump(child_segment);
	debug_assert( childjump );
	int const parentjump = find_jump(parent_segment);
	if ( parentjump == childjump ) {
		TR.Warning << "Parent jump and child jump are the same while sliding jump for " << child_segment << " to have parent " << parent_segment << " --  not doing anything." << std::endl;
		return;
	}
	core::kinematics::Edge const jedge = pose().fold_tree().jump_edge(childjump);
	pose_nonconst().fold_tree( slide_jump( pose().fold_tree(), childjump, segment(parent_segment).safe(), jedge.stop() ) );
}

core::kinematics::FoldTree
StructureDataWithPose::slide_jump(
	core::kinematics::FoldTree const & ft_orig,
	core::Size const jump_idx,
	core::Size const new_start,
	core::Size const new_stop ) const
{
	core::kinematics::FoldTree ft;
	TR << "Sliding jump " << jump_idx << " to " << new_start << " --> " << new_stop << " in " << ft_orig << std::endl;
	debug_assert( jump_idx > 0 );
	debug_assert( jump_idx <= ft_orig.num_jump() );
	core::kinematics::Edge jedge = ft_orig.jump_edge( jump_idx );
	core::kinematics::FoldTree::EdgeList::const_iterator begin_edge = ft_orig.begin();
	debug_assert( begin_edge != ft_orig.end() );

	if ( ft_orig.root() == jedge.start() ) {
		for ( core::kinematics::FoldTree::EdgeList::const_iterator e=ft_orig.begin(); e!=ft_orig.end(); ++e ) {
			if ( ( e->start() == ft_orig.root() ) && ( *e != jedge ) &&
					( static_cast< core::Size >( e->stop() ) != new_start ) &&
					( static_cast< core::Size >( e->stop() ) != new_stop ) ) {
				begin_edge = e;
				break;
			}
		}
	}
	TR << "Selected new begin edge: " << *begin_edge << std::endl;

	ft.add_edge( *begin_edge );
	for ( core::kinematics::FoldTree::EdgeList::const_iterator e=ft_orig.begin(); e!=ft_orig.end(); ++e ) {
		if ( (e != begin_edge) && (*e != jedge) ) {
			ft.add_edge( *e );
		}
	}
	ft.add_edge( jedge );
	TR.Debug << "Sliding jump " << jump_idx << " to " << new_start << " --> " << new_stop << " in " << ft << std::endl;
	ft.slide_jump( jump_idx, new_start, new_stop );
	debug_assert( ft.check_fold_tree() );
	return ft;
}

/// @brief re-arranges residues such that segment 2 follows segment 1 in sequence
/// @details copies segment 2 to immediately after segment 1 and adds a jump
/// original segment 2 is then deleted.
void
StructureDataWithPose::move_segment(
	std::string const & segment1_n,
	std::string const & segment1_c,
	std::string const & segment2_n,
	std::string const & segment2_c )
{
	// do nothing if the segments are already in the correct order
	if ( segment(segment1_c).upper()+1 == segment(segment2_n).lower() ) {
		return;
	}

	core::Size s1 = segment(segment1_n).lower();
	core::Size e1 = segment(segment1_c).upper();
	core::Size s2 = segment(segment2_n).lower();
	core::Size e2 = segment(segment2_c).upper();
	TR << "moving segment " << segment2_n << "(" << s2 << ") <--> " << segment2_c << "(" << e2 << ")" << " to after " << segment1_c << "(" << e1 << ")" << std::endl;
	TR << "Before moving SD: " << *this << std::endl;
	try {
		check_consistency();
	} catch ( EXCN_PoseInconsistent const & e ) {
		TR.Error << "Input to move_segment() is not consistent with pose -- dumping to predump.pdb" << std::endl;
		pose().dump_pdb( "predump.pdb" );
		throw e;
	}

	move_jumps_to_safety();
	move_segment_in_pose( s1, e1, s2, e2 );

	StructureData::move_segments( segment1_c, segment2_n, segment2_c );

	try {
		check_consistency();
	} catch ( EXCN_PoseInconsistent const & e ) {
		TR.Error << "SD inconsistent with pose. Dumping to dump.pdb. SD=" << *this << std::endl;
		pose().dump_pdb( "dump.pdb" );
		throw e;
	}
	move_jumps_to_safety();
}

void
StructureDataWithPose::add_segment(
	std::string const & id_val,
	Segment const & resis,
	StringList::iterator insert_pos,
	core::pose::Pose const & residues )
{
	StructureData::add_segment( id_val, resis, insert_pos );

	core::Size insert_resid = pose().size() + 1;
	if ( insert_pos != segments_end_nonconst() ) {
		insert_resid = segment( *insert_pos ).lower();
	}
	TR << "Inserting new residues before position " << insert_resid << std::endl;

	pose_nonconst().conformation().insert_conformation_by_jump(
		residues.conformation(),
		insert_resid,
		pose_->num_jump()+2,
		1,
		pose_->num_jump()+1 );
}

/// @brief removes all traces of the given segment from the object
void
StructureDataWithPose::delete_segment( std::string const & seg_val )
{
	TR.Debug << "delete_segment(" << seg_val << ")" << std::endl;
	Segment const & seg = segment( seg_val );

	// collect pose information
	core::Size start_del = seg.lower();
	core::Size stop_del = seg.upper();
	std::string const segment1_c = seg.lower_segment();
	std::string const segment2_n = seg.upper_segment();

	if ( !seg.has_free_lower_terminus() ) {
		if ( start_del <= stop_del ) {
			start_del += 1;
		}
	}

	if ( !seg.has_free_upper_terminus() ) {
		if ( start_del <= stop_del ) {
			debug_assert( stop_del > 1 );
			stop_del -= 1;
		}
	}

	// delete the segment from StructureData
	StructureData::delete_segment( seg_val );

	// delete the segment from the pose
	if ( start_del <= stop_del ) {
		core::Size const resi_count = pose().size();
		delete_residues_in_pose( start_del, stop_del );
		debug_assert( resi_count - (stop_del-start_del+1) == pose().size() );
	}

	if ( !segment1_c.empty() ) {
		add_upper_terminus_variant_type( segment( segment1_c ).upper() );
	}
	if ( !segment2_n.empty() ) {
		add_lower_terminus_variant_type( segment( segment2_n ).lower() );
	}
	chains_from_termini();
	utility::vector1< std::string > roots;
	if ( !segment1_c.empty() ) roots.push_back( segment1_c );
	if ( !segment2_n.empty() ) roots.push_back( segment2_n );
	if ( roots.empty() ) {
		move_jumps_to_safety();
	} else {
		consolidate_movable_groups( roots );
	}
	TR << "Deleted segment " << seg_val << std::endl;
}

/// @brief moves a segment of the pose such that the segment from start2<=res<=end2 is moved so that it starts at end1+1
/// returns the jump number that is to be ignored in fold tree searches
void
StructureDataWithPose::move_segment_in_pose(
	core::Size start1,
	core::Size end1,
	core::Size start2,
	core::Size end2 )
{
	TR << "moving segment " << start2 << "->" << end2 << " to after residues " << start1 << "->" << end1 << std::endl;

	debug_assert( pose_ );
	int const len = end2 - start2 + 1;

	insert_after_residue_in_pose( start1, end1, start2, end2 );
	debug_assert( len >= 0 );
	if ( end1 < start2 ) {
		start2 += len;
		end2 += len;
	}
	delete_residues_in_pose( start2, end2 );
}

/// @brief copies and inserts a segment into the permutation after the given segment
void
StructureDataWithPose::insert_after_residue_in_pose(
	core::Size segment1_start,
	core::Size segment1_end,
	core::Size segment2_start,
	core::Size segment2_end )
{
	debug_assert( pose_ );
	TR.Debug << "Inserting segment -- original segment2: " << segment2_start << "-->" << segment2_end << ", segment1_end=" << segment1_end << std::endl;

	// add segment2 residues
	for ( core::Size cseg_res=segment2_start; cseg_res<=segment2_end; ++cseg_res ) {
		TR.Debug << "copying " << cseg_res << " to " << segment1_end+1 << std::endl;
		// new residue will be attached by jump
		if ( cseg_res != segment2_start ) {
			if ( pose().residue(cseg_res).is_lower_terminus() ) {
				pose_nonconst().insert_residue_by_jump( pose().residue( cseg_res ), segment1_end+1, pose().fold_tree().root() );
			} else {
				//remove_lower_terminus_variant_type( cseg_res );
				pose_nonconst().append_polymer_residue_after_seqpos( pose().residue(cseg_res), segment1_end, false );
			}
		} else {
			int upstream_res = pose().fold_tree().root();
			TR.Debug << "Adding jump from " << upstream_res << " after " << segment1_end << " to residue copied from " << cseg_res << std::endl;
			core::kinematics::FoldTree const & ft = pose().fold_tree();
			for ( core::Size i=1; i<=pose().size(); ++i ) {
				TR.Debug << "Res " << i << " : " << ft.is_cutpoint(i) << std::endl;
			}
			pose_nonconst().insert_residue_by_jump( pose().residue(cseg_res), segment1_end+1, upstream_res );
			add_lower_terminus_variant_type( segment1_end + 1 );
		}
		if ( segment1_end < segment2_start ) {
			cseg_res += 1;
			segment2_start += 1;
			segment2_end += 1;
		}
		segment1_end += 1;
	}
	TR.Debug << "after adding, segment1_end=" << segment1_end << " segment2_start=" << segment2_start << " segment2_end=" << segment2_end << std::endl;

	// fix the fold tree -- jumps without label "jump" which are in c2 need to be reset
	core::kinematics::FoldTree ft( pose().fold_tree() );
	TR.Debug << "Initial foldtree: " << ft << std::endl;
	utility::vector1< core::kinematics::Edge > new_edges;
	for ( int j=1; j<=static_cast<int>(ft.num_jump()); ++j ) {
		int const u = ft.upstream_jump_residue(j);
		int const d = ft.downstream_jump_residue(j);
		// these are important for the casting below
		debug_assert( u >= 0 );
		debug_assert( d >= 0 );
		TR.Debug << "Looking at jump " << j << " upstream " << u << " downstream " << d << " Comp2=" << segment2_start << "->" << segment2_end << std::endl;
		bool const u_in_segment2 = ( segment2_start <= static_cast< core::Size >(u) ) && ( static_cast< core::Size >(u) <= segment2_end );
		bool const d_in_segment1 = ( segment1_start <= static_cast< core::Size >(d) ) && ( static_cast< core::Size >(d) <= segment1_end );
		core::kinematics::Edge & new_edge = ft.jump_edge(j);
		debug_assert( new_edge.start() == u );
		debug_assert( new_edge.stop() == d );
		// 4 cases
		if ( u_in_segment2 && d_in_segment1 ) {
			// 1. upstream in cut segment, downstream in new segment
			TR.Debug << "reversing jump " << j << " from newup=" << d << " to " << " newdown=" << u << std::endl;
			new_edge.start() = d;
			new_edge.stop() = u;
		} else if ( u_in_segment2 && (!d_in_segment1) ) {
			// 2. upstream in cut segment, downstream not in new segment
			core::Size const target_res = segment1_end - ( segment2_end - ft.upstream_jump_residue(j) );
			TR.Debug << "sliding up jump " << j << " from " << u << " to " << target_res << " down=" << d << std::endl;
			new_edge.start() = target_res;
			insert_peptide_edges( ft, new_edge );
		} else {
			// 4. Neither residue is involved in this addition -- don't touch
			TR.Debug << "ignoring jump " << j << " up=" << u << " down=" << d << std::endl;
		}
	}

	TR.Debug << "fold tree after insert : " << ft << std::endl;
	ft.delete_extra_vertices();
	TR.Debug << "cleaned fold tree      : " << ft << std::endl;
	debug_assert( ft.check_fold_tree() );
	pose_nonconst().fold_tree( ft );
}

/// @brief deletes the given residues from the pose
void
StructureDataWithPose::delete_residues_in_pose(
	core::Size const start,
	core::Size const end )
{
	debug_assert( pose_ );
	pose_->conformation().buffer_signals();
	TR.Debug << "Deleting pose residues " << start << " to " << end << std::endl;
	if ( start > end ) {
		std::stringstream ss;
		ss << id() << ": deleting residue from " << start << " to " << end
			<< ", pose length = " << pose_->size() << std::endl;
		ss << "Bad start and end given to delete_residues_in_pose()" << std::endl;
		throw utility::excn::EXCN_BadInput( ss.str() );
	}
	pose_->conformation().delete_residue_range_slow( start, end );
	pose_->conformation().unblock_signals();
}

/// @brief updates numbering based on the saved order of Segment objects
void
StructureDataWithPose::update_numbering()
{
	StructureData::update_numbering();
	update_covalent_bonds_in_pose();
	save_into_pose();
}

} //protocols
} //denovo_design
} //components

