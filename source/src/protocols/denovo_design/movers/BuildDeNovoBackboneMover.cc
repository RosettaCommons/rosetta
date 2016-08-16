// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/movers/BuildDeNovoBackboneMover.cc
/// @brief Mover that builds and folds a structure via fragment insertion
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMover.hh>
#include <protocols/denovo_design/movers/BuildDeNovoBackboneMoverCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/AddConstraints.hh>
#include <protocols/constraint_generator/CoordinateConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh>
#include <protocols/denovo_design/movers/SealFoldTreeMover.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/filters/CalculatorFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

// TEMPORARY
#include <core/pose/datacache/CacheableObserverType.hh>
#include <core/pose/datacache/ObserverCache.hh>
static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.BuildDeNovoBackboneMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

BuildDeNovoBackboneMover::BuildDeNovoBackboneMover():
	protocols::moves::Mover( BuildDeNovoBackboneMover::class_name() ),
	architect_(),
	builder_( new components::ExtendedPoseBuilder ),
	folder_( new components::RemodelLoopMoverPoseFolder ),
	score_filter_(),
	id_(),
	dry_run_( false ),
	connection_overlap_( 1 )
{
}

BuildDeNovoBackboneMover::~BuildDeNovoBackboneMover()
{}

void
BuildDeNovoBackboneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & ,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	set_id( tag->getOption< std::string >( "name" ) );

	core::Size idx = 1;
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() == "PreFoldMovers" ) parse_prefold_movers( *subtag, movers );
		else if ( (*subtag)->getName() == "PostFoldMovers" ) parse_postfold_movers( *subtag, movers );
		else if ( idx == 1 ) {
			parse_architect( *subtag, data );
			++idx;
		} else if ( idx == 2 ) {
			parse_folder( *subtag, data );
			++idx;
		} else {
			std::stringstream msg;
			msg << get_name() << "::parse_my_tag(): Invalid subtag found:" << **subtag << std::endl;
			msg << "Valid subtags are \"PreFoldMovers\", \"PostFoldMovers\", a DeNovo architect, and a Pose Folder."
				<< std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
	}
	TR << "Finished parsing tag: " << prefold_movers_.size() << " prefold movers found." << std::endl;

}

protocols::moves::MoverOP
BuildDeNovoBackboneMover::clone() const
{
	return protocols::moves::MoverOP( new BuildDeNovoBackboneMover( *this ) );
}


protocols::moves::MoverOP
BuildDeNovoBackboneMover::fresh_instance() const
{
	return protocols::moves::MoverOP( new BuildDeNovoBackboneMover );
}

std::string
BuildDeNovoBackboneMover::get_name() const
{
	return BuildDeNovoBackboneMover::class_name();
}

std::string
BuildDeNovoBackboneMover::class_name()
{
	return "BuildDeNovoBackboneMover";
}

std::string const &
BuildDeNovoBackboneMover::id() const
{
	return id_;
}

void
BuildDeNovoBackboneMover::set_id( std::string const & id_value )
{
	id_ = id_value;
}

protocols::filters::FilterOP
BuildDeNovoBackboneMover::default_score_filter()
{
	protocols::filters::CalculatorFilterOP calc_filt( new protocols::filters::CalculatorFilter( "-SS" ) );
	protocols::filters::FilterOP ss_filt( new protocols::fldsgn::filters::SecondaryStructureFilter );
	calc_filt->add_filter( "SS", ss_filt );
	return calc_filt;
}

void
BuildDeNovoBackboneMover::apply( core::pose::Pose & pose )
{
	if ( !architect_ ) {
		std::stringstream msg;
		msg << get_name() << "::apply(): No architect is set!"
			<< get_name() << " requires an architect to build and modify poses." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	components::StructureDataOP sd = architect_->apply( pose );

	core::pose::PoseOP built = builder_->apply( *sd );
	sd->check_pose_consistency( *built );
	pose = *built;
	//pose.dump_pdb( "prefold_initial.pdb" );

	components::StructureDataFactory::get_instance()->save_into_pose( pose, *sd );

	try {
		fold( pose );
		set_last_move_status( protocols::moves::MS_SUCCESS );
	} catch ( components::EXCN_Fold const & e ) {
		TR << "Folding failed: ";
		e.show( TR );
		TR << std::endl;
		static std::string const fold_fail_pdb = "fold_failed.pdb";
		TR << "PDB dumped to " << fold_fail_pdb << std::endl;
		pose.dump_pdb( fold_fail_pdb );
		set_last_move_status( protocols::moves::FAIL_RETRY );
	} catch ( EXCN_NothingToFold const & e ) {
		e.show( TR );
		TR.flush();
		TR.Warning << "WARNING: not folding anything!" << std::endl;
		set_last_move_status( protocols::moves::MS_SUCCESS );
	}
}

/// @brief handles acceptance/reject
void
BuildDeNovoBackboneMover::check_and_accept(
	core::pose::Pose const & orig,
	core::pose::Pose & pose,
	FoldScore & bestscore ) const
{
	if ( !score_filter_ ) {
		TR << "Accepting because no score filter is set" << std::endl;
		return;
	}

	FoldScore const myscore = folding_score( pose );
	if ( bestscore < myscore ) {
		// reject
		TR << "New pose (" << myscore << ") is worse than the old one ("
			<< bestscore << ") in score... reverting." << std::endl;
		pose = orig;
	} else {
		// accept
		TR << "Accepting pose with score " << myscore << ", previous best " << bestscore << std::endl;
		bestscore = myscore;
	}
}

SegmentNames
BuildDeNovoBackboneMover::segments_to_fold( components::StructureData const & sd ) const
{
	SegmentNames segs;
	for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		if ( !sd.segment( *s ).template_pose() ) segs.push_back( *s );
	}
	return segs;
}

SegmentNames
BuildDeNovoBackboneMover::find_roots(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops ) const
{
	SegmentNames roots;

	components::StructureData const & sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	for ( protocols::loops::Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
		core::Size const startres =
			loop_start_without_overlap( pose, l->start(), connection_overlap_ );
		std::string const & start_comp = sd.segment_name( startres );
		std::string const & prevcomp = sd.segment( start_comp ).lower_segment();
		if ( prevcomp != "" ) {
			roots.push_back( prevcomp );
		}
		core::Size const stopres =
			loop_stop_without_overlap( pose, l->stop(), connection_overlap_ );
		std::string const & stop_comp = sd.segment_name( stopres );
		std::string const & nextcomp = sd.segment(stop_comp).upper_segment();
		if ( nextcomp != "" ) {
			roots.push_back( nextcomp );
		}
		TR.Debug << id() << ": Loop: " << startres << " " << prevcomp << "__" << start_comp << " --> " << stopres << " " << stop_comp << "__" << nextcomp << std::endl;
	}
	if ( roots.empty() ) {
		roots.push_back( *sd.segments_begin() );
	}
	return roots;
}

void
BuildDeNovoBackboneMover::parse_prefold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & moversmap )
{
	MoverOPs const & movers = parse_movers( tag, moversmap );
	for ( MoverOPs::const_iterator m=movers.begin(); m!=movers.end(); ++m ) {
		prefold_movers_.push_back( *m );
	}
}

void
BuildDeNovoBackboneMover::parse_postfold_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & moversmap )
{
	MoverOPs const & movers = parse_movers( tag, moversmap );
	for ( MoverOPs::const_iterator m=movers.begin(); m!=movers.end(); ++m ) {
		postfold_movers_.push_back( *m );
	}
}

BuildDeNovoBackboneMover::MoverOPs
BuildDeNovoBackboneMover::parse_movers( utility::tag::TagCOP tag, protocols::moves::Movers_map const & moversmap ) const
{
	MoverOPs movers;
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() != "Add" ) {
			std::stringstream msg;
			msg << type() << ": Invalid xml tag name (" << (*subtag)->getName() << ") found in tag " << *tag << std::endl;
			msg << "Valid tags are: \"Add\"." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		std::string const mover_name = (*subtag)->getOption< std::string >( "mover" );
		protocols::moves::Movers_map::const_iterator find_mover = moversmap.find( mover_name );
		if ( find_mover == moversmap.end() ) {
			std::stringstream msg;
			msg << type() << "::parse_movers(): ERROR !! mover not found in map: \n" << **subtag << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		movers.push_back( find_mover->second->clone() );
		TR.Debug << "found mover " << mover_name << std::endl;
	}
	return movers;
}

BuildDeNovoBackboneMover::MoverOPs
BuildDeNovoBackboneMover::prefold_movers(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	ConstraintGeneratorCOPs const & generators ) const
{
	MoverOPs movers;

	// 1. switch to centroid
	protocols::moves::MoverOP switch_cen( new protocols::simple_moves::SwitchResidueTypeSetMover( "centroid" ) );
	movers.push_back( switch_cen );

	// 2. Run user-provided movers
	for ( MoverOPs::const_iterator m=prefold_movers_.begin(); m!=prefold_movers_.end(); ++m ) {
		movers.push_back( *m );
	}

	// 3. Add coordinate constraints for overlapping regions
	if ( !generators.empty() ) {
		protocols::moves::MoverOP add_csts( new protocols::constraint_generator::AddConstraints( generators ) );
		if ( add_csts ) movers.push_back( add_csts );
	}

	// 4. Build proper fold tree for remodel
	SegmentNames const roots = find_roots( pose, loops );
	protocols::moves::MoverOP set_ft( new FoldTreeFromFoldGraphMover( roots, loops ) );
	movers.push_back( set_ft );

	return movers;
}

BuildDeNovoBackboneMover::MoverOPs
BuildDeNovoBackboneMover::postfold_movers(
	core::pose::Pose const & pose,
	protocols::loops::Loops const & loops,
	ConstraintGeneratorCOPs const & generators ) const
{
	MoverOPs movers;

	components::StructureData const sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	// 1. Remove coodinate constraints for overlapping regions
	//    if they were added
	if ( !generators.empty() ) {
		protocols::moves::MoverOP rm_csts( new protocols::constraint_generator::RemoveConstraints( generators ) );
		if ( rm_csts ) movers.push_back( rm_csts );
	}

	// 2. Run user-provided movers
	for ( MoverOPs::const_iterator m=postfold_movers_.begin(); m!=postfold_movers_.end(); ++m ) {
		movers.push_back( *m );
	}

	// 3. Switch to full atom
	protocols::moves::MoverOP switch_fa( new protocols::simple_moves::SwitchResidueTypeSetMover( "fa_standard" ) );
	movers.push_back( switch_fa );

	// 4. Replace sidechains
	// TL: SaveAndRetriveSidechains copies the pose, which gets rid of the observer
	// TODO: fix this in Pose, this is horrible behavior!
	protocols::moves::MoverOP restore_sc_mover(
		new simple_moves::ReturnSidechainMover( pose, 1, pose.total_residue() ) );
	//protocols::moves::MoverOP sars(
	// new simple_moves::SaveAndRetrieveSidechains(
	// pose,
	// true, // allsc -- replace all sidechains, as opposed to only ALAs
	// false, // ensure_variant_matching
	// 0 ) ); // jumpid -- no idea why this needs to be set
	movers.push_back( restore_sc_mover );

	// 5. Seal fold tree
	protocols::moves::MoverOP set_ft( new SealFoldTreeMover( sd, loops ) );
	movers.push_back( set_ft );

	return movers;
}

void
BuildDeNovoBackboneMover::clear_prefold_movers()
{
	prefold_movers_.clear();
}

void
BuildDeNovoBackboneMover::add_prefold_mover( protocols::moves::Mover const & prefold_mover )
{
	prefold_movers_.push_back( prefold_mover.clone() );
}

void
BuildDeNovoBackboneMover::clear_postfold_movers()
{
	postfold_movers_.clear();
}

void
BuildDeNovoBackboneMover::add_postfold_mover( protocols::moves::Mover const & postfold_mover )
{
	postfold_movers_.push_back( postfold_mover.clone() );
}

protocols::loops::Loops
BuildDeNovoBackboneMover::create_loops(
	core::pose::Pose const & pose,
	ResidueVector & overlap_residues,
	SegmentNames & movable_segments ) const
{
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	components::StructureData const & sd = factory.get_from_const_pose( pose );
	components::FoldGraph fg( sd );
	movable_segments = segments_to_fold( sd );
	protocols::loops::LoopsOP loops = fg.create_loops( movable_segments );
	if ( !loops ) {
		std::stringstream msg;
		msg << type() << "::create_loops(): Loops object could not be created! Movable segments="
			<< movable_segments << " SD=" << sd << std::endl;
		throw EXCN_NothingToFold( msg.str() );
	}

	// add requested connection overlap
	overlap_residues = add_overlap_to_loops( *loops, connection_overlap_, pose );
	return *loops;
}

void
BuildDeNovoBackboneMover::fold( core::pose::Pose & pose ) const
{
	TR << "Attempting to fold structure" << std::endl;

	// make Loops object to pass to PoseFolder
	ResidueVector overlap_residues;
	SegmentNames movable_segments;
	protocols::loops::Loops const loops = create_loops( pose, overlap_residues, movable_segments );
	if ( loops.empty() ) {
		TR.Warning << "Loops object is empty.  Not doing anything." << std::endl;
		return;
	}

	// Make movable residue subset for PoseFolder
	core::select::residue_selector::ResidueSubset const movable =
		select_movable_residues( pose, movable_segments, overlap_residues );

	// Create pre-fold and post-fold movers
	ConstraintGeneratorCOPs const generators = create_constraint_generators( overlap_residues );

	// Get StructureData object from pose
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	components::StructureData const sd = factory.get_from_pose( pose );

	// apply user-provided movers
	MoverOPs const prefold = prefold_movers( pose, loops, generators );
	TR << "Before postfold pose is_attached: " << pose.observer_cache().is_attached( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER ) << std::endl;
	MoverOPs const postfold = postfold_movers( pose, loops, generators );

	// TL: RestoreSidechainsMover copies the pose, which nukes my observer. Re-add it here
	// TODO: fix this in pose
	factory.save_into_pose( pose, sd );

	apply_movers( prefold, pose );

	// initial scoring
	core::pose::Pose const orig( pose );

	// TL: Just created a copy of pose, so need to readd the observer
	// TODO: fix in pose
	factory.save_into_pose( pose, sd );

	TR << "After copying pose is_attached: " << pose.observer_cache().is_attached( core::pose::datacache::CacheableObserverType::STRUCTUREDATA_OBSERVER ) << std::endl;
	FoldScore bestscore = folding_score( pose );

	// Store copy of SD with cutpoints removed
	components::StructureData newsd = factory.get_from_pose( pose );
	remove_cutpoints( newsd, loops );

	// fold entire structure first
	TR << id() << ": Folding entire chain" << std::endl;
	folder_->apply( pose, movable, loops );
	factory.save_into_pose( pose, newsd );
	check_and_accept( orig, pose, bestscore );

	// apply user-provided movers
	apply_movers( postfold, pose );
}

void
BuildDeNovoBackboneMover::remove_cutpoints( components::StructureData & sd, protocols::loops::Loops const & loops ) const
{
	for ( protocols::loops::Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
		if ( ! l->cut() ) continue;
		TR.Debug << "Removing cutpoint from StructureData: " << l->cut() << std::endl;
		sd.set_cutpoint( sd.segment_name( l->cut() ), 0 );
	}
}

void
BuildDeNovoBackboneMover::apply_movers( MoverOPs const & movers, core::pose::Pose & pose ) const
{
	//pose.dump_pdb( "prefold_0.pdb" );
	core::Size count = 1;
	for ( MoverOPs::const_iterator m=movers.begin(); m!=movers.end(); ++m, ++count ) {
		TR.Debug << "Running mover " << (*m)->get_name() << std::endl;
		(*m)->apply( pose );
		//pose.dump_pdb( "prefold_" + (*m)->get_name() + boost::lexical_cast< std::string >( count ) + ".pdb" );
	}
}

core::select::residue_selector::ResidueSubset
BuildDeNovoBackboneMover::select_movable_residues(
	core::pose::Pose const & pose,
	SegmentNames const & movable_segments,
	ResidueVector const & overlap_residues ) const
{
	components::StructureData const & sd =
		components::StructureDataFactory::get_instance()->get_from_const_pose( pose );

	core::select::residue_selector::ResidueSubset subset( sd.pose_length(), false );
	for ( SegmentNames::const_iterator s=movable_segments.begin(); s!=movable_segments.end(); ++s ) {
		for ( core::Size resid=sd.segment( *s ).lower(); resid<=sd.segment( *s ).upper(); ++resid ) {
			subset[ resid ] = true;
		}
	}
	for ( ResidueVector::const_iterator r=overlap_residues.begin(); r!=overlap_residues.end(); ++r ) {
		debug_assert( *r > 0 );
		debug_assert( *r <= pose.total_residue() );
		subset[ *r ] = true;
	}
	return subset;
}

BuildDeNovoBackboneMover::ConstraintGeneratorCOPs
BuildDeNovoBackboneMover::create_constraint_generators( ResidueVector const & residues ) const
{
	if ( residues.empty() ) {
		return BuildDeNovoBackboneMover::ConstraintGeneratorCOPs();
	}
	return boost::assign::list_of (create_coordinate_constraint_generator( residues ));
}

BuildDeNovoBackboneMover::ConstraintGeneratorOP
BuildDeNovoBackboneMover::create_coordinate_constraint_generator( ResidueVector const & residues ) const
{
	protocols::constraint_generator::CoordinateConstraintGeneratorOP coord_csts(
		new protocols::constraint_generator::CoordinateConstraintGenerator );

	coord_csts->set_id( "BuildDeNovoBackboneMover_" + boost::lexical_cast< std::string >( numeric::random::rg().uniform() ) );

	std::stringstream index_ss;
	for ( ResidueVector::const_iterator r=residues.begin(); r!=residues.end(); ++r ) {
		if ( r != residues.begin() ) index_ss << ',';
		index_ss << *r;
	}

	core::select::residue_selector::ResidueIndexSelectorOP selector(
		new core::select::residue_selector::ResidueIndexSelector( index_ss.str() ) );

	coord_csts->set_ca_only( true );
	coord_csts->set_residue_selector( selector );

	return coord_csts;
}

BuildDeNovoBackboneMover::FoldScore
BuildDeNovoBackboneMover::folding_score( core::pose::Pose const & pose ) const
{
	if ( !score_filter_ ) return 0.0;

	TR << "Computing score" << std::endl;
	core::Real const score = score_filter_->report_sm( pose );
	TR << "Score is " << score << std::endl;
	return score;
}

void
BuildDeNovoBackboneMover::set_score_filter( protocols::filters::Filter const & filter )
{
	score_filter_ = filter.clone();
}

void
BuildDeNovoBackboneMover::parse_architect( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	architect_ = architects::DeNovoArchitectFactory::get_instance()->create_from_tag( tag, data );
}

void
BuildDeNovoBackboneMover::parse_folder( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	components::RemodelLoopMoverPoseFolderOP folder( new components::RemodelLoopMoverPoseFolder );
	folder->parse_my_tag( tag, data );
	folder_ = folder;
}

void
BuildDeNovoBackboneMover::set_architect( architects::DeNovoArchitect const & architect )
{
	architect_ = architect.clone();
}

void
BuildDeNovoBackboneMover::set_folder( components::PoseFolder const & folder )
{
	folder_ = folder.clone();
}

void
BuildDeNovoBackboneMover::set_connection_overlap( core::Size const overlap_val )
{
	connection_overlap_ = overlap_val;
}

/////////////// Creator ///////////////

protocols::moves::MoverOP
BuildDeNovoBackboneMoverCreator::create_mover() const
{
	return protocols::moves::MoverOP( new BuildDeNovoBackboneMover );
}

std::string
BuildDeNovoBackboneMoverCreator::keyname() const
{
	return BuildDeNovoBackboneMover::class_name();
}

////////////// Helper functions /////////////////

/// @brief goes through loops and adds overlapping residues
BuildDeNovoBackboneMover::ResidueVector
add_overlap_to_loops(
	protocols::loops::Loops & loops,
	core::Size const overlap,
	core::pose::Pose const & pose )
{
	BuildDeNovoBackboneMover::ResidueVector residues;
	if ( !overlap ) return residues;

	// go through loops and add proper overlap
	for ( protocols::loops::Loops::iterator l=loops.v_begin(); l!=loops.v_end(); ++l ) {
		for ( core::Size i=1; i<=overlap; ++i ) {
			// see if this loop is lower-terminal
			if ( core::pose::is_lower_terminus( pose, l->start() ) ) {
				break;
			}
			l->set_start( l->start()-1 );
			residues.push_back( l->start() );
		}
		for ( core::Size i=1; i<=overlap; ++i ) {
			// check if this loop is upper-terminal
			if ( core::pose::is_upper_terminus( pose, l->stop() ) ) {
				break;
			}
			l->set_stop( l->stop()+1 );
			residues.push_back( l->stop() );
		}
		TR << "LOOP is now " << l->start() << " " << l->stop() << " " << l->cut() << std::endl;
	}
	return residues;
}

} //protocols
} //denovo_design
} //movers

