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
#include <protocols/constraint_generator/AtomPairConstraintGenerator.hh>
#include <protocols/constraint_generator/RemoveConstraints.hh>
#include <protocols/denovo_design/architects/CompoundArchitect.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/DivideAndConqueror.hh>
#include <protocols/denovo_design/components/ExtendedPoseBuilder.hh>
#include <protocols/denovo_design/components/FoldGraph.hh>
#include <protocols/denovo_design/components/RemodelLoopMoverPoseFolder.hh>
#include <protocols/denovo_design/components/RandomTorsionPoseFolder.hh>
#include <protocols/denovo_design/components/NullPoseFolder.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/SegmentPairing.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/StructureDataPerturber.hh>
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>
#include <protocols/denovo_design/movers/FoldTreeFromFoldGraphMover.hh>
#include <protocols/denovo_design/movers/SealFoldTreeMover.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/filters/CalculatorFilter.hh>
#include <protocols/fldsgn/filters/SecondaryStructureFilter.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

// Core headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/select/residue_selector/ResidueIndexSelector.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <numeric/random/random.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

// C++ headers
#include <ctime>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.movers.BuildDeNovoBackboneMover" );

namespace protocols {
namespace denovo_design {
namespace movers {

BuildDeNovoBackboneMover::BuildDeNovoBackboneMover():
	protocols::moves::Mover( BuildDeNovoBackboneMover::class_name() ),
	architect_(),
	builder_( new components::ExtendedPoseBuilder ),
	folder_( new components::RemodelLoopMoverPoseFolder ),
	perturber_( new components::NullPerturber ),
	prefold_movers_(),
	postfold_movers_(),
	filters_(),
	score_filter_(),
	id_(),
	dry_run_( false ),
	dump_pdbs_( false ),
	build_overlap_( 1 ),
	iterations_per_phase_( 0 ),
	start_segments_(),
	stop_segments_()
{
}

BuildDeNovoBackboneMover::~BuildDeNovoBackboneMover()
{}

void
BuildDeNovoBackboneMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const & filters,
	protocols::moves::Movers_map const & movers,
	core::pose::Pose const & )
{
	set_id( tag->getOption< std::string >( "name" ) );

	build_overlap_ = tag->getOption< core::Size >( "build_overlap", build_overlap_ );

	std::string const segments_csv = tag->getOption< std::string >( "start_segments", "" );
	if ( tag->hasOption( "start_segments" ) ) set_start_segments( segments_csv );

	std::string const stop_segments_csv = tag->getOption< std::string >( "stop_segments", "" );
	if ( tag->hasOption( "stop_segments" ) ) set_stop_segments( stop_segments_csv );

	iterations_per_phase_ = tag->getOption< core::Size >( "iterations_per_phase", iterations_per_phase_ );
	dump_pdbs_ = tag->getOption< bool >( "dump_pdbs", dump_pdbs_ );

	core::Size idx = 1;
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() == "PreFoldMovers" ) parse_prefold_movers( *subtag, movers );
		else if ( (*subtag)->getName() == "PostFoldMovers" ) parse_postfold_movers( *subtag, movers );
		else if ( (*subtag)->getName() == "Filters" ) parse_filters( *subtag, filters );
		else if ( (*subtag)->getName() == "Perturbers" ) parse_perturbers( *subtag, data );
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

	// create StructureData "blueprint"
	components::StructureDataOP sd = architect_->apply( pose );
	if ( !sd ) {
		std::stringstream msg;
		msg << "BuildDeNovoBackboneMover::architect " << architect_->id()
			<< " failed to generate anything." << std::endl;
		utility_exit_with_message( msg.str() );
	}

	// Break up into managable pieces
	components::DivideAndConqueror divider;
	divider.set_start_segments( start_segments_ );
	divider.set_stop_segments( stop_segments_ );

	components::BuildPhases const build_phases = divider.divide_and_conquer( *sd );

	TR << "Trying to build the following structure:" << std::endl;
	TR << *sd << std::endl;

	core::pose::PoseOP pose_ptr = build_in_phases( *sd, build_phases, SegmentNameSet(), 1 );
	if ( pose_ptr ) pose = *pose_ptr;
	else set_last_move_status( protocols::moves::FAIL_RETRY );
}

/// @brief given an build phase number, returns a pdb filename for checkpointing
///        Example filename: iter_01_20160718133800.pdb
/// @param[in] phase_num  Build phase number
std::string
pdb_filename( core::Size const phase_num )
{
	std::stringstream filename;
	filename << time( NULL ) << "_iter_" << std::setw( 2 )
		<< std::setfill('0') << phase_num << ".pdb";
	return filename.str();
}

/// @brief builds/folds pose in phases using recursive algorithm
/// @throws EXCN_Fold if we couldn't fold the pose
core::pose::PoseOP
BuildDeNovoBackboneMover::build_in_phases(
	components::StructureData const & full_sd,
	components::BuildPhases const & phases,
	SegmentNameSet const & finished,
	core::Size const phase_num ) const
{
	core::Size const max_iter = iterations_per_phase_ == 0 ? 50 : iterations_per_phase_;

	SegmentNames const & phase = *phases.begin();
	TR << "Building phase " << phase_num << " " << phase << std::endl;

	// compute segments to slice out of StructureData
	SegmentNameSet all_phase_segments( finished.begin(), finished.end() );
	all_phase_segments.insert( phase.begin(), phase.end() );

	components::StructureData working_sd = full_sd;
	components::StructureDataPerturberOP phase_perturber = perturber_->clone();

	// compute segments not in this phase
	SegmentNameSet const in_this_phase( phase.begin(), phase.end() );
	SegmentNameSet const in_sd( full_sd.segments_begin(), full_sd.segments_end() );
	SegmentNames const ignore_vec = set_difference< SegmentNames, SegmentName >( in_sd.begin(), in_sd.end(), in_this_phase.begin(), in_this_phase.end() );
	SegmentNameSet const ignore_set( ignore_vec.begin(), ignore_vec.end() );

	phase_perturber->set_ignore_segments( ignore_set );
	TR << "Set ignored segments to " << ignore_set << std::endl;

	for ( core::Size iter=1; iter<=max_iter; ++iter ) {
		TR << "Attempting to fold phase " << phase_num << " " << phase << " "
			<< " iteration " << iter << " / " << max_iter << std::endl;

		// perturb
		phase_perturber->apply( working_sd );

		// slice
		TR << "Slicing " << all_phase_segments << " out of StructureData" << std::endl;
		components::StructureData const phase_sd = working_sd.slice( all_phase_segments, true );
		TR << "Done slicing. New SD=" << phase_sd << std::endl;

		// build pose
		core::pose::PoseOP built = builder_->apply( phase_sd );
		try {
			fold_attempt( *built );

			/// dump pdb
			if ( dump_pdbs_ ) built->dump_pdb( pdb_filename( phase_num ) );

			/// add template information based on the folded segments, so that the next phase can use it
			/// this will never be reached if check_pose() throws an exception
			components::StructureData new_full_sd = working_sd;
			components::StructureData const & built_sd = components::StructureDataFactory::get_instance()->get_from_const_pose( *built );
			for ( SegmentNameList::const_iterator s=built_sd.segments_begin(); s!=built_sd.segments_end(); ++s ) {
				new_full_sd.set_template_pose( *s, *built, built_sd.segment( *s ).start(), built_sd.segment( *s ).stop() );
			}

			if ( ++phases.begin() == phases.end() ) {
				return built;
			} else {
				return build_in_phases(
					new_full_sd,
					components::BuildPhases( ++phases.begin(), phases.end() ),
					SegmentNameSet( all_phase_segments.begin(), all_phase_segments.end() ),
					phase_num + 1 );
			}
		} catch ( components::EXCN_Fold const & e ) {
			e.show( TR );
			TR << std::endl << "Failing phase " << phase_num << " " << phase << std::endl;
		} catch ( EXCN_NothingToFold const & e ) {
			e.show( TR );
			TR.flush();
			TR.Warning << "WARNING: not folding anything!" << std::endl;
		} catch ( EXCN_FilterFailed const & e ) {
			e.show( TR );
			TR.flush();
		}
	}

	// if we reach here without returning something, folding failed.
	std::stringstream msg;
	msg << class_name() << "build_in_phases():  Failed to fold anything passing filters after "
		<< max_iter << " attempts." << std::endl;
	throw components::EXCN_Fold( msg.str() );
	return core::pose::PoseOP();
}

/// @brief checks the pose for whether it folded correctly, according to user's filters
void
BuildDeNovoBackboneMover::check_pose( core::pose::Pose const & pose ) const
{
	core::Size filter_num = 1;
	for ( FilterCOPs::const_iterator f=filters_.begin(); f!=filters_.end(); ++f, ++filter_num ) {
		if ( !(*f)->apply( pose ) ) throw EXCN_FilterFailed( (*f)->get_type(), filter_num );
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
			loop_start_without_overlap( pose, l->start(), build_overlap_ );
		std::string const & start_comp = sd.segment_name( startres );
		std::string const & prevcomp = sd.segment( start_comp ).lower_segment();
		if ( prevcomp != "" ) {
			roots.push_back( prevcomp );
		}
		core::Size const stopres =
			loop_stop_without_overlap( pose, l->stop(), build_overlap_ );
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
BuildDeNovoBackboneMover::parse_perturbers( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	if ( tag->getTags().size() > 1 ) {
		std::stringstream msg;
		msg << get_name() << ": Only one perturber can be used at a time. You have specified "
			<< tag->getTags().size() << " : " << *tag << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		perturber_ = components::StructureDataPerturber::create( **subtag, data );
	}
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

void
BuildDeNovoBackboneMover::parse_filters( utility::tag::TagCOP tag, protocols::filters::Filters_map const & filter_map )
{
	filters_.clear();
	for ( utility::tag::Tag::tags_t::const_iterator subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		if ( (*subtag)->getName() != "Add" ) {
			std::stringstream msg;
			msg << type() << ": Invalid xml tag name (" << (*subtag)->getName() << ") found in tag " << *tag << std::endl;
			msg << "Valid tags are: \"Add\"." << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		std::string const filter_name = (*subtag)->getOption< std::string >( "filter" );
		protocols::filters::Filters_map::const_iterator find_filter = filter_map.find( filter_name );
		if ( find_filter == filter_map.end() ) {
			std::stringstream msg;
			msg << type() << "::parse_filters(): ERROR !! filter not found in map: \n" << **subtag << std::endl;
			throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
		}
		filters_.push_back( find_filter->second->clone() );
		TR.Debug << "found filter " << filter_name << std::endl;
	}
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

	// 2. Set secondary structure in pose
	SetPoseSecstructFromStructureDataMover set_ss;
	movers.push_back( set_ss.clone() );

	// 3. Run user-provided movers
	for ( MoverOPs::const_iterator m=prefold_movers_.begin(); m!=prefold_movers_.end(); ++m ) {
		movers.push_back( *m );
	}

	// 4. Add coordinate constraints for overlapping regions
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

	// 2. Seal fold tree
	protocols::moves::MoverOP seal_ft( new SealFoldTreeMover( sd, loops ) );
	movers.push_back( seal_ft );

	// 2. Run user-provided movers
	for ( MoverOPs::const_iterator m=postfold_movers_.begin(); m!=postfold_movers_.end(); ++m ) {
		movers.push_back( *m );
	}

	// 3. Switch to full atom
	protocols::moves::MoverOP switch_fa( new protocols::simple_moves::SwitchResidueTypeSetMover( "fa_standard" ) );
	movers.push_back( switch_fa );

	// 4. Replace sidechains
	protocols::moves::MoverOP restore_sc_mover(
		new simple_moves::ReturnSidechainMover( pose, 1, pose.size() ) );
	//protocols::moves::MoverOP sars(
	// new simple_moves::SaveAndRetrieveSidechains(
	// pose,
	// true, // allsc -- replace all sidechains, as opposed to only ALAs
	// false, // ensure_variant_matching
	// 0 ) ); // jumpid -- no idea why this needs to be set
	movers.push_back( restore_sc_mover );

	// 5. Re-seal fold tree for full-atom mode
	//    This is necessary if there are sidechain-sidechain covalent bonds
	movers.push_back( seal_ft );

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
	overlap_residues = add_overlap_to_loops( *loops, build_overlap_, pose );
	return *loops;
}

core::scoring::ScoreFunctionOP
get_score_function( core::pose::Pose const & pose )
{
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::pose::symmetry::make_score_function_consistent_with_symmetric_state_of_pose( pose, scorefxn );
	return scorefxn;
}

void
BuildDeNovoBackboneMover::fold_attempt( core::pose::Pose & pose ) const
{
	// Get StructureData object from pose
	components::StructureDataFactory const & factory = *components::StructureDataFactory::get_instance();
	if ( !factory.has_cached_data( pose ) ) {
		utility_exit_with_message( "No StructureData at start of fold attempt" );
	}

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

	components::StructureData const sd = factory.get_from_pose( pose );

	// apply user-provided movers
	MoverOPs const prefold = prefold_movers( pose, loops, generators );
	MoverOPs const postfold = postfold_movers( pose, loops, generators );

	// apply user-provided movers
	if ( !factory.has_cached_data( pose ) ) {
		utility_exit_with_message( "No data before applying premovers" );
	}
	apply_movers( prefold, pose );

	// fold structure
	// folders need to keep sd and pose in sync
	folder_->apply( pose, movable, loops );

	try {
		factory.get_from_pose( pose ).check_pose_consistency( pose );
	} catch ( components::EXCN_PoseInconsistent const & e ) {
		std::stringstream msg;
		msg << "BuildDeNovoBackboneMover::fold_attempt(): After pose folder, StructureData no longer agrees with pose."
			<< " PoseFolders need to ensure that StructureData and pose are synced.  Error =" << std::endl;
		e.show( msg );
		utility_exit_with_message( msg.str() );
	}

// apply user-provided movers
	apply_movers( postfold, pose );

	// modify structuredata so this portion of the structure can be checked independently
	components::StructureData const saved_sd = factory.get_from_pose( pose );
	//components::StructureData check_sd = saved_sd;
	//modify_for_check( check_sd );
	//TR.Debug << "SD for check = " << check_sd << std::endl;
	//factory.save_into_pose( pose, check_sd );

	// save modified secstruct into pose for checking
	//SetPoseSecstructFromStructureDataMover().apply( pose );

	// set Energies object -- sets hbonds info in the Energies object in the pose
	pose.energies().clear();
	(*get_score_function(pose))( pose );

	// check pose against user-provided filters
	check_pose( pose );

	// save new sd into pose
	//TR.Debug << "Saving into pose: " << saved_sd << std::endl;
	//components::StructureDataFactory::get_instance()->save_into_pose( pose, saved_sd );
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
		debug_assert( *r <= pose.size() );
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
	return boost::assign::list_of (create_overlap_constraint_generator( residues ));
}

/// @brief builds a constraint generator for residues that overlap between phases
/// @param[in] residues Vector of resids corresponding to the overlap residues
BuildDeNovoBackboneMover::ConstraintGeneratorOP
BuildDeNovoBackboneMover::create_overlap_constraint_generator( ResidueVector const & residues ) const
{
	using core::select::residue_selector::ResidueIndexSelector;

	protocols::constraint_generator::AtomPairConstraintGeneratorOP dist_csts(
		new protocols::constraint_generator::AtomPairConstraintGenerator );

	dist_csts->set_id( "BuildDeNovoBackboneMover_constrain_overlap_residues" );
	dist_csts->set_ca_only( false );
	dist_csts->set_min_seq_sep( 0 );
	dist_csts->set_sd( 1.0 );

	std::stringstream index_ss;
	for ( ResidueVector::const_iterator r=residues.begin(); r!=residues.end(); ++r ) {
		if ( r != residues.begin() ) index_ss << ',';
		index_ss << *r;
	}

	TR << "Creating index selector from " << index_ss.str() << std::endl;
	ResidueIndexSelector const selector( index_ss.str() );
	dist_csts->set_residue_selector( selector );
	return dist_csts;
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
	components::PoseFolderOP folder;
	if ( tag->getName() == components::RemodelLoopMoverPoseFolder::class_name() ) {
		folder = components::PoseFolderOP( new components::RemodelLoopMoverPoseFolder );
	} else if ( tag->getName() == components::RandomTorsionPoseFolder::class_name() ) {
		folder = components::PoseFolderOP( new components::RandomTorsionPoseFolder );
	} else if ( tag->getName() == components::NullPoseFolder::class_name() ) {
		folder = components::PoseFolderOP( new components::NullPoseFolder );
	} else {
		std::stringstream msg;
		msg << "BuildDeNovoBackboneMover::parse_folder(): Unknown folder type: "
			<< tag->getName() << std::endl;
		utility_exit_with_message( msg.str() );
	}

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
BuildDeNovoBackboneMover::set_build_overlap( core::Size const overlap_val )
{
	build_overlap_ = overlap_val;
}

/// @brief sets names of segments to be included in the starting build phase
/// @param[in] segments_csv Comma-separated string containing segment names
void
BuildDeNovoBackboneMover::set_start_segments( std::string const & segments_csv )
{
	set_start_segments( csv_to_container< SegmentNameSet >( segments_csv ) );
}

/// @brief sets names of segments to be included in the starting build phase
/// @param[in] segments Set of segment names
void
BuildDeNovoBackboneMover::set_start_segments( SegmentNameSet const & segments )
{
	start_segments_ = segments;
}

/// @brief sets names of segments to be included in the final build phase
/// @param[in] segments_csv Comma-separated string containing segment names
void
BuildDeNovoBackboneMover::set_stop_segments( std::string const & segments_csv )
{
	set_stop_segments( csv_to_container< SegmentNameSet >( segments_csv ) );
}

/// @brief sets names of segments to be included in the final build phase
/// @param[in] segments Set of segment names
void
BuildDeNovoBackboneMover::set_stop_segments( SegmentNameSet const & segments )
{
	stop_segments_ = segments;
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

////////////// Helper movers /////////////////
std::string
SetPoseSecstructFromStructureDataMover::get_name() const
{
	return class_name();
}

protocols::moves::MoverOP
SetPoseSecstructFromStructureDataMover::clone() const
{
	return protocols::moves::MoverOP( new SetPoseSecstructFromStructureDataMover(*this) );
}

void
SetPoseSecstructFromStructureDataMover::apply( core::pose::Pose & pose )
{
	components::StructureData const & sd = components::StructureDataFactory::get_instance()->get_from_const_pose( pose );
	std::string target_ss = "";

	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		target_ss = symmetric_secstruct( pose, sd.ss() );
	} else {
		target_ss = sd.ss();
	}

	if ( target_ss.size() != pose.size() ) {
		std::stringstream msg;
		msg << class_name() << "::apply(): StructureData ss size ("
			<< target_ss.size() << ") does not match pose size ("
			<< pose.size() << ") " << std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size resid = 1;
	for ( std::string::const_iterator sschar=target_ss.begin(); sschar!=target_ss.end(); ++sschar, ++resid ) {
		pose.set_secstruct( resid, *sschar );
	}
	TR << "Set Pose secstruct to " << target_ss << std::endl;

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

architects::RegisterShift
compute_nobu_register_shift(
	components::StructureData const & sd,
	components::SegmentPair const & pair,
	architects::StrandOrientation const & o1,
	architects::StrandOrientation const & o2,
	architects::RegisterShift const & shift )
{
	bool const order_reversed = ( sd.segment( pair.first ).safe() > sd.segment( pair.second ).safe() );
	architects::StrandOrientation const ref_orientation = order_reversed ? o2 : o1;

	if ( ref_orientation == architects::UP ) {
		if ( order_reversed ) {
			return -shift;
		} else {
			return shift;
		}
	} else if ( ref_orientation == architects::DOWN ) {
		architects::RegisterShift const inverted_shift =
			sd.segment( pair.first ).elem_length() - ( sd.segment( pair.second ).elem_length() + shift );
		if ( order_reversed ) {
			return -inverted_shift;
		} else {
			return inverted_shift;
		}
	} else {
		std::stringstream msg;
		msg << "compute_nobu_register_shift(): Invalid orientation for strand named " << pair.first << " : " << o1 << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return 0;
}

protocols::fldsgn::topology::StrandPairingSet
compute_strand_pairings( components::StructureData const & sd, components::SegmentPairSet const & pairs )
{
	using protocols::fldsgn::topology::StrandPairing;
	using protocols::fldsgn::topology::StrandPairingOP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::SS_Info2;

	SS_Info2 const ss_info( sd.ss() );
	StrandPairingSet pairset;

	for ( components::SegmentPairSet::const_iterator pair=pairs.begin(); pair!=pairs.end(); ++pair ) {
		if ( !sd.has_segment( pair->first ) ) continue;
		if ( !sd.has_segment( pair->second ) ) continue;
		components::Segment const & seg1 = sd.segment( pair->first );
		components::Segment const & seg2 = sd.segment( pair->second );
		architects::StrandOrientation const o1 =
			architects::StrandArchitect::int_to_orientation( sd.get_data_int( pair->first, architects::StrandArchitect::orientation_keyname() ) );
		architects::StrandOrientation const o2 =
			architects::StrandArchitect::int_to_orientation( sd.get_data_int( pair->second, architects::StrandArchitect::orientation_keyname() ) );
		architects::RegisterShift const tomp_shift = sd.get_data_int( pair->second, architects::StrandArchitect::register_shift_keyname() );

		core::Size const length = seg1.elem_length() <= seg2.elem_length() ? seg1.elem_length() : seg2.elem_length();
		architects::RegisterShift const shift = compute_nobu_register_shift( sd, *pair, o1, o2, tomp_shift );
		char const orient = ( o1 == o2 ) ? 'P' : 'A';
		StrandPairingOP sp( new StrandPairing(
			ss_info.strand_id( seg1.safe() ),  // Strand 1
			ss_info.strand_id( seg2.safe() ),  // Strand 2
			seg1.start(),                     // Strand 1 start res
			seg2.start(),                     // Strand 2 start res
			length,                        // pairing length
			shift,                        // register shift
			orient ) );                      // orientation
		TR << "Created strand pairing for *pair " << *sp << std::endl;
		pairset.push_back( sp );
	}
	pairset.finalize();
	return pairset;
}

architects::RegisterShift
retrieve_shift( components::StructureData const & sd, SegmentName const & segment_name )
{
	return sd.get_data_int( segment_name, architects::StrandArchitect::register_shift_keyname() );
}

architects::StrandOrientation
retrieve_orientation( components::StructureData const & sd, SegmentName const & segment_name )
{
	return architects::StrandArchitect::int_to_orientation(
		sd.get_data_int( segment_name, architects::StrandArchitect::orientation_keyname() ) );
}

/*
/// @brief returns sorted vector of bulge residues
core::select::residue_selector::ResidueVector
get_bulges( components::StructureData const & sd )
{
core::select::residue_selector::ResidueVector bulges;
std::string::const_iterator ss = sd.ss().begin();
std::string::const_iterator ab = sd.abego().begin();
for ( core::Size resid=1; resid<=sd.pose_length(); ++resid, ++ss, ++ab ) {
if ( ( *ss == 'E' ) && ( *ab == 'A' ) ) bulges.push_back( resid );
}
return bulges;
}

void
modify_for_check( components::StructureData & sd )
{
using core::select::residue_selector::ResidueVector;
using core::select::residue_selector::ResidueSubset;
using protocols::fldsgn::topology::StrandPairingSet;
using components::ResiduePairs;

ResidueSubset paired( sd.pose_length(), false );

ResiduePairs const pairs = components::SegmentPairing::get_strand_residue_pairs( sd );
for ( ResiduePairs::const_iterator p=pairs.begin(); p!=pairs.end(); ++p ) {
paired[p->first] = true;
paired[p->second] = true;
}

ResidueVector const bulges = get_bulges( sd );
TR.Debug << "Bulges = " << bulges << std::endl;
for ( ResidueVector::const_iterator r=bulges.begin(); r!=bulges.end(); ++r ) {
paired[*r] = true;
}
TR.Debug << "Pair set = " << paired << std::endl;

core::Size resid = 1;
for ( ResidueSubset::const_iterator p=paired.begin(); p!=paired.end(); ++p, ++resid ) {
if ( !*p && ( sd.ss()[ resid - 1 ] == 'E' ) ) {
TR << "Setting ss for residue " << resid << " to L" << std::endl;
sd.set_ss( resid, 'L' );
}
}
}
*/

} //protocols
} //denovo_design
} //movers

