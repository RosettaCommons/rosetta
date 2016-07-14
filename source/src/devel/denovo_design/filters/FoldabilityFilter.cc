// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
/// @file src/devel/denovo_design/filters/FoldabilityFilter.cc
/// @brief Tom's Denovo design protocol
/// @details
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <devel/denovo_design/filters/FoldabilityFilter.hh>
#include <devel/denovo_design/filters/FoldabilityFilterCreator.hh>

// Project Headers
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/components/Picker.hh>

// Protocol Headers
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/GrowLeft.hh>
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/forge/remodel/RemodelLoopMover.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/rosetta_scripts/util.hh>

// Core Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/select/residue_selector/SecondaryStructureSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/pose/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/util/SwitchResidueTypeSet.hh>

// Basic Headers
#include <basic/datacache/DataMap.hh>
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// ObjexxFCL Headers

// C++ Headers


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif


static THREAD_LOCAL basic::Tracer TR( "devel.denovo_design.filters.FoldabilityFilter" );

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace devel {
namespace denovo_design {
namespace filters {

std::string
FoldabilityFilterCreator::keyname() const
{
	return FoldabilityFilterCreator::filter_name();
}

protocols::filters::FilterOP
FoldabilityFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new FoldabilityFilter() );
}

std::string
FoldabilityFilterCreator::filter_name()
{
	return "Foldability";
}

/*
NOTES
test case is a long helix, 12 residues fixed helical
plus 9 residues variable for fragment insertion
plus 12 residues fixed helical

Architecture: Javier writes code to make a fragment insertion move
I write code that evaluates the move
Possible metrics for study :
independent variable fragment length and abego
angle of H-axis bend
dihedral about H-axis
could go with 6-D translation/rotation IF we are changing insert lengths.
if we are not, project the a parametric helix forward or back to find the closest point
to a predefined end point in euclidean space.
*/

///  ---------------------------------------------------------------------------------
///  FoldabilityFilter main code:
///  ---------------------------------------------------------------------------------
FoldabilityFilter::FoldabilityFilter() :
	Filter( "FoldabilityFilter" ),
	motif_( "" ),
	tries_( 100 ),
	distance_threshold_( 4.0 ),
	ignore_pose_abego_( false ),
	use_sequence_( false ),
	output_poses_( false ),
	scorefxn_(),
	selector_(),
	picker_( PickerOP( new Picker() ) ),
	vlb_( protocols::forge::components::VarLengthBuildOP( new protocols::forge::components::VarLengthBuild() ) )
{
	segments_.clear();
}

/// @brief destructor - this class has no dynamic allocation, so
//// nothing needs to be cleaned. C++ will take care of that for us.
FoldabilityFilter::~FoldabilityFilter()
{}


/// Return a copy of ourselves
protocols::filters::FilterOP
FoldabilityFilter::clone() const
{
	return protocols::filters::FilterOP( new FoldabilityFilter(*this) );
}

protocols::filters::FilterOP
FoldabilityFilter::fresh_instance() const
{
	return protocols::filters::FilterOP( new FoldabilityFilter() );
}

void
FoldabilityFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	tries_ = tag->getOption< core::Size >( "tries", tries_ );
	if ( tag->hasOption("selector") ) {
		if ( tag->hasOption("start_res") || tag->hasOption("end_res" ) ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "\"start_res\" and \"end_res\" cannot be used at the same time as \"selector\" in FoldabilityFilter." );
		}
		std::string const selectorname = tag->getOption< std::string >("selector");
		try {
			selector_ = data.get_ptr< core::select::residue_selector::ResidueSelector const >( "ResidueSelector", selectorname );
		} catch ( utility::excn::EXCN_Msg_Exception & e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' from the Datamap from DisulfidizeMover.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		assert( selector_ );
		TR << "Using residue selector " << selectorname << std::endl;
	} else {
		segments_.clear();
		add_segment( tag->getOption< core::Size >( "start_res", 1 ),
			tag->getOption< core::Size >( "end_res", 1 ) );
	}
	motif_ = tag->getOption< std::string >( "motif", motif_ );
	ignore_pose_abego_ = tag->getOption< bool >( "ignore_pose_abego", ignore_pose_abego_ );
	use_sequence_ = tag->getOption< bool >( "use_sequence", use_sequence_ );
	output_poses_ = tag->getOption< bool >( "output_poses", output_poses_ );
	distance_threshold_ = tag->getOption< core::Real >( "distance_threshold", distance_threshold_ );
	if ( ignore_pose_abego_ && ( motif_ == "" ) ) {
		utility_exit_with_message( "You need to specify a motif if you are ignoring pose abego values." );
	}
	if ( tag->hasOption( "scorefxn" ) ) {
		scorefxn_ = protocols::rosetta_scripts::parse_score_function( tag, data )->clone();
		TR << "score function " << tag->getOption< std::string >( "scorefxn" ) << " is used. " << std::endl;
	}
}

void
FoldabilityFilter::add_segment( core::Size const startval, core::Size const endval )
{
	segments_.push_back( core::select::residue_selector::ResidueRange( startval, endval ) );
}

std::string
FoldabilityFilter::get_name() const
{
	return "Foldability";
}

void
FoldabilityFilter::report( std::ostream & out, core::pose::Pose const & ) const
{
	out << " reporting" << std::endl;
}

core::Real
FoldabilityFilter::report_sm( core::pose::Pose const & pose ) const
{
	return compute( pose );
}

core::Real
FoldabilityFilter::compute( core::pose::Pose const & pose ) const
{
	using core::select::residue_selector::ResidueRanges;

	ResidueRanges segments;
	if ( selector_ ) {
		segments = ResidueRanges( selector_->apply( pose ) );
	} else {
		segments = segments_;
	}

	core::Real score = 0.0;
	for ( ResidueRanges::const_iterator range=segments.begin(); range!=segments.end(); ++range ) {
		score += compute_segment( pose, *range );
	}
	return score / core::Real( segments.size() );
}

core::Real
FoldabilityFilter::compute_segment(
	core::pose::Pose const & pose,
	core::select::residue_selector::ResidueRange const & segment ) const
{
	using protocols::denovo_design::components::StructureDataOP;
	using protocols::denovo_design::components::StructureDataFactory;

	// work on pose copy
	core::pose::PoseOP posecopy = generate_pose( pose );

	core::Size start = segment.start();
	core::Size end = segment.stop();
	runtime_assert( end >= start );

	if ( StructureDataFactory::get_instance()->has_cached_string( *posecopy ) ) {
		StructureDataOP perm = StructureDataFactory::get_instance()->create_from_pose( *posecopy, "foldability" );
		std::string const ss = perm->ss();
		core::Size res = 1;
		for ( std::string::const_iterator c=ss.begin(), endc=ss.end(); c!=endc; ++c, ++res ) {
			posecopy->set_secstruct( res, *c );
		}
		TR << "Pose SS set from permutation info: " << posecopy->secstruct() << std::endl;
	} else {
		protocols::moves::DsspMover dssp;
		dssp.apply( *posecopy );
		TR << "pose dssp is : " << posecopy->secstruct() << std::endl;
	}

	std::string ss( "" );
	std::string aa( "" );
	utility::vector1< std::string > abego;
	get_aa_ss_abego( aa, ss, abego, start, end, *posecopy );

	// save end residue for later comparision
	core::conformation::Residue const end_res( posecopy->residue( end ) );

	// get the pose ready for fragment insertion
	prepare_pose( *posecopy, start, end );

	// safety, clear the energies object
	posecopy->energies().clear();

	TR << "ss to build starting with residue " << start << " is " << ss.substr(start-1,end-start+1) << " with aa " << aa.substr(start-1,end-start+1) << std::endl;

	protocols::moves::MoverOP insert_fragments;
	if ( use_sequence_ ) {
		insert_fragments = create_fragment_insertion_mover( aa, ss, abego, posecopy->conformation().chain_endings(), start, end );
	} else {
		insert_fragments = create_fragment_insertion_mover( "", ss, abego, posecopy->conformation().chain_endings(), start, end );
	}

	runtime_assert( insert_fragments );

	core::Size const good_count = fragment_insertion( *posecopy, *insert_fragments, end, end_res );

	return core::Real( good_count ) / core::Real( tries_ );
}

/// @brief Does the Foldability Filtering
bool
FoldabilityFilter::apply( core::pose::Pose const & pose ) const
{
	report_sm( pose );
	return true;
}

/// @brief gets non-const version the pose for the filter to work on
core::pose::PoseOP
FoldabilityFilter::generate_pose( core::pose::Pose const & pose ) const
{
	core::pose::PoseOP posecopy;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		posecopy = core::pose::PoseOP( new core::pose::Pose );
		core::pose::symmetry::extract_asymmetric_unit(pose, *posecopy);
		for ( core::Size i=1, endi=posecopy->total_residue(); i<=endi; ++i ) {
			if ( posecopy->residue_type(i).name() == "VRT" ) {
				posecopy->conformation().delete_residue_slow( posecopy->total_residue() );
			}
		}
	} else {
		posecopy = pose.clone();
	}
	debug_assert( posecopy );
	return posecopy;
}

/// @brief gets aa string, ss string, and abego vector for the area to rebuild
void
FoldabilityFilter::get_aa_ss_abego(
	std::string & aa,
	std::string & ss,
	utility::vector1< std::string > & abego,
	core::Size const start,
	core::Size & end,
	core::pose::Pose const & pose ) const
{
	utility::vector1< std::string > abego_insert;
	// if a motif is not specified, determine it from the structure
	if ( motif_ == "" ) {
		for ( core::Size i=1, endi=pose.total_residue(); i<=endi; ++i ) {
			if ( i == end+1 ) {
				continue;
			}
			ss += pose.secstruct(i);
			aa += pose.residue(i).name1();
		}
	} else {
		// parse motif string
		core::Size insert_length( 0 );
		for ( core::Size i=1, endi=start-1; i<=endi; ++i ) {
			ss += pose.secstruct(i);
			aa += pose.residue(i).name1();
		}
		utility::vector1< std::string > const & motifs( utility::string_split( motif_, '-' ) );
		for ( core::Size i=1; i<=motifs.size(); ++i ) {
			std::string const ss_type( motifs[i].substr( motifs[i].size()-2, 1 ) );
			std::string const abego_type( motifs[i].substr( motifs[i].size()-1, 1 ) );
			core::Size const len( (core::Size)utility::string2int( motifs[i].substr( 0, motifs[i].size()-2 ) ) );
			TR << "motif" << i << " = " << motifs[i] << " " << ss_type << " " << abego_type << " " << len << std::endl;
			insert_length += len;
			for ( core::Size j=1; j<=len; ++j ) {
				ss += ss_type;
				aa += "V";
				abego_insert.push_back( abego_type );
			}
		}
		for ( core::Size i=end+2, endi=pose.total_residue(); i<=endi; ++i ) {
			ss += pose.secstruct(i);
			aa += pose.residue(i).name1();
		}
		// there is a chance that the motif is a different length -- we need to update end
		end = insert_length + start;
	}

	// create full abego vector
	utility::vector1< std::string > abego_all( core::sequence::get_abego( pose, 2 ) );
	TR << "abego from full input pose is " << abego_all << std::endl;
	if ( ignore_pose_abego_ ) {
		runtime_assert( abego_insert.size() );
		for ( core::Size i=1; i<start; ++i ) {
			abego.push_back( abego_all[i] );
		}
		for ( core::Size i=1, endi=abego_insert.size(); i<=endi; ++i ) {
			abego.push_back( abego_insert[i] );
		}
		for ( core::Size i=end+2, endi=abego.size(); i<=endi; ++i ) {
			abego.push_back( abego_all[i] );
		}
		TR << "abego from motif is: " << abego << std::endl;
	} else {
		for ( core::Size i=1, endi=end; i<=endi; ++i ) {
			abego.push_back( abego_all[i] );
		}
		for ( core::Size i=end+2, endi=abego_all.size(); i<=endi; ++i ) {
			abego.push_back( abego_all[i] );
		}
		TR << "abego from pose is: " << abego << std::endl;
	}
}

/// @brief prepares the pose/segment from start to end for insertion
void
FoldabilityFilter::prepare_pose(
	core::pose::Pose & pose,
	core::Size const start,
	core::Size const end ) const
{
	// cut the residue after end to avoid clashes
	if ( end+1 <= pose.total_residue() ) {
		pose.conformation().delete_residue_slow( end+1 );
	}
	// insert a jump from start-1 to end+1 (which is actually residue end+2 from the original)
	core::Size jumpend = 0;
	if ( end+1 <= pose.total_residue() ) {
		jumpend = end+1;
	}
	core::Size jumpstart = 0;
	if ( start > 1 ) {
		jumpstart = start - 1;
	}
	if ( jumpend && jumpstart ) {
		core::kinematics::FoldTree ft = pose.fold_tree();
		ft.new_jump( jumpstart, jumpend, end );
		debug_assert( ft.check_fold_tree() );
		pose.fold_tree(ft);
	}
	core::pose::add_upper_terminus_type_to_pose_residue( pose, end );
	core::pose::add_lower_terminus_type_to_pose_residue( pose, end+1 );

	// switch to centroid
	if ( !pose.is_centroid() ) {
		core::util::switch_to_residue_type_set( pose, "centroid" );
	}
}

/// @brief performs fragment picking/other preparations for building
protocols::moves::MoverOP
FoldabilityFilter::create_fragment_insertion_mover(
	std::string const & complete_aa,
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	utility::vector1< core::Size > const & chain_endings,
	core::Size const start,
	core::Size const end ) const
{
	debug_assert( picker_ );

	// create loops
	protocols::loops::LoopsOP loops( new protocols::loops::Loops() );
	loops->push_back( protocols::loops::Loop( start, end, start ) );

	core::fragment::FragSetOP frag9 = picker_->pick_and_cache_fragments( complete_aa, complete_ss, complete_abego, chain_endings, start, end, 9 );
	core::fragment::FragSetOP frag3 = picker_->pick_and_cache_fragments( complete_aa, complete_ss, complete_abego, chain_endings, start, end, 3 );

	protocols::forge::remodel::RemodelLoopMoverOP remodel( new protocols::forge::remodel::RemodelLoopMover( loops ) );
	assert( remodel );
	remodel->add_fragments( frag9 );
	remodel->add_fragments( frag3 );

	// setup movemap for remodel so that only the target loops can move
	core::kinematics::MoveMap mm;
	// initialize explicitly to all false
	for ( core::Size res=1, endr=complete_abego.size(); res<=endr; ++res ) {
		mm.set_bb( res, false );
		mm.set_chi( res, false );
	}
	for ( core::Size l=start; l<=end; ++l ) {
		mm.set_bb( l, true );
		mm.set_chi( l, true );
	}
	remodel->false_movemap( mm );
	remodel->set_keep_input_foldtree( true );
	if ( scorefxn_ ) {
		remodel->scorefunction( *scorefxn_ );
	}

	return remodel;
}

/// @brief performs fragment insertion and returns number of successful builds. Assumes setup_vlb() has already been called
core::Size
FoldabilityFilter::fragment_insertion(
	core::pose::Pose const & pose,
	protocols::moves::Mover & fragment_mover,
	core::Size const end,
	core::conformation::Residue const & end_res ) const
{
	core::Size good_count = 0;
	for ( core::Size i=1; i<=tries_; ++i ) {
		core::pose::Pose posecopy( pose );
		fragment_mover.apply( posecopy );
		if ( output_poses_ ) {
			posecopy.dump_pdb( "foldability" + boost::lexical_cast< std::string >( i ) + ".pdb" );
		}
		core::Real const distance = end_res.xyz( "N" ).distance( posecopy.residue( end ).xyz( "N" ) );
		TR << "Trial " << i << " : distance is " << distance << ", success threshold=" << distance_threshold_ << std::endl;
		if ( distance < distance_threshold_ ) {
			++good_count;
		}
	}
	return good_count;
}

} // namespace filters
} // namespace denovo_design
} // namespace devel

