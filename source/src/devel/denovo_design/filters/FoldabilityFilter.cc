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

// Protocol Headers
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/GrowLeft.hh>
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/toolbox/SelectResiduesByLayer.hh>

// Core Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pack/task/residue_selector/ResidueSelector.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/sequence/ABEGOManager.hh>

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


static thread_local basic::Tracer TR( "devel.denovo_design.filters.FoldabilityFilter" );

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
	start_res_( 1 ),
	end_res_( 1 ),
	distance_threshold_( 4.0 ),
	ignore_pose_abego_( false ),
	output_poses_( false ),
	selector_(),
	vlb_(	protocols::forge::components::VarLengthBuildOP( new protocols::forge::components::VarLengthBuild() ) ),
	cached_aa_( "" ),
	cached_ss_( "" ),
	cached_start_( 0 ),
	cached_end_( 0 )
{
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
			selector_ = data.get_ptr< core::pack::task::residue_selector::ResidueSelector const >( "ResidueSelector", selectorname );
		} catch ( utility::excn::EXCN_Msg_Exception e ) {
			std::stringstream error_msg;
			error_msg << "Failed to find ResidueSelector named '" << selectorname << "' from the Datamap from DisulfidizeMover.\n";
			error_msg << e.msg();
			throw utility::excn::EXCN_Msg_Exception( error_msg.str() );
		}
		assert( selector_ );
		TR << "Using residue selector " << selectorname << std::endl;
	} else {
		set_start_res( tag->getOption< core::Size >( "start_res", start_res_ ) );
		set_end_res( tag->getOption< core::Size >( "end_res", end_res_ ) );
	}
	motif_ = tag->getOption< std::string >( "motif", motif_ );
	ignore_pose_abego_ = tag->getOption< bool >( "ignore_pose_abego", ignore_pose_abego_ );
	output_poses_ = tag->getOption< bool >( "output_poses", output_poses_ );
	distance_threshold_ = tag->getOption< core::Real >( "distance_threshold", distance_threshold_ );
	if ( ignore_pose_abego_ && ( motif_ == "" ) ) {
		utility_exit_with_message( "You need to specify a motif if you are ignoring pose abego values." );
	}
}

void
FoldabilityFilter::set_start_res( core::Size const startval )
{
	start_res_ = startval;
}

void
FoldabilityFilter::set_end_res( core::Size const endval )
{
	end_res_ = endval;
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
	// work on pose copy
	core::pose::PoseOP posecopy = generate_pose( pose );

	core::Size start( 0 );
	core::Size end( 0 );
	choose_start_and_end( start, end, *posecopy );
	runtime_assert( end >= start );

	// if motif isn't specified, we will need dssp information
	if ( motif_ == "" ) {
		protocols::moves::DsspMover dssp;
		dssp.apply( *posecopy );
		TR << "pose dssp is : " << pose.secstruct() << std::endl;
	}

	std::string ss( "" );
	std::string aa( "" );
	utility::vector1< std::string > abego;
	get_aa_ss_abego( aa, ss, abego, start, end, *posecopy );

	// save end residue for later comparision
	core::conformation::Residue const end_res( posecopy->residue( end ) );

	// eliminate the segment to be tested.
	delete_segment( *posecopy, start, end );

	// safety, clear the energies object
	posecopy->energies().clear();

	setup_vlb( aa, ss, abego, start, end );
	runtime_assert( vlb_ );

	TR << "ss to build starting with residue " << start << " is " << ss << " with aa " << aa << std::endl;
	core::Size const good_count = fragment_insertion( *posecopy, end, end_res );

	return (core::Real)good_count / (core::Real)tries_;
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

/// @brief determine start and end window, from selector if necessary
void
FoldabilityFilter::choose_start_and_end(
		core::Size & start,
		core::Size & end,
		core::pose::Pose const & pose ) const
{
	start = 0;
	end = 0;
	if ( selector_ ) {
		core::pack::task::residue_selector::ResidueSubset subset = selector_->apply(pose);
		bool finished = false;
		for ( core::Size i=1, endi=subset.size(); i<=endi; ++i ) {
			if ( subset[i] ) {
				if ( finished ) {
					throw utility::excn::EXCN_BadInput( "ResidueSelector does not select one continuous set of residues. FoldabilityFilter requires a continuous segment to be specified." );
				}
				if ( !start )
					start = i;
				end = i;
			} else {
				if ( start && end ) {
					finished = true;
				}
			}
		}
		TR << "Residue selector chose residues " << start << " --> " << end << std::endl;
	} else {
		start = start_res_;
		end = end_res_;
		TR << "Start and end residues are: " << start << " --> " << end << std::endl;
	}
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
		for ( core::Size i=start; i<=end; ++i ) {
			ss += pose.secstruct(i);
			aa += pose.residue(i).name1();
		}
	} else {
		// parse motif string
		core::Size insert_length( 0 );
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

/// @brief deletes the segment from start to end (inclusive) from the pose
void
FoldabilityFilter::delete_segment(
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
	pose.conformation().delete_residue_range_slow( start, end );
}

/// @brief performs setup on the fragment insertion machinery
void
FoldabilityFilter::setup_vlb(
		std::string const & aa,
		std::string const & ss,
		utility::vector1< std::string > const & abego,
		core::Size const start,
		core::Size const end ) const
{
	// only set the manager if something has changed... otherwise it will pick new fragments each time the mover is called
	if ( ( ss != cached_ss_ ) || ( aa != cached_aa_ ) || ( start != cached_start_ ) || ( end != cached_end_ ) ) {
		// the build manager
		protocols::forge::build::BuildManager manager;
		// add the instruction to rebuild this segment
		if ( start == 1 ) {
			manager.add( protocols::forge::build::BuildInstructionOP( new protocols::forge::build::GrowLeft( start, ss, aa ) ) );
		} else {
			manager.add( protocols::forge::build::BuildInstructionOP( new protocols::forge::build::GrowRight( start, ss, aa ) ) );
		}
		// clear fragment cache and set buildmanager
		vlb_->manager( manager );
		// set cache
		cached_ss_ = ss;
		cached_aa_ = aa;
		cached_start_ = start;
		cached_end_ = end;
	}
	//vlb_->scorefunction( scorefxn_->clone() );
	vlb_->vall_memory_usage( protocols::forge::components::VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->loop_mover_str( "RemodelLoopMover" );
	vlb_->set_abego( abego );
}

/// @brief performs fragment insertion and returns number of successful builds. Assumes setup_vlb() has already been called
core::Size
FoldabilityFilter::fragment_insertion(
		core::pose::Pose const & pose,
		core::Size const end,
		core::conformation::Residue const & end_res ) const
{
	debug_assert( vlb_ );
	core::Size good_count = 0;
	for ( core::Size i=1; i<=tries_; ++i ) {
		core::pose::Pose posecopy( pose );
		vlb_->apply( posecopy );
		if ( output_poses_ ) {
			posecopy.dump_pdb( "foldability" + boost::lexical_cast< std::string >( i ) + ".pdb" );
		}
		core::Real const distance = end_res.xyz( "N" ).distance( posecopy.residue( end ).xyz( "N" ) );
		TR << "Distance is " << distance << ", success threshold=" << distance_threshold_ << std::endl;
		if ( distance < distance_threshold_ ) {
			++good_count;
		}
	}
	return good_count;
}

} // namespace filters
} // namespace denovo_design
} // namespace devel

