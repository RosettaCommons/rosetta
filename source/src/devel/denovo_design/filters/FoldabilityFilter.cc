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
/// @detailed
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <devel/denovo_design/filters/FoldabilityFilter.hh>
#include <devel/denovo_design/filters/FoldabilityFilterCreator.hh>

// Protocol Headers

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/sc/ShapeComplementarityCalculator.hh>
#include <core/util/ABEGOManager.hh>
#include <protocols/fldsgn/topology/HelixPairing.hh>
#include <protocols/fldsgn/topology/HSSTriplet.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/forge/build/BuildManager.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/forge/build/GrowRight.hh>
#include <protocols/forge/components/VarLengthBuild.hh>
#include <protocols/jd2/parser/BluePrint.hh>
#include <protocols/moves/DsspMover.hh>
#include <protocols/toolbox/SelectResiduesByLayer.hh>

// Basic Headers
#include <basic/MetricValue.hh>
#include <basic/Tracer.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <utility/string_util.hh>

// ObjexxFCL Headers

//C++ Headers


#ifdef GL_GRAPHICS
#include <protocols/viewer/viewers.hh>
#endif

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

#ifdef BOINC_GRAPHICS
#include <protocols/boinc/boinc.hh>
#endif


static basic::Tracer TR("devel.denovo_design.filters.FoldabilityFilter");

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
	return new FoldabilityFilter();
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
	ignore_pose_abego_( false ),
	vlb_(	new protocols::forge::components::VarLengthBuild() ),
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
	return new FoldabilityFilter(*this);
}

protocols::filters::FilterOP
FoldabilityFilter::fresh_instance() const
{
	return new FoldabilityFilter();
}

void
FoldabilityFilter::parse_my_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
	tries_ = tag->getOption< core::Size >( "tries", tries_ );
	start_res_ = tag->getOption< core::Size >( "start_res", start_res_ );
	end_res_ = tag->getOption< core::Size >( "end_res", end_res_ );
	motif_ = tag->getOption< std::string >( "motif", motif_ );
	ignore_pose_abego_ = tag->getOption< bool >( "ignore_pose_abego", ignore_pose_abego_ );
	if ( ignore_pose_abego_ && ( motif_ == "" ) ) {
		utility_exit_with_message( "You need to specify a motif if you are ignoring pose abego values." );
	}
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
	core::pose::Pose posecopy( pose );
	core::Size good_count( 0 );
	core::Size end( 0 );
	std::string ss( "" );
	std::string aa( "" );
	utility::vector1< std::string > abego_insert;
	// if a motif is not specified, determine it from the structure
	if ( motif_ == "" ) {
		protocols::moves::DsspMover dssp;
		dssp.apply( posecopy );
		for ( core::Size i=start_res_; i<=end_res_; ++i ) {
			ss += posecopy.secstruct( i );
			aa += pose.residue(i).name1();
		}
		end = end_res_;
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
		end = insert_length + start_res_;
	}

	// create abego vector
	utility::vector1< std::string > abego_all( core::util::get_abego( pose, 2 ) );
	TR << "abego from full input pose is " << abego_all << std::endl;
	utility::vector1< std::string > abego;
	if ( ignore_pose_abego_ ) {
		runtime_assert( abego_insert.size() );
		for ( core::Size i=1; i<start_res_; ++i ) {
			abego.push_back( abego_all[i] );
		}
		for ( core::Size i=1; i<=abego_insert.size(); ++i ) {
			abego.push_back( abego_insert[i] );
		}
		TR << "abego from motif is: " << abego << std::endl;
	} else {
		for ( core::Size i=1; i<=end; ++i ) {
			abego.push_back( abego_all[i] );
		}
		TR << "abego from pose is: " << abego << std::endl;
	}

	// eliminate the segment to be tested.
	core::conformation::Residue const end_res( pose.residue( end ) );
	runtime_assert( end >= start_res_ );
	posecopy.conformation().delete_residue_range_slow( start_res_+1, pose.total_residue() );

	// safety, clear the energies object
	posecopy.energies().clear();
	core::pose::Pose saved_pose( posecopy );

	runtime_assert( vlb_ );
	// only set the manager if something has changed... otherwise it will pick new fragments each time the mover is called
	if ( ( ss != cached_ss_ ) || ( aa != cached_aa_ ) || ( start_res_ != cached_start_ ) || ( end != cached_end_ ) ) {
		// the build manager
		protocols::forge::build::BuildManager manager;
		// add the instruction to rebuild this segment
		manager.add( new protocols::forge::build::GrowRight( start_res_, ss, aa ) );
		// clear fragment cache and set buildmanager
		vlb_->manager( manager );
		// set cache
		cached_ss_ = ss;
		cached_aa_ = aa;
		cached_start_ = start_res_;
		cached_end_ = end;
	}
	//vlb_->scorefunction( scorefxn_->clone() );
	vlb_->vall_memory_usage( protocols::forge::components::VLB_VallMemoryUsage::CLEAR_IF_CACHING_FRAGMENTS );
	vlb_->loop_mover_str( "RemodelLoopMover" );
	vlb_->set_abego( abego );

	TR << "ss to build starting with residue " << start_res_ << " is " << ss << " with aa " << aa << std::endl;
	for ( core::Size i=1; i<=tries_; ++i ) {
		vlb_->apply( posecopy );
		posecopy.dump_pdb( "foldability" + boost::lexical_cast< std::string >( i ) + ".pdb" );
		core::Real const distance = end_res.xyz( "N" ).distance( posecopy.residue( end ).xyz( "N" ) );
		TR << "DIstance is " << distance << std::endl;
		if ( distance < 4.0 ) {
			++good_count;
		}
		posecopy = saved_pose;
	}
	return (core::Real)good_count / (core::Real)tries_;
}

/// @brief Does the Foldability Filtering
bool
FoldabilityFilter::apply( core::pose::Pose const & pose ) const
{
	report_sm( pose );
	return true;
}


} // namespace filters
} // namespace denovo_design
} // namespace devel


