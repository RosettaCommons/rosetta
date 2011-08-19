// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/protein_interface_design/filters/RmsdFilter.cc
/// @brief rmsd filtering
/// @author Jacob Corn (jecorn@u.washington.edu)
#include <protocols/protein_interface_design/filters/RmsdFilter.hh>
#include <protocols/protein_interface_design/filters/RmsdFilterCreator.hh>
#include <protocols/filters/Filter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/cacheable_observers.hh>
#include <core/conformation/Conformation.hh>
#include <utility/tag/Tag.hh>
#include <protocols/moves/DataMap.hh>
#include <protocols/moves/Mover.fwd.hh> //Movers_map
#include <core/pose/PDBInfo.hh>

#include <core/scoring/rms_util.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>

#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <algorithm>
#include <list>

//Auto Headers

#include <basic/Tracer.hh>

namespace protocols {
namespace protein_interface_design {
namespace filters {

RmsdFilter::RmsdFilter() :
	protocols::filters::Filter( "Rmsd" ),
	superimpose_( true ),
	symmetry_( false ),
	threshold_( 5.0 ),
	reference_pose_( NULL ),
	selection_from_segment_cache_(false)
{
	selection_.clear();
}

RmsdFilter::RmsdFilter(
	std::list<core::Size> const selection,
	bool const superimpose,
	core::Real const threshold,
	core::pose::PoseOP reference_pose
) : protocols::filters::Filter( "Rmsd" ),
		selection_(selection),
		superimpose_(superimpose),
		threshold_(threshold),
		reference_pose_(reference_pose),
		selection_from_segment_cache_(false)
{}


RmsdFilter::~RmsdFilter() {}

protocols::filters::FilterOP
RmsdFilter::clone() const {
	return new RmsdFilter( *this );
}

static basic::Tracer TR( "protocols.protein_interface_design.filters.RmsdFilter" );
core::Real
RmsdFilter::compute( core::pose::Pose const & pose ) const
{
	using namespace core;
	using namespace core::scoring;
	core::pose::Pose copy_pose = pose;
	core::pose::Pose native = *reference_pose_;
	core::Real rmsd( 0 );

	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native, basic::options::option[ basic::options::OptionKeys::in::file::native ] );
	}

	if ( !symmetry_ )
		runtime_assert( copy_pose.total_residue() == native.total_residue() );

	// generate temporary FArray
	FArray1D_bool selection_array( pose.total_residue(), false );
	if( selection_from_segment_cache_ ) core::pose::datacache::SpecialSegmentsObserver::set_farray_from_sso( selection_array, pose, true );
	else {
		for( std::list<core::Size>::const_iterator it = selection_.begin(); it!=selection_.end(); ++it ) {
			selection_array[*it-1] = true; // FArray1D is 0 indexed
		}
	}

	if ( reference_pose_->total_residue() == pose.total_residue() ) {
		if( superimpose_ ) {
			if ( symmetry_ )
				rmsd = core::scoring::sym_rmsd_with_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
			else
				rmsd = core::scoring::rmsd_with_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
		}
		else{
			rmsd = core::scoring::rmsd_no_super_subset( copy_pose, native, selection_array, core::scoring::is_protein_CA );
		}
	}
	return rmsd;
}

bool
RmsdFilter::apply( core::pose::Pose const & pose ) const {

	core::Real const rmsd( compute( pose ));
	TR << "RMSD over selected residues: " << rmsd ;
	if( rmsd <= threshold_ )
	{
		TR<<" passing."<<std::endl;
		return( true );
	}
	else TR<<" failing." << std::endl;
	return( false );
}

void
RmsdFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	out<<"RMSD: " << rmsd<<'\n';
}

core::Real
RmsdFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const rmsd( compute( pose ));
	return( (core::Real) rmsd );
}

void
RmsdFilter::parse_my_tag( utility::tag::TagPtr const tag, protocols::moves::DataMap & data_map, protocols::filters::Filters_map const &, protocols::moves::Movers_map const &, core::pose::Pose const & reference_pose )
{
	/// @details
	///if the save pose mover has been instantiated, this filter can calculate the rms
	///against the ref pose
	if( tag->hasOption("reference_name") ){
		reference_pose_ = protocols::rosetta_scripts::saved_reference_pose(tag,data_map );
	}
	else reference_pose_ = new core::pose::Pose( reference_pose );

	symmetry_ = tag->getOption<bool>( "symmetry", 0 );
	std::string chains = tag->getOption<std::string>( "chains", "" );
	if( chains != "" ) {
		core::Size chain_start( 0 );
		core::Size chain_end( 0 );

		utility::vector1<core::Size> chain_ends = reference_pose_->conformation().chain_endings();
		chain_ends.push_back( reference_pose_->total_residue() ); // chain_endings() doesn't count last residue as a chain ending (sigh)
		for( std::string::const_iterator chain_it = chains.begin(); chain_it != chains.end(); ++chain_it ) { // for each chain letter
			char const chain = *chain_it;
			TR.Debug << "Chain " << chain << " selected" << std::endl;
			for( utility::vector1<core::Size>::const_iterator it = chain_ends.begin(); it != chain_ends.end(); ++it ) {
				chain_end = *it;
				core::Size const chainid = reference_pose.residue( chain_end ).chain();
				if( reference_pose_->pdb_info()->chain( chain_end ) == chain ) { // if chain letter of chain_end == chain specified
					if( chain_end == reference_pose_->total_residue() ) { // all of this because the last residue doesn't count as a chain ending. why, god, why!?
						//core::Size const chainid = reference_pose_.residue( chain_end ).chain();
						for( core::Size i = 1; i <= reference_pose_->total_residue(); ++i ) {
							if( (core::Size)reference_pose_->residue(i).chain() == chainid ) { // first time we hit this, we're at the start of the chain in question
								chain_start = i;
								break;
							}
						}
					}
					else chain_start = reference_pose_->conformation().chain_begin( chainid );
					for( core::Size i = chain_start; i <= chain_end; ++i ) { // populate selection_ list
						selection_.push_back( i );
					}
				}
			}
		}
	}

	utility::vector0< utility::tag::TagPtr > const rmsd_tags( tag->getTags() );
	for( utility::vector0< utility::tag::TagPtr >::const_iterator it=rmsd_tags.begin(); it!=rmsd_tags.end(); ++it ) {
		utility::tag::TagPtr const rmsd_tag = *it;
		if( rmsd_tag->getName() == "residue" ) {
			core::Size const resnum( protocols::rosetta_scripts::get_resnum( rmsd_tag, *reference_pose_ ) );
			selection_.push_back( resnum );
		}
		if( rmsd_tag->getName() == "span" ) {
			core::Size const begin( protocols::rosetta_scripts::get_resnum( rmsd_tag, *reference_pose_, "begin_" ) );
			core::Size const end( protocols::rosetta_scripts::get_resnum( rmsd_tag, *reference_pose_, "end_" ) );
			runtime_assert( end > begin );
			runtime_assert( begin>=1);
			runtime_assert( end<=reference_pose_->total_residue() );
			for( core::Size i=begin; i<=end; ++i ) selection_.push_back( i );
		}
	}

	if( selection_.size() > 0 ) selection_.unique(); // make sure our selection list doesn't have redudancies
	TR.Debug << "RMSD selected residues: ";
	if( selection_.size() == 0 ) TR.Debug << "ALL" << std::endl;
	else {
		for( std::list<core::Size>::const_iterator it = selection_.begin(); it != selection_.end(); ++it ) {
			TR.Debug << *it << " ";
		}
	}
	TR.Debug << std::endl;

	superimpose_ = tag->getOption<bool>( "superimpose", 1 );
	threshold_ = tag->getOption<core::Size>( "threshold", 5 );
	if( tag->hasOption("rms_residues_from_pose_cache") ){
		selection_from_segment_cache_ = tag->getOption<bool>( "rms_residues_from_pose_cache", 1 );
		if( selection_.size() != 0 ) std::cerr << "Warning: in rmsd filter tag, both a span selection and the instruction to set the residues from the pose cache is given. Incompatible, defined span will be ignored." << std::endl;
	}

	TR<<"RMSD filter with superimpose=" << superimpose_ << " and threshold="<< threshold_ << " over residues ";
	if( selection_from_segment_cache_ ) TR << " that are in pose segment observer cache at apply time." << std::endl;
	else if( selection_.size() == 0 ) {
		TR << "ALL" << std::endl;
		for( core::Size i=1; i<=reference_pose.total_residue(); ++i ) selection_.push_back( i );
	}
	else {
		for( std::list<core::Size>::const_iterator it=selection_.begin(); it != selection_.end(); ++it ) TR << *it << " ";
		TR << std::endl;
	}
}

protocols::filters::FilterOP
RmsdFilterCreator::create_filter() const { return new RmsdFilter; }

std::string
RmsdFilterCreator::keyname() const { return "Rmsd"; }


} // filters
} // protein_interface_design
} // devel


