// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/enzdes/RemoveLigandMover.cc
/// @brief
/// @author Javier Castellanos (javiercv@uw.edu)

// unit headers
#include <protocols/enzdes/RemoveLigandFilter.hh>
#include <protocols/enzdes/RemoveLigandFilterCreator.hh>

// package headers

// Project Headers
#include <core/pose/Pose.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/types.hh>

#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/protein_interface_design/filters/RmsdFilter.hh>
#include <protocols/toolbox/pose_manipulation/pose_manipulation.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
//utility headers
#include <utility/tag/Tag.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

using namespace core;
using namespace core::scoring;

static thread_local basic::Tracer TR( "protocols.enzdes.RemoveLigandFilter" );
namespace protocols {
namespace enzdes {

using namespace protocols::filters;
using namespace protocols::moves;
using namespace utility::tag;
using protocols::simple_moves::MinMover;
using protocols::protein_interface_design::filters::RmsdFilter;
using core::kinematics::MoveMapOP;
using core::pose::PoseOP;
using core::pose::Pose;
using core::Real;

std::string
RemoveLigandFilterCreator::keyname() const
{
    return "RemoveLigandFilter";
}

FilterOP
RemoveLigandFilterCreator::create_filter() const
{
    return new RemoveLigandFilter;
}
   


RemoveLigandFilter::RemoveLigandFilter( ):
		Filter( "RemoveLigandFilter" ),
    threshold_( 99.99 ),
    mover_( new MinMover ),
    filter_( new RmsdFilter )
{
    MoveMapOP movemap = new core::kinematics::MoveMap;
    movemap->set_bb(1);
    movemap->set_chi(1);
    
    MinMover* min_mover = dynamic_cast< MinMover* >( mover_.get() );
    min_mover->movemap( movemap );
    core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
    min_mover->score_function( scorefxn );
}
 
    
RemoveLigandFilter::RemoveLigandFilter( Real threshold ):
		Filter( "RemoveLigandFilter" ),
    threshold_( threshold ),
    mover_( new MinMover ),
    filter_( new RmsdFilter )
{ }
 
RemoveLigandFilter::RemoveLigandFilter( RemoveLigandFilter const & rval ):
		Filter( "RemoveLigandFilter" ),
    threshold_( rval.threshold_ ),
    mover_( rval.mover_ ),
    filter_( rval.filter_ )
{
    
}
    
Real
RemoveLigandFilter::report_sm( Pose const & pose ) const
{
    
    Pose no_lig_pose( pose );
    protocols::toolbox::pose_manipulation::remove_non_protein_residues( no_lig_pose );
    if( filter_->get_type() == "Rmsd" )
    {
        PoseOP init_pose = new Pose( no_lig_pose );
        RmsdFilter* rmsd_filter = dynamic_cast< RmsdFilter* >( filter_.get() );
        rmsd_filter->reference_pose( init_pose );
        rmsd_filter->superimpose( true );
        std::list< core::Size > selection;
        for(Size i = 1; i <= init_pose->n_residue(); i++ )
            selection.push_back( i );
        rmsd_filter->selection( selection );
        
        mover_->apply( no_lig_pose );
        return rmsd_filter->report_sm( no_lig_pose );
        
    } else {
        Real start_score = filter_->report_sm( no_lig_pose );
        mover_->apply( no_lig_pose );
        Real end_score = filter_->report_sm( no_lig_pose );
        return end_score - start_score;
    }
}
    
bool
RemoveLigandFilter::apply(Pose const & pose) const
{
    return report_sm( pose ) < threshold_;
}

void
RemoveLigandFilter::parse_my_tag( TagCOP const tag, basic::datacache::DataMap &, Filters_map const & filters, Movers_map const & movers, Pose const &)
{
    threshold_ = tag->getOption< Real >( "threshold", 3.0 );
    const std::string mover_name = tag->getOption< std::string >( "mover", "");
    const std::string filter_name = tag->getOption< std::string >( "filter", "");
    
    Movers_map::const_iterator  find_mover ( movers.find( mover_name ));
	Filters_map::const_iterator find_filter( filters.find( filter_name ));
	if( find_mover == movers.end() && mover_name != "" ) {
		TR.Error << "ERROR !! mover not found in map: \n" << tag << std::endl;
		runtime_assert( find_mover != movers.end() );
	}
	if( find_filter == filters.end() && filter_name != "" ) {
		TR.Error << "ERROR !! filter not found in map: \n" << tag << std::endl;
		runtime_assert( find_filter != filters.end() );
	}

    if( mover_name != "" )
		mover_ = find_mover->second;
    if( filter_name != "" )
        filter_ = find_filter->second;
}

} // enzdes
} // protocols
