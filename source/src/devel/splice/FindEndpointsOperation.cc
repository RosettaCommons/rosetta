// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/toolbox/task_operations/FindEndpointsOperation.cc
/// @brief
/// @author Sarelf Fleishman sarelf@uw.edu

// Unit Headers
#include <devel/splice/FindEndpointsOperation.hh>
#include <devel/splice/FindEndpointsOperationCreator.hh>
#include <core/conformation/Conformation.hh>
// Project Headers
// Utility Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pose/Pose.hh>
#include <protocols/toolbox/task_operations/DesignAroundOperation.hh>
#include <numeric/xyzVector.hh>

#include <boost/foreach.hpp>
#define foreach BOOST_FOREACH

// C++ Headers

using basic::Error;
using basic::Warning;
static basic::Tracer TR( "devel.splice.FindEndpointsOperation" );

namespace devel {
namespace splice {

using namespace core::pack::task::operation;
using namespace std;

FindEndpointsOperation::FindEndpointsOperation() :
    Cterm_offset_( 0 ),
    Nterm_offset_( 0 ),
    even_( true ),
		odd_( true ),
    distance_cutoff_( 18.0 ),
    neighbors_( 6 ),
    point_inside_( true )
{
}

FindEndpointsOperation::~FindEndpointsOperation() {}

core::pack::task::operation::TaskOperationOP
FindEndpointsOperationCreator::create_task_operation() const
{
	return new FindEndpointsOperation;
}

core::pack::task::operation::TaskOperationOP FindEndpointsOperation::clone() const
{
	return new FindEndpointsOperation( *this );
}

///@brief compute the number of residue neighbours target_res has within neighbours vector
core::Size
    neighbors_in_vector( core::pose::Pose const & pose, core::Size const target_res, utility::vector1< core::Size > const & neighbors, core::Real const dist_threshold, core::scoring::dssp::Dssp & dssp ){

    core::Size neighbor_count( 0 );
    for(utility::vector1<core::Size>::const_iterator res_jt = neighbors.begin(); res_jt != neighbors.end(); ++res_jt){
        if( target_res == *res_jt ) // don't count self as neighbour
            continue;
        if( std::abs( (int)target_res - (int)*res_jt ) <= 15 ) // make sure sequence separation is >=15
            continue;
        bool intervening_helix( false );
        for( core::Size i=std::min( target_res, *res_jt ); i<=std::max( target_res, *res_jt ); ++i ) /// make sure there is an intervening helix
            if( dssp.get_dssp_secstruct( i ) == 'H' )
                intervening_helix = true;
        if( !intervening_helix )
            continue;
        core::Real const distance(pose.residue( target_res ).xyz( "CA" ).distance( pose.residue( *res_jt ).xyz("CA") ) );

        if( distance<=dist_threshold )
            ++neighbor_count;
    }
    return( neighbor_count );
}

///@brief restricts to repacking all residues outside of design_shell_ around each residue
void
FindEndpointsOperation::apply( core::pose::Pose const & pose, core::pack::task::PackerTask & task ) const
{
	using namespace core::pack::task;
    using namespace protocols::toolbox::task_operations;

    core::scoring::dssp::Dssp dssp( pose );
    dssp.dssp_reduced(); // switch to simplified H E L notation
    utility::vector1< core::Size > strand_positions;
    strand_positions.clear();
    TR<<"DSSP: ";
    for( core::Size i = 1; i <= pose.total_residue() - 2; ++i ){
        if( dssp.get_dssp_secstruct( i ) == 'E' && dssp.get_dssp_secstruct( i + 1) == 'E' && dssp.get_dssp_secstruct( i + 2 ) == 'E' )
            strand_positions.push_back( i );
        TR<<dssp.get_dssp_secstruct(i);
    }
    TR<<std::endl;
    utility::vector1< core::Size > strand_ntermini;
    strand_ntermini.clear();
    core::Size prev_resnum( 99999 );
    foreach( core::Size const resi, strand_positions ){
        if( resi != prev_resnum + 1 )
            strand_ntermini.push_back( resi );
        prev_resnum = resi;
    }
    utility::vector1< core::Size > ntermini_w_neighbors;
    ntermini_w_neighbors.clear();
    for( utility::vector1< core::Size >::const_iterator res_it = strand_ntermini.begin(); res_it!=strand_ntermini.end();++res_it){
        core::Size const resi_neighbors( neighbors_in_vector( pose, *res_it, strand_ntermini, distance_cutoff(), dssp ));
        if( resi_neighbors >= neighbors() ){
            ntermini_w_neighbors.push_back( *res_it );
            TR<<*res_it<<" has "<<resi_neighbors<<" neighbors"<<std::endl;
        }//fi resi_neighbors
    }//for res_it

    utility::vector1< core::Size > ntermini_w_close_neighbors;
    ntermini_w_close_neighbors.clear();
    for( utility::vector1< core::Size >::const_iterator res_it = ntermini_w_neighbors.begin(); res_it!=ntermini_w_neighbors.end();++res_it){
        core::Size const close_neighbors( neighbors_in_vector( pose, *res_it, strand_ntermini, 10.0, dssp ));
        if( close_neighbors >= 2 )
            ntermini_w_close_neighbors.push_back( *res_it );
        else{
            core::Size const close_neighbors_one_down( neighbors_in_vector( pose, *res_it+1, strand_ntermini, 10.0, dssp ));
            if( close_neighbors_one_down >= 2 )
                ntermini_w_close_neighbors.push_back( *res_it + 1);
        }
        TR<<*res_it<<" has "<<close_neighbors<<" close neighbors"<<std::endl;
    }

    using namespace numeric;
    xyzVector< core::Real > radius(0.0,0.0,0.0);
    foreach( core::Size const resi, ntermini_w_close_neighbors ){
        radius+=pose.residue( resi ).xyz( "CA" );
    }
    radius/=ntermini_w_close_neighbors.size();
 //   TR<<"Center of circle at: "<<radius<<std::endl;
    utility::vector1< core::Size > residues_pointing_in;
    residues_pointing_in.clear();
    foreach( core::Size const res, ntermini_w_close_neighbors ){
        core::Size resi = res;
        if( pose.residue( res ).name3() == "GLY" ){
            resi = res+1;
        }
        core::Real const Ca_dist( pose.residue( resi ).xyz("CA").distance(radius));
        core::Real const Cb_dist( pose.residue( resi ).xyz("CB").distance(radius));
        TR<<"Ca_dist/Cb_dist: "<<Ca_dist<<'/'<<Cb_dist<<std::endl;
        if( Ca_dist + 0.5 <= Cb_dist ){
            TR<<"Residue "<<resi<<" pointing out from barrel. Choosing instead ";
            if( pose.residue( res ).name3() == "GLY" ){
                TR<<resi-1<<std::endl;
                residues_pointing_in.push_back( res );
            }
            else{
                TR<<resi<<std::endl;
                residues_pointing_in.push_back( resi+1 );
            }
        }
        else
            residues_pointing_in.push_back( resi );
    }

    protocols::toolbox::task_operations::DesignAroundOperation dao;
    dao.repack_shell( 0.01 );
    dao.design_shell( 0.01 );
    TR<<"found "<<residues_pointing_in.size()<<" nterminal residues with enough neighbors"<<std::endl;
		bool curr_odd( true );
    foreach( core::Size const res, residues_pointing_in ){
			if( ( curr_odd && odd() ) || ( !curr_odd && even() ) )
        dao.include_residue( res );
			curr_odd = !curr_odd;
    }

    dao.apply( pose, task);
}

void
FindEndpointsOperation::parse_tag( TagCOP tag , DataMap & )
{
    Cterm_offset( tag->getOption< core::Size >( "Cterm_offset", 0 ) );
    Nterm_offset( tag->getOption< core::Size >( "Nterm_offset", 0 ));
    even( tag->getOption< bool >( "even", true ));
    odd( tag->getOption< bool >( "odd", true ));
    neighbors( tag->getOption< core::Size >( "neighbors", 6 ));
    distance_cutoff( tag->getOption< core::Real >( "distance_cutoff", 18.0 ));
    point_inside( tag->getOption< bool >( "point_inside", true ));
    TR<<"Cterm_offset: "<<Cterm_offset()<<" Nterm_offset: "<<Nterm_offset()<<" even: "<<even()<<" odd: "<<odd()<<" neighbors: "<<neighbors()<<" distance_cutoff: "<<distance_cutoff()<<" point_inside: "<<point_inside()<<std::endl;
}

} //namespace splice
} //namespace devel
