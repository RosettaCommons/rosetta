// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file RBOutMover.cc
/// @brief

// Unit headers
#include <devel/splice/RBOutMover.hh>
#include <devel/splice/RBOutMoverCreator.hh>
#include <core/kinematics/FoldTree.hh>
#include <basic/Tracer.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
using basic::T;
using basic::Error;
using basic::Warning;
static thread_local basic::Tracer TR( "protocols.simple_moves.RBOutMover" );
#include <utility/tag/Tag.hh>

// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <utility/vector1.hh>
#include <boost/foreach.hpp>
#include <core/pose/Pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/VariantType.hh>
#include <core/scoring/constraints/SequenceProfileConstraint.hh>
#include <core/sequence/SequenceProfile.hh>
#include <core/conformation/Residue.hh>
#include <protocols/simple_moves/CutChainMover.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <fstream>
#include <algorithm>
#include <vector>


namespace devel {
namespace splice {

//
using namespace::protocols;

std::string
RBOutMoverCreator::keyname() const
{
	return RBOutMoverCreator::mover_name();
}

protocols::moves::MoverOP
RBOutMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new RBOutMover );
}

std::string
RBOutMoverCreator::mover_name()
{
	return "RBOut";
}

RBOutMover::RBOutMover(): moves::Mover("RBOut"),
    template_pdb_fname_(""),
    jump_dbase_fname_(""),
		jump_from_foldtree_( false )
{
}

RBOutMover::~RBOutMover()
{
}

bool compare(std::pair<core::Size,core::Size> p1, std::pair<core::Size, core::Size> p2) {return p1.first == p2.first;}

/// find disulfide pairs. Simply look for pairs of SG atoms on residues that are disulfide variants and within 3.0A
utility::vector1< std::pair< core::Size, core::Size > >
find_disulfs_in_range( core::pose::Pose const & pose, core::Size const start, core::Size const end ){
    utility::vector1< std::pair< core::Size, core::Size > > disulfs;
    disulfs.clear();

    runtime_assert( end <= pose.total_residue() && end >= 1 && start <= end && start >= 1 );
    for( core::Size resi = start; resi < end; ++resi ){
        for( core::Size resj = resi + 1 ; resj <= end; ++resj ){
            if( pose.residue( resi ).has_variant_type( core::chemical::DISULFIDE ) && pose.residue( resj ).has_variant_type( core::chemical::DISULFIDE ) && pose.residue( resi ).xyz( "SG" ).distance( pose.residue( resj ).xyz("SG")) <= 3.0 )
                    disulfs.push_back( std::pair< core::Size, core::Size >( resi, resj )) ;
        }
    }

    utility::vector1< std::pair<core::Size,core::Size> >::iterator it;
    it = std::unique(disulfs.begin(),disulfs.end(), compare);
    disulfs.resize(std::distance(disulfs.begin(),it) );

    return(disulfs);
}


core::kinematics::Jump
RBOutMover::get_disulf_jump( Pose & pose, core::pose::Pose const & template_pose )
{
    utility::vector1< std::pair< core::Size, core::Size > > const disulfs_in_template = find_disulfs_in_range( template_pose, 1, template_pose.total_residue() );
    runtime_assert( disulfs_in_template.size() == 2 );

    utility::vector1< std::pair< core::Size, core::Size > > const disulfs_in_pose = find_disulfs_in_range( pose, 1, pose.total_residue() );
    runtime_assert( disulfs_in_pose.size() >= 2 );
    utility::vector1< std::pair< core::Size, core::Size > > relevant_disulfs;
    relevant_disulfs.clear();
    core::Real const max_allowed_distance( 5.0 );

    for( utility::vector1< std::pair< core::Size, core::Size > >::const_iterator disulf_templ_it = disulfs_in_template.begin(); disulf_templ_it != disulfs_in_template.end(); ++disulf_templ_it ){
        for( utility::vector1< std::pair< core::Size, core::Size > >::const_iterator disulf_pose_it = disulfs_in_pose.begin(); disulf_pose_it != disulfs_in_pose.end(); ++disulf_pose_it ){
            core::Real const dist1( template_pose.residue( disulf_templ_it->first ).xyz( "CA" ).distance( pose.residue( disulf_pose_it->first ).xyz( "CA" ) ) );
            core::Real const dist2( template_pose.residue( disulf_templ_it->second ).xyz( "CA" ).distance( pose.residue( disulf_pose_it->second ).xyz( "CA" ) ) );
            TR<<"Template disulf: "<<disulf_templ_it->first<<"/"<<disulf_templ_it->second<<" pose disulf: "<<disulf_pose_it->first<<"/"<<disulf_pose_it->second<<" distances: "<<dist1<<"/"<<dist2<<std::endl;
            if( dist1 <= max_allowed_distance && dist2 <= max_allowed_distance ){
                TR<<"Found neighbour"<<std::endl;
                relevant_disulfs.push_back( std::pair< core::Size, core::Size >( disulf_pose_it->first, disulf_pose_it->second ) );
                break;
            }
            else
                TR<<"Not a neighbour"<<std::endl;
        }
    }
    runtime_assert( relevant_disulfs.size() == 2 );

    core::kinematics::FoldTree ft;
    ft.clear();

    /// We don't know what order the chains will come in. Is the first disulfide in the template also first on the pose? Not always the case...
    core::Size const first_disulf( std::min( relevant_disulfs[ 1 ].second, relevant_disulfs[ 2 ].second ));
    core::Size const second_disulf( std::max( relevant_disulfs[ 1 ].second, relevant_disulfs[ 2 ].second ));

    /// Following is an ugly foldtree just used to define the jump between the 2nd and 1st disulfides in the order they're observed in template pose
    using namespace core::kinematics;
    ft.add_edge( 1, first_disulf, Edge::PEPTIDE );
    ft.add_edge( second_disulf, pose.total_residue(), Edge::PEPTIDE );
    ft.add_edge( relevant_disulfs[ 1 ].second, relevant_disulfs[ 2 ].second, 1 );
    if( relevant_disulfs[ 1 ].second > relevant_disulfs[ 2 ].second )
        ft.add_edge( relevant_disulfs[ 1 ].second, relevant_disulfs[ 2 ].second - 1, Edge::PEPTIDE );
    else
        ft.add_edge( relevant_disulfs[ 2 ].second, relevant_disulfs[ 1 ].second + 1, Edge::PEPTIDE );
    ft.reorder(1);
    TR<<"New foldtree:" << ft <<std::endl;
    pose.fold_tree( ft );
    core::kinematics::Jump const pose_disulf_jump( pose.jump( 1 ) );

    TR << "The jump in pose is " << pose_disulf_jump << std::endl;
    return pose_disulf_jump;
}

void
RBOutMover::apply( Pose & pose )
{
	core::pose::Pose template_pose;
  core::import_pose::pose_from_pdb( template_pose, template_pdb_fname() );

  core::kinematics::Jump const pose_disulf_jump = ( jump_from_foldtree() ? pose.jump( 1 ) : get_disulf_jump( pose, template_pose ) );
//write to file
	std::ofstream data;
	data.open( jump_dbase_fname().c_str(), std::ios::out );
	if( !data.good() )
		utility_exit_with_message( "Unable to open jump dbase file for writing: " + jump_dbase_fname() + "\n" );
  data<<pose_disulf_jump<<std::endl;
  TR<<"Wrote jump info: "<<pose_disulf_jump<<std::endl;
}

std::string
RBOutMover::get_name() const {
	return RBOutMoverCreator::mover_name();
}

moves::MoverOP
RBOutMover::clone() const
{
	return moves::MoverOP( new RBOutMover( *this ) );
}

moves::MoverOP
RBOutMover::fresh_instance() const
{
	return moves::MoverOP( new RBOutMover );
}

void
RBOutMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap &,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const & )
{
    template_pdb_fname( tag->getOption< std::string >( "template_fname" ) );
    jump_dbase_fname( tag->getOption< std::string >( "jump_dbase_fname" ));
		jump_from_foldtree( tag->getOption< bool >( "jump_from_foldtree", false ) );
    TR<<"Template pdb fname: "<<template_pdb_fname()<<" jump_dbase_fname: "<<jump_dbase_fname()<<" jump_from_foldtree: "<<jump_from_foldtree()<<std::endl;
}

} // splice
} // devel

