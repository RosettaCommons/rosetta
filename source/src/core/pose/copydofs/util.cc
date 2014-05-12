// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/pose/copydofs/util.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <core/pose/copydofs/util.hh>
#include <core/pose/copydofs/CopyDofs.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/Pose.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/chemical/ResidueType.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "core.pose.copydofs.util" );

namespace core {
namespace pose {
namespace copydofs {

///////////////////////////////////////////////
///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief A very useful function that copies degrees of freedom from one pose to another.
/// @details res_map defines how to map residue numbers from the large pose to the smaller "scratch" pose.
/// @author rhiju, 2009.
/////////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					MiniPose const & scratch_pose,
					core::pose::ResMap const & res_map )
{

	std::map < id::AtomID , id::AtomID > atom_id_map;
	setup_atom_id_map( atom_id_map, res_map, pose ); // note that this is *not* careful about atom names, etc.
	copy_dofs( pose, scratch_pose, atom_id_map );

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs_match_atom_names(
					pose::Pose & pose,
					Pose const & scratch_pose )
{

	// Assumes the poses have the same number of residues
	runtime_assert( pose.total_residue() == scratch_pose.total_residue() );
	std::map< Size, Size > res_map;
	for ( Size n = 1; n <= pose.total_residue(); ++n ) res_map[n] = n;
	copy_dofs_match_atom_names( pose, scratch_pose, res_map );

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs_match_atom_names( //Parin Sripakdeevong Dec 27, 2011.
					pose::Pose & pose,
					MiniPose const & chunk_pose,
					core::pose::ResMap const & res_map )
{

	std::map < id::AtomID , id::AtomID > atom_id_map;
	setup_atom_id_map_match_atom_names( atom_id_map, res_map, pose, chunk_pose ); // note that this is CAREFUL about atom names, etc.
	copy_dofs( pose, chunk_pose, atom_id_map );

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					Pose const & scratch_pose )
{
	// Assumes the poses have the same number of residues
	runtime_assert( pose.total_residue() == scratch_pose.total_residue() );
	std::map< Size, Size > res_map;
	for ( Size n = 1; n <= pose.total_residue(); ++n ) res_map[n] = n;
	copy_dofs( pose, scratch_pose, res_map );
}
////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					Pose const & scratch_pose,
					core::pose::ResMap const & res_map )
{
	// Need to map atom numbers from big pose to scratch pose -- following assumes that
	// variant types are exactly the same.
	std::map < id::AtomID , id::AtomID > atom_id_map;
	setup_atom_id_map( atom_id_map, res_map, pose ); // note that this is *not* careful about atom names, etc.
	copy_dofs( pose, scratch_pose, atom_id_map );
}


////////////////////////////////////////////////////////////////////////////////////////////////
// slower!!!
void
copy_dofs_match_atom_names(
													 pose::Pose & pose,
													 Pose const & scratch_pose,
													 core::pose::ResMap const & res_map,
													 bool const backbone_only /* = false */,
													 bool const side_chain_only /* = false */,
													 bool const ignore_virtual /* = true */ )
{
	// Need to map atom numbers from big pose to scratch pose --
	std::map < id::AtomID , id::AtomID > atom_id_map;
	setup_atom_id_map_match_atom_names( atom_id_map, res_map, pose, scratch_pose, backbone_only, side_chain_only, ignore_virtual );
	copy_dofs( pose, scratch_pose, atom_id_map );
}


////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					Pose const & scratch_pose,
					std::map < id::AtomID , id::AtomID > const & atom_id_map )
{
	copy_dofs( pose, MiniPose( scratch_pose ), atom_id_map );
}

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					MiniPose const & scratch_pose,
					std::map < id::AtomID , id::AtomID > const & atom_id_map ){

	std::map< id::AtomID, Size > atom_id_domain_map = copydofs::blank_atom_id_domain_map( pose );
	copy_dofs( pose, scratch_pose, atom_id_map, atom_id_domain_map );

}

////////////////////////////////////////////////////////////////////////////////////////////////
void
copy_dofs(
					pose::Pose & pose,
					MiniPose const & scratch_pose,
					std::map < id::AtomID , id::AtomID > const & atom_id_map,
					std::map< id::AtomID, Size > const & atom_id_domain_map )
{

	copydofs::CopyDofs copy_dofs( scratch_pose, atom_id_map, atom_id_domain_map );
	copy_dofs.apply( pose );

}


///////////////////////////////////////////////////////////////////
void
setup_atom_id_map(
									std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
									ResMap const & res_map,
									core::pose::Pose const & pose )
{
	using namespace core::id;

	for ( ResMap::const_iterator
					it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

		Size const i = it->first; //Index in big pose.
		Size const i_scratch_pose = it->second; // Index in the little "chunk" or "scratch" pose

		//		std::cout << "setting up atom_id map " << i << " " << i_scratch_pose << std::endl;
		chemical::ResidueType const & rsd_type( pose.residue_type( i ) );
		Size count( 0 );

		/////////////////////////////////////////////////////////////////////////
		// Might be better to figure out correspondence based on atom name!
		//   Current code assumes numbering is similar -- not as robust but fast.
		/////////////////////////////////////////////////////////////////////////
		for ( Size j = 1; j <= rsd_type.natoms(); j++ ) {
			// HEY NEED TO FIX THIS LATER. DO WE NEED TO BE CAREFUL ABOUT VIRT?
			// MUCH BETTER TO MAKE VARIANTS MATCH *BEFORE* CALLING COPY_DOFS();
			//			if ( rsd_type.is_virtual( j ) ) continue;
			count++;
			atom_id_map[  AtomID( j, i ) ] = AtomID( count, i_scratch_pose );
		}

	}
}

///////////////////////////////////////////////////////////////////
void
setup_atom_id_map_match_atom_names(
									std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
									ResMap const & res_map,
									core::pose::Pose const & pose,
									core::pose::Pose const & reference_pose,
									bool const backbone_only /* = false */,
									bool const side_chain_only /* = false */,
									bool const ignore_virtual /* = true */ )
{
	using namespace core::id;

	for ( ResMap::const_iterator
					it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

		Size const i1 = it->first; //Index in big pose.
		Size const i2 = it->second; // Index in the little "chunk" or "scratch" pose

		chemical::ResidueType const & rsd_type1( pose.residue_type( i1 ) );
		chemical::ResidueType const & rsd_type2( reference_pose.residue_type( i2 ) );

		for ( Size j1 = 1; j1 <= rsd_type1.natoms(); j1++ ) {

			// Hey do we need this?
			if ( ignore_virtual && rsd_type1.is_virtual( j1 ) ) continue;

			std::string const & atom_name1 = rsd_type1.atom_name( j1 );

			if ( ! rsd_type2.has( atom_name1 ) ) continue;

			if ( backbone_only &&
					 !( j1 <= rsd_type1.last_backbone_atom() ) &&
					 !( j1 > rsd_type1.nheavyatoms()   && j1 < rsd_type1.first_sidechain_hydrogen() )	 ) continue;

			if ( side_chain_only &&
					 !( j1 > rsd_type1.last_backbone_atom() && j1 <= rsd_type1.nheavyatoms() ) &&
					 !( j1 >= rsd_type1.first_sidechain_hydrogen() )	 ) continue;

			Size const j2 = rsd_type2.atom_index( atom_name1 );

			// Hey do we need this?
			if ( ignore_virtual && rsd_type2.is_virtual( j2 ) ) continue;

			// this is new (Dec. 2010) -- be careful!
			//			if ( rsd_type1.is_virtual( j1 ) && ! rsd_type2.is_virtual( j2 ) ) continue;
			//			if ( ! rsd_type1.is_virtual( j1 ) && rsd_type2.is_virtual( j2 ) ) continue;

			atom_id_map[  AtomID( j1, i1 ) ] = AtomID( j2, i2 );
		}

	}
}

////////////////////////////////////////////////////////////////////////////////////////////////

void
setup_atom_id_map_match_atom_names( //June 16, 2011 Parin Sripakdeevong
									std::map < core::id::AtomID , core::id::AtomID > & atom_id_map,
									ResMap const & res_map,
									core::pose::Pose const & pose,
									MiniPose const & chunk_pose ){
	using namespace core::id;

	for ( ResMap::const_iterator it=res_map.begin(), it_end = res_map.end(); it != it_end; ++it ) {

		Size const full_seq_num = it->first; //Index in full pose.
		Size const chunk_seq_num = it->second; // Index in the little "chunk" or "scratch" pose

		chemical::ResidueType const & rsd_type1( pose.residue_type( full_seq_num ) );

		utility::vector1< utility::vector1< std::string > > const & chunk_atom_names_list = chunk_pose.atom_names_list();

		for(Size j1 = 1; j1 <= rsd_type1.natoms(); j1++ ) {
			for(Size j2=1; j2<=chunk_atom_names_list[chunk_seq_num].size(); j2++){

				std::string const & atom_name1 = rsd_type1.atom_name( j1 );
				std::string const & atom_name2 = chunk_atom_names_list[chunk_seq_num][j2];

				if(atom_name1==atom_name2){ //found matching atom_name!
					atom_id_map[  AtomID( j1, full_seq_num ) ] = AtomID( j2, chunk_seq_num );
					break;
				}
			}
		}

	}
}


/////////////////////////////////////////////////////////////////////
// specify dof_tolerance for speed -- changing dofs (even to the same
// value) triggers pose refold which can take some time.
void
apply_dofs( pose::Pose & pose, CopyDofsInfo const & copy_dofs_info,
						core::Real const dof_tolerance /* = 1.0e-5*/ ){
	for ( Size n = 1; n <= copy_dofs_info.size(); n++ ){

		if ( dof_tolerance > 0.0 ) {
			Real const dof_value_original = pose.dof( copy_dofs_info[n].first );
			Real const dof_value_new      = copy_dofs_info[n].second;
			if ( std::abs( dof_value_original - dof_value_new ) < dof_tolerance ) continue;
		}

		pose.set_dof( copy_dofs_info[n].first, copy_dofs_info[n].second );
	}
}

/////////////////////////////////////////////////////////////////////
std::map< id::AtomID, Size >
blank_atom_id_domain_map( pose::Pose const & pose ) {
	std::map< id::AtomID, Size > atom_id_domain_map;
	for ( Size i = 1; i <= pose.total_residue(); i++ ){
		for ( Size j = 1; j <= pose.residue_type( i ).natoms(); j++ ){
				atom_id_domain_map[ id::AtomID( j, i ) ] = 0;
		}
	}
	return atom_id_domain_map;
}


} //copydofs
} //pose
} //core
