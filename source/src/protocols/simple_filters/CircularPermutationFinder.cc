// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/CircularPermutationFinderFilter.cc
/// @brief
/// @author Jonathan Weinstein & Sarel Fleishman


//Unit Headers
#include <protocols/simple_filters/CircularPermutationFinder.hh>
#include <protocols/simple_filters/CircularPermutationFinderCreator.hh>
#include <protocols/simple_filters/SSMotifFinderFilter.hh>
#include <utility/tag/Tag.hh>
//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <core/pose/selection.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <fstream>
#include <utility/io/izstream.hh>
#include <sstream>
#include <boost/foreach.hpp>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/PDBInfo.hh>
#include <protocols/toolbox/task_operations/ProteinInterfaceDesignOperation.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <protocols/toolbox/superimpose.hh>
#include <core/pack/task/TaskFactory.hh>

namespace protocols{
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.CircularPermutationFinder" );

protocols::filters::FilterOP
CircularPermutationFinderFilterCreator::create_filter() const { return new CircularPermutationFinder; }

std::string
CircularPermutationFinderFilterCreator::keyname() const { return "CircularPermutationFinder"; }

//default ctor
CircularPermutationFinder::CircularPermutationFinder() :
protocols::filters::Filter( "CircularPermutationFinder" ),
rmsd_( 0.0 ),
N_C_distance_( 0.0 ),
filename_( "" ),
align_only_interface_( true )
{
}

CircularPermutationFinder::~CircularPermutationFinder() {}

void
CircularPermutationFinder::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	rmsd( tag->getOption< core::Real >( "rmsd", 5.0 ) );
	filename( tag->getOption< std::string >( "filename" ) );
	N_C_distance( tag->getOption< core::Real >( "N_C_distance", 10.0 ) );
	align_only_interface( tag->getOption< bool >( "align_only_interface", true ) );

	TR<<"rmsd: "<<rmsd()<<" filename: "<<filename_<<" N_C_distance: "<<N_C_distance()<<" align_only_interfacce: "<<align_only_interface()<<std::endl;
}

void
write_to_file( std::string const filename, core::pose::Pose const & pose, core::Size const interface_begin, core::Size const cut, core::Size const interface_end, core::Real const rmsd ){
	core::pose::PDBInfoCOP pdb_info( pose.pdb_info() );

	char const chain( pdb_info->chain( interface_begin ) );
	int const resnum_begin( pdb_info->number( interface_begin )), resnum_cut( pdb_info->number( cut ) ), resnum_end( pdb_info->number( interface_end ) );

  std::ofstream file;
	file.open( filename.c_str(), std::ios::app );
	runtime_assert( file );
	file<<resnum_begin<<chain<<' '<<resnum_cut<<chain<<' '<<resnum_end<<chain<<' '<<rmsd<<'\n';
}

utility::vector1< numeric::xyzVector< core::Real > >
bb_coords( core::pose::Pose const & pose, utility::vector1< core::Size > const & positions ){
  utility::vector1< numeric::xyzVector< core::Real > > coords;

  coords.clear();
  BOOST_FOREACH( core::Size const pos, positions ){
    coords.push_back( pose.residue( pos ).xyz( "N" ) );
    coords.push_back( pose.residue( pos ).xyz( "CA" ) );
    coords.push_back( pose.residue( pos ).xyz( "C" ) );
    coords.push_back( pose.residue( pos ).xyz( "O" ) );
  }
  return coords;
}

core::Real
coord_rmsd( utility::vector1< numeric::xyzVector< core::Real > > const & c1, utility::vector1< numeric::xyzVector< core::Real > > const & c2 ){
	runtime_assert( c1.size() == c2.size() );

	core::Real sum_dist2( 0.0 );
	for( core::Size i = 1; i <= c1.size(); ++i )
		sum_dist2 += c1[ i ].distance_squared( c2[ i ] );

	sum_dist2 /= ( core::Real ) c1.size();
	sum_dist2 = std::sqrt( sum_dist2 );
	return( sum_dist2 );
}

core::Real
superimpose_parts( utility::vector1< numeric::xyzVector< core::Real > > & ref_coords, core::Size const cut ){
	core::Size const aa_atom_num( 4 );
	utility::vector1< numeric::xyzVector< core::Real > > rotated_coords( ref_coords );

	std::rotate( rotated_coords.begin(), rotated_coords.begin() + ( cut - 1 ) * aa_atom_num, rotated_coords.end());

  numeric::xyzMatrix< core::Real > rotation;
  numeric::xyzVector< core::Real > to_init_center, to_fit_center;

  protocols::toolbox::superposition_transform( rotated_coords, ref_coords, rotation, to_init_center, to_fit_center );

	return( coord_rmsd( rotated_coords, ref_coords ) );
 // apply_superposition_transform( pose, rotation, to_init_center, to_fit_center );
}

void
CircularPermutationFinder::circular_permutation( core::pose::Pose const & pose, core::Size const chainA, core::Size const chainB ) const{
	using namespace std;
	using core::Size;
	using core::Real;
	using utility::vector1;
	using namespace core::pose;
	using namespace protocols::toolbox::task_operations;
	using namespace core::pack::task;

	if( !pose.conformation().residue( pose.conformation().chain_begin( chainA ) ).is_protein() ||
			!pose.conformation().residue( pose.conformation().chain_begin( chainB ) ).is_protein() )
		return;
	TR<<"Testing chainA: "<<chainA<<" chainB: "<<chainB<<std::endl;
	vector1< PoseOP > pose_chains( pose.split_by_chain() );
	pose_chains[ chainA ]->append_pose_by_jump( *pose_chains[ chainB ], 1 );
	Pose const rewired_pose( *pose_chains[ chainA ] );

	Real const n_c_dist( atom_distance( rewired_pose.conformation().residue( 1 ), "CA",
																			rewired_pose.conformation().residue( rewired_pose.conformation().chain_end( 1 ) ), "CA" ) );
	TR<<"Distance = "<<n_c_dist<<std::endl;
	if( n_c_dist >= N_C_distance() ){
		TR<<" failed distance cutoff"<<std::endl;
		return;
	}

	ProteinInterfaceDesignOperationOP pido = new ProteinInterfaceDesignOperation;
	pido->repack_chain1( true );
	pido->design_chain1( true );
	pido->repack_chain2( false );
	pido->interface_distance_cutoff( 8.0 );
	TaskFactoryOP tf_interface( new TaskFactory);
	tf_interface->push_back(pido);
  vector1< Size > const chainA_interface( protocols::rosetta_scripts::residue_packer_states( rewired_pose, tf_interface, true, true));
	TR<<"Number of binding residues: "<<chainA_interface.size()<<std::endl;
	if( chainA_interface.size() <= 5 ){
		TR<<"Number too small; failing"<<std::endl;
		return;
	}

	map< Size/*residue number*/, Real/*rmsd*/ > rmsds;
	rmsds.clear();
	Size const min_res( chainA_interface[ 1 ] ), max_res( chainA_interface[ chainA_interface.size() ] );
  vector1< Size > aligned_positions;
	aligned_positions.clear();
	if( align_only_interface() )
		aligned_positions = chainA_interface;
	else{
		for( Size i = min_res; i <= max_res; ++i )
			aligned_positions.push_back( i );
	}
	vector1< numeric::xyzVector< core::Real > > aligned_coords( bb_coords( rewired_pose, aligned_positions ) );
	for( Size cut = 2; cut < aligned_positions.size(); ++cut ){
		Real const internal_rmsd = superimpose_parts( aligned_coords, cut );
		TR<<"At cut position: "<<aligned_positions[ cut ]<<" rmsd: "<<internal_rmsd<<std::endl;
		rmsds[ aligned_positions[ cut ] ] = internal_rmsd;
	}
	Real min_rmsd( 100000.0 );
	Size min_cut( 0 );
	for( map< Size, Real >::const_iterator rmsd_it = rmsds.begin(); rmsd_it != rmsds.end(); ++rmsd_it ){
		if( rmsd_it->second <= min_rmsd ){
			min_cut = rmsd_it->first;
			min_rmsd = rmsd_it->second;
		}
	}
	if( min_rmsd <= rmsd() )
		write_to_file( filename_, rewired_pose, chainA_interface[ 1 ], min_cut, chainA_interface[ chainA_interface.size() ], min_rmsd );
}


bool
CircularPermutationFinder::apply( core::pose::Pose const & pose ) const {
	using core::Size;
	using core::Real;

	Size const chain_num( pose.conformation().num_chains() );
	for( Size chaini = 1; chaini < chain_num; ++chaini ){
		for( Size chainj = chaini + 1; chainj <= chain_num; ++chainj ){
			circular_permutation( pose, chaini, chainj );
			circular_permutation( pose, chainj, chaini );
		}
	}

	return true;
}

void
CircularPermutationFinder::report( std::ostream &, core::pose::Pose const & ) const {
}

core::Real
CircularPermutationFinder::report_sm( core::pose::Pose const & ) const {
	return 0.0;
}

}
}
