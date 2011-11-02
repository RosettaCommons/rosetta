// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file  core/conformation/symmetry/util.hh
/// @brief utility functions for handling of symmetric conformations
/// @author Ingemar Andre

// Unit headers
#include <core/conformation/symmetry/util.hh>

// Package headers
#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
//#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

//#include <core/pose/Pose.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Edge.hh>
#include <core/id/AtomID.hh>

// Basic Headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/symmetry.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>

// Numeric headers
#include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/xyzTriple.hh>


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>

#include <utility/vector1.hh>



namespace core {
namespace conformation {
namespace symmetry {

static basic::Tracer TR("core.conformation.symmetry.util");
static numeric::random::RandomGenerator RG(408529); // <- Magic number, do not change it!!!

/// @details  Attempt to detect whether a conformation is symmetric
bool
is_symmetric( conformation::Conformation const & conf )
{
	return ( dynamic_cast< conformation::symmetry::SymmetricConformation const * >( &conf ) );
}

bool
is_symmetric( conformation::symmetry::SymmetryInfo const & symminfo )
{
	return symminfo.get_use_symmetry();
}

/// @details Generate a symmetric conformation from a monomeric pose and symmetry information
/// stored in the SymmData object
conformation::symmetry::SymmetricConformationOP
setup_symmetric_conformation(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata
)
{
	Size njump_orig = src_conformation.fold_tree().num_jump();

	// maybe a little inefficient: first build a standard conformation, then construct symmetric conformation at the end
  conformation::Conformation conf( src_conformation );

	// if anchor is 0 then use com
	if ( symmdata.get_anchor_residue() == 0 ) {
		core::Size new_anchor ( core::conformation::symmetry::residue_center_of_mass( conf, 1, conf.size() ) );
		symmdata.set_anchor_residue( new_anchor );
	}

// recenter the pose to the origin?
	if ( symmdata.get_recenter() ) {
		core::conformation::symmetry::recenter( conf, symmdata );
	}

	// Setup temporary variables
  Size const nres_monomer( conf.size() );
  Size const njump_monomer( conf.fold_tree().num_jump() );

	//runtime_assert( nres_monomer == symmdata.get_nres_subunit() );
	//runtime_assert( njump_monomer == symmdata.get_njump_subunit() );

  // setup the pseudo residues
	std::map< std::string, conformation::symmetry::VirtualCoordinate > coords ( symmdata.get_virtual_coordinates() );
	std::map< std::string, std::pair< std::string, std::string > > virtual_connects( symmdata.get_virtual_connects() );
	std::map< Size, std::string > virtual_num_to_id (symmdata.get_virtual_num_to_id() );

	// Setup virtual residues
	{
    // create the new residue
    chemical::ResidueTypeSet const & rsd_set( src_conformation.residue(1).residue_type_set() );
		conformation::ResidueOP rsd( conformation::ResidueFactory::create_residue( rsd_set.name_map( "VRT" ) ) );

    // root the fold_tree at this pseudo residue
//    {
//      kinematics::FoldTree f( conf.fold_tree() );
//      f.reorder( conf.size() );
//      conf.fold_tree( f );
//    }

		for ( Size i =1; i <= virtual_num_to_id.size(); ++i ) {
			if ( virtual_num_to_id.find( i ) == virtual_num_to_id.end() ) {
				utility_exit_with_message( "[ERROR] Cannot find jump number..." );
			}
			std::string tag ( virtual_num_to_id.find( i )->second );
			if ( coords.find( tag ) == coords.end() ) {
        utility_exit_with_message( "[ERROR] Cannot find tag " + tag );
      }
			conformation::symmetry::VirtualCoordinate virt_coord( coords.find( tag )->second );
      rsd->set_xyz( "ORIG", virt_coord.get_origin() );
      rsd->set_xyz( "X", virt_coord.get_x().normalized() + virt_coord.get_origin() );
      rsd->set_xyz( "Y", virt_coord.get_y().normalized() + virt_coord.get_origin() );
      conf.append_residue_by_jump( *rsd, conf.size() ); //append it to the end of the monomeri
		}
	}

	// now insert the other conformations. This step will create a fold_tree that will be discarded at a later stage.
	// Optimally we want to set the fold tree directly but that turns out to be difficult. TODO change
  kinematics::Jump const monomer_jump( conf.jump( njump_monomer+1 ) );
  for ( Size i=1; i< symmdata.get_subunits(); ++i ) {
		 Size const insert_seqpos( nres_monomer * i + 1 ); // desired sequence position of residue 1 in inserted conformation
    Size const insert_jumppos( njump_monomer * i + 1 ); // desired jump number of 1st jump in inserted conformation
    Size const anchor_pseudo( nres_monomer * i + i + 1 ); // in current numbering
    Size const new_jump_number( njump_monomer*( i + 1 ) + i + 1 );
    conf.insert_conformation_by_jump( src_conformation, insert_seqpos, insert_jumppos, anchor_pseudo, new_jump_number);
//		conf.set_jump( new_jump_number, monomer_jump );
  }
	// Now generate the fold_tree from the symmdata and set it into the conformation. In the future we want to be able to use a
	// atom tree for this.
	kinematics::FoldTree f ( core::conformation::symmetry::set_fold_tree_from_symm_data( src_conformation, symmdata ) );
	// set the root of the fold tree
	Size new_root ( conf.size() - coords.size() + symmdata.get_root() );
	f.reorder( new_root );
	TR.Debug << f << std::endl;
	conf.fold_tree( f );

  // now build the symmetry info object
	conformation::symmetry::SymmetryInfo const symm_info( symmdata, nres_monomer, njump_monomer );

  // now create the symmetric conformation
	conformation::symmetry::SymmetricConformationOP symm_conf = new conformation::symmetry::SymmetricConformation( conf, symm_info );

	// apply independent jumps so the structure is created symmetrically
	for ( Size i=1; i<= f.num_jump(); ++i ) {
		if ( symm_info.jump_is_independent( i ) )
			symm_conf->set_jump( i, conf.jump( i ) );
	}

	// renumber the dof information to take the internal jumps into consideration
	core::conformation::symmetry::shift_jump_numbers_in_dofs( *symm_conf, njump_orig * symmdata.get_subunits() );

	return symm_conf;

}

// @details Make a fold tree from data stored in the SymmData object. This would have to make sense.
// this function would probably enjoy more checks...
kinematics::FoldTree
set_fold_tree_from_symm_data(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata
)
{
	using namespace kinematics;
	FoldTree f, f_orig = src_conformation.fold_tree();
	f.clear();

	// Save number of monomer residues, subunits and virtuals
	Size num_res_subunit( src_conformation.size() );
	Size subunits( symmdata.get_subunits() );
	Size num_cuts_subunit( f_orig.num_cutpoint() );
	Size num_jumps_subunit( f_orig.num_jump() );
	Size num_virtuals( symmdata.get_virtual_coordinates().size() );
	Size num_virtual_connects ( symmdata.get_virtual_connects().size() );
	Size num_res_real( num_res_subunit*subunits );
	Size anchor ( symmdata.get_anchor_residue() );

	// Check that input data makes sense
	if ( anchor > num_res_subunit ) {
			utility_exit_with_message( "Anchor outside the subunit..." );
	}

	// Store number of jumps and cuts
	Size njumps( 0 );
	Size total_num_cuts( num_cuts_subunit*subunits + subunits + num_virtuals - 1 );
	Size total_num_jumps( num_jumps_subunit*subunits + num_virtual_connects );

	// store information from subunit foldtree
	utility::vector1< int > cuts_subunit( f_orig.cutpoints() );
	ObjexxFCL::FArray1D_int cuts( total_num_cuts );
	ObjexxFCL::FArray2D_int jump_points_subunit( 2, num_jumps_subunit ),
													jumps( 2, total_num_jumps );

	for ( Size i = 1; i<= f_orig.num_jump(); ++i ) {
		jump_points_subunit(1,i) = f_orig.upstream_jump_residue(i);
		jump_points_subunit(2,i) = f_orig.downstream_jump_residue(i);
	}

	// Now add jumps and cuts for the other subunits
	for ( Size i = 0; i < subunits; ++i ) {
    for ( Size j = 1; j <= num_jumps_subunit; ++j ) {
      ++njumps;
			int new_cut_pos( i*num_res_subunit + cuts_subunit.at( j ) );
			cuts( njumps ) = new_cut_pos;
      	int new_jump_pos1( i*num_res_subunit + jump_points_subunit(1,j) );
			int new_jump_pos2( i*num_res_subunit + jump_points_subunit(2,j) );
      	jumps(1, njumps ) = std::min(new_jump_pos1,new_jump_pos2);
			jumps(2, njumps ) = std::max(new_jump_pos1,new_jump_pos2);
    }
  }

	// why copy the maps?
	std::map< std::string, std::pair< std::string, std::string > > const virtual_connect( symmdata.get_virtual_connects() );
	std::map< std::string, Size > virtual_id_to_num (symmdata.get_virtual_id_to_num() );
	std::map< Size, std::string > virtual_num_to_id (symmdata.get_virtual_num_to_id() );
	std::map< std::string, Size > virtual_id_to_subunit (symmdata.get_virt_id_to_subunit_num() );

	std::map< std::string, std::pair< std::string, std::string > >::const_iterator it_start = virtual_connect.begin();
	std::map< std::string, std::pair< std::string, std::string > >::const_iterator it_end = virtual_connect.end();
	for ( std::map< std::string, std::pair< std::string, std::string > >::const_iterator it = it_start; it != it_end; ++it ) {
		std::pair< std::string, std::string > connect( it->second );
		std::string pos_id1( connect.first );
		std::string pos_id2( connect.second );
		if ( pos_id2 == "SUBUNIT" ) {
			int subunit = virtual_id_to_subunit.find( pos_id1 )->second;
			++njumps;
			cuts( njumps ) = num_res_subunit*subunit;
			int jump_pos1( ( subunit-1 )*num_res_subunit + anchor );
			int jump_pos2( num_res_real + virtual_id_to_num.find( pos_id1 )->second );
			jumps(1, njumps ) = std::min(jump_pos1,jump_pos2);
			jumps(2, njumps ) = std::max(jump_pos1,jump_pos2);
		}
	}

	// Read jump structure from symmdata

	// Read the virtual connect data and add the new connections. They are stored in virtual_connects. We also need to
	// know the mapping between mapping between name of virtual residues and the jump number
	Size pos( 0 );
	for ( std::map< std::string, std::pair< std::string, std::string > >::const_iterator it = it_start; it != it_end; ++it ) {
				std::pair< std::string, std::string > connect( it->second );
		std::string pos_id1( connect.first );
		std::string pos_id2( connect.second );
		// We have already added the jumps from virtual residues to their corresponding subunits
		if ( pos_id2 == "SUBUNIT" ) continue;
		++njumps;
		++pos;
		assert( virtual_id_to_num.find( pos_id1 ) != virtual_id_to_num.end() && virtual_id_to_num.find( pos_id2 ) != virtual_id_to_num.end() );
		Size pos1 ( virtual_id_to_num.find( pos_id1 )->second );
		Size pos2 ( virtual_id_to_num.find( pos_id2 )->second );
     	cuts( njumps ) = num_res_real + pos;
		jumps(1, njumps ) = num_res_real + std::min(pos1,pos2);
		jumps(2, njumps ) = num_res_real + std::max(pos1,pos2);
	}

	// Now create foldtree
	f.tree_from_jumps_and_cuts( num_res_real + num_virtuals, njumps, jumps, cuts );
	f.reorder( num_res_real + 1 );
	return f;
}

// this function is STILL under construction...
kinematics::FoldTree
replaced_symmetric_foldtree_with_new_monomer(
	kinematics::FoldTree symm_f,
	conformation::symmetry::SymmetryInfo symmetry_info,
	kinematics::FoldTree monomer_f
)
{
	using namespace kinematics;
	FoldTree f=symm_f, f_new=monomer_f;

	// Save number of monomer residues, subunits and virtuals
	Size num_res_subunit( symmetry_info.num_independent_residues() );
	Size subunits( symmetry_info.subunits() );
	//Size anchor ( symmdata.get_anchor_residue() );


	// Now add jumps and cuts for the other subunits
	for ( Size i = 0; i < subunits; ++i ) {
		int begin ( i*num_res_subunit+1 );
		int end ( (i+1)*num_res_subunit );
		f.delete_segment( begin, end );
		//f.insert_fold_tree_by_jump(f, insert_seqpos, insert_jumppos, anchor_pos, anchor_jump_number, anchor_atom, root_atom);

  }

	return f;
}


// @details center the pose at the origin. Use the CA of the anchor residue
// as the anchor point
void
recenter(
	conformation::Conformation & src_conformation,
	conformation::symmetry::SymmData & symmdata
)
{
	conformation::Residue anchor ( src_conformation.residue( symmdata.get_anchor_residue() ) );
	Vector trans( anchor.xyz( "CA" ) );
	for ( Size i = 1; i <= src_conformation.size(); ++i ) {
    for ( Size j = 1; j <= src_conformation.residue_type(i).natoms(); ++j ) {
      id::AtomID id( j, i );
      src_conformation.set_xyz( id, src_conformation.xyz(id) - trans );
    }
  }
}

// @details shift jump numbers in dof
void
shift_jump_numbers_in_dofs(
	conformation::Conformation & conformation,
	Size shift
)
{
	using namespace core::conformation::symmetry;

	runtime_assert( is_symmetric( conformation ) );
	SymmetricConformation & symm_conf (
	      dynamic_cast<SymmetricConformation & > ( conformation ) );
	SymmetryInfoOP symm_info( symm_conf.Symmetry_Info() );

	std::map< Size, SymDof > dofs ( symm_conf.Symmetry_Info()->get_dofs() );
	std::map< Size, SymDof > dofs_shifted;
	std::map< Size, SymDof >::iterator it;
	std::map< Size, SymDof >::iterator it_begin = dofs.begin();
	std::map< Size, SymDof >::iterator it_end = dofs.end();

	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first + shift );
		dofs_shifted.insert( std::make_pair( jump_nbr, (*it).second ) );
	}
	symm_info->set_dofs( dofs_shifted );
}


kinematics::FoldTree
get_asymm_unit_fold_tree( core::conformation::Conformation const &conf ) {
  if( !is_symmetric( conf ) ) {
		return conf.fold_tree();
	}

	// basic idea: only take edges with at least one end in the 1st subunit
	kinematics::FoldTree const &f = conf.fold_tree();
	kinematics::FoldTree f_new;

	conformation::symmetry::SymmetricConformation const & symm_conf ( dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( conf ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	Size nres_subunit ( symm_info->num_independent_residues() );

	for ( core::kinematics::FoldTree::const_iterator it=f.begin(), eit=f.end(); it != eit; ++it) {
		if ( it->start() <=	(int)nres_subunit && it->stop() <=(int)nres_subunit ) {
			f_new.add_edge( *it );
		} else if ( it->stop() <= (int)nres_subunit ) { // this is the jump to the subunit
			core::kinematics::Edge e_new = *it;
			e_new.start() = nres_subunit + 1;
			f_new.add_edge( e_new );
		}
	}
	f_new.renumber_jumps();

	//TR.Error << "get_asymm_unit_fold_tree() called with " << f << std::endl;
	//TR.Error << "get_asymm_unit_fold_tree() returning " << f_new << std::endl;
	return f_new;
}

// this function is directly stolen from the docking code in protocols. The code duplication is introduced
// to avoid dependencies of core functionality from protocols
int
residue_center_of_mass(
	conformation::Conformation const & conformation,
	int const start,
	int const stop
)
{
	Vector center( 0.0 );
	for ( int i=start; i<=stop; ++i ) {
		//if( !pose.residue( i ).is_protein() ) continue;
		if( !conformation.residue( i ).is_protein()) {
			Vector ca_pos( conformation.residue( i ).nbr_atom_xyz() );
			center += ca_pos;
		} else {
			Vector ca_pos( conformation.residue( i ).atom( "CA" ).xyz() );
			center += ca_pos;
			}
	}
	center /= (stop-start+1);

	return core::conformation::symmetry::return_nearest_residue( conformation, start, stop, center );
}

// this function is directly stolen from the docking code in protocols. The code duplication is introduced
// to avoid dependencies of protocols on core functionality
int
return_nearest_residue(
	conformation::Conformation const & conformation,
	int const begin,
	int const end,
	Vector center
)
{
	Real min_dist = 9999.9;
	int res = 0;
	for ( int i=begin; i<=end; ++i )
	{
		Vector ca_pos;
		if( !conformation.residue( i ).is_protein() ){
			ca_pos = conformation.residue( i ).nbr_atom_xyz();
			} else {
			ca_pos = conformation.residue( i ).atom( "CA" ).xyz() ;
			}

		ca_pos -= center;
		Real tmp_dist( ca_pos.length_squared() );
		if ( tmp_dist < min_dist ) {
			res = i;
			min_dist = tmp_dist;
		}
	}
	return res;
}


} // symmetry
} // pose
} // core
