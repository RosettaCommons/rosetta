// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file  core/pose/symmetry/util.hh
/// @brief utility functions for handling of symmetric conformations
/// @author Ingemar Andre

// Unit headers
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/util.hh>

#include <core/conformation/symmetry/SymmData.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/MirrorSymmetricConformation.hh>
#include <core/conformation/carbohydrates/GlycanTreeSet.hh>
#include <core/conformation/util.hh>
#include <core/scoring/symmetry/SymmetricEnergies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/chains_util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/chemical/AA.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/conformation/symmetry/SymmetryInfo.hh>

// Utility functions
#include <basic/options/option.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/fold_and_dock.OptionKeys.gen.hh>
#include <core/id/AtomID.hh>
#include <numeric/random/random.hh>

// Package Headers
#include <core/kinematics/Edge.hh>
#include <core/io/HeaderInformation.hh>
#include <core/pose/PDBInfo.hh>

#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>

#include <utility/vector1.hh>

#include <core/conformation/symmetry/SymmetryTransform.hh> // IWYU MANUAL
#include <core/kinematics/Jump.hh> // AUTO IWYU For Jump

namespace core {
namespace pose {
namespace symmetry {

static basic::Tracer TR( "core.pose.symmetry.util" );

bool
is_symmetric( scoring::Energies const & energies )
{
	return ( dynamic_cast< scoring::symmetry::SymmetricEnergies const * >( &energies ) );
}

/// @details  Attempt to detect whether a pose is symmetric
bool
is_symmetric( pose::Pose const & pose )
{
	return conformation::symmetry::is_symmetric( pose.conformation() );
}

/// @details  Attempt to detect whether a pose is symmetric
bool
is_mirror_symmetric( pose::Pose const & pose )
{
	return conformation::symmetry::is_mirror_symmetric( pose.conformation() );
}

/// @brief Convenience function for the number of residues in a subunit.
/// Will return the total size for a asymmetric pose
core::Size
get_nres_asymmetric_unit( pose::Pose const & pose ) {
	if ( is_symmetric( pose ) ) {
		core::conformation::symmetry::SymmetricConformation conf ( dynamic_cast<core::conformation::symmetry::SymmetricConformation const &>( pose.conformation()) );
		return conf.Symmetry_Info()->get_nres_subunit();
	} else {
		return pose.size();
	}
}

/// @details Attempts to grab the symmetry info if the pose is symmetric
conformation::symmetry::SymmetryInfoCOP symmetry_info( pose::Pose const & pose )
{
	runtime_assert( is_symmetric( pose ) );
	auto const & SymmConf (
		dynamic_cast<conformation::symmetry::SymmetricConformation const &> ( pose.conformation()) );
	return SymmConf.Symmetry_Info();
}

/// @details  Attempt to detect whether a conformation is symmetric
bool
scorefxn_is_symmetric( conformation::Conformation const & conf )
{
	return ( dynamic_cast< conformation::symmetry::SymmetricConformation const * >( &conf ) );
}

/// @details  Attempt to detect whether a pose is symmetric
bool
scorefxn_is_symmetric( pose::Pose const & pose )
{
	return conformation::symmetry::is_symmetric( pose.conformation() );
}

/// @details constructs a symmetric pose with a symmetric conformation and energies object
///    from a monomeric pose and symmetryinfo object.
/// Unlike the version of make_symmetric_pose from symmdata, this does not expand the
///    pose; it assumes the symmetric fold tree and residues are already present
/// For example, this is used to reconstruct a symm pose from a silent file
/// All data in the pose datacache will be cleared
void
make_symmetric_pose(
	pose::Pose & pose,
	conformation::symmetry::SymmetryInfo symmetry_info
)
{
	using namespace basic::options;

	// if the pose contains mirrored virtuals, we need to create a mirrored symmetric conformation
	bool has_mirror=false;
	for ( core::Size i=1; i<=pose.size() && !has_mirror; ++i ) {
		if ( pose.residue_type(i).is_inverted_virtual_residue() ) {
			has_mirror=true;
		}
	}

	conformation::symmetry::SymmetricConformationOP symm_conf;
	if ( has_mirror ) {
		symm_conf = conformation::symmetry::SymmetricConformationOP (
			new core::conformation::symmetry::MirrorSymmetricConformation( pose.conformation(), symmetry_info ) );
	} else {
		symm_conf = conformation::symmetry::SymmetricConformationOP (
			new core::conformation::symmetry::SymmetricConformation( pose.conformation(), symmetry_info ) );
	}

	scoring::symmetry::SymmetricEnergiesOP symm_energy( new scoring::symmetry::SymmetricEnergies( pose.energies()) );

	pose.set_new_conformation( symm_conf );
	pose.set_new_energies_object( symm_energy );

	pose::PDBInfoOP pdb_info( new pose::PDBInfo( pose, true ) );

	//fpd if the input pdb info is valid copy it
	if ( pose.pdb_info() && pose.pdb_info()->nres() == pose.size() ) {
		pdb_info = utility::pointer::make_shared< pose::PDBInfo >( *(pose.pdb_info()) );
	}
	pose.pdb_info( pdb_info );

	debug_assert( is_symmetric( pose ) );

	if ( basic::options::option[ basic::options::OptionKeys::symmetry::detect_bonds ] ) {
		pose.conformation().detect_bonds();
	}

}

/// @details constructs a symmetric pose with a symmetric conformation and energies object from a monomeric pose
/// and symmData object
void
make_symmetric_pose(
	pose::Pose & pose,
	conformation::symmetry::SymmData & symmdata,
	bool keep_pdb_info_labels /* false */
)
{

	using namespace basic::options;

	pose::PDBInfoOP pdb_info_src( pose.pdb_info() );
	if ( !pose.pdb_info() ) {
		pdb_info_src = utility::pointer::make_shared< pose::PDBInfo >( pose, true );
	}

	conformation::symmetry::SymmetricConformationOP symm_conf
		( setup_symmetric_conformation( pose.conformation(), symmdata, conf2pdb_chain(pose) ) );

	scoring::symmetry::SymmetricEnergiesOP symm_energy( new scoring::symmetry::SymmetricEnergies( pose.energies()) );

	pose.set_new_conformation( symm_conf );

	pose.set_new_energies_object( symm_energy );

	pose::PDBInfoOP pdb_info( new pose::PDBInfo( pose, true ) );
	core::pose::symmetry::make_symmetric_pdb_info( pose, pdb_info_src, pdb_info, keep_pdb_info_labels /* false */ );
	pose.pdb_info( pdb_info );

	debug_assert( is_symmetric( pose ) );

	//pose.conformation().detect_disulfides();
	core::pose::initialize_disulfide_bonds( pose );

	if ( option[ OptionKeys::symmetry::detect_bonds ] ) {
		pose.conformation().detect_bonds();
	}

	///Make GlycanTreeSet have symmetric info.
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		pose.conformation().clear_glycan_trees();
		pose.conformation().setup_glycan_trees();
	}
}

/// @details constructs a symmetric pose with a symmetric conformation and energies object from a monomeric pose
/// and symmetry definition file on command line. Requires the presence of a symmetry_definition file
void
make_symmetric_pose(
	pose::Pose & pose,
	std::string symmdef_file
)
{
	using namespace basic::options;

	conformation::symmetry::SymmData symmdata( pose.size(), pose.num_jump() );
	std::string symm_def = symmdef_file.length()==0 ? option[ OptionKeys::symmetry::symmetry_definition ] : symmdef_file;
	symmdata.read_symmetry_data_from_file(symm_def);
	make_symmetric_pose( pose, symmdata );
}

// @details make a symmetric pose assymmetric by instantiating new conformation and energies objects
void
make_asymmetric_pose(
	pose::Pose & pose
)
{
	conformation::ConformationOP conformation (
		new conformation::Conformation(static_cast< conformation::Conformation const & >( ( pose.conformation() ) ) ) );
	pose.set_new_conformation( conformation );

	scoring::EnergiesOP energies( new scoring::Energies( static_cast< scoring::Energies const & >( pose.energies() ) ) );

	energies->clear_energies();
	pose.set_new_energies_object( energies );

	///Make GlycanTreeSet have non-symmetric info
	if ( pose.conformation().contains_carbohydrate_residues() ) {
		pose.conformation().clear_glycan_trees();
		pose.conformation().setup_glycan_trees();
	}

	debug_assert( !is_symmetric( pose ) );
}

/// @brief extract the asu from a pose... unlike previous function symmetric clones are thrown out
/// @param[in]  pose_in            Symmetric input pose containing the asymmetric subunit of interest
/// @param[out] pose_out           Asymmetric subunit will be placed into this object
/// @param[in]  with_virtual_atoms If true, virtual atoms related to symmetry will be kept with the asymmetric subunit.
///                                If false, virtual atoms will be removed (default=true)
void extract_asymmetric_unit(
	core::pose::Pose const & pose_in,
	core::pose::Pose & pose_out,
	bool const with_virtual_atoms
)
{
	using core::conformation::Residue;
	using core::chemical::aa_vrt;
	using core::chemical::aa_h2o;
	using core::chemical::aa_unk;
	using namespace core::conformation::symmetry;

	if ( !is_symmetric( pose_in ) ) {
		pose_out = pose_in;
		return;
	}

	auto const & symm_conf (
		dynamic_cast<SymmetricConformation const & > ( pose_in.conformation() ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	for ( Size i=1; i<=symm_info->num_total_residues_without_pseudo(); i++ ) {
		if ( !symm_info->bb_is_independent(i) ) {
			continue;
		}
		Residue residue( pose_in.residue( i ) );

		//Special case carbohydrates - this is almost special-cased for branches, but we need the parent somehow.
		if ( residue.is_carbohydrate() && (i != 1) && (symm_conf.glycan_tree_set()->get_parent(i) != 0 ) ) {
			Size const parent_resnum = symm_conf.glycan_tree_set()->get_parent(i);
			Residue const & parent_res = symm_conf.residue( parent_resnum );

			//Connection between child (i) and parent is a branched connection.
			if ( parent_res.connected_residue_at_upper() != i ) {
				Size const anchor_atom_num = parent_res.connect_atom(residue);

				Size const child_atom_num = symm_conf.residue_type(i).lower_connect_atom();


				// Determine both connections...
				uint anchor_connection( 0 );
				uint connection( 0 );
				Size n_connections( parent_res.n_possible_residue_connections() );
				for ( Size ii = 1; ii <= n_connections; ++ii ) {
					if ( parent_res.residue_connect_atom_index( ii ) == anchor_atom_num ) {
						anchor_connection = ii;
						break;
					}
				}
				n_connections = residue.n_possible_residue_connections();

				for ( Size ii = 1; ii <= n_connections; ++ii ) {
					if ( residue.residue_connect_atom_index( ii ) == child_atom_num ) {
						connection = ii;
						break;
					}
				}

				if ( ( ! anchor_connection ) || ( ! connection ) ) {
					utility_exit_with_message( "Can't append by these atoms, "
						"since they do not correspond to connections on the Residues in question!" );

				}
				pose_out.append_residue_by_bond( residue, false, connection, parent_resnum, anchor_connection);
			} else {
				//Connection is normal.
				pose_out.append_residue_by_bond( residue, false );
			}
		} else if ( i>1 && residue.is_polymer_bonded( i-1 ) ) {
			pose_out.append_residue_by_bond( residue );
		} else {
			pose_out.append_residue_by_jump( residue, 1);
		}

		if ( i>1 ) {
			conformation::Residue const & prev_rsd( pose_in.residue( i-1 ) );
			if ( prev_rsd.is_upper_terminus() || residue.is_lower_terminus() || prev_rsd.chain() != residue.chain() ) {
				debug_assert( pose_out.size() == i );
				pose_out.conformation().insert_chain_ending( i-1 );
			}
		}
	}

	// disulfide topology has to be reconstructed
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	if ( pose_out.is_fullatom() ) {
		pose_out.conformation().detect_disulfides();
	}

	// add the parent VRT, set foldtree
	if ( with_virtual_atoms ) {
		core::pose::addVirtualResAsRoot(pose_out);
		core::kinematics::FoldTree const &f_in = core::conformation::symmetry::get_asymm_unit_fold_tree( pose_in.conformation() );
		pose_out.fold_tree( f_in );
	}

	pose::PDBInfoOP pdb_info( new pose::PDBInfo( pose_out, true ) );

	pose::PDBInfoCOP pdb_info_src ( pose_in.pdb_info() );
	if ( !pose_in.pdb_info() ) {
		pdb_info_src = utility::pointer::make_shared< pose::PDBInfo >( pose_in, true );
	}

	core::pose::symmetry::extract_asymmetric_unit_pdb_info( pose_in, pdb_info_src, pdb_info );
	pose_out.pdb_info( pdb_info );

	// reconstruct secondary structure
	debug_assert( pose_out.size() <= pose_in.size() );
	for ( core::Size resid=1; resid<=pose_out.size(); ++resid ) {
		pose_out.set_secstruct( resid, pose_in.secstruct( resid ) );
	}

	if ( pose_out.conformation().contains_carbohydrate_residues() ) {
		pose_out.conformation().clear_glycan_trees();
		pose_out.conformation().setup_glycan_trees();
	}

	runtime_assert( !core::pose::symmetry::is_symmetric( pose_out ) );
}


// @details make a asymmetric pose copying residues of original pose
core::pose::Pose
get_asymmetric_pose_copy_from_symmetric_pose(
	pose::Pose const & pose
)
{
	using core::conformation::Residue;
	using core::chemical::aa_vrt;
	using core::chemical::aa_unk;

	core::pose::Pose new_pose;

	bool jump_to_next = false;
	for ( Size i=1; i<=pose.size(); i++ ) {

		Residue residue( pose.residue( i ) );

		if ( residue.type().is_lower_terminus() ||
				residue.aa() == aa_unk || residue.aa() == aa_vrt || jump_to_next ) {

			if ( residue.aa() == aa_unk || residue.aa() == aa_vrt ) {
				jump_to_next = true;
			} else if ( jump_to_next ) {
				jump_to_next = false;
				///fpd ^^^ the problem is that the residue following the X should be connected by a jump as well.
				///     it should be of LOWER_TERMINUS variant type, but if not, we'll recover & spit out a warning for now.
				///     same thing for ligands???
				if ( ! residue.is_lower_terminus() ) {
					TR.Warning << "Residue following X, Z, or an upper terminus is _not_ a lower terminus type!  Continuing ..." << std::endl;
				}
			}
			new_pose.append_residue_by_jump( residue, 1, "", "", true ); // each time this happens, a new chain should be started

		} else {

			new_pose.append_residue_by_bond( residue, false );

			//fpd If res i is an upper terminus but (i+1) is not a lower terminus, the code exits on a failed assertion
			//fpd Don't let this happen; always jump in these cases
			if ( residue.type().is_upper_terminus(  ) ) jump_to_next = true;
		}
	}

	if ( new_pose.conformation().contains_carbohydrate_residues() ) {
		new_pose.conformation().clear_glycan_trees();
		new_pose.conformation().setup_glycan_trees();
	}

	return new_pose;
}

// @details make symmetric PDBIinfo
void
make_symmetric_pdb_info(
	pose::Pose const & pose,
	pose::PDBInfoOP pdb_info_src,
	pose::PDBInfoOP pdb_info_target,
	bool keep_pdb_info_labels /* false */
)
{
	using namespace core::conformation::symmetry;

	runtime_assert( is_symmetric( pose ) );
	runtime_assert( pdb_info_target->nres() == pose.size() );
	auto const & symm_conf( dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	// setup chainID mapping
	std::set< std::string > asymm_chains;
	std::set< std::string > used_chainIDs;

	std::map< std::pair< std::string,Size >, std::string > symmChainIDMap;

	// 1: find chainIDs in base unit, keep these the same
	for ( Size res=1; res <= pdb_info_src->nres(); ++res ) {
		std::string chn_id = pdb_info_src->chain( res );
		asymm_chains.insert( chn_id );
		used_chainIDs.insert( chn_id );
	}

	// 2: map (base chain id,clone#) => new chain id
	Size nclones = symm_info->num_bb_clones();
	for ( uint clone_i = 1; clone_i <= nclones; ++clone_i ) {
		for ( Size res=1; res <= pdb_info_src->nres(); ++res ) {
			std::string chn_id = pdb_info_src->chain( res );

			if ( symmChainIDMap.find( std::make_pair(chn_id,clone_i) ) == symmChainIDMap.end() ) {
				uint newchainidx = 0;
				std::string new_chain = "A";
				do {
					++newchainidx;
					new_chain = core::conformation::canonical_chain_letter_for_chain_number(newchainidx, /*extended*/ true );
				} while ( used_chainIDs.count(new_chain) >= 1 );

				symmChainIDMap[ std::make_pair(chn_id,clone_i) ] = new_chain;
			}
		}
	}

	// 3: now build the symmetrized pdbinfo
	for ( Size res=1; res <= pdb_info_src->nres(); ++res ) {
		// resids in scoring subunit
		int src_res = res;
		if ( !symm_info->bb_is_independent(res) ) {
			src_res = symm_info->bb_follows(res);
		}

		int tgt_res_id = pdb_info_src->number( res );
		pdb_info_target->number( src_res, tgt_res_id );

		// insertion codes
		char tgt_icode = pdb_info_src->icode( res );
		pdb_info_target->icode( src_res, tgt_icode );

		// chnids in scoring subunit
		std::string tgt_chn_id = pdb_info_src->chain( res );
		pdb_info_target->chain( src_res, tgt_chn_id );

		// symmetrize B's
		for ( Size atm=1; atm <= pdb_info_src->natoms(res); ++atm ) {
			core::Real tgt_b_atm = pdb_info_src->temperature( res, atm );
			pdb_info_target->temperature( src_res, atm, tgt_b_atm );
		}

		int clone_counter = 1;
		utility::vector1< Size > const & clones_i = symm_info->bb_clones( src_res );
		for ( core::Size i=1; i<=clones_i.size(); ++i ) {
			int clone_res = clones_i[i];
			pdb_info_target->number( clone_res, tgt_res_id );
			pdb_info_target->icode( clone_res, tgt_icode );

			std::string newchn_idx = symmChainIDMap[ std::make_pair(tgt_chn_id,clone_counter) ];

			if ( newchn_idx.empty() ) {
				pdb_info_target->chain( clone_res, " " ); // Didn't assign for some reason.
			} else {
				pdb_info_target->chain( clone_res, newchn_idx );
			}

			for ( Size atm=1; atm <= pdb_info_src->natoms(res) /*pose.residue(res).natoms()*/; ++atm ) {
				pdb_info_target->temperature( clone_res, atm, pdb_info_src->temperature( res, atm ) );
			}

			clone_counter++;
		}

		// keep info labels if desired
		if ( keep_pdb_info_labels ) {
			utility::vector1 < std::string > reslabels_for_res = pdb_info_src->get_reslabels(res);
			for ( auto const &reslabel : reslabels_for_res ) {
				pdb_info_target->add_reslabel( src_res, reslabel );
			}
		}
	}

	// rebuild pdb2pose
	pdb_info_target->rebuild_pdb2pose();

	// copy crystinfo
	pdb_info_target->set_crystinfo( pdb_info_src->crystinfo() );

	// copy header and remark lines
	if ( pdb_info_src->header_information() ) {
		pdb_info_target->header_information( utility::pointer::make_shared< io::HeaderInformation >(*pdb_info_src->header_information()));
	}
	pdb_info_target->remarks( pdb_info_src->remarks() );
}

void
extract_asymmetric_unit_pdb_info(
	pose::Pose const & pose,
	pose::PDBInfoCOP pdb_info_src,
	pose::PDBInfoOP pdb_info_target
)
{
	using namespace core::conformation::symmetry;

	if ( !is_symmetric(pose) ) {
		pdb_info_target = utility::pointer::make_shared< pose::PDBInfo >( *pdb_info_src );
		return;
	}

	auto const & symm_conf (
		dynamic_cast<SymmetricConformation const & > ( pose.conformation() ) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	Size nres = symm_info->get_nres_subunit();
	for ( Size res=1; res <= nres; ++res ) {
		int res_id = pdb_info_src->number( res );
		pdb_info_target->number( res, res_id );

		// insertion codes
		char icode = pdb_info_src->icode( res );
		pdb_info_target->icode( res, icode );

		std::string chn_id = pdb_info_src->chain( res );
		pdb_info_target->chain( res, chn_id );

		//std::cout << "Remap: " << res << " to " << res_id << chn_id << std::endl;

		// keep info labels
		utility::vector1 < std::string > reslabels_for_res = pdb_info_src->get_reslabels(res);
		for ( auto const &reslabel : reslabels_for_res ) {
			pdb_info_target->add_reslabel( res, reslabel );
		}

		// symmetrize B's
		// there may be a CYD<->CYS switch
		// account for this
		Size natoms = std::min( pdb_info_src->natoms(res), pdb_info_target->natoms(res) );
		for ( Size atm=1; atm <= natoms; ++atm ) {
			pdb_info_target->temperature( res, atm, pdb_info_src->temperature( res, atm ) );
		}
		for ( Size atm=natoms+1; atm <= pdb_info_target->natoms(res); ++atm ) {
			pdb_info_target->temperature( res, atm, pdb_info_src->temperature( res, natoms ) );
		}
	}
	// vrt
	if ( pdb_info_target->nres() > nres ) {
		pdb_info_target->number( nres+1, 1 );
		pdb_info_target->chain( nres+1, "z" );  //fpd  is this a problem???? should this be an "illegal" chainID instead?
	}

	// rebuild pdb2pose
	pdb_info_target->rebuild_pdb2pose();

	// copy crystinfo
	pdb_info_target->set_crystinfo( pdb_info_src->crystinfo() );

	// copy header and remark lines
	if ( pdb_info_src->header_information() ) {
		pdb_info_target->header_information( utility::pointer::make_shared< io::HeaderInformation >(*pdb_info_src->header_information() ));
	}
	pdb_info_target->remarks( pdb_info_src->remarks() );
}


// @details setting the movemap to only allow for symmetrical dofs.
// Jumps information is discarded and dof information in symmetry_info
// object is used instead
void
make_symmetric_movemap(
	pose::Pose const & pose,
	kinematics::MoveMap & movemap
)
{
	using namespace core::conformation::symmetry;

	runtime_assert( is_symmetric( pose ) );
	auto const & symm_conf (
		dynamic_cast<SymmetricConformation const & > ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	// allow only one subunit to change chi and torsions
	for ( Size i=1; i<= pose.conformation().size(); ++i ) {
		if ( !symm_info->bb_is_independent( i ) ) {
			movemap.set_bb ( i,false );
			movemap.set_chi( i, false );
			movemap.set_branches(i, false );
		}
	}

	// get number of non-VRT reses
	int numNonVrt = symm_info->num_total_residues_without_pseudo();

	// use the dof information in the symmetry_info object set allowed rb
	// dofs
	std::map< Size, SymDof > dofs ( symm_info->get_dofs() );

	std::map< Size, SymDof >::iterator it;
	auto it_begin = dofs.begin();
	auto it_end = dofs.end();

	// first find out whether we have any allowed jumps. we just need to
	// find one jump that is true to know that set_jump have been set to
	// true
	kinematics::MoveMapOP movemap_in( new kinematics::MoveMap( movemap ) );

	//fpd  only let THETA and D move  in master subunit
	if ( movemap.get( core::id::THETA ) ) {
		movemap.set( core::id::THETA , false );
		for ( Size i=1; i<= pose.conformation().size(); ++i ) {
			if ( symm_info->bb_is_independent( i ) ) {
				int const n_atoms( pose.residue(i).natoms() );
				for ( int j=1; j<=n_atoms; ++j ) {
					id::DOF_ID id( id::AtomID(j,i) ,id::THETA );
					if ( id.valid() && pose.has_dof(id) ) {
						movemap.set( id, true );
					}
				}
			}
		}
	}
	if ( movemap.get( core::id::D ) ) {
		movemap.set( core::id::D , false );
		for ( Size i=1; i<= pose.conformation().size(); ++i ) {
			if ( symm_info->bb_is_independent( i ) ) {
				int const n_atoms( pose.residue(i).natoms() );
				for ( int j=1; j<=n_atoms; ++j ) {
					id::DOF_ID id( id::AtomID(j,i) ,id::D );
					if ( id.valid() && pose.has_dof(id) ) {
						movemap.set( id, true );
					}
				}
			}
		}
	}

	///Switch any ON atoms that are not part of the master subunit to OFF for cartesian minimization
	if ( movemap.get_atoms().size() != 0 ) {
		std::map<id::AtomID, bool > const & atom_settings = movemap.get_atoms();
		for ( Size i=1; i<= pose.conformation().size(); ++i ) {
			if ( symm_info->bb_is_independent( i ) ) {
				continue;
			} else {
				int const n_atoms( pose.residue(i).natoms() );
				for ( int j=1; j<=n_atoms; ++j ) {
					id::AtomID atm = id::AtomID(j,i);
					if ( atom_settings.count(atm) ) {
						movemap.set_atom( atm, false );
					}
				}
			}
		}
	}

	// First set all jump to false
	movemap.set_jump(false);

	// Allow internal jumps to move (if allowed in movemap_in
	// ... allow it to move if it is in the controlling subunit
	// ... otherwise, it can't move
	for ( int jump_nbr = 1; jump_nbr <= (int)pose.num_jump(); ++jump_nbr ) {
		int upstream_resid = pose.fold_tree().upstream_jump_residue (jump_nbr);
		int downstream_resid = pose.fold_tree().downstream_jump_residue (jump_nbr);
		if ( upstream_resid <= numNonVrt && downstream_resid <= numNonVrt &&
				symm_info->bb_is_independent(upstream_resid) && symm_info->bb_is_independent(downstream_resid) ) {
			bool jump_move( movemap_in->get_jump( jump_nbr ) );
			for ( Size i = X_DOF; i <= Z_ANGLE_DOF; ++i ) {
				id::DOF_ID const & id( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,i)));
				if ( jump_move || movemap_in->get( id ) ) {
					movemap.set( id, true );
				}
			}
			//TR << "Setting internal jump (" << upstream_resid << "->" << downstream_resid
			//         << ") to " << movemap_in->get_jump(jump_nbr) << std::endl;
		}
	}

	// Then allow only dofs according to the dof information
	for ( it = it_begin; it != it_end; ++it ) {
		int jump_nbr ( (*it).first );
		SymDof dof( (*it).second );

		//bool jump_move( false );
		bool jump_move( movemap_in->get_jump( jump_nbr ) );
		for ( Size i = X_DOF; i <= Z_ANGLE_DOF; ++i ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,i)));
			jump_move |= movemap_in->get( id );
		}
		if ( !jump_move ) {
			continue;
		}

		if ( dof.allow_dof( X_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,1)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " X_DOF jump is allowed" << std::endl;

		}
		if ( dof.allow_dof( Y_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,2)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " Y_DOF jump is allowed" << std::endl;
		}
		if ( dof.allow_dof( Z_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,3)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " Z_DOF jump is allowed" << std::endl;
		}
		if ( dof.allow_dof( X_ANGLE_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,4)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " X_ANGLE_DOF jump is allowed" << std::endl;
		}
		if ( dof.allow_dof( Y_ANGLE_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,5)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " Y_ANGLE_DOF jump is allowed" << std::endl;
		}
		if ( dof.allow_dof( Z_ANGLE_DOF ) ) {
			id::DOF_ID const & id
				( pose.conformation().dof_id_from_torsion_id(id::TorsionID(jump_nbr,id::JUMP,6)));
			movemap.set( id, true );
			//std::cout << "make symmetric movemap: " << jump_nbr << " Z_ANGLE_DOF jump is allowed" << std::endl;
		}
	}
}

// @details find the anchor residue in the first subunit. This function assumes that the
// fold tree is rooted at the first jump from a virtual to its subunit
/*int
find_symmetric_basejump_anchor( pose::Pose & pose )
{
kinematics::FoldTree const & f  = pose.conformation().fold_tree();
kinematics::FoldTree::const_iterator it=f.begin();
debug_assert ( it->is_jump() );
return it->stop();
}*/
// @details find the anchor residue in the first subunit. This function assumes that
// there is only one one jump between a subunit and a virtual and that the subunit is
// folded downstream to the virtual
int
find_symmetric_basejump_anchor( pose::Pose & pose )
{
	kinematics::FoldTree const & f  = pose.conformation().fold_tree();
	for ( int i = 1; f.num_jump(); ++i ) {
		if ( pose.residue( f.downstream_jump_residue(i) ).name() != "VRT" &&
				pose.residue( f.upstream_jump_residue(i) ).name() == "VRT"  ) {
			return f.downstream_jump_residue(i);
		}
	}
	utility_exit_with_message( "No anchor residue is found..." );
	return 0;
}

// @ details move the anchor residue of a symmetric system. This function moves the anchors from
// virtual residues that define the coordinate system to the subunits. Useful during fragment insertion
// make sure to call rotate_to_x after this step so that the anchor moves correctly in the coordinate
// system defined by the virtual. The new anchor point is selected randomly or is selected as the point
// where the chains are closest together
void
find_new_symmetric_jump_residues( core::pose::Pose & pose )
{
	using namespace core::conformation::symmetry;

	runtime_assert( is_symmetric( pose ) );
	auto & symm_conf (
		dynamic_cast<SymmetricConformation & > ( pose.conformation()) );
	SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	Size nres_subunit ( symm_info->num_independent_residues() );
	Size anchor_start ( pose::symmetry::find_symmetric_basejump_anchor( pose ) );

	//Looking in central half of each segment -- pick a random point.
	Size anchor = static_cast<Size>( numeric::random::rg().uniform() * core::Size(nres_subunit/2) ) +
		(nres_subunit/4) + 1;

	if ( basic::options::option[ basic::options::OptionKeys::fold_and_dock::set_anchor_at_closest_point ] ) {
		//Find Closest point of contact.
		core::Real mindist = pose.residue(anchor).xyz("CEN").distance( pose.residue(anchor + nres_subunit).xyz("CEN") );
		for ( Size i = 1; i <= nres_subunit; i++ ) {
			Size const j = i + nres_subunit;
			core::Real dist = pose.residue(i).xyz("CEN").distance( pose.residue(j).xyz("CEN") );
			if ( dist < mindist ) {
				mindist = dist;
				anchor = i;
			}
		}
	}

	// Update the fold tree with the new jump points
	kinematics::FoldTree f ( pose.conformation().fold_tree() );

	// Setyp the lists of jumps and cuts
	Size num_jumps( f.num_jump() );
	Size num_cuts( f.num_cutpoint() );

	utility::vector1< Size > cuts_vector( f.cutpoints() );
	ObjexxFCL::FArray1D< Size > cuts( num_cuts );
	ObjexxFCL::FArray2D< Size > jumps( 2, num_jumps );

	// Initialize jumps
	for ( Size i = 1; i<= num_jumps; ++i ) {
		Size down ( f.downstream_jump_residue(i) );
		Size up ( f.upstream_jump_residue(i) );
		if ( down < up ) {
			jumps(1,i) = down;
			jumps(2,i) = up;
		} else {
			jumps(1,i) = up;
			jumps(2,i) = down;
		}
	}

	for ( Size i = 1; i<= num_cuts; ++i ) {
		cuts(i) = cuts_vector[i];
	}

	// anchor is always in first subunit. We need to convert it to the controlling subunit.
	if ( symm_info->bb_follows( anchor ) != 0 ) {
		//Size diff = symm_info->bb_follows( anchor_start ) - anchor_start;
		//anchor_start = symm_info->bb_follows( anchor_start );
		//anchor += diff;
		anchor = symm_info->bb_follows( anchor );
	}

	// This is the basejump
	Size root ( f.root() );
	Size const jump_number ( f.get_jump_that_builds_residue( anchor_start ) );
	Size residue_that_builds_anchor( f.upstream_jump_residue( jump_number ) );

	jumps(1, jump_number ) = anchor;
	jumps(2, jump_number ) = residue_that_builds_anchor;

	bool try_assert ( true );
	if ( symm_info->bb_follows( anchor_start ) != 0 ) {
		anchor_start = symm_info->bb_follows( anchor_start );
		try_assert = false;
	}

	for ( auto
			clone     = symm_info->bb_clones( anchor_start ).begin(),
			clone_end = symm_info->bb_clones( anchor_start ).end();
			clone != clone_end; ++clone ) {
		Size jump_clone ( f.get_jump_that_builds_residue( *clone ) );
		Size takeoff_pos ( f.upstream_jump_residue( jump_clone ) );
		Size new_anchor ( anchor - anchor_start + *clone );
		if ( try_assert ) runtime_assert( jumps(1,jump_clone) == *clone && jumps(2,jump_clone) == takeoff_pos );
		jumps(1, jump_clone ) = new_anchor;
		jumps(2, jump_clone ) = takeoff_pos;
		//std::cout<<"new_anchor "<<new_anchor<< " " <<anchor<<" "<<anchor_start<<" "<< takeoff_pos <<std::endl;
	}
	/* debug
	std::cout<<"cuts ";
	for ( Size i = 1; i<= num_cuts; ++i ) {
	std::cout<< cuts(i) << ' ';
	}
	std::cout<<std::endl;
	std::cout<<"jumps ";
	for ( Size i = 1; i<= num_jumps; ++i ) {
	std::cout<< " ( "<<jumps(1,i) << " , " << jumps(2,i) << " ) ";
	}
	std::cout<<std::endl;
	*/
	f.tree_from_jumps_and_cuts( pose.conformation().size(), num_jumps, jumps, cuts );
	f.reorder( root );
	pose.conformation().fold_tree( f );
}

// @details this function rotates a anchor residue to the x-axis defined by the virtual residue
// at the root of the fold tree. It is necessary that the anchor is placed in the coordinate system
// of the virtual residue so that rigid body movers and minimization code moves in the correct direction
void
rotate_anchor_to_x_axis( core::pose::Pose & pose ){

	//first anchor point -- assume its the first jump.
	kinematics::FoldTree f( pose.fold_tree() );
	int anchor ( find_symmetric_basejump_anchor( pose ) );
	int const jump_number ( f.get_jump_that_builds_residue( anchor ) );
	int coordsys_residue( f.upstream_jump_residue( jump_number ) );
	//Where is this stupid anchor point now?
	Vector anchor_pos ( pose.residue( anchor ).xyz("CA") );
	float theta = std::atan2(  anchor_pos(2), anchor_pos(1) );

	// we should make sure that the residue type is VRT
	runtime_assert( pose.residue( coordsys_residue ).name() == "VRT" );
	Vector x_axis( pose.residue( coordsys_residue ).xyz("X") - pose.residue( coordsys_residue ).xyz("ORIG") );
	Vector y_axis( pose.residue( coordsys_residue ).xyz("Y") - pose.residue( coordsys_residue ).xyz("ORIG") );

	Vector z_axis( x_axis );
	z_axis.normalize();
	Vector y_axis_nrml( y_axis );
	y_axis_nrml.normalize();
	z_axis = z_axis.cross( y_axis_nrml );
	z_axis.normalize();


	numeric::xyzMatrix< Real > coordsys_rot ( numeric::xyzMatrix< Real >::rows( x_axis, y_axis, z_axis ) );
	// make sure it is a rotation matrix
	runtime_assert( std::abs( coordsys_rot.det() - 1.0 ) < 1e-6 );

	// initialize a rotation matrix that rotates the anchor to cartesian x
	numeric::xyzMatrix< Real > z_rot = numeric::z_rotation_matrix_degrees(
		( -1.0 * numeric::conversions::to_degrees( theta ) ) );

	kinematics::Jump base_jump( pose.jump( jump_number ) );
	// rotate around absolute x
	Vector center(0,0,0);
	// Find the stub
	core::kinematics::Stub upstream_stub = pose.conformation().upstream_jump_stub( jump_number );
	// rotate to cartesian x, then rotate to the coordinate system defined by the virtual residue
	base_jump.rotation_by_matrix( upstream_stub, center, coordsys_rot.transpose()*z_rot );
	// set the jump
	pose.set_jump( jump_number, base_jump );
	// The convention is to place the subunit so that it has positive x values for slide moves to
	// go in the right direction. Silly!!!
	Vector new_anchor_pos ( pose.residue( anchor ).xyz("CA") );
	// std::cout << "before: " << anchor_pos(1) << " " << anchor_pos(2) << " " << anchor_pos(3) << std::endl
	//  <<  "after: " << new_anchor_pos(1) << " " << new_anchor_pos(2) << " " << new_anchor_pos(3) << std::endl;
	if ( new_anchor_pos(1) < 0 ) {
		upstream_stub = pose.conformation().upstream_jump_stub( jump_number );
		numeric::xyzMatrix< Real > twofold = numeric::z_rotation_matrix_degrees( 180.0 );
		base_jump.rotation_by_matrix( upstream_stub, center, twofold );
		pose.set_jump( jump_number, base_jump );
	}
}

/// @brief    Converts an asymmetric foldtree (f) with virtual root into a
///           symmetric foldtree compatible with symmetric pose (p)
/// @param    p - A symmetric pose
/// @param    f - An asymmetric foldtree. This foldtree MUST have a virtual root
/// @details  This function does not require the symm data
void
symmetrize_fold_tree( core::pose::Pose const &p, kinematics::FoldTree &f ) {
	core::conformation::symmetry::symmetrize_fold_tree( p.conformation(), f );
}

void
set_asymm_unit_fold_tree( core::pose::Pose &p, kinematics::FoldTree const &f) {
	core::conformation::symmetry::set_asymm_unit_fold_tree( p.conformation(), f );
}

// symmetry-aware version of FoldTree::partition_by_jump().  Accepts multiple jumps.
//  "floodfills" from root, not crossing any of the input jumps
void
partition_by_symm_jumps(
	utility::vector1< int > jump_numbers,
	core::kinematics::FoldTree const & ft,
	conformation::symmetry::SymmetryInfoCOP symm_info,
	ObjexxFCL::FArray1D_bool & partner1
) {
	using namespace core::kinematics;

	// expand jumps to include jump clones
	Size njumps = jump_numbers.size();
	for ( int i=1; i<=(int)njumps; ++i ) {
		utility::vector1< Size > clones_i = symm_info->jump_clones( jump_numbers[i] );
		for ( int j=1; j<= (int)clones_i.size(); ++j ) {
			jump_numbers.push_back( clones_i[j] );
		}
	}

	//int const pos1( ft.root() );
	int pos1=-1;
	for ( int i=1; i<= (int)symm_info->num_total_residues_without_pseudo(); ++i ) {
		if ( symm_info->bb_is_independent(i) ) {
			pos1 = i;
			break;
		}
	}
	debug_assert(pos1>0);

	partner1 = true;
	partner1( pos1 ) = false;

	bool new_member ( true );
	auto it_begin( ft.begin() );
	auto it_end  ( ft.end() );

	while ( new_member ) {     // keep adding new members
		new_member = false;
		for ( auto it = it_begin; it != it_end; ++it ) {
			if ( std::find ( jump_numbers.begin(), jump_numbers.end(), it->label() ) != jump_numbers.end() ) continue;

			int const start( std::min( it->start(), it->stop() ) );
			int const stop ( std::max( it->start(), it->stop() ) );
			if ( (partner1( start ) && !partner1( stop )) ||
					(partner1( stop ) && !partner1( start )) ) {
				new_member = true;
				if ( it->is_polymer() ) {
					// all the residues
					for ( int i=start; i<= stop; ++i ) {
						partner1( i ) = false;
					}
				} else {
					// just the vertices
					partner1( start ) = false;
					partner1( stop ) = false;
				}
			}
		}
	}
}


// find symmetry axis
// returns <0,0,0> if:
//    system is nonsymmetric
//    dimer symmetry
numeric::xyzVector< core::Real >
get_symm_axis( core::pose::Pose & pose ) {
	if ( !is_symmetric( pose ) ) return numeric::xyzVector< core::Real >(0,0,0);

	auto const & symm_conf (
		dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );

	// find first independent residue
	core::Size base = 1;
	for ( Size i=1; i <=symm_info->num_total_residues_with_pseudo(); ++i ) {
		if ( symm_info->bb_is_independent(i) ) {
			base = i;
			break;
		}
	}

	// find clones
	utility::vector1< core::Size > base_clones = symm_info->bb_clones( base );
	if ( base_clones.size() <= 1 ) {
		return numeric::xyzVector< core::Real >(0,0,0);
	}

	// base_clones >= 2
	numeric::xyzVector< core::Real > X = pose.residue(base_clones[1]).xyz(1) - pose.residue(base).xyz(1);
	numeric::xyzVector< core::Real > Y = pose.residue(base_clones[2]).xyz(1) - pose.residue(base).xyz(1);
	return X.cross( Y ).normalize();
}

// Figure out which chains touch chain A, and return those chains
core::pose::Pose
get_buildingblock_and_neighbor_subs (Pose const &pose_in, utility::vector1<Size> intra_subs) {

	if ( is_multicomponent(pose_in) ) {
		utility_exit_with_message("core::pose::symmetry::get_buildingblock_and_neighbor_subs is only for use with singlecomponent symmetries.");
	}

	//fpd we need to first score the pose
	Pose pose = pose_in;
	core::scoring::ScoreFunction sc_rep;
	sc_rep.set_weight( core::scoring::fa_atr, 1.0 );
	sc_rep( pose );

	Pose sub_pose;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	Size nres_monomer = symm_info->num_independent_residues();

	// add all base subunits
	for ( Size i=1; i<=symm_info->subunits(); ++i ) {
		if ( std::find(intra_subs.begin(), intra_subs.end(), i) == intra_subs.end() ) continue;
		Size start = (i-1)*nres_monomer;
		sub_pose.append_residue_by_jump(pose.residue(start+1),sub_pose.size());
		for ( Size ir=2; ir<=nres_monomer; ir++ ) {
			sub_pose.append_residue_by_bond(pose.residue(ir+start));
		}
	}

	for ( Size i=1; i<=symm_info->subunits(); ++i ) {
		if ( std::find(intra_subs.begin(), intra_subs.end(), i) != intra_subs.end() ) continue;
		bool contact = false;
		Size start = (i-1)*nres_monomer;
		for ( Size ir=1; ir<=nres_monomer; ir++ ) {
			if ( pose.energies().residue_total_energies(ir+start)[core::scoring::fa_atr] < 0 ) {
				contact = true;
				break;
			}
		}
		if ( contact ) {
			sub_pose.append_residue_by_jump(pose.residue(start+1),sub_pose.size());
			for ( Size ir=2; ir<=nres_monomer; ir++ ) {
				sub_pose.append_residue_by_bond(pose.residue(ir+start));
			}
		}
	}

	return sub_pose;
}

core::pose::Pose
get_subpose(Pose const &pose, utility::vector1<std::pair<Size,std::string>> const & subs) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	Pose sub_pose;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	std::pair<Size,std::string> prev_sub;
	bool start = true;
	for ( Size i=1; i<=symm_info->num_total_residues_without_pseudo(); i++ ) {
		if ( find(subs.begin(),subs.end(), get_resnum_to_subunit_component(pose,i)) == subs.end() ) continue;
		if ( start || (get_resnum_to_subunit_component(pose,i) != prev_sub) ) {
			sub_pose.append_residue_by_jump(pose.residue(i),sub_pose.size());
		} else {
			sub_pose.append_residue_by_bond(pose.residue(i));
		}
		start = false;
		prev_sub = get_resnum_to_subunit_component(pose,i);
	}
	return sub_pose;
}

utility::vector1<Size>
get_resis(Pose const &pose, utility::vector1<std::pair<Size,std::string>> const & subs ) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1<Size> resis;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	for ( Size i=1; i<=symm_info->num_total_residues_without_pseudo(); i++ ) {
		if ( find(subs.begin(),subs.end(), get_resnum_to_subunit_component(pose,i)) == subs.end() ) continue;
		resis.push_back(i);
	}
	return resis;
}

// Contacting intracomponent neighboring subunits
utility::vector1< std::pair<Size,std::string> >
get_full_intracomponent_and_neighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	if ( is_singlecomponent(pose) ) {
		utility_exit_with_message("core::pose::symmetry::get_full_intracomponent_and_neighbor_subs is only for use with multicomponent symmetries.");
	}

	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	int jump = sym_dof_jump_num( pose, sym_dof_name );
	Real const contact_dist_sq = contact_dist * contact_dist;

	ObjexxFCL::FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump, is_upstream );
	for ( Size i=1; i<=symm_info->num_independent_residues(); ++i ) {
		if ( is_upstream(i) ) continue;
	}
	utility::vector1<std::pair<Size,std::string>> subs = get_full_intracomponent_subs(pose,sym_dof_name);
	for ( Size i=1; i<=symm_info->num_total_residues_without_pseudo(); ++i ) {
		if ( is_upstream(i) ) continue;
		for ( Size j=1; j<=symm_info->num_total_residues_without_pseudo(); j++ ) {
			if ( !is_upstream(j) ) continue;
			if ( find(subs.begin(),subs.end(), get_resnum_to_subunit_component(pose,j) )!=subs.end() ) continue;
			std::string atom_i = (pose.residue(i).name3() == "GLY") ? "CA" : "CB";
			std::string atom_j = (pose.residue(j).name3() == "GLY") ? "CA" : "CB";
			if ( pose.residue(i).xyz(atom_i).distance_squared(pose.residue(j).xyz(atom_j)) <= contact_dist_sq ) {
				subs.push_back(get_resnum_to_subunit_component(pose,j));
			}
		}
	}
	return subs;
}

core::pose::Pose
get_full_intracomponent_and_neighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose, get_full_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist));
}

utility::vector1<Size>
get_full_intracomponent_and_neighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose, get_full_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist));
}

core::pose::Pose
get_full_intracomponent_subpose(Pose const &pose, std::string const & sym_dof_name) {
	return get_subpose(pose, get_full_intracomponent_subs(pose, sym_dof_name));
}

utility::vector1<Size>
get_full_intracomponent_resis(Pose const &pose, std::string const & sym_dof_name) {
	return get_resis(pose, get_full_intracomponent_subs(pose, sym_dof_name));
}

// Contacting neighbor subunits of the full component building block controlled by the symdof
utility::vector1<std::pair<Size,std::string>>
get_full_intracomponent_neighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1<std::pair<Size,std::string>> subs = get_full_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist);
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	// Loop through subunit chains and only keep those that are not the primary subunits of the component(s) controlled by the specified sym_dof.
	utility::vector1<std::pair<Size,std::string>> full_intracomponent_subs = get_full_intracomponent_subs(pose,sym_dof_name);
	utility::vector1<std::pair<Size,std::string>> sub_subs;
	for ( Size i=1; i<=subs.size(); i++ ) {
		if ( find(full_intracomponent_subs.begin(),full_intracomponent_subs.end(),subs[i])!=full_intracomponent_subs.end() ) continue;
		sub_subs.push_back(subs[i]);
	}
	return sub_subs;
}

core::pose::Pose
get_full_intracomponent_neighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose,get_full_intracomponent_neighbor_subs(pose,sym_dof_name,contact_dist));
}

utility::vector1<Size>
get_full_intracomponent_neighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose,get_full_intracomponent_neighbor_subs(pose,sym_dof_name,contact_dist));
}


// Are there intracomponent contacts for the specified symdof?
bool
intracomponent_contact(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	if ( is_singlecomponent(pose) ) {
		utility_exit_with_message("core::pose::symmetry::get_intracomponent_and_neighbor_subs is only for use with multicomponent symmetries.");
	}

	Size monomer_lower_bound = 0;
	Size monomer_upper_bound = 0;
	utility::vector1<std::string> subs;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	int jump = sym_dof_jump_num( pose, sym_dof_name );
	Real const contact_dist_sq = contact_dist * contact_dist;

	ObjexxFCL::FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump, is_upstream );
	bool start = true;
	for ( Size i=1; i<=symm_info->num_independent_residues(); ++i ) {
		if ( is_upstream(i) ) continue;
		if ( start ) monomer_lower_bound = i;
		start = false;
		monomer_upper_bound = i;
	}
	bool contact = false;
	utility::vector1< std::pair<Size,std::string>> full_intracomponent_subs(get_full_intracomponent_subs(pose, sym_dof_name));
	for ( Size i=monomer_lower_bound; i<=monomer_upper_bound; ++i ) {
		for ( Size j=1; j<=symm_info->num_total_residues_without_pseudo(); j++ ) {
			if ( find(full_intracomponent_subs.begin(),full_intracomponent_subs.end(), get_resnum_to_subunit_component(pose,j) )!=full_intracomponent_subs.end() ) continue;
			if ( get_component_of_residue(pose,i) != get_component_of_residue(pose,j) ) continue;
			std::string atom_i = (pose.residue(i).name3() == "GLY") ? "CA" : "CB";
			std::string atom_j = (pose.residue(j).name3() == "GLY") ? "CA" : "CB";
			if ( pose.residue(i).xyz(atom_i).distance_squared(pose.residue(j).xyz(atom_j)) <= contact_dist_sq ) {
				contact = true;
			}
			if ( contact ) break; // j loop
		}
		if ( contact ) {
			TR << "Intracontact detected for the component controlled by sym_dof: " << sym_dof_name << std::endl;
			break; // i loop
		}
	}
	return contact;
}

// Contacting intracomponent and neighboring subunits
utility::vector1< std::pair<Size,std::string> >
get_intracomponent_and_neighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	if ( is_singlecomponent(pose) ) {
		utility_exit_with_message("core::pose::symmetry::get_intracomponent_and_neighbor_subs is only for use with multicomponent symmetries.");
	}

	Size monomer_lower_bound = 0;
	Size monomer_upper_bound = 0;
	utility::vector1< std::pair<Size,std::string> > subs;
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);
	int jump = sym_dof_jump_num( pose, sym_dof_name );
	Real const contact_dist_sq = contact_dist * contact_dist;

	ObjexxFCL::FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump, is_upstream );
	bool start = true;
	for ( Size i=1; i<=symm_info->num_independent_residues(); ++i ) {
		if ( is_upstream(i) ) continue;
		if ( start ) monomer_lower_bound = i;
		start = false;
		monomer_upper_bound = i;
	}
	subs.push_back(get_resnum_to_subunit_component(pose,monomer_lower_bound));
	for ( Size i=monomer_lower_bound; i<=monomer_upper_bound; ++i ) {
		for ( Size j=1; j<=symm_info->num_total_residues_without_pseudo(); j++ ) {
			if ( find(subs.begin(),subs.end(), get_resnum_to_subunit_component(pose,j) )!=subs.end() ) continue;
			std::string atom_i = (pose.residue(i).name3() == "GLY") ? "CA" : "CB";
			std::string atom_j = (pose.residue(j).name3() == "GLY") ? "CA" : "CB";
			if ( pose.residue(i).xyz(atom_i).distance_squared(pose.residue(j).xyz(atom_j)) <= contact_dist_sq ) {
				subs.push_back(get_resnum_to_subunit_component(pose,j));
			}
		}
	}
	return subs;
}

core::pose::Pose
get_intracomponent_and_neighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose, get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist));
}

utility::vector1<Size>
get_intracomponent_and_neighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose, get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist));
}

// Contacting intracomponent subunits only
utility::vector1<std::pair<Size,std::string>>
get_intracomponent_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1<std::pair<Size,std::string>> subs = get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist);
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	// Loop through subunit chains and only keep those that are primary subunits of the component(s) controlled by the specified sym_dof.
	utility::vector1<std::pair<Size,std::string>> full_intracomponent_subs = get_full_intracomponent_subs(pose,sym_dof_name);
	utility::vector1<std::pair<Size,std::string>> sub_subs;
	for ( Size i=1; i<=subs.size(); i++ ) {
		if ( find(full_intracomponent_subs.begin(),full_intracomponent_subs.end(),subs[i])==full_intracomponent_subs.end() ) continue;
		sub_subs.push_back(subs[i]);
	}

	return sub_subs;
}

core::pose::Pose
get_intracomponent_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose,get_intracomponent_subs(pose,sym_dof_name,contact_dist));
}

utility::vector1<Size>
get_intracomponent_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose,get_intracomponent_subs(pose,sym_dof_name,contact_dist));
}

// Contacting neighbor subunits only
utility::vector1<std::pair<Size,std::string>>
get_neighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1<std::pair<Size,std::string>> subs = get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist);
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	// Loop through subunit chains and only keep those that are not the primary subunits of the component(s) controlled by the specified sym_dof.
	utility::vector1<std::pair<Size,std::string>> full_intracomponent_subs = get_full_intracomponent_subs(pose,sym_dof_name);
	utility::vector1<std::pair<Size,std::string>> sub_subs;
	for ( Size i=1; i<=subs.size(); i++ ) {
		if ( find(full_intracomponent_subs.begin(),full_intracomponent_subs.end(),subs[i])!=full_intracomponent_subs.end() ) continue;
		sub_subs.push_back(subs[i]);
	}
	return sub_subs;
}

core::pose::Pose
get_neighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose,get_neighbor_subs(pose,sym_dof_name,contact_dist));
}

utility::vector1<Size>
get_neighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose,get_neighbor_subs(pose,sym_dof_name,contact_dist));
}

// Contacting intracomponents and intraneighbors
utility::vector1< std::pair<Size,std::string> >
get_intracomponent_and_intraneighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1< std::pair<Size,std::string> > subs = get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist);
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	// Loop through subunit chains and only keep those that are the primary subunits of the component(s) controlled by the specified sym_dof or neighbors that of the same component.
	utility::vector1< std::pair<Size,std::string> > sub_subs;
	utility::vector1<std::string> components = get_jump_name_to_components(pose,sym_dof_name);
	for ( Size i=1; i<=subs.size(); i++ ) {
		if ( find(components.begin(),components.end(),subs[i].second) ==components.end() ) continue;
		sub_subs.push_back(subs[i]);
	}
	return sub_subs;
}

core::pose::Pose
get_intracomponent_and_intraneighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose,get_intracomponent_and_intraneighbor_subs(pose,sym_dof_name,contact_dist));
}

utility::vector1<Size>
get_intracomponent_and_intraneighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose,get_intracomponent_and_intraneighbor_subs(pose,sym_dof_name,contact_dist));
}

// Contacting intracomponents and interneighbors
utility::vector1<std::pair<Size,std::string>>
get_intracomponent_and_interneighbor_subs(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	using namespace basic;
	using namespace core::conformation::symmetry;
	using namespace core::pose::symmetry;

	utility::vector1<std::pair<Size,std::string>> subs = get_intracomponent_and_neighbor_subs(pose, sym_dof_name, contact_dist);
	core::conformation::symmetry::SymmetryInfoCOP symm_info = core::pose::symmetry::symmetry_info(pose);

	// Loop through subunit chains and only keep those that are the primary subunits of the component(s) controlled by the specified sym_dof or neighbors that of the not of the same component.
	utility::vector1<std::pair<Size,std::string>> full_intracomponent_subs = get_full_intracomponent_subs(pose,sym_dof_name);
	utility::vector1<std::pair<Size,std::string>> sub_subs;
	utility::vector1<std::string> components = get_jump_name_to_components(pose,sym_dof_name);
	for ( Size i=1; i<=subs.size(); i++ ) {
		if ( find(full_intracomponent_subs.begin(),full_intracomponent_subs.end(),subs[i])==full_intracomponent_subs.end() && find(components.begin(),components.end(),subs[i].second)!=components.end() ) continue;
		sub_subs.push_back(subs[i]);
	}
	return sub_subs;
}

core::pose::Pose
get_intracomponent_and_interneighbor_subpose(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_subpose(pose,get_intracomponent_and_interneighbor_subs(pose,sym_dof_name,contact_dist));
}

utility::vector1<Size>
get_intracomponent_and_interneighbor_resis(Pose const &pose, std::string const & sym_dof_name, Real contact_dist) {
	return get_resis(pose,get_intracomponent_and_interneighbor_subs(pose,sym_dof_name,contact_dist));
}

void
make_residue_mask_symmetric( core::pose::Pose const &p, utility::vector1< bool > & msk ) {
	if ( !is_symmetric( p ) ) return;

	auto const & symm_conf (
		dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( p.conformation() ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	// Size nres_subunit ( symm_info->num_independent_residues() );
	// Size nsubunits ( symm_info->subunits() );

	runtime_assert( symm_info->num_total_residues_with_pseudo() == msk.size() );

	for ( core::Size i=1; i<= msk.size(); ++i ) {
		bool msk_i=msk[i];
		if ( msk_i && !symm_info->fa_is_independent( i ) ) {
			msk[i] = false;
			msk[ symm_info->bb_follows( i ) ] = true;
		}
	}
}

// Remove internal subunit jumps and cuts and create a fold tree that only has the symmetric
// jump framework
kinematics::FoldTree
sealed_symmetric_fold_tree( core::pose::Pose & pose ) {
	if ( !is_symmetric( pose ) ) {
		TR.Error << "sealed_symmetric_fold_tree called with assymetric fold tree. Return FoldTree" << std::endl;
		return pose.fold_tree();
	}

	kinematics::FoldTree f_orig = pose.fold_tree();
	kinematics::FoldTree f = f_orig;
	//TR.Error << "sealed_symmetric_fold_tree called with " << f_orig << std::endl;

	auto const & symm_conf (
		dynamic_cast<conformation::symmetry::SymmetricConformation const & > ( pose.conformation() ) );
	conformation::symmetry::SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
	//Size nres_subunit ( symm_info->num_independent_residues() );

	// mjo commenting out 'nsubunits' because it is not used and casus a warning
	//Size nsubunits ( symm_info->subunits() );
	//Size num_nonvrt = symm_info->num_total_residues_without_pseudo();

	// 1 - Get monomer jumps, cuts, and anchor

	// mjo commenting out 'new_anchor' because it is not used and causes a warning
	//Size new_anchor = 0;
	utility::vector1< Size > new_cuts;
	utility::vector1< std::pair<Size,Size> > new_jumps;

	utility::vector1< int > cuts_vector( f_orig.cutpoints() );

	// 2 - Get symmetic jumps cuts and anchor
	// inter-VRT jumps
	for ( Size i = 1; i<= f_orig.num_jump(); ++i ) {
		Size down ( f_orig.downstream_jump_residue(i) );
		Size up ( f_orig.upstream_jump_residue(i) );
		/*
		// connections between VRTs are unchanged
		if ( up > num_nonvrt && down > num_nonvrt) {
		new_jumps.push_back( std::pair<Size,Size>(up,down) );
		}
		// jumps to anchor
		if ( up > num_nonvrt && down <= num_nonvrt) {
		new_jumps.push_back( std::pair<Size,Size>( up, down ) );
		}
		// jumps from anchor
		if ( up <= num_nonvrt && down > num_nonvrt) {
		new_jumps.push_back( std::pair<Size,Size>( up, down ) );
		}
		*/
		// all virts are of type ligand, so just use it to catch all the jumps
		if ( pose.conformation().residue(up).is_ligand() == true || pose.conformation().residue(down).is_ligand() == true ) {
			new_jumps.push_back( std::pair<Size,Size>( up, down) );
		}

	}

	// cuts
	cuts_vector = f_orig.cutpoints();
	for ( Size i = 1; i<=cuts_vector.size(); ++i ) {
		// cuts between vrts and jumps between chains
		//  if ( cuts_vector[i] >= (int) num_nonvrt || cuts_vector[i]%nres_subunit == 0 ) // not good if ligand is present.

		//if terminus
		if ( pose.conformation().residue(cuts_vector[i]).is_terminus() ) {
			new_cuts.push_back( cuts_vector[i] );
		}

		//if virtural or ligand
		if ( pose.conformation().residue(cuts_vector[i]).is_ligand() ) {
			new_cuts.push_back( cuts_vector[i] );
		}
	}


	// scan for ligands

	// 3 - put them in FArrays
	ObjexxFCL::FArray1D< Size > cuts( new_cuts.size() );
	ObjexxFCL::FArray2D< Size > jumps( 2, new_jumps.size() );
	// Initialize jumps
	for ( Size i = 1; i<= new_jumps.size(); ++i ) {
		jumps(1,i) = std::min( (int)new_jumps[i].first, (int)new_jumps[i].second);
		jumps(2,i) = std::max( (int)new_jumps[i].first, (int)new_jumps[i].second);
		// DEBUG -- PRINT JUMPS AND CUTS
		//TR.Error << " jump " << i << " : " << jumps(1,i) << " , " << jumps(2,i) << std::endl;
	}
	for ( Size i = 1; i<= new_cuts.size(); ++i ) {
		cuts(i) = (int)new_cuts[i];
		// TR.Error << " cut " << i << " : " << cuts(i) << std::endl;
	}

	// 4 make the sealed foldtree
	f.clear();
	f.tree_from_jumps_and_cuts( pose.size(), new_jumps.size(), jumps, cuts, f_orig.root(), false );
	return f;
}

// @brief given a symmetric pose and a jump number, get_sym_aware_jump_num
//   translates the jump number into the corresponding SymDof so that the
//   symmetric pose can moved moved appropriately.
int
get_sym_aware_jump_num ( core::pose::Pose const & pose, core::Size jump_num ) {
	using namespace core::conformation::symmetry;
	Size sym_jump = jump_num;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		SymmetryInfoCOP sym_info = core::pose::symmetry::symmetry_info(pose);
		std::map<Size,SymDof> dofs = sym_info->get_dofs();
		sym_jump = 0;
		Size jump_counter = 0;

		for ( auto & dof : dofs ) {
			//fpd  if slide moves are not allowed on this jump, then skip it
			if ( !dof.second.allow_dof(1) && !dof.second.allow_dof(2) && !dof.second.allow_dof(3) ) continue;
			if ( ++jump_counter == (Size)jump_num ) {
				sym_jump = dof.first;
				break;
			}
		}
		if ( sym_jump == 0 ) {
			TR.Error << "Failed to find sym_dof with index " << jump_num << std::endl;
			utility_exit_with_message("Symmetric slide DOF "" not found!");
		}
	}
	return sym_jump;
} // get_symdof_from_jump_num


utility::vector1<std::string>
sym_dof_names(core::pose::Pose const & pose) {
	using namespace core::conformation::symmetry;
	utility::vector1<std::string> names;
	SymmetryInfoCOP syminfo = core::pose::symmetry::symmetry_info(pose);
	std::map<Size,SymDof> dofs = syminfo->get_dofs();
	for ( auto & dof : dofs ) {
		names.push_back(syminfo->get_jump_name(dof.first));
	}
	return names;
}

int
sym_dof_jump_num(core::pose::Pose const & pose, std::string const & jname){
	return core::pose::symmetry::symmetry_info(pose)->get_jump_num(jname);
}

std::string
jump_num_sym_dof(core::pose::Pose const & pose, Size const & jnum){
	return core::pose::symmetry::symmetry_info(pose)->get_jump_name(jnum);
}

utility::vector1<Size>
get_symdof_subunits(core::pose::Pose const & pose, std::string const & jname){
	using namespace core::conformation::symmetry;
	utility::vector1<Size> subs;
	if ( core::pose::symmetry::is_symmetric(pose) ) {
		SymmetryInfo const & sym_info = *core::pose::symmetry::symmetry_info(pose);
		int jnum = sym_dof_jump_num(pose,jname);
		ObjexxFCL::FArray1D_bool is_upstream( pose.size(), false );
		pose.fold_tree().partition_by_jump( jnum, is_upstream );
		for ( Size i = 1; i <= sym_info.num_total_residues_without_pseudo(); i+=sym_info.get_nres_subunit() ) {
			if ( is_upstream(i) ) {
				subs.push_back( sym_info.subunit_index(i) );
				std::cout << subs.back() << std::endl;
			}
		}
	} else {
		utility_exit_with_message("pose not symmetric!");
	}
	return subs;
}

utility::vector1<std::string> symmetric_components(core::pose::Pose const & pose){
	return symmetry_info(pose)->get_components();
}

bool is_singlecomponent(core::pose::Pose const & pose){
	return symmetry_info(pose)->get_num_components()==1;
}
bool is_multicomponent(core::pose::Pose const & pose){
	return symmetry_info(pose)->get_num_components()>=2;
}

Size get_component_lower_bound(core::pose::Pose const & pose, std::string const & c){
	return symmetry_info(pose)->get_component_lower_bound(c);
}
Size get_component_upper_bound(core::pose::Pose const & pose, std::string const & c){
	return symmetry_info(pose)->get_component_upper_bound(c);
}
std::string get_component_of_residue(core::pose::Pose const & pose, Size ir){
	return symmetry_info(pose)->get_component_of_residue(ir);
}
std::string get_subunit_name_to_component(core::pose::Pose const & pose, std::string const & vname){
	return symmetry_info(pose)->get_subunit_name_to_component(vname);
}
utility::vector1<std::string> const & get_jump_name_to_components(core::pose::Pose const & pose, std::string const & jname){
	return symmetry_info(pose)->get_jump_name_to_components(jname);
}
utility::vector1<Size> const & get_jump_name_to_subunits(core::pose::Pose const & pose, std::string const & jname){
	return symmetry_info(pose)->get_jump_name_to_subunits(jname);
}

std::pair< Size, std::string> get_resnum_to_subunit_component(core::pose::Pose const & pose, Size const & resnum) {
	return std::make_pair( symmetry_info(pose)->subunit_index(resnum), symmetry_info(pose)->get_component_of_residue(resnum) );
}

utility::vector1< std::pair<Size,std::string> >
get_full_intracomponent_subs(core::pose::Pose const & pose, std::string const & jname){
	utility::vector1<std::string> components = get_jump_name_to_components(pose,jname);
	utility::vector1<Size> intrasubs = get_jump_name_to_subunits(pose,jname);
	utility::vector1< std::pair<Size,std::string> > full_intracomponent_subs;
	for ( Size c=1; c<=components.size(); c++ ) {
		for ( Size s=1; s<=intrasubs.size(); s++ ) {
			full_intracomponent_subs.push_back( std::make_pair( intrasubs[s], components[c] ) );
		}
	}
	return full_intracomponent_subs;
}

} // symmetry
} // pose
} // core
