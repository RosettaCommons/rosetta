// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/SWM_RMSD_Energy.cc
/// @brief  Hack to force SWM to generate native-like conformations
/// @author Arvind Kannan


// Unit headers
#include <core/scoring/methods/SWM_RMSD_Energy.hh>
#include <core/scoring/methods/SWM_RMSD_EnergyCreator.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/id/AtomID.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/pose/rna/RNA_Util.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/io/pdb/file_data.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

// C++
#include <basic/Tracer.hh>
#include <map>

static basic::Tracer TR("core.scoring.methods.SWM_RMSD_Energy");

/////////////////////////////////////////////////////////////////////////////////////
//
// Created in order to generate native-like conformations for score comparisons with
// the results of SWM calculations.
//
/////////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {
namespace methods {


/// @details This must return a fresh instance of the SWM_RMSD_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
SWM_RMSD_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new SWM_RMSD_Energy;
}

ScoreTypes
SWM_RMSD_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( swm_rmsd );
	return sts;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
SWM_RMSD_Energy::add_to_atom_id_map_after_checks( std::map< id::AtomID, id::AtomID> & atom_id_map,
								std::string const & atom_name,
								Size const & n1, Size const & n2,
								pose::Pose const & pose1, pose::Pose const & pose2 ) const {
	
	using namespace core::id;
	
	runtime_assert ( n1 >= 1 && n1 <= pose1.total_residue() );
	runtime_assert ( n2 >= 1 && n2 <= pose2.total_residue() );
	runtime_assert( pose1.residue_type( n1 ).aa() == pose2.residue_type( n2 ).aa() );
	
	if ( ! pose1.residue_type( n1 ).has( atom_name ) ) return;
	if ( ! pose2.residue_type( n2 ).has( atom_name ) ) return;
	
	Size const idx1 = pose1.residue_type( n1 ).atom_index( atom_name );
	Size const idx2 = pose2.residue_type( n2 ).atom_index( atom_name );
	
	if ( pose1.residue_type( n1 ).is_virtual( idx1 ) ) return;
	if ( pose2.residue_type( n2 ).is_virtual( idx2 ) ) return;
	
	atom_id_map[ AtomID( idx1, n1 ) ] = AtomID( idx2, n2 );
	
}
	
bool
SWM_RMSD_Energy::mutate_position( pose::Pose & pose, Size const i, char const & new_seq ) const {
	
	using namespace core::conformation;
	using namespace core::chemical;
	
	if ( new_seq == pose.sequence()[i-1] ) return false;
	
	ResidueTypeSet const & rsd_set = pose.residue( i ).residue_type_set();
	
	ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( new_seq ).match_variants( pose.residue(i).type() ).select( rsd_set )[1] );
	ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );
	
	Real const save_chi = pose.chi(i);
	pose.replace_residue( i, *new_rsd, false );
	pose.set_chi( i, save_chi );
	
	return true;
}
	
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
SWM_RMSD_Energy::superimpose_at_fixed_res( pose::Pose & pose, pose::Pose const & native_pose,
						 Real & rmsd, Size & natoms_rmsd ) const {
	
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::pose;
	using namespace core::pose::full_model_info;
	using namespace core::scoring;
	
	Pose native_pose_local = native_pose; // local working copy, mutated in cases where nucleotides have been designed ('n')
	
	// first need to slice up native_pose to match residues in actual pose.
	// define atoms over which to compute RMSD, using rmsd_res.
	FullModelInfo const & full_model_info = const_full_model_info( pose );
	utility::vector1< Size > const & res_list = get_res_list_from_full_model_info_const( pose );
	utility::vector1< Size > const & fixed_domain_map = full_model_info.fixed_domain_map();
	std::string const full_sequence = full_model_info.full_sequence();
	
	// following needs to be updated.
	utility::vector1< Size > const rmsd_res = full_model_info.moving_res_in_full_model();
	
	utility::vector1< Size > calc_rms_res;
	for ( Size n = 1; n <= pose.total_residue(); n++ ){
		if ( rmsd_res.has_value( res_list[ n ] ) ) {
			calc_rms_res.push_back( n );
			
			char const pose_nt = pose.sequence()[ n-1 ];
			if ( full_sequence[ res_list[ n ] - 1 ] == 'n' ){
				mutate_position( native_pose_local, res_list[ n ], pose_nt );
			} else {
				runtime_assert( full_sequence[ res_list[ n ] - 1 ] == pose_nt);
			}
			runtime_assert( native_pose_local.sequence()[ res_list[ n ] - 1] == pose_nt );
		}
	}
	
	std::map< AtomID, AtomID > calc_rms_atom_id_map;
	
	for ( Size k = 1; k <= calc_rms_res.size(); k++ ){
		Size const n = calc_rms_res[ k ];
		for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
			add_to_atom_id_map_after_checks( calc_rms_atom_id_map,
											pose.residue_type( n ).atom_name( q ),
											n, res_list[ n ],
											pose, native_pose_local );
		}
	}
	
	utility::vector1< Size > calc_rms_suites;
	// additional RNA suites over which to calculate RMSD
	for ( Size n = 1; n < pose.total_residue(); n++ ){
		
		if ( !pose.residue_type( n ).is_RNA() || !pose.residue_type( n + 1 ).is_RNA() ) continue;
		if ( calc_rms_res.has_value( n+1 ) ) continue;
		
		// Atoms at ends of rebuilt loops:
		if ( calc_rms_res.has_value( n ) &&
			( !pose.fold_tree().is_cutpoint( n ) || pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) ) ) {
			calc_rms_suites.push_back( n ); continue;
		}
		
		// Domain boundaries:
		if ( (res_list[ n+1 ] == res_list[ n ] + 1) &&
			fixed_domain_map[ res_list[ n ] ] != 0 &&
			fixed_domain_map[ res_list[ n+1 ] ] != 0 &&
			fixed_domain_map[ res_list[ n ] ] != fixed_domain_map[ res_list[ n+1 ] ] ){
			calc_rms_suites.push_back( n );
		}
	}
	
	utility::vector1< std::string > const extra_suite_atoms = utility::tools::make_vector1( " P  ", " OP1", " OP2", " O5'" );
	for ( Size k = 1; k <= calc_rms_suites.size(); k++ ){
		Size const n = calc_rms_suites[ k ];
		for ( Size q = 1; q <= extra_suite_atoms.size(); q++ ){
			add_to_atom_id_map_after_checks( calc_rms_atom_id_map, extra_suite_atoms[ q ],
											n+1, res_list[ n+1 ],
											pose, native_pose_local );
		}
	}
	
	//rms_map_.insert( calc_rms_atom_id_map.begin(), calc_rms_atom_id_map.end() );
	
	//		for ( std::map < AtomID, AtomID >::const_iterator it = calc_rms_atom_id_map.begin();
	//					it != calc_rms_atom_id_map.end(); it++ ){
	//			TR << it->first << " mapped to " << it->second << std::endl;
	//		}
	
	// define superposition atoms. Should be over atoms in any fixed domains. This should be
	// the 'inverse' of calc_rms atoms.
	std::map< AtomID, AtomID > superimpose_atom_id_map;
	for ( Size n = 1; n < pose.total_residue(); n++ ){
		for ( Size q = 1; q <= pose.residue_type( n ).nheavyatoms(); q++ ){
			if ( calc_rms_atom_id_map.find( AtomID( q, n ) ) == calc_rms_atom_id_map.end() ){
				add_to_atom_id_map_after_checks( superimpose_atom_id_map,
												pose.residue_type( n ).atom_name( q ),
												n, res_list[ n ],
												pose, native_pose_local );
			}
		}
	}
	
	// What if there weren't any fixed atoms? superimpose over everything.
	if ( superimpose_atom_id_map.size() == 0 ) superimpose_atom_id_map = calc_rms_atom_id_map;
	
	rmsd = 0.0;
	natoms_rmsd = calc_rms_atom_id_map.size();
	if ( natoms_rmsd > 0 && superimpose_atom_id_map.size() > 0 ) {
		//		Real const rmsd0 = rms_at_corresponding_atoms( pose, native_pose, atom_id_map );
		scoring::superimpose_pose( pose, native_pose, superimpose_atom_id_map );
		rmsd = rms_at_corresponding_atoms_no_super( pose, native_pose, calc_rms_atom_id_map );
	}
	TR << "Pose " << make_tag_with_dashes(res_list) << ": RMSD " << rmsd << " over " << natoms_rmsd << " atoms, superimposing on " << superimpose_atom_id_map.size() << " atoms. " << std::endl;
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
SWM_RMSD_Energy::superimpose_recursively( pose::Pose & pose, pose::Pose const & native_pose, Real & rmsd, Size & natoms ) const {
	
	using namespace core::pose;
	using namespace core::pose::full_model_info;
	
	Real rmsd_pose;
	Size natoms_pose;
	superimpose_at_fixed_res( pose, native_pose, rmsd_pose, natoms_pose );
	
	Real const total_sd = ( rmsd * rmsd * natoms) + (rmsd_pose * rmsd_pose * natoms_pose );
	natoms += natoms_pose;
	if ( natoms > 0 ) {
		rmsd = std::sqrt( total_sd / Real( natoms ) );
	} else {
		runtime_assert( std::abs( rmsd ) < 1e-5 );
	}
	
	utility::vector1< PoseOP > const & other_pose_list = nonconst_full_model_info( pose ).other_pose_list();
	for ( Size n = 1; n <= other_pose_list.size(); n++ ){
		superimpose_recursively( *( other_pose_list[ n ] ), native_pose, rmsd, natoms );
	}
	
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Real
SWM_RMSD_Energy::superimpose_at_fixed_res_and_get_all_atom_rmsd( pose::Pose & pose, pose::Pose const & native_pose ) const {
	Real rmsd( 0.0 );
	Size natoms( 0 );
	superimpose_recursively( pose, native_pose, rmsd, natoms );
	return rmsd;
}
	
/// c-tor
SWM_RMSD_Energy::SWM_RMSD_Energy() :
	parent( new SWM_RMSD_EnergyCreator )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace basic::options::OptionKeys::rna;
	using namespace core::chemical;
	using namespace core::io;
	
	if ( option[ in::file::native ].user() ) {
		core::io::pdb::build_pose_from_pdb_as_is(native_pose_, option[ in::file::native ]() );
	} else {
		utility_exit_with_message( "Error: must provide native pose when scoring with SWM_RMSD_Energy!\n" );
	}
}

/// clone
methods::EnergyMethodOP
SWM_RMSD_Energy::clone() const
{
	return new SWM_RMSD_Energy;
}


/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
	
///////////////////////////////////////////////////////////////////////////////
void
SWM_RMSD_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	core::pose::Pose temp_pose = pose;
	core::pose::Pose temp_native = native_pose_;
	Real rms = superimpose_at_fixed_res_and_get_all_atom_rmsd( temp_pose, temp_native );
	
	if ( rms < 1.0 ) {
		totals[ swm_rmsd ] = 0.0;
	} else {
		totals[ swm_rmsd ] = rms - 1.0;
	}

} // finalize_total_energy


///////////////////////////////////////////////////////////////////////////////
void
SWM_RMSD_Energy::eval_atom_derivative(
	id::AtomID const & id,
	pose::Pose const & pose,
	kinematics::DomainMap const &,
	ScoreFunction const &,
	EnergyMap const &,
	Vector & F1,
	Vector & F2
 	) const
{
} // eval atom derivative

core::Size
SWM_RMSD_Energy::version() const
{
	return 1; // Initial versioning
}



} // methods
} // scoring
} // core
