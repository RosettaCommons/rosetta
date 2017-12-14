// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/magnesium/util.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/magnesium/util.hh>
#include <protocols/magnesium/MgWaterHydrogenPacker.hh>
#include <protocols/magnesium/MgOrbitalFrameFinder.hh>
#include <protocols/magnesium/MgHydrater.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/variant_util.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/subpose_manipulation_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/tools/make_vector1.hh>
#include <basic/Tracer.hh>
#include <numeric/UniformRotationSampler.hh>
#include <numeric/xyzMatrix.io.hh>
#include <ObjexxFCL/format.hh>

static basic::Tracer TR( "protocols.magnesium.util" );

using namespace core;
using namespace ObjexxFCL::format;
using core::Real;
using core::Size;
using utility::vector1;
using Matrix = numeric::xyzMatrix<Real>;

namespace protocols {
namespace magnesium {

///////////////////////////////////////////////////////////////////////////////
void
fixup_magnesiums( pose::Pose & pose ) {
	// set frames -- soon will update to 'modern' version.
	MgOrbitalFrameFinder mg_orbital_frame_finder;
	mg_orbital_frame_finder.apply( pose );

	MgWaterHydrogenPacker mg_water_hydrogen_packer;
	mg_water_hydrogen_packer.apply( pose );
}

///////////////////////////////////////////////////////////////////////////////
void
hydrate_magnesiums( pose::Pose & pose,
	bool use_virtual_waters_as_placeholders /* = true */,
	bool test_all_mg_hydration_frames /* = true */ ) {
	utility::vector1< Size > const mg_res = get_mg_res( pose );
	if ( mg_res.size() == 0 ) return;
	if ( !use_virtual_waters_as_placeholders ) remove_mg_bound_waters( pose, mg_res, false /*leave_other_waters*/ );
	MgHydrater mg_hydrater( mg_res );
	mg_hydrater.set_use_fast_frame_heuristic( !test_all_mg_hydration_frames );
	mg_hydrater.set_use_virtual_waters_as_placeholders( use_virtual_waters_as_placeholders );
	mg_hydrater.set_excise_mini_pose( true );
	mg_hydrater.apply( pose );
}

///////////////////////////////////////////////////////////////////////////////
core::scoring::ScoreFunctionOP
get_mg_scorefxn() {
	using namespace scoring;
	core::scoring::ScoreFunctionOP scorefxn( new ScoreFunction );
	methods::EnergyMethodOptions options;
	options.hbond_options().use_hb_env_dep( false );
	scorefxn->set_energy_method_options( options );
	scorefxn->set_weight( fa_atr, 0.21 );
	scorefxn->set_weight( fa_rep, 0.05 ); // used for sampling, so downweight clashes a bit.
	scorefxn->set_weight( geom_sol_fast, 0.17 );
	scorefxn->set_weight( hbond_sc, 1.0 ); // this is slowest due to calculation over entire pose.

	// scorefxn->set_weight( mg, 1.0 ); // new score function.
	scorefxn->set_weight( mg_lig, 1.0 ); // new score function.
	scorefxn->set_weight( mg_sol, 1.0 ); // new score function.
	scorefxn->set_weight( mg_ref, 1.0 ); // new score function.
	scorefxn->set_weight( hoh_ref, 1.0 ); // new score function.

	return scorefxn;
}


/////////////////////////////////////////////////////////////////////////////////////////////
numeric::UniformRotationSamplerCOP
get_water_uniform_rotation_sampler() {

	using namespace numeric;

	core::Real rotstep = 15.0; //option[score::wat_rot_sampling]();
	static numeric::UniformRotationSamplerOP urs( new UniformRotationSampler( rotstep ) );

	static bool init( false );
	if ( !init ) {
		std::cout << "About to prepare rotations " << std::endl;
		// system is now invariant in rotation of 180 about X.  Eliminate 1/2 of rotations
		numeric::xyzMatrix<core::Real> flipZ = numeric::xyzMatrix<core::Real>::rows(  1,0,0,
			0,-1,0,
			0,0,-1);
		clock_t start_time = clock();
		urs->remove_redundant( flipZ );
		clock_t end_time = clock();
		init = true;
		std::cout << "remove_redundant finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

		// start_time = clock();
		// numeric::UniformRotationSampler urs_fast( rotstep );
		// urs_fast.remove_redundant_quaternion( flipZ );
		// end_time = clock();
		// std::cout << "remove_redundant_fast finished in " << double(end_time - start_time) / CLOCKS_PER_SEC << " seconds." << std::endl;

		std::cout << "Done preparing rotations " << std::endl;
	}

	return urs;
}

///////////////////////////////////////////////////////////////////////////////
numeric::UniformRotationSamplerCOP
get_octahedral_uniform_rotation_sampler( core::Real const rotstep /* = 2.5 */,
	bool const remove_redundant /*= false */ ) {
	using namespace numeric;
	using conversions::radians;
	using core::Real;
	using core::Size;

	static UniformRotationSamplerOP urs( new UniformRotationSampler( rotstep ) );

	static bool init( false );
	if ( !init ) {
		std::cout << "About to prepare rotations " << std::endl;
		if ( remove_redundant ) {

			// I'm pretty sure this closes the octahedral sub-group -- keeps a cube invariant.
			utility::vector1< xyzMatrix< core::Real > > rotateSquare( 4 ), rotateFaceToFace( 6 );

			// sub-group of 4-fold rotations
			rotateSquare[ 1 ] = xyzMatrix< core::Real >::identity();
			for ( Size i = 2; i <= 4; i++ ) {
				rotateSquare[ i ] = rotation_matrix( Vector( 1.0, 0.0, 0.0 ), radians( 90.0 * ( i - 1 ) ) );
			}

			// co-set -- put x-axis onto: -x, +y, -y, +z, -z
			rotateFaceToFace[ 1 ] = xyzMatrix< core::Real >::identity();
			rotateFaceToFace[ 2 ] = rotation_matrix( Vector( 0.0, 1.0, 0.0 ), radians( 180.0 ) );
			rotateFaceToFace[ 3 ] = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), radians( +90.0 ) );
			rotateFaceToFace[ 4 ] = rotation_matrix( Vector( 0.0, 0.0, 1.0 ), radians( -90.0 ) );
			rotateFaceToFace[ 5 ] = rotation_matrix( Vector( 0.0, 1.0, 0.0 ), radians( -90.0 ) );
			rotateFaceToFace[ 6 ] = rotation_matrix( Vector( 0.0, 1.0, 0.0 ), radians( +90.0 ) );

			for ( Size i = 1; i <= 4; i++ ) {
				for ( Size j = 1; j <= 6; j++ ) {
					if ( i == 1 && j == 1 ) continue; // identity.
					xyzMatrix< core::Real > R = rotateFaceToFace[ j ] * rotateSquare[ i ];
					urs->remove_redundant( R ); // argh, takes too long -- N^2 comparison...
				}
			}
		}
		std::cout << "Done preparing rotations " << std::endl;
		init = true;
	}

	return urs;
}

core::conformation::ResidueOP
get_useful_HOH_coords( Vector & Oc, Vector & OH1c, Vector & OH2c,
	core::chemical::ResidueTypeSet const & residue_set ) {

	core::conformation::ResidueOP rsd_canonic = conformation::ResidueFactory::create_residue( residue_set.name_map("HOH") );

	Oc = rsd_canonic->xyz("O");
	runtime_assert( Oc == 0.0 );
	OH1c = rsd_canonic->xyz("H1");
	OH2c = rsd_canonic->xyz("H2");

	// rotate so H's are symm about [1,0,0]
	core::Real cos_rot_angle = std::sqrt( 0.5* ( 1+ (OH1c.dot(OH2c)/(OH1c.length()*OH2c.length())) ) );
	core::Real sin_rot_angle = std::sqrt( 1-cos_rot_angle*cos_rot_angle );
	numeric::xyzMatrix<core::Real> rotH = numeric::xyzMatrix<core::Real>::rows(
		cos_rot_angle,sin_rot_angle, 0,
		-sin_rot_angle,cos_rot_angle, 0,
		0,            0, 1);
	OH1c = rotH*OH1c;
	OH2c = rotH*OH2c;
	rsd_canonic->set_xyz( "H1", OH1c );
	rsd_canonic->set_xyz( "H2", OH2c );
	return rsd_canonic;
}


///////////////////////////////////////////////////////////////////////////////
// to determine orbital frame, should only use ligands that really are acceptors.
bool
is_ligand( core::pose::Pose const & pose, core::id::AtomID const & ligand )
{
	return pose.residue( ligand.rsd() ).atom_type( ligand.atomno() ).is_acceptor();
}

// handles Mg res + ligand AtomID
utility::vector1< std::pair< Size, core::id::AtomID > >
filter_acceptor_ligands( pose::Pose const & pose, utility::vector1< std::pair< Size, core::id::AtomID > > const & ligands ) {
	using namespace core::id;
	utility::vector1< std::pair< Size, core::id::AtomID > > acceptor_ligands;
	for ( Size n = 1; n <= ligands.size(); n++ ) {
		AtomID const & ligand = ligands[ n ].second;
		if ( is_ligand( pose, ligand ) ) acceptor_ligands.push_back( ligands[ n ] );
	}
	return acceptor_ligands;
}

// handles ligand AtomID
utility::vector1< core::id::AtomID >
filter_acceptor_ligands( core::pose::Pose const & pose, utility::vector1< core::id::AtomID > const & ligands ) {
	using namespace core::id;
	utility::vector1< core::id::AtomID  > acceptor_ligands;
	for ( Size n = 1; n <= ligands.size(); n++ ) {
		AtomID const & ligand = ligands[ n ];
		if ( pose.residue( ligand.rsd() ).atom_type( ligand.atomno() ).is_acceptor() ) acceptor_ligands.push_back( ligands[ n ] );
	}
	return acceptor_ligands;
}

///////////////////////////////////////////////////////////////////////////////
// to determine orbital frame, should only use ligands that really are acceptors.
utility::vector1< std::pair< Size, core::id::AtomID > >
filter_water_ligands( pose::Pose const & pose, utility::vector1< std::pair< Size, core::id::AtomID > > const & ligands ) {
	using namespace core::id;
	utility::vector1< std::pair< Size, core::id::AtomID > > water_ligands;
	for ( Size n = 1; n <= ligands.size(); n++ ) {
		AtomID const & ligand = ligands[ n ].second;
		if ( pose.residue( ligand.rsd() ).name3() == "HOH" &&
				pose.residue( ligand.rsd() ).atom_name( ligand.atomno() ) == " O  " ) {
			water_ligands.push_back( ligands[ n ] );
		}
	}
	return water_ligands;
}

///////////////////////////////////////////////////////////////////////////////
// The idea:
//   Go through each mg(2+)-water that could be a pair.
//   But then some waters could have two possible parent Mg(2+) ions.
//   Decide which Mg(2+) gets the ion based on distance, and don't allow
//    more than 6 waters per Mg(2+).
//
//  This is heuristic and may fail -- some waters may not get assigned, in which case error out.
//
std::map< Size, vector1< Size > >
define_mg_water_map( pose::Pose const & pose )
{
	std::map< Size, vector1< Size > > mg_water_map; // each mg mapped to a list of up to 6 waters
	std::map< Size, Size > water_to_magnesium; // each water assigned to one magnesium
	utility::vector1< std::pair< Size, Size > > mg_water_pairs = get_mg_water_pairs( pose, false /*exclude_virtual_water*/ );

	TR << "NUM WATER PAIRS " << mg_water_pairs.size() << std::endl;

	for ( Size n = 1; n <= mg_water_pairs.size(); n++ ) {
		Size const mg_res    = mg_water_pairs[ n ].first;
		Size const water_res = mg_water_pairs[ n ].second;
		if ( water_to_magnesium[ water_res ] > 0 ) continue; // already assigned.
		if ( mg_water_map[ mg_res ].size() == 6 ) continue; // mg is full.
		mg_water_map[ mg_res ].push_back( water_res );
		water_to_magnesium[ water_res ] = mg_res;
	}

	for ( std::map< Size, Size >::const_iterator it = water_to_magnesium.begin(); it != water_to_magnesium.end(); it++ ) {
		if ( it->second == 0 ) {
			utility_exit_with_message( "Problem -- did not assign water " + I( 3, it->first )
				+ " to a magnesium because the magnesiums have full shells" );
		}
	}
	return mg_water_map;
}

///////////////////////////////////////////////////////////////////////////////
utility::vector1< std::pair< Size, core::id::AtomID > >
get_mg_ligands( pose::Pose const & pose, utility::vector1< Size > const mg_res,
	bool const filter_for_acceptors /* true */,
	bool const exclude_virtual_waters /* true */,
	bool const just_waters = false )
{
	using namespace core::conformation;
	using namespace core::id;

	// keep track of distances so that later we can sort by closeness.
	utility::vector1< std::pair< Distance, std::pair< Size, core::id::AtomID > > > distance_ligand_pairs;
	for ( Size n = 1; n <= mg_res.size(); n++ ) {
		Size const i = mg_res[ n ];
		Vector const xyz_mg =  pose.xyz( NamedAtomID( "MG  ", i ) );

		// would be much faster to only look over neighbors
		for ( Size j = 1; j <= pose.size(); j++ ) {
			if ( i == j ) continue;

			Residue const & rsd_j = pose.residue( j );
			// look through all non-hydrogen atoms
			for ( Size jj = 1; jj <= rsd_j.nheavyatoms(); jj++ ) {

				if ( rsd_j.is_virtual( jj ) &&
						( exclude_virtual_waters || rsd_j.aa() != core::chemical::aa_h2o ) ) continue;

				Vector const & xyz_jj = rsd_j.xyz( jj );
				core::Real distance = (xyz_jj - xyz_mg).length();
				if ( distance < MG_LIGAND_DISTANCE_CUTOFF ) {
					distance_ligand_pairs.push_back( std::make_pair( distance, std::make_pair( i, core::id::AtomID( jj, j ) ) ) );
				}
			}

		} // loop over RNA residues to find distance_ligand_pairs.
	} // loop over Mg(2+) residues

	std::sort( distance_ligand_pairs.begin(), distance_ligand_pairs.end() );

	utility::vector1< std::pair< Size, core::id::AtomID > > ligands;
	for ( Size n = 1; n <= distance_ligand_pairs.size(); n++ ) {
		ligands.push_back( distance_ligand_pairs[ n ].second );
	}
	if ( filter_for_acceptors ) {
		ligands = filter_acceptor_ligands( pose, ligands );
	}
	if ( just_waters ) {
		ligands = filter_water_ligands( pose, ligands );
	}
	return ligands;
}

///////////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::AtomID >
get_mg_ligands( pose::Pose const & pose, Size const i,
	bool const filter_for_acceptors /* true */,
	bool const exclude_virtual_waters /* true */ )
{
	vector1< std::pair< Size, core::id::AtomID > > mg_and_ligands = get_mg_ligands( pose, utility::tools::make_vector1( i ),
		filter_for_acceptors, exclude_virtual_waters );
	utility::vector1< core::id::AtomID > ligands;
	for ( Size n = 1; n <= mg_and_ligands.size(); n++ ) ligands.push_back( mg_and_ligands[ n ].second );
	return ligands;
}

///////////////////////////////////////////
utility::vector1< std::pair< Size, Size > >
get_mg_water_pairs( pose::Pose const & pose,
	bool const exclude_virtual_waters /* = true */ ) {
	return get_mg_water_pairs( pose, get_mg_res( pose ), exclude_virtual_waters );
}

///////////////////////////////////////////
utility::vector1< std::pair< Size, Size > >
get_mg_water_pairs( pose::Pose const & pose,
	vector1< Size > const & mg_res,
	bool const exclude_virtual_waters /* = true */ ) {

	vector1< std::pair< Size, id::AtomID > > const ligands = get_mg_ligands( pose, mg_res, false /*filter_for_acceptors*/,
		exclude_virtual_waters, true /*just_waters*/ );
	utility::vector1< std::pair< Size, Size > > mg_water_pairs;

	for ( Size n = 1; n <= ligands.size(); n++ ) {
		mg_water_pairs.push_back( std::make_pair( ligands[n].first, ligands[n].second.rsd() ) );
	}
	return mg_water_pairs;
}

///////////////////////////////////////////////////////////////////////////////
// If fold tree is set up correctly, expect to find 6 waters as daughters of
// each Mg(2+). Some might be virtual, we don't care.
utility::vector1< core::Size >
find_bound_waters_that_are_daughters_in_fold_tree( pose::Pose const & pose, Size const mg_res )
{
	vector1< Size > water_res;
	for ( Size i = 1; i <= pose.size(); i++ ) {
		if ( pose.residue_type( i ).aa() != core::chemical::aa_h2o ) continue;
		if ( pose.fold_tree().get_parent_residue( i ) != mg_res ) continue;
		if ( ( pose.residue( mg_res ).xyz( "MG  " ) - pose.residue( i ).xyz(" O  " ) ).length() > MG_LIGAND_DISTANCE_CUTOFF ) continue;
		water_res.push_back( i );
	}
	return water_res;
}

///////////////////////////////////////////////////////////////////////////////
core::Size
get_bound_water_res( pose::Pose const & pose, Size const mg_res, Size const n /* = 1..6 */) {
	runtime_assert( n >= 1);
	runtime_assert( n <= 6);
	vector1< Size > water_res = find_bound_waters_that_are_daughters_in_fold_tree( pose, mg_res );
	if ( water_res.size() != 6 ) {
		TR << TR.Red << pose.annotated_sequence() << std::endl;
		TR << TR.Red << pose.fold_tree() << std::endl;
		TR << TR.Red << "Problem: " << mg_res << " has " << water_res.size() << " waters, not 6." << std::endl;
	}
	runtime_assert( water_res.size() == 6 );
	return water_res[ n ];
}

//////////////////////////////////////////////////////////
core::Size
instantiate_water_at_octahedral_vertex( pose::Pose & pose,
	Size const mg_res,
	Size const n /* 1 ... 6*/,
	Distance const hoh_distance /*= 2.1 */,
	bool const replace_residue /* = false */,
	bool const virtual_water /* = false */ ) {
	using namespace core::conformation;
	runtime_assert( n >= 1 && n <= 6 );
	Vector const & xyz_mg( pose.residue( mg_res ).xyz( 1 ) );
	utility::vector1< Vector > vec( 6 ); // unit vectors pointing to V1, ... V6 'orbital' virtual atoms.
	for ( Size k = 1; k <= 6; k++ ) {
		vec[ k ] = ( pose.residue( mg_res ).xyz( k+1 ) - xyz_mg ).normalized();
	}
	Vector Oc, OH1c, OH2c;
	ResidueOP water_rsd = get_useful_HOH_coords( Oc, OH1c, OH2c, *pose.residue_type_set_for_pose( pose.residue_type( mg_res ).mode() ) );
	// set appropriate frame, pointing towards water, and with axes set by
	// O_h geometry.
	Vector x = vec[ n ];
	Vector y = ( n == 6 ) ? vec[ 1 ] : vec[ n+1 ];
	Vector z = cross( x, y );

	Vector const xyz_water = xyz_mg + hoh_distance * vec[ n ];
	core::kinematics::Stub stub( Matrix::cols( x, y, z ), xyz_water );
	water_rsd->set_xyz( " O  ", stub.local2global( Oc )  );
	water_rsd->set_xyz( " H1 ", stub.local2global( OH1c ) );
	water_rsd->set_xyz( " H2 ", stub.local2global( OH2c ) );

	Size water_res( 0 );
	if ( replace_residue ) {
		water_res = get_bound_water_res( pose, mg_res, n );
		pose.replace_residue( water_res, *water_rsd, false /*orient_backbone*/ );
	} else {
		water_res = append_mg_bound_water( pose, *water_rsd, mg_res );
		update_numbers_in_pdb_info( pose );
	}
	runtime_assert( water_res > 0 );

	if ( virtual_water ) {
		add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RESIDUE_VARIANT, water_res );
	} else {
		remove_variant_type_from_pose_residue( pose, core::chemical::VIRTUAL_RESIDUE_VARIANT, water_res );
	}
	return water_res;
}


///////////////////////////////////////////////////////////////////////////////
void
update_jump_atoms_for_mg_bound_water( core::pose::Pose & pose, core::Size const n )
{
	core::kinematics::FoldTree f = pose.fold_tree();
	Size const njump( f.get_jump_that_builds_residue( n ) );
	runtime_assert( pose.residue_type( f.upstream_jump_residue(   njump ) ).name3() == " MG" );
	runtime_assert( pose.residue_type( f.downstream_jump_residue( njump ) ).name3() == "HOH" );
	f.set_jump_atoms( njump, "MG  ", " O  ", true );
	pose.fold_tree( f );
}

///////////////////////////////////////////////////////////////////////////////
Size
append_mg_bound_water(  core::pose::Pose & pose,
	core::conformation::Residue const & rsd,
	core::Size const mg_res )
{
	pose.append_residue_by_jump( rsd, mg_res );
	update_jump_atoms_for_mg_bound_water( pose, pose.size() );
	return ( pose.size() );
}


///////////////////////////////////////////////////////////////////////////////
core::conformation::ResidueOP
get_mg_rsd() {
	using namespace core::chemical;
	// stick following in a get_single_mg_pose function.
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	return core::conformation::ResidueFactory::create_residue ( *( rsd_set->get_representative_type_name3( " MG" ) ) );
}

///////////////////////////////////////////////////////////////////////////////
void
add_single_magnesium( pose::Pose & pose )
{
	pose.append_residue_by_jump ( *get_mg_rsd(), pose.size() );
}

////////////////////////////////////////////////////
void
strip_out_magnesiums( pose::Pose & pose ){
	utility::vector1< Size > slice_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue(n).name3() == " MG" ) continue;
		slice_res.push_back( n );
	}

	core::pose::pdbslice( pose, slice_res );
}

///////////////////////////////////////////
utility::vector1< Size >
get_res_with_name( pose::Pose const & pose, std::string const & name ) {
	utility::vector1< Size > mg_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue( n ).name3() == name ) mg_res.push_back( n );
	}
	return mg_res;
}

///////////////////////////////////////////
utility::vector1< Size >
get_mg_res( pose::Pose const & pose ) {
	return get_res_with_name( pose, " MG" );
}

///////////////////////////////////////////
utility::vector1< Size >
get_water_res( pose::Pose const & pose ) {
	return get_res_with_name( pose, "HOH" );
}

///////////////////////////////////////////
void
remove_waters_except_mg_bound( pose::Pose & pose,
	utility::vector1< std::pair< Size, Size > > const & mg_water_pairs ) {
	utility::vector1< Size > slice_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue( n ).name3() == "HOH" ) continue;
		slice_res.push_back( n );
	}
	for ( Size i = 1; i <= mg_water_pairs.size(); i++ ) {
		Size const n = mg_water_pairs[ i ].second;
		if ( slice_res.has_value( n ) ) continue;
		slice_res.push_back( n );
	}
	std::sort( slice_res.begin(), slice_res.end() );
	pdbslice( pose, slice_res);
}

///////////////////////////////////////////
void
remove_mg_bound_waters( pose::Pose & pose, utility::vector1< Size > const & mg_res, bool const leave_other_waters /* = false */ ) {

	using namespace core::pose;
	using namespace core::pose::full_model_info;

	vector1< std::pair< Size, Size > > mg_water_pairs = get_mg_water_pairs( pose, mg_res );
	vector1< Size > mg_bound_waters;
	for ( Size i = 1; i <= mg_water_pairs.size(); i++ ) {
		Size const n = mg_water_pairs[ i ].second;
		mg_bound_waters.push_back( n );
	}

	utility::vector1< Size > slice_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue( n ).name3() == "HOH" ) {
			if ( !leave_other_waters ) continue; // no waters in final pose!
			if ( mg_bound_waters.has_value( n ) ) continue;
		}
		slice_res.push_back( n );
	}
	std::sort( slice_res.begin(), slice_res.end() );

	if ( full_model_info_defined( pose ) ) {
		// use slice_out_pose which is smart about fold-tree and full_model_info.

		// water res to remove is everything that was not in slice_res, defined above.
		vector1< Size > const & res_list = const_full_model_info( pose ).res_list();
		vector1< Size > remove_water_res_in_full_model_numbering;
		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( !slice_res.has_value( n ) ) remove_water_res_in_full_model_numbering.push_back( res_list[ n ] );
		}

		// have to use slice_out_pose one at a time -- makes assumption that there is a single connection (jump or bond)
		// between the sliced out segment and the remainder pose
		Pose sliced_out_waters_pose;
		for ( Size n = 1; n <= remove_water_res_in_full_model_numbering.size(); n++ ) {
			Size const water_res = const_full_model_info( pose ).res_list().index( remove_water_res_in_full_model_numbering[ n ] );
			stepwise::modeler::slice_out_pose( pose, sliced_out_waters_pose, utility::tools::make_vector1( water_res ) );
		}
	} else {
		pdbslice( pose, slice_res);
	}
}

//////////////////////////////////////////
void
set_water_numbers_to_zero( pose::Pose & pose ) {
	using namespace core::pose;
	PDBInfoCOP pdb_info = pose.pdb_info();
	vector1< Size > numbering;
	vector1< char > chains;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue( n ).name3() == "HOH" ) {
			numbering.push_back( 0 );
		} else {
			numbering.push_back( pdb_info->number( n ) );
		}
		chains.push_back( pdb_info->chain( n ) );
	}

	PDBInfoOP new_pdb_info( new PDBInfo( pose ) );
	new_pdb_info->set_numbering( numbering );
	new_pdb_info->set_chains( chains );
	pose.pdb_info( new_pdb_info );
}

//////////////////////////////////////////
void
update_numbers_in_pdb_info( pose::Pose & pose, bool const reset_waters /* = false */ ) {
	using namespace core::pose;
	if ( reset_waters ) set_water_numbers_to_zero( pose );
	PDBInfoCOP pdb_info = pose.pdb_info();
	int max_number( 0 );
	vector1< Size > numbering;
	vector1< char > chains;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pdb_info->number( n ) > max_number ) max_number = pdb_info->number( n );
		if ( pdb_info->number( n ) == 0 ) {
			max_number++;
			numbering.push_back( max_number );
			chains.push_back( ' ' ); // otherwise get a warning about '^' placeholder
		} else {
			numbering.push_back( pdb_info->number( n ) );
			chains.push_back( pdb_info->chain( n ) );
		}
	}

	PDBInfoOP new_pdb_info( new PDBInfo( pose ) );
	new_pdb_info->set_numbering( numbering );
	new_pdb_info->set_chains( chains );
	pose.pdb_info( new_pdb_info );

}

utility::vector1< core::Size >
pdbslice( core::pose::Pose & pose, core::Size const center_res, core::Distance distance_cutoff /* = 12.0 */ ) {
	vector1< Size > slice_res;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( ( pose.residue( n ).nbr_atom_xyz() - pose.residue( center_res ).nbr_atom_xyz() ).length() < distance_cutoff ) {
			slice_res.push_back( n );
		}
	}
	runtime_assert( slice_res.has_value( center_res ) );
	pdbslice( pose, slice_res);
	return slice_res;
}


utility::vector1< numeric::xyzVector< core::Real > >
get_hoh_xyz( pose::Pose const & pose, Size const pdb_mg_res, Size & nlig ) {
	using namespace core::conformation;
	using namespace core::id;
	using namespace protocols::magnesium;
	utility::vector1< numeric::xyzVector< core::Real > > hoh_xyz;
	Size pose_mg_res = pdb_to_pose( pose, pdb_mg_res );
	vector1< AtomID > mg_ligands = get_mg_ligands( pose, pose_mg_res );
	nlig = mg_ligands.size();
	for ( Size i = 1; i <= mg_ligands.size(); i++ ) {
		AtomID const & ligand = mg_ligands[ i ];
		Residue const & rsd = pose.residue( ligand.rsd() );
		if ( rsd.name3() != "HOH" ) continue;
		hoh_xyz.push_back( rsd.xyz( ligand.atomno() ) );
	}
	return hoh_xyz;
}

//
// what stats would be useful?
//
// (1) coordination number - ref
// (2) number of direct contacts - ref
// (3) number of hohs - ref
// (4) number of hohs that make contacts - ref
// (5) score of mg(2+) (subtract with and without)
//
// and *again* for model
//
// and then RMSD of waters
//
void
get_hydration_stats( pose::Pose const & pose,
	pose::Pose const & reference_pose,
	utility::vector1< Size > const & pdb_mg_res_list_in,
	std::string const & outfile ){

	using namespace protocols::magnesium;

	utility::io::ozstream out( outfile, std::ios_base::app );
	vector1< Size> pdb_mg_res_list( pdb_mg_res_list_in );
	if ( pdb_mg_res_list.size() == 0 ) pdb_mg_res_list = pose_to_pdb( pose, get_mg_res( pose ) );
	for ( Size n = 1; n <= pdb_mg_res_list.size(); n++ ) {
		Size pdb_mg_res = pdb_mg_res_list[ n ];
		//  runtime_assert( );
		Size nlig_ref( 0 ), nlig_model( 0 );
		vector1< Vector > ref_hoh_xyz   = get_hoh_xyz( reference_pose, pdb_mg_res, nlig_ref );
		vector1< Vector > model_hoh_xyz = get_hoh_xyz( pose, pdb_mg_res, nlig_model );
		Size const num_ref_hoh(   ref_hoh_xyz.size() );
		Size const num_model_hoh( model_hoh_xyz.size() );
		core::Real rmsd( 0.0 );
		for ( Size i = 1; i <= ref_hoh_xyz.size(); i++ ) {
			Distance best_distance( 10.0 ); // placeholder
			for ( Size j = 1; j <= model_hoh_xyz.size(); j++ ) {
				Distance current_distance( ( model_hoh_xyz[ j ] - ref_hoh_xyz[ i ] ).length() );
				if ( current_distance < best_distance ) best_distance = current_distance;
			}
			rmsd += best_distance*best_distance;
		}
		if ( num_ref_hoh > 0 ) rmsd = std::sqrt( rmsd / core::Real( num_ref_hoh ) );
		core::Real const Bfactor = reference_pose.pdb_info()->temperature( pdb_to_pose( reference_pose, pdb_mg_res ), 1 );
		//  TR << I( 5, pdb_mg_res ) << " " << I( 2, num_ref_hoh ) << " " << I( 2, num_model_hoh ) << " " << Bfactor << F( 8, 3, rmsd ) << std::endl;
		out << I( 5, pdb_mg_res )
			<< " " << I( 2, nlig_ref )   << " " <<  I( 2, num_ref_hoh )
			<< " " << I( 2, nlig_model ) << " " << I( 2, num_model_hoh )
			<< " " << Bfactor << F( 8, 3, rmsd ) << std::endl;
	}
	out.close();
}


} //magnesium
} //protocols
