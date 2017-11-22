// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/rna/movers/ErraserMinimizerMover.cc
/// @brief
/// @author Andy Watkins, amw579@stanford.edu

#include <protocols/rna/movers/ErraserMinimizerMover.hh>
#include <protocols/rna/movers/ErraserMinimizerMoverCreator.hh>

#include <core/pose/rna/RNA_IdealCoord.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/util.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/sequence/util.hh>
#include <core/sequence/Sequence.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <basic/Tracer.hh>
#include <basic/database/open.hh>

#include <protocols/viewer/viewers.hh>

#include <core/io/pdb/build_pose_as_is.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/rna/RNA_TorsionPotential.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/rna/util.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/annotated_sequence.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pack/rotamer_trials.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/moves/Mover.hh>

#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/constraints/util.hh>

//////////////////////////////////////////////////
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
///////////////////////////////////////////////////

#include <core/pose/PDBInfo.hh>
#include <core/chemical/rna/RNA_FittedTorsionInfo.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/types.hh>

#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/options/OptionCollection.hh>
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

#include <fstream>
#include <algorithm>

#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.rna.movers.ErraserMinimizerMover" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::rna::erraser;
using namespace core;
using namespace core::chemical;
using namespace core::chemical::rna;
using namespace core::conformation;
using namespace core::id;
using namespace core::kinematics;
using namespace core::optimization;
using namespace core::pose;
using namespace core::pose::rna;
using namespace core::scoring;
using namespace core::scoring::constraints;
using namespace core::scoring::func;
using namespace numeric::conversions;
using namespace protocols::stepwise::modeler::rna;
using namespace protocols::moves;
using namespace protocols::filters;
using utility::vector1;

namespace protocols {
namespace rna {
namespace movers {

// AMW: TODO
// 1. Setting up NCNT ideal poses for ideal coordinate constraints is expensive
// to do at construction. What if we only did it for the NCNTs actually found in
// the pose (usually 1-2, not hundreds)? Maybe we could cache constraints and/or
// poses for repeated inputs... Poses for sure.

utility::vector1< core::Size >
string_to_size_vector( std::string const & sv ) {
	utility::vector1< core::Size > vec;

	std::stringstream ss( sv );
	Size s;
	while ( ss >> s ) {
		vec.push_back( s );
	}

	return vec;
}

typedef utility::pointer::shared_ptr< ErraserMinimizerMover > ErraserMinimizerMoverOP;
typedef utility::pointer::shared_ptr< ErraserMinimizerMover const > ErraserMinimizerMoverCOP;

ErraserMinimizerMover::ErraserMinimizerMover():
	Mover("ErraserMinimizerMover"),
	vary_bond_geometry_( true ),
	constrain_phosphate_( true ),
	ready_set_only_( false ),
	skip_minimize_( false ),
	scorefxn_( new ScoreFunction() ),
	edens_scorefxn_( new ScoreFunction() )
{
	initialize_from_options( option );
}

void ErraserMinimizerMover::initialize_from_options( utility::options::OptionCollection const & options ) {
	// We're going to keep this default true -- if you have to turn it on with this option, then
	// that turns it on for stepwise steps too. That's terrible.
	//vary_bond_geometry_ = options[ basic::options::OptionKeys::rna::vary_geometry ].value();
	constrain_phosphate_ = options[ constrain_P ].value();
	ready_set_only_ =     options[ ready_set_only ].value();
	skip_minimize_ =      options[ skip_minimize ].value();
	attempt_pyrimidine_flip_ = options[ attempt_pyrimidine_flip ].value();
	std::copy( options[ fixed_res ].value().begin(), options[ fixed_res ].value().end(), std::inserter( fixed_res_list_, fixed_res_list_.end() ) );
	cutpoint_list_ = options[ full_model::cutpoint_open ].value();

	if ( !ready_set_only_ ) {
		//Setup score function.
		std::string score_weight_file = "stepwise/rna/rna_res_level_energy4_elec_dens";
		if ( options[ basic::options::OptionKeys::score::weights ].user() ) {
			score_weight_file = options[ basic::options::OptionKeys::score::weights ]();
			TR << "User passed in score:weight option: " << score_weight_file << std::endl;
		}
		scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
		//edens_scorefxn_->set_weight( elec_dens_atomwise, scorefxn_->get_weight( elec_dens_atomwise ) );
	}
}

void ErraserMinimizerMover::parse_my_tag( TagCOP tag,
	basic::datacache::DataMap & /*data*/,
	Filters_map const & /*filters*/,
	Movers_map const & /*movers*/,
	Pose const & /*pose*/
) {
	vary_bond_geometry_ = tag->getOption< bool >( "vary_bond_geometry", vary_bond_geometry_ );
	constrain_phosphate_ = tag->getOption< bool >( "constrain_phosphate", constrain_phosphate_ );
	ready_set_only_ = tag->getOption< bool >( "ready_set_only", ready_set_only_ );
	skip_minimize_ = tag->getOption< bool >( "skip_minimize", skip_minimize_ );
	attempt_pyrimidine_flip_ = tag->getOption< bool >( "attempt_pyrimidine_flip", attempt_pyrimidine_flip_ );

	utility::vector1< Size > temp_vector = string_to_size_vector( tag->getOption< std::string >( "fixed_res_list", "" ) );
	std::copy( temp_vector.begin(), temp_vector.end(), std::inserter( fixed_res_list_, fixed_res_list_.end() ) );

	cutpoint_list_ = string_to_size_vector( tag->getOption< std::string >( "cutpoint_list", "" ) );
	output_pdb_name_ = tag->getOption< std::string >( "output_pdb_name", output_pdb_name_ );

	if ( !ready_set_only_ ) {
		//Setup score function.
		std::string score_weight_file = "stepwise/rna/rna_res_level_energy4_elec_dens";
		if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
			score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
			TR << "User passed in score:weight option: " << score_weight_file << std::endl;
		}
		scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
		//edens_scorefxn_->set_weight( elec_dens_atomwise, scorefxn_->get_weight( elec_dens_atomwise ) );
	}
}

// Virts and sidechain atoms that aren't the first base atom should not move
bool
ErraserMinimizerMover::i_want_this_atom_to_move(
	chemical::ResidueType const & residue_type,
	Size const & atomno
) {
	// First trial. I just want to fix this nonplanar ring, but maybe this will make tiny improvements in general.
	runtime_assert( scorefxn_->get_weight( cart_bonded ) > 0 );

	if ( residue_type.is_virtual( atomno ) ) {
		//  TR << "Is this virtual? " << residue2.atom_name( atomno ) << std::endl;
		return false;
	}

	return true;
}

//////////////////////////////////////////////////////////////////////////////
// Following has not (yet) been carefully debugged.
void
ErraserMinimizerMover::vary_bond_geometry(
	core::kinematics::MoveMap & mm,
	pose::Pose & pose,
	ObjexxFCL::FArray1D< bool > & allow_insert, // Operationally: not fixed, cutpoint, virt
	std::set< Size > const & chunk
) {
	Size const nres( pose.size() );
	TR << "Enter vary_bond_geometry....." << std::endl;

	TR << "My impression of what residues should maybe move: " << std::endl;
	TR << "[ ";
	for ( auto const elem : chunk ) TR << elem << " ";
	TR << "]" << std::endl;
	// First, set appropriate DOFs to move in the movemap, mm
	for ( Size i = 1; i <= nres; ++i )  {

		// Don't do anything for protein residues, because we don't have them as ideals.
		// In the future, apply cart_bonded.
		if ( pose.residue_type( i ).is_protein() ) continue;

		if ( chunk.find( i ) == chunk.end() ) continue;
		if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue; //FCC
		// can't skip just for allow_insert false: carbohydrates can misleadingly have this
		if ( !allow_insert( i ) && !pose.residue_type( i ).is_carbohydrate() ) continue;

		chemical::ResidueType const & residue_type( pose.residue_type( i ) );

		for ( Size j = 1; j <= residue_type.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move( residue_type, j ) ) continue;

			tree::AtomCOP current_atom( pose.atom_tree().atom( AtomID( j, i ) ).get_self_ptr() );
			if ( current_atom->is_jump() ) continue;

			// STEP ONE: DISTANCES
			tree::AtomCOP j_distatom( current_atom->input_stub_atom1() );

			if ( !j_distatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_distatom->id() ) ) continue;

			mm.set( DOF_ID( AtomID( j, i ), D ), true );

			if ( j_distatom->is_jump() ) continue;

			// STEP TWO: ANGLES
			core::kinematics::tree::AtomCOP j_angleatom( current_atom->input_stub_atom2() );

			if ( !j_angleatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_angleatom->id() ) ) continue;
			if ( j_angleatom == current_atom ) continue;

			mm.set( DOF_ID( AtomID( j, i ), THETA ), true );

			if ( j_angleatom->is_jump() ) continue;

			// STEP THREE: DIHEDRALS
			core::kinematics::tree::AtomCOP j_dihedralatom( current_atom->input_stub_atom3() );

			if ( !j_dihedralatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_dihedralatom->id() ) ) continue;
			if ( j_dihedralatom == current_atom ) continue;

			mm.set( DOF_ID( AtomID( j, i ), PHI ), true );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
//FCC: Adding Virtual res
int
ErraserMinimizerMover::add_virtual_res( core::pose::Pose & pose ) {
	int nres = pose.size();

	// if already rooted on virtual residue , return
	if ( pose.residue_type( pose.fold_tree().root() ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "add_virtual_res() called but pose is already rooted on a VRT residue ... continuing." << std::endl;
		return pose.fold_tree().root();
	}

	if ( pose.residue_type( nres ).aa() == core::chemical::aa_vrt ) {
		TR.Warning << "add_virtual_res() called but pose already has a VRT residue ... rerooting." << std::endl;
		kinematics::FoldTree newF( pose.fold_tree() );
		newF.reorder( nres );
		pose.fold_tree( newF );
		return nres;
	}

	// attach virt res there
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *core::pose::virtual_type_for_pose(pose) ) );

	// OK, what we need is to save the PDBInfo, then add it back for every residue.
	PDBInfo info = *pose.pdb_info();
	pose.append_residue_by_jump( *new_res , nres );
	pose.pdb_info()->copy( info, 1, nres, 1 );

	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres + 1 );
	pose.fold_tree( newF );

	core::pose::full_model_info::FullModelInfoOP full_model_info( new core::pose::full_model_info::FullModelInfo( pose ) );
	set_full_model_info( pose, full_model_info );
	// This is fine
	pose.pdb_info()->obsolete( false );

	return nres + 1;
}

///////////////////////////////////////////////////////////
void
ErraserMinimizerMover::setup_fold_tree( pose::Pose & pose ) {
	Size const nres( pose.size() );
	Size const num_jumps( cutpoint_list_.size() );
	ObjexxFCL::FArray2D< Size > jump_points( 2, num_jumps );
	ObjexxFCL::FArray1D< Size > cuts( num_jumps );

	for ( Size n = 1; n <= cutpoint_list_.size(); n++ ) {
		jump_points( 1, n ) = cutpoint_list_[n];
		jump_points( 2, n ) = cutpoint_list_[n] + 1;
		cuts( n ) = cutpoint_list_[n];
	}

	FoldTree f( nres );
	f.tree_from_jumps_and_cuts( nres, num_jumps, jump_points, cuts, 1, false );
	pose.fold_tree( f );
}

///////////////////////////////////////////
void
ErraserMinimizerMover::pyrimidine_flip_trial( pose::Pose & pose )
{
	Size const total_res = pose.size();
	Pose screen_pose = pose;
	Real orig_score, new_score;
	orig_score = ( *scorefxn_ )( pose );
	new_score = ( *scorefxn_ )( screen_pose );
	TR << "Start pyrimidine_flip_trial. Flip residue :";
	for ( Size i = 1; i <= total_res; ++i ) {
		if ( fixed_res_list_.find( i ) != fixed_res_list_.end() ) continue;

		if ( !pose.residue_type( i ).is_pyrimidine() ) continue;

		// Pyrimidine flipping?
		Real const orig_chi = pose.torsion( TorsionID( i, id::CHI, 1 ) );
		Real const new_chi = orig_chi + 180.0;
		screen_pose.set_torsion( TorsionID( i, id::CHI, 1 ), new_chi );
		new_score = ( *scorefxn_ )( screen_pose );
		if ( new_score < orig_score ) { //Flip the chi!
			pose.set_torsion( TorsionID( i, id::CHI, 1 ), new_chi );
			orig_score = new_score;
			TR << ' ' << i;
		} else { //Keep the original chi
			screen_pose.set_torsion( TorsionID( i, id::CHI, 1 ), orig_chi );
		}
	}
	TR << std::endl;
}

///////////////////////////////////////////

void
ErraserMinimizerMover::add_fixed_res_constraints(
	pose::Pose & pose,
	Size const fixed_res_num,
	Size const my_anchor
) {
	Residue const & rsd( pose.residue( fixed_res_num ) );
	ConstraintSetOP cst_set = pose.constraint_set()->clone();

	if ( !rsd.has( "P" ) ) return;

	Real const coord_sdev( 0.1 );
	Size const atm_indexP   = rsd.atom_index( "P" );
	Size const atm_indexO3  = rsd.atom_index( "O3'" );
	Size const atm_indexOP2 = rsd.atom_index( "OP2" );
	Size const atm_indexC6  = rsd.atom_index( "C6" );
	// Note that this is safe: if we can ask for these phosphate names,
	// we can ask for chi1 atom 3 and be confident it's a nucleobases's
	// first base atom.
	Size const atm_indexBase = chemical::rna::first_base_atom_index( rsd.type() );

	cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
		AtomID( atm_indexP, fixed_res_num ),
		AtomID( 1, my_anchor ), rsd.xyz( atm_indexP ),
		FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) );
	cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
		AtomID( atm_indexO3, fixed_res_num ),
		AtomID( 1, my_anchor ), rsd.xyz( atm_indexO3 ),
		FuncOP( new HarmonicFunc ( 0.0, coord_sdev ) ) ) ) );
	cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
		AtomID( atm_indexBase, fixed_res_num ),
		AtomID( 1, my_anchor ), rsd.xyz ( atm_indexBase ),
		FuncOP( new HarmonicFunc ( 0.0, coord_sdev ) ) ) ) );
	cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
		AtomID( atm_indexC6, fixed_res_num ),
		AtomID( 1, my_anchor ), rsd.xyz( atm_indexC6 ),
		FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) );
	cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
		AtomID( atm_indexOP2, fixed_res_num ),
		AtomID( 1, my_anchor ), rsd.xyz( atm_indexOP2 ),
		FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) );

	pose.constraint_set( cst_set );
}

template< typename T >
void remove_set1_elements_from_set2(
	std::set< T > const & set1,
	std::set< T > & set2
) {
	for ( auto it = set1.begin(); it != set1.end(); ++it ) {
		set2.erase( set2.find( *it ) );
	}
}

void
fill_gaps_and_remove_isolated_res(
	std::set< Size > & res_list,
	Size const total_res,
	std::set< Size > & res_remove
) {

	// res_list is sorted
	if ( *std::next(res_list.begin()) - *res_list.begin() != 1 ) {
		res_remove.insert( *res_list.begin() );
	}

	for ( auto it = std::next(res_list.begin()); it != std::prev(res_list.end()); ++it ) {
		if ( *std::next(it) - *it != 1 && *it - *std::prev(it) != 1 ) {
			res_remove.insert( *it );
		}
	}
	if ( *std::prev(res_list.end()) - *std::prev(res_list.end(),2) != 1 ) {
		res_remove.insert( *std::prev(res_list.end()) );
	}

	// remove res in res_remove from res_list
	remove_set1_elements_from_set2( res_remove, res_list );

	// add some new residues
	std::set< Size > new_res;
	for ( auto it = res_list.begin(); it != std::prev(res_list.end()); ++it ) {
		Size const gap = *std::next(it) - *it;
		// fill in gaps of length less than four
		if ( gap > 1 && gap <= 4 ) {
			for ( Size n = *it+1; n <= *std::next(it)+1; ++n ) {
				new_res.insert( n );
			}
		}
	}

	Size const gap = total_res - *std::prev(res_list.end());
	if ( gap > 0 && gap <= 4 ) {
		for ( Size n = *std::prev(res_list.end()) + 1; n <= total_res; ++n ) {
			new_res.insert( n );
		}
	}
	res_list.insert( new_res.begin(), new_res.end() );
}

std::set< Size >
find_nearby_res(
	Pose const & pose,
	std::set< Size > res_list_current,
	Real const dist_cutoff
) {
	std::set< Size > res_list;
	for ( auto const i : res_list_current ) {
		//TR << "Evaluating residue " << i << " from res_list_current " << std::endl;
		std::string const dist_atom_i = pose.residue_type( i ).is_protein() ? "CA" : "C1'";
		if ( !pose.residue_type( i ).has( dist_atom_i ) ) {
			// defect
			continue;
		}
		for ( Size j = 1; j <= pose.total_residue(); ++j ) {
			if ( res_list_current.find(j) != res_list_current.end() || res_list.find(j) != res_list.end() ) continue;
			std::string const dist_atom_j = pose.residue_type( j ).is_protein() ? "CA" : "C1'";
			if ( !pose.residue_type( j ).has( dist_atom_j ) ) {
				// defect
				continue;
			}
			Real const dist_C1 = pose.residue( i ).xyz( dist_atom_i ).distance_squared(
				pose.residue( j ).xyz( dist_atom_j ) );

			if ( dist_C1 > ( dist_cutoff + 8 ) * ( dist_cutoff + 8 ) ) continue;

			// Now that we've passed the neighborhood filter, check all vs all.
			for ( Size atom_i = 1; atom_i <= pose.residue_type(i).natoms(); ++atom_i ) {
				bool found_qualifying_atom = false;
				for ( Size atom_j = 1; atom_j <= pose.residue_type(j).natoms(); ++atom_j ) {
					Real const dist_sq = pose.residue( i ).xyz( atom_i ).distance_squared(
						pose.residue( j ).xyz( atom_j ) );
					if ( dist_sq < dist_cutoff * dist_cutoff ) {
						res_list.insert(j);
						//TR << " Now res_list: [ ";
						//for ( auto const elem : res_list ) TR << elem << " ";
						//TR << "]" << std::endl;
						found_qualifying_atom = true;
						break;
					}
				}
				if ( found_qualifying_atom ) break;
			}

		}
	}
	return res_list;
}

void
erase_if_in_any_slice( utility::vector1< std::set< Size > > const & res_list_sliced, Size const res, std::set< Size > & res_list_new
) {
	for ( auto & res_set : res_list_sliced ) {
		if ( res_set.find( res ) != res_set.end() ) {
			res_list_new.erase( res_list_new.find( res ) );
		}
	}
}

void clean_res_list ( std::set< Size > & res_list_new, utility::vector1< std::set< Size > > const & res_list_sliced ) {
	std::set< Size > clean_res_list_new;
	for ( auto const & res_list : res_list_sliced ) {
		std::set_difference( res_list_new.begin(), res_list_new.end(), res_list.begin(), res_list.end(), inserter(clean_res_list_new,clean_res_list_new.begin()) );
		TR << "cleanres_list_new: [ ";
		for ( auto const elem : clean_res_list_new ) TR << elem << " ";
		TR << "]" << std::endl;
		res_list_new = clean_res_list_new;
		clean_res_list_new.clear();
	}
}


void
identify_chunks(
	Pose const & pose,
	utility::vector1< std::set< Size > > & sliced_list_final,
	Size const virtual_res_pos
) {
	TR << "Identifying chunks..." << std::endl;
	Size const total_res = pose.size() - 1;
	if ( total_res <= 150 ) {
		// All in one.
		std::set< Size > res_set;
		for ( Size i = 1; i <= total_res; ++i ) {
			if ( i != virtual_res_pos ) res_set.insert( i );
		}
		sliced_list_final.push_back( res_set );
		return;
	}

	TR << "Found " << total_res << " residues." << std::endl;
	Size const n_chunk = pose.size() / 100;
	Size n_chunk_left = n_chunk;
	Size chunk_size = 0;
	std::set< Size > res_list_unsliced;
	for ( Size i = 1; i <= total_res; ++i ) res_list_unsliced.insert(i);

	utility::vector1< std::set< Size > > res_list_sliced;
	std::set< Size > res_list_current;
	while ( res_list_unsliced.size() != 0 ) {
		TR.Trace << "Unsliced: " << res_list_unsliced.size() << std::endl;
		TR.Trace << "Current:  " << res_list_current.size() << std::endl;
		for ( auto const & sl : sliced_list_final ) {
			TR.Trace << "One slice:  " << sl.size() << std::endl;
		}

		if ( res_list_current.size() == 0 ) {
			// python pop 0 to res_list_current
			Size res = *res_list_unsliced.begin();
			TR.Trace << "Starting new chunk from scratch with res " << res << std::endl;
			chunk_size = res_list_unsliced.size() / n_chunk_left;
			res_list_unsliced.erase(res_list_unsliced.begin());
			res_list_current.insert( res );
			TR.Trace << "Objective is to get " << chunk_size << " res in this chunk." << std::endl;
			n_chunk_left -= 1;
		}

		std::set< Size > res_list_new = find_nearby_res(pose, res_list_current, 3.5 );
		TR.Trace << "Adding " << res_list_new.size() << " new residues near res_list_current" << std::endl;
		TR.Trace << "To wit, res_list_new: [ ";
		for ( auto const elem : res_list_new ) TR.Trace << elem << " ";
		TR.Trace << "]" << std::endl;

		// Remove all residues previously in a res list from res_list_new
		clean_res_list( res_list_new, res_list_sliced );

		if ( res_list_new.size() == 0 && res_list_current.size() < chunk_size * 0.7 ) {
			TR.Trace << "No nearby new residues identified, but the current res list size " <<  res_list_current.size();
			TR.Trace << " is still less than chunk_size * 0.7 " << chunk_size * 0.7 << std::endl;
			while ( true ) {
				// "pop 0th element to res"
				Size res = *res_list_unsliced.begin();
				TR.Trace << "So, we pop " << res << std::endl;
				res_list_unsliced.erase(res_list_unsliced.begin());
				if ( res_list_current.find( res ) == res_list_current.end() ) {
					res_list_new.insert(res);
					break;
				}
			}
		}

		TR.Trace << "Inserting res_list_new ( " << res_list_new.size() << " residues ) into current " << std::endl;
		res_list_current.insert( res_list_new.begin(), res_list_new.end() );

		TR.Trace << "We have obtained a residue list of " << res_list_current.size() << " and our target chunk size is " << chunk_size << std::endl;
		if ( res_list_current.size() >= chunk_size || res_list_new.size() == 0 ) {
			TR.Trace << "We have exceeded chunk size with res_list_current or failed at res_list_new" << std::endl;

			TR.Trace << "res_list_current: [ ";
			for ( auto const elem : res_list_current ) TR << elem << " ";
			TR.Trace << "]" << std::endl;

			std::set< Size > unused;
			fill_gaps_and_remove_isolated_res( res_list_current, total_res, unused );

			TR.Trace << "Gonna remove res_list_current from res_list_unsliced" << std::endl;
			for ( auto const & res : res_list_current ) {
				//TR << "Trying to remove " << res << std::endl;
				if ( res_list_unsliced.find(res) != res_list_unsliced.end() ) {
					//TR << "Found" << res << std::endl;
					res_list_unsliced.erase(res_list_unsliced.find(res));
					//TR << "So the unsliced length is now " << res_list_unsliced.size()  << std::endl;
				}
			}

			TR.Trace << "Gonna add this chunk to my lists of chunks" << std::endl;
			res_list_sliced.push_back( res_list_current );
			sliced_list_final.push_back( res_list_current );
			if ( n_chunk_left == 0 ) {
				return;
			} else if ( n_chunk_left == 1 ) {
				TR.Trace << "What remains is: ";
				for ( Size const res : res_list_unsliced ) TR << res << " ";
				//if ( res_list_unsliced.size() == 0 ) break;
				TR.Trace << std::endl;
				res_list_current = res_list_unsliced;
				std::set< Size > removed_res;
				fill_gaps_and_remove_isolated_res(res_list_current, total_res, removed_res);
				sliced_list_final.push_back( res_list_current );
				TR.Trace << "Adding res_list_current to the list of sliced_list_final" << std::endl;
				while ( removed_res.size() != 0 ) {
					TR.Trace << "Trying to remove res from each sliced list, to some end" << std::endl;
					std::set< Size > res_remove;
					for ( Size const res : removed_res ) {
						std::set< Size > just_res;
						just_res.insert( res );
						std::set< Size > res_list_near = find_nearby_res( pose, just_res, 2.0 );
						for ( auto & sliced_list : sliced_list_final  ) {
							bool can_break = false;
							for ( Size const near_res : res_list_near ) {
								if ( sliced_list.find( near_res ) != sliced_list.end() ) {
									sliced_list.insert( res );
									res_remove.insert( res );
									can_break = true;
									break;
								}
							}
							if ( can_break ) break;
						}
					}
					for ( auto const res : res_remove ) {
						removed_res.erase(removed_res.find(res));
					}
				}
				break;
			} else {
				TR.Trace << "Ready for new beginnings." << std::endl;
				res_list_current.clear();
				res_list_current.insert(*res_list_unsliced.begin());
				res_list_unsliced.erase(res_list_unsliced.begin());
				chunk_size = Size( res_list_unsliced.size() / n_chunk_left * 1.1 );
				n_chunk_left -= 1;
			}
		}
	}
}

void
ErraserMinimizerMover::turn_off_for_chunks(
	MoveMap & mm,
	Pose const & pose,
	std::set< Size > const & chunk
) {
	// Turn off all DOFs that ARE NOT in chunk_i
	for ( Size i = 1; i <= pose.total_residue(); ++i ) {
		if ( chunk.find( i ) == chunk.end() ) {
			mm.set_chi( i, false );
			mm.set_bb(  i, false );
			for ( Size j = 1; j <= pose.residue_type( i ).natoms(); ++j ) {
				mm.set( DOF_ID( AtomID( j, i ), D ),     false );
				mm.set( DOF_ID( AtomID( j, i ), THETA ), false );
				mm.set( DOF_ID( AtomID( j, i ), PHI ),   false );
			}
		}
	}

	// AMW TODO: Turn ON all DOFs that cross between chunk_i and the universe
}

core::Vector com_of_true_residues( kinematics::MoveMap const & mm, Pose const & pose ) {
	// Just average heavyatom position
	core::Vector avg( 0 );
	Size natoms = 0;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( !mm.get_bb( ii ) && !mm.get_chi( ii ) ) continue;

		for ( Size jj = 1; jj <= pose.residue_type( ii ).nheavyatoms(); ++jj ) {
			avg += pose.residue( ii ).xyz( jj );
			natoms += 1;
		}
	}
	return avg / Real(natoms);
}


///////////////////////////////////////////

// Main workhorse function
void
ErraserMinimizerMover::apply(
	Pose & pose
) {
	core::pose::rna::make_phosphate_nomenclature_matches_mini( pose );

	if ( ready_set_only_ ) return;

	// If this pose 'comes with' a particularly opinionated FT, don't change it.
	if ( pose.fold_tree().num_cutpoint() == 0 ) {
		// Setup fold tree using user input or using Rhiju's function
		if ( cutpoint_list_.size() == 0 ) {
			core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
		} else {
			setup_fold_tree( pose );
		}
	}
	//TR << pose.fold_tree();

	// Add a virtual residue for density scoring
	Size const virtual_res_pos = add_virtual_res( pose );
	pose::Pose const pose_full = pose;
	Size const nres( pose.size() );
	Size const nres_moving( nres - fixed_res_list_.size() );

	// Output the sequence
	std::string const & working_sequence = pose.sequence();
	TR << "Pose sequence = " << working_sequence << std::endl;

	// Try flipping the pyrimidines
	if ( attempt_pyrimidine_flip_ ) {
		TR << "Do we have to flip pyrimidines?" << std::endl;
		pyrimidine_flip_trial( pose );
	}

	// If we have more than 150 residues, we need to "slice." This involves
	// designing ~100-residue chunks of the pose and making only those DOFs
	// visible to the minimizer. All residues of the pose will be present.
	// Each minimization gets its own little output file. Notably, while we
	// must figure out the chunks every time, we can skip chunks for which
	// an appropriate output file exists.
	utility::vector1< std::set< Size > > chunks;
	// AMW TODO: read chunks from temp file if exists
	identify_chunks( pose, chunks, virtual_res_pos );
	Size const n_chunk = chunks.size();
	TR << "Identified " << n_chunk << " chunks" << std::endl;
	for ( Size ii = 1; ii <= n_chunk; ++ii ) {
		TR << "[";
		for ( core::Size jj : chunks[ii] ) {
			TR << " " << jj;
		}
		TR << "]" << std::endl;
	}

	Size first_chunk = 1;
	for ( ; first_chunk <= n_chunk; ++first_chunk ) {
		std::stringstream inname;
		inname << "full_minimize_temp_" << first_chunk << ".out";
		std::ifstream f( inname.str() );
		if ( !f.good() ) break;
		bool found_token = false;
		while ( f.good() ) {
			std::string tok;
			f >> tok;
			if ( tok == "DONE!" ) {
				found_token = true;
			}
		}
		if ( !found_token ) {
			// this out is bad
			//first_chunk--;
			break;
		}
	}

	//Set the MoveMap, avoiding moving the virtual residue
	TR << "Setting up movemap ..." << std::endl;
	kinematics::MoveMap mm;
	ObjexxFCL::FArray1D< bool > allow_insert( nres, false );
	mm.set_bb( false );
	mm.set_chi( false );
	mm.set_jump( false );

	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		if ( pose.residue_type( ii ).aa() == core::chemical::aa_vrt ) continue;

		allow_insert( ii ) = true;
		mm.set_bb(  ii, true );
		mm.set_chi( ii, true );
	}

	kinematics::FoldTree fold_tree( pose.fold_tree() );

	std::set< core::Size > cut_upper, cut_lower;
	for ( Size i = 1; i <= fold_tree.num_jump(); ++i ) {
		Size const ustream   = fold_tree.upstream_jump_residue( i );
		Size const dstream = fold_tree.downstream_jump_residue( i );
		TR << "Considering jump number " << i << " from " << ustream << " to " << dstream << std::endl;
		cut_lower.insert( ustream );
		cut_upper.insert( dstream );

		if ( pose.residue_type( ustream ).aa() == core::chemical::aa_vrt ) continue;
		if ( pose.residue_type( dstream ).aa() == core::chemical::aa_vrt ) continue;

		// Skip MG-water jumps. We do these poorly.
		if ( ( pose.residue_type( ustream ).aa() == core::chemical::aa_h2o &&
				pose.residue_type( dstream ).name() == "MG" ) ||
				( pose.residue_type( dstream ).aa() == core::chemical::aa_h2o &&
				pose.residue_type( ustream ).name() == "MG" ) ) continue;

		TR << "not all virts..." << std::endl;
		if ( fixed_res_list_.empty() || ( fixed_res_list_.size() != 0 &&
				fixed_res_list_.find( ustream ) == fixed_res_list_.end() &&
				fixed_res_list_.find( dstream ) == fixed_res_list_.end() ) ) {
			mm.set_jump( i, true );

			TR << "not fixed_res! " << fixed_res_list_ << std::endl;
		}
	}

	//Fixed res mode
	if ( fixed_res_list_.size() != 0 ) {
		scorefxn_->set_weight( coordinate_constraint, 10 );

		TR << "fixed res: ";
	}

	for ( Size const fixed_res_num : fixed_res_list_ ) {
		TR << fixed_res_num << " ";

		add_fixed_res_constraints( pose, fixed_res_num, virtual_res_pos );

		mm.set_chi( fixed_res_num, false );
		mm.set_bb(  fixed_res_num, false );

		allow_insert( fixed_res_num ) = false;

		if ( fixed_res_num - 1 > 0 &&
				fixed_res_list_.find( fixed_res_num ) == fixed_res_list_.end() &&
				cut_lower.find( fixed_res_num ) == cut_lower.end() ) {
			allow_insert( fixed_res_num ) = true;
		}
		if ( fixed_res_num + 1 <= nres &&
				fixed_res_list_.find( fixed_res_num + 1 ) == fixed_res_list_.end() &&
				cut_upper.find( fixed_res_num ) == cut_upper.end() ) {
			allow_insert( fixed_res_num ) = true;
		}
	}
	TR << std::endl;

	// Handle phosphate constraints
	if ( constrain_phosphate_ ) {
		scorefxn_->set_weight( coordinate_constraint, 10 );
		ConstraintSetOP cst_set = pose.constraint_set()->clone();

		for ( Size i = 1; i <= nres; ++i ) {
			if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue;
			// Fixed res phosphates can't move anyway, so don't bother.
			if ( fixed_res_list_.find( i ) != fixed_res_list_.end() ) continue;

			Real const coord_sdev( 0.3 );
			Size const my_anchor( virtual_res_pos ); //anchor on virtual residue
			Residue const & rsd( pose.residue( i ) );
			// P for RNA, CA for protein (if they sneak in)
			if ( !rsd.has( "P" ) && !rsd.has( "CA" ) ) continue;
			Size const atm_indexP = rsd.has( "P" ) ? rsd.atom_index( "P" ) : rsd.atom_index( "CA" );
			cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
				AtomID( atm_indexP, i ),
				AtomID( 1, my_anchor ), rsd.xyz( atm_indexP ),
				FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) );
		}
		pose.constraint_set( cst_set );
	}
	ConstraintSetOP saved_cst_set = pose.constraint_set()->clone();

	for ( Size chunk_i = first_chunk; chunk_i <= n_chunk; ++chunk_i ) {
		time_t chunk_start = time(0);

		// Don't retain those dirty constraints from chunk n-1 -- they don't
		// matter, since those residues can't move!
		pose.constraint_set( saved_cst_set );
		std::stringstream outname;
		outname << "full_minimize_temp_" << chunk_i << ".out";
		std::ofstream out( outname.str() );
		out << "Starting chunk " << chunk_i << "..." << std::endl;

		// Constrain bonded atom sets to the starting geometry
		kinematics::MoveMap chunk_mm( mm );
		if ( vary_bond_geometry_ ) {
			vary_bond_geometry( chunk_mm, pose, allow_insert, chunks[ chunk_i ] );
		}
		//Ensure OFF for chunks - probably not needed
		if ( n_chunk != 1 ) {
			turn_off_for_chunks( chunk_mm, pose, chunks[chunk_i] );
		}
		protocols::stepwise::modeler::output_movemap( chunk_mm, pose );
		scorefxn_->show( TR, pose );
		Real const score_before = ( ( *scorefxn_ )( pose ) );
		Real const edens_score_before = ( ( *edens_scorefxn_ )( pose ) );

		core::Vector chunk_com = com_of_true_residues( chunk_mm, pose );
		protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400, false, true, chunk_com );

		// Start Minimizing the Full Structure
		Pose const start_pose = pose;
		AtomTreeMinimizer minimizer;
		float const dummy_tol( 0.0001 );
		Size const min_iter = std::min( 5000, std::max( 2000, int( nres_moving * 24 ) ) );

		TR << "Minimize using dfpmin with use_nb_list=true .." << std::endl;
		MinimizerOptions min_options_dfpmin( "lbfgs_armijo_nonmonotone" /*option[ run::min_type ]*/, dummy_tol, true, false, false );
		min_options_dfpmin.max_iter( min_iter );
		minimizer.run( pose, chunk_mm, *scorefxn_, min_options_dfpmin );

		scorefxn_->show( TR, pose );
		Real const score = ( ( *scorefxn_ )( pose ) );
		Real const edens_score = ( ( *edens_scorefxn_ )( pose ) );
		TR << "current_score = " << score << ", start_score = " << score_before << std::endl;
		TR << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
		/*
		if ( score > score_before + 5 || edens_score > edens_score_before * 0.9 ) {
		TR << "current_score = " << score << ", start_score = " << score_before << std::endl;
		TR << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
		TR << "The minimization went wild!!! Try alternative minimization using dfpmin with use_nb_list=false .." << std::endl;

		pose = start_pose;

		MinimizerOptions min_options_dfpmin_no_nb( "lbfgs_armijo_nonmonotone", dummy_tol, false, false, false );
		min_options_dfpmin_no_nb.max_iter( min_iter );
		minimizer.run( pose, chunk_mm, *scorefxn_, min_options_dfpmin_no_nb );
		scorefxn_->show( std::cout, pose );
		Real const score = ( ( *scorefxn_ ) ( pose ) );
		Real const edens_score = ( ( *edens_scorefxn_ ) ( pose ) );
		if ( score > score_before + 5 || edens_score > edens_score_before * 0.9 ) {
		TR << "current_score = " << score << ", start_score = " << score_before << std::endl;
		TR << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
		pose = start_pose;
		TR << "The minimization went wild again!!! Skip the minimization!!!!!" << std::endl;
		}
		}
		*/
		time_t chunk_end = time(0);
		out << "CPU time for chunk: " << (chunk_end-chunk_start) << " seconds." << std::endl;
		out << "DONE!" << std::endl;
		std::stringstream fn;
		fn << "debug_chunk_" << chunk_i << ".pdb";
		pose.pdb_info()->obsolete( false );

		pose.dump_pdb( fn.str() );

		protocols::viewer::clear_conformation_viewers();
	}
	pose.pdb_info()->obsolete( false );

	TR << "Job completed sucessfully." << std::endl;

	// Remove slice output files
	for ( Size chunk_i = 1; chunk_i <= n_chunk; ++chunk_i ) {
		std::stringstream outname;
		outname << "full_minimize_temp_" << chunk_i << ".out";
		utility::file::file_delete( outname.str() );
	}
}

// XRW TEMP std::string
// XRW TEMP ErraserMinimizerMoverCreator::keyname() const
// XRW TEMP {
// XRW TEMP  return ErraserMinimizerMover::mover_name();
// XRW TEMP }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP ErraserMinimizerMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new ErraserMinimizerMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP ErraserMinimizerMover::mover_name()
// XRW TEMP {
// XRW TEMP  return "ErraserMinimizerMover";
// XRW TEMP }

std::string ErraserMinimizerMover::get_name() const {
	return mover_name();
}

std::string ErraserMinimizerMover::mover_name() {
	return "ErraserMinimizerMover";
}

void ErraserMinimizerMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "vary_bond_geometry", xsct_rosetta_bool, "Vary bond lengths and angles, constrained to ideal values." )
		+ XMLSchemaAttribute( "constrain_phosphate", xsct_rosetta_bool, "Constrain phosphates to their initial positions." )
		+ XMLSchemaAttribute( "ready_set_only", xsct_rosetta_bool, "Do nothing but ensure phosphate nomenclature is correct: useful for PHENIX integration and little else." )
		+ XMLSchemaAttribute( "skip_minimize", xsct_rosetta_bool, "Skip the minimization step -- only do pyrimidine rotamer trials (if enabled)." )
		+ XMLSchemaAttribute( "attempt_pyrimidine_flip", xsct_rosetta_bool, "Rotamer trials on pyrimidine bases." )
		+ XMLSchemaAttribute( "fixed_res_list", xsct_int_wsslist, "Whitespace-separated list of sequence positions in Rosetta numbering." )
		+ XMLSchemaAttribute( "cutpoint_list", xsct_int_wsslist, "Whitespace-separated list of sequence positions in Rosetta numbering indicating cutpoints." )
		+ XMLSchemaAttribute( "output_pdb_name", xs_string, "Output a PDB to this file, subverting JD2." );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Optimize an RNA in the presence of electron density from X-ray or cryo-EM: minimization phase", attlist );
}

std::string ErraserMinimizerMoverCreator::keyname() const {
	return ErraserMinimizerMover::mover_name();
}

protocols::moves::MoverOP
ErraserMinimizerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ErraserMinimizerMover );
}

void ErraserMinimizerMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ErraserMinimizerMover::provide_xml_schema( xsd );
}


} //movers
} //rna
} //protocols
