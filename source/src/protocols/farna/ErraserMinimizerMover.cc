// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/farna/ErraserMinimizerMover.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/farna/ErraserMinimizerMover.hh>
#include <protocols/farna/ErraserMinimizerMoverCreator.hh>

#include <core/pose/rna/RNA_IdealCoord.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/output_util.hh>
// libRosetta headers
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
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

#include <protocols/viewer/viewers.hh>

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
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
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
#include <numeric/conversions.hh>
#include <numeric/xyz.functions.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.farna.ErraserMinimizerMover" );

using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace basic::options::OptionKeys::rna::farna::erraser;
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
using namespace protocols::farna;
using namespace protocols::filters;
using utility::vector1;

namespace protocols {
namespace farna {


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
{}

void ErraserMinimizerMover::initialize_from_options() {
	vary_bond_geometry_ = option[ basic::options::OptionKeys::rna::vary_geometry ].value();
	constrain_phosphate_ = option[ constrain_P ].value();
	ready_set_only_ =     option[ ready_set_only ].value();
	skip_minimize_ =      option[ skip_minimize ].value();
	attempt_pyrimidine_flip_ = option[ attempt_pyrimidine_flip ].value();
	std::copy( option[ fixed_res ].value().begin(), option[ fixed_res ].value().end(), std::inserter( fixed_res_list_, fixed_res_list_.end() ) );
	cutpoint_list_ = option[ full_model::cutpoint_open ].value();
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
}

void
ErraserMinimizerMover::create_pose_reference( Pose const & pose ) {
	core::chemical::ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	make_pose_from_sequence( pose_reference_, pose.annotated_sequence(), *rsd_set );
	apply_ideal_coordinates( pose );
}

bool
ErraserMinimizerMover::check_in_bonded_list(
	core::id::AtomID const & atom_id1,
	core::id::AtomID const & atom_id2
) {
	for ( auto & bond : bonded_atom_list_ ) {
		if ( atom_id1 == bond.first && atom_id2 == bond.second ) return true;
		if ( atom_id2 == bond.first && atom_id1 == bond.second ) return true;
	}

	return false;
}

bool
ErraserMinimizerMover::check_in_bond_angle_list(
	core::id::AtomID const & central_atom,
	core::id::AtomID const & side_one,
	core::id::AtomID const & side_two
) {
	for ( auto & angle : bond_angle_list_ ) {
		if ( central_atom != angle.first ) continue;

		if ( side_one == angle.second.first && side_two == angle.second.second ) return true;
		if ( side_two == angle.second.first && side_one == angle.second.second ) return true;
	}

	return false;
}

void
ErraserMinimizerMover::apply_ideal_coordinates( pose::Pose const & pose ) {
	RNA_FittedTorsionInfo const rna_fitted_torsion_info;
	Real const DELTA_CUTOFF( rna_fitted_torsion_info.delta_cutoff() );
	bool const use_phenix_geo = option[ basic::options::OptionKeys::rna::corrected_geo ];
	utility::vector1< PuckerState > pucker_conformation( pose_reference_.size(), NO_PUCKER );

	RNA_IdealCoord ideal_coord;
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( pose.residue_type( n ).is_virtual_residue() ) continue;

		Real const delta = pose.residue( n ).mainchain_torsion( DELTA );

		if ( delta > DELTA_CUTOFF ) { //south
			apply_ideal_c2endo_sugar_coords( pose_reference_, n );
			pucker_conformation[n] = SOUTH;
		} else { //north
			pucker_conformation[n] = NORTH;
		}
	}
	if ( use_phenix_geo ) {
		ideal_coord.apply( pose_reference_, pucker_conformation, false );
	}
}

void
ErraserMinimizerMover::add_bond_constraint(
	AtomID const & atom_id1,
	AtomID const & atom_id2,
	core::pose::Pose const & pose,
	core::scoring::constraints::ConstraintSetOP & cst_set
) {
	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );

	if ( !pose_reference_.residue_type( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference_.residue_type( atom_id2.rsd() ).has( atom_name2 ) ) return;

	// already addressed, just return
	if ( check_in_bonded_list( atom_id1, atom_id2 ) ) return;

	bonded_atom_list_.push_back( std::make_pair( atom_id1, atom_id2 ) );
	Real const bond_length_sd_( 0.05 );
	Real const bond_length = ( pose_reference_.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
		pose_reference_.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();
	func::FuncOP dist_harm_func_( new core::scoring::func::HarmonicFunc( bond_length, bond_length_sd_ ) );
	cst_set->add_constraint( ConstraintCOP( ConstraintOP( new AtomPairConstraint( atom_id1, atom_id2, dist_harm_func_, rna_bond_geometry ) ) ) );

	TR.Trace << "PUTTING CONSTRAINT ON DISTANCE: " <<
		atom_id2.rsd() << " " << atom_name1 << "; "  <<
		atom_id1.rsd() << " " << atom_name2 << " "  <<
		bond_length << std::endl;
}

void
ErraserMinimizerMover::add_bond_angle_constraint(
	AtomID const & atom_id1,
	AtomID const & atom_id2,
	AtomID const & atom_id3,
	core::pose::Pose const & pose,
	core::scoring::constraints::ConstraintSetOP & cst_set
) {
	if ( atom_id2 == atom_id3 ) return;
	if ( atom_id1 == atom_id3 ) return;
	if ( atom_id1 == atom_id2 ) return;

	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue_type( atom_id3.rsd() ).atom_name( atom_id3.atomno() );

	if ( !pose_reference_.residue_type( atom_id1.rsd() ).has( atom_name1 ) ) return;
	if ( !pose_reference_.residue_type( atom_id2.rsd() ).has( atom_name2 ) ) return;
	if ( !pose_reference_.residue_type( atom_id3.rsd() ).has( atom_name3 ) ) return;

	// if already handled, return
	if ( check_in_bond_angle_list( atom_id1, atom_id2, atom_id3 ) ) return;

	bond_angle_list_.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );
	Real const bond_angle_sd_( radians ( 3.0 ) );
	Real const bond_angle = angle_radians(
		pose_reference_.residue( atom_id2.rsd() ).xyz( atom_name2 ),
		pose_reference_.residue( atom_id1.rsd() ).xyz( atom_name1 ),
		pose_reference_.residue( atom_id3.rsd() ).xyz( atom_name3 ) );

	if ( bond_angle < 0.001 ) TR.Warning << "WHAT THE HELL????????? " << std::endl;

	cst_set->add_constraint( ConstraintCOP( new AngleConstraint(
		atom_id2, atom_id1, atom_id3,
		FuncOP( new HarmonicFunc( bond_angle, bond_angle_sd_ ) ),
		rna_bond_geometry ) ) );

	TR.Trace << "PUTTING CONSTRAINT ON ANGLE: "
		<< atom_id2.rsd() << " " << atom_name2 << "; "
		<< atom_id1.rsd() << " " << atom_name1 << "; "
		<< atom_id3.rsd() << " " << atom_name3 << " ==> "
		<< degrees( bond_angle ) << " " << degrees( bond_angle_sd_ ) << std::endl;
}

bool
ErraserMinimizerMover::check_if_connected_in_atom_tree(
	core::pose::Pose const & pose,
	AtomID const & atom_id1,
	AtomID const & atom_id2
) {
	if ( atom_id1.rsd() == atom_id2.rsd() ) return true;

	core::kinematics::tree::AtomCOP atom1( pose.atom_tree().atom( atom_id1 ).get_self_ptr() );
	core::kinematics::tree::AtomCOP atom2( pose.atom_tree().atom( atom_id2 ).get_self_ptr() );

	if ( atom1->parent() == atom2 ) return true;
	if ( atom2->parent() == atom1 ) return true;

	return false;
}

// Virts and sidechain atoms that aren't the first base atom should not move
bool
ErraserMinimizerMover::i_want_this_atom_to_move(
	chemical::ResidueType const & residue_type,
	Size const & atomno
) {
	if ( atomno > residue_type.first_sidechain_atom() &&
			atomno != chemical::rna::first_base_atom_index( residue_type ) ) return false;

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
	ObjexxFCL::FArray1D< bool > & allow_insert // Operationally: not fixed, cutpoint, virt
) {
	ConstraintSetOP cst_set = pose.constraint_set()->clone();

	Size const nres( pose.size() );
	TR << "Enter vary_bond_geometry....." << std::endl;

	// First, set appropriate DOFs to move in the movemap, mm
	for ( Size i = 1; i <= nres; ++i )  {
		if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue; //FCC
		if ( !allow_insert( i ) ) continue;

		chemical::ResidueType const & residue_type( pose.residue_type( i ) );

		for ( Size j = 1; j <= residue_type.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move( residue_type, j ) ) continue;
			if ( !does_atom_exist_in_reference( pose, AtomID( j, i ) ) )  continue;

			tree::AtomCOP current_atom( pose.atom_tree().atom( AtomID( j, i ) ).get_self_ptr() );
			if ( current_atom->is_jump() ) continue;

			// STEP ONE: DISTANCES
			tree::AtomCOP j_distatom( current_atom->input_stub_atom1() );

			if ( !j_distatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_distatom->id() ) ) continue;
			if ( !does_atom_exist_in_reference( pose, j_distatom->id() ) ) continue;

			mm.set( DOF_ID( AtomID( j, i ), D ), true );

			if ( j_distatom->is_jump() ) continue;

			// STEP TWO: ANGLES
			core::kinematics::tree::AtomCOP j_angleatom( current_atom->input_stub_atom2() );

			if ( !j_angleatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_angleatom->id() ) ) continue;
			if ( !does_atom_exist_in_reference( pose, j_angleatom->id() ) ) continue;
			if ( j_angleatom == current_atom ) continue;

			mm.set( DOF_ID( AtomID( j, i ), THETA ), true );

			if ( j_angleatom->is_jump() ) continue;

			// STEP THREE: DIHEDRALS
			core::kinematics::tree::AtomCOP j_dihedralatom( current_atom->input_stub_atom3() );

			if ( !j_dihedralatom ) continue;
			if ( !i_want_this_atom_to_move( pose, j_dihedralatom->id() ) ) continue;
			if ( !does_atom_exist_in_reference( pose, j_dihedralatom->id() ) ) continue;
			if ( j_dihedralatom == current_atom ) continue;

			mm.set( DOF_ID( AtomID( j, i ), PHI ), true );
		}
	}

	// Next, build lists on bonds and bond angles and restrain them to present values.
	for ( Size i = 1; i <= nres; i++ )  {
		if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue; //FCC
		if ( !allow_insert( i ) ) continue;

		chemical::ResidueType const & residue_type( pose.residue_type( i ) );


		for ( Size j = 1; j <= residue_type.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move( residue_type, j ) ) continue;

			AtomID j_atomid( j, i );
			utility::vector1< AtomID > nbrs( pose.conformation().bonded_neighbor_all_res( j_atomid ) );

			// Constrain all mobile bond lengths.
			for ( auto & nbr : nbrs ) {
				if ( nbr.rsd() > nres || nbr.rsd() < 1 ) continue;
				chemical::ResidueType const & residue_type2( pose.residue_type( nbr.rsd() ) );
				Size const & k( nbr.atomno() );

				if ( ! check_if_connected_in_atom_tree( pose, j_atomid, nbr ) ) continue;

				if ( i_want_this_atom_to_move( residue_type2, k ) )  {
					add_bond_constraint( j_atomid, nbr, pose, cst_set );
				}
			}

			// Bond angles
			for ( auto nbr = nbrs.begin(); nbr != nbrs.end(); ++nbr ) {
				if ( nbr->rsd() > nres || nbr->rsd() < 1 ) continue;

				chemical::ResidueType const & residue_type2( pose.residue_type( nbr->rsd() ) );

				if ( ! check_if_connected_in_atom_tree( pose, j_atomid, *nbr ) ) continue;

				Size const & k( nbr->atomno() ) ;

				for ( auto & ang_nbr : nbrs ) {
					if ( ang_nbr.rsd() > nres || ang_nbr.rsd() < 1 ) continue;
					chemical::ResidueType const & residue_type3( pose.residue_type( ang_nbr.rsd() ) );

					if ( !check_if_connected_in_atom_tree( pose, j_atomid, ang_nbr ) ) continue;

					Size const & q( ang_nbr.atomno() ) ;

					if ( i_want_this_atom_to_move( residue_type2, k ) &&
							i_want_this_atom_to_move( residue_type3, q ) )  {
						add_bond_angle_constraint( j_atomid, *nbr, ang_nbr,
							pose, cst_set );
					}
				}

				utility::vector1< AtomID > nbrs2( pose.conformation().bonded_neighbor_all_res( *nbr ) );

				for ( auto & nbr2 : nbrs2 ) {
					if ( nbr2.rsd() > nres || nbr2.rsd() < 1 ) break;
					chemical::ResidueType const & residue_type3( pose.residue_type( nbr2.rsd() ) );

					if ( ! check_if_connected_in_atom_tree( pose, *nbr, nbr2 ) ) continue;

					Size const & q( nbr2.atomno() ) ;

					if ( i_want_this_atom_to_move( residue_type2, k ) &&
							i_want_this_atom_to_move( residue_type3, q ) )  {
						add_bond_angle_constraint( j_atomid, *nbr, nbr2,
							pose, cst_set );
					}
				}
			}
		}
	}

	pose.constraint_set( cst_set );
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

	// attach virt res there
	core::chemical::ResidueTypeSet const & residue_set = *pose.residue_type( 1 ).residue_type_set();
	core::conformation::ResidueOP new_res( core::conformation::ResidueFactory::create_residue( *( residue_set.get_representative_type_name3( "VRT" ) ) ) );
	pose.append_residue_by_jump( *new_res , nres );
	// make the virt atom the root
	kinematics::FoldTree newF( pose.fold_tree() );
	newF.reorder( nres + 1 );
	pose.fold_tree( newF );

	core::pose::full_model_info::FullModelInfoOP full_model_info( new core::pose::full_model_info::FullModelInfo( pose ) );
	set_full_model_info( pose, full_model_info );

	return nres + 1;
}

bool
ErraserMinimizerMover::does_atom_exist_in_reference(
	pose::Pose const & pose,
	core::id::AtomID const & atom_id
) {
	std::string const & atom_name = pose.residue_type( atom_id.rsd() ).atom_name( atom_id.atomno() );

	if ( pose_reference_.residue_type( atom_id.rsd() ).has( atom_name ) ) {
		return true;
	} else {
		TR << atom_name << std::endl;
		return false;
	}
}

///////////////////////////////////////////////////////////
void
ErraserMinimizerMover::setup_fold_tree( pose::Pose & pose ) {
	Size const nres( pose.size() );
	Size const num_jumps( cutpoint_list_.size() );
	ObjexxFCL::FArray2D< int > jump_points( 2, num_jumps );
	ObjexxFCL::FArray1D< int > cuts( num_jumps );

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

///////////////////////////////////////////

// Main workhorse function
void
ErraserMinimizerMover::apply( Pose & pose ) {
	core::pose::rna::make_phosphate_nomenclature_matches_mini( pose );

	if ( ready_set_only_ ) return;

	//Setup score function.
	std::string score_weight_file = "stepwise/rna/rna_hires_elec_dens";
	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
		TR << "User passed in score:weight option: " << score_weight_file << std::endl;
	}
	scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
	edens_scorefxn_->set_weight( elec_dens_atomwise, scorefxn_->get_weight( elec_dens_atomwise ) );

	// Setup fold tree using user input or using Rhiju's function
	if ( cutpoint_list_.size() == 0 ) {
		core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );
	} else {
		setup_fold_tree( pose );
	}

	// Add a virtual residue for density scoring
	Size const virtual_res_pos = add_virtual_res( pose );
	pose::Pose const pose_full = pose;
	Size const nres( pose.size() );
	Size const nres_moving( nres - fixed_res_list_.size() );

	// Output the sequence
	std::string working_sequence = pose.sequence();
	TR << "Pose sequence = " << working_sequence << std::endl;

	// Try flipping the pyrimidines
	if ( attempt_pyrimidine_flip_ ) {
		pyrimidine_flip_trial( pose );
	}

	if ( skip_minimize_ ) return;

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
		cut_lower.insert( ustream );
		cut_upper.insert( dstream );

		if ( pose.residue_type( ustream ).aa() == core::chemical::aa_vrt ) continue;
		if ( pose.residue_type( dstream ).aa() == core::chemical::aa_vrt ) continue;

		if ( fixed_res_list_.size() != 0 &&
				fixed_res_list_.find( ustream ) == fixed_res_list_.end() &&
				fixed_res_list_.find( dstream ) == fixed_res_list_.end() ) {
			mm.set_jump ( i, true );
		}
	}

	//Fixed res mode
	if ( fixed_res_list_.size() != 0 ) {
		scorefxn_->set_weight( coordinate_constraint, 10 );

		TR << "fixed res: ";
	}

	for ( auto fixed_res_num = fixed_res_list_.begin(); fixed_res_num != fixed_res_list_.end(); ++fixed_res_num ) {

		TR << *fixed_res_num << " ";

		add_fixed_res_constraints( pose, *fixed_res_num, virtual_res_pos );

		mm.set_chi( *fixed_res_num, false );
		mm.set_bb(  *fixed_res_num, false );

		allow_insert( *fixed_res_num ) = false;

		if ( *fixed_res_num - 1 > 0 &&
				fixed_res_list_.find( *fixed_res_num ) == fixed_res_list_.end() &&
				cut_lower.find( *fixed_res_num ) == cut_lower.end() ) {
			allow_insert( *fixed_res_num ) = true;
		}
		if ( *fixed_res_num + 1 <= nres &&
				fixed_res_list_.find( *fixed_res_num + 1 ) == fixed_res_list_.end() &&
				cut_upper.find( *fixed_res_num ) == cut_upper.end() ) {
			allow_insert( *fixed_res_num ) = true;
		}
	}
	TR << std::endl;

	// Handle phosphate constraints
	if ( constrain_phosphate_ ) {
		scorefxn_->set_weight( coordinate_constraint, 10 );

		for ( Size i = 1; i <= nres; ++i ) {
			if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue;
			// Fixed res phosphates can't move anyway, so don't bother.
			if ( fixed_res_list_.find( i ) != fixed_res_list_.end() ) continue;

			Real const coord_sdev( 0.3 );
			Size const my_anchor( virtual_res_pos ); //anchor on virtual residue
			ConstraintSetOP cst_set = pose.constraint_set()->clone();
			Residue const & rsd( pose.residue( i ) );
			Size const atm_indexP = rsd.atom_index( "P" );
			cst_set->add_constraint( ConstraintCOP( new CoordinateConstraint(
				AtomID( atm_indexP, i ),
				AtomID( 1, my_anchor ), rsd.xyz( atm_indexP ),
				FuncOP( new HarmonicFunc( 0.0, coord_sdev ) ) ) ) );
			pose.constraint_set( cst_set );
		}
	}

	// Create a starting reference pose and constrain bonded atom sets to the starting geometry
	if ( vary_bond_geometry_ ) {
		TR << "Setup vary_bond_geometry" << std::endl;
		create_pose_reference( pose_full );
		vary_bond_geometry( mm, pose, allow_insert );
	}

	protocols::stepwise::modeler::output_movemap( mm, pose );
	scorefxn_->show( TR, pose );
	Real const score_before = ( ( *scorefxn_ )( pose ) );
	Real const edens_score_before = ( ( *edens_scorefxn_ )( pose ) );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	// Start Minimizing the Full Structure
	Pose const start_pose = pose;
	AtomTreeMinimizer minimizer;
	float const dummy_tol( 0.00000001 );

	TR << "Minimize using dfpmin with use_nb_list=true .." << std::endl;
	MinimizerOptions min_options_dfpmin( "lbfgs_armijo_nonmonotone", dummy_tol, true, false, false );
	min_options_dfpmin.max_iter( std::min( 3000, std::max( 1000, int(nres_moving * 12) ) ) );
	minimizer.run( pose, mm, *scorefxn_, min_options_dfpmin );

	scorefxn_->show( TR, pose );
	Real const score = ( ( *scorefxn_ )( pose ) );
	Real const edens_score = ( ( *edens_scorefxn_ )( pose ) );
	if ( score > score_before + 5 || edens_score > edens_score_before * 0.9 ) {
		TR << "current_score = " << score << ", start_score = " << score_before << std::endl;
		TR << "current_edens_score = " << edens_score << ", start_edens_score = " << edens_score_before << std::endl;
		TR << "The minimization went wild!!! Try alternative minimization using dfpmin with use_nb_list=false .." << std::endl;

		pose = start_pose;

		MinimizerOptions min_options_dfpmin_no_nb( "lbfgs_armijo_nonmonotone", dummy_tol, false, false, false );
		min_options_dfpmin_no_nb.max_iter( std::min( 3000, std::max( 1000, int(nres_moving * 12) ) ) );
		minimizer.run( pose, mm, *scorefxn_, min_options_dfpmin_no_nb );
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

	TR << "Job completed sucessfully." << std::endl;
}

std::string
ErraserMinimizerMoverCreator::keyname() const
{
	return ErraserMinimizerMoverCreator::mover_name();
}

protocols::moves::MoverOP
ErraserMinimizerMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new ErraserMinimizerMover );
}

std::string
ErraserMinimizerMoverCreator::mover_name()
{
	return "ErraserMinimizerMover";
}

} //farna
} //protocols
