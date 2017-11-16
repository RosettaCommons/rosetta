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
	// Initialize our ideal residues.
	//Names of the pdb files
	std::string const path( basic::database::full_name("chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_phenix/ideal_geometry/") );
	utility::vector1< std::string > pdb_file_list;
	pdb_file_list.push_back( path + "/A_n_std.pdb" );
	pdb_file_list.push_back( path + "/A_s_std.pdb" );
	pdb_file_list.push_back( path + "/G_n_std.pdb" );
	pdb_file_list.push_back( path + "/G_s_std.pdb" );
	pdb_file_list.push_back( path + "/C_n_std.pdb" );
	pdb_file_list.push_back( path + "/C_s_std.pdb" );
	pdb_file_list.push_back( path + "/U_n_std.pdb" );
	pdb_file_list.push_back( path + "/U_s_std.pdb" );

	chemical::ResidueTypeSetCOP rsd_set = chemical::ChemicalManager::get_instance()->residue_type_set(chemical::FA_STANDARD);
	Pose ref_pose;
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[1] );
	ideal_poses_[ "north" ][ "  A" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[2] );
	ideal_poses_[ "south" ][ "  A" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[3] );
	ideal_poses_[ "north" ][ "  G" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[4] );
	ideal_poses_[ "south" ][ "  G" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[5] );
	ideal_poses_[ "north" ][ "  C" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[6] );
	ideal_poses_[ "south" ][ "  C" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[7] );
	ideal_poses_[ "north" ][ "  U" ] = ref_pose;
	ref_pose.clear();
	io::pdb::build_pose_from_pdb_as_is( ref_pose, *rsd_set, pdb_file_list[8] );
	ideal_poses_[ "south" ][ "  U" ] = ref_pose;
	ref_pose.clear();

	// Get all nonnatural RNA RTs.
	// We don't have phenix ideal coordinates for NCNTs, but we have QM.
	utility::vector1< ResidueTypeCOP > ncnts = ResidueTypeFinder( *rsd_set ).name1('X').base_property( RNA ).get_possible_base_residue_types();
	for ( Size ii = 1; ii <= ncnts.size(); ++ii ) {
		std::stringstream ss;
		ss << "X[" << ncnts[ii]->name() << "]X[" << ncnts[ii]->name() << "]X[" << ncnts[ii]->name() << "]";
		make_pose_from_sequence( ref_pose, ss.str(), rsd_set );
		ideal_ncnt_poses_[ ncnts[ii]->name3() ] = ref_pose;
	}

	initialize_from_options();
}

void ErraserMinimizerMover::initialize_from_options() {
	vary_bond_geometry_ = option[ basic::options::OptionKeys::rna::vary_geometry ].value();
	constrain_phosphate_ = option[ constrain_P ].value();
	ready_set_only_ =     option[ ready_set_only ].value();
	skip_minimize_ =      option[ skip_minimize ].value();
	attempt_pyrimidine_flip_ = option[ attempt_pyrimidine_flip ].value();
	std::copy( option[ fixed_res ].value().begin(), option[ fixed_res ].value().end(), std::inserter( fixed_res_list_, fixed_res_list_.end() ) );
	cutpoint_list_ = option[ full_model::cutpoint_open ].value();

	if ( !ready_set_only_ ) {
		//Setup score function.
		std::string score_weight_file = "stepwise/rna/rna_hires_elec_dens";
		if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
			score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
			TR << "User passed in score:weight option: " << score_weight_file << std::endl;
		}
		scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
		edens_scorefxn_->set_weight( elec_dens_atomwise, scorefxn_->get_weight( elec_dens_atomwise ) );
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
		std::string score_weight_file = "stepwise/rna/rna_hires_elec_dens";
		if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
			score_weight_file = option[ basic::options::OptionKeys::score::weights ]();
			TR << "User passed in score:weight option: " << score_weight_file << std::endl;
		}
		scorefxn_ = ScoreFunctionFactory::create_score_function( score_weight_file );
		edens_scorefxn_->set_weight( elec_dens_atomwise, scorefxn_->get_weight( elec_dens_atomwise ) );
	}
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

Real
ErraserMinimizerMover::ideal_length(
	std::string const & pucker,
	std::string const & name,
	std::string const & name1,
	std::string const & name2
) {
	// We could work with this for inter-residue measures too...
	// OK, here's the deal. We have a map of ideal residues to draw on.
	return ideal_poses_[ pucker ][ name ].residue(2).xyz( name1 ).distance(
		ideal_poses_[ pucker ][ name ].residue(2).xyz( name2 ) );
}

Real
ErraserMinimizerMover::ideal_length_ncnt(
	std::string const & name,
	std::string const & name1,
	std::string const & name2
) {
	// We could work with this for inter-residue measures too...
	// OK, here's the deal. We have a map of ideal residues to draw on.
	return ideal_ncnt_poses_[ name ].residue(2).xyz( name1 ).distance(
		ideal_ncnt_poses_[ name ].residue(2).xyz( name2 ) );
}


Real
ErraserMinimizerMover::ideal_angle(
	std::string const & pucker,
	std::string const & name,
	std::string const & name1,
	Size const no1,
	std::string const & name2,
	Size const no2,
	std::string const & name3,
	Size const no3
) {
	// We could work with this for inter-residue measures too...
	// OK, here's the deal. We have a map of ideal residues to draw on.
	return angle_radians(
		ideal_poses_[ pucker ][ name ].residue(no1).xyz( name1 ),
		ideal_poses_[ pucker ][ name ].residue(no2).xyz( name2 ),
		ideal_poses_[ pucker ][ name ].residue(no3).xyz( name3 )
	);
}

Real
ErraserMinimizerMover::ideal_angle_ncnt(
	std::string const & name,
	std::string const & name1,
	Size const no1,
	std::string const & name2,
	Size const no2,
	std::string const & name3,
	Size const no3
) {
	// We could work with this for inter-residue measures too...
	// OK, here's the deal. We have a map of ideal residues to draw on.
	return angle_radians(
		ideal_ncnt_poses_[ name ].residue(no1).xyz( name1 ),
		ideal_ncnt_poses_[ name ].residue(no2).xyz( name2 ),
		ideal_ncnt_poses_[ name ].residue(no3).xyz( name3 )
	);
}

bool
ErraserMinimizerMover::ideal_has_atom(
	ResidueType const & rt,
	std::string const & an
) {
	// Look in ncnt or normal?
	if ( rt.name1() == 'X' ) {
		return ideal_ncnt_poses_[ rt.name3() ].residue_type( 2 ).has( an );
	} else {
		return ideal_poses_[ "north" ][ rt.name3() ].residue_type( 2 ).has( an );
	}
}

void
ErraserMinimizerMover::add_bond_constraint(
	AtomID const & atom_id1,
	AtomID const & atom_id2,
	core::pose::Pose const & pose,
	core::scoring::constraints::ConstraintSetOP & cst_set
) {
	#define IDEAL_PHOS 1.60708825
	Real const delta_cutoff = 115.0;
	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );

	if ( !ideal_has_atom( pose.residue_type( atom_id1.rsd() ), atom_name1 ) ) return;
	if ( !ideal_has_atom( pose.residue_type( atom_id2.rsd() ), atom_name2 ) ) return;

	// already addressed, just return
	if ( check_in_bonded_list( atom_id1, atom_id2 ) ) return;
	bonded_atom_list_.push_back( std::make_pair( atom_id1, atom_id2 ) );

	Real const bond_length_sd_( 0.05 );
	Real bond_length = 0;//( pose_reference_.residue( atom_id1.rsd() ).xyz( atom_name1 ) -
	//pose_reference_.residue( atom_id2.rsd() ).xyz( atom_name2 ) ).length();
	// either same or different residue
	if ( atom_id1.rsd() == atom_id2.rsd() ) {
		if ( pose.residue_type( atom_id1.rsd() ).name1() == 'X' ) {
			bond_length = ideal_length_ncnt( pose.residue_type( atom_id1.rsd() ).name3(), atom_name1, atom_name2 );
		} else {
			std::string pucker = pose.torsion( TorsionID( atom_id1.rsd(), id::BB, DELTA) ) > delta_cutoff ? "south" : "north";
			bond_length = ideal_length( pucker, pose.residue_type( atom_id1.rsd() ).name3(), atom_name1, atom_name2 );
		}
	} else {
		// polymeric connection for them to be bonded in atom tree--so we are
		// assuming either O5' to LOWER or O3' to UPPER?
		bond_length = IDEAL_PHOS; // ( ideal_length( pucker, pose.residue_type( atom_id1.rsd() ).name()
	}
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
	Real const delta_cutoff = 115.0;
	if ( atom_id2 == atom_id3 ) return;
	if ( atom_id1 == atom_id3 ) return;
	if ( atom_id1 == atom_id2 ) return;

	std::string const & atom_name1 = pose.residue_type( atom_id1.rsd() ).atom_name( atom_id1.atomno() );
	std::string const & atom_name2 = pose.residue_type( atom_id2.rsd() ).atom_name( atom_id2.atomno() );
	std::string const & atom_name3 = pose.residue_type( atom_id3.rsd() ).atom_name( atom_id3.atomno() );

	if ( !ideal_has_atom( pose.residue_type( atom_id1.rsd() ), atom_name1 ) ) return;
	if ( !ideal_has_atom( pose.residue_type( atom_id2.rsd() ), atom_name2 ) ) return;
	if ( !ideal_has_atom( pose.residue_type( atom_id3.rsd() ), atom_name3 ) ) return;

	// if already handled, return
	if ( check_in_bond_angle_list( atom_id1, atom_id2, atom_id3 ) ) return;

	bond_angle_list_.push_back( std::make_pair( atom_id1, std::make_pair( atom_id2, atom_id3 ) ) );
	Real const bond_angle_sd_( radians( 3.0 ) );
	// pass rn - min + 1
	Size smallest_rn = std::min( atom_id1.rsd(), std::min( atom_id2.rsd(), atom_id3.rsd() ) );
	std::string pucker = pose.torsion( TorsionID( atom_id1.rsd(), id::BB, DELTA) ) > delta_cutoff ? "south" : "north";
	Real bond_angle = 0;
	if ( pose.residue_type( atom_id1.rsd() ).name1() == 'X' ) {
		bond_angle =
			ideal_angle_ncnt(
			pose.residue_type( atom_id1.rsd() ).name3(),
			atom_name2,
			atom_id2.rsd() - smallest_rn + 2,
			atom_name1,
			atom_id1.rsd() - smallest_rn + 2,
			atom_name3,
			atom_id3.rsd() - smallest_rn + 2
		);
	} else {
		bond_angle =
			ideal_angle( pucker,
			pose.residue_type( atom_id1.rsd() ).name3(),
			atom_name2,
			atom_id2.rsd() - smallest_rn + 2,
			atom_name1,
			atom_id1.rsd() - smallest_rn + 2,
			atom_name3,
			atom_id3.rsd() - smallest_rn + 2
		);
	}

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
	ObjexxFCL::FArray1D< bool > & allow_insert, // Operationally: not fixed, cutpoint, virt
	std::set< Size > const & chunk
) {
	ConstraintSetOP cst_set = pose.constraint_set()->clone();

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
		if ( !allow_insert( i ) ) continue;

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

	// Next, build lists on bonds and bond angles and restrain them to present values.
	for ( Size i = 1; i <= nres; i++ )  {
		// Don't do anything for protein residues, because we don't have them as ideals.
		// In the future, apply cart_bonded.
		if ( pose.residue_type( i ).is_protein() ) continue;

		if ( chunk.find( i ) == chunk.end() ) continue;
		if ( pose.residue_type( i ).aa() == core::chemical::aa_vrt ) continue; //FCC
		if ( !allow_insert( i ) ) continue;

		chemical::ResidueType const & residue_type( pose.residue_type( i ) );

		for ( Size j = 1; j <= residue_type.natoms(); j++ ) {
			if ( !i_want_this_atom_to_move( residue_type, j ) ) continue;

			AtomID j_atomid( j, i );
			utility::vector1< AtomID > nbrs( pose.conformation().bonded_neighbor_all_res( j_atomid ) );

			// Constrain all mobile bond lengths.
			for ( auto const & nbr : nbrs ) {
				if ( nbr.rsd() > nres || nbr.rsd() < 1 ) continue;
				// Don't do anything for protein residues, because we don't have them as ideals.
				// In the future, apply cart_bonded.
				if ( pose.residue_type( nbr.rsd() ).is_protein() ) continue;

				chemical::ResidueType const & residue_type2( pose.residue_type( nbr.rsd() ) );
				Size const k( nbr.atomno() );

				if ( ! check_if_connected_in_atom_tree( pose, j_atomid, nbr ) ) continue;

				if ( i_want_this_atom_to_move( residue_type2, k ) )  {
					add_bond_constraint( j_atomid, nbr, pose, cst_set );
				}
			}

			// Bond angles
			for ( auto const & nbr : nbrs ) {
				if ( nbr.rsd() > nres || nbr.rsd() < 1 ) continue;
				// Don't do anything for protein residues, because we don't have them as ideals.
				// In the future, apply cart_bonded.
				if ( pose.residue_type( nbr.rsd() ).is_protein() ) continue;

				chemical::ResidueType const & residue_type2( pose.residue_type( nbr.rsd() ) );

				if ( ! check_if_connected_in_atom_tree( pose, j_atomid, nbr ) ) continue;
				Size const k( nbr.atomno() ) ;

				for ( auto const & ang_nbr : nbrs ) {
					// Don't do anything for protein residues, because we don't have them as ideals.
					// In the future, apply cart_bonded.
					if ( pose.residue_type( ang_nbr.rsd() ).is_protein() ) continue;

					if ( ang_nbr.rsd() > nres || ang_nbr.rsd() < 1 ) continue;
					chemical::ResidueType const & residue_type3( pose.residue_type( ang_nbr.rsd() ) );

					if ( !check_if_connected_in_atom_tree( pose, j_atomid, ang_nbr ) ) continue;
					Size const q( ang_nbr.atomno() ) ;

					if ( i_want_this_atom_to_move( residue_type2, k ) &&
							i_want_this_atom_to_move( residue_type3, q ) )  {
						add_bond_angle_constraint( j_atomid, nbr, ang_nbr, pose, cst_set );
					}
				}

				utility::vector1< AtomID > nbrs2( pose.conformation().bonded_neighbor_all_res( nbr ) );

				for ( auto const & nbr2 : nbrs2 ) {
					// Don't do anything for protein residues, because we don't have them as ideals.
					// In the future, apply cart_bonded.
					if ( pose.residue_type( nbr2.rsd() ).is_protein() ) continue;

					if ( nbr2.rsd() > nres || nbr2.rsd() < 1 ) break;
					chemical::ResidueType const & residue_type3( pose.residue_type( nbr2.rsd() ) );

					if ( ! check_if_connected_in_atom_tree( pose, nbr, nbr2 ) ) continue;
					Size const q( nbr2.atomno() ) ;

					if ( i_want_this_atom_to_move( residue_type2, k ) &&
							i_want_this_atom_to_move( residue_type3, q ) )  {
						add_bond_angle_constraint( j_atomid, nbr, nbr2, pose, cst_set );
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
		TR << "Unsliced: " << res_list_unsliced.size() << std::endl;
		TR << "Current:  " << res_list_current.size() << std::endl;
		for ( auto const & sl : sliced_list_final ) {
			TR << "One slice:  " << sl.size() << std::endl;
		}

		if ( res_list_current.size() == 0 ) {
			// python pop 0 to res_list_current
			Size res = *res_list_unsliced.begin();
			TR << "Starting new chunk from scratch with res " << res << std::endl;
			chunk_size = res_list_unsliced.size() / n_chunk_left;
			res_list_unsliced.erase(res_list_unsliced.begin());
			res_list_current.insert( res );
			TR << "Objective is to get " << chunk_size << " res in this chunk." << std::endl;
			n_chunk_left -= 1;
		}

		std::set< Size > res_list_new = find_nearby_res(pose, res_list_current, 3.5 );
		TR << "Adding " << res_list_new.size() << " new residues near res_list_current" << std::endl;
		TR << "To wit, res_list_new: [ ";
		for ( auto const elem : res_list_new ) TR << elem << " ";
		TR << "]" << std::endl;

		// Remove all residues previously in a res list from res_list_new
		clean_res_list( res_list_new, res_list_sliced );

		if ( res_list_new.size() == 0 && res_list_current.size() < chunk_size * 0.7 ) {
			TR << "No nearby new residues identified, but the current res list size " <<  res_list_current.size();
			TR << " is still less than chunk_size * 0.7 " << chunk_size * 0.7 << std::endl;
			while ( true ) {
				// "pop 0th element to res"
				Size res = *res_list_unsliced.begin();
				TR << "So, we pop " << res << std::endl;
				res_list_unsliced.erase(res_list_unsliced.begin());
				if ( res_list_current.find( res ) == res_list_current.end() ) {
					res_list_new.insert(res);
					break;
				}
			}
		}

		TR << "Inserting res_list_new ( " << res_list_new.size() << " residues ) into current " << std::endl;
		res_list_current.insert( res_list_new.begin(), res_list_new.end() );

		TR << "We have obtained a residue list of " << res_list_current.size() << " and our target chunk size is " << chunk_size << std::endl;
		if ( res_list_current.size() >= chunk_size || res_list_new.size() == 0 ) {
			TR << "We have exceeded chunk size with res_list_current or failed at res_list_new" << std::endl;

			TR << "res_list_current: [ ";
			for ( auto const elem : res_list_current ) TR << elem << " ";
			TR << "]" << std::endl;

			std::set< Size > unused;
			fill_gaps_and_remove_isolated_res( res_list_current, total_res, unused );

			TR << "Gonna remove res_list_current from res_list_unsliced" << std::endl;
			for ( auto const & res : res_list_current ) {
				//TR << "Trying to remove " << res << std::endl;
				if ( res_list_unsliced.find(res) != res_list_unsliced.end() ) {
					//TR << "Found" << res << std::endl;
					res_list_unsliced.erase(res_list_unsliced.find(res));
					//TR << "So the unsliced length is now " << res_list_unsliced.size()  << std::endl;
				}
			}

			TR << "Gonna add this chunk to my lists of chunks" << std::endl;
			res_list_sliced.push_back( res_list_current );
			sliced_list_final.push_back( res_list_current );
			if ( n_chunk_left == 0 ) {
				return;
			} else if ( n_chunk_left == 1 ) {
				TR << "What remains is: ";
				for ( Size const res : res_list_unsliced ) TR << res << " ";
				//if ( res_list_unsliced.size() == 0 ) break;
				TR << std::endl;
				res_list_current = res_list_unsliced;
				std::set< Size > removed_res;
				fill_gaps_and_remove_isolated_res(res_list_current, total_res, removed_res);
				sliced_list_final.push_back( res_list_current );
				TR << "Adding res_list_current to the list of sliced_list_final" << std::endl;
				while ( removed_res.size() != 0 ) {
					TR << "Trying to remove res from each sliced list, to some end" << std::endl;
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
				TR << "Ready for new beginnings." << std::endl;
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
	std::string const & working_sequence = pose.sequence();
	TR << "Pose sequence = " << working_sequence << std::endl;

	TR << "Do we have to flip pyrimidines?" << std::endl;
	// Try flipping the pyrimidines
	if ( attempt_pyrimidine_flip_ ) {
		pyrimidine_flip_trial( pose );
	}
	TR << "Back from a possibly pyrimidine flip trial." << std::endl;

	TR << "Do we skip minimization?" << std::endl;
	if ( skip_minimize_ ) return;
	TR << "Apparently not." << std::endl;

	// If we have more than 150 residues, we need to "slice." This involves
	// designing ~100-residue chunks of the pose and making only those DOFs
	// visible to the minimizer. All residues of the pose will be present.
	// Each minimization gets its own little output file. Notably, while we
	// must figure out the chunks every time, we can skip chunks for which
	// an appropriate output file exists.
	utility::vector1< std::set< Size > > chunks;
	// AMW TODO: read chunks from temp file if exists
	TR << "About to try to identify chunks if needed." << std::endl;
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
		cut_lower.insert( ustream );
		cut_upper.insert( dstream );

		if ( pose.residue_type( ustream ).aa() == core::chemical::aa_vrt ) continue;
		if ( pose.residue_type( dstream ).aa() == core::chemical::aa_vrt ) continue;

		if ( fixed_res_list_.size() != 0 &&
				fixed_res_list_.find( ustream ) == fixed_res_list_.end() &&
				fixed_res_list_.find( dstream ) == fixed_res_list_.end() ) {
			mm.set_jump( i, true );
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
		float const dummy_tol( 0.00000001 );

		TR << "Minimize using dfpmin with use_nb_list=true .." << std::endl;
		MinimizerOptions min_options_dfpmin( option[ run::min_type ], dummy_tol, true, false, false );
		min_options_dfpmin.max_iter( std::min( 3000, std::max( 1000, int(nres_moving * 12) ) ) );
		minimizer.run( pose, chunk_mm, *scorefxn_, min_options_dfpmin );

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
		time_t chunk_end = time(0);
		out << "CPU time for chunk: " << (chunk_end-chunk_start) << " seconds." << std::endl;
		out << "DONE!" << std::endl;
		std::stringstream fn;
		fn << "debug_chunk_" << chunk_i << ".pdb";
		pose.pdb_info()->obsolete( false );

		pose.dump_pdb( fn.str() );
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
