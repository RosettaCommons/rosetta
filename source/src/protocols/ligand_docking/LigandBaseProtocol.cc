// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/LigandBaseProtocol.cc
///
/// @brief
/// @author Ian W. Davis


#include <protocols/ligand_docking/LigandBaseProtocol.hh>

#include <core/graph/Graph.hh>

#include <core/chemical/automorphism.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeFinder.hh>
#include <core/chemical/VariantType.hh>

#include <core/conformation/Residue.hh>

#include <core/id/AtomID.hh>
#include <basic/database/open.hh>

#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>

#include <basic/options/option.hh>
//#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/UnboundRotamersOperation.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
//#include <core/scoring/func/Func.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <basic/Tracer.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/rigid/RB_geometry.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray1.io.hh>
#include <numeric/random/random.hh>

#include <algorithm>
#include <cmath>


// option key includes

#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/AtomType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/tree/Atom.hh>
#include <core/scoring/EnergyGraph.hh>


namespace protocols {
namespace ligand_docking {


using namespace ObjexxFCL;


static THREAD_LOCAL basic::Tracer TR( "protocols.ligand_docking.LigandBaseProtocol" );


LigandBaseProtocol::LigandBaseProtocol():
	Mover(),
	sc_interface_padding_(0),
	bb_interface_cutoff_(0)
{
	Mover::type( "LigandBaseProtocol" );

	unboundrot_ = core::pack::rotamer_set::UnboundRotamersOperationOP( new core::pack::rotamer_set::UnboundRotamersOperation() );
	unboundrot_->initialize_from_command_line();

	using namespace basic::options;
	use_soft_rep_ = option[ OptionKeys::docking::ligand::soft_rep ];
	bool const rosetta_electrostatics = option[ OptionKeys::docking::ligand::old_estat ];
	bool const hbonds_downweight = true;

	// Set up scoring function
	// Meiler & Baker 06 uses the soft-repulsive weights, but David wants me to use the hard ones...

	if ( option[ OptionKeys::docking::ligand::tweak_sxfn ] ) {
		TR.Debug << "Using regular tweaked scorefunctions" << std::endl;
		hard_scorefxn_ = make_tweaked_scorefxn("ligand", rosetta_electrostatics, rosetta_electrostatics, hbonds_downweight);
		soft_scorefxn_ = make_tweaked_scorefxn("ligand_soft_rep", rosetta_electrostatics, rosetta_electrostatics, hbonds_downweight);
	} else {
		// Use plain scorefunctions - user specified ones if given, else the regular ones (but non-tweaked).
		// Note that you should now be able to specify exclude_protein_protein_fa_elec in the weights files.
		if ( option[ OptionKeys::score::weights ].user() || option[ OptionKeys::score::patch ].user() ) {
			TR.Debug << "Using untweaked command-line specified (hard) scorefunction. " << std::endl;
			hard_scorefxn_ = core::scoring::get_score_function();
		} else {
			TR.Debug << "Using untweaked ligand.wts hard scorefunction." << std::endl;
			hard_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("ligand.wts");
		}
		if ( option[ OptionKeys::score::soft_wts ].user() ) {
			TR.Debug << "Using untweaked -score:soft_wts specified soft scorefunction." << std::endl;
			soft_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function( option[ OptionKeys::score::soft_wts ] );
		} else {
			TR.Debug << "Using untweaked ligand_soft_rep.wts soft scorefunction." << std::endl;
			soft_scorefxn_ = core::scoring::ScoreFunctionFactory::create_score_function("ligand_soft_rep.wts");
		}
	}

	// "Default" score is always hard;  use soft-rep only inside of apply()
	scorefxn_ = hard_scorefxn_; //( use_soft_rep_ ? soft_scorefxn_ : hard_scorefxn_ );
}


LigandBaseProtocol::~LigandBaseProtocol() = default;

core::scoring::ScoreFunctionOP LigandBaseProtocol::get_scorefxn() {
	return scorefxn_;
}
core::scoring::ScoreFunctionCOP LigandBaseProtocol::get_scorefxn() const {
	return scorefxn_;
}

void LigandBaseProtocol::apply( core::pose::Pose & /*pose*/ ) { utility_exit_with_message("Not intended to actually be used!");}


std::string
LigandBaseProtocol::get_name() const {
	return "LigandBaseProtocol";
}

core::scoring::ScoreFunctionOP
LigandBaseProtocol::make_tweaked_scorefxn(
	std::string const & weights_tag,
	bool estat_exclude_protein,
	bool estat_upweight,
	bool hbonds_downweight
)
{
	using namespace core::scoring;

	ScoreFunctionOP sfxn( new ScoreFunction() );
	sfxn->reset();

	// manipulate EnergyMethodOptions here
	methods::EnergyMethodOptions options( sfxn->energy_method_options() );
	options.exclude_protein_protein_fa_elec( estat_exclude_protein );
	sfxn->set_energy_method_options( options );

	sfxn->add_weights_from_file( basic::database::full_name( "scoring/weights/"+weights_tag+".wts" ) );

	// Tiny weight here (like standard.wts) presumably eliminates the worst intra-ligand clashes...
	// Weight increased because I was still getting significant overlaps between ligand atoms,
	// enough to get new "bonds" in PyMol (though no more interpenetrating rings even at 0.004).
	// However, using 0.04 meant total fa_intra_rep ~ 20 in most structures, which I'm concerned
	// biased the selection of ligand conformers too much.  Needs more rigorous testing.
	if ( sfxn->has_zero_weight( fa_intra_rep ) ) sfxn->set_weight( fa_intra_rep, 0.004 ); // from standard.wts

	// For some reason, electrostatics is not in the .wts files...
	// fa_elec has a different dielectric constant than Rosetta++ (10r vs. 6r in ++)
	// It also includes all atom pairs instead of only ligand-protein interactions.
	if ( sfxn->has_zero_weight( fa_elec ) ) sfxn->set_weight( fa_elec, 0.25 ); // from Meiler & Baker 2006

	if ( estat_upweight ) sfxn->set_weight( fa_elec, (10./6.) * sfxn->get_weight( fa_elec ) ); // make like Rosetta++

	if ( hbonds_downweight ) {
		sfxn->set_weight( hbond_sc, 1.30 ); // from Lin Jiang
		sfxn->set_weight( hbond_bb_sc, 1.30 ); // from Lin Jiang
	}

	// You can adjust later, but at least make sure it's enabled!
	// Shouldn't cause problems when no constraints are present.
	// Need atom_pair_constraint and angle_constraint to work with EnzDes constraints.
	if ( sfxn->has_zero_weight( coordinate_constraint ) ) sfxn->set_weight( coordinate_constraint, 1.0 );
	if ( sfxn->has_zero_weight( atom_pair_constraint ) )  sfxn->set_weight( atom_pair_constraint, 1.0 );
	if ( sfxn->has_zero_weight( angle_constraint ) )      sfxn->set_weight( angle_constraint, 1.0 );
	if ( sfxn->has_zero_weight( dihedral_constraint ) )   sfxn->set_weight( dihedral_constraint, 1.0 );
	if ( sfxn->has_zero_weight( chainbreak ) )            sfxn->set_weight( chainbreak, 1.0 );
	if ( sfxn->has_zero_weight( omega ) )                 sfxn->set_weight( omega, 0.5 );
	if ( sfxn->has_zero_weight( rama ) )                  sfxn->set_weight( rama, 0.2 );

	return sfxn;
}

/// @details First discards ligands that aren't touching, then takes the top 5% by total_score.
/// (Take given number of poses if to_keep > 1.0).
void select_best_poses(
	core::import_pose::atom_tree_diffs::ScoresPairList const & scores_in,
	core::import_pose::atom_tree_diffs::ScoresPairList & scores_out,
	core::Real to_keep /* = 0.05 */
)
{
	// Keep only the top 5% by total score
	//scores_out.reserve( scores_in.size() );
	//std::cout << "scores_in.size() = " << scores_in.size() << "\n";
	//std::cout << "scores_out.size() = " << scores_out.size() << " (0)\n";
	for ( core::Size ii = 1; ii <= scores_in.size(); ++ii ) {
		// Drop out cases where the ligand isn't touching the protein
		core::import_pose::atom_tree_diffs::Scores scores = scores_in[ii].second;
		if ( scores.find("ligand_is_touching") == scores.end()
				||  scores["ligand_is_touching"] != 0 ) {
			scores_out.push_back( scores_in[ii] );
		}
	}
	//std::cout << "scores_out.size() = " << scores_out.size() << " (1)\n";
	core::import_pose::atom_tree_diffs::AtomTreeDiff::sort_by("total_score", scores_out);
	//std::cout << "scores_out.size() = " << scores_out.size() << " (2)\n";
	if ( to_keep <= 1.0  && to_keep >= 0.0 ) {
		scores_out.resize( (core::Size) std::ceil(to_keep * scores_out.size()) );
	} else if ( to_keep > 1.0 && scores_out.size() > to_keep ) {
		scores_out.resize( (core::Size) std::ceil(to_keep) );
	} else {
		utility_exit_with_message("Cannot select a negative quantity of poses.");
	}
	//std::cout << "scores_out.size() = " << scores_out.size() << " (3)\n";
	// Although not everyone will care, this is useful for some applications.
	core::import_pose::atom_tree_diffs::AtomTreeDiff::sort_by("interface_delta", scores_out);
	//for(core::Size jj = 1; jj <= 3; ++jj)
	// std::cout << "Best structure: " << scores_out[jj].first << '\n';
}

void select_best_poses(
	core::import_pose::atom_tree_diffs::AtomTreeDiff const & atdiff,
	core::import_pose::atom_tree_diffs::ScoresPairList & scores_out,
	core::Real to_keep /* = 0.05*/
)
{
	core::import_pose::atom_tree_diffs::ScoresPairList const & scores_list = atdiff.scores();
	select_best_poses(scores_list, scores_out, to_keep);
}

void select_best_poses(
	core::import_pose::atom_tree_diffs::AtomTreeDiff const & atdiff,
	std::set< std::string > & tags_out
)
{
	core::import_pose::atom_tree_diffs::ScoresPairList scores_list2;
	scores_list2.reserve( atdiff.scores().size() );
	select_best_poses(atdiff, scores_list2);
	for ( core::Size ii = 1; ii <= scores_list2.size(); ++ii ) {
		tags_out.insert( scores_list2[ii].first );
	}
}

/// @details Currently considers only heavy atoms, not hydrogens.
void frac_atoms_within(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< core::Real > const & cutoffs,
	utility::vector1< core::Real > & fractions_out
)
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::conformation;
	// name() and total number of atoms may actually be different, if we're comparing e.g. tautomers
	if ( rsd1.type().name3() != rsd2.type().name3() ) utility_exit_with_message("Residue type name3 mismatch");
	if ( rsd1.nheavyatoms()  != rsd2.nheavyatoms()  ) utility_exit_with_message("Residue number-of-heavy-atoms mismatch");
	fractions_out.resize(cutoffs.size(), 0);
	int counter = 0;
	// Make atom-number translation table
	AutomorphismIterator ai( rsd1.type() );
	AtomIndices old2new( ai.next() );
	// For each permutation of automorphisms...
	while ( old2new.size() > 0 ) {
		counter++;
		//if( counter%10000 == 0 ) tr.Info << counter << " so far..." << std::endl;

		// Print out translation table for debugging
		//std::cout << "[";
		//for(Size i = 1; i <= old2new.size(); ++i) std::cout << " " << old2new[i];
		//std::cout << " ]\n";
		//for(Size j = 1; j <= old2new.size(); ++j) std::cout << "  " << j << " --> " << old2new[j] << "  /  " << rsd1.type().atom_name(j) << " --> " << rsd1.type().atom_name(old2new[j]) << "\n";

		// Each cutoff might find its maximum fraction from a different automorphism (?)
		// I'm not sure what this means, if anything, for the validity of this measure...
		for ( core::Size i = 1; i <= cutoffs.size(); ++i ) {
			core::Real const cutoff2 = cutoffs[i] * cutoffs[i];
			core::Real nwithin( 0 );
			core::Real natoms( 0 );
			for ( core::Size j = 1; j <= rsd1.type().natoms(); ++j ) {
				if ( !rsd1.atom_type(j).is_hydrogen() ) {
					natoms += 1;
					// This is the step where we effectively re-assign atom names
					// in hopes of reducing RMS (e.g. by "flipping" a phenyl ring).
					core::Vector diff = rsd1.xyz( j ) - rsd2.xyz( old2new[j] );
					if ( diff.length_squared() <= cutoff2 ) nwithin += 1;
				}
			}
			core::Real fracwithin = (natoms > 0 ? nwithin / natoms : 0);
			if ( fracwithin > fractions_out[i] ) {
				//tr.Debug << "New rms of " << curr_rms << " beats previous best of " << best_rms << std::endl;
				fractions_out[i] = fracwithin;
			}
		}
		old2new = ai.next();
	} // done checking all automorphisms
	//tr.Info << counter << " automorphisms from iterator; best rms is " << best_rms << std::endl;
}

core::Size
LigandBaseProtocol::get_ligand_jump_id(
	core::pose::Pose const & pose
)const{
	int jump_id = pose.num_jump(); // assume ligand attached by last jump
	if ( jump_id == 0 ) {
		utility_exit_with_message("Pose has no jumps!");
	}
	return jump_id;
}

/// @brief Return the residue sequence number for our ligand.
/// @details Only works with single residue ligands, really.
/// Reconsider this in the future.
core::Size
LigandBaseProtocol::get_ligand_id(
	core::pose::Pose const & pose
) const
{
	return get_ligand_id(pose, get_ligand_jump_id(pose));
}

/// @brief Return the residue sequence number for our ligand.
/// @details Only works with single residue ligands, really.
/// Reconsider this in the future.
core::Size
LigandBaseProtocol::get_ligand_id(
	core::pose::Pose const & pose,
	core::Size jump_id
) const
{
	core::Size const lig_id = (core::Size) pose.fold_tree().downstream_jump_residue(jump_id);

	// Safety checks...
	FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );
	core::Size num_downstream = 0;
	for ( core::Size i = 1; i <= pose.size(); ++i ) if ( !is_upstream(i) ) num_downstream += 1;
	if ( lig_id != pose.size() || num_downstream != 1 || is_upstream(lig_id) || pose.residue(lig_id).is_polymer() ) {
		utility_exit_with_message("Expected ligand to be last residue in pose and only one downstream of the jump");
	}

	return lig_id;
}

core::Vector LigandBaseProtocol::choose_desired_centroid(
	core::pose::Pose const & pose,
	core::Size jump_id,
	utility::vector1< core::Vector >  start_from_pts
){
	// Choose desired centroid:  either -start_from or the current position.
	core::Vector desired_centroid;
	if ( !start_from_pts.empty() ) {
		int const which_triple = numeric::random::rg().random_range(1, start_from_pts.size());
		desired_centroid = start_from_pts[ which_triple ];
	} else {
		core::Vector dummy;
		protocols::geometry::centroids_by_jump(pose, jump_id, dummy, desired_centroid);
	}
	//std::cout << "Desired centroid: " << desired_centroid << std::endl;
	return desired_centroid;
}
void LigandBaseProtocol::move_ligand_to_desired_centroid(
	core::pose::Pose & pose,
	core::Size jump_id,
	utility::vector1< core::Vector >  start_from_pts
){
	core::Vector desired_centroid= choose_desired_centroid(pose, jump_id, start_from_pts);
	move_ligand_to_desired_centroid(pose, jump_id, desired_centroid);
}

void LigandBaseProtocol::move_ligand_to_desired_centroid(
	core::pose::Pose & pose,
	core::Size jump_id,
	core::Vector desired_centroid
){//(This could be a no-op if desired == current)
	core::Vector ligand_centroid = protocols::geometry::downstream_centroid_by_jump(pose, jump_id);
	//std::cout << "Centroid before: " << ligand_centroid << std::endl;
	core::Vector const trans_vec = desired_centroid - ligand_centroid;
	core::Real const trans_len = trans_vec.length();
	if ( trans_len > 1e-3 ) { // otherwise we get NaNs
		protocols::rigid::RigidBodyTransMover mover( pose, jump_id);
		mover.step_size(trans_len);
		mover.trans_axis(trans_vec);
		mover.apply(pose);
	}
	//std::cout << "Centroid after: " << protocols::geometry::downstream_centroid_by_jump(pose, jump_id) << std::endl;
}

core::kinematics::MoveMapOP
LigandBaseProtocol::make_movemap(
	core::pose::Pose const & pose,
	core::Size jump_id,
	core::Real sc_padding,
	bool include_all_rsds,
	bool include_backbone,
	bool include_ligands,
	bool include_water
) const
{
	// All DOF start false (frozen)
	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap() );
	movemap->set_jump(jump_id, true);
	//if( include_backbone ) movemap->set_bb(true); // held in check by restraints (elsewhere)

	FArray1D_bool allow_min( pose.size(), true );
	utility::vector1< bool > dont_care( pose.size(), false );
	utility::vector1< bool > allow_min_bb( pose.size(), false );
	if ( !include_all_rsds ) {
		find_interface_rsds(pose, jump_id, sc_padding, allow_min);
		if ( include_backbone ) {
			find_interface_backbone(pose, jump_id, bb_interface_cutoff_, dont_care, allow_min_bb);
		}
	}
	if ( !include_ligands ) {
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( !pose.residue(i).is_polymer() ) allow_min(i) = false;
		}
	}
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( allow_min(i) ) {
			//std::cout << "Allow residue to minimize: " << i << "\n";
			movemap->set_chi(i, true);
			if ( include_backbone && allow_min_bb[i] ) movemap->set_bb(i, true);
		}
	}
	if ( include_water ) {
		for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
			if ( ! pose.residue(i).has_property("WATER") ) continue;
			core::kinematics::Edge const & e = pose.fold_tree().get_residue_edge(i);
			if ( ! e.is_jump() ) continue;
			movemap->set_jump( e.label(), true );
			TR << "Minimize water jump " << e.label() << " to residue " << i << " " << pose.residue_type(i).name3() << std::endl;
		}
	}

	return movemap;
}

/// @details If ligand_protonation is true, ligand will be allowed to "mutate"
/// to any residue type with the same name3 (3-letter PDB name).
/// If these residues are all protonation / tautomer states, then effectively
/// it is protonation / tautomer states that will be sampled.
/// Some care is required to ensure superposition works -- the nbr_atom must have
/// 2+ heavy atom neighbors, and heavy atoms must be named consistently.
core::pack::task::PackerTaskOP
LigandBaseProtocol::make_packer_task(
	core::pose::Pose const & pose,
	FArray1D_bool const & allow_repack,
	bool ligand_protonation
) const
{
	using namespace core::pack::task;
	PackerTaskOP pack_task = TaskFactory::create_packer_task(pose);
	pack_task->initialize_from_command_line(); // -ex1 -ex2  etc.
	pack_task->append_rotamerset_operation( unboundrot_ );
	//pack_task->restrict_to_repacking(); // all residues -- now set individually below

	// Can be accomplished from the command line with -packing:use_input_sc
	// Meiler and Baker points out that in some cases you want this off for honesty's sake.
	//pack_task->or_include_current(true);

	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		core::conformation::Residue const & this_rsd = pose.residue(i);
		if ( !this_rsd.is_polymer() && ligand_protonation ) {
			using namespace core::chemical;
			ResidueTypeSetCOP rsd_type_set = this_rsd.residue_type_set();
			ResidueTypeCOPs allowed_types = core::chemical::ResidueTypeFinder( *rsd_type_set ).name3( this_rsd.name3() ).get_all_possible_residue_types(); // a vector1
			for ( core::Size j = 1; j <= allowed_types.size(); ++j ) {
				if ( allowed_types[j]->name() == this_rsd.name() ) continue; // already in the task's list
				pack_task->nonconst_residue_task( i ).allow_noncanonical_aa( allowed_types[j]->name() );
			}
			TR << "Allowed residues at position " << i << ":" << std::endl;
			for ( auto rt = pack_task->nonconst_residue_task( i ).allowed_residue_types_begin();
					rt != pack_task->nonconst_residue_task( i ).allowed_residue_types_end(); ++rt ) {
				TR << "  " << (*rt)->name() << std::endl;
			}
		} else {
			pack_task->nonconst_residue_task( i ).restrict_to_repacking();
		}
		if ( !allow_repack(i) /*|| !pose.residue(i).is_polymer()*/ ) {
			pack_task->nonconst_residue_task( i ).prevent_repacking();
		}
	}

	return pack_task;
}

core::pack::task::PackerTaskOP
LigandBaseProtocol::make_packer_task(
	core::pose::Pose const & pose,
	int jump_id,
	core::Real sc_padding,
	bool include_all_rsds,
	bool ligand_protonation
) const
{
	FArray1D_bool allow_repack( pose.size(), true );
	// Disable packing for residues that are too far from the ligand.
	if ( !include_all_rsds ) find_interface_rsds(pose, jump_id, sc_padding, allow_repack);
	return make_packer_task(pose, allow_repack, ligand_protonation);
}

core::pack::task::PackerTaskOP
LigandBaseProtocol::make_packer_task_ligand_only(
	core::pose::Pose const & pose,
	int /*jump_id*/,
	bool ligand_protonation
) const
{
	FArray1D_bool allow_repack( pose.size(), false );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( !pose.residue(i).is_polymer() ) allow_repack(i) = true;
	}
	return make_packer_task(pose, allow_repack, ligand_protonation);
}

void
LigandBaseProtocol::find_interface_rsds(
	core::pose::Pose const & pose,
	int jump_id,
	core::Real padding,
	FArray1D_bool & is_interface //< output
) const
{
	// The Rosetta++ criterion for which sidechains repack/minimize in docking:
	//  ligand heavy atom within paircutoff(aa,GLY)+1 of aa's CB
	// See
	//  docking_minimize.cc   docking_MCM_pack_side_chains()
	//  docking_movement.cc   docking_repack()
	//  docking_scoring.cc    docking_interf_residues()
	//  ligand.cc             detect_ligand_interface[_res]()
	//                        hetero_atom_amino_acid_distance()
	// Ian's criterion to approximate this in Mini:
	//  ligand heavy atom within rsd.nbr_radius()+6 of rsd.nbr_atom()
	// 6A is an eyeballed magic number to get ~ agreement w/ Rosetta++ paircutoffs+1

	int num_in_interface = 0;
	is_interface.dimension( pose.size(), false ); // init all positions to false
	FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		// all residues on ligand side can move
		if ( ! is_upstream(i) ) {
			is_interface(i) = true;
			num_in_interface += 1;
			continue;
		}
		// on protein side, have to do distance check
		core::conformation::Residue const & prot_rsd = pose.residue(i);
		for ( core::Size j = 1, j_end = pose.size(); j <= j_end; ++j ) {
			if ( is_upstream(j) ) continue; // compare against only ligand residues
			core::conformation::Residue const & lig_rsd = pose.residue(j);
			for ( core::Size k = 1, k_end = lig_rsd.nheavyatoms(); k <= k_end; ++k ) {
				double dist2 = lig_rsd.xyz(k).distance_squared( prot_rsd.xyz(prot_rsd.nbr_atom()) );
				double cutoff = prot_rsd.nbr_radius() + 6.0 + padding;
				if ( dist2 <= cutoff * cutoff ) {
					is_interface(i) = true;
					num_in_interface += 1;
					goto END_LIGRES_LOOP; // C++ lacks multi-level break  :(
				}
			}
		}
		END_LIGRES_LOOP: ; // compiler needs ; as a no-op before end of loop
	}
	TR << "Interface is " << num_in_interface << " / " << pose.size()
		<< " residues (" << 100*(double(num_in_interface) / double(pose.size())) << "%)" << std::endl;
}

/// Find residues that would most benefit docking by backbone movement.
/// Based on absolute distance to CA/CB -- being near the end of a long sc
/// is no reason its backbone should move.
/// This is a new invention; not taken from Rosetta++.
void
LigandBaseProtocol::find_interface_backbone(
	core::pose::Pose const & pose,
	int jump_id,
	core::Real cutoff_dist,
	utility::vector1< bool > & is_interface, //< output
	utility::vector1< bool > & is_around_interface //< output
) const
{
	runtime_assert( cutoff_dist > 0 );
	double const cutoff2 = cutoff_dist * cutoff_dist;

	int num_in_interface = 0;
	is_interface.resize( pose.size(), false ); // init all positions to false
	is_around_interface.resize( pose.size(), false ); // init all positions to false

	FArray1D_bool is_upstream ( pose.size(), false );
	pose.fold_tree().partition_by_jump( jump_id, is_upstream );
	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		// all residues on ligand side can move
		if ( ! is_upstream(i) ) {
			is_interface[i] = true;
			num_in_interface += 1;
			continue;
		}
		// on protein side, have to do distance check
		core::conformation::Residue const & prot_rsd = pose.residue(i);
		if ( !prot_rsd.is_protein() ) continue; // only minimize protein backbone, not other stuff.  Maybe also DNA/RNA in the future?
		core::Vector prot_cb;
		if ( prot_rsd.has("CB") ) prot_cb = prot_rsd.xyz("CB");
		else if ( prot_rsd.has("CA") ) prot_cb = prot_rsd.xyz("CA"); // GLY
		else {
			TR << "Can't find CA/CB for residue " << i << std::endl;
			continue; // non-protein residues not part of the interface
		}
		for ( core::Size j = 1, j_end = pose.size(); j <= j_end; ++j ) {
			if ( is_upstream(j) ) continue; // compare against only ligand residues
			core::conformation::Residue const & lig_rsd = pose.residue(j);
			for ( core::Size k = 1, k_end = lig_rsd.nheavyatoms(); k <= k_end; ++k ) {
				double dist2 = lig_rsd.xyz(k).distance_squared( prot_cb );
				if ( dist2 <= cutoff2 ) {
					is_interface[i] = true;
					//std::cout << "{bb iface " << prot_rsd.name3() << " " << i << "} "
					// << prot_cb.x() << " " << prot_cb.y() << " " << prot_cb.z() << "\n";
					num_in_interface += 1;
					goto END_LIGRES_LOOP; // C++ lacks multi-level break  :(
				}
			}
		}
		END_LIGRES_LOOP: ; // compiler needs ; as a no-op before end of loop
	}
	TR << "Backbone interface is " << num_in_interface << " / " << pose.size()
		<< " residues (" << 100*(double(num_in_interface) / double(pose.size())) << "%)" << std::endl;

	int const window = 3; // how many residues on either side of the truly mobile ones?
	for ( Size i = 1, nres = pose.size(); i <= nres; ++i ) {
		if ( pose.residue_type(i).is_polymer() ) {
			// Track backwards
			for ( Size j = i; j >= Size(std::max(1,int(i)-window)); --j ) {
				if ( !pose.residue_type(j).is_polymer() || pose.residue(j).is_upper_terminus() ) break;
				is_around_interface[i] = is_around_interface[i] | is_interface[j];
				if ( pose.residue(j).is_lower_terminus() ) break;
			}
			// Track forwards
			for ( Size j = i; j <= std::min(nres,i+window); ++j ) {
				if ( !pose.residue_type(j).is_polymer() || pose.residue(j).is_lower_terminus() ) break;
				is_around_interface[i] = is_around_interface[i] | is_interface[j];
				if ( pose.residue(j).is_upper_terminus() ) break;
			}
		}
		//std::cout << i << " unrestr " << is_interface[i] << " mobile " << is_around_interface[i] << std::endl;
	}
}

/// @details Depends on atomtree topology;  atomtree setup should be finished before calling this.
void
LigandBaseProtocol::restrain_protein_Calphas(
	core::pose::Pose & pose,
	utility::vector1< bool > const & is_restrained,
	//core::Real stddev_Angstroms,
	core::scoring::func::FuncOP restr_func
) const
{
	using namespace core::scoring::constraints;
	using core::chemical::ResidueType;
	using core::chemical::AtomIndices;
	using core::conformation::Residue;
	using core::id::AtomID;

	// Can use same function for all b/c target is always 0 distance to orig coords.
	//FuncOP restr_func = new HarmonicFunc(0, stddev_Angstroms);
	// An atom that should never move in terms of absolute coordinates.
	// Needed as a proxy for the origin, b/c Rosetta assumes all energies are
	// translation-invariant.  So it's a "two-body" energy with this fixed atom.
	AtomID fixed_pt( pose.atom_tree().root()->atom_id() );

	for ( core::Size i = 1, i_end = pose.size(); i <= i_end; ++i ) {
		if ( pose.residue(i).is_protein() && is_restrained[i] ) { // protein residues
			Residue const & rsd = pose.residue(i);
			ResidueType const & rsd_type = pose.residue_type(i);
			ConstraintOP constraint( new CoordinateConstraint(
				AtomID(rsd_type.atom_index("CA"), i),
				fixed_pt,
				rsd.xyz("CA"),
				restr_func
				) );
			TR.Debug << "Restraining C-alpha of residue " << i << std::endl;
			pose.add_constraint( constraint );
		}
	}

	// You can adjust later, but at least make sure it's enabled!
	if ( scorefxn_->has_zero_weight(core::scoring::coordinate_constraint) ) {
		scorefxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	}
}

/// @details Depends on atomtree topology;  atomtree setup should be finished before calling this.
/// Returns the newly created Constraint for future reference, but also adds it to the Pose.
core::scoring::constraints::ConstraintOP
LigandBaseProtocol::restrain_ligand_nbr_atom(
	core::pose::Pose & pose,
	core::Size lig_id,
	core::Real stddev_Angstroms
) const
{
	using namespace core::scoring::constraints;
	using core::chemical::ResidueType;
	using core::chemical::AtomIndices;
	using core::conformation::Residue;
	using core::id::AtomID;

	core::scoring::func::FuncOP restr_func( new core::scoring::func::HarmonicFunc(0, stddev_Angstroms) );
	// An atom that should never move in terms of absolute coordinates.
	// Needed as a proxy for the origin, b/c Rosetta assumes all energies are
	// translation-invariant.  So it's a "two-body" energy with this fixed atom.
	AtomID fixed_pt( pose.atom_tree().root()->atom_id() );

	Residue const & rsd = pose.residue(lig_id);
	ConstraintOP constraint( new CoordinateConstraint(
		AtomID(rsd.nbr_atom(), lig_id),
		fixed_pt,
		rsd.nbr_atom_xyz(),
		restr_func
		) );
	TR << "Restraining ligand residue " << lig_id << std::endl;
	pose.add_constraint( constraint );

	// You can adjust later, but at least make sure it's enabled!
	if ( scorefxn_->has_zero_weight(core::scoring::coordinate_constraint) ) {
		scorefxn_->set_weight( core::scoring::coordinate_constraint, 1.0 );
	}
	return constraint;
}

/// @details This function actually alters both the pose, and the jump_id to be docked.
/// The pose will have chainbreak variant types and distance constraints added.
void
LigandBaseProtocol::setup_bbmin_foldtree(
	core::pose::Pose & pose,
	core::Size const & jump_id,
	core::Real cutoff_dist,
	core::Real stddev_Angstroms
)
{
	core::Size const nres = pose.size();
	core::Size const lig_id = get_ligand_id(pose, jump_id);

	// Residues whose backbone can move freely
	// This includes the ligand residue for some reason... why did I do that?
	utility::vector1< bool > unrestr_bb( nres, false );
	// Residues whose backbone can move, but is held in place by restraints
	utility::vector1< bool > mobile_bb( nres, false );

	find_interface_backbone(pose, jump_id, cutoff_dist, unrestr_bb, mobile_bb);

	//TR<< "unrestr_bb"<< unrestr_bb << std::endl;
	//TR<< "mobile_bb"<< mobile_bb << std::endl;

	reorder_foldtree_around_mobile_regions( pose, jump_id, mobile_bb, lig_id);

	// Set allow_move_bb to only restrain mobile residues
	utility::vector1< bool >  allow_move_bb( nres, true );
	for ( core::Size i = 1; i <= nres; ++i ) {
		//allow_move_bb[i] = mobile_bb[i] & !unrestr_bb[i];
		allow_move_bb[i] = mobile_bb[i];
	}
	// Add constraints
	core::scoring::func::FuncOP restr_func( new core::scoring::func::HarmonicFunc(0, stddev_Angstroms) );
	restrain_protein_Calphas(pose, allow_move_bb, restr_func);

}//setup_bbmin_foldtree


/// @brief reorders a fold tree such that movement in the mobile regions will
/// @brief have zero effect on the non-mobile regions. for every contiguous string
/// @brief of mobile regions, the downstream non-mobile stretch will be connected to the
/// @brief upstream non-mobile stretch through a jump. further, a cutpoint will be
/// @brief introduced into every mobile stretch. the ligand will be connected to the
/// @brief nearest non-mobile residue.
/// @brief Note: every contiguous stretch in the mobile regions must be at least 4 residues long
void
LigandBaseProtocol::reorder_foldtree_around_mobile_regions(
	core::pose::Pose & pose,
	core::Size const & jump_id,
	utility::vector1< bool > const & mobile_bb,
	core::Size const & lig_id
) const
{

	using namespace core::kinematics;
	FoldTree f = pose.fold_tree(); // a copy
	core::Size const nres = pose.size();

	//sanity check
	if ( mobile_bb.size() != nres ) {
		std::cerr << "Error: vector containing mobile residue information doesn't have the same number of residues as the pose." << std::endl;
		utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// Attach ligand to nearest CA that is not mobile!
	// This used to be activated by the -shorten_jump flag, but now is the default
	if ( lig_id != 0 ) { //flo sept '12this function now can also  be used to modify fold trees for poses that contain new ligands
		Size attach_pt = 1; // for now, attach ligand to residue 1
		core::Real shortest_dist2 = 1e99;
		core::Vector lig_nbr_atom = pose.residue(lig_id).nbr_atom_xyz();
		for ( core::Size i = 1; i <= nres; ++i ) {
			if ( !pose.residue(i).is_polymer() ) continue;
			if ( i == lig_id ) continue; // should be unneccessary
			if ( !pose.residue(i).has("CA") ) continue;
			if ( mobile_bb[i] ) continue;
			core::Real const new_dist2 = lig_nbr_atom.distance_squared( pose.residue(i).xyz("CA") );
			if ( new_dist2 < shortest_dist2 ) {
				shortest_dist2 = new_dist2;
				attach_pt = i;
			}
		}
		Size old_nres = f.nres(), old_njump = f.num_jump();
		// Deleting edges is a dicey proposition ... better to add them to a new foldtree
		FoldTree f_new;
		// Without the stupid cast-to-const, GCC tries to use the private, non-const version.
		FoldTree const & f_const = f;
		 for ( auto const & e : f_const ) {
			bool contains_attach_pt = (( e.start() < int(attach_pt) && int(attach_pt) < e.stop() )
				|| ( e.stop() < int(attach_pt) && int(attach_pt) < e.start() ));
			if ( e.label() == (int) jump_id ) {
				f_new.add_edge( attach_pt, lig_id, jump_id );
			} else if ( e.label() == Edge::PEPTIDE && contains_attach_pt ) {
				f_new.add_edge( e.start(), attach_pt, Edge::PEPTIDE );
				f_new.add_edge( attach_pt,  e.stop(), Edge::PEPTIDE );
			} else {
				f_new.add_edge( e.start(), e.stop(), e.label() );
			}
		}
		f = f_new;
		TR << "Moved ligand " << f << std::endl;
		if ( !f.connected() ) utility_exit_with_message("Fold tree not connected?!");
		if ( !f.check_fold_tree() ) utility_exit_with_message("Fold tree did not pass check!");
		if ( old_njump != f.num_jump() ) utility_exit_with_message("Number of jumps changed?!");
		if ( old_nres != f.nres() ) utility_exit_with_message("Number of residues changed?!");
		if ( lig_id != (Size) f.downstream_jump_residue(jump_id) ) utility_exit_with_message("Ligand no longer at end of last jump!");
	}

	{
		// Deleting edges is a dicey proposition ... better to add them to a new foldtree
		//for(Size k = 1; k <= nres; ++k) std::cerr << k << ' ' << mobile_bb[k] << std::endl;
		FoldTree f_new;
		int new_jump = f.num_jump();
		// Without the stupid cast-to-const, GCC tries to use the private, non-const version.
		FoldTree const & f_const = f;
		 for ( auto const & edge_itr : f_const ) {
			Size const e_start = edge_itr.start();
			Size const e_stop = edge_itr.stop();
			if ( e_stop < e_start ) utility_exit_with_message("Not prepared to deal with backwards fold tree edges!");
			//std::cout << "Considering edge from " << e_start << " to " << e_stop << ", label " << edge_itr->label() << std::endl;
			if ( !edge_itr.is_polymer() ) {
				f_new.add_edge( edge_itr );
				continue; // only subdivide "peptide" edges of the fold tree
			}
			// else:
			using std::max; using std::min;
			using namespace protocols::loops;
			utility::vector1< Loop > loops;
			for ( Size i = e_start; i <= e_stop; ++i ) {
				while ( i <= e_stop && !mobile_bb[i] ) ++i; // find the start of a mobile region
				if ( i > e_stop ) break; // no more mobile regions
				Size const start = i; // first mobile residue
				//std::cout << "  Start from " << start << std::endl;
				while ( i <= e_stop && mobile_bb[i] ) ++i; // find the end of the mobile region
				if ( i > e_stop ) i = e_stop;
				if ( !mobile_bb[i] ) --i; // back up one unless the end is flexible
				Size const stop = i; // last mobile residue
				//std::cout << "  Stop at " << stop << std::endl;

				if ( (stop-start+1) < 4 ) {
					// This kind of thing can come up at chain terminii, and should not cause a fatal error.
					//for(Size k = 1; k <= nres; ++k) std::cerr << k << ' ' << mobile_bb[k] << std::endl;
					TR.Warning << "WARNING: for backbone minimization to work properly, a stretch of at least 4 residues needs to be allowed to move. Stretch between " << start << " and " << stop << " is too short." << std::endl;
					//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
					// But it will probably break the cutpoint logic below, so we have to skip this group of residues.
					continue;
				}

				// Cutpoint should fall between start and stop, but not if the rsd on either side is a terminus.
				// Because this happens within one peptide edge, no "internal" residue should ever be a terminus.
				Size const cut_start = ( pose.residue(start).is_terminus() ? start+1 : start );
				Size const cut_end = ( pose.residue(stop).is_terminus() ? stop-2 : stop-1 );
				runtime_assert( cut_start <= cut_end );
				Size cutpt = Size( numeric::random::rg().random_range(cut_start, cut_end) ); // cut is made between cutpt and cutpt+1
				// Can't use this function while iterating -- invalidates the iterators!
				//f.new_jump( max(e_start,start-1), min(e_stop,stop+1), cutpt );
				//loops.push_back( Loop( max(e_start,start-1), min(e_stop,stop+1), cutpt ) );
				loops.push_back( Loop( start, stop, cutpt ) );
				// also need to set up residue variants so chainbreak score works correctly!
				core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, cutpt );
				core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, cutpt+1 );
			}
			int last_rigid = e_start;
			for ( Size i = 1; i <= loops.size(); ++i ) {
				Loop const & l = loops[i];
				int first = max(int(e_start), int(l.start()-1) );
				int last = min(int(e_stop), int(l.stop()+1) );
				if ( first != last_rigid ) {
					f_new.add_edge( last_rigid, first, Edge::PEPTIDE );
				}
				f_new.add_edge( first, l.cut(), Edge::PEPTIDE );
				//f_new.add_edge( first, last, ++new_jump );
				// Using CA instead of N may help keep downstream res from moving as much?
				f_new.add_edge( Edge( first, last, ++new_jump, "CA", "CA", false /* default val but req'd by compiler */ ) );
				f_new.add_edge( last, l.cut()+1, Edge::PEPTIDE );
				last_rigid = last;
			}
			if ( last_rigid != int(e_stop) ) {
				f_new.add_edge( last_rigid, e_stop, Edge::PEPTIDE );
			}
		}
		f = f_new;
	}

	f.delete_extra_vertices();
	// Find the first non-mobile residue to be the root
	for ( Size i = 1; i <= nres; ++i ) {
		if ( pose.residue(i).is_polymer() && !mobile_bb[i] ) {
			f.reorder( i );
			break;
		}
	}

	TR << "Final loops foldtree " << f << std::endl;
	if ( !f.check_fold_tree() ) {
		utility_exit_with_message("Invalid fold tree after trying to set up for minimization!");
	}

	// The new jump is the assigned the next number in sequence!
	//jump_id = foldtree.num_jump();
	pose.fold_tree( f );

	//core::kinematics::dump_pose_kinemage("final_foldtree.kin", pose);
	if ( lig_id != 0 ) runtime_assert( lig_id == get_ligand_id(pose, jump_id) ); // runtime_assert ligand ID hasn't changed


} // reorder_foldtree_around_mobile_regions


void
LigandBaseProtocol::get_non_bb_clashing_rotamers(
	core::pose::Pose const & pose,
	core::Size seqpos,
	core::scoring::ScoreFunctionCOP scofx,
	utility::vector1< core::conformation::ResidueCOP > & accepted_rotamers
) const {
	using namespace core;

	utility::vector1< conformation::ResidueOP > suggested_rotamers;

	Real bb_bump_cutoff = basic::options::option[ basic::options::OptionKeys::enzdes::bb_bump_cutoff ];

	pack::task::PackerTaskOP help_task = pack::task::TaskFactory::create_packer_task( pose );
	for ( core::Size i = 1; i <= pose.size(); ++i ) {
		if ( i == seqpos ) help_task->nonconst_residue_task( i ).restrict_to_repacking();
		else help_task->nonconst_residue_task( i ).prevent_repacking();
	}

	//let's see if this works
	//graph::GraphOP neighbor_graph = pack::create_packer_graph( pose, *scofx, help_task );
	graph::GraphOP neighbor_graph( new graph::Graph( pose.energies().energy_graph() ) );

	chemical::ResidueTypeCOP res_type = pose.residue_type( seqpos ).get_self_ptr();
	conformation::Residue const & existing_residue( pose.residue( seqpos ) );

	utility::vector1< utility::vector1< Real > > extra_chi_steps( res_type->nchi() );

	pack::rotamers::SingleResidueRotamerLibraryCOP rotlib( core::pack::rotamers::SingleResidueRotamerLibraryFactory::get_instance()->get( *res_type ) );

	if ( rotlib ) {
		rotlib->fill_rotamer_vector( pose, *scofx, *help_task, neighbor_graph, res_type, existing_residue, extra_chi_steps, true /*buried*/, suggested_rotamers );
	} else return;


	//now that we have the suggested rotamers, let's do the bump_check
	utility::vector1< core::conformation::ResidueOP > temp_accepted_rotamers; //buffer to allow sorting
	core::conformation::ResidueOP best_rot;
	core::Real bestE(1000000.0);

	for ( auto & suggested_rotamer : suggested_rotamers ) {

		scoring::EnergyMap emap;
		for ( graph::Graph::EdgeListConstIter
				ir  = neighbor_graph->get_node( seqpos )->const_edge_list_begin(),
				ire = neighbor_graph->get_node( seqpos )->const_edge_list_end();
				ir != ire; ++ir ) {

			int const neighbor_id( (*ir)->get_other_ind( seqpos ) );
			conformation::Residue const & neighbor( pose.residue( neighbor_id ) );

			scofx->bump_check_backbone( *suggested_rotamer, neighbor, pose, emap );
		}
		core::Real thisrotE = scofx->weights().dot( emap ) ;

		if (  thisrotE < bb_bump_cutoff ) {
			temp_accepted_rotamers.push_back( suggested_rotamer );
			if ( thisrotE < bestE ) {
				best_rot = suggested_rotamer;
				bestE = thisrotE;
			}
		}
		//debug
		//std::cerr << "rotamer has energy " << (scofx->weights().dot( emap )) << " ... " << std::endl;

	}

	//we half-sort the output array such that the lowest energy rotamer is the first one in the output array
	if ( temp_accepted_rotamers.size() > 0 ) {
		accepted_rotamers.push_back( best_rot );
		for ( auto & temp_accepted_rotamer : temp_accepted_rotamers ) {

			if ( temp_accepted_rotamer != best_rot ) accepted_rotamers.push_back( temp_accepted_rotamer );
		}

	}


} //get_non_bb_clashing_rotamers


} // namespace ligand_docking
} // namespace protocols
