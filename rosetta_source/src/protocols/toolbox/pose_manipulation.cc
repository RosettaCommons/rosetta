// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/toolbox/pose_manipulation.cc
/// @brief some general functions to manipulate poses. feel free to add your own
/// @brief if you add your own, please mention your name:)
/// @author Florian Richter, floric@u.washington.edu
/// @author Steven Lewis, smlewi@gmail.com (insert_pose_into_pose) domain insertion code

// Unit headers
#include <protocols/toolbox/pose_manipulation.hh>

// Project headers
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
//#include <core/chemical/VariantType.hh>
 //needed for adding variant types

#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/id/AtomID_Map.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh> //typeset swapping

#include <protocols/moves/MonteCarlo.hh>

#include <protocols/loops/loop_closure/ccd/CcdLoopClosureMover.hh>
#include <protocols/loops/Loops.hh>


// Utility Headers
// AUTO-REMOVED #include <basic/options/option.hh>
#include <basic/Tracer.hh>
#include <core/types.hh>

#include <core/chemical/VariantType.hh>
#include <core/pose/util.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/util.tmpl.hh>

// C++ Headers

namespace protocols {
namespace toolbox {
namespace pose_manipulation{

static basic::Tracer TR("protocols.toolbox.pose_manipulation");
static basic::Tracer TR_DI("protocols.toolbox.pose_manipulation.insert_pose_into_pose");
using basic::T;
using basic::Error;
using basic::Warning;
using core::chemical::ResidueType;

void
construct_poly_ala_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	construct_poly_uniq_restype_pose( pose, positions, core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map("ALA"), keep_pro, keep_gly, keep_disulfide_cys);
}


void
construct_poly_uniq_restype_pose(
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	ResidueType const & restype,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	using namespace core;

	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	conformation::Residue const replace_res( restype, true );

	for( utility::vector1< Size >::const_iterator pos_it = positions.begin();
			 pos_it != positions.end(); ++pos_it )
		{

			chemical::ResidueTypeCOP cur_restype = & pose.residue_type( *pos_it );


			if( ( keep_pro && ( cur_restype->aa() == chemical::aa_pro ) )
				||( keep_gly && ( cur_restype->aa() == chemical::aa_gly ) )
				||( keep_disulfide_cys && ( cur_restype->aa() == chemical::aa_cys ) && cur_restype->has_variant_type( chemical::DISULFIDE ) ) )
				{
					continue;
				}

			utility::vector1< std::string > current_variants;

			bool variants_match = cur_restype->variants_match( replace_res.type() );
			//TR<< "replacing: " << *pos_it << std::endl;

			if( !variants_match ){

				current_variants = cur_restype->variant_types();
				chemical::ResidueType var_replace_type = replace_res.type();

				for(core::Size var = 1; var <= current_variants.size(); var++){
					var_replace_type = restype_set->get_residue_type_with_variant_added( var_replace_type, current_variants[ var ] );
				}

				//runtime_assert( var_replace_type.name3() == "ALA" );
				conformation::Residue const var_replace_res( restype_set->name_map( var_replace_type.name() ), true );

				pose.replace_residue( *pos_it, var_replace_res, true );
			}

			else	pose.replace_residue( *pos_it, replace_res, true);

		} //iterator over positions to replace

} // construct_poly_ala_pose function

void
construct_poly_XXX_pose(
	std::string const & aa,
	core::pose::Pose & pose,
	utility::vector1< core::Size > const & positions,
	bool keep_pro,
	bool keep_gly,
	bool keep_disulfide_cys
)
{
	using namespace core;

	chemical::ResidueTypeSetCAP restype_set = chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );

	conformation::Residue const replace_res( restype_set->name_map( aa ), true );

	for( utility::vector1< Size >::const_iterator pos_it = positions.begin();
			 pos_it != positions.end(); ++pos_it )
		{

			chemical::ResidueTypeCOP cur_restype = & pose.residue_type( *pos_it );


			if( ( keep_pro && ( cur_restype->aa() == chemical::aa_pro ) )
				||( keep_gly && ( cur_restype->aa() == chemical::aa_gly ) )
				||( keep_disulfide_cys && ( cur_restype->aa() == chemical::aa_cys ) && cur_restype->has_variant_type( chemical::DISULFIDE ) ) )
				{
					continue;
				}

			utility::vector1< std::string > current_variants;

			bool variants_match = cur_restype->variants_match( replace_res.type() );

			if( !variants_match ){

				current_variants = cur_restype->variant_types();
				chemical::ResidueType var_replace_type = replace_res.type();

				for(core::Size var = 1; var <= current_variants.size(); var++){
					var_replace_type = restype_set->get_residue_type_with_variant_added( var_replace_type, current_variants[ var ] );
				}

				runtime_assert( var_replace_type.name3() == aa );
				conformation::Residue const var_replace_res( restype_set->name_map( var_replace_type.name() ), true );

				pose.replace_residue( *pos_it, var_replace_res, true );
			}

			else	pose.replace_residue( *pos_it, replace_res, true);

		} //iterator over positions to replace

} // construct_poly_XXX_pose function


void
remove_non_protein_residues(
	core::pose::Pose & pose
)
{

	bool residues_deleted(false);

	for( core::Size i = pose.total_residue(); i > 0 ; --i){

		if( ! pose.residue_type( i ).is_protein() ){
			pose.conformation().delete_residue_slow( i );
			residues_deleted = true;
		}
	}

	if( residues_deleted ) pose.energies().clear();

} //remove_non_protein_residues


void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS ) ) continue;

		if ( !pose.residue_type( this_cutpoint ).is_protein() ) continue;
		if ( !pose.residue_type( this_cutpoint +1 ).is_protein() ) continue;

		if( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}

void
add_chainbreaks_according_to_jumps( core::pose::Pose & pose, utility::vector1< core::Size > const& no_cutpoint_residues  )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );
		//exclude residue numbers in array
		if ( find( no_cutpoint_residues.begin(), no_cutpoint_residues.end(), this_cutpoint ) != no_cutpoint_residues.end() ) {
			continue;
		}

		if ( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::UPPER_TERMINUS ) ) continue;

		if ( pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::LOWER_TERMINUS ) ) continue;

		if( !pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( !pose.residue_type( this_cutpoint +1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::add_variant_type_to_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint +1 );
		}
	}
}


void
remove_chainbreaks_according_to_jumps( core::pose::Pose & pose )
{
	for( core::Size i =1; i <= pose.fold_tree().num_jump(); ++i ){

		core::Size this_cutpoint( pose.fold_tree().cutpoint( i ) );

		if( pose.residue_type( this_cutpoint ).has_variant_type( core::chemical::CUTPOINT_LOWER ) ){
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_LOWER, this_cutpoint );
		}
		if( pose.residue_type( this_cutpoint + 1).has_variant_type( core::chemical::CUTPOINT_UPPER ) ){
			core::pose::remove_variant_type_from_pose_residue( pose, core::chemical::CUTPOINT_UPPER, this_cutpoint+1 );
		}
	}
}

core::Real
superimpose_pose_on_subset_CA(
	core::pose::Pose & pose,
	core::pose::Pose const & ref_pose,
	utility::vector1< core::Size > const & positions,
	int const offset
)
{
	using namespace core::id;

	AtomID_Map< AtomID > atom_map;

	core::pose::initialize_atomid_map( atom_map, pose, BOGUS_ATOM_ID );

	for( utility::vector1< core::Size >::const_iterator res_it = positions.begin(); res_it != positions.end(); ++res_it){

		AtomID id1( pose.residue( *res_it + offset).atom_index("CA"), *res_it + offset );
		AtomID id2( ref_pose.residue( *res_it ).atom_index("CA"), *res_it );

		atom_map.set( id1, id2);

	} //iterator over residues

	return core::scoring::superimpose_pose( pose, ref_pose, atom_map );
} //superimpose_pose_on_subset_CA

///@author Steven Lewis smlewi@gmail.com
///@details brief inserts one pose into another pose, returning the product as a new value.  This is basically a seed for a domain insertion application.  The three core::Size arguments define a flexible surface loop on the scaffold, the insert pose will be added immediately after insert_point.  The insert will be left unchanged in internal-coordinate space except for the phi on the first residue, and the psi/omega on the last residue, and atoms whose bonding partners change as a result of the insertion.  Note that insert_loop_end is pass-by-reference: this field will be updated to reflect its numbering in the result pose (post-insertion).  Internally, this function performs the insertion, idealizes the loop residues (omegas to 180, peptide bonds idealized) and the newly made polymer connections at the insert point, and then attempts to close the loop (chainbreak is next to the end of the loop).  It is intended, but not guarunteed, to produce a loop with good rama, omega, and chainbreak/peptide_bond scores.  If there is a bad chainbreak, it will be at the return value in (insert_loop_end-1).  It does NOT attempt to give a loop with good sidechains (it does not repack at all) or worry overmuch about van der Waals.
core::pose::Pose
insert_pose_into_pose(
											core::pose::Pose const & scaffold_pose,
											core::pose::Pose const & insert_pose,
											core::Size const insert_loop_start,
											core::Size const insert_point,
											core::Size & insert_loop_end,
											core::Size const cycles
){
	//local copies
	core::pose::Pose scaffold(scaffold_pose);
	core::pose::Pose insert(insert_pose);
	core::Size const insertlength = insert.total_residue();
	//strip termini variants from insert if necessary
	using core::pose::remove_variant_type_from_pose_residue;
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::LOWER_TERMINUS, 1);
	core::pose::remove_variant_type_from_pose_residue(insert, core::chemical::UPPER_TERMINUS, insertlength);

	//calculate relevant portions of the loop
	//change name of insert_loop_start
	core::Size const loop_start(insert_loop_start); //this position is the first flexible member of the loop
	core::Size const loop_start_foldtree_anchor(loop_start-1); //this is the N-terminal jump anchor for the loop
	//insert_point passed into function
	//core::Size const insert_point; //this is the last scaffold residue before the insert
	core::Size const insert_start(insert_point+1); //this is the first residue of the inflexible insert
	core::Size const insert_end(insert_point+insertlength); //the last residue of the inflexible insert
	core::Size const loop_end(insert_loop_end+insertlength); //last flexible loop residue
	//echo back out to argument insert_loop_end
	insert_loop_end = loop_end;
	core::Size const cutpoint_lower(loop_end-1); //cutpoint at the end of the loop
	core::Size const cutpoint_upper(loop_end); //cutpoint at the end of the loop
	core::Size const loop_end_foldtree_anchor(loop_end+1); //C-terminal jump anchor

	TR_DI << "loop_start" << " " << "insert_start" << " " << "loop_end" << std::endl;
	TR_DI << "scaffold numbering " <<  insert_loop_start  << " " <<  insert_point  << " " << insert_loop_end  << std::endl;
	TR_DI << "insertlength " << insertlength << "\nmodified numbering" << std::endl;
	TR_DI << "loop_start " << loop_start << std::endl;
	TR_DI << "loop_start_foldtree_anchor " << loop_start_foldtree_anchor << std::endl;
	TR_DI << "insert_point " << insert_point << std::endl;
	TR_DI << "insert_start " << insert_start << std::endl;
	TR_DI << "insert_end " << insert_end << std::endl;
	TR_DI << "loop_end " << loop_end << std::endl;
	TR_DI << "cutpoint_lower " << cutpoint_lower << std::endl;
	TR_DI << "cutpoint_upper " << cutpoint_upper << std::endl;
	TR_DI << "loop_end_foldtree_anchor " << loop_end_foldtree_anchor << std::endl;

	//Fold tree allows insertion into scaffold (all via jump)
	using core::kinematics::Edge;
	core::pose::Pose combined(scaffold);
	core::kinematics::FoldTree inserting_tree(combined.total_residue());
	inserting_tree.clear();
	inserting_tree.add_edge(Edge(1, insert_point, Edge::PEPTIDE));
	inserting_tree.add_edge(Edge(insert_point+1, combined.total_residue(), Edge::PEPTIDE));
	inserting_tree.add_edge(Edge(insert_point, insert_point+1, 1));
	inserting_tree.reorder(1);
	TR_DI << inserting_tree << std::endl;
	combined.fold_tree(inserting_tree);

	//insert insert residues by jump (just to get them in)
	for (core::Size i = 1; i <= insertlength; ++i){
		combined.insert_residue_by_jump(insert.residue(i), insert_point+i, loop_start);
	}

	//This fold tree is the proper fold tree for loop remodeling
	core::kinematics::FoldTree remodeling_tree(combined.total_residue());
	remodeling_tree.clear();
	remodeling_tree.add_edge(Edge(1, loop_start_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, cutpoint_lower, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(cutpoint_upper, loop_end_foldtree_anchor, Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_end_foldtree_anchor, combined.total_residue(), Edge::PEPTIDE));
	remodeling_tree.add_edge(Edge(loop_start_foldtree_anchor, loop_end_foldtree_anchor, 1));
	remodeling_tree.reorder(1);
	TR_DI << remodeling_tree << std::endl;
	combined.fold_tree(remodeling_tree);

	//label cutpoints properly - necessary for CCD
	using core::pose::add_variant_type_to_pose_residue;
	core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_LOWER, cutpoint_lower );
	core::pose::add_variant_type_to_pose_residue(combined, core::chemical::CUTPOINT_UPPER, cutpoint_upper );

	protocols::loops::Loops loops;
	loops.add_loop(loop_start, loop_end, cutpoint_lower);
	TR_DI << loops << std::endl;

	//create movemap: loops mobile, insert not mobile
	//this code also resets conformation variables: omegas to 180, newly made connections phi or psi to reasonable
	//edges of insert will be somewhat mobile inside minimization (small and CCD moves will ignore it)
	using namespace core::id;
	core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap();
	core::Size i(loop_start);
	for(; i<=insert_point; ++i) {
		movemap->set_bb(i, true);
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
		combined.set_omega(i, 180);
		TR_DI << "mobile " << i << std::endl;
	}
	TR_DI << "ideal " << insert_start << std::endl;
 	movemap->set( TorsionID(insert_start, BB, phi_torsion), true);
	combined.set_phi(insert_start, -60);

 	movemap->set( TorsionID(insert_end, BB, psi_torsion), true);
	combined.conformation().insert_ideal_geometry_at_polymer_bond(insert_end);
	TR_DI << "ideal " << insert_end << std::endl;
	combined.set_omega(insert_end, 180);
	combined.set_psi(insert_end, -40);
	i = insert_end+1;
	for(; i<=loop_end; ++i){
		movemap->set_bb(i, true);
		//movemap->set( TorsionID(i, BB, omega_torsion), false ); //fixes omega angle
		combined.conformation().insert_ideal_geometry_at_polymer_bond(i);
		combined.set_omega(i, 180);
		TR_DI << "mobile " << i << std::endl;
	}

	//centroidize the pose before we do stuff to it - sidechains are expensive and unnecessary
	protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap(core::chemical::CENTROID);
	protocols::simple_moves::ReturnSidechainMover return_sidechains( combined );
	typeset_swap.apply( combined );

	//combined.dump_pdb("combined_preclose_cen.pdb");

	//create CCD mover, smallmover
	protocols::loops::loop_closure::ccd::CcdLoopClosureMover close( *(loops.begin()), movemap );
	protocols::simple_moves::SmallMover small(movemap, 10, 200); //huge moves for sampling
	small.angle_max( 'H', 180.0 );
	small.angle_max( 'E', 180.0 );
	small.angle_max( 'L', 180.0 );

	//this scorefunction is pretty arbitrary - feel free to hack around with it
	using namespace core::scoring;
	ScoreFunctionOP scfxn = new ScoreFunction;
	(*scfxn).set_weight( chainbreak, 20.0 );
	(*scfxn).set_weight( cbeta,       1.0 );
	(*scfxn).set_weight( vdw,         1.0 );
	(*scfxn).set_weight( pair,        1.0 );
	(*scfxn).set_weight( cenpack,     1.0 );
	(*scfxn).set_weight( rama,        5.0 );
	(*scfxn).set_weight( hbond_lr_bb, 1.0 );
	(*scfxn).set_weight( hbond_sr_bb, 1.0 );
	(*scfxn).set_weight( omega,       5.0 );

	protocols::simple_moves::MinMover min_mover(movemap, scfxn, "dfpmin_armijo", 0.01, true /*use_nblist*/ );

	/////////////////////////Monte Carlo//////////////////////////////////////////////////////////
	//make the monte carlo object
	using protocols::moves::MonteCarlo;
	using protocols::moves::MonteCarloOP;
	MonteCarlo mc(combined, (*scfxn), 0.8);

	TR_DI << "start " << ((*scfxn))(combined) << std::endl;
	for( core::Size i(1); i<=cycles; ++i){

		small.apply(combined);
		close.apply(combined);
		combined.conformation().insert_ideal_geometry_at_polymer_bond(cutpoint_lower);
		min_mover.apply(combined);
		combined.conformation().insert_ideal_geometry_at_polymer_bond(cutpoint_lower);

		if(mc.boltzmann(combined)) TR_DI << i << " " << ((*scfxn))(combined) << std::endl;
	}
	mc.recover_low(combined);
	TR_DI << "finish " << ((*scfxn))(combined) << std::endl;
	//combined.conformation().insert_ideal_geometry_at_polymer_bond(cutpoint_lower);
	//TR_DI << "finish " << ((*scfxn))(combined) << std::endl;

	//combined.dump_pdb("combined_closed_cen.pdb");
	return_sidechains.apply( combined );
	//combined.dump_pdb("combined_closed.pdb");

	return combined;
}


} // namespace pose_manipulation
} //toolbox
} //protocols
