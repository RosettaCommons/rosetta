// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    protocols/grafting/CCDEndsGraftMover.cc
/// @brief   Method definitions for CCDEndsGraftMover
/// @author  Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/grafting/CCDEndsGraftMover.hh>
#include <protocols/grafting/CCDEndsGraftMoverCreator.hh>
#include <protocols/grafting/AnchoredGraftMover.hh>
#include <protocols/grafting/util.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/chemical/VariantType.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/selection.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>

#include <protocols/loops/FoldTreeFromLoopsWrapper.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>

#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.grafting.CCDEndsGraftMover" );


namespace protocols {
namespace grafting {

using namespace core::scoring;
using namespace core::pose;
using namespace core::pack::task;
using namespace protocols::loops;
using namespace protocols::simple_moves;
using namespace core::kinematics;

using core::pose::Pose;
using core::kinematics::MoveMapCOP;
using core::kinematics::MoveMapOP;
using core::scoring::ScoreFunctionCOP;
using core::scoring::ScoreFunctionOP;
using protocols::simple_moves::MinMoverOP;
using protocols::simple_moves::SmallMoverOP;
using core::Size;

CCDEndsGraftMover::CCDEndsGraftMover():
	AnchoredGraftMover(0u, 0u) // 0u to avoid spurious rewrite
{
	Nter_overhang_length(2);
	Cter_overhang_length(2);
}

CCDEndsGraftMover::CCDEndsGraftMover(const Size start, Size const end, bool copy_pdb_info /*false*/)
:AnchoredGraftMover(start, end)
{
	Nter_overhang_length(2);
	Cter_overhang_length(2);

	copy_pdbinfo(copy_pdb_info);
}

CCDEndsGraftMover::CCDEndsGraftMover(
	const Size start,
	const Size end,
	core::pose::Pose const & piece,
	Size Nter_overhang_length,
	Size Cter_overhang_length,
	bool copy_pdb_info /*false*/)
:AnchoredGraftMover(start, end, piece, Nter_overhang_length, Cter_overhang_length)
{
	copy_pdbinfo(copy_pdb_info);
}

CCDEndsGraftMover::CCDEndsGraftMover(const CCDEndsGraftMover& src ):
	AnchoredGraftMover( src )
{

}

CCDEndsGraftMover::~CCDEndsGraftMover() = default;

// XRW TEMP std::string
// XRW TEMP CCDEndsGraftMover::get_name() const { return "CCDEndsGraftMover"; }

// XRW TEMP protocols::moves::MoverOP
// XRW TEMP CCDEndsGraftMoverCreator::create_mover() const {
// XRW TEMP  return protocols::moves::MoverOP( new CCDEndsGraftMover );
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CCDEndsGraftMoverCreator::keyname() const {
// XRW TEMP  return CCDEndsGraftMover::mover_name();
// XRW TEMP }

// XRW TEMP std::string
// XRW TEMP CCDEndsGraftMover::mover_name(){
// XRW TEMP  return "CCDEndsGraftMover";
// XRW TEMP }

protocols::moves::MoverOP
CCDEndsGraftMover::clone() const{
	return protocols::moves::MoverOP( new CCDEndsGraftMover(*this) );
}

protocols::moves::MoverOP
CCDEndsGraftMover::fresh_instance() const{
	return protocols::moves::MoverOP( new CCDEndsGraftMover() );
}

SmallMoverOP
CCDEndsGraftMover::setup_default_small_mover(){

	using namespace core::kinematics;
	using namespace protocols::simple_moves;
	//Smaller moves for sampling than AnchoredGraftMover - Just to help CCD a bit
	//Will check whether this does indeed help or hinder later.

	//200 seems extraordinarily high.
	// Lets take nmoves as twice the number of
	// movable residues in the movemap instead.

	core::Size nmoves = 0;

	for ( auto it=movemap()->movemap_torsion_id_begin(), it_end=movemap()->movemap_torsion_id_end(); it !=it_end; ++it ) {
		//Scaffold to new MM
		if ( it->second && it->first.second == core::id::BB ) {
			nmoves+=1;
		}
	}

	SmallMoverOP small( new SmallMover(movemap(), 2, nmoves*2) );
	small->angle_max( 'H', 10.0 );
	small->angle_max( 'E', 10.0 );
	small->angle_max( 'L', 30.0 );
	return small;
}

void
CCDEndsGraftMover::apply(Pose & pose){
	using namespace core::id;
	using namespace protocols::loops;
	using core::pose::add_variant_type_to_pose_residue;
	using protocols::moves::MonteCarlo;
	using protocols::moves::MonteCarloOP;

	//TR <<"Beginning of anchored graft mover" <<std::endl;
	//pose.constraint_set()->show(TR);

	if ( !scaffold_start().empty() && ! scaffold_end().empty() ) {
		core::Size scaff_start = core::pose::parse_resnum(scaffold_start(), pose);
		core::Size scaff_end = core::pose::parse_resnum(scaffold_end(), pose);
		set_insert_region(scaff_start, scaff_end);
	}

	TR << "Start: " << start() << " End: " << end() << " NterO: "<<Nter_overhang_length() <<" CterO: " << Cter_overhang_length()<<std::endl;

	protocols::grafting::superimpose_overhangs_heavy(pose, *piece(), true, start(), end(), Nter_overhang_length(), Cter_overhang_length() );


	Pose combined = insert_piece(pose);

	setup_movemap_and_regions(pose);

	TR << "Start: "<<start() << std::endl;
	TR << "Original End: " << original_end() << std::endl;
	TR << "End: " << end() << std::endl;
	TR << "Insert Length: " << insertion_length() << std::endl;

	core::kinematics::FoldTree original_ft = combined.fold_tree();
	//Setup for the remodeling
	core::Size const insert_start(start()+1); //this will be the first residue of the insert
	core::Size const insert_end(start()+insertion_length()); //this will be the last residue of the insert


	setup_scorefunction();

	///Add variants, create the loops and set the foldtree that will be used for CCD.


	Loop Nter_loop;
	Loop Cter_loop;
	LoopsOP loop_set( new Loops() );
	std::map< Loop, loop_closure::ccd::CCDLoopClosureMoverOP > loop_set_map; //Would not work without owning pointer.

	//Cut should always be at the position of the insert.

	//Nter_loop = Loop(Nter_loop_start(), Nter_loop_end()+1, Nter_loop_end());//(LEFT LOOP)
	//Cter_loop = Loop(Cter_loop_start()-1, Cter_loop_end(), Cter_loop_start()-1);//(RIGHT LOOP)


	//We need the loops to be one off from the cutpoint.  However, if we are insert flexibility, we are already OK.

	core::Size Nter_end;
	core::Size Cter_start;

	if ( get_nterm_insert_flexibility() == 0 ) {
		Nter_end = Nter_loop_end() +1;
	} else {
		Nter_end = Nter_loop_end();
	}

	if ( get_cterm_insert_flexibility() == 0 ) {
		Cter_start = Cter_loop_start() -1;
	} else {
		Cter_start = Cter_loop_start();
	}

	Nter_loop = Loop(Nter_loop_start(), Nter_end , start());//(LEFT LOOP)
	Cter_loop = Loop(Cter_start, Cter_loop_end(), end()-1);//(RIGHT LOOP)

	TR << Nter_loop << std::endl;
	TR << Cter_loop << std::endl;

	loop_set->add_loop(Nter_loop);
	loop_set->add_loop(Cter_loop);

	//movemap()->show(TR);
	//TR << "Loops: " << *loop_set << std::endl;

	if ( Cter_start - Nter_end < 1 ) {
		utility_exit_with_message("Cannot setup correct loops for CCDEndsGraftMover as the loop FT is not proper.\n"
			"The loop is short two loops on either side cannot be set.  Decrease Insert movement.  \n");
	}
	//combined.dump_pdb("combined_pose.pdb");

	FoldTreeFromLoops ft_loop = FoldTreeFromLoops();
	ft_loop.loops(loop_set);
	ft_loop.apply(combined);

	set_loops(loop_set);

	add_cutpoint_variants_for_ccd(combined, *loop_set);

	loop_set_map[Nter_loop] = protocols::loops::loop_closure::ccd::CCDLoopClosureMoverOP( new loop_closure::ccd::CCDLoopClosureMover(Nter_loop, movemap()) );
	loop_set_map[Cter_loop] = protocols::loops::loop_closure::ccd::CCDLoopClosureMoverOP( new loop_closure::ccd::CCDLoopClosureMover(Cter_loop, movemap()) );

	//combined.dump_pdb("before_idealize.pdb");

	idealize_combined_pose(combined, movemap(), start(), insert_start, insert_end, Nter_loop_start(), Cter_loop_end(), idealize_insert());
	movemap()->set( TorsionID(insert_start, BB, phi_torsion), true);
	movemap()->set( TorsionID(insert_end, BB, psi_torsion), true);

	//combined.dump_pdb("after_idealize.pdb");

	//centroidize the pose before we do stuff to it - sidechains are expensive and unnecessary
	protocols::simple_moves::SwitchResidueTypeSetMover typeset_swap(core::chemical::CENTROID);
	protocols::simple_moves::ReturnSidechainMoverOP return_sidechains( new  protocols::simple_moves::ReturnSidechainMover(combined ) );
	typeset_swap.apply( combined );

	//TR <<"After type swap" <<std::endl;

	//combined.dump_pdb("combined_preclose_cen.pdb");

	MinMoverOP min_mover = setup_default_min_mover();
	SmallMoverOP small = setup_default_small_mover();

	//Testing
	if ( test_control_mode() ) { perturb_backbone_for_test(combined, movemap());}

	MonteCarlo mc(combined, (*cen_scorefxn()), 0.8);


	/////////////////////////Protocol//////////////////////////////////////////////////////////
	TR << "start " << ((*cen_scorefxn()))(combined) << std::endl;

	bool skip_mc = false;
	for ( core::Size i(1); i<=cycles(); ++i ) {
		TR <<"round "<<i <<std::endl;
		if ( !skip_sampling() ) { small->apply(combined);}

		for ( auto const & it : *loop_set ) {

			loop_set_map[it]->apply(combined);
			combined.conformation().insert_ideal_geometry_at_polymer_bond(it.cut());
			min_mover->apply(combined);
			combined.conformation().insert_ideal_geometry_at_polymer_bond(it.cut());

		}
		if ( stop_at_closure() && graft_closed(combined, *loop_set) ) {
			TR << "Graft Closed early - returning" << std::endl;
			skip_mc = true;
			TR << i << " " << ((*cen_scorefxn()))(combined) << std::endl;
			break;
		}
		if ( mc.boltzmann(combined) ) TR << i << " " << ((*cen_scorefxn()))(combined) << std::endl;
	}

	if ( ! skip_mc ) mc.recover_low(combined);
	TR << "finish " << ((*cen_scorefxn()))(combined) << std::endl;
	//combined.conformation().insert_ideal_geometry_at_polymer_bond(Cter_loop.cut());

	return_sidechains->apply( combined );

	//Remove cutpoints that were required for CCD.
	remove_cutpoint_variants_for_ccd(combined, *loop_set);

	//Give back foldtree from pose_into_pose.
	combined.fold_tree(original_ft);
	if ( final_repack() ) {
		repack_connection_and_residues_in_movemap_and_piece_and_neighbors( pose, fa_scorefxn(),
			start(), end(), movemap(), neighbor_dis());
	}
	TR <<"Graft meets ideal geometry " << std::boolalpha << graft_closed(combined, *loop_set) << std::endl;
	TR << "Complete"<<std::endl;
	pose = combined;

}

std::string CCDEndsGraftMover::get_name() const {
	return mover_name();
}

std::string CCDEndsGraftMover::mover_name() {
	return "CCDEndsGraftMover";
}

void CCDEndsGraftMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	XMLSchemaComplexTypeGeneratorOP ct_gen = AnchoredGraftMover::complex_type_generator_for_anchored_graft_mover( xsd );
	ct_gen->element_name( mover_name() )
		.complex_type_naming_func( & moves::complex_type_name_for_mover )
		.description( "AnchoredGraftMover that uses CCD to close loops" )
		.write_complex_type_to_schema( xsd );
}

std::string CCDEndsGraftMoverCreator::keyname() const {
	return CCDEndsGraftMover::mover_name();
}

protocols::moves::MoverOP
CCDEndsGraftMoverCreator::create_mover() const {
	return protocols::moves::MoverOP( new CCDEndsGraftMover );
}

void CCDEndsGraftMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	CCDEndsGraftMover::provide_xml_schema( xsd );
}


} //grafting
} //protocols
