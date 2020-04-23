// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/rna/movers/RNAIdealizeMover.cc
/// @brief Slowly accomodate movement from non-ideal to ideal bond lengths and angles by repeated minimization
/// @author Andy Watkins (amw579@nyu.edu)

// Unit headers
#include <protocols/rna/movers/RNAIdealizeMover.hh>
#include <protocols/rna/movers/RNAIdealizeMoverCreator.hh>

// Core headers
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/RNA_SuiteName.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/CoordinateConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID.hh>

//#include <protocols/recces/sampler/rna/MC_RNA_KIC_Sampler.hh>
#include <protocols/stepwise/sampler/rna/RNA_KIC_Sampler.hh>

#include <protocols/stepwise/monte_carlo/mover/StepWiseMove.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>

#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.hh>
#include <protocols/stepwise/monte_carlo/mover/TransientCutpointHandler.hh>

#include <protocols/minimization_packing/MinMover.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/angle.functions.hh>
#include <numeric/conversions.hh>
#include <numeric/random/random.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static basic::Tracer TR( "protocols.rna.movers.RNAIdealizeMover" );

namespace protocols {
namespace rna {
namespace movers {

using namespace core;
using namespace core::id;
using namespace core::scoring;
using namespace core::conformation;
using namespace scoring::constraints;
using namespace scoring::func;

RNAIdealizeMover::RNAIdealizeMover():
	protocols::moves::Mover( RNAIdealizeMover::mover_name() )
{}

RNAIdealizeMover::~RNAIdealizeMover()= default;

void
RNAIdealizeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap&
)
{
	iterations_ = tag->getOption< Size >( "iterations", iterations_ );
	noise_ = tag->getOption< bool >( "noise", noise_ );
	final_minimize_ = tag->getOption< bool >( "final_minimize", final_minimize_ );
	ang_significance_threshold_ = tag->getOption< Real >( "ang_significance_threshold", ang_significance_threshold_ );
}

protocols::moves::MoverOP
RNAIdealizeMover::clone() const
{
	return utility::pointer::make_shared< RNAIdealizeMover >( *this );
}

protocols::moves::MoverOP
RNAIdealizeMover::fresh_instance() const
{
	return utility::pointer::make_shared< RNAIdealizeMover >();
}



void
RNAIdealizeMover::show( std::ostream & output ) const
{
	protocols::moves::Mover::show( output );
}

std::ostream &
operator<<( std::ostream & os, RNAIdealizeMover const & mover )
{
	mover.show(os);
	return os;
}

/// @brief Add Gaussian noise to the xyz coordinates of each atom.
/// @details Intentionally uses set_xyz because we don't want changes to propagate
void
RNAIdealizeMover::perturb_pose( pose::Pose & pose ) const
{
	using namespace numeric;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		for ( Size jj = 1; jj <= pose.residue_type( ii ).natoms(); ++jj ) {
			pose.conformation().set_xyz(
				AtomID( jj, ii ),
				pose.conformation().xyz( AtomID( jj, ii ) )
				+ xyzVector< core::Real >(
				random::rg().gaussian() * 0.02,
				random::rg().gaussian() * 0.02,
				random::rg().gaussian() * 0.02 ) );
		}
	}
}

void
RNAIdealizeMover::apply( pose::Pose & pose )
{
	using namespace core::chemical::rna;

	TR << "Idealizing pose with " << iterations_ << " iterations." << std::endl;
	TR << "A final round of minimization " << ( final_minimize_ ? "will" : "won't" ) << " follow." << std::endl;
	TR << "Noise " << ( noise_ ? "is" : "isn't" ) << " added to starting coords." << std::endl;

	ScoreFunctionOP scorefxn = get_score_function();

	scorefxn->set_weight( core::scoring::coordinate_constraint, 1.0 );

	// Perturb input pose with Gaussian xyz noise
	if ( noise_ ) perturb_pose( pose );

	// Store original cst set and foldtree
	auto const cst_set = pose.constraint_set()->clone();
	auto const orig_ft = pose.fold_tree();

	Pose ideal_pose;
	// Pose copy_pose;// = pose;
	for ( Size imax = 1 ; imax <= pose.size(); ++imax ) {
		if ( imax == 1 || !pose.residue( imax ).is_bonded( imax-1 ) ) {
			ideal_pose.append_residue_by_jump( core::conformation::Residue( pose.residue_type( imax ), true ), 1 );
		} else {
			ideal_pose.append_residue_by_bond( core::conformation::Residue( pose.residue_type( imax ), true ), true );
		}

		// if ( imax == 1 || !pose.residue( imax ).is_bonded( imax-1 ) ) {
		//  copy_pose.append_residue_by_jump( pose.residue( imax ), 1 );
		// } else {
		//  copy_pose.append_residue_by_bond( pose.residue( imax ), true );
		// }

		// add constraints
		core::pose::addVirtualResAsRoot( ideal_pose );
		// core::pose::addVirtualResAsRoot( copy_pose );
		Size const my_anchor( ideal_pose.fold_tree().root() );

		FuncOP constraint_func = utility::pointer::make_shared< FlatHarmonicFunc >( 0, 0.25, 1 );
		// FuncOP constraint_func( new HarmonicFunc( 0, 1 ) );

		kinematics::MoveMapOP suite_mm = utility::pointer::make_shared< kinematics::MoveMap >();
		for ( Size kk = 0; kk <= 4; ++kk ) {
			suite_mm->set_bb( ideal_pose.size() - kk, true );
			suite_mm->set_chi( ideal_pose.size() - kk, true );
		}
		// suite_mm->set_bb( true );
		// suite_mm->set_chi( true );
		suite_mm->set_jump( true );
		protocols::minimization_packing::MinMoverOP minm = utility::pointer::make_shared< protocols::minimization_packing::MinMover >( suite_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );

		// Clone dihedrals, just for imax and imax - 1
		if ( imax > 1 && pose.residue( imax ).is_bonded( imax-1 ) ) {
			for ( Size jj = 1; jj <= ideal_pose.residue_type( imax - 1 ).mainchain_atoms().size(); ++jj ) {
				ideal_pose.set_torsion( core::id::TorsionID( imax - 1, core::id::BB, jj ), pose.torsion( core::id::TorsionID( imax - 1, core::id::BB, jj ) ) );
			}
		}
		for ( Size jj = 1; jj <= ideal_pose.residue_type( imax ).mainchain_atoms().size(); ++jj ) {
			ideal_pose.set_torsion( core::id::TorsionID( imax, core::id::BB, jj ), pose.torsion( core::id::TorsionID( imax, core::id::BB, jj ) ) );
		}

		for ( Size jj = 1; jj <= ideal_pose.residue_type( imax ).nchi(); ++jj ) {
			ideal_pose.set_torsion( core::id::TorsionID( imax, core::id::CHI, jj ), pose.torsion( core::id::TorsionID( imax, core::id::CHI, jj ) ) );
		}


		// for ( Size ii = 1; ii < ideal_pose.size(); ++ii ) {

		for ( Size jj = 1; jj <= ideal_pose.residue_type( ideal_pose.size()-1 ).natoms(); ++jj ) {
			ConstraintOP constraint = utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( jj, ideal_pose.size()-1 ), core::id::AtomID( 1, my_anchor ),
				pose.residue( ideal_pose.size()-1 ).xyz( jj ),
				constraint_func );
			ideal_pose.add_constraint( constraint );
			// copy_pose.add_constraint( constraint );
		}

		// minimize ideal pose, I guess.
		// }
		TR << "Minimized ideal pose " << imax << " from " << ( *scorefxn )( ideal_pose ) << " to ";
		minm->apply( ideal_pose );
		TR << ( *scorefxn )( ideal_pose ) << "." << std::endl;

		// TR << "Minimized orig pose " << imax << " from " << ( *scorefxn )( copy_pose ) << " to ";
		// minm->apply( copy_pose );
		// TR << ( *scorefxn )( copy_pose ) << "." << std::endl;

		// std::stringstream ss;
		// ss << "min_ideal_" << imax << ".pdb";
		// ideal_pose.dump_pdb( ss.str() );
		core::pose::remove_virtual_residues( ideal_pose );

		// ss.clear();
		// ss << "min_orig_" << imax << ".pdb";
		// copy_pose.dump_pdb( ss.str() );
		// core::pose::remove_virtual_residues( copy_pose );
	}

	core::pose::addVirtualResAsRoot( ideal_pose );
	// core::pose::addVirtualResAsRoot( copy_pose );
	Size const my_anchor( ideal_pose.fold_tree().root() );

	FuncOP constraint_func = utility::pointer::make_shared< FlatHarmonicFunc >( 0, 1, 1 );

	kinematics::MoveMapOP suite_mm = utility::pointer::make_shared< kinematics::MoveMap >();
	suite_mm->set_bb( true );
	suite_mm->set_chi( true );
	suite_mm->set_jump( true );
	protocols::minimization_packing::MinMoverOP minm = utility::pointer::make_shared< protocols::minimization_packing::MinMover >( suite_mm, scorefxn, "lbfgs_armijo_nonmonotone", 0.001, true );

	for ( Size ii = 1; ii < ideal_pose.size(); ++ii ) {

		for ( Size jj = 1; jj <= ideal_pose.residue_type( ii ).natoms(); ++jj ) {
			ConstraintOP constraint = utility::pointer::make_shared< CoordinateConstraint >( core::id::AtomID( jj, ii ), core::id::AtomID( 1, my_anchor ),
				pose.residue( ii ).xyz( jj ),
				constraint_func );
			ideal_pose.add_constraint( constraint );
		}

	}
	minm->apply( ideal_pose );
	ideal_pose.dump_pdb( "final_ideal_1.pdb" );

	// for ( Size ii = 1; ii < copy_pose.size(); ++ii ) {

	//  for ( Size jj = 1; jj <= copy_pose.residue_type( ii ).natoms(); ++jj ) {
	//   ConstraintOP constraint( new CoordinateConstraint( core::id::AtomID( jj, ii ), core::id::AtomID( 1, my_anchor ),
	//    pose.residue( ii ).xyz( jj ),
	//    constraint_func ) );
	//   copy_pose.add_constraint( constraint );
	//  }

	// }
	// minm->apply( copy_pose );
	// copy_pose.dump_pdb( "final_orig_1.pdb" );

	// OK, reminimize but with just C4' constraints.

	constraint_func = utility::pointer::make_shared< FlatHarmonicFunc >( 0, 1.5, 1 ); //ideal_pose.constraint_set( cst_set );
	minm->apply( ideal_pose );
	// minm->apply( copy_pose );
	// ideal_pose.dump_pdb( "final_2.pdb" );

	// now, looser constraints.
	constraint_func = utility::pointer::make_shared< FlatHarmonicFunc >( 0, 2, 2 );
	minm->apply( ideal_pose );
	// minm->apply( copy_pose );
	// ideal_pose.dump_pdb( "final_3.pdb" );

	core::pose::remove_virtual_residues( ideal_pose );
	// core::pose::remove_virtual_residues( copy_pose );

	pose = ideal_pose;
	// Restore original constraints
	pose.constraint_set( cst_set );
	pose.fold_tree( orig_ft );
}

/////////////// Creator ///////////////



std::string RNAIdealizeMover::get_name() const {
	return mover_name();
}

std::string RNAIdealizeMover::mover_name() {
	return "RNAIdealizeMover";
}

void RNAIdealizeMover::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute( "iterations", xsct_positive_integer, "Number of iterations over which to spread idealization" )
		+ XMLSchemaAttribute( "noise", xsct_rosetta_bool, "Salt initial pose with some Gaussian noise" )
		+ XMLSchemaAttribute( "final_minimize", xsct_rosetta_bool, "Minimize the whole pose at once at the end" )
		+ XMLSchemaAttribute( "ang_significance_threshold", xsct_real, "Size of angle deviation to correct (degrees)" );

	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(), "Idealize a pose, moving its bond lengths and angles towards ideal values by repeated relaxation", attlist );
}

std::string RNAIdealizeMoverCreator::keyname() const {
	return RNAIdealizeMover::mover_name();
}

protocols::moves::MoverOP
RNAIdealizeMoverCreator::create_mover() const {
	return utility::pointer::make_shared< RNAIdealizeMover >();
}

void RNAIdealizeMoverCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	RNAIdealizeMover::provide_xml_schema( xsd );
}


} //movers
} //rna
} //protocols

