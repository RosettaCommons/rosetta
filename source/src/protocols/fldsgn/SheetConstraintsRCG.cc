// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/SheetConstraintsRCG.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu), Nov 2012

// Unit header
#include <protocols/fldsgn/SheetConstraintsRCG.hh>
#include <protocols/fldsgn/SheetCstGeneratorCreator.hh>

// Package headers
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ headers
#include <cmath>

static basic::Tracer TR( "protocols.fldsgn.SheetConstraintsRCG" );

namespace protocols{
namespace fldsgn{

std::string
SheetCstGeneratorCreator::keyname() const
{
	return SheetCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
SheetCstGeneratorCreator::create_mover() const {
	return new SheetConstraintsRCG();
}

std::string
SheetCstGeneratorCreator::mover_name()
{
	return "SheetCstGenerator";
}

/// @brief
SheetConstraintsRCG::SheetConstraintsRCG():
	RemodelConstraintGenerator(),
	weight_( 1.0 ),
	dist_( 5.5 ),
	angle_tolerance_( 0.35 ),
	cacb_dihedral_tolerance_( 0.9 ),
	bb_dihedral_tolerance_( 0.52 ),
	constrain_dist_only_( false ),
	blueprint_( NULL )
{}

/// @ copy constructor
SheetConstraintsRCG::SheetConstraintsRCG( SheetConstraintsRCG const & rval ):
	RemodelConstraintGenerator( rval ),
	weight_( rval.weight_ ),
	dist_( rval.dist_ ),
	angle_tolerance_( rval.angle_tolerance_ ),
	cacb_dihedral_tolerance_( rval.cacb_dihedral_tolerance_ ),
	bb_dihedral_tolerance_( rval.bb_dihedral_tolerance_ ),
	constrain_dist_only_( rval.constrain_dist_only_ ),
	blueprint_( rval.blueprint_ )
{}

/// @brief
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue ):
	RemodelConstraintGenerator(),
	weight_( 1.0 ),
	dist_( 5.5 ),
	angle_tolerance_( 0.35 ),
	cacb_dihedral_tolerance_( 0.9 ),
	bb_dihedral_tolerance_( 0.52 ),
	constrain_dist_only_( false ),
	blueprint_( blue )
{}

/// @brief value constructor
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue, Real const coef ):
	RemodelConstraintGenerator(),
	weight_( coef ),
	dist_( 5.5 ),
	angle_tolerance_( 0.35 ),
	cacb_dihedral_tolerance_( 0.9 ),
	bb_dihedral_tolerance_( 0.52 ),
	constrain_dist_only_( false ),
	blueprint_( blue )
{}

/// @brief value constructor
SheetConstraintsRCG::SheetConstraintsRCG( BluePrintOP const & blue, Real const coef, Real const dist ):
	RemodelConstraintGenerator(),
	weight_( coef ),
	dist_( dist ),
	angle_tolerance_( 0.35 ),
	cacb_dihedral_tolerance_( 0.9 ),
	bb_dihedral_tolerance_( 0.52 ),
	constrain_dist_only_( false ),
	blueprint_( blue )
{}

/// @brief
SheetConstraintsRCG::~SheetConstraintsRCG() {}

void
SheetConstraintsRCG::parse_my_tag( TagCOP const tag,
												basic::datacache::DataMap & data,
												protocols::filters::Filters_map const & filters,
												protocols::moves::Movers_map const & movers,
												core::pose::Pose const & pose )
{
	RemodelConstraintGenerator::parse_my_tag( tag, data, filters, movers, pose );
	// for all of these, the default is to not change what is already  present in the class
	set_blueprint( tag->getOption<std::string>( "blueprint", "" ) );
	set_distance( tag->getOption< core::Real >( "dist", dist_ ) );
	set_angle_tolerance( tag->getOption< core::Real >( "angle_tolerance", angle_tolerance_ ) );
	set_cacb_dihedral_tolerance( tag->getOption< core::Real >( "cacb_dihedral_tolerance", cacb_dihedral_tolerance_ ) );
	set_bb_dihedral_tolerance( tag->getOption< core::Real >( "bb_dihedral_tolerance", bb_dihedral_tolerance_ ) );
	set_weight( tag->getOption< core::Real >( "weight", weight_ ) );
	set_constrain_dist_only( tag->getOption< bool >( "constrain_dist_only", constrain_dist_only_ ) );
}

std::string
SheetConstraintsRCG::get_name() const
{
	return SheetCstGeneratorCreator::mover_name();
}

protocols::moves::MoverOP
SheetConstraintsRCG::fresh_instance() const
{
	return new SheetConstraintsRCG();
}

protocols::moves::MoverOP
SheetConstraintsRCG::clone() const
{
	return new SheetConstraintsRCG( *this );
}


/// @brief sets teh blueprint file used for determining proper sheet pairing
/// This function will create the blueprint object for you from the file and use it
void
SheetConstraintsRCG::set_blueprint( std::string const & blueprint_file )
{
 	if ( blueprint_file == "" ) {
		utility_exit_with_message( "SheetCstGenerator requires a blueprint file" );
  }
  set_blueprint( new protocols::jd2::parser::BluePrint( blueprint_file ) );
	if ( ! blueprint_ ) {
		utility_exit_with_message( "SheetCstGenerator tried to read a blueprint file, but failed to create the proper object." );
	}
}

/// @brief
void
SheetConstraintsRCG::set_blueprint( BluePrintOP const & blue )
{
	blueprint_ = blue;
}

/// @brief set weight
void
SheetConstraintsRCG::set_weight( Real const coef )
{
	weight_ = coef;
}

/// @brief set distance of constraint
void
SheetConstraintsRCG::set_distance( Real const dist )
{
	dist_ = dist;
}


/// @brief set the flat-bottom tolerance for the backbone angle between strands for each pair
/// This is N1-C1-C2 and N2-C2-C1 for parallel sheets, and N1-C1-N2/N2-C2-N1 for antiparallel.
void
SheetConstraintsRCG::set_angle_tolerance( Real const angle_tolerance )
{
	angle_tolerance_ = angle_tolerance;
}

/// @brief set the flat-bottom tolerance for the Cb1-Ca1-Ca2-Cb2 dihedral angle (0 = optimal)
void
SheetConstraintsRCG::set_cacb_dihedral_tolerance( Real const dihedral_tolerance )
{
	cacb_dihedral_tolerance_ = dihedral_tolerance;
}

/// @brief set the flat-bottom tolerance for the backbone dihedrals (0=optimal)
/// Dihedral 1 = O1-N1-C1-C2, Dihedral 2 = O2-N2-C2-C1
void
SheetConstraintsRCG::set_bb_dihedral_tolerance( Real const dihedral_tolerance )
{
	bb_dihedral_tolerance_ = dihedral_tolerance;
}

/// @brief sets whether we should constrain distance only, and not generate dihedral and angle constraints
void
SheetConstraintsRCG::set_constrain_dist_only( bool const constrain_dist_only )
{
	constrain_dist_only_ = constrain_dist_only;
}

/// @brief
void
SheetConstraintsRCG::generate_remodel_constraints( Pose const & pose )
{
	using core::scoring::constraints::AtomPairConstraint;
	using core::scoring::constraints::BoundFunc;
	using core::scoring::constraints::ConstraintOPs;
	using core::scoring::constraints::OffsetPeriodicBoundFunc;
	using core::scoring::constraints::ScalarWeightedFunc;
	using core::scoring::constraints::ScalarWeightedFuncOP;
	using protocols::fldsgn::topology::StrandPairing;
	using protocols::fldsgn::topology::StrandPairings;
	using protocols::fldsgn::topology::StrandPairingOP;
	using protocols::fldsgn::topology::StrandPairingSet;
	using protocols::fldsgn::topology::SS_Info2;
	using protocols::fldsgn::topology::SS_Info2_OP;

	// returned data
	ConstraintOPs csts;

	// set constraint func
	core::Real lb( 0.0 );
	core::Real ub( dist_ );
	core::Real sd( 1.0 );
	std::string tag( "constraints_in_beta_sheet" );
	ScalarWeightedFuncOP cstfunc = new ScalarWeightedFunc( weight_, new BoundFunc( lb, ub, sd, tag ) );
	// TL Oct '12 add weights to angle/dihedral constraints
	ScalarWeightedFuncOP cacb_dihedral_func = new ScalarWeightedFunc( weight_, new OffsetPeriodicBoundFunc(-cacb_dihedral_tolerance_,cacb_dihedral_tolerance_, std::sqrt(1.0/42.0), "dihed_cacb", 6.28, 0.0 ) );
	ScalarWeightedFuncOP bb_dihedral_func = new ScalarWeightedFunc( weight_, new OffsetPeriodicBoundFunc(-bb_dihedral_tolerance_,bb_dihedral_tolerance_, std::sqrt(1.0/42.0), "dihed_bb", 3.14, 0.0 ) );
	ScalarWeightedFuncOP bb_angle_func = new ScalarWeightedFunc( weight_, new BoundFunc(1.57-angle_tolerance_,1.57+angle_tolerance_, sqrt(1.0/42.0), "angle_bb") );

	// set constraints to csts
	core::Size nres( pose.total_residue() );
	//flo sep '12 in case we have ligands in the pose, don't count them
	for( core::Size i = nres; i != 0; i-- ){
		if( pose.residue_type(i).is_ligand() ) {
			nres--;
			TR << pose.residue( i ).name3() << i << " = ligand" << std::endl;
		}	
	}
	runtime_assert( blueprint_ );
	TR << "Blueprint num res=" << blueprint_->total_residue() << std::endl;
	TR << "Pose nres protein=" << nres << std::endl;
	runtime_assert( nres == blueprint_->total_residue() );

	TR << "Blueprint file is used for determining constrained residue pairs.  " << std::endl;
	TR << "Constrains between CA-CA atoms in sheet are applied for the following residues. " << std::endl;
	TR << "dist=" << dist_ << ", weight_=" << weight_ << ", cacb_dihedral=" << cacb_dihedral_tolerance_ << ", bb_dihedral=" << bb_dihedral_tolerance_ << ", angle=" << angle_tolerance_ << std::endl;

	SS_Info2_OP ssinfo = new SS_Info2( pose, blueprint_->secstruct() );
	StrandPairingSet spairset( blueprint_->strand_pairings(), ssinfo );
	StrandPairings spairs = spairset.strand_pairings();
	for ( utility::vector1< StrandPairingOP >::const_iterator it=spairs.begin(); it!=spairs.end(); ++it ) {

		StrandPairing spair=**it;
		for( core::Size iaa=spair.begin1(); iaa<=spair.end1(); iaa++ ) {
			core::Size jaa( spair.residue_pair( iaa ) );
			core::id::AtomID atom1( pose.residue_type( iaa ).atom_index( "CA" ), iaa );
			core::id::AtomID atom2( pose.residue_type( jaa ).atom_index( "CA" ), jaa );
			csts.push_back( new AtomPairConstraint( atom1, atom2, cstfunc ) );
			//flo sep '12: constrain dihedral, might be more accurate
			if( ( !constrain_dist_only_ ) ||
					( basic::options::option[ basic::options::OptionKeys::flxbb::constraints_sheet_include_cacb_pseudotorsion ].value() ) ){
				core::id::AtomID resi_n( pose.residue_type( iaa ).atom_index( "N" ), iaa );
				core::id::AtomID resi_c( pose.residue_type( iaa ).atom_index( "C" ), iaa );
				core::id::AtomID resi_o( pose.residue_type( iaa ).atom_index( "O" ), iaa );
        core::id::AtomID resj_n( pose.residue_type( jaa ).atom_index( "N" ), jaa );
        core::id::AtomID resj_c( pose.residue_type( jaa ).atom_index( "C" ), jaa );
        core::id::AtomID resj_o( pose.residue_type( jaa ).atom_index( "O" ), jaa );
				csts.push_back( new core::scoring::constraints::DihedralConstraint( resi_o, resi_n, resi_c, resj_c, bb_dihedral_func ) );
				csts.push_back( new core::scoring::constraints::DihedralConstraint( resj_o, resj_n, resj_c, resi_c, bb_dihedral_func ) );
				if( spair.orient() == 'P' ){
					csts.push_back( new core::scoring::constraints::AngleConstraint( resi_n, resi_c, resj_c, bb_angle_func ) );
					csts.push_back( new core::scoring::constraints::AngleConstraint( resj_n, resj_c, resi_c, bb_angle_func ) );
				}
				else if( spair.orient() == 'A' ){
          csts.push_back( new core::scoring::constraints::AngleConstraint( resi_n, resi_c, resj_n, bb_angle_func ) );
          csts.push_back( new core::scoring::constraints::AngleConstraint( resj_n, resj_c, resi_n, bb_angle_func ) );
				}
				if( (pose.residue_type( iaa ).name3() == "GLY") || (pose.residue_type( jaa ).name3() == "GLY" ) ) continue; // don't bother restraining cacb dihedral with gly
				core::id::AtomID resi_cb( pose.residue_type( iaa ).atom_index( "CB" ), iaa );
      	core::id::AtomID resj_cb( pose.residue_type( jaa ).atom_index( "CB" ), jaa );
				csts.push_back( new core::scoring::constraints::DihedralConstraint( resi_cb, atom1, atom2, resj_cb, cacb_dihedral_func ) );
				TR << "Added dihedral constraint between residues " << iaa << " and " << jaa << std::endl;
			}
			// flo sep '12 over
		} // for( Size i=1 )

	} // StrandPairingOP
	this->add_constraints( csts );
} //generate constraints


} //namespace fldsgn
} //namespace protocols
