// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraint_generator/denovo_design/HydrogenBondConstraintGenerator.cc
/// @brief
/// @author Tom Linsky ( tlinsky at uw dot edu )

// Unit headers
#include <protocols/constraint_generator/HydrogenBondConstraintGenerator.hh>
#include <protocols/constraint_generator/HydrogenBondConstraintGeneratorCreator.hh>

// Protocol headers
#include <protocols/constraint_generator/util.hh>

// Core headers
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/select/residue_selector/NeighborhoodResidueSelector.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Numeric headers

// Boost headers
#include <boost/assign.hpp>

// C++ headers

namespace protocols {
namespace constraint_generator {

static THREAD_LOCAL basic::Tracer TR( "protocols.constraint_generator.HydrogenBondConstraintGenerator" );

protocols::constraint_generator::ConstraintGeneratorOP
HydrogenBondConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new HydrogenBondConstraintGenerator );
}

std::string
HydrogenBondConstraintGeneratorCreator::keyname() const
{
	return HydrogenBondConstraintGeneratorCreator::constraint_generator_name() ;
}

std::string
HydrogenBondConstraintGeneratorCreator::constraint_generator_name()
{
	return "HydrogenBondConstraintGenerator";
}

HydrogenBondConstraintGenerator::HydrogenBondConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( HydrogenBondConstraintGeneratorCreator::constraint_generator_name() ),
	selector1_( new core::select::residue_selector::TrueResidueSelector ),
	selector2_( new core::select::residue_selector::TrueResidueSelector ),
	atoms1_(),
	atoms2_(),
	atom_pair_func_( new core::scoring::func::FlatHarmonicFunc( 2.0, 0.5, 1.5 ) ),
	angle1_func_( new core::scoring::func::FlatHarmonicFunc( 2.0, 0.5, 0.4 ) ),
	angle2_func_( new core::scoring::func::FlatHarmonicFunc( 2.0, 0.5, 0.4 ) )
{
}

HydrogenBondConstraintGenerator::~HydrogenBondConstraintGenerator()
{
}

protocols::constraint_generator::ConstraintGeneratorOP
HydrogenBondConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new HydrogenBondConstraintGenerator( *this ) );
}

void
HydrogenBondConstraintGenerator::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & data )
{
	std::string const & selector1_name = tag->getOption< std::string >( "residue_selector1", "" );
	if ( !selector1_name.empty() ) {
		core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::get_residue_selector( selector1_name, data );
		if ( !selector ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "'residue_selector1' residue selector was not found in the data map!\n" );
		}
		set_residue_selector1( selector );
	}

	std::string const & selector2_name = tag->getOption< std::string >( "residue_selector2", "" );
	if ( !selector2_name.empty() ) {
		core::select::residue_selector::ResidueSelectorCOP selector = core::select::residue_selector::get_residue_selector( selector2_name, data );
		if ( !selector ) {
			throw utility::excn::EXCN_RosettaScriptsOption( "'residue_selector2' residue selector was not found in the data map!\n" );
		}
		set_residue_selector2( selector );
	}

	set_atoms1( tag->getOption< std::string >( "atoms1", "" ) );
	set_atoms2( tag->getOption< std::string >( "atoms2", "" ) );

	std::string const & func_def = tag->getOption< std::string >( "atom_pair_func", "" );
	if ( !func_def.empty() ) set_atom_pair_func( func_def );

	std::string const & ang1_func_def = tag->getOption< std::string >( "angle1_func", "" );
	if ( !ang1_func_def.empty() ) set_angle1_func( ang1_func_def );

	std::string const & ang2_func_def = tag->getOption< std::string >( "angle2_func", "" );
	if ( !ang2_func_def.empty() ) set_angle2_func( ang2_func_def );

}

void
HydrogenBondConstraintGenerator::set_residue_selector1( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector1_ = selector;
}

void
HydrogenBondConstraintGenerator::set_residue_selector2( core::select::residue_selector::ResidueSelectorCOP selector )
{
	selector2_ = selector;
}

void
HydrogenBondConstraintGenerator::set_atoms1( std::string const & atoms_str )
{
	set_atoms1( utility::string_split( atoms_str, ',' ) );
}

void
HydrogenBondConstraintGenerator::set_atoms1( utility::vector1< std::string > const & atoms )
{
	atoms1_ = atoms;
}

void
HydrogenBondConstraintGenerator::set_atoms2( std::string const & atoms_str )
{
	set_atoms2( utility::string_split( atoms_str, ',' ) );
}

void
HydrogenBondConstraintGenerator::set_atoms2( utility::vector1< std::string > const & atoms )
{
	atoms2_ = atoms;
}

void
HydrogenBondConstraintGenerator::set_atom_pair_func( std::string const & func_def )
{
	atom_pair_func_ = create_func( func_def );
}

void
HydrogenBondConstraintGenerator::set_atom_pair_func( core::scoring::func::FuncOP atompairfunc )
{
	atom_pair_func_ = atompairfunc;
}

void
HydrogenBondConstraintGenerator::set_angle1_func( std::string const & func_def )
{
	angle1_func_ = create_func( func_def );
}

void
HydrogenBondConstraintGenerator::set_angle1_func( core::scoring::func::FuncOP angle_func )
{
	angle1_func_ = angle_func;
}

void
HydrogenBondConstraintGenerator::set_angle2_func( std::string const & func_def )
{
	angle2_func_ = create_func( func_def );
}

void
HydrogenBondConstraintGenerator::set_angle2_func( core::scoring::func::FuncOP angle_func )
{
	angle2_func_ = angle_func;
}

utility::vector1< core::Size >
existing_atoms( core::conformation::Residue const & rsd, utility::vector1< std::string > const & atoms )
{
	utility::vector1< core::Size > atom_idxs;
	for ( utility::vector1< std::string >::const_iterator aname=atoms.begin(); aname!=atoms.end(); ++aname ) {
		if ( rsd.has( *aname ) ) atom_idxs.push_back( rsd.type().atom_index( *aname ) );
	}
	return atom_idxs;
}

core::scoring::constraints::ConstraintCOPs
HydrogenBondConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	if ( !selector1_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "\"residue_selector1\" must be specified to HydrogenBondConstraintGenerator!" );
	}

	if ( !selector2_ ) {
		throw utility::excn::EXCN_RosettaScriptsOption( "\"residue_selector2\" must be specified to HydrogenBondConstraintGenerator!" );
	}

	debug_assert( selector1_ );
	debug_assert( selector2_ );
	core::select::residue_selector::ResidueSubset const subset1 = selector1_->apply( pose );
	core::select::residue_selector::ResidueSubset const subset2 = selector2_->apply( pose );

	core::scoring::constraints::ConstraintCOPs csts;

	for ( core::Size resid1=1; resid1<=pose.total_residue(); ++resid1 ) {
		if ( !subset1[ resid1 ] ) continue;
		if ( ! pose.residue( resid1 ).is_protein() ) continue;

		for ( core::Size resid2=1; resid2<=pose.total_residue(); ++resid2 ) {
			if ( !subset2[ resid2 ] ) continue;
			if ( ! pose.residue( resid2 ).is_protein() ) continue;
			if ( resid1 == resid2 ) continue;

			core::scoring::constraints::ConstraintOP cst =
				create_residue_constraint( pose.residue( resid1 ), pose.residue( resid2 ) );

			if ( !cst ) {
				TR.Debug << "NO valid atoms found to generator Hbond constraints between residues " << resid1
					<< " and " << resid2 << " - skipping" <<  std::endl;
				continue;
			}

			csts.push_back( cst );
		}
	}

	return csts;
}

utility::vector1< core::Size >
HydrogenBondConstraintGenerator::choose_atoms(
	core::conformation::Residue const & rsd,
	utility::vector1< std::string > const & atoms ) const
{
	// if no atoms are specified, look through the residue type and try to find some
	if ( atoms.empty() ) {
		utility::vector1< core::Size > atom_idxs;
		core::chemical::AtomIndices const & sc_atoms = rsd.type().all_sc_atoms();
		for ( core::chemical::AtomIndices::const_iterator at=sc_atoms.begin(); at!=sc_atoms.end(); ++at ) {
			if ( ( rsd.type().atom_type( *at ).element() == "O" ) ||
					( rsd.type().atom_type( *at ).element() == "N" ) ||
					( rsd.type().atom_type( *at ).element() == "S" ) ) {
				atom_idxs.push_back( *at );
			}
		}
		return atom_idxs;
	}
	return existing_atoms( rsd, atoms );
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2 ) const
{
	return create_residue_constraint( rsd1, rsd2, choose_atoms( rsd1, atoms1_ ), choose_atoms( rsd2, atoms2_ ) );
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	utility::vector1< core::Size > const & atoms1,
	utility::vector1< core::Size > const & atoms2 ) const
{
	core::scoring::constraints::ConstraintCOPs csts;
	for ( utility::vector1< core::Size >::const_iterator a1=atoms1.begin(); a1!=atoms1.end(); ++a1 ) {
		core::id::AtomID const parent_atomid1( rsd1.type().atom_base( *a1 ), rsd1.seqpos() );
		core::id::AtomID const atomid1( *a1, rsd1.seqpos() );
		for ( utility::vector1< core::Size >::const_iterator a2=atoms2.begin(); a2!=atoms2.end(); ++a2 ) {
			core::id::AtomID const parent_atomid2( rsd2.type().atom_base( *a2 ), rsd2.seqpos() );
			core::id::AtomID const atomid2( *a2, rsd2.seqpos() );
			csts.push_back( create_residue_constraint( atomid1, parent_atomid1, atomid2, parent_atomid2 ) );
		}
	}
	if ( csts.empty() ) {
		return core::scoring::constraints::ConstraintOP();
	} else {
		return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint( csts ) );
	}
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::id::AtomID const & atomid1,
	core::id::AtomID const & parent_atomid1,
	core::id::AtomID const & atomid2,
	core::id::AtomID const & parent_atomid2 ) const
{
	using namespace core::scoring::constraints;

	TR.Debug << "Creating cst: " << atomid1 << " " << parent_atomid1 << " "
		<< atomid2 << " " << parent_atomid2 << std::endl;

	ConstraintOPs hb_csts;
	hb_csts.push_back( ConstraintOP( new AtomPairConstraint( atomid1, atomid2, atom_pair_func_ ) ) );
	hb_csts.push_back( ConstraintOP( new AngleConstraint( parent_atomid1, atomid1, atomid2, angle1_func_ ) ) );
	hb_csts.push_back( ConstraintOP( new AngleConstraint( atomid1, atomid2, parent_atomid2, angle2_func_ ) ) );
	return ConstraintOP( new MultiConstraint( hb_csts ) );
}

} // namespace constraint_generator
} // namespace protocols
