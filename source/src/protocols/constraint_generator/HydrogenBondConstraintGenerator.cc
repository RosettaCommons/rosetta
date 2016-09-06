// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/constraints/MultiConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <core/scoring/func/FlatHarmonicFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/select/residue_selector/TrueResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/stream_util.hh>
#include <utility/tag/Tag.hh>

// Numeric headers
#include <numeric/constants.hh>

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
	return HydrogenBondConstraintGenerator::class_name();
}

HydrogenBondConstraintGenerator::HydrogenBondConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( HydrogenBondConstraintGenerator::class_name() ),
	selector1_( new core::select::residue_selector::TrueResidueSelector ),
	selector2_( new core::select::residue_selector::TrueResidueSelector ),
	atoms1_(),
	atoms2_(),
	atom_pair_func_(),
	angle1_func_(),
	angle2_func_(),
	bounded_( false ),
	atom_pair_sd_( 0.5 ),
	angle_sd_( 0.5 )
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

	std::string const hbond_def = tag->getOption< std::string >( "atom_definitions", "" );
	if ( !hbond_def.empty() ) add_atom_definitions( hbond_def );

	std::string const & func_def = tag->getOption< std::string >( "atom_pair_func", "" );
	if ( !func_def.empty() ) set_atom_pair_func( func_def );

	std::string const & ang1_func_def = tag->getOption< std::string >( "angle1_func", "" );
	if ( !ang1_func_def.empty() ) set_angle1_func( ang1_func_def );

	std::string const & ang2_func_def = tag->getOption< std::string >( "angle2_func", "" );
	if ( !ang2_func_def.empty() ) set_angle2_func( ang2_func_def );

	set_bounded( tag->getOption< bool >( "bounded", bounded_ ) );
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
	if ( atoms_str.empty() ) {
		set_atoms2( AtomNameSet() );
	} else {
		utility::vector1< std::string > const atom_strs = utility::string_split( atoms_str, ',' );
		set_atoms1( AtomNameSet( atom_strs.begin(), atom_strs.end() ) );
	}
}

void
HydrogenBondConstraintGenerator::set_atoms1( AtomNameSet const & atoms )
{
	atoms1_ = atoms;
}

void
HydrogenBondConstraintGenerator::set_atoms2( std::string const & atoms_str )
{
	if ( atoms_str.empty() ) {
		set_atoms2( AtomNameSet() );
	} else {
		utility::vector1< std::string > const atom_strs = utility::string_split( atoms_str, ',' );
		set_atoms2( AtomNameSet( atom_strs.begin(), atom_strs.end() ) );
	}
}

void
HydrogenBondConstraintGenerator::set_atoms2( AtomNameSet const & atoms )
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

void
HydrogenBondConstraintGenerator::set_bounded( bool const bounded )
{
	bounded_ = bounded;
}

void
HydrogenBondConstraintGenerator::set_atom_pair_sd( core::Real const sd )
{
	atom_pair_sd_ = sd;
}

void
HydrogenBondConstraintGenerator::set_angle_sd( core::Real const sd )
{
	angle_sd_ = sd;
}

void
HydrogenBondConstraintGenerator::add_atom_definitions( std::string const & atom_def )
{
	HydrogenBondInfo & info = *HydrogenBondInfo::get_instance();

	utility::vector1< std::string > const definitions = utility::string_split( atom_def, ';' );
	for ( utility::vector1< std::string >::const_iterator def=definitions.begin(); def!=definitions.end(); ++def ) {
		TR.Debug << "Adding atom definition from string: " << *def << std::endl;
		TR.Debug << "Added " << info.add_atoms_from_string( *def ) << std::endl;
	}
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

		for ( core::Size resid2=1; resid2<=pose.total_residue(); ++resid2 ) {
			if ( !subset2[ resid2 ] ) continue;
			if ( resid1 == resid2 ) continue;

			core::scoring::constraints::ConstraintOP cst =
				create_residue_constraint( pose.residue( resid1 ), pose.residue( resid2 ) );

			if ( !cst ) {
				TR.Debug << "No valid atoms found to generator Hbond constraints between residues " << resid1
					<< " and " << resid2 << " - skipping" <<  std::endl;
				continue;
			}

			csts.push_back( cst );
		}
	}

	return csts;
}

bool
can_hydrogen_bond( std::string const & element_string )
{
	if ( element_string == "O" ) return true;
	if ( element_string == "N" ) return true;
	if ( element_string == "S" ) return true;
	return false;
}

void
HydrogenBondConstraintGenerator::compute_valid_atoms(
	HydrogenBondingAtoms & ,
	core::conformation::Residue const & rsd ) const
{
	core::chemical::AtomIndices const & sc_atoms = rsd.type().all_sc_atoms();
	for ( core::chemical::AtomIndices::const_iterator at=sc_atoms.begin(); at!=sc_atoms.end(); ++at ) {
		if ( ! can_hydrogen_bond( rsd.type().atom_type( *at ).element() ) ) {
			continue;
		}
	}
}

HydrogenBondingAtoms
HydrogenBondConstraintGenerator::choose_atoms(
	core::conformation::Residue const & rsd,
	AtomNameSet const & allowed_atoms ) const
{
	HydrogenBondingAtoms const hb_atoms = HydrogenBondInfo::get_instance()->atoms( rsd.name3() );
	HydrogenBondingAtoms valid_atoms;

	if ( hb_atoms.empty() ) {
		compute_valid_atoms( valid_atoms, rsd );
		TR.Warning << "No hydrogen bonding atoms found in residue " << rsd.name3()
			<< rsd.seqpos() << " - not generating any constraints for this residue."
			<< " Full residue name is " << rsd.name3() << std::endl;
	} else {
		// if atoms are found for this residue, look for existing ones in the allowed set
		for ( HydrogenBondingAtoms::const_iterator a=hb_atoms.begin(); a!=hb_atoms.end(); ++a ) {
			if ( ! rsd.has( a->hb_atom() ) ) {
				TR << "Skipping HBond atom " << a->hb_atom() << " because it does not exist in residue "
					<< rsd.name() << rsd.seqpos() << std::endl;
				continue;
			}
			if ( ! rsd.has( a->atom2() ) ) {
				TR << "Skipping HBond atom " << a->hb_atom() << " because atom " << a->atom2()
					<< " does not exist in residue " << rsd.name() << rsd.seqpos() << std::endl;
				continue;
			}
			if ( ! rsd.has( a->atom3() ) ) {
				TR << "Skipping HBond atom " << a->hb_atom() << " because atom " << a->atom3()
					<< " does not exist in residue " << rsd.name() << rsd.seqpos() << std::endl;
				continue;
			}
			// skip if allowed_atoms has things, but doesn't have this atom
			if ( ( !allowed_atoms.empty() ) && ( allowed_atoms.find( a->hb_atom() ) == allowed_atoms.end() ) ) {
				TR.Debug << "Skipping HBond atom " << a->hb_atom() << " because it is not in the set of allowed atoms" << std::endl
					<< "Allowed atoms = [";
				for ( AtomNameSet::const_iterator a=allowed_atoms.begin(); a!=allowed_atoms.end(); ++a ) {
					TR.Debug << " " << *a;
				}
				TR.Debug << " ], empty=" << allowed_atoms.empty() << std::endl;
				continue;
			}
			valid_atoms.push_back( *a );
		}
	}

	return valid_atoms;
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2 ) const
{
	return create_residue_constraint( rsd1, rsd2, choose_atoms( rsd1, atoms1_ ), choose_atoms( rsd2, atoms2_ ) );
}

core::scoring::constraints::ConstraintOP
ambiguous_constraint_wrap( core::scoring::constraints::ConstraintOPs const & csts )
{
	if ( csts.size() == 1 ) {
		// no sense to create ambiguous constraint for 1 cst
		return *( csts.begin() );
	} else {
		return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::AmbiguousConstraint( csts ) );
	}
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::conformation::Residue const & rsd1,
	core::conformation::Residue const & rsd2,
	HydrogenBondingAtoms const & atoms1,
	HydrogenBondingAtoms const & atoms2 ) const
{
	core::scoring::constraints::ConstraintOPs csts;
	TR.Debug << "Going to create constraints for res1=" << rsd1.seqpos() << " atoms1=" << atoms1
		<< " res2=" << rsd2.seqpos() << " atoms2=" << atoms2 << std::endl;
	for ( HydrogenBondingAtoms::const_iterator a1=atoms1.begin(); a1!=atoms1.end(); ++a1 ) {
		core::id::AtomID const atomid1( rsd1.type().atom_index( a1->hb_atom() ), rsd1.seqpos() );
		core::id::AtomID const parent_atomid1( rsd1.type().atom_index( a1->atom2() ), rsd1.seqpos() );
		core::id::AtomID const parent2_atomid1( rsd1.type().atom_index( a1->atom3() ), rsd1.seqpos() );
		for ( HydrogenBondingAtoms::const_iterator a2=atoms2.begin(); a2!=atoms2.end(); ++a2 ) {
			core::id::AtomID const atomid2( rsd2.type().atom_index( a2->hb_atom() ), rsd2.seqpos() );
			core::id::AtomID const parent_atomid2( rsd2.type().atom_index( a2->atom2() ), rsd2.seqpos() );
			core::id::AtomID const parent2_atomid2( rsd2.type().atom_index( a2->atom3() ), rsd2.seqpos() );
			csts.push_back( create_residue_constraint(
				atomid1, parent_atomid1, parent2_atomid1,
				atomid2, parent_atomid2, parent2_atomid2, *a1, *a2 ) );
		}
	}
	if ( csts.empty() ) {
		return core::scoring::constraints::ConstraintOP();
	}
	return ambiguous_constraint_wrap( csts );
}

core::scoring::func::FuncOP
HydrogenBondConstraintGenerator::atom_pair_func(
	HydrogenBondingAtom const & a1,
	HydrogenBondingAtom const & a2 ) const
{
	if ( atom_pair_func_ ) {
		return atom_pair_func_;
	}
	if ( bounded_ ) {
		return core::scoring::func::FuncOP(
			new core::scoring::func::FlatHarmonicFunc( a1.distance() + a2.distance(), atom_pair_sd_, 1.5 ) );
	} else {
		return core::scoring::func::FuncOP(
			new core::scoring::func::HarmonicFunc( a1.distance() + a2.distance(), atom_pair_sd_ ) );
	}
}

core::scoring::func::FuncOP
HydrogenBondConstraintGenerator::angle1_func( HydrogenBondingAtom const & a1 ) const
{
	if ( angle1_func_ ) {
		return angle1_func_;
	}
	if ( bounded_ ) {
		return core::scoring::func::FuncOP(
			new core::scoring::func::FlatHarmonicFunc( a1.angle(), angle_sd_, 0.4 ) );
	} else {
		return core::scoring::func::FuncOP(
			new core::scoring::func::CircularHarmonicFunc( a1.angle(), angle_sd_ ) );
	}
}

core::scoring::func::FuncOP
HydrogenBondConstraintGenerator::angle2_func( HydrogenBondingAtom const & a2 ) const
{
	if ( angle2_func_ ) {
		return angle2_func_;
	}
	if ( bounded_ ) {
		return core::scoring::func::FuncOP(
			new core::scoring::func::FlatHarmonicFunc( a2.angle(), angle_sd_, 0.4 ) );
	} else {
		return core::scoring::func::FuncOP(
			new core::scoring::func::CircularHarmonicFunc( a2.angle(), angle_sd_ ) );
	}
}

core::scoring::func::FuncOP
HydrogenBondConstraintGenerator::dihedral_func( core::Real const dihedral_value ) const
{
	if ( bounded_ ) {
		return core::scoring::func::FuncOP(
			new core::scoring::func::FlatHarmonicFunc( dihedral_value, angle_sd_, 0.4 ) );
	} else {
		return core::scoring::func::FuncOP(
			new core::scoring::func::CircularHarmonicFunc( dihedral_value, angle_sd_ ) );
	}
}

core::scoring::constraints::ConstraintOP
HydrogenBondConstraintGenerator::create_residue_constraint(
	core::id::AtomID const & atomid1,
	core::id::AtomID const & parent_atomid1,
	core::id::AtomID const & parent2_atomid1,
	core::id::AtomID const & atomid2,
	core::id::AtomID const & parent_atomid2,
	core::id::AtomID const & parent2_atomid2,
	HydrogenBondingAtom const & a1,
	HydrogenBondingAtom const & a2 ) const
{
	using namespace core::scoring::constraints;

	TR.Debug << "Creating cst: " << atomid1 << " " << parent_atomid1 << " "
		<< atomid2 << " " << parent_atomid2 << std::endl;

	ConstraintOPs hb_csts;
	hb_csts.push_back( ConstraintOP( new AtomPairConstraint( atomid1, atomid2, atom_pair_func( a1, a2 ) ) ) );
	hb_csts.push_back( ConstraintOP( new AngleConstraint( parent_atomid1, atomid1, atomid2, angle1_func( a1 ) ) ) );
	hb_csts.push_back( ConstraintOP( new AngleConstraint( atomid1, atomid2, parent_atomid2, angle2_func( a2 ) ) ) );

	if ( !a1.dihedrals().empty() ) {
		ConstraintOPs dihedral_csts;
		for ( HydrogenBondingAtom::Dihedrals::const_iterator d=a1.dihedrals().begin(); d!=a1.dihedrals().end(); ++d ) {
			dihedral_csts.push_back( ConstraintOP( new DihedralConstraint(
				parent2_atomid1, parent_atomid1, atomid1, atomid2, dihedral_func( *d ) ) ) );
		}
		hb_csts.push_back( ambiguous_constraint_wrap( dihedral_csts ) );
	}
	if ( !a2.dihedrals().empty() ) {
		ConstraintOPs dihedral_csts;
		for ( HydrogenBondingAtom::Dihedrals::const_iterator d=a2.dihedrals().begin(); d!=a2.dihedrals().end(); ++d ) {
			dihedral_csts.push_back( ConstraintOP( new DihedralConstraint(
				atomid1, atomid2, parent_atomid2, parent2_atomid2, dihedral_func( *d ) ) ) );
		}
		hb_csts.push_back( ambiguous_constraint_wrap( dihedral_csts ) );
	}
	return ConstraintOP( new MultiConstraint( hb_csts ) );
}

HydrogenBondingAtom::HydrogenBondingAtom(
	std::string const & atom1,
	std::string const & atom2,
	std::string const & atom3,
	core::Real const ideal_distance,
	core::Real const ideal_angle,
	Dihedrals const & ideal_dihedral ):
	utility::pointer::ReferenceCount(),
	atom_( atom1 ),
	atom2_( atom2 ),
	atom3_( atom3 ),
	distance_( ideal_distance ),
	angle_( numeric::constants::f::pi / 180.0 * ideal_angle ),
	dihedrals_()
{
	for ( Dihedrals::const_iterator d=ideal_dihedral.begin(); d!=ideal_dihedral.end(); ++d ) {
		dihedrals_.push_back( numeric::constants::f::pi / 180.0 * *d );
	}
}

std::ostream &
operator<<( std::ostream & os, HydrogenBondingAtom const & atom )
{
	os << atom.atom_ << ", " << atom.atom2_ << ", " << atom.atom3_ << ", "
		<< atom.distance_ << ", " << atom.angle_ << ", " << atom.dihedrals_;
	return os;
}


HydrogenBondInfo *
HydrogenBondInfo::create_singleton_instance()
{
	return new HydrogenBondInfo;
}

HydrogenBondInfo::HydrogenBondInfo():
	utility::SingletonBase< HydrogenBondInfo >(),
	atoms_()
{
	HydrogenBondingAtom::Dihedrals const zero_and_180 = boost::assign::list_of (0.0) (180.0);
	HydrogenBondingAtom::Dihedrals const extended = boost::assign::list_of (180.0);
	HydrogenBondingAtom::Dihedrals const none;

	HydrogenBondingAtoms & arg = create_residue( "ARG" );
	arg.push_back( HydrogenBondingAtom( "NE", "CZ", "NH1", 1.4, 120.0, zero_and_180 ) );
	arg.push_back( HydrogenBondingAtom( "NH1", "CZ", "NE", 1.4, 120.0, zero_and_180 ) );
	arg.push_back( HydrogenBondingAtom( "NH2", "CZ", "NE", 1.4, 120.0, zero_and_180 ) );

	HydrogenBondingAtoms & asn = create_residue( "ASN" );
	asn.push_back( HydrogenBondingAtom( "OD1", "CG", "ND2", 1.4, 120.0, zero_and_180 ) );
	asn.push_back( HydrogenBondingAtom( "ND2", "CG", "OD1", 1.4, 120.0, zero_and_180 ) );

	HydrogenBondingAtoms & asp = create_residue( "ASP" );
	asp.push_back( HydrogenBondingAtom( "OD1", "CG", "OD2", 1.4, 120.0, zero_and_180 ) );
	asp.push_back( HydrogenBondingAtom( "OD2", "CG", "OD1", 1.4, 120.0, zero_and_180 ) );

	HydrogenBondingAtoms & cys = create_residue( "CYS" );
	cys.push_back( HydrogenBondingAtom( "SG", "CB", "CA", 2.0, 109.5, none ) );

	HydrogenBondingAtoms & gln = create_residue( "GLN" );
	gln.push_back( HydrogenBondingAtom( "OE1", "CD", "OE2", 1.4, 120.0, zero_and_180 ) );
	gln.push_back( HydrogenBondingAtom( "OE2", "CD", "OE1", 1.4, 120.0, zero_and_180 ) );

	HydrogenBondingAtoms & glu = create_residue( "GLU" );
	glu.push_back( HydrogenBondingAtom( "OE1", "CD", "OE2", 1.4, 120.0, zero_and_180 ) );
	glu.push_back( HydrogenBondingAtom( "OE2", "CD", "OE1", 1.4, 120.0, zero_and_180 ) );

	HydrogenBondingAtoms & his = create_residue( "HIS" );
	his.push_back( HydrogenBondingAtom( "ND1", "CE1", "NE2", 1.4, 120.0, extended ) );
	his.push_back( HydrogenBondingAtom( "NE2", "CE1", "ND1", 1.4, 120.0, extended ) );

	HydrogenBondingAtoms & his_d = create_residue( "HIS_D" );
	his_d.push_back( HydrogenBondingAtom( "ND1", "CE1", "NE2", 1.4, 120.0, extended ) );
	his_d.push_back( HydrogenBondingAtom( "NE2", "CE1", "ND1", 1.4, 120.0, extended ) );

	HydrogenBondingAtoms & lys = create_residue( "LYS" );
	lys.push_back( HydrogenBondingAtom( "NZ", "CE", "CD", 1.4, 109.5, none ) );

	HydrogenBondingAtoms & ser = create_residue( "SER" );
	ser.push_back( HydrogenBondingAtom( "OG", "CB", "CA", 1.4, 109.5, none ) );

	HydrogenBondingAtoms & thr = create_residue( "THR" );
	thr.push_back( HydrogenBondingAtom( "OG1", "CB", "CA", 1.4, 109.5, none ) );

	HydrogenBondingAtoms & trp = create_residue( "TRP" );
	trp.push_back( HydrogenBondingAtom( "NE1", "CD1", "CG", 1.4, 120.0, extended ) );

	HydrogenBondingAtoms & tyr = create_residue( "TYR" );
	tyr.push_back( HydrogenBondingAtom( "OH", "CZ", "CE1", 1.4, 109.5, none ) );
}

HydrogenBondingAtoms const HydrogenBondInfo::empty_;

HydrogenBondingAtoms &
HydrogenBondInfo::create_residue( std::string const & rsd_name )
{
	atoms_[ rsd_name ] = HydrogenBondingAtoms();
	return atoms_[ rsd_name ];
}

HydrogenBondingAtoms &
HydrogenBondInfo::retrieve_residue( std::string const & rsd_name )
{
	AtomNameMap::iterator hb_atoms = atoms_.find( rsd_name );
	if ( hb_atoms == atoms_.end() ) {
		hb_atoms = atoms_.insert( std::make_pair( rsd_name, HydrogenBondingAtoms() ) ).first;
	}
	debug_assert( hb_atoms != atoms_.end() );
	return hb_atoms->second;
}

HydrogenBondingAtoms &
HydrogenBondInfo::add_atoms_from_string( std::string const & description_str )
{
	utility::vector1< std::string > const fields = utility::string_split( description_str, ',' );
	if ( fields.size() < 6 ) {
		std::stringstream msg;
		msg << "HydrogenBondingAtoms::from_string(): Bad description string, not enough fields ("
			<< fields.size() << " fields found, >=7 fields required)!" << std::endl;
		msg << "Description string: " << description_str << std::endl;
		utility_exit_with_message( msg.str() );
	}

	utility::vector1< std::string >::const_iterator cur_field = fields.begin();
	std::string const & res_name = *cur_field;
	++cur_field;
	std::string const & atom1 = *cur_field;
	++cur_field;
	std::string const & atom2 = *cur_field;
	++cur_field;
	std::string const & atom3 = *cur_field;
	++cur_field;
	core::Real const dist = boost::lexical_cast< core::Real >( *cur_field );
	++cur_field;
	core::Real const angle = boost::lexical_cast< core::Real >( *cur_field );
	++cur_field;
	HydrogenBondingAtom::Dihedrals dihedrals;
	for ( ; cur_field!=fields.end(); ++cur_field ) {
		dihedrals.push_back( boost::lexical_cast< core::Real >( *cur_field ) );
	}

	HydrogenBondingAtoms & hb_atoms = retrieve_residue( res_name );
	hb_atoms.push_back( HydrogenBondingAtom( atom1, atom2, atom3, dist, angle, dihedrals ) );
	return hb_atoms;
}

HydrogenBondingAtoms const &
HydrogenBondInfo::atoms( std::string const & rsd_name ) const
{
	AtomNameMap::const_iterator a = atoms_.find( rsd_name );
	if ( a == atoms_.end() ) return empty_;
	else return a->second;
}

} // namespace constraint_generator
} // namespace protocols

namespace utility {
using namespace protocols::constraint_generator;
#ifdef MULTI_THREADED
template<> std::mutex SingletonBase< HydrogenBondInfo >::singleton_mutex_{};
template<> std::atomic< HydrogenBondInfo * > SingletonBase< HydrogenBondInfo >::instance_( NULL );
#else
template<> HydrogenBondInfo * SingletonBase< HydrogenBondInfo >::instance_ = NULL;
#endif
}

