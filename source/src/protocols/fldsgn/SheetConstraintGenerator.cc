// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file  protocols/forge/constraints/SheetConstraintGenerator.cc
///
/// @brief
/// @author Nobuyasu Koga( nobuyasu@uw.edu ) , October 2009
/// @modified Tom Linsky (tlinsky@uw.edu)

// Unit header
#include <protocols/fldsgn/SheetConstraintGenerator.hh>
#include <protocols/fldsgn/SheetConstraintGeneratorCreator.hh>

// Package headers
#include <protocols/constraint_generator/util.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <protocols/fldsgn/topology/StrandPairing.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueType.hh>
#include <core/scoring/constraints/AmbiguousConstraint.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/ScalarWeightedFunc.hh>
#include <core/scoring/func/HarmonicFunc.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Utility headers
#include <basic/options/option.hh>
#include <basic/options/keys/flxbb.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <numeric/constants.hh>
#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>

// C++ headers
#include <boost/assign.hpp>
#include <cmath>

static THREAD_LOCAL basic::Tracer TR( "protocols.fldsgn.SheetConstraintGenerator" );

namespace protocols {
namespace fldsgn {

std::string
SheetConstraintGeneratorCreator::keyname() const
{
	return SheetConstraintGeneratorCreator::constraint_generator_name();
}

protocols::constraint_generator::ConstraintGeneratorOP
SheetConstraintGeneratorCreator::create_constraint_generator() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new SheetConstraintGenerator );
}

std::string
SheetConstraintGeneratorCreator::constraint_generator_name()
{
	return "SheetConstraintGenerator";
}

/// @brief
SheetConstraintGenerator::SheetConstraintGenerator():
	protocols::constraint_generator::ConstraintGenerator( SheetConstraintGeneratorCreator::constraint_generator_name() ),
	weight_( 1.0 ),
	dist_( 5.5 ),
	dist_tolerance_( 1.0 ),
	angle_tolerance_( 0.35 ),
	cacb_dihedral_tolerance_( 0.9 ),
	bb_dihedral_tolerance_( 0.52 ),
	constrain_ca_ca_dist_( true ),
	constrain_bb_cacb_dihedral_( true ),
	constrain_bb_dihedral_( true ),
	constrain_bb_angle_( true ),
	flat_bottom_constraints_( true ),
	secstruct_( "" ),
	spairs_( "" )
{}

/// @brief
SheetConstraintGenerator::~SheetConstraintGenerator() {}

protocols::constraint_generator::ConstraintGeneratorOP
SheetConstraintGenerator::clone() const
{
	return protocols::constraint_generator::ConstraintGeneratorOP( new SheetConstraintGenerator( *this ) );
}

/// @brief
core::scoring::constraints::ConstraintCOPs
SheetConstraintGenerator::apply( core::pose::Pose const & pose ) const
{
	using core::scoring::func::FuncOP;

	// returned data
	core::scoring::constraints::ConstraintOPs csts;

	// set constraint func
	FuncOP caca_atom_pair_func = create_caca_atom_pair_func( dist_ );
	FuncOP cacb_dihedral_func = create_cacb_dihedral_func( 0.0 );
	FuncOP bb_shift_dihedral_func1 = create_bb_dihedral_func( 0.0 );
	FuncOP bb_shift_dihedral_func2 = create_bb_dihedral_func( numeric::constants::f::pi );
	FuncOP bb_angle_func = create_bb_angle_func( numeric::constants::f::pi / 2 );

	std::pair< std::string, std::string > ss_spair = get_secstruct_and_strandpairings( pose );
	std::string const & secstruct = ss_spair.first;
	std::string const & spair_str = ss_spair.second;

	TR << "Found SS = " << secstruct << std::endl;
	TR << "Found strand pairings " << spair_str << std::endl;
	TR << "Constrains between CA-CA atoms in sheet are applied for the following residues. " << std::endl;
	TR << "dist=" << dist_ << ", weight_=" << weight_ << ", cacb_dihedral=" << cacb_dihedral_tolerance_ << ", bb_dihedral=" << bb_dihedral_tolerance_ << ", angle=" << angle_tolerance_ << std::endl;

	protocols::fldsgn::topology::SS_Info2_OP ssinfo( new protocols::fldsgn::topology::SS_Info2( pose, secstruct ) );
	protocols::fldsgn::topology::StrandPairingSet spairset( spair_str, ssinfo );
	ResiduePairs const respairs = compute_residue_pairs( spairset.strand_pairings() );
	for ( ResiduePairs::const_iterator pair=respairs.begin(); pair!=respairs.end(); ++pair ) {
		core::Size const iaa = pair->first;
		core::Size const jaa = pair->second;

		// Residues probably need to be protein to get 'E' secondary structure.
		// However, these checks are fast and could prevent issues
		if ( !pose.residue( iaa ).is_protein() ) continue;
		if ( !pose.residue( jaa ).is_protein() ) continue;

		// atom definitions
		core::id::AtomID const atom1( pose.residue_type( iaa ).atom_index( "CA" ), iaa );
		core::id::AtomID const atom2( pose.residue_type( jaa ).atom_index( "CA" ), jaa );
		core::id::AtomID const resi_n( pose.residue_type( iaa ).atom_index( "N" ), iaa );
		core::id::AtomID const resi_c( pose.residue_type( iaa ).atom_index( "C" ), iaa );
		core::id::AtomID const resi_o( pose.residue_type( iaa ).atom_index( "O" ), iaa );
		core::id::AtomID const resj_n( pose.residue_type( jaa ).atom_index( "N" ), jaa );
		core::id::AtomID const resj_c( pose.residue_type( jaa ).atom_index( "C" ), jaa );
		core::id::AtomID const resj_o( pose.residue_type( jaa ).atom_index( "O" ), jaa );

		// distance
		if ( constrain_ca_ca_dist_ ) {
			csts.push_back( create_ca_ca_atom_pair_constraint( atom1, atom2, caca_atom_pair_func ) );
		}

		// bb_dihedral
		if ( constrain_bb_dihedral_ ) {
			csts.push_back( create_bb_dihedral_constraint( pose, iaa, jaa, bb_shift_dihedral_func1, bb_shift_dihedral_func2 ) );
			csts.push_back( create_bb_dihedral_constraint( pose, jaa, iaa, bb_shift_dihedral_func1, bb_shift_dihedral_func2 ) );
		}

		// angle
		if ( constrain_bb_angle_ ) {
			if ( pair->orientation == 'P' ) {
				csts.push_back( create_bb_angle_constraint( resi_n, resi_c, resj_c, bb_angle_func ) );
				csts.push_back( create_bb_angle_constraint( resj_n, resj_c, resi_c, bb_angle_func ) );
			} else if ( pair->orientation == 'A' ) {
				csts.push_back( create_bb_angle_constraint( resi_n, resi_c, resj_n, bb_angle_func ) );
				csts.push_back( create_bb_angle_constraint( resj_n, resj_c, resi_n, bb_angle_func ) );
			}
		}

		if ( !pose.residue_type( iaa ).has( "CB" ) || !pose.residue_type( jaa ).has( "CB" ) ) {
			continue; // don't bother restraining cacb dihedral when there is no "CB" (e.g. gly)
		}

		if ( constrain_bb_cacb_dihedral_ || basic::options::option[ basic::options::OptionKeys::flxbb::constraints_sheet_include_cacb_pseudotorsion ].value() ) {
			core::id::AtomID const resi_cb( pose.residue_type( iaa ).atom_index( "CB" ), iaa );
			core::id::AtomID const resj_cb( pose.residue_type( jaa ).atom_index( "CB" ), jaa );
			csts.push_back( create_bb_cacb_dihedral_constraint( resi_cb, atom1, atom2, resj_cb, cacb_dihedral_func ) );
			TR << "Added dihedral constraint between residues " << iaa << " and " << jaa << std::endl;
		}
	} // residue pairs
	return csts;
} //apply

void
SheetConstraintGenerator::parse_tag(
	utility::tag::TagCOP const tag,
	basic::datacache::DataMap & )
{
	// for all of these, the default is to not change what is already  present in the class
	if ( tag->hasOption( "blueprint" ) ) {
		initialize_from_blueprint( tag->getOption< std::string >( "blueprint" ) );
	}

	set_secstruct( tag->getOption< std::string >( "secstruct", secstruct_ ) );
	set_strand_pairs( tag->getOption< std::string >( "spairs", spairs_ ) );
	set_flat_bottom_constraints( tag->getOption< bool >( "flat_bottom_constraints", flat_bottom_constraints_ ) );
	set_distance( tag->getOption< core::Real >( "dist", dist_ ) );
	set_distance_tolerance( tag->getOption< core::Real >( "dist_tolerance", dist_tolerance_ ) );
	set_angle_tolerance( tag->getOption< core::Real >( "angle_tolerance", angle_tolerance_ ) );
	set_cacb_dihedral_tolerance( tag->getOption< core::Real >( "cacb_dihedral_tolerance", cacb_dihedral_tolerance_ ) );
	set_bb_dihedral_tolerance( tag->getOption< core::Real >( "bb_dihedral_tolerance", bb_dihedral_tolerance_ ) );
	set_weight( tag->getOption< core::Real >( "weight", weight_ ) );
	set_constrain_ca_ca_dist( tag->getOption< bool >( "constrain_ca_ca_dist", constrain_ca_ca_dist_ ) );
	set_constrain_bb_cacb_dihedral( tag->getOption< bool >( "constrain_cacb_dihedral", constrain_bb_cacb_dihedral_ ) );
	set_constrain_bb_dihedral( tag->getOption< bool >( "constrain_bb_dihedral", constrain_bb_dihedral_ ) );
	set_constrain_bb_angle( tag->getOption< bool >( "constrain_bb_angle", constrain_bb_angle_ ) );
}

/// @brief sets the blueprint file used for determining proper sheet pairing
/// This function will create the blueprint object for you from the file and use it
void
SheetConstraintGenerator::initialize_from_blueprint( std::string const & blueprint_file )
{
	if ( blueprint_file == "" ) {
		utility_exit_with_message( "SheetConstraintGenerator requires a blueprint file" );
	}
	BluePrintOP bp( new protocols::jd2::parser::BluePrint( blueprint_file ) );
	if ( ! bp ) {
		utility_exit_with_message( "SheetConstraintGenerator tried to read a blueprint file, but failed to create the proper object." );
	}
	initialize_from_blueprint( bp );
	TR << "Using blueprint " << blueprint_file << " for constraint generation." << std::endl;
}

void
SheetConstraintGenerator::initialize_from_blueprint( protocols::jd2::parser::BluePrintCOP bp )
{
	set_secstruct( bp->secstruct() );
	set_strand_pairs( bp->strand_pairings() );
}

void
SheetConstraintGenerator::set_flat_bottom_constraints( bool const flat_csts )
{
	flat_bottom_constraints_ = flat_csts;
}

void
SheetConstraintGenerator::set_secstruct( std::string const & ss )
{
	secstruct_ = ss;
}

void
SheetConstraintGenerator::set_strand_pairs( std::string const & spairs )
{
	spairs_ = spairs;
}

/// @brief set weight
void
SheetConstraintGenerator::set_weight( Real const coef )
{
	weight_ = coef;
}

/// @brief set distance of constraint
void
SheetConstraintGenerator::set_distance( Real const dist )
{
	dist_ = dist;
}

void
SheetConstraintGenerator::set_distance_tolerance( Real const dist_tol )
{
	dist_tolerance_ = dist_tol;
}

/// @brief set the flat-bottom tolerance for the backbone angle between strands for each pair
/// This is N1-C1-C2 and N2-C2-C1 for parallel sheets, and N1-C1-N2/N2-C2-N1 for antiparallel.
void
SheetConstraintGenerator::set_angle_tolerance( Real const angle_tolerance )
{
	angle_tolerance_ = angle_tolerance;
}

/// @brief set the flat-bottom tolerance for the Cb1-Ca1-Ca2-Cb2 dihedral angle (0 = optimal)
void
SheetConstraintGenerator::set_cacb_dihedral_tolerance( Real const dihedral_tolerance )
{
	cacb_dihedral_tolerance_ = dihedral_tolerance;
}

/// @brief set the flat-bottom tolerance for the backbone dihedrals (0=optimal)
/// Dihedral 1 = O1-N1-C1-C2, Dihedral 2 = O2-N2-C2-C1
void
SheetConstraintGenerator::set_bb_dihedral_tolerance( Real const dihedral_tolerance )
{
	bb_dihedral_tolerance_ = dihedral_tolerance;
}

/// @brief sets whether we should constrain distance only, and not generate dihedral and angle constraints
void
SheetConstraintGenerator::set_constrain_ca_ca_dist( bool const constrain_dist_only )
{
	constrain_ca_ca_dist_ = constrain_dist_only;
}

/// @brief sets whether we should constrain cb_ca_ca_cb dihedrals
void
SheetConstraintGenerator::set_constrain_bb_cacb_dihedral( bool const constrain_dihedral )
{
	constrain_bb_cacb_dihedral_ = constrain_dihedral;
}

/// @brief sets whether we should constrain bb dihedrals
void
SheetConstraintGenerator::set_constrain_bb_dihedral( bool const constrain_dihedral )
{
	constrain_bb_dihedral_ = constrain_dihedral;
}

/// @brief sets whether we should constrain bb angle
void
SheetConstraintGenerator::set_constrain_bb_angle( bool const constrain_angle )
{
	constrain_bb_angle_ = constrain_angle;
}

std::pair< std::string, std::string >
SheetConstraintGenerator::get_secstruct_and_strandpairings( core::pose::Pose const & pose ) const
{
	std::string secstruct = secstruct_;
	std::string spairs = spairs_;

	// fill in gaps with structuredata information, if necessary
	denovo_design::components::StructureDataOP sd =
		denovo_design::components::StructureDataFactory::get_instance()->create_from_pose( pose, id() );

	if ( secstruct.empty() || spairs.empty() ) {
		if ( secstruct.empty() ) {
			secstruct = sd->ss();
			TR << "StructureData information is used for determinining secondary structure: " << secstruct << std::endl;
		}
		if ( spairs.empty() ) {
			spairs = protocols::denovo_design::get_strandpairings( *sd, true );
			TR << "StructureData information is used for determinining constraint residue pairs: " << spairs << std::endl;
		}
	}

	// catch the case where ligands aren't present, e.g. when initialized from blueprint
	if ( secstruct.size() != pose.total_residue() ) {
		std::string ss_with_ligands = "";
		for ( core::Size res=1; res<=pose.total_residue(); ++res ) {
			if ( pose.residue( res ).is_protein() ) {
				ss_with_ligands += secstruct[ res - 1 ];
			} else {
				ss_with_ligands += 'L';
			}
		}
		if ( ss_with_ligands.size() == pose.total_residue() ) {
			secstruct = ss_with_ligands;
		} else {
			std::stringstream msg;
			msg << "SheetConstraintGenerator: given secondary structure " << secstruct
				<< "(length " << secstruct.size() << ") does not match pose size ("
				<< pose.total_residue() << ")" << std::endl;
			throw utility::excn::EXCN_BadInput( msg.str() );
		}
	}
	debug_assert( secstruct.size() == pose.total_residue() );
	if ( secstruct.empty() ) {
		throw utility::excn::EXCN_Msg_Exception( "SheetConstraintGenerator: secstruct is empty!" );
	}
	if ( spairs.empty() ) {
		throw utility::excn::EXCN_Msg_Exception( "SheetConstraintGenerator: strand pairs list is empty!" );
	}
	debug_assert( !secstruct.empty() );
	debug_assert( !spairs.empty() );
	return std::make_pair( secstruct, spairs );
}

ResiduePairs
SheetConstraintGenerator::compute_residue_pairs( topology::StrandPairings const & spairs ) const
{
	ResiduePairs res_pairs;
	for ( topology::StrandPairings::const_iterator sp=spairs.begin(); sp!=spairs.end(); ++sp ) {
		TR << "Pair " << **sp << " from " << (*sp)->begin1() << " to " << (*sp)->end1() << " has bulge? " << (*sp)->has_bulge() << std::endl;
		for ( core::Size resid=(*sp)->begin1(); resid<=(*sp)->end1(); ++resid ) {
			if ( (*sp)->is_bulge( resid ) ) {
				TR << "Skipping residue " << resid << " because it is a bulge" << std::endl;
				continue;
			}
			core::Size const paired_resid = (*sp)->residue_pair( resid );
			res_pairs.push_back( ResiduePair( resid, paired_resid, (*sp)->orient() ) );
			TR.Debug << "Added paired residues " << resid << " <--> " << paired_resid << std::endl;
		}
	}
	return res_pairs;
}

core::scoring::func::FuncOP
SheetConstraintGenerator::weighted_func( core::scoring::func::FuncOP func ) const
{
	return protocols::constraint_generator::scalar_weighted( func, weight_ );
}

core::scoring::func::FuncOP
SheetConstraintGenerator::create_caca_atom_pair_func( core::Real const ideal_dist ) const
{
	using namespace core::scoring::func;
	core::Real const lb = 0.0;
	std::string const tag( "constraints_in_beta_sheet" );
	if ( flat_bottom_constraints_ ) {
		return weighted_func( FuncOP( new core::scoring::constraints::BoundFunc( lb, ideal_dist, dist_tolerance_, tag ) ) );
	} else {
		return weighted_func( FuncOP( new HarmonicFunc( ideal_dist, dist_tolerance_ ) ) );
	}
}

core::scoring::func::FuncOP
SheetConstraintGenerator::create_bb_angle_func( core::Real const ideal_angle ) const
{
	using core::scoring::constraints::BoundFunc;
	using namespace core::scoring::func;
	if ( flat_bottom_constraints_ ) {
		return weighted_func( FuncOP( new BoundFunc(
			ideal_angle-angle_tolerance_,
			ideal_angle+angle_tolerance_, sqrt(1.0/42.0), "angle_bb") ) );
	} else {
		return weighted_func( FuncOP( new CircularHarmonicFunc( ideal_angle, angle_tolerance_ ) ) );
	}
}

core::scoring::func::FuncOP
SheetConstraintGenerator::create_bb_dihedral_func( core::Real const ideal_dihedral ) const
{
	using core::scoring::constraints::OffsetPeriodicBoundFunc;
	using namespace core::scoring::func;
	core::Real const periodicity = numeric::constants::f::pi;
	if ( flat_bottom_constraints_ ) {
		return weighted_func( FuncOP( new OffsetPeriodicBoundFunc(
			ideal_dihedral-bb_dihedral_tolerance_,
			ideal_dihedral+bb_dihedral_tolerance_,
			std::sqrt(1.0/42.0), "dihed_bb", periodicity, 0.0 ) ) );
	} else {
		return weighted_func( FuncOP( new CircularHarmonicFunc( ideal_dihedral, bb_dihedral_tolerance_ ) ) );
	}
}

core::scoring::func::FuncOP
SheetConstraintGenerator::create_cacb_dihedral_func( core::Real const ideal_dihedral ) const
{
	using core::scoring::constraints::OffsetPeriodicBoundFunc;
	using namespace core::scoring::func;
	core::Real const periodicity = 2 * numeric::constants::f::pi;
	if ( flat_bottom_constraints_ ) {
		return weighted_func( FuncOP( new OffsetPeriodicBoundFunc(
			ideal_dihedral-cacb_dihedral_tolerance_,
			ideal_dihedral+cacb_dihedral_tolerance_,
			std::sqrt(1.0/42.0), "dihed_cacb", periodicity, 0.0 ) ) );
	} else {
		return weighted_func( FuncOP( new CircularHarmonicFunc( ideal_dihedral, cacb_dihedral_tolerance_ ) ) );
	}
}

core::scoring::constraints::ConstraintOP
SheetConstraintGenerator::create_bb_dihedral_constraint(
	core::pose::Pose const & pose,
	core::Size const res1,
	core::Size const res2,
	core::scoring::func::FuncOP func1,
	core::scoring::func::FuncOP func2 ) const
{
	using namespace core::scoring::constraints;
	core::id::AtomID const resi_n( pose.residue_type( res1 ).atom_index( "N" ), res1 );
	core::id::AtomID const resi_c( pose.residue_type( res1 ).atom_index( "C" ), res1 );
	core::id::AtomID const resi_o( pose.residue_type( res1 ).atom_index( "O" ), res1 );
	core::id::AtomID const resj_c( pose.residue_type( res2 ).atom_index( "C" ), res2 );
	ConstraintCOPs const csts = boost::assign::list_of
		( ConstraintOP( new DihedralConstraint( resi_o, resi_n, resi_c, resj_c, func1 ) ) )
		( ConstraintOP( new DihedralConstraint( resi_o, resi_n, resi_c, resj_c, func2 ) ) );
	return ConstraintOP( new AmbiguousConstraint( csts ) );
}

core::scoring::constraints::ConstraintOP
SheetConstraintGenerator::create_ca_ca_atom_pair_constraint(
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	core::scoring::func::FuncOP func ) const
{
	using core::scoring::constraints::AtomPairConstraint;
	return  core::scoring::constraints::ConstraintOP( new AtomPairConstraint( atom1, atom2, func ) );
}

core::scoring::constraints::ConstraintOP
SheetConstraintGenerator::create_bb_angle_constraint(
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	core::id::AtomID const & atom3,
	core::scoring::func::FuncOP func ) const
{
	using core::scoring::constraints::AngleConstraint;
	return core::scoring::constraints::ConstraintOP( new AngleConstraint( atom1, atom2, atom3, func ) );
}

core::scoring::constraints::ConstraintOP
SheetConstraintGenerator::create_bb_cacb_dihedral_constraint(
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	core::id::AtomID const & atom3,
	core::id::AtomID const & atom4,
	core::scoring::func::FuncOP func ) const
{
	return core::scoring::constraints::ConstraintOP( new core::scoring::constraints::DihedralConstraint( atom1, atom2, atom3, atom4, func ) );
}

} //namespace fldsgn
} //namespace protocols
