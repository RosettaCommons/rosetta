// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/ncbb/SecStructMinimizeMover.cc
/// @brief
/// @detailed
/// @author Andy Watkins, andy.watkins2@gmail.com


#include <protocols/ncbb/SecStructMinimizeMover.hh>
#include <protocols/ncbb/SecStructMinimizeMoverCreator.hh>

#include <protocols/ncbb/SecStructMinimizeMultiFunc.hh>
#include <protocols/ncbb/util.hh>

#include <protocols/rosetta_scripts/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>

#include <protocols/simple_moves/MinMover.hh>

#include <numeric/conversions.hh>
#include <core/scoring/constraints/DihedralConstraint.hh>
#include <core/scoring/func/CircularHarmonicFunc.hh>

#include <basic/datacache/DataMap.hh>
#include <utility/tag/Tag.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

static basic::Tracer TR( "protocols.ncbb.SecStructMinimizeMover" );

using namespace core;

namespace protocols {
namespace ncbb {

//Destructor
SecStructMinimizeMover::~SecStructMinimizeMover()
{}

void
SecStructMinimizeMover::apply( core::pose::Pose & pose )
{
	using namespace core::optimization;

	Size number_dihedral_sets;
	utility::vector1< char > uniqs;
	count_uniq_char( dihedral_pattern_, number_dihedral_sets, uniqs );
	Size number_dihedrals = get_number_dihedrals( uniqs, dihedral_pattern_, alpha_beta_pattern_ );
	dihedrals_.resize( number_dihedrals );

	// Initialize dihedrals from pose.
	for ( Size i = 1, d = 1, m = 1; i <= pose.total_residue() && d <= number_dihedrals; ++d, ++m ) {
		TR.Trace << "Taking res " << i << " mc torsion " << m << " for dihedral " << d << "/" << number_dihedrals << std::endl;
		dihedrals_[ d ] = pose.residue( i ).mainchain_torsions()[ m ];

		// Assume no omega manipulation, ever.
		if ( pose.residue_type( i ).is_beta_aa() ) {
			if ( m == 3 ) {
				++i;
				m = 0;
			}
		} else if ( pose.residue_type( i ).is_alpha_aa() || pose.residue_type( i ).is_peptoid() ) {
			if ( m == 2 ) {
				++i;
				m = 0;
			}
		} else {
			utility_exit_with_message( "Non-alpha, non-beta, non-peptoid residues are not yet supported." );
		}
	}

	if ( constrain_ ) {
		add_dihedral_constraints_to_pose( pose, number_dihedral_sets, uniqs );
	}

	TR.Trace << "Dihedrals values: " << std::endl;
	for ( Size ii = 1; ii <= number_dihedrals; ++ii ) {
		TR.Trace << dihedrals_[ ii ] << "\t";
	}
	TR.Trace << std::endl;

	kinematics::MoveMap min_mm;
	min_mm.set_bb( true );
	min_mm.set_chi( true );

	MinimizerMap min_map;
	min_map.setup( pose, min_mm );

	( *score_fxn_ ) ( pose );
	SecStructMinimizeMultiFunc ssmmf( pose, *score_fxn_, min_map, alpha_beta_pattern_, dihedral_pattern_ );

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	MinimizerOptions minoptions( option[ run::min_type ], 0.0001, true, false, false ); // investigate the final bools?
	pose.energies().set_use_nblist( pose, min_map.domain_map(), false );
	Minimizer minimizer( ssmmf, minoptions );

	utility::vector1< Real > dihedrals_for_minimization;
	for ( Size index = 1; index <= dihedrals_.size() - 1; ++index ) {
		dihedrals_for_minimization.push_back( dihedrals_[ index ] );
	}
	minimizer.run( dihedrals_ );

	pose.energies().reset_nblist();
}

void
SecStructMinimizeMover::add_dihedral_constraints_to_pose(
	Pose & pose,
	Size number_dihedral_sets,
	utility::vector1< char > uniqs
) {
	using namespace core::scoring::constraints;
	using namespace core::scoring::func;
	using namespace numeric::conversions;

	for ( Size resi = 1; resi <= pose.n_residue(); ++resi ) {

		Size vec_index = 1;
		for ( Size i = 2; i <= number_dihedral_sets; ++i ) {
			if ( dihedral_pattern_[ resi-1 ] == uniqs[ i ] ) vec_index = i;
		}

		Size index = give_dihedral_index( vec_index, uniqs, dihedral_pattern_, alpha_beta_pattern_ );

		if ( pose.residue( resi ).type().is_beta_aa() ) {

			CircularHarmonicFuncOP dih_func_phi( new CircularHarmonicFunc( radians( dihedrals_[ index ] ), radians( 10 ) ) );
			CircularHarmonicFuncOP dih_func_tht( new CircularHarmonicFunc( radians( dihedrals_[ index+1 ] ), radians( 10 ) ) );
			CircularHarmonicFuncOP dih_func_psi( new CircularHarmonicFunc( radians( dihedrals_[ index+2 ] ), radians( 10 ) ) );

			AtomID aidC1( ( resi == 1 ) ? pose.residue( resi ).atom_index( "CO" )
				: pose.residue(resi-1).atom_index( "C" ),
				( resi == 1 ) ? resi : resi - 1 );


			AtomID aidN1( pose.residue( resi ).atom_index( "N" ), resi );
			AtomID aidCA( pose.residue( resi ).atom_index( "CA" ), resi );
			AtomID aidCM( pose.residue( resi ).atom_index( "CM" ), resi );
			AtomID aidC2( pose.residue( resi ).atom_index( "C" ), resi );
			AtomID aidN2( ( resi == pose.n_residue() ) ? pose.residue( resi ).atom_index( "NM" )
				: pose.residue(resi+1).atom_index( "N" ),
				( resi == pose.n_residue() ) ? resi : resi + 1 );


			ConstraintCOP phiconstraint( new DihedralConstraint( aidC1, aidN1, aidCA, aidCM, dih_func_phi ) );
			ConstraintCOP thtconstraint( new DihedralConstraint( aidN1, aidCA, aidCM, aidC2, dih_func_tht ) );
			ConstraintCOP psiconstraint( new DihedralConstraint( aidCA, aidCM, aidC2, aidN2, dih_func_psi ) );

			pose.add_constraint( phiconstraint );
			pose.add_constraint( thtconstraint );
			pose.add_constraint( psiconstraint );

		} else if ( pose.residue( resi ).type().is_alpha_aa() ) {

			CircularHarmonicFuncOP dih_func_phi( new CircularHarmonicFunc( radians( dihedrals_[ index ] ), radians( 10 ) ) );
			CircularHarmonicFuncOP dih_func_psi( new CircularHarmonicFunc( radians( dihedrals_[ index+1 ] ), radians( 10 ) ) );

			AtomID aidC1( ( resi == 1 ) ? pose.residue( resi ).atom_index( "CO" )
				: pose.residue(resi-1).atom_index( "C" ),
				( resi == 1 ) ? resi : resi - 1 );


			AtomID aidN1( pose.residue( resi ).atom_index( "N" ), resi );
			AtomID aidCA( pose.residue( resi ).atom_index( "CA" ), resi );
			AtomID aidC2( pose.residue( resi ).atom_index( "C" ), resi );
			AtomID aidN2( ( resi == pose.n_residue() ) ? pose.residue( resi ).atom_index( "NM" )
				: pose.residue(resi+1).atom_index( "N" ),
				( resi == pose.n_residue() ) ? resi : resi + 1 );

			ConstraintCOP phiconstraint( new DihedralConstraint( aidC1, aidN1, aidCA, aidC2, dih_func_phi ) );
			ConstraintCOP psiconstraint( new DihedralConstraint( aidN1, aidCA, aidC2, aidN2, dih_func_psi ) );

			pose.add_constraint( phiconstraint );
			pose.add_constraint( psiconstraint );
		}
	}
}

void
SecStructMinimizeMover::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	protocols::filters::Filters_map const &,
	protocols::moves::Movers_map const &,
	core::pose::Pose const &
) {
	alpha_beta_pattern_ = tag->getOption<std::string>("alpha_beta_pattern", "A" );
	constrain_ = tag->getOption<bool>("constrain", false );
	dihedral_pattern_ = tag->getOption<std::string>("dihedral_pattern", "A" );

	try {
		score_fxn_ = protocols::rosetta_scripts::parse_score_function( tag, "score_fxn", data )->clone();
	} catch ( utility::excn::EXCN_RosettaScriptsOption const & e ) {
		TR << "kill me " << e.msg() << std::endl;
	}
}

// MoverCreator
moves::MoverOP SecStructMinimizeMoverCreator::create_mover() const {
	return moves::MoverOP( new SecStructMinimizeMover() );
}

std::string SecStructMinimizeMoverCreator::keyname() const {
	return SecStructMinimizeMoverCreator::mover_name();
}

std::string SecStructMinimizeMoverCreator::mover_name() {
	return "SecStructMinimizeMover";
}

} //ncbb
} //protocols
