// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief contains declarations for a type of constraint that holds a number of constraints
/// @brief where only the lowest one is evaluated
/// @author Florian Richter (floric@u.washington.edu, march 2008)


#include <core/scoring/constraints/AmbiguousNMRConstraint.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/AmbiguousNMRDistanceConstraint.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <utility/excn/Exceptions.hh>
#include <basic/prof.hh>

//Auto Headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace scoring {
namespace constraints {

////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousNMRConstraint::AmbiguousNMRConstraint( func::FuncOP func ):
	MultiConstraint( atom_pair_constraint ),
	func_( func )
{
	//	init_cst_score_types();
debug_assert ( member_constraints().size() == 0 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief Constructor
AmbiguousNMRConstraint::AmbiguousNMRConstraint( ConstraintCOPs const& cst_in, func::FuncOP func ):
	MultiConstraint( cst_in, atom_pair_constraint ),
	func_( func )
{
	//	init_cst_score_types();
debug_assert ( member_constraints().size() > 0 );
}
////////////////////////////////////////////////////////////////////////////////////////////////////
// void
// AmbiguousNMRConstraint::init_cst_score_types()
// {
// 	cst_score_types_.clear();
// 	cst_score_types_.push_back(atom_pair_constraint);
// }
////////////////////////////////////////////////////////////////////////////////////////////////////
/// @brief ScoreFunction, scores all member constraints but only reports the lowest one
void
AmbiguousNMRConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & /*weights*/, EnergyMap & emap ) const
{

	core::Real cum_invdist6 = 0;

	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
		AmbiguousNMRDistanceConstraintCOP cst_in_casted;
		cst_in_casted = utility::pointer::dynamic_pointer_cast< AmbiguousNMRDistanceConstraint const >( *member_it );
		if ( cst_in_casted ) cum_invdist6 += cst_in_casted->inv_dist6( xyz_func );
		if ( !cst_in_casted ) {
			AtomPairConstraintCOP cst_in_casted;
			cst_in_casted = utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( *member_it );
			if ( cst_in_casted ) {
				Real dist = cst_in_casted->dist( xyz_func );
				Real inv_dist = 1.0/dist;
				Real inv_dist2 = inv_dist*inv_dist;
				cum_invdist6 += inv_dist2*inv_dist2*inv_dist2;
			} else {
				runtime_assert( 0 == 1 );
			}
		}
	}
	//add lowest score to the actual emap
	Real eff_dist = pow( cum_invdist6, -1.0/6 );
	emap[ score_type() ] +=  get_func().func( eff_dist );
} //score


core::Real
AmbiguousNMRConstraint::dist( core::pose::Pose const& pose ) const
{
	return dist( func::ConformationXYZ( pose.conformation() ) );
}

core::Real
AmbiguousNMRConstraint::dist( func::XYZ_Func const & xyz ) const
{
	core::Real cum_invdist6 = 0;

	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
		AmbiguousNMRDistanceConstraintCOP cst_in_casted;
		PROF_START( basic::NOESY_ASSIGN_DIST_CST_CAST );
		cst_in_casted = utility::pointer::dynamic_pointer_cast< AmbiguousNMRDistanceConstraint const >( *member_it );
		PROF_STOP( basic::NOESY_ASSIGN_DIST_CST_CAST );
		if ( cst_in_casted ) cum_invdist6 += cst_in_casted->inv_dist6( xyz );
		if ( !cst_in_casted ) {
			AtomPairConstraintCOP cst_in_casted;
			PROF_START( basic::NOESY_ASSIGN_DIST_CST_CAST );
			cst_in_casted = utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( *member_it );
			PROF_STOP( basic::NOESY_ASSIGN_DIST_CST_CAST );
			if ( cst_in_casted ) {
				Real dist = cst_in_casted->dist( xyz );
				Real inv_dist = 1.0/dist;
				Real inv_dist2 = inv_dist*inv_dist;
				cum_invdist6 += inv_dist2*inv_dist2*inv_dist2;
			} else {
				runtime_assert( 0 == 1 );
			}
		}
	}
	//add lowest score to the actual emap
	Real eff_dist = pow( cum_invdist6, -1.0/6 );
	return eff_dist;
}

ConstraintOP
AmbiguousNMRConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	ConstraintCOPs new_csts;
	for( ConstraintCOPs::const_iterator cst_it = member_constraints_.begin(); cst_it != member_constraints_.end(); ++cst_it ){
		ConstraintOP new_cst = (*cst_it)->remap_resid( seqmap );
		if( new_cst ) new_csts.push_back( new_cst );
	}
	if( new_csts.size() > 0 ){
		return ConstraintOP( new AmbiguousNMRConstraint( new_csts, get_func().clone() ) );
	}
	else return NULL;
}


/// @brief function to minimize lowest scoring member constraint
void
AmbiguousNMRConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	//active_constraint_->fill_f1_f2(atom, conformation, F1, F2, weights);
	Real eff_dist = dist( xyz );
	Real out_wderiv( weights[ score_type() ] * get_func().dfunc( eff_dist ));
	Real in_deriv = -1.0/6.0 * pow( eff_dist, 7.0 );
	//	tr.Trace << "deriv for atom " << atom << eff_dist << " " << out_wderiv << " " << in_deriv << std::endl;

	//	tr.Trace << "the_other_atoms: " << the_other_atoms.size() << " " << the_other_atoms.front() << std::endl;
	//Vector f1(0.0), f2(0.0);
	for( ConstraintCOPs::const_iterator member_it = member_constraints().begin(); member_it != member_constraints().end(); member_it++){
              Vector f1(0.0), f2(0.0);

               //fpd hack to get at vector from atom i->j
                (*member_it)->fill_f1_f2( atom, xyz, f1, f2, weights );
    core::Real member_cst_is_scaled_by = weights[ (*member_it)->score_type() ] * (*member_it)->get_func().dfunc( eff_dist
);

               core::Real scale_i = -6.0*pow((*member_it)->dist( xyz ),-7.0);

               if (std::fabs(member_cst_is_scaled_by) > 1e-14) scale_i /= member_cst_is_scaled_by;

               F1 += 1.0 * scale_i * out_wderiv * in_deriv * f1;
               F2 += 1.0 * scale_i * out_wderiv * in_deriv * f2;
        }
        //              tr.Trace << "wderiv " << wderiv << std::endl;


	//		tr.Trace << "wderiv " << wderiv << std::endl;
//	F1 += out_wderiv * in_deriv * f1;
//	F2 += out_wderiv * in_deriv * f2;
}


// void
// AmbiguousNMRConstraint::show( std::ostream& out) const
// {
// 	out << "AmbiguousNMRConstraint Active constraint:" << std::endl;
// 	out << "AmbiguousNMRConstraint containing the following " << member_constraints().size() << " constraints: " << std::endl;
// 	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
// 		(*cst_it)->show(out);
// 	}

// 	out << " ...all member constraints of this AmbiguousNMRConstraint shown." << std::endl;
// }

////////////////////////////////////////////////////////////////////////////////////////////////////
///@details read definition of a multiconstraint. Since a MultiConstraint is essentially a vector of
void
AmbiguousNMRConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	func::FuncFactory const & func_factory
)
{
	std::string func_type;
	data >> func_type;

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
	//chu skip the rest of line since this is a single line defintion.
		while( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}
	MultiConstraint::read_def( data, pose, func_factory );
}

void AmbiguousNMRConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " ";
	if ( func_ ) func_->show_definition( out );
	else out << std::endl;
  for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
    (*cst_it)->show_def( out, pose );
		//		out<<std::endl;
  }
	out << "End_"<< type() << std::endl;
}

Size
AmbiguousNMRConstraint::show_violations( std::ostream& out, pose::Pose const& pose, Size verbose_level, Real threshold ) const
{
	Size total_viol = 0;
	bool passed( false );
	for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
		Size viol = (*cst_it)->show_violations( out, pose, verbose_level, threshold);
		if ( viol == 0 && verbose_level > 70 ) {
			utility::vector1< int > pos_list( (*cst_it)->residues() );
			passed = true;
			out << "\nResiduePairConstraints (" <<  pos_list[1] << ", " << pos_list[pos_list.size()] << " ) . of total: 1  0 violated" << std::endl;
		}
		total_viol += viol;
	}
	if ( !passed && verbose_level > 70 ) {
		for( ConstraintCOPs::const_iterator cst_it = member_constraints().begin(); cst_it != member_constraints().end(); cst_it++){
			utility::vector1< int > pos_list( (*cst_it)->residues() );
			out << "\nResiduePairConstraints (" <<  pos_list[1] << ", " << pos_list[pos_list.size()] << " ) + of total: 1  1 violated" << std::endl;
		}
	}
	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

bool cst_eq( Constraint const& cst1, Constraint const& cst2 ) {
	if ( cst1.natoms() != cst2.natoms() ) return false;
	for ( Size i=1; i<=cst1.natoms(); ++i ) {
		if ( cst1.atom( i ) != cst2.atom( i ) ) return false;
	}
	return true;
}

void
AmbiguousNMRConstraint::add_individual_constraint( ConstraintCOP cst_in )
{
	ConstraintCOP cst_in_casted;
	cst_in_casted = utility::pointer::dynamic_pointer_cast< AmbiguousNMRDistanceConstraint const >( cst_in );
	if ( !cst_in_casted ) {
		cst_in_casted = utility::pointer::dynamic_pointer_cast< AtomPairConstraint const >( cst_in );
	}
	if ( !cst_in_casted ) {
		throw utility::excn::EXCN_BadInput( "failed attempt to add " + cst_in->type() + " to AmbiguousNMRConstraint. Can only add AmbiguousNMRDistanceConstraint and AtomPairConstraint");
	}
	//is it unique ? -- otherwise reject constraint
	for ( ConstraintCOPs::const_iterator it=member_constraints().begin(); it != member_constraints().end(); ++it ) {
		if ( cst_eq( **it, *cst_in ) ) return; //we got it already...
	}

	MultiConstraint::add_individual_constraint( cst_in );
}


} //constraints
} //scoring
} //core
