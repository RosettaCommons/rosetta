// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/BackboneStubLinearConstraint.cc
///
/// @brief
/// @author Lei Shi


#include <core/scoring/constraints/BackboneStubLinearConstraint.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
// AUTO-REMOVED #include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/func/PeriodicFunc.hh>
#include <basic/Tracer.hh>

// used to make temporary alanines for gly cst's
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

//Auto Headers
#include <core/id/SequenceMapping.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/Jump.hh>

namespace core {
namespace scoring {
namespace constraints {


static thread_local basic::Tracer TR( "core.scoring.constraints.BackboneStubLinearConstraint" );

utility::pointer::owning_ptr< AngleConstraint > BackboneStubLinearConstraint::ang_cst_(0);

BackboneStubLinearConstraint::BackboneStubLinearConstraint(
	pose::Pose const & pose,
	Size const seqpos,
	AtomID const & fixed_atom_id,
	conformation::Residue const & target_rsd,
	core::Real const & superposition_bonus,
	core::Real const & CB_force_constant
):
	Constraint( core::scoring::backbone_stub_linear_constraint ),
	superposition_bonus_( superposition_bonus ),
	CB_force_constant_( CB_force_constant ),
	seqpos_( seqpos ),
	fixed_atom_id_( fixed_atom_id )
{
	// store info about the target residue
	assert( target_rsd.is_protein() );
	if ( (target_rsd.aa() == chemical::aa_gly ) ) {
		TR << "ERROR - Gly residues cannot be used in BackboneStubLinearConstraints." << std::endl;
		return;
	}
	CB_target_ = target_rsd.xyz("CB");
	CA_target_ = target_rsd.xyz("CA");
	C_target_ = target_rsd.xyz("C");
	N_target_ = target_rsd.xyz("N");
	CB_CA_target_ = CB_target_ - CA_target_;
	C_N_target_ = C_target_ - N_target_;

	// constraint depends on CB, CA, C coordinates
	conformation::Residue const & rsd( pose.residue(seqpos_) );
	assert( rsd.is_protein() );
	if ( (rsd.aa() == chemical::aa_gly) ) {
		TR << "ERROR - Gly residues cannot be used in BackboneStubLinearConstraints." << std::endl;
		return;
	}
	CB_atom_id_ = AtomID( rsd.atom_index("CB"), seqpos_ );
	atom_ids_.push_back(CB_atom_id_);
	CA_atom_id_ = AtomID( rsd.atom_index("CA"), seqpos_ );
	atom_ids_.push_back(CA_atom_id_);
	C_atom_id_ = AtomID( rsd.atom_index("C"), seqpos_ );
	atom_ids_.push_back(C_atom_id_);
	N_atom_id_ = AtomID( rsd.atom_index("N"), seqpos_ );
	atom_ids_.push_back( N_atom_id_ );

	// need a fixed reference point on the pose
	atom_ids_.push_back(fixed_atom_id_);
	// we don't actually need to save the coors of the reference point,
	// but will do so to ensure that it doesn't move
	fixed_reference_point_ = pose.xyz(fixed_atom_id_);

	// to get access to AngleConstraint derivatives
	if ( ang_cst_ == 0 ) {
		// note: PeriodicFunc has functional form y = ( k * cos(n * (x - x0) ) ) + C
		func::FuncOP cos_func = new func::PeriodicFunc(0., 1., 1., 0.);
		ang_cst_ = AngleConstraintOP( new AngleConstraint( cos_func ) );
	}
}

core::Size BackboneStubLinearConstraint::seqpos() const
{ return seqpos_; }

void BackboneStubLinearConstraint::show( std::ostream& out ) const
{
	out << "BackboneStubCst Seqpos: " << seqpos_ << "    bonus: " << superposition_bonus_ << "    CB force constant: " << CB_force_constant_ << std::endl;
}

bool BackboneStubLinearConstraint::operator == ( Constraint const & other_cst ) const
{

	if( !dynamic_cast< BackboneStubLinearConstraint const * > ( &other_cst ) ) return false;

	BackboneStubLinearConstraint const & other( static_cast< BackboneStubLinearConstraint const & > (other_cst) );

	if( superposition_bonus_ != other.superposition_bonus_ ) return false;
	if( CB_force_constant_ != other.CB_force_constant_ ) return false;
	if( seqpos_ != other.seqpos_ ) return false;

	if( CB_atom_id_ != other.CB_atom_id_ ) return false;
	if( CA_atom_id_ != other.CA_atom_id_ ) return false;
	if( C_atom_id_ != other.C_atom_id_ ) return false;
	if( N_atom_id_ != other.N_atom_id_ ) return false;

	if( CB_target_ != other.CB_target_ ) return false;
	if( CA_target_ != other.CA_target_ ) return false;
	if( C_target_ != other.C_target_ ) return false;
	if( N_target_ != other.N_target_ ) return false;
	if( CB_CA_target_ != other.CB_CA_target_ ) return false;
	if( C_N_target_ != other.C_N_target_ ) return false;

	if( fixed_atom_id_ != other.fixed_atom_id_ ) return false;
	if( fixed_reference_point_ != other.fixed_reference_point_ ) return false;

	return true;
}

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
// the constraint score, C, is given by C = (b + kx^2) * z, where b is the score bonus for the stub, k is the user-defined CB force constant
// z is the dot product between the Ca-Cb vectors of the stub and the current
// position on the scaffold times the C-N vectors of the stub and scaffold position.
// The intuition in this is that a stub that lands
// perfectly on the scaffold's bacbkone position AND for which the two vectors
// match the scaffold's perfectly would give maximal bonus. This choice of vectors is
// orthogonal and so provides unbiased resolution for conformation space.
// Note that we don't know in advance which phi/psi angles a scaffold position would have, and
// that may change the C_N vector by up to 25o from the stub's. Since the cosine of small angular
// differences is close to 1, this uncertainty would not cause large problems.
// Since the constraint does not produce positive values, the effect of a stub
// on the energy is only if the distance between the Cb's of stub and scaffold is
// sqrt( -b / k ). As a rule of thumb, b=-4, so k can be decided in a way that will
// determine the radius of the effect of a stub on pulling a scaffold towards it: k = 4/dist^2, where dist is the preferred radius

//Lei Shi
// constraint score, C, is given by C = b + k*(1-exp(-x^2/2)+z_Cab^2+z_Ccn^2
// this should help better the ranking of patchdock solutions and constrainted guided docking
void BackboneStubLinearConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const {

	if ( weights[ this->score_type() ] == 0 ) return;

	core::conformation::Residue const & curr_rsd = xyz_func.residue(seqpos_);
	assert( curr_rsd.is_protein() );

	// verify that the fixed reference point is in the same place
	core::Vector curr_ref_location = xyz_func(fixed_atom_id_);
	core::Real ref_dist = curr_ref_location.distance_squared( fixed_reference_point_ );
	if ( ref_dist > 1E-8 ) {
		TR << "ERROR - BackboneStubLinearConstraint requires a fixed reference atom, but this atom has moved!!" << std::endl;
		TR << "Reference location was " << fixed_reference_point_.x() << ", " << fixed_reference_point_.y() << ", " << fixed_reference_point_.z() << ", now it's " << curr_ref_location.x() << ", " << curr_ref_location.y() << ", " << curr_ref_location.z() << std::endl;
		std::exit(1);
	}

	// return a value between superposition_bonus_ (-ve) and zero
  // This bonus really does nothing. Just used to keep the same interface
	core::Real cst_val(0.);
	//core::Real cst_val(superposition_bonus_);
	//assert( cst_val < 0. );

	// apply a harmonic constraint on the CB's
	core::Vector CB_curr,CA_curr,C_curr,N_curr;
	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
			TR << "ERROR - Gly residues cannot be used in BackboneStubLinearConstraints." << std::endl;
			return;
	}

	CB_curr = curr_rsd.xyz("CB");
	core::Real CB_d2 = CB_curr.distance_squared( CB_target_ );
  cst_val += CB_force_constant_ * CB_d2;
  //TR << "CB_CB_dist curr: "<< CB_d2 << "," << CB_curr[0] << "," << CB_curr[1] << ","<< CB_curr[2] << std::endl;
  //TR << "CB_CB_dist target: "<< CB_d2 << "," << CB_target_[0] << "," << CB_target_[1] << ","<< CB_target_[2] << std::endl;
	//cst_val += CB_force_constant_*(1-std::exp(-CB_d2/2));

	CA_curr = curr_rsd.xyz("CA");
	core::Real CA_d2 = CA_curr.distance_squared( CA_target_ );
  cst_val += CA_d2;

	C_curr = curr_rsd.xyz("C");
	core::Real C_d2 = C_curr.distance_squared( C_target_ );
  cst_val += CB_force_constant_ * C_d2;

	N_curr = curr_rsd.xyz("N");
	core::Real N_d2 = N_curr.distance_squared( N_target_ );
  cst_val += N_d2;

/*
	// plus by the (cos-1)^2 of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
	//cst_val *= (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)*(ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)/2;
	cst_val += 1.0-ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_);
	//cst_val += (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)*(ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)/2;
*/

/*
	// plus by the (cos-1)^2 of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
	//cst_val *= (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)*(ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)/2;
	cst_val += 1.0-ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ );
	//cst_val += (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)*(ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)/2;
*/

  //TR << "cst_val dist: " << CB_force_constant_*(1-std::exp(-CB_d2/2)) << std::endl;
  //TR << "cst_val angle1: " << (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)*(ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_)-1)/2 << std::endl;
  //TR << "cst_val angle2: " << (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)*(ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)/2 << std::endl;
  //TR << "cst_val bonus: " << superposition_bonus_ << std::endl;
  //TR << "cst_val geometry: " << cst_val-superposition_bonus_ << std::endl;
  //TR << "cst_val total: " << cst_val  << std::endl;

	emap[ this->score_type() ] += cst_val;
}

void
BackboneStubLinearConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	if ( weights[ this->score_type() ] == 0 ) return;

	core::conformation::Residue const & curr_rsd = xyz.residue( seqpos_ );

	assert( curr_rsd.is_protein() );

	if ( ( atom != CB_atom_id_ ) && ( atom != CA_atom_id_ ) && ( atom != C_atom_id_ ) && ( atom != N_atom_id_ ) ) {
		return;
	}

	// return a value between superposition_bonus_ (-ve) and zero
	core::Real cst_val(0.);
	//core::Real cst_val(superposition_bonus_);
	//assert( cst_val < 0. );
  //TR << " cst_val1: "<< cst_val << std::endl;

	// start by computing the cst value normally, collecting score components along the way. if cst_val is non-negative no derivative is computed

	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
			TR << "ERROR - Gly residues cannot be used in BackboneStubLinearConstraints." << std::endl;
			return;
	}

	// apply a harmonic constraint on the CB's
	core::Vector const CB_curr( curr_rsd.xyz("CB") );
	core::Real const CB_d2 = CB_curr.distance_squared( CB_target_ );
	core::Real const CB_pos_term = CB_force_constant_ * CB_d2;
	cst_val += CB_pos_term;

	core::Vector const CA_curr( curr_rsd.xyz("CA") );
	core::Real const CA_d2 = CA_curr.distance_squared( CA_target_ );
	core::Real const CA_pos_term = CA_d2;
	cst_val += CA_pos_term;

	core::Vector const C_curr( curr_rsd.xyz("C") );
	core::Real const C_d2 = C_curr.distance_squared( C_target_ );
	core::Real const C_pos_term = CB_force_constant_ * C_d2;
	cst_val += C_pos_term;

	core::Vector const N_curr( curr_rsd.xyz("N") );
	core::Real const N_d2 = N_curr.distance_squared( N_target_ );
	core::Real const N_pos_term = N_d2;
	cst_val += N_pos_term;


  //derivatives
	core::Real const constant_dist_term( weights[ this->score_type() ] );

  if ( atom == CB_atom_id_ ) {
    Vector const CB_f2( CB_curr - CB_target_ );
    core::Real const CB_dist( CB_f2.length() );
    core::Real const CB_deriv = 2 * CB_force_constant_ * CB_dist;
    if ( CB_deriv != 0.0 && CB_dist != 0.0 ) {
      Vector const CB_f1( CB_curr.cross( CB_target_ ) );
      F1 += ( ( CB_deriv / CB_dist ) * CB_f1 ) * constant_dist_term;
      F2 += ( ( CB_deriv / CB_dist ) * CB_f2 ) * constant_dist_term;
    }
    return;
  }

  if ( atom == CA_atom_id_ ) {
    Vector const CA_f2( CA_curr - CA_target_ );
    core::Real const CA_dist( CA_f2.length() );
    core::Real const CA_deriv = 2 * CA_dist;
    if ( CA_deriv != 0.0 && CA_dist != 0.0 ) {
      Vector const CA_f1( CA_curr.cross( CA_target_ ) );
      F1 += ( ( CA_deriv / CA_dist ) * CA_f1 ) * constant_dist_term;
      F2 += ( ( CA_deriv / CA_dist ) * CA_f2 ) * constant_dist_term;
    }
    return;
  }

  if ( atom == C_atom_id_ ) {
    Vector const C_f2( C_curr - C_target_ );
    core::Real const C_dist( C_f2.length() );
    core::Real const C_deriv = 2 * CB_force_constant_ * C_dist;
    if ( C_deriv != 0.0 && C_dist != 0.0 ) {
      Vector const C_f1( C_curr.cross( C_target_ ) );
      F1 += ( ( C_deriv / C_dist ) * C_f1 ) * constant_dist_term;
      F2 += ( ( C_deriv / C_dist ) * C_f2 ) * constant_dist_term;
    }
    return;
  }

  if ( atom == N_atom_id_ ) {
    Vector const N_f2( N_curr - N_target_ );
    core::Real const N_dist( N_f2.length() );
    core::Real const N_deriv = 2 * N_dist;
    if ( N_deriv != 0.0 && N_dist != 0.0 ) {
      Vector const N_f1( N_curr.cross( N_target_ ) );
      F1 += ( ( N_deriv / N_dist ) * N_f1 ) * constant_dist_term;
      F2 += ( ( N_deriv / N_dist ) * N_f2 ) * constant_dist_term;
    }
    return;
  }

/*
	// multiply by the cos of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
	core::Vector const CB_CA_curr = CB_curr - CA_curr;
	core::Real const CB_CA_angle_term( 1.0-ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ ));
	//core::Real const CB_CA_angle_term( (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1)*(ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1)/2 );
  cst_val += CB_CA_angle_term;
  //cst_val *= CB_CA_angle_term;
  //TR << " CB_CA_angle: "<< ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ ) << std::endl;
  //TR << " CB_CA_angle_term: "<< CB_CA_angle_term << std::endl;

	// multiply by the cos of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
	core::Vector const C_N_curr = C_curr - N_curr;
	core::Real const C_N_angle_term( 1.0-ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ ));
	//core::Real const C_N_angle_term( (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)*(ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1)/2 );
  cst_val += C_N_angle_term;
  //cst_val *= C_N_angle_term;
  //TR << " C_N_angle: "<< ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ ) << std::endl;
  //TR << " C_N_angle_term: "<< C_N_angle_term << std::endl;

	core::Real const constant_dist_term( weights[ this->score_type() ] );
	core::Real const constant_ang_term( constant_dist_term  * ( 0. + CB_pos_term ) );
	//core::Real const constant_ang_term( constant_dist_term  * ( superposition_bonus_ + CB_pos_term ) );

	// contribution from differentiating CB_dist ( * CA_angle_term * C_angle_term from product rule)
	// and then adding the effect of the angular constraint on cb using the chain rule
	if ( atom == CB_atom_id_ ) {
// the effects of the coordinate constraint
		Vector const CB_f2( CB_curr - CB_target_ );
		core::Real const CB_dist( CB_f2.length() );
		//core::Real const CB_deriv = CB_force_constant_ ;
		core::Real const CB_deriv = 2 * CB_force_constant_ * CB_dist;
		//core::Real const CB_deriv = CB_force_constant_ * (std::exp(-CB_d2/2)) * CB_dist;
		if ( CB_deriv != 0.0 && CB_dist != 0.0 ) {
			Vector const CB_f1( CB_curr.cross( CB_target_ ) );
			//F1 += ( ( CB_deriv / CB_dist ) * CB_f1 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
			//F2 += ( ( CB_deriv / CB_dist ) * CB_f2 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
			F1 += ( ( CB_deriv / CB_dist ) * CB_f1 ) * constant_dist_term;
			F2 += ( ( CB_deriv / CB_dist ) * CB_f2 ) * constant_dist_term;
		}

// the angular constraint on cb
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CB_curr, CA_curr, CA_curr + CB_CA_target_, partial_F1, partial_F2 );
		//F1 += partial_F1 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * C_N_angle_term * constant_ang_term;
		//F2 += partial_F2 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * C_N_angle_term * constant_ang_term;
		F1 += partial_F1 * (-1) * constant_dist_term;
		F2 += partial_F2 * (-1) * constant_dist_term;
		//F1 += partial_F1 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * constant_dist_term;
		//F2 += partial_F2 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * constant_dist_term;
		return;
	}

// the angular constraint on ca
// See N atom for explanation on the derivs
	if ( atom == CA_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CA_curr, CA_curr - CB_CA_curr, CA_curr + CB_CA_target_ - CB_CA_curr, partial_F1, partial_F2 );
		//F1 += -partial_F1 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * C_N_angle_term * constant_ang_term;
		//F2 += -partial_F2 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * C_N_angle_term * constant_ang_term;
		F1 += -partial_F1 * (-1) * constant_dist_term;
		F2 += -partial_F2 * (-1) * constant_dist_term;
		//F1 += -partial_F1 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * constant_dist_term;
		//F2 += -partial_F2 * (ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ )-1) * constant_dist_term;
		return;
	}

// the angular constraint on c
	if ( atom == C_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( C_curr, N_curr, N_curr + C_N_target_, partial_F1, partial_F2 );
		//F1 += partial_F1 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * CB_CA_angle_term * constant_ang_term;
		//F2 += partial_F2 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * CB_CA_angle_term * constant_ang_term;
		F1 += partial_F1 * (-1) * constant_dist_term;
		F2 += partial_F2 * (-1) * constant_dist_term;
		//F1 += partial_F1 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * constant_dist_term;
		//F2 += partial_F2 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * constant_dist_term;
		return;
	}

// the angular constraint on n
// This is tricky, and thanks to Frank for coming up with this!
// Ang_cst uses the actual coordinates of p1 in computing the derivatives, requiring it to be constant.
// But in our case, p1 changes with respect to p2->p3. Solution: translate the vectors by -C_N_curr.
// Now, p1 is at the site of the nonchanging vector, and so the coordinates are safe for ang_cst.
// However, the derivative for the angular constraint is reversed, so we multiply by -1
	if( atom == N_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( N_curr, N_curr - C_N_curr, N_curr + C_N_target_ - C_N_curr, partial_F1, partial_F2 );
		//F1 += -partial_F1 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * CB_CA_angle_term * constant_ang_term;
		//F2 += -partial_F2 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * CB_CA_angle_term * constant_ang_term;
		F1 += -partial_F1 * (-1) * constant_dist_term;
		F2 += -partial_F2 * (-1) * constant_dist_term;
		//F1 += -partial_F1 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * constant_dist_term;
		//F2 += -partial_F2 * (ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ )-1) * constant_dist_term;
	}
*/

	return;
}

ConstraintOP BackboneStubLinearConstraint::clone() const
{
	return ConstraintOP( new BackboneStubLinearConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP BackboneStubLinearConstraint::remapped_clone( pose::Pose const& /*src*/, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {

	core::Size new_seqpos = seqpos_;
	AtomID new_fixed_atom_id = fixed_atom_id_;

	if ( smap ) {
		new_seqpos = (*smap)[ seqpos_ ];
		new_fixed_atom_id.rsd() = (*smap)[ fixed_atom_id_.rsd() ];
		if( new_seqpos == 0 ) return NULL;
	}

	// make an alanine
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( alatype );
	ala->set_xyz("CB",CB_target_);
	ala->set_xyz("CA",CA_target_);
	ala->set_xyz("C",C_target_);
	ala->set_xyz("N",N_target_);

	return new BackboneStubLinearConstraint(dest, new_seqpos, new_fixed_atom_id, *ala, superposition_bonus_, CB_force_constant_ );
}

} // namespace constraints
} // namespace scoring
} // namespace core
