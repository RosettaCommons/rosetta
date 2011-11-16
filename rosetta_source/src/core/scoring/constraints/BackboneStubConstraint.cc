// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/BackboneStubConstraint.cc
///
/// @brief
/// @author John Karanicolas, Sarel Fleishman


#include <core/scoring/constraints/BackboneStubConstraint.hh>

// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
#include <core/conformation/Conformation.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/PeriodicFunc.hh>
#include <basic/Tracer.hh>

// used to make temporary alanines for gly cst's
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>

//Auto Headers
#include <platform/types.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/Adduct.fwd.hh>
#include <core/chemical/Adduct.hh>
#include <core/chemical/AtomICoor.fwd.hh>
#include <core/chemical/AtomICoor.hh>
#include <core/chemical/AtomType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.fwd.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/chemical/ElementSet.fwd.hh>
#include <core/chemical/MMAtomType.fwd.hh>
#include <core/chemical/MMAtomTypeSet.fwd.hh>
#include <core/chemical/ResConnID.fwd.hh>
#include <core/chemical/ResConnID.hh>
#include <core/chemical/ResidueConnection.fwd.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/chemical/ResidueSelector.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/chemical/VariantType.fwd.hh>
// AUTO-REMOVED #include <core/chemical/types.hh>
#include <core/chemical/orbitals/ICoorOrbitalData.hh>
#include <core/chemical/orbitals/OrbitalType.fwd.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.fwd.hh>
#include <core/chemical/sdf/MolData.fwd.hh>
#include <core/chemical/sdf/MolData.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Atom.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/PseudoBond.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/conformation/orbitals/OrbitalXYZCoords.hh>
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.fwd.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/DOF_ID_Map.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>
#include <core/id/NamedAtomID.fwd.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
#include <core/id/NamedStubID.fwd.hh>
#include <core/id/SequenceMapping.fwd.hh>
#include <core/id/SequenceMapping.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/types.hh>
#include <core/kinematics/AtomPointer.fwd.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/AtomWithDOFChange.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/Edge.fwd.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MinimizerMapBase.fwd.hh>
#include <core/kinematics/RT.fwd.hh>
#include <core/kinematics/RT.hh>
#include <core/kinematics/ResidueCoordinateChangeList.fwd.hh>
#include <core/kinematics/ShortestPathInFoldTree.fwd.hh>
#include <core/kinematics/Stub.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/tree/Atom.fwd.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/pose/PDBInfo.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/ObserverCache.fwd.hh>
#include <core/pose/metrics/PoseMetricContainer.fwd.hh>
#include <core/pose/signals/ConformationEvent.fwd.hh>
#include <core/pose/signals/DestructionEvent.fwd.hh>
#include <core/pose/signals/EnergyEvent.fwd.hh>
#include <core/pose/signals/GeneralEvent.fwd.hh>
#include <core/scoring/Energies.fwd.hh>
#include <core/scoring/EnergyMap.fwd.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/AngleConstraint.hh>
#include <core/scoring/constraints/BackboneStubConstraint.fwd.hh>
#include <core/scoring/constraints/Constraint.fwd.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/Func.hh>
#include <core/scoring/constraints/FuncFactory.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.fwd.hh>
#include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/constraints/XYZ_Func.fwd.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/down_cast.hh>
#include <utility/exit.hh>
#include <utility/vector0_bool.hh>
#include <utility/vector1.fwd.hh>
#include <utility/vector1.hh>
#include <utility/vector1_bool.hh>
#include <utility/vectorL.fwd.hh>
#include <utility/vectorL.hh>
#include <utility/vectorL_Selector.hh>
#include <utility/vectorL_bool.hh>
#include <utility/keys/Key2Tuple.fwd.hh>
#include <utility/keys/Key2Tuple.hh>
#include <utility/keys/Key3Tuple.fwd.hh>
#include <utility/keys/Key3Tuple.hh>
#include <utility/keys/Key4Tuple.fwd.hh>
#include <utility/keys/Key4Tuple.hh>
#include <utility/pointer/ReferenceCount.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/access_ptr.fwd.hh>
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.functions.hh>
#include <utility/pointer/owning_ptr.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/signals/BufferedSignalHub.fwd.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.fwd.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/LinkUnit.fwd.hh>
#include <utility/signals/LinkUnit.hh>
#include <utility/signals/PausableSignalHub.fwd.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/signals/SignalHub.fwd.hh>
#include <utility/signals/SignalHub.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/sphericalVector.fwd.hh>
#include <numeric/trig.functions.hh>
#include <numeric/xyz.functions.fwd.hh>
#include <numeric/xyzMatrix.fwd.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.fwd.hh>
#include <numeric/xyzVector.hh>
#include <numeric/internal/ColPointers.hh>
#include <numeric/internal/ColVectors.hh>
#include <numeric/internal/ColsPointer.hh>
#include <numeric/internal/RowPointers.hh>
#include <numeric/internal/RowVectors.hh>
#include <numeric/internal/RowsPointer.hh>
#include <numeric/random/random.fwd.hh>
#include <ObjexxFCL/Dimension.fwd.hh>
#include <ObjexxFCL/Dimension.hh>
#include <ObjexxFCL/DimensionExpression.hh>
#include <ObjexxFCL/DynamicIndexRange.fwd.hh>
#include <ObjexxFCL/DynamicIndexRange.hh>
#include <ObjexxFCL/FArray.fwd.hh>
#include <ObjexxFCL/FArray.hh>
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/FArray1.hh>
#include <ObjexxFCL/FArray1D.fwd.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2.fwd.hh>
#include <ObjexxFCL/FArray2.hh>
#include <ObjexxFCL/FArray2D.fwd.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArrayInitializer.fwd.hh>
#include <ObjexxFCL/FArrayInitializer.hh>
#include <ObjexxFCL/FArraySection.fwd.hh>
#include <ObjexxFCL/FArraySection.hh>
#include <ObjexxFCL/FArrayTraits.fwd.hh>
#include <ObjexxFCL/FArrayTraits.hh>
#include <ObjexxFCL/IndexRange.fwd.hh>
#include <ObjexxFCL/IndexRange.hh>
#include <ObjexxFCL/InitializerSentinel.hh>
#include <ObjexxFCL/Observer.fwd.hh>
#include <ObjexxFCL/Observer.hh>
#include <ObjexxFCL/ObserverMulti.hh>
#include <ObjexxFCL/ObserverSingle.hh>
#include <ObjexxFCL/ProxySentinel.hh>
#include <ObjexxFCL/SetWrapper.fwd.hh>
#include <ObjexxFCL/Star.fwd.hh>
#include <ObjexxFCL/Star.hh>
#include <ObjexxFCL/proxy_const_assert.hh>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iosfwd>
#include <iostream>
#include <limits>
#include <list>
#include <map>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>
#include <basic/MetricValue.fwd.hh>
#include <basic/Tracer.fwd.hh>
#include <basic/datacache/BasicDataCache.fwd.hh>
#include <boost/bind.hpp>
#include <boost/function.hpp>
#include <boost/functional/hash.hpp>

//Auto Headers


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR("core.scoring.constraints.BackboneStubConstraint");

utility::pointer::owning_ptr< AngleConstraint > BackboneStubConstraint::ang_cst_(0);

BackboneStubConstraint::BackboneStubConstraint(
	pose::Pose const & pose,
	Size const seqpos,
	AtomID const & fixed_atom_id,
	conformation::Residue const & target_rsd,
	core::Real const & superposition_bonus,
	core::Real const & CB_force_constant
):
	Constraint( core::scoring::backbone_stub_constraint ),
	superposition_bonus_( superposition_bonus ),
	CB_force_constant_( CB_force_constant ),
	seqpos_( seqpos ),
	fixed_atom_id_( fixed_atom_id )
{

	// store info about the target residue
	assert( target_rsd.is_protein() );
	if ( (target_rsd.aa() == chemical::aa_gly ) ) {
		TR << "ERROR - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
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
		TR << "ERROR - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
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
		FuncOP cos_func = new PeriodicFunc(0., 1., 1., 0.);
		ang_cst_ = new AngleConstraint( cos_func );
	}
}

core::Size
BackboneStubConstraint::seqpos() const
{ return seqpos_; }

void BackboneStubConstraint::show( std::ostream& out ) const
{
	out << "BackboneStubCst Seqpos: " << seqpos_ << "    bonus: " << superposition_bonus_ << "    CB force constant: " << CB_force_constant_ << std::endl;
}


bool
BackboneStubConstraint::operator == ( Constraint const & other_cst ) const
{

	if( !dynamic_cast< BackboneStubConstraint const * > ( &other_cst ) ) return false;

	BackboneStubConstraint const & other( static_cast< BackboneStubConstraint const & > (other_cst) );

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
void
BackboneStubConstraint::score( XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const {

	if ( weights[ this->score_type() ] == 0 ) return;

	core::conformation::Residue const & curr_rsd = xyz_func.residue(seqpos_);
	assert( curr_rsd.is_protein() );

	// verify that the fixed reference point is in the same place
	core::Vector curr_ref_location = xyz_func(fixed_atom_id_);
	core::Real ref_dist = curr_ref_location.distance_squared( fixed_reference_point_ );
	if ( ref_dist > 1E-8 ) {
		TR << "ERROR - BackboneStubConstraint requires a fixed reference atom, but this atom has moved!!" << std::endl;
		std::exit(1);
	}

	// return a value between superposition_bonus_ (-ve) and zero
	core::Real cst_val(superposition_bonus_);
	assert( cst_val < 0. );

	// apply a harmonic constraint on the CB's
	core::Vector CB_curr;
	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
/// SJF I'm disabling the option for glycine backbone stub constraints b/c it hardly seems useful and is something of a pain
/// to deal with since changing curr_rsd to a Residue const &
/*		if( basic::options::option[basic::options::OptionKeys::hotspot::allow_gly]() ) {
			// make an alanine
			core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
			core::conformation::ResidueOP ala_copy = core::conformation::ResidueFactory::create_residue( alatype );
			// move ala on top of gly
			ala_copy->orient_onto_residue( curr_rsd );
			curr_rsd = *(ala_copy); // gly => ala
		}
		else { */
			TR << "ERROR - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
			return;
//		}
	}
	CB_curr = curr_rsd.xyz("CB");

	core::Real CB_d2 = CB_curr.distance_squared( CB_target_ );
	cst_val += CB_force_constant_ * CB_d2;
 	if ( cst_val > 0. ) return;

///// SJF the following are the dot product results (or cosines) of the various
///// vectors in an amino-acid:
///// ca-cb x ca_c = 0.34
///// ca-cb x cb_c = 0.82
///// ca-cb x n-c = 0.02
///// So, by far the most orthogonal choice is to use ca-b and n-c

	// multiply by the cos of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
//	core::Vector const CB_CA_curr = CB_curr - CA_curr;
	cst_val *= ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_);
	if ( cst_val > -1E-10 ) return;

	// multiply by the cos of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
//	core::Vector const C_N_curr = C_curr - N_curr;

	cst_val *= ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ );
	if ( cst_val > -1E-10 ) return;

	emap[ this->score_type() ] += cst_val;
}

void
BackboneStubConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const & xyz,
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
	core::Real cst_val(superposition_bonus_);
	assert( cst_val < 0. );

	// start by computing the cst value normally, collecting score components along the way. if cst_val is non-negative no derivative is computed

	if ( (curr_rsd.aa() == chemical::aa_gly) ) {
/*		if( basic::options::option[basic::options::OptionKeys::hotspot::allow_gly] ) {
			// make an alanine
			core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
			core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
			core::conformation::ResidueOP ala_copy = core::conformation::ResidueFactory::create_residue( alatype );
			// move ala on top of gly
			ala_copy->orient_onto_residue( *curr_rsd );
			curr_rsd = *ala_copy; // gly => ala
		}
		else { */
			TR << "ERROR - Gly residues cannot be used in BackboneStubConstraints." << std::endl;
			return;
//		}
	}

	// apply a harmonic constraint on the CB's
	core::Vector const CB_curr( curr_rsd.xyz("CB") );

	core::Real const CB_d2 = CB_curr.distance_squared( CB_target_ );
	core::Real const CB_pos_term = CB_force_constant_ * CB_d2;
	cst_val += CB_pos_term;
	if ( cst_val > -1E-10 ) return;

	// multiply by the cos of the angle between the CB-CA and CB'-CA' vectors
	core::Vector const CA_curr = curr_rsd.xyz("CA");
	core::Vector const CB_CA_curr = CB_curr - CA_curr;
	core::Real const CB_CA_angle_term( ang_cst_->score( CB_curr, CA_curr, CA_curr + CB_CA_target_ ) );
	if( CB_CA_angle_term <= 1E-10 ) return;

	// multiply by the cos of the angle between the C-N and C'-N' vectors
	core::Vector const C_curr = curr_rsd.xyz("C");
	core::Vector const N_curr = curr_rsd.xyz("N");
	core::Vector const C_N_curr = C_curr - N_curr;
	core::Real const C_N_angle_term( ang_cst_->score( C_curr, N_curr, N_curr + C_N_target_ ) );
	if ( C_N_angle_term <= 1E-10 ) return;

	core::Real const constant_dist_term( weights[ this->score_type() ] );
	core::Real const constant_ang_term( constant_dist_term  * ( superposition_bonus_ + CB_pos_term ) );

	// contribution from differentiating CB_dist ( * CA_angle_term * C_angle_term from product rule)
	// and then adding the effect of the angular constraint on cb using the chain rule
	if ( atom == CB_atom_id_ ) {
// the effects of the coordinate constraint
		Vector const CB_f2( CB_curr - CB_target_ );
		core::Real const CB_dist( CB_f2.length() );
		core::Real const CB_deriv = 2. * CB_force_constant_ * CB_dist;
		if ( CB_deriv != 0.0 && CB_dist != 0.0 ) {
			Vector const CB_f1( CB_curr.cross( CB_target_ ) );
			F1 += ( ( CB_deriv / CB_dist ) * CB_f1 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
			F2 += ( ( CB_deriv / CB_dist ) * CB_f2 ) * CB_CA_angle_term * C_N_angle_term * constant_dist_term;
		}

// the angular constraint on cb
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CB_curr/*p1*/, CA_curr/*p2*/, CA_curr + CB_CA_target_/*p3*/, partial_F1, partial_F2 );
		F1 += partial_F1 * C_N_angle_term * constant_ang_term;
		F2 += partial_F2 * C_N_angle_term * constant_ang_term;
		return;
	}

// the angular constraint on ca
// See N atom for explanation on the derivs
	if ( atom == CA_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( CA_curr, CA_curr - CB_CA_curr, CA_curr + CB_CA_target_ - CB_CA_curr, partial_F1, partial_F2 );
		F1 += -partial_F1 * C_N_angle_term * constant_ang_term;
		F2 += -partial_F2 * C_N_angle_term * constant_ang_term;
		return;
	}

// the angular constraint on c
	if ( atom == C_atom_id_ ) {
		Vector partial_F1(0.), partial_F2(0.);
		ang_cst_->p1_deriv( C_curr/*p1*/, N_curr/*p2*/, N_curr + C_N_target_/*p3*/, partial_F1, partial_F2 );
		F1 += partial_F1 * CB_CA_angle_term * constant_ang_term;
		F2 += partial_F2 * CB_CA_angle_term * constant_ang_term;
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
		F1 += -partial_F1 * CB_CA_angle_term * constant_ang_term;
		F2 += -partial_F2 * CB_CA_angle_term * constant_ang_term;
	}
	return;
}

ConstraintOP BackboneStubConstraint::clone() const
{
	return ConstraintOP( new BackboneStubConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP BackboneStubConstraint::remapped_clone( pose::Pose const& /*src*/, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {

	core::Size new_seqpos = seqpos_;
	AtomID new_fixed_atom_id = fixed_atom_id_;

	if ( smap ) {
		new_seqpos = (*smap)[ seqpos_ ];
		new_fixed_atom_id.rsd() = (*smap)[ fixed_atom_id_.rsd() ];
		if( new_seqpos == 0 ) return NULL;
	}

	// make an alanine
	core::chemical::ResidueTypeSetCAP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::chemical::ResidueType const & alatype( residue_set->name_map( "ALA" ) );
	core::conformation::ResidueOP ala = core::conformation::ResidueFactory::create_residue( alatype );
	ala->set_xyz("CB",CB_target_);
	ala->set_xyz("CA",CA_target_);
	ala->set_xyz("C",C_target_);
	ala->set_xyz("N",N_target_);

	return new BackboneStubConstraint(dest, new_seqpos, new_fixed_atom_id, *ala, superposition_bonus_, CB_force_constant_ );
}


/*ConstraintOP
BackboneStubConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
	for(
  if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 && seqmap[atom3_.rsd()] != 0 && seqmap[atom4_.rsd()] != 0 ) {
    AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
      remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] ),
			remap_a3( atom3_.atomno(), seqmap[atom3_.rsd()] ),
			remap_a4( atom4_.atomno(), seqmap[atom4_.rsd()] );
    return ConstraintOP( new DihedralConstraint( remap_a1, remap_a2, remap_a3, remap_a4, this->func_ ) );
  } else {
    return NULL;
  }
}
*/


} // namespace constraints
} // namespace scoring
} // namespace core
