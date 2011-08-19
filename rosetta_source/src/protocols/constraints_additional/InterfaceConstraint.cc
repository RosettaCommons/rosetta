// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/constraints_additional/InterfaceConstraint.cc
///
/// @brief
/// @author Monica Berrondo

// Unit Headers
#include <protocols/constraints_additional/InterfaceConstraint.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/FuncFactory.hh>
#include <protocols/scoring/InterfaceInfo.hh>
#include <protocols/scoring/InterchainPotential.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <core/pose/util.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
// Utility Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>


static basic::Tracer tr("protocols.constraints_additional");


namespace protocols {
namespace constraints_additional {

using namespace core;
using namespace core::scoring;
using namespace core::scoring::constraints;

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP InterfaceConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( atom_id_to_named_atom_id( atom(1), src ) );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1_.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( named_atom_id_to_atom_id( atom1, dest ) );
	if ( id1.valid() ) {
		return new InterfaceConstraint( id1, func_, score_type() );
	} else {
		return NULL;
	}
}

void
InterfaceConstraint::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	// calculate the interface and set the second atom to the correct residue
	conformation::ResidueCOP src_rsd = new conformation::Residue( pose.residue( atom1_.rsd() ) );
	conformation::ResidueCOP int_rsd;

	protocols::scoring::InterfaceInfo & interface( nonconst_interface_from_pose( pose ) );
	interface.initialize();
	interface.calculate( pose );

	int_rsd = new conformation::Residue( pose.residue( interface.closest_interface_residue( pose, src_rsd->seqpos(), dist_ ) ) );

	atom2_ = id::AtomID( int_rsd->nbr_atom(), int_rsd->seqpos() );

	f1_ = src_rsd->xyz( src_rsd->nbr_atom() ).cross( int_rsd->xyz( int_rsd->nbr_atom() ) );
	f2_ = int_rsd->xyz( int_rsd->nbr_atom() ) - src_rsd->xyz( src_rsd->nbr_atom() );

}

void InterfaceConstraint::show( std::ostream& out ) const {
	out << "InterfaceConstraint ("
			<< atom1_.atomno() << "," << atom1_.rsd() << ")" << std::endl;
	func_->show( out );
}

void InterfaceConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " " << atom_id_to_named_atom_id( atom1_, pose )  << " ";
	func_->show_definition( out );
}

// atom deriv
void
InterfaceConstraint::fill_f1_f2(
	AtomID const & atom,
	XYZ_Func const &,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	AtomID other_atom;
	if ( atom == atom1_ ) other_atom = atom2_;
	else if ( atom == atom2_ ) other_atom = atom1_;
	else {
		// std::cout << "Error in AtomPairConstraint::fill_f1_f2()" << std::endl;
		// std::cout << "Bad AtomID: (" << atom.rsd() << ", " << atom.atomno() << ") -- options: (";
		// std::cout << atom1_.rsd() << ", " << atom1_.atomno() << ") and (";
		// std::cout << atom2_.rsd() << ", " << atom2_.atomno() << ");" << std::endl;
		return;
	}
	Real wderiv( weights[ this->score_type() ] * func_->dfunc( dist_ ));
	F1 += wderiv * f1_;
	F2 += wderiv * f2_;

}

ConstraintOP
InterfaceConstraint::remap_resid( id::SequenceMapping const &seqmap ) const
{
  if ( seqmap[atom1_.rsd()] != 0 ) {
    AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] );
    return ConstraintOP( new InterfaceConstraint( remap_a1, this->func_ ) );
  } else {
    return NULL;
  }
}


///@details one line definition "Interfaces atom1 res1 partners function_type function_definition"
void
InterfaceConstraint::read_def(
	std::istream & data,
	pose::Pose const & pose,
	constraints::FuncFactory const & func_factory
) {
	Size res1, chainid;
	std::string tempres1;
	std::string name1;
	std::string partners;
	std::string func_type;
	std::string type;

	tr << "reading interface constraints" << std::endl;
	data
		>> name1 >> tempres1
		>> partners
		>> func_type;

	constraints::ConstraintIO::parse_residue( pose, tempres1, res1 );

	tr.Debug << "read: " << name1 << " " << res1 << " func: " << func_type << std::endl;
	if ( res1 > pose.total_residue() ) {
		tr.Warning 	<< "ignored constraint (no such atom in pose!)"
			<< name1 << " "  << res1 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	chainid = pose.residue( res1 ).chain();

	if ( name1 == "H" && res1 == 1 && pose.is_fullatom() ) name1 = "1H";

	atom1_ = pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose );
	// Make a dummy atom2 that will be modified during scoring to the opposite interface atom
	// must be different from atom1_, otherwise ConstraintSet will think it's a one-body energy
	atom2_ = pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, 1 ), pose );

	if ( atom1_.atomno() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "), "
			<< "and found AtomIDs (" << atom1_ << ")" << std::endl;
			data.setstate( std::ios_base::failbit );
			runtime_assert( false );
			return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
	//chu skip the rest of line since this is a single line defintion.
		while( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( tr.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}

	runtime_assert( atom1_.valid() );
	runtime_assert( atom2_.valid() );
} // read_def



/// @details Pose must already contain a Interface object or this method will fail
protocols::scoring::InterfaceInfo const &
InterfaceConstraint::interface_from_pose( core::pose::Pose const & pose ) const
{
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;
	return *( static_cast< protocols::scoring::InterfaceInfo const * >(pose.data().get_const_ptr( INTERFACE_INFO )() ) );
}

/// @details Either returns a non-const reference to the Interface object that already exists
/// in the pose, or creates a new Interface object, places it in the pose, and then returns
/// a non-const reference to it
protocols::scoring::InterfaceInfo &
InterfaceConstraint::nonconst_interface_from_pose( core::pose::Pose & pose ) const
{
	//using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

	if ( pose.data().has( INTERFACE_INFO ) ) {
		return *( static_cast< protocols::scoring::InterfaceInfo * >(pose.data().get_ptr( INTERFACE_INFO )() ) );
	}
	// else
	protocols::scoring::InterfaceInfoOP interface = new protocols::scoring::InterfaceInfo;
	pose.data().set( INTERFACE_INFO, interface );
	return *interface;
}

} // constraints
} // protocols
