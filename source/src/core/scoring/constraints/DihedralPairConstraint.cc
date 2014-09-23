// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/DihedralPairConstraint.cc
/// @brief Restrain a pair of residues to take the same torsion angles
/// @author Frank DiMaio

#include <core/scoring/constraints/DihedralPairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FuncFactory.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/dihedral_deriv.hh>

#include <utility/exit.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {
namespace constraints {

static thread_local basic::Tracer TR( "core.io.constraints" );


/////////////////////////////////////////////////////////////////////////////

ConstraintOP
DihedralPairConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
  if (   seqmap[atomA1_.rsd()] != 0 && seqmap[atomA2_.rsd()] != 0
	    && seqmap[atomA3_.rsd()] != 0 && seqmap[atomA4_.rsd()] != 0
	    && seqmap[atomB1_.rsd()] != 0 && seqmap[atomB2_.rsd()] != 0
	    && seqmap[atomB3_.rsd()] != 0 && seqmap[atomB4_.rsd()] != 0 ) {
    AtomID remap_a1( atomA1_.atomno(), seqmap[atomA1_.rsd()] ),
           remap_a2( atomA2_.atomno(), seqmap[atomA2_.rsd()] ),
		       remap_a3( atomA3_.atomno(), seqmap[atomA3_.rsd()] ),
		       remap_a4( atomA4_.atomno(), seqmap[atomA4_.rsd()] );
    AtomID remap_b1( atomB1_.atomno(), seqmap[atomB1_.rsd()] ),
           remap_b2( atomB2_.atomno(), seqmap[atomB2_.rsd()] ),
		       remap_b3( atomB3_.atomno(), seqmap[atomB3_.rsd()] ),
		       remap_b4( atomB4_.atomno(), seqmap[atomB4_.rsd()] );
    return ConstraintOP( new DihedralPairConstraint(
				remap_a1, remap_a2, remap_a3, remap_a4,
				remap_b1, remap_b2, remap_b3, remap_b4,
				this->func_ ) );
  } else {
    return NULL;
  }
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP DihedralPairConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
  id::NamedAtomID atomA1( core::pose::atom_id_to_named_atom_id(atom(1), src ) );
  id::NamedAtomID atomA2( core::pose::atom_id_to_named_atom_id(atom(2), src ) );
  id::NamedAtomID atomA3( core::pose::atom_id_to_named_atom_id(atom(3), src ) );
  id::NamedAtomID atomA4( core::pose::atom_id_to_named_atom_id(atom(4), src ) );
  id::NamedAtomID atomB1( core::pose::atom_id_to_named_atom_id(atom(5), src ) );
  id::NamedAtomID atomB2( core::pose::atom_id_to_named_atom_id(atom(6), src ) );
  id::NamedAtomID atomB3( core::pose::atom_id_to_named_atom_id(atom(7), src ) );
  id::NamedAtomID atomB4( core::pose::atom_id_to_named_atom_id(atom(8), src ) );
  if ( smap ) {
    atomA1.rsd() = (*smap)[ atomA1_.rsd() ];
    atomA2.rsd() = (*smap)[ atomA2_.rsd() ];
    atomA3.rsd() = (*smap)[ atomA3_.rsd() ];
    atomA4.rsd() = (*smap)[ atomA4_.rsd() ];
    atomB1.rsd() = (*smap)[ atomB1_.rsd() ];
    atomB2.rsd() = (*smap)[ atomB2_.rsd() ];
    atomB3.rsd() = (*smap)[ atomB3_.rsd() ];
    atomB4.rsd() = (*smap)[ atomB4_.rsd() ];
  }

  //get AtomIDs for target pose
  id::AtomID id1( core::pose::named_atom_id_to_atom_id(atomA1, dest ));
  id::AtomID id2( core::pose::named_atom_id_to_atom_id(atomA2, dest ));
  id::AtomID id3( core::pose::named_atom_id_to_atom_id(atomA3, dest ));
  id::AtomID id4( core::pose::named_atom_id_to_atom_id(atomA4, dest ));
  id::AtomID id5( core::pose::named_atom_id_to_atom_id(atomB1, dest ));
  id::AtomID id6( core::pose::named_atom_id_to_atom_id(atomB2, dest ));
  id::AtomID id7( core::pose::named_atom_id_to_atom_id(atomB3, dest ));
  id::AtomID id8( core::pose::named_atom_id_to_atom_id(atomB4, dest ));
  if (    id1.valid() && id2.valid() &&  id3.valid() && id4.valid()
	     && id5.valid() && id6.valid() &&  id7.valid() && id8.valid() ) {
    return ConstraintOP( new DihedralPairConstraint( id1, id2, id3, id4, id5, id6, id7, id8, func_, score_type() ) );
  } else {
    return NULL;
  }
}


id::AtomID const &
DihedralPairConstraint::atom( Size const n ) const {
	switch( n ) {
	case 1: return atomA1_;
	case 2: return atomA2_;
	case 3: return atomA3_;
	case 4: return atomA4_;
	case 5: return atomB1_;
	case 6: return atomB2_;
	case 7: return atomB3_;
	case 8: return atomB4_;
	default:
		utility_exit_with_message( "DihedralPairConstraint::atom() bad argument" );
	}
	return atomA1_;
}



////////////////////////////////////////////////////////////////////////////////////////////////////
///@details one line definition "Dihedral atom1 res1 atom2 res2 atom3 res3 atom4 res4 function_type function_definition"
void
DihedralPairConstraint::read_def(
	std::istream & in,
	pose::Pose const & pose,
	func::FuncFactory const & func_factory
) {
	Size res1, res2, res3, res4, res5, res6, res7, res8;
	std::string tempres1, tempres2, tempres3, tempres4, tempres5, tempres6, tempres7, tempres8;
	std::string name1, name2, name3, name4, name5, name6, name7, name8;
	std::string func_type;
	std::string type;

	in
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> name4 >> tempres4
		>> name5 >> tempres5
		>> name6 >> tempres6
		>> name7 >> tempres7
		>> name8 >> tempres8
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );
	ConstraintIO::parse_residue( pose, tempres4, res4 );
	ConstraintIO::parse_residue( pose, tempres5, res5 );
	ConstraintIO::parse_residue( pose, tempres6, res6 );
	ConstraintIO::parse_residue( pose, tempres7, res7 );
	ConstraintIO::parse_residue( pose, tempres8, res8 );

	TR.Debug 	<< "read: " << name1 << " " << name2 << " " << name3 << " " << name4 << " "
						<< res1 << " " << res2 << " " << res3 << " " << res4 << " func: " << func_type
						<< std::endl;
	if (    res1 > pose.total_residue() || res2 > pose.total_residue()
	     || res3 > pose.total_residue() || res4 > pose.total_residue()
	     || res5 > pose.total_residue() || res6 > pose.total_residue()
	     || res7 > pose.total_residue() || res8 > pose.total_residue() ) {
		TR.Warning 	<< "ignored constraint (no such atom in pose!)"
								<< name1 << " " << name2 << " " << name3 << " " << name4 << " "
													<< res1 << " " << res2 << " " << res3 << " " << res4 << " func: " << func_type
													<< std::endl;
		in.setstate( std::ios_base::failbit );
		return;
	}

	atomA1_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose ) );
	atomA2_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose ) );
	atomA3_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name3, res3 ), pose ) );
	atomA4_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name4, res4 ), pose ) );
	atomB1_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name5, res5 ), pose ) );
	atomB2_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name6, res6 ), pose ) );
	atomB3_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name7, res7 ), pose ) );
	atomB4_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name8, res8 ), pose ) );

	if (    atomA1_.atomno() == 0 || atomA2_.atomno() == 0 || atomA3_.atomno() == 0 || atomA4_.atomno() == 0
	     || atomB1_.atomno() == 0 || atomB2_.atomno() == 0 || atomB3_.atomno() == 0 || atomB4_.atomno() == 0 ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "," << name4 << ","
			<< name5 << "," << name6 << "," << name7 << "," << name8 << "), "
			<< "and found AtomIDs (" << atomA1_ << "," << atomA2_ << "," << atomA3_ << "," << atomA4_ << " -- "
			<< atomB1_ << "," << atomB2_ << "," << atomB3_ << "," << atomB4_ << ")" << std::endl;
		in.setstate( std::ios_base::failbit );
		return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( in );

	while( in.good() && (in.get() != '\n') ) {}

	if ( TR.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////

Real
DihedralPairConstraint::score(
	conformation::Conformation const & conformation
) const {
	return score(
		conformation.xyz( atomA1_ ), conformation.xyz( atomA2_ ),
		conformation.xyz( atomA3_ ), conformation.xyz( atomA4_ ),
		conformation.xyz( atomB1_ ), conformation.xyz( atomB2_ ),
		conformation.xyz( atomB3_ ), conformation.xyz( atomB4_ ));
}

void
DihedralPairConstraint::score(
	func::XYZ_Func const & xyz,
	EnergyMap const &,
	EnergyMap & emap
) const {
	emap[ this->score_type() ] += score(
		xyz( atomA1_ ), xyz( atomA2_ ), xyz( atomA3_ ), xyz( atomA4_ ),
		xyz( atomB1_ ), xyz( atomB2_ ), xyz( atomB3_ ), xyz( atomB4_ )
	);
}


Real
DihedralPairConstraint::score(
	Vector const & p1, Vector const & p2, Vector const & p3, Vector const & p4,
	Vector const & p5, Vector const & p6, Vector const & p7, Vector const & p8
) const {
	core::Real difference =  dihedral_degrees( p1, p2, p3, p4 ) - dihedral_degrees( p5, p6, p7, p8 ) ;
	if (difference > 180)  difference -=360;
	if (difference < -180) difference +=360;
	TR.Debug << "difference of " << difference << " degrees" << std::endl;
	return func(difference);
}


void
DihedralPairConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
 	Vector & F2,
 	EnergyMap const & weights
) const {
	using namespace numeric::deriv;

	Vector f1(0.0) ,f2(0.0);

	Real theta(0.0), theta0;
	if ( atom == atomA1_ ) {
		dihedral_p1_cosine_deriv( xyz( atomA1_ ),xyz( atomA2_ ), xyz( atomA3_ ), xyz( atomA4_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomB1_ ),xyz( atomB2_ ), xyz( atomB3_ ), xyz( atomB4_ ) );
	} else if ( atom == atomA2_ ) {
		dihedral_p2_cosine_deriv( xyz( atomA1_ ),xyz( atomA2_ ), xyz( atomA3_ ), xyz( atomA4_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomB1_ ),xyz( atomB2_ ), xyz( atomB3_ ), xyz( atomB4_ ) );
	} else if ( atom == atomA3_ ) {
		dihedral_p2_cosine_deriv( xyz( atomA4_ ), xyz( atomA3_ ), xyz( atomA2_ ), xyz( atomA1_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomB4_ ), xyz( atomB3_ ), xyz( atomB2_ ), xyz( atomB1_ ) );
	} else if ( atom == atomA4_ ) {
		dihedral_p1_cosine_deriv( xyz( atomA4_ ), xyz( atomA3_ ), xyz( atomA2_ ), xyz( atomA1_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomB4_ ), xyz( atomB3_ ), xyz( atomB2_ ), xyz( atomB1_ ) );
	} else if ( atom == atomB1_ ) {
		dihedral_p1_cosine_deriv( xyz( atomB1_ ),xyz( atomB2_ ), xyz( atomB3_ ), xyz( atomB4_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomA1_ ),xyz( atomA2_ ), xyz( atomA3_ ), xyz( atomA4_ ) );
	} else if ( atom == atomB2_ ) {
		dihedral_p2_cosine_deriv( xyz( atomB1_ ),xyz( atomB2_ ), xyz( atomB3_ ), xyz( atomB4_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomA1_ ),xyz( atomA2_ ), xyz( atomA3_ ), xyz( atomA4_ ) );
	} else if ( atom == atomB3_ ) {
		dihedral_p2_cosine_deriv( xyz( atomB4_ ), xyz( atomB3_ ), xyz( atomB2_ ), xyz( atomB1_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomA4_ ), xyz( atomA3_ ), xyz( atomA2_ ), xyz( atomA1_ ) );
	} else if ( atom == atomB4_ ) {
		dihedral_p1_cosine_deriv( xyz( atomB4_ ), xyz( atomB3_ ), xyz( atomB2_ ), xyz( atomB1_ ), theta, f1, f2 );
		theta0 = dihedral_degrees( xyz( atomA4_ ), xyz( atomA3_ ), xyz( atomA2_ ), xyz( atomA1_ ) );
	} else {
		return;
	}

	core::Real difference =  numeric::conversions::degrees( theta ) - theta0;
	if (difference > 180)  difference -=360;
	if (difference < -180) difference +=360;
	Real const dE_dtheta( numeric::conversions::degrees( dfunc( difference ) ) );

	F1 += dE_dtheta * weights[ this->score_type() ] * f1;
	F2 += dE_dtheta * weights[ this->score_type() ] * f2;
}


bool
DihedralPairConstraint::operator == ( Constraint const & other_cst ) const {
	if( !dynamic_cast< DihedralPairConstraint const * > ( &other_cst ) ) return false;

	DihedralPairConstraint const & other( static_cast< DihedralPairConstraint const & > (other_cst) );

	if( atomA1_ != other.atomA1_ ) return false;
	if( atomA2_ != other.atomA2_ ) return false;
	if( atomA3_ != other.atomA3_ ) return false;
	if( atomA4_ != other.atomA4_ ) return false;
	if( atomB1_ != other.atomB1_ ) return false;
	if( atomB2_ != other.atomB2_ ) return false;
	if( atomB3_ != other.atomB3_ ) return false;
	if( atomB4_ != other.atomB4_ ) return false;
	if( func_ != other.func_ ) return false;
	if( this->score_type() != other.score_type() ) return false;

	return true;
}


void DihedralPairConstraint::show( std::ostream & out ) const {
	out << "DihedralPairConstraint";
	for ( Size i = 1; i <= natoms(); ++i ) {
		AtomID const & id = atom(i);
		out << ' ' << id.rsd() << ' ' << id.atomno();
	}
	out << ' ';
	func_->show_definition(out);
}


Size DihedralPairConstraint::show_violations(
        std::ostream& out,
        pose::Pose const& pose,
        Size verbose_level,
        Real threshold
) const {
	     if ( verbose_level > 80 ) {
                out << "DihedralPair ("
                        << pose.residue_type(atomA2_.rsd() ).atom_name( atomA2_.atomno() ) << ":"
                        << atomA2_.atomno() << "," << atomA2_.rsd() << "-"
                        << pose.residue_type(atomA3_.rsd() ).atom_name( atomA3_.atomno() ) << ":"
                        << atomA3_.atomno() << "," << atomA3_.rsd() << " === "
                        << pose.residue_type(atomB2_.rsd() ).atom_name( atomB2_.atomno() ) << ":"
                        << atomB2_.atomno() << "," << atomB2_.rsd() << "-"
                        << pose.residue_type(atomB3_.rsd() ).atom_name( atomB3_.atomno() ) << ":"
                        << atomB3_.atomno() << "," << atomB3_.rsd() << "-";
        }
        return func_->show_violations( out, 0.0, verbose_level, threshold );
}


Real
DihedralPairConstraint::func( Real const theta ) const {
	return func_->func( theta );
}

Real
DihedralPairConstraint::dfunc( Real const theta ) const {
	return func_->dfunc( theta );
}

} // constraints
} // scoring
} // core
