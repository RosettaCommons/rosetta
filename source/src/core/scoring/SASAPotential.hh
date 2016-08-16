// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/SASAPotential.hh
/// @brief
/// @author


#ifndef INCLUDED_core_scoring_SASAPotential_hh
#define INCLUDED_core_scoring_SASAPotential_hh

#include <core/scoring/SASAPotential.fwd.hh>

#include <core/scoring/power_diagram/PowerDiagram.hh>

#include <core/types.hh>

#include <core/scoring/DerivVectorPair.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/DomainMap.fwd.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#if defined(WIN32)  // required, because Windows uses rad1 and rad2 as keywords
#undef rad1
#undef rad2
#endif

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

namespace core {
namespace scoring {

typedef numeric::xyzVector< core::Real > Vector;
typedef numeric::xyzMatrix< core::Real > Matrix;

////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////



class PDatom : public utility::pointer::ReferenceCount {

public:
	PDatom() : res_(1), atom_(1), rad_(0.0), rad2_(0.0) {}
	Size const & res() const { return res_; }
	Size & nonconst_res() { return res_; }
	Size const & atom() const { return atom_; }
	Size & nonconst_atom() { return atom_; }
	Vector const & xyz() const { return xyz_; }
	Vector & nonconst_xyz() { return xyz_; }
	Real const & rad() const { return rad_; }
	Real & nonconst_rad() { return rad_; }
	Real const & rad2() const { return rad2_; }
	Real & nonconst_rad2() { return rad2_; }
	std::list< PDvertexOP > & vertices(){ return cell_vertices_; }

private:
	Size res_;
	Size atom_;
	Real rad_;
	Real rad2_;
	Vector xyz_;
	std::list< PDvertexOP > cell_vertices_;

};


class PDvertex : public utility::pointer::ReferenceCount {

public:
	PDvertex() : id_(1), finite_(false), live_(false), power_(0.0) {}
	Size & nonconst_id() { return id_; }
	Size const & id() const { return id_; }
	bool live() const { return live_; }
	void set_live( bool in_val ) { live_ = in_val; }
	bool finite() const { return finite_; }
	void set_finite( bool in_val ) { finite_ = in_val; }
	Vector const & xyz() const { return xyz_; }
	Vector & nonconst_xyz() { return xyz_; }
	Vector const & direction() const { return xyz_; }
	Vector & nonconst_direction() { return xyz_; }
	Real power() { return power_; }
	Real & nonconst_power() { return power_; }
	utility::vector1< PDvertexOP >  const & partners() const { return partners_; }
	utility::vector1< PDvertexOP > & nonconst_partners() { return partners_; }
	utility::vector1< PDatomOP >  const & generators() const { return generators_; }
	utility::vector1< PDatomOP > & nonconst_generators() { return generators_; }

private:
	Size id_;
	bool finite_;
	bool live_;
	Vector xyz_;
	Real power_;
	utility::vector1< PDvertexOP > partners_;
	utility::vector1< PDatomOP > generators_;

};


// Class for the 'vertices' that correspond to triple atom
// intersections.
class PDinter : public utility::pointer::ReferenceCount {

public:
	PDinter( Vector & xyz, PDvertexCOP & v1, PDvertexCOP & v2 ) : xyz_(xyz), v1_(v1), v2_(v2), circle_(false) {}

	Vector const & xyz() const { return xyz_; }
	PDvertexCOP const & vrt1() const { return v1_; }
	PDvertexCOP const & vrt2() const { return v2_; }
	utility::vector1< PDatomCOP > const & atoms() const { return atoms_; }
	void add_atom( PDatomCOP pa ) { atoms_.push_back( pa ); }
	bool & circle() { return circle_; }

private:
	Vector xyz_;
	PDvertexCOP v1_;
	PDvertexCOP v2_;
	utility::vector1< PDatomCOP > atoms_;
	bool circle_;

};


// Class for each node in a contour around a buried surface area patch
class SAnode : public utility::pointer::ReferenceCount {

public:
	SAnode( PDinterOP & inter, PDatomCOP & other_atom, Real ct, Real phi ) : inter_(inter), other_atom_(other_atom), cos_theta_(ct), phi_(phi) {}
	SAnode( PDinterOP & inter, PDatomCOP & other_atom ) : inter_(inter), other_atom_(other_atom), cos_theta_(0.0), phi_(0.0) {}

	Vector const & xyz() const { return inter_->xyz(); }
	PDinterOP const & inter() const { return inter_; }
	PDatomCOP const & other_atom() const { return other_atom_; }
	Real cos_theta() const { return cos_theta_; }
	Real phi() const { return phi_; }
	void set_cos_theta( Real ct ) { cos_theta_ = ct; }
	void set_phi( Real p ) { phi_ = p; }

private:
	PDinterOP inter_;
	PDatomCOP other_atom_;
	Real cos_theta_;
	Real phi_;

};

class SASAPotential : public utility::pointer::ReferenceCount {
public:
	typedef core::conformation::Residue Residue;

public:
	/// ctor
	SASAPotential()
	{ pd_ = power_diagram::PowerDiagramOP( new power_diagram::PowerDiagram ); }

	///
	void
	setup_for_scoring(
		pose::Pose & pose
	) const;

	///
	Real
	get_res_res_sasa(
		Residue const & rsd1,
		Residue const & rsd2
	) const;

	void
	eval_residue_pair_derivatives(
		conformation::Residue const & rsd1,
		conformation::Residue const & rsd2,
		pose::Pose const & pose, // provides context
		Real const & factor,
		utility::vector1< DerivVectorPair > & r1_atom_derivs,
		utility::vector1< DerivVectorPair > & r2_atom_derivs
	) const;

	Size & vertex_count() { return vertex_count_; }

public:
private:
	Size vertex_count_;
	power_diagram::PowerDiagramOP pd_;

};



} // scoring
} // core

#endif
