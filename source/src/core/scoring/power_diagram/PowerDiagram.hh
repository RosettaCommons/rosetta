// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   scoring/power_diagram/PowerDiagram.hh
/// @brief  Class for PowerDiagram, including consruction and access
/// @author Jim Havranek

#ifndef INCLUDED_core_scoring_power_diagram_PowerDiagram_hh
#define INCLUDED_core_scoring_power_diagram_PowerDiagram_hh

#include <core/scoring/power_diagram/PowerDiagram.fwd.hh>

#include <numeric/numeric.functions.hh>
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// unit headers
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <list>

#if defined(WIN32)  // required, because Windows uses rad1 and rad2 as keywords
#undef rad1
#undef rad2
#endif

namespace core {
namespace scoring {
namespace power_diagram {

typedef numeric::xyzVector< core::Real > Vector;
typedef numeric::xyzMatrix< core::Real > Matrix;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Class for the 'vertices' that correspond to triple atom
// intersections.

class PDinter : public utility::pointer::ReferenceCount {

public:
	PDinter( Vector & xyz, PDvertex const * v1, PDvertex const * v2 ) : xyz_(xyz), v1_(v1), v2_(v2), circle_(false) {}

	Vector const & xyz() const { return xyz_; }
	PDvertex const * vrt1() const { return v1_; }
	PDvertex const * vrt2() const { return v2_; }
	utility::vector1< PDsphere const * > const & atoms() const { return atoms_; }
	utility::vector1< PDsphere const * > & nonconst_atoms() { return atoms_; }
	void add_atom( PDsphere const * pa ) { atoms_.push_back( pa ); }
	bool const & circle() const { return circle_; }
	bool & nonconst_circle() { return circle_; }

private:
	Vector xyz_;
	PDvertex const * v1_;
	PDvertex const * v2_;
	utility::vector1< PDsphere const * > atoms_;
	bool circle_;
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Class for each node in a contour around a buried surface area patch
class SAnode : public utility::pointer::ReferenceCount {

public:
	SAnode( PDinterOP & inter, PDsphere const *other_atom, core::Real ct, core::Real phi ) : inter_(inter), other_atom_(other_atom), cos_theta_(ct), phi_(phi) {}
	SAnode( PDinterOP & inter, PDsphere const * other_atom ) : inter_(inter), other_atom_(other_atom), cos_theta_(0.0), phi_(0.0) {}

	Vector const & xyz() const { return inter_->xyz(); }
	PDinterOP const & inter() const { return inter_; }
	PDsphere const * other_atom() const { return other_atom_; }
	core::Real cos_theta() const { return cos_theta_; }
	core::Real phi() const { return phi_; }
	void set_cos_theta( core::Real ct ) { cos_theta_ = ct; }
	void set_phi( core::Real p ) { phi_ = p; }

private:
	PDinterOP inter_;
	PDsphere const * other_atom_;
	core::Real cos_theta_;
	core::Real phi_;
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Class for the basic sphere (aka an atom)

class PDsphere : public utility::pointer::ReferenceCount {

public:
	explicit PDsphere() : res_(1), atom_(1), rad_(0.0), rad2_(0.0)
	{
		xyz_.clear();
		cell_vertices_.clear();
		cycles_.clear();
	}
	~PDsphere() { cell_vertices_.clear(); }
	Size const & res() const { return res_; }
	Size & nonconst_res() { return res_; }
	Size const & atom() const { return atom_; }
	Size & nonconst_atom() { return atom_; }
	Vector const & xyz() const { return xyz_; }
	Vector & nonconst_xyz() { return xyz_; }
	core::Real const & rad() const { return rad_; }
	core::Real & nonconst_rad() { return rad_; }
	core::Real const & rad2() const { return rad2_; }
	core::Real & nonconst_rad2() { return rad2_; }
	std::list< PDvertex * > & vertices(){ return cell_vertices_; }
	utility::vector1< utility::vector1< SAnode > > & cycles() { return cycles_; }

private:
	Size res_;
	Size atom_;
	core::Real rad_;
	core::Real rad2_;
	Vector xyz_;
	std::list< PDvertex * > cell_vertices_;
	utility::vector1< utility::vector1< SAnode > > cycles_;
};

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////

// Class for each vertex in the power diagram

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
	core::Real power() { return power_; }
	core::Real & nonconst_power() { return power_; }
	utility::vector1< PDvertex * >  const & partners() const { return partners_; }
	utility::vector1< PDvertex * > & nonconst_partners() { return partners_; }
	utility::vector1< PDsphere * > const & generators() const { return generators_; }
	utility::vector1< PDsphere * > & nonconst_generators() { return generators_; }
	std::list< PDvertexOP >::iterator & my_itr() { return my_itr_; }

private:
	Size id_;
	bool finite_;
	bool live_;
	Vector xyz_;
	core::Real power_;
	utility::vector1< PDvertex * > partners_;
	utility::vector1< PDsphere * > generators_;
	std::list< PDvertexOP >::iterator my_itr_;
};


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////


class PowerDiagram {

public: // construct/destruct

	/// @brief default constructor
	/// @warning no initialization of corners for speed, make sure you reset()
	/// @warning or otherwise set corners before adding points
	inline
	PowerDiagram() {}

	/// @brief power diagram from pose
	PowerDiagram( core::pose::Pose & pose );

	/// @brief copy constructor
	inline
	PowerDiagram(
		PowerDiagram const & // pd
	)
	{
		std::cout << "PD copy constructor called?" << std::endl;
		spheres_.clear();
		finite_vertices_.clear();
		infinite_vertices_.clear();
		sphere_lookup_.clear();
	}

	/// @brief default destructor
	inline
	~PowerDiagram() {}


public: // assignment

	/// @brief copy assignment
	inline
	PowerDiagram &
	operator =(
		PowerDiagram const & pd
	)
	{
		// Copy these over
		this->spheres_ = pd.spheres_;
		this->finite_vertices_ = pd.finite_vertices_;
		this->infinite_vertices_ = pd.infinite_vertices_;
		this->sphere_lookup_ = pd.sphere_lookup_;
		return *this;
	}


public:

	PDsphereOP
	make_new_sphere( core::pose::Pose & p, Size ires, Size iatm );

	PDsphereOP
	make_new_sphere( Vector & pos, Real rad );

	void
	construct_from_pose( core::pose::Pose & p );

	void
	make_initial_power_diagram();

	void
	add_single_atom_to_power_diagram( PDsphereOP & new_atm );

	void
	store_new_sphere( PDsphereOP new_sph ) {
		spheres_.push_back( new_sph );
		sphere_lookup_[ new_sph->res() ][ new_sph->atom() ] = new_sph.get();
	}

	std::list< PDinterOP > get_intersections_for_atom( Size ires, Size iatm );

	std::list< PDinterOP > get_intersections_for_atom( PDsphere* patom );

	PDsphere* sphere_lookup( Size ires, Size iatm  )
	{ return sphere_lookup_[ ires ][ iatm ]; }

	core::Real extract_sasa_for_atom( Size ires, Size iatm );

	/// @brief Clear all data from power diagram
	void clear();

public: // query constructed power diagram
	Size const & vertex_count() { return vertex_count_; };

private: // data
	Size vertex_count_;
	utility::vector1< PDsphereOP > spheres_;
	std::list< PDvertexOP > finite_vertices_;
	std::list< PDvertexOP > infinite_vertices_;
	utility::vector1< utility::vector1< PDsphere* > > sphere_lookup_;

};

// Utility functions
core::Real power_distance( Vector const & pt, PDsphere const * sph );

inline void
link_vertex_to_generators( PDvertexOP vrt )
{
	//TR << "New Vertex has generators: ";
	for ( auto pvert : vrt->nonconst_generators() ) {
		pvert->vertices().push_back( vrt.get() );
	}
	return;
}

Vector vertex_xyz_from_generators( PDsphere const * a1, PDsphere const * a2, PDsphere const * a3, PDsphere const * a4 );
Vector vertex_xyz_from_generators( utility::vector1< PDsphere* > const & gen );

PDvertex * find_next_vertex_with_smallest_dist( PDvertex * srch_vrt, PDsphereOP & new_sph, Real & this_dist );

///////////////////////////////////////////////////////////////////
// Start Surface Area Calculation Functions ///////////////////////
///////////////////////////////////////////////////////////////////


void find_intersections(
	PDsphere const * patom,
	PDvertex const * vrt1,
	PDvertex const * vrt2,
	Vector const & start_pt,
	Vector const & dir,
	core::Real const max_extent,
	std::list< PDinterOP > & intersections
);

bool
find_intersections(
	PDsphereCOP patom,
	Vector const & start_pt,
	Vector const & dir,
	core::Real const max_extent
);

void
find_common_intersection_atoms( PDinterOP inter );

utility::vector1< utility::vector1< SAnode > >
get_cycles_from_intersections( std::list< PDinterOP > & intersections, PDsphere const * this_atom );

core::Real
get_sasa_from_cycles( utility::vector1< utility::vector1< SAnode > > & cycles, PDsphere* this_atom );

void
get_derivs_from_cycles( utility::vector1< utility::vector1< SAnode > > & cycles,
	PDsphereOP & this_atom, PDsphereOP & check_atom, Vector & f1, Vector & f2 );

void
get_derivs_from_cycle( utility::vector1< SAnode > & cycle,
	PDsphereOP & this_atom, PDsphereOP & check_atom, Vector & f1, Vector & f2 );

#ifdef NOTDEF
void
check_deriv_cycles(
		utility::vector1< utility::vector1< SAnode > > & cycles1,
		utility::vector1< utility::vector1< SAnode > > & cycles2,
		PDsphereOP & this_atom, PDsphereOP & check_atom );

void
check_deriv_cycle(
		utility::vector1< SAnode > & cycle1,
		utility::vector1< SAnode > & cycle2,
		PDsphereOP & this_atom, PDsphereOP & check_atom );
#endif

bool share_axis_atoms( PDinterCOP v1, PDsphere const * a1, PDsphere const * a2 );

core::Real get_area_from_cycle( PDsphere* this_atom, utility::vector1< SAnode > & cycle );


///////////////////////////////////////////////////////////////////
// End Surface Area Calculation Functions /////////////////////////
///////////////////////////////////////////////////////////////////


// Diagnostic output functions

// Prints out a list of points as dummy PDB lines for visualization
void print_points( std::list< PDinterOP > & inters );

// Prints out the partner vertices for the input vertex
void print_partners( PDvertexCOP vrt );

// Prints out the generator atoms for the input vertex
void print_generators( PDvertexCOP vrt );

// Prints out finite vertices as dummy PDB lines, and infinite vertices along
// the line that connects them to their single finite partner
void print_vertices( std::list< PDvertexOP > & fv, std::list< PDvertexOP > & iv );

// Prints out finite vertices, their partners and generators as dummy PDB lines
void print_vertices( std::list< PDvertexOP > & v );

}
}
}


#endif /*INCLUDED_core_scoring_power_diagram_PowerDiagram_HH*/
