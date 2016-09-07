// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/vardist_solaccess/VarSolDRotamerDots.hh
/// @brief  VarSolDRotamerDots classes header file
/// @author Andrew Leaver-Fay
/// @author Ron Jacak

#ifndef INCLUDED_devel_vardist_solaccess_VarSolDRotamerDots_HH
#define INCLUDED_devel_vardist_solaccess_VarSolDRotamerDots_HH

// Unit Headers
#include <devel/vardist_solaccess/VarSolDRotamerDots.fwd.hh>

// Project headers
#include <core/id/AtomID_Map.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pack/interaction_graph/RotamerDots.fwd.hh>
#include <core/types.hh>

#include <basic/MetricValue.fwd.hh>

//Utilitiy Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

//C++ Headers
#include <algorithm>
#include <cmath>
#include <vector>

#include <ObjexxFCL/FArray2D.fwd.hh>


#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/access.fwd.hpp>
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace devel {
namespace vardist_solaccess {

/// @brief Handles sphere-sphere overlap calculations
class VarSolDRotamerDots : public utility::pointer::ReferenceCount {

public:
	/*
	VarSolDRotamerDots( core::conformation::ResidueCOP rotamer, bool allatoms = false,
	core::Real probe_radius=1.4, core::Real wobble=0.0);
	*/

	VarSolDRotamerDots(
		core::conformation::ResidueCOP rotamer,
		VarSolDistSasaCalculatorCOP vsasa_calc
	);

	~VarSolDRotamerDots() override;
	VarSolDRotamerDots(const VarSolDRotamerDots& rhs);

	void copy( VarSolDRotamerDots const & rhs );
	VarSolDRotamerDots const & operator= ( VarSolDRotamerDots const & rhs );

	// do two rotamers have any sphere overlaps? doesn't actually determine which dots are covered if there is overlap.
	bool overlaps( VarSolDRotamerDots const & other ) const;

	core::conformation::ResidueCOP
	rotamer() const;

	//bool state_unassigned() const;
	core::Size get_num_atoms() const;
	core::Vector get_atom_coords_xyz( core::Size atom_index ) const;
	core::Real get_atom_collision_radius( core::Size atom_index ) const;
	core::Real get_atom_interaction_radius( core::Size atom_index ) const;

	core::Size nshells_for_atom( core::Size atom_index ) const;
	core::Real shell_radius_for_atom( core::Size atom_index, core::Size shell_index ) const;
	core::Size ndots() const;

	bool get_dot_covered( core::Size atom_index, core::Size shell_index, core::Size dot_index ) const;

	// add rotamer coverage counts produced by self overlap.
	void increment_self_overlap();

	void intersect_residues( VarSolDRotamerDots & other_res );

	void
	write_dot_kinemage( std::ostream & kinfile ) const;

	bool
	any_exposed_dots( core::Size atom ) const;

	// @brief The area of the molecular surface accessible to solvent.  Computes an "AND" for each
	// dot on atoms with a variable solvent radius.  MSAS radii taken from the "SASA_RADIUS_LEGACY" extra
	// property taken from the database.
	core::Real
	msas_for_atom( core::Size atom_index ) const;

	//core::Real
	//get_probe_radius() const;

	//core::Real
	//get_wobble() const;

private:

	bool
	get_atom_overlap_masks(
		VarSolDRotamerDots const & other,
		core::Size at_this,
		core::Size at_other,
		core::Real & distance,
		core::Size & closest_dot1,
		core::Size & closest_dot2
	) const;

	bool
	overlap_atoms(
		VarSolDRotamerDots const & other,
		core::Size at_this,
		core::Size at_other,
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & at_this_coverage,
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & at_other_coverage
	) const;

	void initialize_sasa_arrays();
	core::Real interaction_radii_squared( core::Size attype1, core::Size attype2 ) const;

	void
	write_dot(
		std::ostream & kinfile,
		core::Size atom,
		core::Size dot,
		core::Real radius
	) const;

	void write_dot(
		std::ostream & kinfile,
		core::Vector const & coord,
		std::string const & atname
	) const;

	void
	write_dotlist_header(
		std::ostream & kinfile,
		std::string const & master_name,
		std::string const & color
	) const;


	VarSolDistSasaCalculatorCAP owner_;
	core::conformation::ResidueCOP rotamer_;

	//const core::Real& probe_radius_;
	//const core::Real& wobble_;
	const utility::vector1< utility::vector1< core::Real > >& radii_;
	const utility::vector1< core::Real >& msas_radii_;
	const utility::vector1< core::Real >& coll_radii_;
	const utility::vector1< core::Real >& int_radii_;
	const utility::vector1< utility::vector1< core::Real > >& int_radii_sum_;
	const utility::vector1< utility::vector1< core::Real > >& int_radii_sum2_;
	const ObjexxFCL::FArray2D_int*   lg_angles_;
	const ObjexxFCL::FArray2D_ubyte* lg_masks_;
	const core::Size& num_bytes_;
	const core::Real& polar_expansion_radius_;

	core::Size num_atoms_;
	utility::vector1< utility::vector1< core::pack::interaction_graph::DotSphere > > atom_coverage_;
	utility::vector1< core::Vector > dot_coords_;

	//VarSolDistSasaCalculatorCOP vsasa_calc_;
#ifdef    SERIALIZATION
protected:
	friend class cereal::access;

public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > static void load_and_construct( Archive & arc, cereal::construct< VarSolDRotamerDots > & construct );
#endif // SERIALIZATION

};


class VarSolDistSasaCalculator :
	public core::pose::metrics::StructureDependentCalculator,
	public utility::pointer::enable_shared_from_this< VarSolDistSasaCalculator >
{
	friend class VarSolDRotamerDots;
public:
	VarSolDistSasaCalculator();

	inline VarSolDistSasaCalculatorCOP get_self_ptr() const      { return shared_from_this(); }
	inline VarSolDistSasaCalculatorOP  get_self_ptr()            { return shared_from_this(); }
	inline VarSolDistSasaCalculatorCAP get_self_weak_ptr() const { return VarSolDistSasaCalculatorCAP( shared_from_this() ); }
	inline VarSolDistSasaCalculatorAP  get_self_weak_ptr()       { return VarSolDistSasaCalculatorAP( shared_from_this() ); }



	core::pose::metrics::PoseMetricCalculatorOP clone() const override;

	void set_element_radii(std::string atype_name, core::Real coll_radius, core::Real int_radius, Size nshells);
	void set_atom_type_radii(std::string atype_name, core::Real coll_radius, core::Real int_radius, Size nshells);

	core::id::AtomID_Map< core::Real > calculate(const core::pose::Pose& pose);
protected:
	void lookup(const std::string& key, basic::MetricValueBase* valptr ) const override;
	std::string print(const std::string& key ) const override;
	void recompute(const core::pose::Pose & this_pose ) override;

private:
	core::Real probe_radius_;
	core::Real wobble_;
	core::Real total_sasa_;
	core::id::AtomID_Map< core::Real > atom_sasa_;
	utility::vector1< core::Real > residue_sasa_;
	utility::vector1< VarSolDRotamerDotsOP > rotamer_dots_vec_;

	void initialize_sasa_arrays();
	//static core::Real interaction_radii_squared( core::Size attype1, core::Size attype2 );
	core::Real interaction_radii_squared( core::Size attype1, core::Size attype2 ) const;

	utility::vector1< utility::vector1< core::Real > > radii_;
	utility::vector1< core::Real > msas_radii_;
	utility::vector1< core::Real > coll_radii_; // collision radii
	utility::vector1< core::Real > int_radii_;  // interaction radii -- larger than coll radii for N and O
	utility::vector1< utility::vector1< core::Real > > int_radii_sum_; // atom-type pair interaction radii sums
	utility::vector1< utility::vector1< core::Real > > int_radii_sum2_; // atom-type pair interaction radii sums, squared
	core::Size const num_bytes_;
	const ObjexxFCL::FArray2D_ubyte* lg_masks_;
	const ObjexxFCL::FArray2D_int*   lg_angles_;

	core::Real polar_expansion_radius_;

	bool up_to_date;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // vardist_solaccess
} // core


#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( devel_vardist_solaccess_VarSolDRotamerDots )
#endif // SERIALIZATION


#endif // INCLUDED_devel_vardist_sollaccess_VarSolDRotamerDots_HH
