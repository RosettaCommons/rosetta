// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
// AUTO-REMOVED #include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/metrics/PoseMetricCalculatorBase.hh>
#include <core/pack/interaction_graph/RotamerDots.fwd.hh>
#include <core/types.hh>

#include <basic/MetricValue.fwd.hh>

// Objexx Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
// AUTO-REMOVED #include <ObjexxFCL/ubyte.hh>

// Numeric headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>

//Utilitiy Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

//C++ Headers
#include <vector>

#include <ObjexxFCL/FArray2D.fwd.hh>



namespace devel {
namespace vardist_solaccess {

/// @brief Handles sphere-sphere overlap calculations
class VarSolDRotamerDots : public utility::pointer::ReferenceCount {

public:
	VarSolDRotamerDots();
	VarSolDRotamerDots( core::conformation::ResidueCOP rotamer, bool allatoms = false );
	virtual ~VarSolDRotamerDots();
	VarSolDRotamerDots( VarSolDRotamerDots const & rhs );

	void copy( VarSolDRotamerDots const & rhs );
	VarSolDRotamerDots const & operator= ( VarSolDRotamerDots const & rhs );
	//bool operator!= ( VarSolDRotamerDots const & rhs );

	//void zero();

	// do two rotamers have any sphere overlaps? doesn't actually determine which dots are covered if there is overlap.
	bool overlaps( VarSolDRotamerDots const & other ) const;

	core::conformation::ResidueCOP
	rotamer() const;

	//bool state_unassigned() const;
	core::Size get_num_atoms() const;
	core::Vector get_atom_coords_xyz( core::Size atom_index ) const;
	core::Real get_atom_collision_radius( core::Size atom_index ) const;
	core::Real get_atom_interaction_radius( core::Size atom_index ) const;
	//core::Real radius_for_attype( Size const attype_index );
	//core::Real max_atom_radius();
	//utility::vector1< Real >* get_radii() const;

	// add rotamer coverage counts produced by self overlap.
	void increment_self_overlap();

	// increments the atom counts on this only, and saves the overlap into the passed in rotamerdotscache
	//void increment_this_and_cache( RotamerDots const & other, RotamerDotsCache & this_overlap_on_other,
	//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache );

	void intersect_residues( VarSolDRotamerDots & other_res );

	// rotamerdots frequently have to calculate how they overlap with another rotamerdots object. this method figures it all out.
	//void get_res_res_overlap( RotamerDots const & other_res,
	//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res1_covered_by_res2,
	//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res2_covered_by_res1,
	//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache ) const;

	// given two atoms, figures out which dots are covered on both atoms and sets the appropriate bits in the passed in vectors
	//bool get_atom_atom_coverage( Vector const & at1_xyz, Real at1_base_radius, Vector const & at2_xyz, Real at2_base_radius,
	//	utility::vector1< ObjexxFCL::ubyte > & at1_sphere_covered, utility::vector1< ObjexxFCL::ubyte > & at2_sphere_covered, Real dist_sq ) const;

	//void increment_from_cached( RotamerDotsCache const & cached_dot_overlap );
	//void decrement_from_cached( RotamerDotsCache const & cached_dot_overlap );

	// add rotamer coverage counts for dots on both this and other.
	//void increment_both( RotamerDots & other );

	//void increment_both_and_cache( RotamerDots & other_rotamer, RotamerDotsCache & others_dots_covered_by_this,
	//	RotamerDotsCache & this_dots_covered_by_other, utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache );

	//core::Real get_sasa() const;
	//core::Real get_atom_sasa( Size atom_index ) const;

	//core::Size get_num_uncovered( core::Size atom ) const;
	//core::Size get_num_covered_total() const;

	//void write_dot_kinemage( std::ofstream & kinfile );

	//void print( std::ostream & os ) const;
	//std::string name3() const;
	//core::Size seqpos() const;

	// used by the unit tests only
	//inline utility::vector1< DotSphere > const & get_atom_counts() { return atom_counts_; }

	//void invert_to_boolmasks( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots ) const;

	//void invert_to_boolmasks(
	//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots,
	//	utility::vector1< Size > const & ats_to_update
	//) const;

	void
	write_dot_kinemage( std::ostream & kinfile ) const;

	bool
	any_exposed_dots( core::Size atom ) const;

	// @brief The area of the molecular surface accessible to solvent.  Computes an "AND" for each
	// dot on atoms with a variable solvent radius.  MSAS radii taken from the "SASA_RADIUS_LEGACY" extra
	// property taken from the database.
	core::Real
	msas_for_atom( core::Size atom_index ) const;

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

	static void initialize_sasa_arrays();
	static core::Real interaction_radii_squared( core::Size attype1, core::Size attype2 );

	//void get_overlap_cache(
	//	RotamerDots const & other,
	//	RotamerDotsCache & others_dots_covered_by_this,
	//	RotamerDotsCache & this_dots_covered_by_other,
	//	utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
	//) const;

	//void write_dotlist_header( std::ofstream & kinfile, std::string master_name, std::string color );
	//void write_dot( std::ofstream & kinfile, core::Size atom, core::Size dot, core::Real radius );
	//void write_dot( std::ofstream & kinfile, numeric::xyzVector< core::Real > const & coord, std::string atname );

	//numeric::xyzVector< core::Real > const & get_dot_coord( core::Size dot_id );
	//void initialize_dot_sphere_coordinates_from_file();

public:
	static core::Size const num_bytes_;
	static core::Real probe_radius_;
	static core::Real polar_expansion_radius_;

private:

	void
	write_dotlist_header(
		std::ostream & kinfile,
		std::string const & master_name,
		std::string const & color
	) const;

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

private:
	core::conformation::ResidueCOP rotamer_;
	core::Size num_atoms_;
	utility::vector1< utility::vector1< core::pack::interaction_graph::DotSphere > > atom_coverage_;

	static bool sasa_arrays_initialized_;
	static utility::vector1< utility::vector1< core::Real > > radii_;
	// radii for calculating the molecular surface area accessible to solvent; taken from the SASA_RADIUS
	// parameters specified in the database
	static utility::vector1< core::Real > msas_radii_;


	static utility::vector1< core::Real > coll_radii_; // collision radii
	static utility::vector1< core::Real > int_radii_;  // interaction radii -- larger than coll radii for N and O
	static utility::vector1< utility::vector1< core::Real > > int_radii_sum_; // atom-type pair interaction radii sums
	static utility::vector1< utility::vector1< core::Real > > int_radii_sum2_; // atom-type pair interaction radii sums, squared

	static utility::vector1< core::Vector > dot_coords_;


public:	/// TEMP!
	static ObjexxFCL::FArray2D_int const *   lg_angles_;
	static ObjexxFCL::FArray2D_ubyte const * lg_masks_;

};


class VarSolDistSasaCalculator : public core::pose::metrics::StructureDependentCalculator {
public:

	VarSolDistSasaCalculator();

	core::pose::metrics::PoseMetricCalculatorOP clone() const;

protected:
	virtual void lookup( std::string const & key, basic::MetricValueBase * valptr ) const;
	virtual std::string print( std::string const & key ) const;
	virtual void recompute( core::pose::Pose const & this_pose );

private:

	core::Real total_sasa_;
	core::id::AtomID_Map< core::Real > atom_sasa_;
	utility::vector1< core::Real > residue_sasa_;
	utility::vector1< VarSolDRotamerDotsOP > rotamer_dots_;
};

//std::ostream & operator<< ( std::ostream & os, RotamerDots const & rd );


///
/// @begin RotamerDotsRadiusData
///
/// @brief
/// A singleton class which reads in database SASA radii files and provides accessors for those values to the RotamerDots class.
///
/// @detailed
/// The RotamerDots class keeps track of the SASA of some rotamer using DotSphere objects and lots of get overlap calls
/// on atoms. The SASA of a given atom (or residue) depends on what radii are used for the atoms. Different program use
/// different sets of radii, and for the hpatch score, the polar atom radii are expanded. Previously RotamerDots objects
/// would have a static vector member variable that represented the particular set of radii being used.  It would also have
/// a second static vector member variable that represented the expanded polar version of the radii.  So, no matter how
/// many instances of RotamerDots objects we created, we only maintained two vectors for the radii.
/// Now, we want to make the RotamerDots class only use one set of radii for the hpatch score: the expanded polar radii.
/// The super easy solution would be to remove all the logic that keeps track of SASA when using the standard radii and
/// make the RotamerDots class only calculate that kind of SASA.  But, there will probably be more uses in the future
/// for a RotamerDots class that can calculate SASA with a standard set of radii than with the expanded polar ones. So, that's
/// where this class comes in. It will read in database files, depending on which set of radii are requested and provide
/// access to that data in some way.
///
/// Let's say we have to create RotamerDots objects for two different sets of SASA radii. When those objects are constructed
/// a pointer to the right set of radii will have to be set. Definitely don't want to copy/store the radii to each instance.
/// If I make it a pointer to this class, then it would only work for one set of radii.  So at RotamerDots construct time
/// we have to ask this class to give us a pointer to the set of radii that are requested. Then only one instance of this
/// class exists, and all RotamerDots objects have pointers set to the radii they care about. What would this change in
/// the RotamerDots API? The get_radius() function would dereference the pointer and index into the vector to get the
/// right radius. This structure is similar to how the ChemicalManager works. But, in that class, the functions return
/// AtomTypeSet objects.
///
//class RotamerDotsRadiusData {
//
//public:
//	static RotamerDotsRadiusData * get_instance();
//
//	/// @brief return a pointer to the standard Rosetta SASA radii
//	utility::vector1< Real >* get_ROSETTA_SASA_radii();
//
//	/// @brief return a pointer to the SASA radii used by NACCESS
//	utility::vector1< Real >* get_NACCESS_SASA_radii();
//
//	/// @brief return a pointer to the SASA radii used by NACCESS, with polar atom radii expanded
//	utility::vector1< Real >* get_NACCESS_SASA_radii_with_expanded_polars( Real polar_expansion_radius = 1.0 );
//
//private:
//	/// @brief private constructor
//	RotamerDotsRadiusData();
//
//	/// @brief static data member holding pointer to the singleton class itself
//	static RotamerDotsRadiusData * instance_;
//
//	utility::vector1< Real > ROSETTA_SASA_radii_;
//	utility::vector1< Real > NACCESS_SASA_radii_;
//	utility::vector1< Real > NACCESS_SASA_radii_with_expanded_polars;
//
//};


///
/// @begin RotamerDotsCache
///
/// @brief
/// A lightweight version of the RotamerDots class. Used to cache overlap between interaction graph Nodes and BGNodes.
///
/// @detailed
/// During packing, when a first class node has to respond to another Node changing state, it's faster to decrement the
/// coverage the previous state produced and increment the coverage the new state produces than to completely recalculate
/// how the new state overlaps the node. But instead of holding that coverage information in a RotamerDots object (which
/// does alot of other things), hold it in this cache class instead.
///
//class RotamerDotsCache {
//
//public:
//	friend class RotamerDots;
//
//	RotamerDotsCache();
//	RotamerDotsCache( core::Size num_atoms );
//	RotamerDotsCache( RotamerDotsCache const & rhs );
//	RotamerDotsCache const & operator= ( RotamerDotsCache const & rhs );
//	~RotamerDotsCache();
//
//	void resize( core::Size num_atoms );
//	void zero();
//
//	void increment_count( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & covered );
//	void write_to_compact_array( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & compact ) const;
//
//	void print_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const;
//	void print( std::ostream & os ) const;
//
//	// for unit tests only
//	inline utility::vector1< DotSphere > const & get_atom_counts() const { return atom_counts_; }
//
//private:
//	utility::vector1< DotSphere > atom_counts_;
//
//};


///
/// @begin InvRotamerDots
///
/// @brief
/// Used to determine whether the overlap between two atoms is buried or exposed.
///
//class InvRotamerDots : public utility::pointer::ReferenceCount {
//
//public:
//	InvRotamerDots();
//	InvRotamerDots( InvRotamerDots const & src );
//	~InvRotamerDots();
//
//	InvRotamerDots const & operator = ( InvRotamerDots const & rhs );
//
//	void setup_from_rotamer_dots( RotamerDots const & rdots );
//	void setup_from_rotamer_dots(
//		RotamerDots const & rdots,
//		utility::vector1< Size > const & ats_to_update
//	);
//
//	core::conformation::ResidueCOP
//	rotamer() const;
//
//	/// @brief Is the intersection between two atoms on this inv-rotamer-dots object exposed?
//	bool atom_overlap_is_exposed( Size at1, Size at2 ) const;
//
//	/// @brief Is the intersection between one atom on this inv-rotamer-dots object,
//	/// and one atom on another inv-rotamer-dots object exposed?
//	bool atom_overlap_is_exposed( Size at_this, InvRotamerDots const & other, Size at_other ) const;
//
//	bool dot_exposed( Size atomid, Size dot_index ) const;
//
//	void write_exposed_dots_to_kinemage( std::ostream & ostr, bool group = false ) const;
//
//	void write_circle_intersection_mask_to_kinemage(
//		std::ostream & ostr,
//		Size const atom_this,
//		InvRotamerDots const & invdots_other,
//		Size const atom_other,
//		bool group = false
//	) const;
//
//private:
//
//	bool
//	overlap_exposed(
//		core::conformation::Atom const & at1,
//		utility::vector1< ObjexxFCL::ubyte > const & at1exposed_dots,
//		core::conformation::Atom const & at2,
//		utility::vector1< ObjexxFCL::ubyte > const & at2exposed_dots
//	) const;
//
//	void print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const;
//
//private:
//	core::conformation::ResidueCOP rotamer_;
//	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > inv_dots_;
//	static Real const max_dist_from_dot_to_intersection;
//	utility::vector1< Real >* radii_;
//
//};


} // vardist_solaccess
} // core


#endif // INCLUDED_devel_vardist_sollaccess_VarSolDRotamerDots_HH
