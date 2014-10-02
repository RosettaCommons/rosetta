// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/interaction_graph/RotamerDots.hh
/// @brief  RotamerDots classes header files - ported from rosetta++
/// @author Andrew Leaver-Fay
/// @author Ron Jacak

#ifndef INCLUDED_core_pack_interaction_graph_RotamerDots_hh
#define INCLUDED_core_pack_interaction_graph_RotamerDots_hh

// Unit Headers
// AUTO-REMOVED #include <core/chemical/AtomType.hh>
#include <core/conformation/Atom.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/interaction_graph/RotamerDots.fwd.hh>
#include <core/types.hh>

// Objexx Headers
// AUTO-REMOVED #include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/ubyte.hh>

#ifdef PYROSETTA
#include <ObjexxFCL/FArray2D.hh>
#endif

// Numeric headers
// AUTO-REMOVED #include <numeric/xyzVector.hh>

//Utilitiy Headers
#include <utility/SingletonBase.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>
#include <utility/thread/ReadWriteMutex.hh>

//C++ Headers
#include <vector>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray2D.fwd.hh>

namespace core {
namespace pack {
namespace interaction_graph {

///
/// @begin DotSphere
///
/// @brief
/// Represents the sphere of dots on the vdW surface of an atom, for use in the LeGrand and Merz method of calculating SASA.
///
/// @detailed
/// For every atom in a protein, the vdW surface is patterned with dots. Each dot has to keep track of how many other
/// residues are "covering" this dot. So, that's 1 character for each dot.  Each character is the count of the number of
/// residues overlapping with this dot.  An assumption we're making here is that a single atom (or really, dot) will never
/// be covered by more than 255 residues.
///
/// In this implementation of the LeGrand and Merz algorithm, we're going to be using 162 dots per atom. Turns out that
/// you can distribute 162 dots evenly on the surface of a sphere.
///
/// This class is extremely simple. The RotamerDots class below does all the work of tying a particular residues atoms to
/// DotSpheres. As a matter of fact, DotSphere doesn't even know what atom it's representing. It just has the one C-style
/// array for the coverage count and that's it.
///
class DotSphere {

public:
	DotSphere();
	~DotSphere();
	DotSphere( DotSphere const & rhs );
	DotSphere const & operator= ( DotSphere const & rhs );
	bool operator!= ( DotSphere const & rhs );

	void zero();

	inline core::Size get_total_dots() const { return NUM_DOTS_TOTAL; }
	void increment_count( utility::vector1< ObjexxFCL::ubyte > const & );

	void count_num_covered() const;
	core::Size get_num_uncovered() const;
	core::Size get_num_covered() const;

	DotSphere const & operator -= ( DotSphere const & rhs );
	DotSphere const & operator += ( DotSphere const & rhs );

	void print( std::ostream & os) const;

	bool get_dot_covered( core::Size dot_index ) const;

	void write_to_compact_array( utility::vector1< ObjexxFCL::ubyte > & compact ) const;
	void invert_to_compact_array( utility::vector1< ObjexxFCL::ubyte > & inv_compact ) const;

	static core::Size const NUM_DOTS_TOTAL = 162;

private:

	static core::Size const NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS = 21;
	static core::Size const NUM_COUNTS_TO_ALLOCATE = NUM_BYTES_IN_DOT_SPHERE_OVERLAP_ARRAYS * 8;

	unsigned char dots_coverage_count_[ NUM_COUNTS_TO_ALLOCATE ]; // sized to 168
	mutable core::Size num_covered_;
	mutable bool num_covered_current_;

};

std::ostream & operator<< ( std::ostream & os, DotSphere const & ds );



///
/// @begin RotamerDots
///
/// @brief
/// Handles sphere-sphere overlap calculations for the HPatchInteractionGraph.
///
/// @detailed
/// One big change from the r++ version of this class is that this class now includes all of the information that was
/// previously stored in the RotamerCoords class. Since I'm not storing atoms in trie ordering (perhaps I'll add this
/// later), there is no need to have a separate class for the rotamer coordinates.
///
/// RotamerDots hold DotSphere objects for the atoms of a given Residue (or really, of the current state on some interaction
/// graph Node). Default use of the class will result in RotamerDots objects keeping DotSpheres for every atom of a residue.
/// For the hpatch interaction graph, though, we only care about the SASA of the heavy atoms.  No need to include the
/// hydrogens when looking for hydrophobic patches. By not keeping track of the hydrogens, we save a huge amount of time
/// on computing updates to the SASA score because hydrogen atoms generally make up half of a protein. So I'm adding a
/// boolean flag to the non-default constructor which toggles whether we're tracking SASA of all atoms, or just the heavy
/// atoms.
///
/// Two other big changes being made to this class are that 1) the class will now keep track of two kinds of SASA and
/// 2) it will no longer keep score or score_is_current variables.  The two kinds of SASA the class will keep track of
/// are the standard SASA, and a SASA with polar atom radii extended. The expanded polar SASA will only be kept if a
/// boolean flag is set at construct time.  If not, it will just calculate standard SASA and that's it.  The second
/// change is that this class is now only keeping track of SASA. RotamerDots objects will not be responsible for calculating
/// a score.
///
class RotamerDots : public utility::pointer::ReferenceCount {

public:
	RotamerDots();
	RotamerDots( conformation::ResidueCOP rotamer, bool exclude_hydrogen_atoms = false, bool use_expanded_polar_atom_radii = false );
	virtual ~RotamerDots();
	RotamerDots( RotamerDots const & rhs );

	void copy(RotamerDots const & rhs);
	RotamerDots const & operator= ( RotamerDots const & rhs );
	bool operator!= ( RotamerDots const & rhs );
	// zeros out the dot coverage counts on all atoms. only used to reinit IG's for multiple packing runs.
	void zero();

	// do two rotamers have any sphere overlaps? doesn't actually determine which dots are covered if there is overlap.
	bool overlaps( RotamerDots const & other ) const;

	core::conformation::ResidueCOP
	rotamer() const;

	bool state_unassigned() const;
	core::Size get_num_atoms() const;
	numeric::xyzVector< Real > get_atom_coords_xyz( Size atom_index ) const;
	core::Real get_atom_radius( Size atom_index ) const;
	core::Real radius_for_attype( Size const attype_index );
	core::Real max_atom_radius();
	utility::vector1< Real > const * get_radii() const;

	// add rotamer coverage counts produced by self overlap.
	void increment_self_overlap();

	// increments the atom counts on this only, and saves the overlap into the passed in rotamerdotscache
	void increment_this_and_cache( RotamerDots const & other, RotamerDotsCache & this_overlap_on_other,
		utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache );

	// rotamerdots frequently have to calculate how they overlap with another rotamerdots object. this method figures it all out.
	void get_res_res_overlap( RotamerDots const & other_res,
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res1_covered_by_res2,
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & res2_covered_by_res1,
		utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache ) const;

	// given two atoms, figures out which dots are covered on both atoms and sets the appropriate bits in the passed in vectors
	bool get_atom_atom_coverage( Vector const & at1_xyz, Real at1_base_radius, Vector const & at2_xyz, Real at2_base_radius,
		utility::vector1< ObjexxFCL::ubyte > & at1_sphere_covered, utility::vector1< ObjexxFCL::ubyte > & at2_sphere_covered, Real dist_sq ) const;

	void increment_from_cached( RotamerDotsCache const & cached_dot_overlap );
	void decrement_from_cached( RotamerDotsCache const & cached_dot_overlap );

	// add rotamer coverage counts for dots on both this and other.
	void increment_both( RotamerDots & other );

	void increment_both_and_cache( RotamerDots & other_rotamer, RotamerDotsCache & others_dots_covered_by_this,
		RotamerDotsCache & this_dots_covered_by_other, utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache );

	core::Real get_sasa() const;
	core::Real get_atom_sasa( Size atom_index ) const;

	core::Size get_num_uncovered( core::Size atom ) const;
	core::Size get_num_covered_total() const;

	//void write_dot_kinemage( std::ofstream & kinfile );

	void print( std::ostream & os ) const;
	std::string name3() const;
	core::Size seqpos() const;

	// used by the unit tests only
	inline utility::vector1< DotSphere > const & get_atom_counts() { return atom_counts_; }

	void invert_to_boolmasks( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots ) const;

	void invert_to_boolmasks(
		utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & inv_dots,
		utility::vector1< Size > const & ats_to_update
	) const;

	static core::Vector dot_coord( Size index );

private:
	static void initialize_sasa_arrays();

	void get_overlap_cache(
		RotamerDots const & other,
		RotamerDotsCache & others_dots_covered_by_this,
		RotamerDotsCache & this_dots_covered_by_other,
		utility::vector1< utility::vector1< bool > > & atom_atom_overlaps_cache
	) const;

	//void write_dotlist_header( std::ofstream & kinfile, std::string master_name, std::string color );
	//void write_dot( std::ofstream & kinfile, core::Size atom, core::Size dot, core::Real radius );
	//void write_dot( std::ofstream & kinfile, numeric::xyzVector< core::Real > const & coord, std::string atname );

	//numeric::xyzVector< core::Real > const & get_dot_coord( core::Size dot_id );
	//void initialize_dot_sphere_coordinates_from_file();

public:
	static core::Size const num_bytes_;
	static core::Real probe_radius_;

private:
	core::conformation::ResidueCOP rotamer_;
	Size num_atoms_;
	utility::vector1< DotSphere > atom_counts_;

	static bool sasa_arrays_initialized_;
	mutable core::Real sasa_;
	mutable bool sasa_is_current_;
	mutable utility::vector1< core::Real > atom_sasa_;

	utility::vector1< Real > const * radii_;

	// APL FIX THIS: This belongs in its own singleton
	static utility::vector1< core::Vector > dot_coords_;

	static void initialize_dot_coords( utility::vector1< core::Vector > & dot_coords );

public:	/// TEMP!

	static ObjexxFCL::FArray2D_int const *   lg_angles_;
	static ObjexxFCL::FArray2D_ubyte const * lg_masks_;

};

std::ostream & operator<< ( std::ostream & os, RotamerDots const & rd );


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
class RotamerDotsRadiusData : public utility::SingletonBase< RotamerDotsRadiusData >
{
public:
	friend class utility::SingletonBase< RotamerDotsRadiusData >;

public:
	/// @brief return a pointer to the standard Rosetta SASA radii
	utility::vector1< Real > const * get_ROSETTA_SASA_radii() const;

	/// @brief return a pointer to the SASA radii used by NACCESS
	utility::vector1< Real > const * get_NACCESS_SASA_radii() const;

	/// @brief return a pointer to the SASA radii used by NACCESS, with polar atom radii expanded
	utility::vector1< Real > const * get_NACCESS_SASA_radii_with_expanded_polars() const;

private:
	/// @brief private constructor
	RotamerDotsRadiusData();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static RotamerDotsRadiusData * create_singleton_instance();

private:

	core::Real polar_expansion_radius_;

	utility::vector1< Real > ROSETTA_SASA_radii_;
	utility::vector1< Real > NACCESS_SASA_radii_;
	utility::vector1< Real > NACCESS_SASA_radii_with_expanded_polars_;

};


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
class RotamerDotsCache {

public:
	friend class RotamerDots;

	RotamerDotsCache();
	RotamerDotsCache( core::Size num_atoms );
	RotamerDotsCache( RotamerDotsCache const & rhs );
	RotamerDotsCache const & operator= ( RotamerDotsCache const & rhs );
	~RotamerDotsCache();

	void resize( core::Size num_atoms );
	void zero();

	void increment_count( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & covered );
	void write_to_compact_array( utility::vector1< utility::vector1< ObjexxFCL::ubyte > > & compact ) const;

	void print_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const;
	void print( std::ostream & os ) const;

	// for unit tests only
	inline utility::vector1< DotSphere > const & get_atom_counts() const { return atom_counts_; }

private:
	utility::vector1< DotSphere > atom_counts_;

};


///
/// @begin InvRotamerDots
///
/// @brief
/// Used to determine whether the overlap between two atoms is buried or exposed.
///
class InvRotamerDots : public utility::pointer::ReferenceCount {

public:
	InvRotamerDots();
	InvRotamerDots( InvRotamerDots const & src );
	virtual ~InvRotamerDots();

	InvRotamerDots const & operator = ( InvRotamerDots const & rhs );

	void setup_from_rotamer_dots( RotamerDots const & rdots );
	void setup_from_rotamer_dots(
		RotamerDots const & rdots,
		utility::vector1< Size > const & ats_to_update
	);

	core::conformation::ResidueCOP
	rotamer() const;

	/// @brief Is the intersection between two atoms on this inv-rotamer-dots object exposed?
	bool atom_overlap_is_exposed( Size at1, Size at2 ) const;

	/// @brief Is the intersection between one atom on this inv-rotamer-dots object,
	/// and one atom on another inv-rotamer-dots object exposed?
	bool atom_overlap_is_exposed( Size at_this, InvRotamerDots const & other, Size at_other ) const;

	bool dot_exposed( Size atomid, Size dot_index ) const;

	void write_exposed_dots_to_kinemage( std::ostream & ostr, bool group = false ) const;

	void write_circle_intersection_mask_to_kinemage(
		std::ostream & ostr,
		Size const atom_this,
		InvRotamerDots const & invdots_other,
		Size const atom_other,
		bool group = false
	) const;

private:

	bool
	overlap_exposed(
		core::conformation::Atom const & at1,
		utility::vector1< ObjexxFCL::ubyte > const & at1exposed_dots,
		core::conformation::Atom const & at2,
		utility::vector1< ObjexxFCL::ubyte > const & at2exposed_dots
	) const;

	void print_dot_bit_string( utility::vector1< ObjexxFCL::ubyte > & values ) const;

private:
	core::conformation::ResidueCOP rotamer_;
	utility::vector1< utility::vector1< ObjexxFCL::ubyte > > inv_dots_;
	static Real const max_dist_from_dot_to_intersection;
	utility::vector1< Real > const * radii_;

};


} // interaction_graph
} // pack
} // core


#endif // INCLUDED_core_pack_interaction_graph_RotamerDots_HH
