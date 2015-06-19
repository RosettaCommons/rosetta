// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/MoveMap.hh
/// @brief  Declarations for the MoveMap class
/// @author Phil Bradley
/// @author Christopher Miles (cmiles@uw.edu)

#ifndef INCLUDED_core_kinematics_MoveMap_hh
#define INCLUDED_core_kinematics_MoveMap_hh

// Unit header
#include <core/kinematics/MoveMap.fwd.hh>

// Package headers
#include <core/kinematics/types.hh>

// Project headers
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/JumpID.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/DOF_ID_Mask.fwd.hh>

// Utility headers
#include <utility/py/PyAssert.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>
#include <vector>


namespace core {
namespace kinematics {

/// @brief A class specifying DOFs to be flexible or fixed
///
/// @details Currently there are two groups of data, one is a residue-based Torsion
/// definition, such as BB, CHI, NU, and JUMP; the other is an atom-based DOF
/// definition, such as bond length D, bond angle THETA, and torsion angle PHI,
/// which are used in the AtomTree. MoveMap does not automatically handle
/// conversion from one group to the other, i.e., setting PHI false for
/// DOF_type does not affect setting for BB and CHI torsion though they are
/// PHIs in atom-tree.
///
/// Within each group, there are multiple levels of control
/// (from general/high to specific/lower):
/// @li Torsion-based: TorsionType(BB, CHI, NU, BRANCH, JUMP) -> MoveMapTorsionID
/// (BB, CHI of one residue) -> TorsionID ( BB torsion 2 or CHI torsion 3 of
/// one residue)
/// @li DOF-base: DOF_type( D, THETA, PHI ) -> DOF_ID (D, THETA, PHI of one atom)
///
/// Settings for each level are stored in a map structure and they are only
/// added to the map when setting methods are invoked. As a result, MoveMap
/// does not behave like a "Boolean vector", which always contains setting for
/// each residue or atom in a conformation. Setting for a higher level will
/// override setting for lower levels (remove it from map); Similarly, when
/// querying a lower level finds no setting, it will check setting for its
/// higher level. For example, setting TorsionType BB to be true will remove
/// any data of BB setting for a residue or a specific BB torsion (such as
/// backbone psi) in a residue. And querying the setting for BB torsion 2 of
/// residue 4 will first check if there is any specific setting, if not, it will
/// check if there is a setting for all BB torsions for residue 4, if not
/// again, it will use the setting for BB torsions for all residues.
///
/// Example:
///     movemap = MoveMap()
/// See also:
///     Pose
///     MinMover
///     ShearMover
///     SmallMover
class MoveMap: public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~MoveMap();
	// ids
	typedef id::AtomID AtomID;
	typedef id::AtomID_Mask AtomID_Mask;
	typedef id::DOF_Type DOF_Type;
	typedef id::DOF_ID DOF_ID;
	typedef id::DOF_ID_Mask DOF_ID_Mask;
	typedef id::TorsionType TorsionType;
	typedef id::TorsionID TorsionID;

	// TODO: should probably make this a class/struct
	/// @brief  Our own specific torsion_id
	/// @details TorsionType can be BB, CHI, NU, or JUMP.  Therefore, it doesn't distinguish specific torsions,
	/// e.g., phi/psi/omega in BB torsion.  Useful when setting residue k backbone fixed.
	typedef std::pair< Size, TorsionType > MoveMapTorsionID;

	/// @brief flexible or fixed for this TorsionType (BB, CHI, NU, BRANCH, JUMP), for all residues
	typedef std::map< TorsionType, bool > TorsionTypeMap;

	/// @brief flexible or fixed for this TorsionType (BB, CHI, NU, BRANCH, JUMP), for one residue,
	/// e.g., no distinction between phi/psi/omega
	typedef std::map< MoveMapTorsionID, bool > MoveMapTorsionID_Map;

	/// @brief flexible or fixed for a single torsion, e.g., psi (BB torsion 2) of residue 10
	typedef std::map< TorsionID, bool > TorsionID_Map;

	/// @brief flexible or fixed for this DOF_Type (PHI, THETA, D, RB1-6), for all atoms
	typedef std::map< DOF_Type, bool > DOF_TypeMap;

	/// @brief flexible or fixed for a single DOF, eg, D of atom 5 in residue 10
	typedef std::map< DOF_ID, bool > DOF_ID_Map;

	/// @brief flexible or fixed jumps (fold-tree independent definition via residue pairs )
	typedef std::map< id::JumpID, bool > JumpID_Map;

	/// @brief default constructor
	MoveMap(){}

	MoveMapOP clone() const {
		return MoveMapOP( new MoveMap( *this ) );
	}

	/// @brief clear -- sets all to FALSE
	void
	clear();

	/// @brief Set whether or not BB TorsionType is moveable
	///
	/// example:
	///     movemap.set_bb(True)
	/// See also:
	///     MoveMap
	///     MoveMap.set_chi
	///     MoveMap.set_nu
	///     MoveMap.set_branches
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover
	inline
	void
	set_bb( bool const setting )
	{
		set( id::BB, setting );
	}

	/// @brief Sets whether or not the BB torsions of residue  <seqpos>  are movable
	///
	/// example:
	///     movemap.set_bb(49,True)
	/// See also:
	///     MoveMap
	///     MoveMap.set_chi
	///     MoveMap.set_nu
	///     MoveMap.set_branches
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover
	inline
	void
	set_bb( Size const seqpos, bool const setting )
	{
		PyAssert( (seqpos>0), "MoveMap::set_bb( Size const seqpos , bool const setting ): input variable seqpos has a meaningless value");
		set( MoveMapTorsionID( seqpos, BB ), setting );
	}

	/// @brief Sets BB torsions movable based on input array
	inline
	void
	set_bb( utility::vector1< bool > allow_bb )
	{
		for( Size ii = 1; ii <= allow_bb.size(); ++ii )
			set( MoveMapTorsionID( ii, id::BB ), allow_bb[ii] );
	}

	/// @brief Sets the BB torsions between residues  <begin>  and  <end>
	/// as movable, all other residues are non-movable
	///
	/// example:
	///     movemap.set_bb_true_range(40,60)
	/// See also:
	///     MoveMap
	///     MoveMap.set_bb
	///     MoveMap.set_chi
	///     MoveMap.set_nu
	///     MoveMap.set_branches
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover
	inline
	void
	set_bb_true_range( Size const begin, Size const end  )
	{
		PyAssert( (begin>0), "MoveMap::set_bb_true_range( Size const begin , Size const end): input variable begin has a meaningless value");
		PyAssert( (end>0), "MoveMap::set_bb_true_range( Size const begin , Size const end ): input variable end has a meaningless value");
		PyAssert( (begin <= end), "MoveMap::set_bb_true_range( Size const begin, Size const end): input variable begin < input variable end");
		set_bb(false);
		for( Size ir=begin; ir<=end; ++ir ) set_bb(ir, true);
	}

	/// @brief Prevents backbone torsion modifications to the intervals specified
	/// by <ranges>. Each element of the vector is a pair, which specifies the
	/// begin and end indices. Counting begins with 1.
	void set_ranges_unmodifiable(const std::vector<std::pair<Size, Size> >& ranges);

	/// @brief Sets whether or not CHI TorsionType is movable
	///
	/// example:
	///     movemap.set_chi(True)
	/// See also:
	///     MoveMap
	///     MoveMap.set_bb
	///     MoveMap.set_nu
	///     MoveMap.set_branches
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover
	inline
	void
	set_chi( bool const setting )
	{
		set( id::CHI, setting );
	}

	/// @brief Sets whether or not the CHI torsions of residue  <seqpos>  are movable
	///
	/// example:
	///     movemap.set_chi(49,True)
	/// See also:
	///     MoveMap
	///     MoveMap.set_bb
	///     MoveMap.set_nu
	///     MoveMap.set_branches
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover	inline
	void
	set_chi( Size const seqpos, bool const setting )
	{
		PyAssert( (seqpos>0), "MoveMap::set_chi( Size const seqpos , bool const setting ): input variable seqpos has a meaningless value");
		set( MoveMapTorsionID( seqpos, id::CHI ), setting );
	}

	/// @brief set CHI torsions movable based on input array
	inline
	void
	set_chi( utility::vector1< bool > allow_chi )
	{
		for( Size ii = 1; ii <= allow_chi.size(); ++ii )
			set( MoveMapTorsionID( ii, id::CHI ), allow_chi[ii] );
	}

	/// @brief Sets the chi torsions between residues <begin> and <end>
	/// as movable and all other residues are non-movable.
	///
	/// Example:
	///     movemap.set_chi_true_range(40, 60)
	/// See also:
	///     MoveMap
	///     MoveMap.set_bb
	///     MoveMap.set_bb_true_range
	///     MoveMap.set_chi
	///     MoveMap.set_nu
	///     MoveMap.set_nu_true_range
	///     MoveMap.set_branches
	///     MoveMap.set_branches_true_range
	///     Pose
	///     MinMover
	///     ShearMover
	///     SmallMover
	inline
	void set_chi_true_range(core::uint const begin, core::uint const end)
	{
		PyAssert((begin > 0), "MoveMap::set_chi_true_range(core::uint const begin, core::uint const end): "
				"Input variable <begin> has a meaningless value.");
		PyAssert((end > 0), "MoveMap::set_chi_true_range(core::uint const begin, core::uint const end): "
				"Input variable <end> has a meaningless value.");
		PyAssert((begin <= end), "MoveMap::set_chi_true_range(core::uint const begin, core::uint const end): "
				"Input variable <begin> must be <= input variable <end>.");
		set_chi(false);
		for(Size res = begin; res <= end; ++res) {
			set_chi(res, true);
		}
	}

	/// @brief Set whether or not NU TorsionTypes are movable.
	///
	/// @details Example: movemap.set_nu(True)\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_branches\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	inline void
	set_nu( bool const setting )
	{
		set( id::NU, setting );
	}

	/// @brief Set whether or not the NU torsions of residue <seqpos> are movable.
	/// @details Example: movemap.set_nu(49, True)\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_branches\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	void
	set_nu( core::uint const seqpos, bool const setting )
	{
		PyAssert( ( seqpos > 0 ), "MoveMap::set_nu( core::uint const seqpos, bool const setting ): "
				"Input variable <seqpos> has a meaningless value." );
		set( MoveMapTorsionID( seqpos, id::NU ), setting );
	}

	/// @brief Set which NU torsions are movable based on input array.
	inline void
	set_nus( utility::vector1< bool > const settings )
	{
		Size const n_settings( settings.size() );
		for ( uint i( 1 ); i <= n_settings; ++i ) {
			set( MoveMapTorsionID( i, id::NU ), settings[ i ] );
		}
	}

	/// @brief Set the NU torsions between residues <begin> and <end>
	/// as movable and all other residues are non-movable.
	/// @details Example: movemap.set_nu_true_range(40, 60)\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_bb_true_range\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_chi_true_range\n
	///     MoveMap.set_nu\n
	///     MoveMap.set_branches\n
	///     MoveMap.set_branches_true_range\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	inline void
	set_nu_true_range( core::uint const begin, core::uint const end )
	{
		PyAssert( ( begin > 0 ), "MoveMap::set_nu_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <begin> has a meaningless value." );
		PyAssert( ( end > 0 ), "MoveMap::set_nu_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <end> has a meaningless value." );
		PyAssert( ( begin <= end ), "MoveMap::set_nu_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <begin> must be <= input variable <end>." );
		set_nu( false );
		for ( uint resnum( begin ); resnum <= end; ++resnum ) {
			set_nu( resnum, true );
		}
	}


	/// @brief Set whether or not BRANCH TorsionTypes are movable.
	///
	/// @details Example: movemap.set_branches( True )\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_nu\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	inline
	void
	set_branches( bool const setting )
	{
		set( id::BRANCH, setting );
	}

	/// @brief Set whether or not the BRANCH torsions of residue <seqpos> are movable.
	/// @details Example: movemap.set_branches(49, True)\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_nu\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	void
	set_branches( core::uint const seqpos, bool const setting )
	{
		PyAssert( ( seqpos > 0 ), "MoveMap::set_branches( core::uint const seqpos, bool const setting ): "
				"Input variable <seqpos> has a meaningless value." );
		set( MoveMapTorsionID( seqpos, id::BRANCH ), setting );
	}

	/// @brief Set which BRANCH torsions are movable based on input array.
	inline void
	set_branches( utility::vector1< bool > const settings )
	{
		Size const n_settings( settings.size() );
		for ( uint i( 1 ); i <= n_settings; ++i ) {
			set( MoveMapTorsionID( i, id::BRANCH ), settings[ i ] );
		}
	}

	/// @brief Set the BRANCH torsions between residues <begin> and <end>
	/// as movable and all other residues are non-movable.
	/// @details Example: movemap.set_branches_true_range(40, 60)\n
	/// See also:\n
	///     MoveMap\n
	///     MoveMap.set_bb\n
	///     MoveMap.set_bb_true_range\n
	///     MoveMap.set_chi\n
	///     MoveMap.set_chi_true_range\n
	///     MoveMap.set_nu\n
	///     MoveMap.set_nu_true_range\n
	///     MoveMap.set_branches\n
	///     Pose\n
	///     MinMover\n
	///     ShearMover\n
	///     SmallMover
	inline void
	set_branches_true_range( core::uint const begin, core::uint const end )
	{
		PyAssert( ( begin > 0 ), "MoveMap::set_branches_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <begin> has a meaningless value." );
		PyAssert( ( end > 0 ), "MoveMap::set_branches_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <end> has a meaningless value." );
		PyAssert( ( begin <= end ), "MoveMap::set_branches_true_range( core::uint const begin, core::uint const end ): "
				"Input variable <begin> must be <= input variable <end>." );
		set_branches( false );
		for ( uint resnum( begin ); resnum <= end; ++resnum ) {
			set_branches( resnum, true );
		}
	}


	/// @brief Sets whether or not JUMP TorsionType is moveable
	///
	/// example:
	///     movemap.set_jump(True)
	inline
	void
	set_jump( bool const setting )
	{
		set( id::JUMP, setting );
	}

	/// @brief Sets the movability of JUMP  <jump_number>  to  <setting>
	///
	/// example:
	///     movemap.set_jump(1,True)
	inline
	void
	set_jump( int const jump_number, bool const setting )
	{
		PyAssert( (jump_number>0), "MoveMap::set_jump( int const jump_number , bool const setting ): input variable jump_number has a meaningless value");
		set( MoveMapTorsionID( jump_number, id::JUMP ), setting );
	}

	/// @brief set JUMP moveable or not for one specific residue pair
	/// this mechanism  does not mix with the "jump_nr" mechanism...
	/// i.e., if you set jump_nr 3 = movable ---> movemap does not know that this refers to say residue 23-89
	// but it does mix with the global switch set_jump( true/false ).
	inline
	void
	set_jump( Size const pos1, Size const pos2, bool const setting )
	{
		set_jump( id::JumpID( pos1, pos2 ), setting );
	}

	void
	set_jump( id::JumpID const & jump, bool const setting );

	/// @brief set a specific TorsionType movable: currently BB, CHI, NU, BRANCH, or JUMP
	void
	set( TorsionType const & t, bool const setting );

	/// @brief set TorsionType flexible or fixed for one residue, e.g., BB torsions for residue 10
	void
	set( MoveMapTorsionID const & id, bool const setting );

	/// @brief set an individual Torsion movable for now, e.g., "BB torsion 2 of residue 4"
	void
	set( TorsionID const & id, bool const setting );

	/// @brief set atom tree DOF, e.g., D, PHI, THETA
	void
	set( DOF_Type const & t, bool const setting );

	/// @brief set for an individual DoF, e.g., "PHI of Atom 3 in Residue 5"
	void
	set( DOF_ID const & id, bool const setting );

public: // accessors
	/// @brief Returns if BB torsions are movable or not for residue  <seqpos>
	///
	/// example:
	///     movemap.get_bb(49)
	inline
	bool
	get_bb( Size const seqpos ) const
	{
		PyAssert( (seqpos>0), "MoveMap::get_bb( Size const seqpos ): input variable seqpos has a meaningless value");
		return get( MoveMapTorsionID( seqpos, id::BB ) );
	}

	/// @brief Returns if SC torsions are movable or not for residue  <seqpos>
	///
	/// example:
	///     movemap.get_chi(49)
	inline
	bool
	get_chi( int const seqpos ) const
	{
		PyAssert( (seqpos>0), "MoveMap::get_chi( int const seqpos ): input variable seqpos has a meaningless value");
		return get( MoveMapTorsionID( seqpos, id::CHI ) );
	}

	/// @brief Return if NU torsions are movable or not for residue <seqpos>.
	/// @details: Example: movemap.get_nu(49)
	inline
	bool
	get_nu( core::uint const seqpos ) const
	{
		PyAssert( ( seqpos > 0 ), "MoveMap::get_nu(core::uint const seqpos): "
				"Input variable <seqpos> has a meaningless value." );
		return get( MoveMapTorsionID( seqpos, id::NU ) );
	}

	/// @brief Return if BRANCH torsions are movable or not for residue <seqpos>.
	/// @details: Example: movemap.get_branches(49)
	inline
	bool
	get_branches( core::uint const seqpos ) const
	{
		PyAssert( ( seqpos > 0 ), "MoveMap::get_branches(core::uint const seqpos): "
				"Input variable <seqpos> has a meaningless value." );
		return get( MoveMapTorsionID( seqpos, id::BRANCH ) );
	}

	/// @brief Returns if JUMP  <jump_number>  is movable or not
	///
	/// example:
	///     movemap.get_jump(1)
	bool
	get_jump( int const jump_number ) const
	{
		PyAssert( (jump_number>0), "MoveMap::get_jump( int const jump_number ): input variable jump_number has a meaningless value");
		return get( MoveMapTorsionID( jump_number, id::JUMP ) );
	}

	inline
	bool
	get_jump( Size const pos1, Size const pos2 ) const
	{
		return get_jump( id::JumpID( pos1, pos2 ) );
	}

	bool
	get_jump( id::JumpID const & jump ) const;

	/// @brief get setting for a specific TorsionType, such as "BB"
	bool
	get( TorsionType const & t ) const;

	/// @brief get TorsionType flexible or fixed for one residue, eg BB torsions for residue 10
	bool
	get( MoveMapTorsionID const & id ) const;

	bool
	get( TorsionID const & id ) const;

	/// @brief get the default for this type of DOF, eg "PHI"
	bool
	get( DOF_Type const & t ) const;

	/// @brief get the setting for an individual dof, eg, PHI of Atom 3 in Residue 5
	bool
	get( DOF_ID const & id ) const;

	/// @brief find the explicit setting for the given TorsionType
	/// @return iterator pointing to the TorsionType-bool pair, otherwise
	///  torsion_type_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	TorsionTypeMap::const_iterator
	find( TorsionType const & t ) const;

	/// @brief find the explicit setting for the given MoveMapTorsionID
	/// @return iterator pointing to the MoveMapTorsionID-bool pair, otherwise
	///  movemap_torsion_id_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	MoveMapTorsionID_Map::const_iterator
	find( MoveMapTorsionID const & id ) const;

	/// @brief find the explicit setting for the given TorsionID
	/// @return iterator pointing to the TorsionID-bool pair, otherwise torsion_id_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	TorsionID_Map::const_iterator
	find( TorsionID const & id ) const;

	/// @brief find the explicit setting for the given JumpID
	/// @return iterator pointing to the JumpID-bool pair, otherwise jump_id_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	JumpID_Map::const_iterator
	find( id::JumpID const & jump ) const;

	/// @brief find the explicit setting for the given DOF_Type
	/// @return iterator pointing to the DOF_Type-bool pair, otherwise dof_type_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	DOF_TypeMap::const_iterator
	find( DOF_Type const & t ) const;

	/// @brief find the explicit setting for the given DOF_ID
	/// @return iterator pointing to the DOF_ID-bool pair, otherwise dof_id_end()
	/// @warning Do not use this for general lookup, as it does not take
	///  into account the stringency levels.  Only use this when you need
	///  to check if a setting explicitly exists.
	DOF_ID_Map::const_iterator
	find( DOF_ID const & id ) const;

 public: // iterators
	/// @brief return an iterator pointing just past the last element of the TorsionTypeMap
	inline
	TorsionTypeMap::const_iterator
	torsion_type_begin() const { return torsion_type_map_.begin(); }

	/// @brief return an iterator pointing at the first element of the TorsionTypeMap
	inline
	TorsionTypeMap::const_iterator
	torsion_type_end() const { return torsion_type_map_.end(); }

	/// @brief return an iterator pointing at the first element of the MoveMapTorsionID_Map
	inline
	MoveMapTorsionID_Map::const_iterator
	movemap_torsion_id_begin() const { return move_map_torsion_id_map_.begin(); }

	/// @brief return an iterator pointing just past the last element of the MoveMapTorsionID_Map
	inline
	MoveMapTorsionID_Map::const_iterator
	movemap_torsion_id_end() const { return move_map_torsion_id_map_.end(); }

	/// @brief return an iterator pointing at the first element of the TorsionID_Map
	inline
	TorsionID_Map::const_iterator
	torsion_id_begin() const { return torsion_id_map_.begin(); }

	/// @brief return an iterator pointing just past the last element of the TorsionID_Map
	inline
	TorsionID_Map::const_iterator
	torsion_id_end() const { return torsion_id_map_.end(); }

	/// @brief return an iterator pointing just past the last element of the DOF_TypeMap
	inline
	DOF_TypeMap::const_iterator
	dof_type_begin() const { return dof_type_map_.begin(); }

	/// @brief return an iterator pointing at the first element of the DOF_TypeMap
	inline
	DOF_TypeMap::const_iterator
	dof_type_end() const { return dof_type_map_.end(); }

	/// @brief return an iterator pointing at the first element of the DOF_ID_Map
	inline
	DOF_ID_Map::const_iterator
	dof_id_begin() const { return dof_id_map_.begin(); }

	/// @brief return an iterator pointing just past the last element of the DOF_ID_Map
	inline
	DOF_ID_Map::const_iterator
	dof_id_end() const { return dof_id_map_.end(); }

	/// @brief return an iterator pointing at the first element of the JumpID_Map
	inline
	JumpID_Map::const_iterator
	jump_id_begin() const { return jump_id_map_.begin(); }

	/// @brief return an iterator pointing just past the last element of the JumpID_Map
	inline
	JumpID_Map::const_iterator
	jump_id_end() const { return jump_id_map_.end(); }

public: // movemap-movemap functionality

	/// @brief import settings from another MoveMap
	/// @return The total number of settings imported.
	/// @remarks This function calls set() for each setting that exists in the
	///  'rval' MoveMap in order from lowest to highest stringency.
	inline
	Size import( MoveMap const & rval ) {
		return import( rval, true, true ); // bool: import_true_settings, import_false_settings
	}

	/// @brief import only False settings from another MoveMap
	/// @return The total number of settings imported.
	/// @remarks This function calls set() for each setting that is marked as
	///  False in the 'rval' MoveMap in order from lowest to highest
	///  stringency.
	inline
	Size import_false( MoveMap const & rval ) {
		return import( rval, false, true ); // bool: import_true_settings, import_false_settings
	}

	/// @brief import only True settings from another MoveMap
	/// @return The total number of settings imported.
	/// @remarks This function calls set() for each setting that is marked as
	///  True in the 'rval' MoveMap in order from lowest to highest
	///  stringency.
	inline
	Size import_true( MoveMap const & rval ) {
		return import( rval, true, false ); // bool: import_true_settings, import_false_settings
	}

public: // i/o
	/// @brief Read MoveMap from file.
	void init_from_file( std::string const & filename );

public: // status
	/// @brief Give the TorsionType bool values up to a given residue number.
	void
	show( std::ostream & out, Size i) const;

	/// @brief Give the TorsionType bool values up to a given residue number.
	/// wrapper for PyRosetta
	void
	show( Size i) const { show(std::cout, i);};

	void
	show( std::ostream & out) const;

	void
	show() const { show(std::cout);};

private: // movemap-movemap functionality

	/// @brief import settings from another MoveMap
	/// @param[in] rval The MoveMap to import settings from.
	/// @param[in] import_true_settings Import True settings?
	/// @param[in] import_false_settings Import False settings?
	/// @return The total number of settings imported.
	/// @remarks This function calls set() for each setting that exists in the
	///  'rval' MoveMap in order from lowest to highest stringency.
	Size import(
		MoveMap const & rval,
		bool const import_true_settings,
		bool const import_false_settings
	);

 private:
	// implementation of the data
	// WARNING: if you add something here also add it to ::clear() and ::import()

	/// @brief flexible or fixed for this TorsionType (BB, CHI, NU, BRANCH, JUMP), for all residues
	TorsionTypeMap torsion_type_map_;

	/// @brief flexible or fixed for this TorsionType (BB, CHI, NU, BRANCH, JUMP), for one residue,
	/// e.g., no distinction between phi/psi/omega
	MoveMapTorsionID_Map move_map_torsion_id_map_;

	/// @brief flexible or fixed for a single torsion, eg, psi (BB torsion 2) of residue 10
	TorsionID_Map torsion_id_map_;

	/// @brief flexible or fixed for this DOF_Type (PHI, THETA, D, RB1-6), for all atoms
	DOF_TypeMap dof_type_map_;

	/// @brief flexible or fixed for a single DOF, eg, D of atom 5 in residue 10
	DOF_ID_Map dof_id_map_;

	JumpID_Map jump_id_map_;
};  // MoveMap

inline
std::ostream &
operator <<( std::ostream & out, MoveMap const & mm )
{
	mm.show( out );
	return out;
}

}  // namespace kinematics
}  // namespace core

#endif  // INCLUDED_core_kinematics_DOF_ID_HH
