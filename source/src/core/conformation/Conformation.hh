// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  conformation container
/// @file   core/conformation/Conformation.hh
/// @author Phil Bradley


#ifndef INCLUDED_core_conformation_Conformation_hh
#define INCLUDED_core_conformation_Conformation_hh

// Unit headers
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.hh>

// Package headers
#include <core/conformation/signals/ConnectionEvent.fwd.hh>
#include <core/conformation/signals/GeneralEvent.fwd.hh>
#include <core/conformation/signals/IdentityEvent.fwd.hh>
#include <core/conformation/signals/LengthEvent.fwd.hh>
#include <core/conformation/signals/XYZEvent.fwd.hh>

#ifdef WIN32
#include <core/conformation/signals/ConnectionEvent.hh>
#include <core/conformation/signals/GeneralEvent.hh>
#include <core/conformation/signals/IdentityEvent.hh>
#include <core/conformation/signals/LengthEvent.hh>
#include <core/conformation/signals/XYZEvent.hh>
#endif

// Project headers
#include <core/chemical/ResidueType.fwd.hh>
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/DOF_ID.hh>
#include <core/id/NamedAtomID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/AtomTree.fwd.hh>
#include <core/kinematics/DomainMap.fwd.hh>
#include <core/kinematics/FoldTree.fwd.hh>

// Utility headers
#include <utility/pointer/access_ptr.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/signals/BufferedSignalHub.hh>
#include <utility/signals/Link.hh>
#include <utility/signals/PausableSignalHub.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>


namespace core {
namespace conformation {

/// @brief A container of Residues and the kinematics to manage them
class Conformation : public utility::pointer::ReferenceCount {

public: // typedefs

	typedef kinematics::Jump Jump;
	typedef kinematics::FoldTree   FoldTree;
	typedef kinematics::FoldTreeOP FoldTreeOP;
	typedef kinematics::AtomTree   AtomTree;
	typedef kinematics::AtomTreeOP AtomTreeOP;
	typedef id::AtomID AtomID;

	//typedef id::AtomID_Map AtomID_Map;

	typedef id::AtomID_Mask AtomID_Mask;
	typedef id::DOF_ID DOF_ID;
	typedef id::TorsionID TorsionID;
	typedef kinematics::DomainMap DomainMap;

	typedef core::conformation::signals::ConnectionEvent ConnectionEvent;
	typedef core::conformation::signals::GeneralEvent GeneralEvent;
	typedef core::conformation::signals::IdentityEvent IdentityEvent;
	typedef core::conformation::signals::LengthEvent LengthEvent;
	typedef core::conformation::signals::XYZEvent XYZEvent;

	// Mirroring typedefs in AtomTree.hh to avoid it's #inclusion.
	// for fragment insertions
	typedef std::map< id::AtomID, Vector > FragXYZ;
	typedef std::map< id::StubID, kinematics::RT > FragRT;

public:

	// APL Removing accessor functions that voilate the data integrity guarantees of this class.
	// Conformation forbids non-const access to its residues.  Iterate from 1 to total_residue and
	// request a Residue const & instead of iterating from res_begin to res_end.
	/// @brief HIGHLY HIGHLY ILLEGAL ACCESS GRANTED TO CONFORMATION DATA.
	/// This function will be removed very very shortly.
	/// @brief Returns a random-access iterator that points at the first residue in the Conformation.
	/// ResidueOPs::iterator res_begin() { return residues_.begin(); }
	/// @brief Returns a random-access iterator that points just beyond the last residue in the Conformation.
	/// ResidueOPs::iterator res_end  () { return residues_.end  (); }

	/////////////////////////////////////////////////////////////////////////////

	/// @brief constructor
	/// if you are using PyRosetta, you should NOT BE HERE!
	Conformation();
		//	utility::pointer::ReferenceCount(),
		//residue_coordinates_need_updating_( false ),
		//residue_torsions_need_updating_( false ),
		//structure_moved_( true )
		//{}

	/// @brief default destructor
	virtual
	~Conformation();

	/// @brief copy constructor
	Conformation( Conformation const & src );

	/// @brief operator
	virtual
	Conformation &
	operator=( Conformation const & src );

	/// @brief clone the conformation
	virtual
	ConformationOP
	clone() const;

	///@details determine the type of the ConformationOP
	virtual
	bool
	same_type_as_me( Conformation const & other, bool recurse /* = true */ ) const;

	/// @brief do the names of all residues in this and src match?
	bool
	sequence_matches( Conformation const & other ) const;

	/// @brief Returns the number of residues in the Conformation
	Size
	size() const
	{
		return residues_.size();
	}

	/// @brief Returns true if this conformation does not have any residues
	bool
	empty() const
	{
		return residues_.empty();
	}

	/// @brief Returns the position number of the last residue in  <chain>
	Size
	chain_end( Size const chain ) const
	{
		if ( chain <= chain_endings_.size() ) return chain_endings_[chain];
		else return size();
	}

	/// @brief Returns the position number of the first residue in  <chain>
	Size
	chain_begin( Size const chain ) const
	{
		if ( chain == 1 ) return 1;
		else return chain_endings_[ chain-1 ]+1;
	}

	/// @brief Returns the number of chains
	Size
	num_chains() const
	{
		return chain_endings_.size() + 1; // last residue is not counted as chain ending
	}


	/// @brief Return true if this conformation contains any carbohydrate residues.
	bool
	contains_carbohydrate_residues() const
	{
		return contains_carbohydrate_residues_;
	}

	/// @brief Set whether this conformation contains any carbohydrate residues.
	void
	contains_carbohydrate_residues(bool const setting)
	{
		contains_carbohydrate_residues_ = setting;
	}


	/// @brief Returns the secondary structure the position  <seqpos>
	/// @return character representing secondary structure; returns 'L' if the
	/// requested sequence position is larger than the length in the
	/// secondary structure array
	char secstruct( Size const seqpos ) const
	{
		if ( Size(seqpos) > secstruct_.size() ) return 'L';
		return secstruct_[seqpos];
	}

	/// @brief Sets the secondary structure of the position  <seqpos>  to  <setting>
	/// @details Sets secondary structure character of a sequence position.
	/// Will resize the secondary structure array if the requested sequence
	/// position is larger than the length of the array.
	virtual void
	set_secstruct( Size const seqpos, char const setting )
	{
		if ( secstruct_.size() < Size(seqpos) ) secstruct_.resize( seqpos, 'L' );
		secstruct_[seqpos] = setting;
	}

	/// @brief Returns the conformation's FoldTree
	virtual FoldTree const &
	fold_tree() const
	{
		return *fold_tree_;
	}

	/// @brief Returns the conformation's AtomTree
	AtomTree const &
	atom_tree() const
	{
		return *atom_tree_;
	}

	/// @brief Sets the FoldTree to  <fold_tree_in>
	virtual void
	fold_tree( FoldTree const & fold_tree_in );


	/// @brief Returns the list of chain endings
	utility::vector1< Size > const &
	chain_endings() const;


	/// @brief Sets the list of chain endings
	/// @remarks All positions must be strictly less than the number of
	/// residues in the Conformation, otherwise the routine will fail fast.
	/// Note that the last residue position is not counted as a chain end.
	void
	chain_endings( utility::vector1< Size > const & endings );


	/// @brief Marks  <seqpos>  as the end of a new chain
	/// @remarks The last residue position is not counted as a chain ending.
	/// Also increases the chain ID number by 1 for all residues upstream from seqpos.
	void
	insert_chain_ending( Size const seqpos );

	/// @brief Deletes  <seqpos>  from the list of chain endings
	/// @remarks The last residue position is not counted as a chain ending.
	void
	delete_chain_ending( Size const seqpos );


	/// @brief Resets chain data so that the Conformation is marked as a single chain
	void
	reset_chain_endings();


	/// @brief Rederive the chains from the termini/polymer status
	void
	chains_from_termini();

	/////////////////////////////////////////////////////////////////////////////

	/// @brief Returns the AA enum for position  <seqpos>
	chemical::AA const &
	aa( Size const seqpos ) const {
		assert( seqpos >= 1 );
		assert( seqpos <= size() );
		return residues_[seqpos]->aa();
	}

	/// @brief access one of the residues;  this access is inlined, since otherwise it
	/// shows up in the profiler.  This will call non-inlined refold methods if necessary.
	///
	/// @details update coordinates and torsions for this and all other residues before
	/// allowing read access
	inline
	Residue const &
	residue( Size const seqpos ) const {
		runtime_assert( seqpos >= 1 );
		runtime_assert( seqpos <= size() );
		if ( residue_coordinates_need_updating_ ) update_residue_coordinates();
		if ( residue_torsions_need_updating_ )    update_residue_torsions();
		return *( residues_[ seqpos ] );
	}

	/// @brief access one of the residue's types -- avoids coord/torsion update
	chemical::ResidueType const &
	residue_type( Size const seqpos ) const
	{
		assert( seqpos >=1 );
		assert( seqpos <= size() );
		return residues_[seqpos]->type();
	}

	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////
	// insert/append/delete residues
	///////////////////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////////////////

	/// @brief Append a new residue by a jump.
	virtual
	void
	append_residue_by_jump(
		conformation::Residue const & new_rsd,
		Size const anchor_residue,
		std::string const& anchor_atom = "", // the atom in the anchor_residue
		std::string const& root_atom = "", // the atom in the new residue
		bool const start_new_chain = false
	);


	/// @brief Insert a new residue by jump.  If new_chain is "true", then seqpos must be the last 
	/// residue of one chain (i.e. residue(seqpos).chain() != residue(seqpos+1).chain() )
	void
	insert_residue_by_jump(
		Residue const & new_rsd_in,
		Size const seqpos, // desired seqpos of new_rsd
		Size anchor_pos, // in the current sequence numbering, ie before insertion of seqpos
		std::string const& anchor_atom = "",
		std::string const& root_atom = "",
		bool new_chain = false // insert this residue as a new chain, displacing all downstream chains
	);

	/// @brief Append a new residue by a bond.
	void
	append_residue_by_bond(
		conformation::Residue const & new_rsd,
		bool const build_ideal_geometry = false,
		int connection_index = 0,
		Size anchor_residue = 0,
		int anchor_connection_index = 0,
		bool const start_new_chain = false,
		bool const lookup_bond_length = false
	);

	/// @brief glues to seqpos and perhaps also seqpos+1
	void
	append_polymer_residue_after_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	/// @brief glues to seqpos and perhaps also seqpos+1, removes termini variants if necessary
	void
	safely_append_polymer_residue_after_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	/// @brief glues to seqpos and perhaps also seqpos-1
	void
	prepend_polymer_residue_before_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	/// @brief glues to seqpos and perhaps also seqpos-1, removes termini variants if necessary
	void
	safely_prepend_polymer_residue_before_seqpos(
		Residue const & new_rsd,
		Size const seqpos,
		bool const build_ideal_geometry
	);

	/// @brief replace residue
	virtual void
	replace_residue(
		Size const seqpos,
		Residue const & new_rsd,
		bool const orient_backbone
	);

		/// @brief funtion to replace a residue based on superposition on
	/// @brief the specified input atom pairs
	/// @brief NOTE: at the moment, only superposition on 3 atoms works
	virtual void
	replace_residue(
		Size const seqpos,
		Residue const & new_rsd,
		utility::vector1< std::pair< std::string, std::string > > const & atom_pairs
	);

	/// @brief Delete polymer residue at the given sequence position
	void
	delete_polymer_residue( Size const seqpos );

	/// @brief Slow method that relies on FoldTree::delete_seqpos, rebuilds atomtree, can handle jumps/root residue
	void
	delete_residue_slow( Size const seqpos );

	/// @brief Slow method that relies on FoldTree::delete_seqpos, rebuilds atomtree, can handle jumps/root residue
	void
	delete_residue_range_slow( Size const range_begin, Size const range_end );

	/// @brief Declare that a chemical bond exists between two residues
	void
	declare_chemical_bond(
												Size const seqpos1,
												std::string const & atom_name1,
												Size const seqpos2,
												std::string const & atom_name2
												);


  /// @brief  Insert one conformation into another. See FoldTree::insert_fold_tree_by_jump
	virtual
	void
	insert_conformation_by_jump(
			Conformation const & conf,             // the conformation to be inserted
			Size const insert_seqpos,              // rsd 1 in conf goes here
			Size const insert_jumppos,             // jump#1 in conf goes here, see insert_fold_tree_by_jump
			Size const anchor_pos,                 // in the current sequence numbering, ie before insertion of conf
			Size const anchor_jump_number = 0,     // the desired jump number of the anchoring jump, default=0
			std::string const & anchor_atom = "",  // "" means take default anchor atom
			std::string const & root_atom   = ""   // "" means take default root   atom
	);

	//////////////////////////////////////////////////////////////////////////
	/// @brief copy a stretch of coordinates/torsions from another Conformation
	void
	copy_segment(
		Size const size,
		Conformation const & src,
		Size const begin,
		Size const src_begin
	);

	/// @brief Rebuild the atoms ( like HN(seqpos), OC(seqpos+1) ) that are dependent on the polymer bond between seqpos and seqpos+1
	void
	rebuild_polymer_bond_dependent_atoms( Size const seqpos );

	/// @brief  Set the transform between two stubs -- only works if there's a jump between the two sets of stubatoms
	void
	set_stub_transform(
		id::StubID const & stub_id1,
		id::StubID const & stub_id2,
		kinematics::RT const & target_rt
	);
	//{
	//	set_dof_moved( atom_tree_.set_stub_transform( stub_id1, stub_id2, target_rt ) );
	//}

	/// @brief  get the transform between two stubs
	kinematics::RT
	get_stub_transform(
		id::StubID const & stub_id1,
		id::StubID const & stub_id2
	) const;
	//{
	//	return atom_tree_.get_stub_transform( stub_id1, stub_id2 );
	//}


	void
	set_jump_atom_stub_id( id::StubID const& id );
	//{
	//	atom_tree_.set_jump_atom_stub_id( id );
	//}

	kinematics::Stub
	stub_from_id( id::StubID const& id ) const;
	//{
	//	return atom_tree_.stub_from_id( id );
	//}

	void
	rebuild_residue_connection_dependent_atoms( Size const seqpos, Size const connid );


	/// @brief  Inefficient -- constructs copy of residues_
	ResidueCAPs
	const_residues() const;

	/////////////////////////////////////////////////////////////////////////////
	// access/modify dofs/xyz's

	/// @brief Returns the AtomTree degree of freedom (DOF)  <id>
	Real
	dof( DOF_ID const & id ) const;
	//{
	//	return atom_tree_.dof( id );
	//}

	/// @brief Sets the AtomTree degree of freedom (DOF)  <id>  to  <setting>
	virtual
	void
	set_dof( DOF_ID const & id, Real const setting );
	//{
	//	set_dof_moved( id );
	//	residue_torsions_need_updating_ = true; // might have been a torsion angle
	//	atom_tree_.set_dof( id, setting );
	//}

	/// @brief Returns the torsion angle  <id>
	Real
	torsion( TorsionID const & id ) const;

	/// @brief Sets the AtomTree DOF and the torsion in the corresponding Residue
	virtual
	void
	set_torsion( TorsionID const & id, Real const setting );

	///
	void
	insert_ideal_geometry_at_polymer_bond( Size const seqpos );

	void
	insert_ideal_geometry_at_residue_connection( Size const pos1, Size const connid1 );

	/// @brief Sets the torsion angle defined by  <atom[1-4]>  to  <setting>
	virtual
	void
	set_torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		Real const setting
	);

	/// @brief Sets the bond angle defined by  <atom[1-3]>  to  <setting>
	virtual
	void
	set_bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		Real const setting
	);

	/// @brief Sets the cond length between  <atom1>  and  <atom2>  to  <setting>
	virtual
	void
	set_bond_length(
		AtomID const & atom1,
		AtomID const & atom2,
		Real const setting
	);

	///
	void
	insert_fragment(
		id::StubID const & instub_id,
		FragRT const & outstub_transforms,
		FragXYZ const & frag_xyz
	);

	/// @brief Returns the torsion angle defined by  <atom[1-4]>
	Real
	torsion_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4
	) const;
	//{
	//	return atom_tree_.torsion_angle( atom1, atom2, atom3, atom4 );
	//}

	/// @brief Returns the bond angle defined by  <atom[1-3]>
	Real
	bond_angle(
		AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3
	) const;
	//{
	//	return atom_tree_.bond_angle( atom1, atom2, atom3 );
	//}

	/// @brief Returns the bond length between  <atom1>  and  <atom2>
	Real
	bond_length(
		AtomID const & atom1,
		AtomID const & atom2
	) const;
	//{
	//	return atom_tree_.bond_length( atom1, atom2 );
	//}


	/// @brief Returns the Jump with jump number  <jump_number>
	const Jump &
	jump( int const jump_number ) const;
	//{
	//	return atom_tree_.jump( jump_atom_id( jump_number ) );
	//}

	/// @brief Sets the jump  <jump_number>  to  <new_jump>
	virtual
	void
	set_jump(
		int const jump_number,
		Jump const & new_jump
	);
	//{
	//	assert( new_jump.ortho_check() );
	//	AtomID const id( jump_atom_id( jump_number ) );
	//	atom_tree_.set_jump( id, new_jump );
	//	set_dof_moved( id );
	//}

	/// @brief Sets a jump and forces immediate calculation of affected XYZ coords
	virtual
	void
	set_jump_now(
		int const jump_number,
		Jump const & new_jump
	);
	//{
	//	AtomID const id( jump_atom_id( jump_number ) );
	//	assert( new_jump.ortho_check() );
	//	atom_tree_.set_jump_now( id, new_jump );
	//	set_dof_moved( id );
	//}

	/// @brief access a jump
	const Jump &
	jump( AtomID const & id ) const;
	//{
	//	return atom_tree_.jump( id );
	//}

	/// @brief set a jump
	virtual
	void
	set_jump(
		AtomID const & id,
		Jump const & new_jump
	);
	//{
	//	assert( new_jump.ortho_check() );
	//	atom_tree_.set_jump( id, new_jump );
	//	set_dof_moved( id );
	//}

	/// @brief identify polymeric connections
	void
	set_polymeric_connection(
		Size res_id_lower,
		Size res_id_upper
	);

	/// @brief Update the polymer connection status between lower_seqpos and lower_seqpos+1 based on chainID's and termini
	void
	update_polymeric_connection( Size const lower_seqpos );

	void
	detect_bonds();

	void
	detect_pseudobonds();

	/// @brief Detect existing disulfides from the protein structure.
	/// @note Assumes full atom
	/// @details looks at SG-SG distance. If the SG-SG are about 2.02 A apart, calls
	/// it a disulfide bond.
	virtual
	void
	detect_disulfides();

	/// @brief Assigns disulfide bonds based on a pre-determined list
	/// @note works in centroid and full-atom modes
    void
    fix_disulfides(utility::vector1< std::pair<Size, Size> > disulf_bonds);

    // void find_disulfides();

	/// @brief The upstream and downstream Stubs are the coordinate frames between which this jump is transforming
	kinematics::Stub
	upstream_jump_stub( int const jump_number ) const;
	//{
	//	return atom_tree_.atom( jump_atom_id( jump_number ) ).get_input_stub();
	//}

	/// @brief  The upstream and downstream Stubs are the coordinate frames between which this jump is transforming
	kinematics::Stub
	downstream_jump_stub( int const jump_number ) const;
	//{
	//	return atom_tree_.atom( jump_atom_id( jump_number ) ).get_stub();
	//}

	/// @brief access xyz coordinates of an atom
	PointPosition const &
	xyz( AtomID const & id ) const;
	//{
	//	return atom_tree_.xyz( id );
	//}

	///
	virtual
	void
	set_xyz( AtomID const & id, PointPosition const & position );

	///
	virtual
	void
	batch_set_xyz( utility::vector1<AtomID> const & id, utility::vector1<PointPosition> const & position );

	///
	virtual
	void
	batch_get_xyz( utility::vector1<AtomID> const & id, utility::vector1<PointPosition> & position ) const;

	//
	void
	update_actcoords();

	void
	update_actcoord( Size resid );

	void update_orbital_coords( Residue & rsd) const;

	void
	update_orbital_coords( Size resid );

	/////////////////////////////////////////////////////////////////////////////
	// ID access and conversions

	///
	DOF_ID
	dof_id_from_torsion_id( TorsionID const & id ) const;

	///
	id::AtomID
	jump_atom_id( int const jump_number ) const;

	///@brief get four atoms which defined this torsion
	///@note  Returns TRUE to signal FAILURE
	bool
	get_torsion_angle_atom_ids(
		TorsionID const & tor_id,
		AtomID & id1,
		AtomID & id2,
		AtomID & id3,
		AtomID & id4
	) const;


	/////////////////////////////////////////////////////////////////////////////
	// for tracking changes to the structure

	/// @brief Generate a domain_map from the current dof/xyz moved data
	void
	update_domain_map( DomainMap & domain_map ) const;

	/// @brief has the structure moved since the last call to reset_move_data or reset_structure_moved
	bool
	structure_moved() const
	{
		return structure_moved_;
	}

	/// @brief reset the structure_moved_ bool
	void
	reset_structure_moved() const
	{
		structure_moved_ = false;
	}

	/// @brief forget all the structure modifications
	/**
		 called after domain map information is transferred to the Energies class
	**/
	void
	reset_move_data();
	//{
	//	structure_moved_ = false;
	//	dof_moved_.fill_with( false );
	//	xyz_moved_.fill_with( false );
	//}


	/////////////////////////////////////////////////////////////////////////////
	/// @brief clear data
	void
	clear();


	/////////////
	/// @brief debugging
	void
	debug_residue_torsions( bool verbose = false ) const;

	/// @brief Show residue connections for debugging purposes.
	void
	show_residue_connections() const;

	/// @brief Show residue connections for debugging purposes.
	void
	show_residue_connections(std::ostream &os) const;

	/// @brief  This returns the AtomID of the atom in the other residue to which the "connection_index"-th
	/// @brief  connection of residue seqpos is connected to.

	AtomID
	inter_residue_connection_partner(
		Size const seqpos,
		int const connection_index
	) const;


	/// @brief get all atoms bonded to another
	utility::vector1<core::id::AtomID>
	bonded_neighbor_all_res(
		core::id::AtomID atomid,
		bool virt = false
	) const;

	///
	void
	fill_missing_atoms(
		id::AtomID_Mask missing
	);

	bool
	atom_is_backbone_norefold( Size const pos, Size const atomno ) const;

	/// @brief returns a mask of residues to be used in scoring
	virtual utility::vector1<bool>
	get_residue_mask() const;

	/// @brief returns a residue-pair weight
	virtual Real
	get_residue_weight(core::Size, core::Size) const;

public: // observer management


	/// @brief attach ConnectionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( ConnectionEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection
	/// @remarks ConnectionEvent observers will only be notified upon a change
	///  in the state of the connection with the Conformation, e.g. if the
	///  Conformation is destroyed or if the connection is being transferred.
	///  SUGGESTION: Try to use the Link objects that are returned when attaching
	///  observers instead of attaching via this function and watching for
	///  ConnectionEvents, as it typically makes connection management easier.
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_connection_obs( MemFn fn, Ptr ptr ) const {
		return connection_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach ConnectionEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( ConnectionEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	/// @remarks ConnectionEvent observers will only be notified upon a change
	///  in the state of the connection with the Conformation, e.g. if the
	///  Conformation is destroyed or if the connection is being transferred.
	///  SUGGESTION: Try to use the Link objects that are returned when attaching
	///  observers instead of attaching via this function and watching for
	///  ConnectionEvents, as it typically makes connection management easier.
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_connection_obs( MemFn fn, Ptr ptr ) const {
		return connection_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach GeneralEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( GeneralEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_general_obs( MemFn fn, Ptr ptr ) const {
		return general_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach GeneralEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( GeneralEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_general_obs( MemFn fn, Ptr ptr ) const {
		return general_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach IdentityEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( IdentityEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_identity_obs( MemFn fn, Ptr ptr ) const {
		return identity_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach IdentityEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( IdentityEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_identity_obs( MemFn fn, Ptr ptr ) const {
		return identity_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach LengthEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( LengthEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_length_obs( MemFn fn, Ptr ptr ) const {
		return length_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach LengthEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( LengthEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_length_obs( MemFn fn, Ptr ptr ) const {
		return length_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief attach XYZEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( XYZEvent const & )
	/// @param ptr pointer to observer object
	/// @return Link that can be used to manage the connection
	template< typename MemFn, typename Ptr >
	inline
	utility::signals::Link
	attach_xyz_obs( MemFn fn, Ptr ptr ) const {
		return xyz_obs_hub_.connect( fn, ptr );
	}


	/// @brief detach XYZEvent observer function
	/// @param fn pointer to observer's unary member function with signature void( XYZEvent const & )
	/// @param ptr pointer to observer object
	/// @return true if disconnect successful, false if connection does not exist
	template< typename MemFn, typename Ptr >
	inline
	bool
	detach_xyz_obs( MemFn fn, Ptr ptr ) const {
		return xyz_obs_hub_.disconnect( fn, ptr );
	}


	/// @brief clear all observers
	/// @remarks ConnectionEvent::DISCONNECT will be sent to all observers
	void clear_observers();


	/// @brief fire a ConnectionEvent::TRANSFER to transfer observers from some
	///  source Conformation
	/// @param src Take observers from this source Conformation.
	/// @remarks Only observers that properly honor the TRANSFER event by
	///  re-attaching themselves to 'this' Conformation will be transferred.
	void receive_observers_from( Conformation const & src );


public: // additional observer behavior


	/// @brief wait for stdin after sending a GeneralEvent signal
	void
	debug_pause( bool const flag ) const;


	/// @brief waiting for stdin after sending a GeneralEvent signal?
	bool
	debug_pause() const;


public: // signal management

	///@brief convenience test for residue_type_set ( based on two middle residue -- to avoid hitting on ligands or pseudos )
	///@note this is not a good test --Doug
	bool is_residue_typeset( std::string tag ) const;

	///@brief convenience test for residue_type_set ( based on two middle residue -- to avoid hitting on ligands or pseudos )
	///@note this is not a good test --Doug
	bool is_fullatom() const;

	///@brief convenience test for residue_type_set ( based on two middle residue -- to avoid hitting on ligands or pseudos )
	///@note this is not a good test --Doug
	bool is_centroid() const;


	/// @brief block signals from being sent and buffer them to be
	///  sent after unblocking
	void
	buffer_signals();


	/// @brief block signals from being sent
	/// @warning for safety, ConnectionEvent are never blocked
	void
	block_signals();


	/// @brief allow signals to be sent
	void
	unblock_signals();


	/// @brief are signals being blocked and buffered?
	bool
	buffering_signals() const;


	/// @brief are signals being blocked?
	bool
	blocking_signals() const;


	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// private methods
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

private:

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void save(Archive & ar, const unsigned int version) const {
			int m = residues_.size();
			ar << m;
			Residue * tmp;
			for( int i = 1; i <= residues_.size(); i++ ) {
				tmp = residues_[i].get(); 
				ar << tmp;
			}
			// second passthrough is to turn off serialized flag inside ResidueType
			// This is how we store each new ResType only once, even though there may
			// be many residues using it
			for( int i = 1; i <= residues_.size(); i++ ) {
				residues_[i]->serialized(false);
			}
			ar << secstruct_;
			ar << fold_tree_;
			ar << chain_endings_;
			ar << residue_coordinates_need_updating_;
			ar << residue_torsions_need_updating_;
			ar << dof_moved_;
			ar << xyz_moved_;
			ar << structure_moved_;
	}

	template<class Archive>
	void load(Archive & ar, const unsigned int version) {
			int idx;
			ar >> idx;
			Residue * tmpp = NULL;
			for( int i = 1; i <= idx; i++ ) {
				ar >> tmpp;
				ResidueOP tmp( tmpp );
				residues_.push_back( tmp );
			}
			ar >> secstruct_;
			ar >> fold_tree_;
			ar >> chain_endings_;
			setup_atom_tree();
			ar >> residue_coordinates_need_updating_;
			ar >> residue_torsions_need_updating_;
			ar >> dof_moved_;
			ar >> xyz_moved_;
			ar >> structure_moved_;
	}
	BOOST_SERIALIZATION_SPLIT_MEMBER()
#endif

	/// @brief Returns a residue without triggering coordinate/torsion update
	/// @details Use with care. Useful inside torsion/coordinate setters where we want chemical info
	/// about a given residue but don't want to trigger the coordinate/torsion updates that go along
	/// with a call to residue(seqpos)
	Residue const &
	residue_( Size const seqpos ) const
	{
		assert( seqpos >=1 );
		assert( seqpos <= size() );
		return *residues_[ seqpos ];
	}

	/// @brief remap *_moved arrays, sequence numbering in the residues_ arrays, etc, after insertion or deletion of rsds
	void
	update_sequence_numbering(
		Size const new_size,
		utility::vector1< Size > const & old2new
	);

	/// @brief rebuild atoms in residue seqpos dependent on either the lower (-1) or upper(1) polymer residue
	void
	rebuild_polymer_bond_dependent_atoms( Size const seqpos, int const upper_lower );

	///
	void
	insert_polymer_residue(
		Residue const & new_rsd_in,
		Size const seqpos, // desired seqpos of new_rsd
		bool const join_lower,
		bool const join_upper
	);

	/// @brief  Now a private method
	/// public interface:  append_residue_by_bond or append_residue_by_jump
	///
	void
	append_residue(
		Residue const & new_rsd_in,
		bool const attach_by_jump,
		std::string const& root_atom,
		id::NamedAtomID anchor_id,
		bool const start_new_chain
	);


	////////////////////////////////////////////////////////////////////////////
	/// @brief wrap direct access to the Residues container
	void
	residues_replace(
		Size const seqpos,
		Residue const & new_rsd
	);

	///
	void
	residues_insert(
		Size const seqpos,
		Residue const & new_rsd,
		bool const use_lower_chain = false,
		bool const new_chain = false
	);

	///
	void
	residues_append( Residue const & new_rsd, bool const start_new_chain );

	///
	void
	residues_delete( Size const seqpos );
	//////////////////////////////////////////////////////////////////////////////



	/// @brief (re-)builds the AtomTree using the FoldTree and the Residues
	void
	setup_atom_tree();

	/// @brief access a torsion from the atom_tree
	Real
	atom_tree_torsion( TorsionID const & tor_id ) const;

	///@brief get four backbone atoms which define this backbone torsion
	///@note  Returns TRUE to signal FAILURE
	bool
	backbone_torsion_angle_atoms(
		TorsionID const & id,
		AtomID & id1,
		AtomID & id2,
		AtomID & id3,
		AtomID & id4
	) const;


	/////////////////////////////////////////////////////////////////////////////
	// setting the moved data

	/// @brief notify of xyz-change
	void
	set_xyz_moved( AtomID const & id )
	{
		structure_moved_ = true;
		xyz_moved_[ id ] = true;
		residue_torsions_need_updating_ = true;
	}

	/// @brief notify of mutiple-xyz-change
	void
	set_xyz_moved( utility::vector1<AtomID> const & ids )
	{
		structure_moved_ = true;
		for(core::Size i=1; i<=ids.size(); ++i)
			xyz_moved_[ ids[i] ] = true;
		residue_torsions_need_updating_ = true;
	}


	/// @brief notify of dof-change
	void
	set_dof_moved( AtomID const & id )
	{
		structure_moved_ = true;
		dof_moved_[ id ] = true;
		// Residues xyz not in sync with internal xyz
		residue_coordinates_need_updating_ = true;
	}

	/// @brief notify of dof-change
	void
	set_dof_moved( DOF_ID const & id )
	{
		structure_moved_ = true;
		dof_moved_[ id.atom_id() ] = true;
		// Residues xyz not in sync with internal xyz
		residue_coordinates_need_updating_ = true;
	}


	/// @brief  Will (if necessary) copy the xyz coordinates from the AtomTree to the Residues being managed
	/// Always safe to call. Nothing will happen unless coords_need_updating_ is true.
	void
	update_residue_coordinates() const;

	/// @brief called by above
	void
	update_residue_coordinates( Size const seqpos, bool const fire_signal = true ) const;

	///
	void
	rederive_chain_endings();

	///
	void
	rederive_chain_ids();

	/// @brief  Will (if necessary) copy the torsion angles (mainchain/chi) from the AtomTree to the Residues being managed
	/// Always safe to call. Nothing will happen unless torsions_need_updating_
	/// is true.
	void
	update_residue_torsions() const;

	/// called by above
	void
	update_residue_torsions( Size const seqpos, bool const fire_signal = true ) const;

	void
	add_pseudobond(
		Size lr,
		Size lr_connid,
		Size ur,
		Size ur_connid,
		Size nbonds
	);

	void
	in_place_copy(
		Conformation const & src
	);

	/// @brief The Conformation must transfer lazily accumulated "moved" data from the AtomTree before
	/// changing the number of residues.
	void
	pre_nresidue_change();

	private: // observer notifications


	/// @brief notify ConnectionEvent observers
	/// @remarks called upon a change in the state of connection between
	///  the Conformation and the observer (e.g. destruction of Conformation
	///  or transfer of connection)
	void
	notify_connection_obs( ConnectionEvent const & e ) const;


	/// @brief notify GeneralEvent observers
	/// @remarks should only be called when there are no other suitable event types
	///  since specific event notifications will automatically fire a GeneralEvent signal
	void
	notify_general_obs( GeneralEvent const & e ) const;


	/// @brief notify IdentityEvent observers
	/// @param e the event
	/// @param fire_general fire a GeneralEvent afterwards? default true
	void
	notify_identity_obs( IdentityEvent const & e, bool const fire_general = true ) const;


	/// @brief notify LengthEvent observers
	/// @param e the event
	/// @param fire_general fire a GeneralEvent afterwards? default true
	void
	notify_length_obs( LengthEvent const & e, bool const fire_general = true ) const;


	/// @brief notify XYZEvent observers
	/// @param e the event
	/// @param fire_general fire a GeneralEvent afterwards? default true
	void
	notify_xyz_obs( XYZEvent const & e, bool const fire_general = true ) const;

private:
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////
	// data
	/////////////////////////////////////////////////////////////////////////////
	/////////////////////////////////////////////////////////////////////////////

protected:
	/// @brief container of Residues
	ResidueOPs residues_;

private:
	/// ResidueCOPs const_residues_; // mirrors residues_ allowing const access -- this will be reinstated soon.

	/// @brief chain number for each position
	/**
		 conformation is in charge of making sure that the Residue chain
		 ID's and the chain_endings_ vector stay in sync.
	**/
	utility::vector1< Size > chain_endings_;

	/// @brief fold tree for the kinematics
	FoldTreeOP fold_tree_;
	/// @brief atom tree for the kinematics
	AtomTreeOP atom_tree_;

	// Does this conformation contain any carbohydrate residues at all?
	bool contains_carbohydrate_residues_;

	/// @brief do we need to update the coordinates in the Residues?
	mutable bool residue_coordinates_need_updating_;

	/// @brief do we need to update the torsion angles in the Residues?
	mutable bool residue_torsions_need_updating_;


	/// @brief book-keeping array for energy evaluations
	/**
		 store which DOF's have changed since the last call to reset_move_data
		 note that we are not currently differentiating dof's from the same atom
	**/
	AtomID_Mask dof_moved_;


	/// @brief book-keeping array for energy evaluations
	/**
		 store which xyz's have changed since the last call to reset_move_data
	**/
	AtomID_Mask xyz_moved_;


	/// @brief has the structure moved since the last call to reset_move_data?
	mutable bool structure_moved_;

	///
	utility::vector1< char > secstruct_;

	/// @brief ConnectionEvent observers
	/// @remarks Notification only occurs when there is a change in the state
	///  of the connection between observers and the Conformation object, e.g.
	///  destruction or transfer of the connection.
	mutable utility::signals::BufferedSignalHub< void, ConnectionEvent > connection_obs_hub_;

	/// @brief GeneralEvent observers
	/// @remarks GeneralEvent observers will be notified whenever any signal
	///  derived from GeneralEvent occurs.
	mutable utility::signals::PausableSignalHub< void, GeneralEvent > general_obs_hub_;

	/// @brief IdentityEvent observers
	mutable utility::signals::BufferedSignalHub< void, IdentityEvent > identity_obs_hub_;

	/// @brief LengthEvent observers
	mutable utility::signals::BufferedSignalHub< void, LengthEvent > length_obs_hub_;

	/// @brief LengthEvent observers
	mutable utility::signals::BufferedSignalHub< void, XYZEvent > xyz_obs_hub_;
};

std::ostream &operator<< (std::ostream &os, Conformation const &conf);

} // conformation
} // core

#endif
