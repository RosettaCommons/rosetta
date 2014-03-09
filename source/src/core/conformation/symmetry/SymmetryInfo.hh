// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief
/// @author Ingemar Andre, Phil Bradley

#ifndef INCLUDED_core_conformation_symmetry_SymmetryInfo_hh
#define INCLUDED_core_conformation_symmetry_SymmetryInfo_hh

#include <core/conformation/symmetry/SymmetryInfo.fwd.hh>

// Unit Headers
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>

//core
#include <core/types.hh>
#include <core/id/DOF_ID.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/AtomID.fwd.hh>
#include <core/conformation/Conformation.fwd.hh>

// Utility Headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <string>

#include <utility/vector1_bool.hh>


namespace core {
namespace conformation {
namespace symmetry {


//Symm_info

class SymmetryInfo : public utility::pointer::ReferenceCount {
public:

	typedef utility::vector1< Size > Clones; // NOTE vector1 *not* std::vector
	typedef utility::vector1< std::pair<Size,Real> > WtedClones; // NOTE vector1 *not* std::vector

	typedef id::DOF_ID DOF_ID;
	typedef id::TorsionID TorsionID;
	typedef id::AtomID AtomID;

	/// convenience: these could go somewhere else
	typedef utility::vector1< DOF_ID > DOF_IDs;
	typedef utility::vector1< TorsionID > TorsionIDs;
	typedef utility::vector1< AtomID > AtomIDs;

public:
	SymmetryInfo();
	virtual ~SymmetryInfo();
	//SymmetryInfo(SymmetryInfo const &);
	SymmetryInfo( SymmData const & symmdata, Size const nres_subunit, Size const njump_subunit );

	SymmetryInfo(
		Size const nres_monomer,
		Size const njump_monomer,
		Size const N,
		std::map< Size, SymDof > dofs,
		Size const score_subunit,
		utility::vector1< Size > score_multiply,
		SymSlideInfo slide_info,
		Size const num_interfaces = 1,
		std::string const & type = "simple"
	);

	SymmetryInfoOP clone() const;

	void
	initialize(
		Size const nres_monomer,
		Size const njump_monomer,
		Size const N,
		Size const num_virtual,
		std::map< Size, SymDof > dofs,
		Size const score_subunit,
		utility::vector1< Size > score_multiply,
		SymSlideInfo slide_info,
		Size const num_interfaces = 1,
		std::string const & type = "simple"
	);

	void
	initialize(
		Size const nres_monomer,
		Size const njump_monomer,
		Size const N,
		Size const num_virtual,
		std::map< Size, WtedClones > jump_clones,
		std::map< Size, SymDof > dofs,
		Size const score_subunit,
		utility::vector1< Size > score_multiply,
		SymSlideInfo slide_info,
		Size const num_interfaces = 1,
		std::string const & type = "simple"
	);

	bool operator== ( SymmetryInfo const & s );

	bool operator!= ( SymmetryInfo const & s );

	//fpd  bb_* and chi_* stuff should get merged at some point since they are always equivalent
	Size bb_follows( Size const seqpos ) const;
	Size chi_follows( Size const seqpos ) const;
	Size jump_follows( Size const seqpos ) const;

	std::vector < std::pair < Size, Size > > map_symmetric_res_pairs( Size res1, Size res2 );

	bool bb_is_independent( Size const seqpos ) const;
	bool chi_is_independent( Size const seqpos ) const;
	bool fa_is_independent( Size const seqpos ) const;
	bool jump_is_independent( Size const seqpos ) const;

	bool is_virtual( Size const seqpos ) const;

	//fpd  resize the asymm unit to contain nres_new residues
	void resize_asu( Size nres_new );

	//fpd  resize the asymm unit to contain nmonomer_new internal jumps
	void update_nmonomer_jumps( Size nmonomer_new );

	Size subunits() const;

	bool contiguous_monomers() const { return contiguous_monomers_; }

	// we rely on the user to set this properly, for the time being...
	bool torsion_changes_move_other_monomers() const { return torsion_changes_move_other_monomers_; }
	void torsion_changes_move_other_monomers( bool const setting ) { torsion_changes_move_other_monomers_ = setting;}

	Size subunit_index( Size const seqpos ) const;

	Size score_multiply_factor() const;

	utility::vector1< bool >
	independent_residues() const;

	Size num_bb_clones() const;
	Size num_chi_clones() const;
	Size num_jump_clones() const;

	Size num_independent_residues() const;
	Size num_total_residues() const;
	Size num_total_residues_with_pseudo() const;
	Size num_total_residues_without_pseudo() const;
	Size num_interfaces() const;
	Size num_virtuals() const;
	Size last_independent_residue() const;

	void num_virtuals( Size const setting );

	Size get_nres_subunit() const;  //fpd same as num_independent_residues?
	Size get_njumps_subunit() const;

	// accessors for Torsion/DOF IDs
	bool dof_is_independent( DOF_ID const & id, Conformation const & conf ) const;

	// get a weight for derivative calculations
	// weights are 1 for indep DOFs, 0 for dependent NON-JUMP DOFs
	//    and may be any real for dependent jump dofs
	core::Real get_dof_derivative_weight( DOF_ID const & id, Conformation const & conf ) const;

	bool torsion_is_independent( TorsionID const & id ) const;
	bool atom_is_independent( AtomID const & id ) const;

	/// @brief  Returns a list of dofs that depend on id. Inefficient -- creates list anew each time.
	DOF_IDs dependent_dofs( DOF_ID const & id, Conformation const & conf ) const;

	/// @brief  Returns a list of dofs that depend on id. Inefficient -- creates list anew each time.
	TorsionIDs dependent_torsions( TorsionID const & id ) const;
	AtomIDs dependent_atoms( AtomID const & id ) const;

	// clone list accessors
	Clones const & bb_clones( Size const seqpos ) const;
	Clones const & chi_clones( Size const seqpos ) const;
	Clones const & jump_clones( Size const base_jump ) const;

	void add_bb_clone( Size const base_pos, Size const clone_pos );
	void add_chi_clone( Size const base_pos, Size const clone_pos );
	void add_jump_clone( Size const base_pos, Size const clone_jump, Real const wt );

	std::map< Size, SymDof > const &get_dofs() const;

	void set_dofs( std::map< Size, SymDof > const & dofs );

	Size interface_number( Size const res1, Size const res2 ) const;

	// score multiply factors
	Real score_multiply( Size const res1, Size const res2 ) const;
	Real deriv_multiply( Size const res1, Size const res2 ) const;
	void set_score_multiply_from_subunit_factors( utility::vector1< Size > const & score_multiply_vector_subunit, Size const nres_subunit, Size const n_subunits );
	void set_score_multiply( Size const res, Size const factor );
	void set_flat_score_multiply( Size const nres, Size const factor );
	bool get_use_symmetry() const;

	void set_use_symmetry( bool setting );  //fpd  used in silent file reading(?)

	SymSlideInfo get_slide_info() const;

	//fpd  these should be unnecessary given *_is_independent & get_*_clone functions
	bool is_asymmetric_seqpos( Size const res ) const;
	Size get_asymmetric_seqpos( Size const res ) const;

	void update_score_multiply_factor();

	// io
	friend std::istream& operator>> ( std::istream & s, SymmetryInfo & symminfo );
	friend std::ostream& operator<< ( std::ostream & s, const SymmetryInfo & symminfo );

	bool read_silent_struct( std::string const & filename );
	bool write_silent_struct( std::string const & filename );

	std::string get_jump_name(Size i) const;
	Size get_jump_num(std::string i) const;
	void set_jump_name(Size i, std::string);
	Size num_slidablejumps() const;

	void
	set_multicomponent_info(
		Size const & num_components,
		utility::vector1<char> const & components,
		std::map<char,std::pair<Size,Size> > const & component_bounds,
		std::map<std::string,char> const & name2component,
		std::map<std::string,utility::vector1<char> > const & jname2component,
		std::map<std::string,utility::vector1<Size> > const & jname2subunits
	);

	Size const & get_num_components() const;
	utility::vector1<char> const & get_components() const;
	std::map<char,std::pair<Size,Size> > const & get_component_bounds() const;
	std::map<std::string,char> const & get_subunit_name_to_component() const;
	std::map<std::string,utility::vector1<char> > const & get_jump_name_to_components() const;
	std::map<std::string,utility::vector1<Size> > const & get_jump_name_to_subunits() const;

	std::pair<Size,Size> const & get_component_bounds(char c) const;
	Size get_component_lower_bound(char c) const;
	Size get_component_upper_bound(char c) const;
	char get_component_of_residue(Size ir) const;
	char get_subunit_name_to_component(std::string const & vname) const;
	utility::vector1<char> const & get_jump_name_to_components(std::string const & jname) const;
	utility::vector1<Size> const & get_jump_name_to_subunits  (std::string const & jname) const;

private:
	void
	update_contiguous_monomers();

#ifdef USEBOOSTSERIALIZE
	friend class boost::serialization::access;

	template<class Archive>
	void serialize(Archive & ar, const unsigned int version){
		ar & bb_clones_;
		ar & chi_clones_;
		ar & jump_clones_;
		ar & jump_clone_wts_;
		ar & bb_follows_;
		ar & chi_follows_;
		ar & jump_follows_;
		ar & empty_list;
		ar & nres_monomer_;
		ar & scoring_subunit_;
		ar & npseudo_;
		ar & njump_monomer_;
		ar & type_;
		ar &  interfaces_;
		ar & score_multiply_;
		ar & score_multiply_factor_;
		ar & dofs_;
		ar & use_symmetry_;
		ar & slide_info_;
	}
#endif

	// mapping from each primary jump to it's clones
	std::map< Size, Clones >   bb_clones_;
	std::map< Size, Clones >  chi_clones_;
	std::map< Size, Clones > jump_clones_;

	// a weight applied to the motion of each jump clone
	std::map< Size, Real > jump_clone_wts_;

	// these are derived from the above lists:
	std::map< Size, Size >   bb_follows_;
	std::map< Size, Size >  chi_follows_;
	std::map< Size, Size > jump_follows_;

	// silly! should make static class data
	Clones empty_list;

	Size nres_monomer_;
	Size scoring_subunit_;
	Size npseudo_;
	Size njump_monomer_;
	Size last_indep_residue_;
	std::string type_;

	// store the number of interfaces
	Size interfaces_;

	// score multiplication factors
	utility::vector1< Size > score_multiply_;

	// intra/inter subunit reweighing factors
	Real reweight_symm_interactions_;

	// total number of subunits in the entire symm complex (not just those in the model)
	// read in the 'E =' line in the symm definition file
	core::Size score_multiply_factor_;

	// store the allowed dofs
	std::map< Size, SymDof > dofs_;

	//  Toggle use of symmetry
	bool use_symmetry_;

	// Slide info
	SymSlideInfo slide_info_;

	std::map<Size,std::string> jnum2dofname_;
	std::map<std::string,Size> dofname2jnum_;

	bool contiguous_monomers_; // are the residues in each monomer sequence-contiguous
	bool torsion_changes_move_other_monomers_; // if we change a backbone torsion in one monomer, do other monomers move?

	Size num_components_;
	utility::vector1<char> components_;
	std::map<char,std::pair<Size,Size> > component_bounds_;
	std::map<std::string,char> name2component_;
	std::map<std::string,utility::vector1<char> > jname2components_;
	std::map<std::string,utility::vector1<Size> > jname2subunits_;

	utility::vector1<char> components_moved_by_jump(std::string const & jname) const;
	utility::vector1<Size> subunits_moved_by_jump(std::string const & jname) const;

}; // SymmetryInfo



} // symmetry
} // conformation
} // core

#endif
