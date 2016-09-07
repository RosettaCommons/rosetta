// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/environment/ProtectedConformation.hh
/// @brief A conformation built by the environment to be protected.
///
/// @author Justin R. Porter

#ifndef INCLUDED_protocols_environment_ProtectedConformation_hh
#define INCLUDED_protocols_environment_ProtectedConformation_hh

// Unit Headers
#include <protocols/environment/ProtectedConformation.fwd.hh>

// Package headers
#include <protocols/environment/DofUnlock.fwd.hh>
#include <protocols/environment/Environment.fwd.hh>
#include <core/environment/SequenceAnnotation.hh>

// Project headers
#include <core/conformation/Conformation.hh>

#include <core/kinematics/FoldTree.fwd.hh>

#include <core/environment/DofPassport.hh>

#include <core/id/TorsionID.fwd.hh>

// C++ Headers
#include <stack>

// ObjexxFCL Headers

namespace protocols {
namespace environment {

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core::conformation;

class ProtectedConformation : public core::conformation::Conformation {
	friend class DofUnlock;
	friend class Environment;
	friend class EnvClaimBroker;

	typedef core::conformation::Conformation Parent;
	typedef core::environment::SequenceAnnotationCOP SequenceAnnotationCOP;

public:
	ProtectedConformation( ProtectedConformation const& src );
	ProtectedConformation( EnvironmentCAP, core::conformation::Conformation const& );

	~ProtectedConformation() override;

	//Environment-related functions
	EnvironmentCAP environment() const;

	void env_destruction(){
		environment_exists_ = false;
	}

	bool is_protected() const override { return true; }

	//@brief return a copy of in_rsd with variants matched to seqpos in this conformation
	ResidueOP match_variants( core::Size seqpos, Residue const& in_rsd ) const;

	//Annotation functions:
	SequenceAnnotationCOP resolver() const;

	SequenceAnnotationCOP annotations() const { return resolver(); }

	core::conformation::ConformationOP clone() const override;

	//Security overloads:
	Conformation& operator=( Conformation const& src ) override;

	void set_torsion( TorsionID const & id, core::Real setting ) override;

	void set_jump( AtomID const&, core::kinematics::Jump const& ) override;

	void set_jump( int, core::kinematics::Jump const& ) override;

	void set_secstruct( Size seqpos, char setting ) override;

	void replace_residue( Size seqpos, core::conformation::Residue const & new_rsd,
		utility::vector1< std::pair< std::string, std::string > > const& atom_pairs ) override;

	void replace_residue( Size seqpos, core::conformation::Residue const& new_rsd, bool orient_backbone ) override;

	virtual void set_stub_transform( core::id::StubID const & stub_id1,
		core::id::StubID const & stub_id2,
		core::kinematics::RT const & target_rt );

	void set_dof( DOF_ID const& id, core::Real setting ) override;

	void set_torsion_angle( AtomID const & atom1,
		AtomID const & atom2,
		AtomID const & atom3,
		AtomID const & atom4,
		core::Real setting,
		bool quiet=false) override;

	void set_bond_angle( AtomID const & atom1, AtomID const & atom2, AtomID const & atom3,
		core::Real setting ) override;

	void set_bond_length( AtomID const & atom1, AtomID const & atom2, core::Real setting ) override;

	virtual void insert_fragment( core::id::StubID const& instub_id, FragRT const& outstub_transforms,
		FragXYZ const& frag_xyz );

	// Always-failing Security Overloads
	// The parameters aren't named on purpose. You're not supposed to use these functions.
	void fold_tree( FoldTree const& ) override;

	virtual void chain_endings( utility::vector1< Size > const& );

	virtual void insert_chain_ending( Size );

	virtual void delete_chain_ending( Size );

	virtual void reset_chain_endings();

	virtual void chains_from_termini();

	void append_residue_by_jump( core::conformation::Residue const &, Size ,
		std::string const& = "", std::string const& = "",
		bool _= false ) override;

	virtual void append_polymer_residue_after_seqpos( core::conformation::Residue const&, Size ,
		bool );

	virtual void safely_append_polymer_residue_after_seqpos( core::conformation::Residue const&, Size ,
		bool);

	virtual void prepend_polymer_residue_before_seqpos( core::conformation::Residue const&, Size ,
		bool);

	virtual void safely_prepend_polymer_residue_before_seqpos( core::conformation::Residue const&,
		Size , bool);

	virtual void delete_polymer_residue( Size );

	virtual void delete_residue_slow( Size );

	virtual void delete_residue_range_slow( Size range_begin, Size range_end );

	void declare_chemical_bond( Size, std::string const&, Size, std::string const& ) override;

	void insert_conformation_by_jump( Conformation const&, Size , Size, Size,
		Size, std::string const& = "",  std::string const& = "" ) override;

	virtual void rebuild_polymer_bond_dependent_atoms( Size );

	virtual void insert_ideal_geometry_at_polymer_bond( Size seqpos );

	virtual void insert_ideal_geometry_at_residue_connection( Size pos1, Size connid1 );

	virtual void set_polymeric_connection( Size, Size );

	// TODO: decide what to do with disulfides
	virtual void fix_disulfides( utility::vector1< std::pair<Size, Size> > );

	// TODO: this method can fail in fullatom
	// This just uses the parent method anyway, and we need a default argument.
	//virtual void detect_disulfides();

	// TODO: decide what to do with direct xyz settings.
	void set_xyz( AtomID const & id, core::PointPosition const & position ) override;
	void batch_set_xyz( utility::vector1<AtomID> const & id,
		utility::vector1< core::PointPosition > const & position ) override;


	virtual void reset_move_data();

	virtual void clear();

	virtual void fill_missing_atoms( core::id::AtomID_Mask missing );

	// I decided not to overload the following functions because they seem safe for anyone to call.
	// If that's wrong, uncomment, implement an always-fail (or something smarter!), and make them
	// virtual in Conformation.hh

	// virtual void update_actcoords();
	// virtual void update_actcoord( Size resid );
	// virtual void update_orbital_coords( Size resid );

	// virtual void update_polymeric_connection( Size );
	// virtual void detect_bonds();
	// virtual void detect_pseudobonds();

	//Misc overloads:
	bool same_type_as_me( Conformation const & other, bool recurse /* = true */ ) const override;

	// Verification Helpers:
private:
	inline bool verify( core::id::TorsionID const& );

	inline bool verify( core::id::DOF_ID const& );

	inline bool verify_jump( core::id::AtomID const& );

	inline bool verify_backbone( Size const& seqpos );

	//used to build a bogus DOF_ID, but can include sequence position
	inline void fail_verification( std::string const& mod_type,
		core::Size const& seqpos = 0 ) const;

	inline void fail_verification( core::id::DOF_ID const& id,
		std::string const& mod_type ) const;

	void replace_residue_sandbox( Size seqpos, core::conformation::Residue const& new_rsd, bool );

	// Passport Management:
private:
	void push_passport( core::environment::DofPassportCOP ) override;

	core::environment::DofPassportCOP pop_passport() override;

	bool has_passport() const override;

	void set_environment( EnvironmentCAP );

	void attach_annotation( SequenceAnnotationCOP );

	std::stack<core::environment::DofPassportCOP> unlocks_;
	SequenceAnnotationCOP annotations_;
	EnvironmentCAP env_;

	bool environment_exists_;

}; // end ProtectedConformation base class

} // environment
} // protocols

#endif //INCLUDED_protocols_environment_ProtectedConformation_hh
