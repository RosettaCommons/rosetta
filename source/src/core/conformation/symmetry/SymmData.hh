// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief  Symmetry data container
/// @file   core/conformation/symmetry/SymmData.hh
/// @author Ingemar Andre


#ifndef INCLUDED_core_conformation_symmetry_SymmData_hh
#define INCLUDED_core_conformation_symmetry_SymmData_hh


// Unit headers
#include <core/conformation/symmetry/SymmData.fwd.hh>
#include <core/conformation/symmetry/VirtualCoordinate.hh>
#include <core/conformation/symmetry/SymDof.hh>
#include <core/conformation/symmetry/SymSlideInfo.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/types.hh>

// C++ headers
#include <map>

#include <utility/vector1.hh>
#include <numeric/xyzMatrix.hh>


namespace core {
namespace conformation {
namespace symmetry {

class SymmData : public utility::pointer::ReferenceCount
{

	// typedef utility::vector1< Size > Clones;
	typedef utility::vector1< std::pair<Size,Real> > WtedClones;

public:

	SymmData();
	SymmData( core::Size nres, core::Size njump );
	SymmData( SymmData const &);
	SymmDataOP
	clone() const;
	virtual ~SymmData();

	friend
	bool
	operator==(SymmData const & a, SymmData const & b);

	friend
	bool
	operator!=(SymmData const & a, SymmData const & b);

private:

	std::string symmetry_name_;
	std::string symmetry_type_;
	core::Size subunits_;
	core::Size num_components_;
	core::Size interfaces_;
	core::Size score_subunit_;
	std::string anchor_residue_;
	bool recenter_;
	core::Size root_;
	SymSlideInfo slide_info_;
	std::vector< std::string > slide_order_string_;
	std::vector< std::vector< std::string > > symm_transforms_;
	std::vector< numeric::xyzMatrix< core::Real> > rotation_matrices_;
	std::vector< numeric::xyzMatrix< core::Real> > translation_matrices_;
	std::map< std::string, VirtualCoordinate > virtual_coordinates_;
	std::map< std::string, std::pair< std::string, std::string > > jump_string_to_virtual_pair_;
	std::map< std::string, Size > jump_string_to_jump_num_;
	std::map< std::string, Size > virt_id_to_virt_num_;
	std::map< std::string, Size > virt_id_to_subunit_num_;
	std::map< std::string, char > virt_id_to_subunit_chain_;
	std::map< std::string, std::string > virt_id_to_subunit_residue_;
	std::map< Size, std::string > virt_num_to_virt_id_;
	std::map< Size, std::string > subunit_num_to_virt_id_;
	std::map< Size, WtedClones > jump_clones_;
	std::map< Size, SymDof > dofs_;
	std::vector< Size > allow_virtual_;
	utility::vector1< Size > score_multiply_subunit_;
	utility::vector1< Size > include_subunit_;
	utility::vector1< Size > output_subunit_;

	utility::vector1<char> components_;
	std::map<std::string,char> name2component_;
	std::map<std::string,utility::vector1<char> > jname2components_;
	std::map<std::string,utility::vector1<Size> > jname2subunits_;

	utility::vector1<std::string> leaves_of_jump(std::string const & jname) const;

public:
	typedef numeric::xyzVector< core::Real > Vector;
	typedef numeric::xyzMatrix< core::Real > Matrix;

	void
	read_symmetry_data_from_file(
		std::string filename
	);

	void
	read_symmetry_data_from_stream(
		std::istream & infile
	);

	void

	sanity_check();
	void
	show();

	// Accessor functions
	std::string const &
	get_symmetry_name() const;

	std::string const &
	get_symmetry_type() const;

	core::Size
	get_subunits() const;

	core::Size get_num_components() const;

	core::Size
	get_interfaces() const;

	core::Size
	get_score_subunit() const;

	std::string const &
	get_anchor_residue() const;

	bool
	get_recenter() const;

	core::Size
	get_root() const;

	utility::vector1< Size > const &
	get_score_multiply_subunit() const;

	utility::vector1< Size > const &
	get_include_subunit() const;

	utility::vector1< Size > const &
	get_output_subunit() const;

	std::vector< numeric::xyzMatrix< core::Real > > const &
	get_rotation_matrix() const;

	std::vector< numeric::xyzMatrix< core::Real > > const &
	get_translation_matrix() const;

	std::map< std::string, VirtualCoordinate > const &
	get_virtual_coordinates() const;

	core::Size
	get_num_virtual() const;

	std::map< Size, SymDof > const &
	get_dofs() const;

	std::map< Size, WtedClones > const &
	get_jump_clones() const;

	std::map< std::string, Size > const &
	get_jump_string_to_jump_num() const;

	std::map< std::string, Size > const &
	get_virtual_id_to_num() const;

	std::map< std::string, Size > const &
	get_virt_id_to_subunit_num() const;

	std::map< std::string, char > const &
	get_virt_id_to_subunit_chain() const;

	std::map< std::string, std::string > const &
	get_virt_id_to_subunit_residue() const;

	std::map< Size, std::string > const &
	get_subunit_num_to_virt_id() const;

	std::map< Size, std::string > const &
	get_virtual_num_to_id() const;

	std::map< std::string, std::pair< std::string, std::string > > const &
	get_virtual_connects() const;

	SymSlideInfo const &
	get_slide_info() const;

	void
	set_symmetry_name(
		std::string symm_name );

	void
	set_symmetry_type(
		std::string symm_type );

	void
	set_subunits(
		core::Size num_subunits );

	void
	set_interfaces(
		core::Size interfaces );

	void
	set_anchor_residue(
		std::string anchor
	);

	void
	set_score_multiply_subunit( utility::vector1< Real > & score_multiply_vector );

	void
	set_slide_info( SymSlideInfo slide_info );

	void
	set_rotation_matrix(
		std::vector< numeric::xyzMatrix< core::Real > > rotation_matrices );

	void
	set_translation_matrix(
		std::vector< numeric::xyzMatrix< core::Real > > translation_matrices );

	void
	set_symm_transforms(
		std::vector< std::vector< std::string> > symm_transforms );

	utility::vector1<char> const & get_components() const;
	std::map<std::string,char> const & get_subunit_name_to_component() const;
	std::map<std::string,utility::vector1<char> > const & get_jump_name_to_components() const;
	std::map<std::string,utility::vector1<Size> > const & get_jump_name_to_subunits() const;

	std::string const & get_parent_jump(std::string const & jname) const;
	std::string const & get_parent_virtual(std::string const & vname) const;
	bool is_ancestor_virtual(std::string const & ancestor, std::string const & child) const;
	std::string const & get_least_common_ancertor_jump(std::string const & jname1, std::string const & jname2) const;
	std::string const & get_least_common_ancertor_virtual(std::string const & vname1, std::string const & vname2) const;

	utility::vector1<char> components_moved_by_jump(std::string const & jname) const;
	utility::vector1<Size> subunits_moved_by_jump(std::string const & jname) const;

};

} // symmetry
} // conformation
} // core
#endif
