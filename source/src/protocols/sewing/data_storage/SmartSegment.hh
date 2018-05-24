// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/SmartSegment.hh
/// @brief a neighbor-aware SewSegment version
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_data_storage_SmartSegment_hh
#define INCLUDED_protocols_sewing_data_storage_SmartSegment_hh

#include <protocols/sewing/data_storage/SmartSegment.fwd.hh>
#include <protocols/sewing/data_storage/SmartSewingResidue.hh>
#include <protocols/sewing/data_storage/Basis.hh>
//Core Headers
#include <core/types.hh>
#include <utility/vector1.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/chemical/Atom.hh>
namespace protocols {
namespace sewing {
namespace data_storage {
///@brief a neighbor-aware SewSegment version
class SmartSegment : public utility::pointer::ReferenceCount {

public:

	SmartSegment(); // "everything is empty" default constructor. This gets used for building segments from files
	SmartSegment(bool is_vital); // "everything is empty" default constructor. This gets used for building segments from files
	SmartSegment(bool is_vital, core::Size max_segment_length); // use this constructor to actually make segments in an assembly, since it preallocates everything
	SmartSegment(SmartSegment const & src);

	virtual ~SmartSegment();

	SmartSegmentOP
	clone() const;

	void
	set_segment_id(core::Size new_segment_id);

	core::Size
	get_segment_id() const;

	void
	set_n_terminal_neighbor(SmartSegmentOP new_n_terminal_neighbor);

	SmartSegmentOP
	get_n_terminal_neighbor() const;

	void
	set_c_terminal_neighbor(SmartSegmentOP new_c_terminal_neighbor);

	SmartSegmentOP
	get_c_terminal_neighbor() const;

	void
	set_dssp_code(char new_dssp_code);

	char
	get_dssp_code() const;

	void
	set_chimaeric_status(bool is_chimaera);

	bool
	is_chimaeric() const;

	void
	set_n_terminal_parent(SmartSegmentOP new_n_terminal_parent);

	void
	set_c_terminal_parent(SmartSegmentOP new_c_terminal_parent);

	SmartSegmentOP
	get_n_terminal_parent() const;

	static SmartSegmentOP
	get_far_n_terminal_parent( SmartSegmentOP current_segment );
	//will return itself if segment has no parent.

	SmartSegmentOP
	get_c_terminal_parent() const;

	static SmartSegmentOP
	get_far_c_terminal_parent( SmartSegmentOP current_segment );
	//will return itself if segment has no parent.

	bool
	is_vital() const;

	void
	set_is_vital(bool);

	bool
	is_in_Assembly() const;

	void
	set_is_in_Assembly(bool);

	bool
	is_n_terminus_fixed() const;

	bool
	is_c_terminus_fixed() const;

	static SmartSegmentCOP
	get_n_most_segment(SmartSegmentCOP, bool cross_chimaerae);

	static SmartSegmentOP
	get_n_most_segment(SmartSegmentOP, bool cross_chimaerae);

	static SmartSegmentCOP
	get_c_most_segment(SmartSegmentCOP, bool cross_chimaerae);

	static SmartSegmentOP
	get_c_most_segment(SmartSegmentOP, bool cross_chimaerae);

	bool
	is_hashable() const;

	void
	set_residue_vector(utility::vector1<SmartSewingResidueOP> new_residue_vector);

	utility::vector1<SmartSewingResidueOP> &
	get_residue_vector();

	utility::vector1<SmartSewingResidueOP> const &
	get_const_residue_vector() const;

	SmartSewingResidueOP
	get_residue(core::Size resnum);

	void
	set_length(core::Size new_length);

	core::Size
	get_length() const; // returns number of residues, not the _length data member. Although they may or maynot be equal.

	void
	set_basis_pair(BasisPair basis_pair);

	// BasisPair
	// get_basis_pair();

	BasisPair
	get_basis_pair() const;

	std::set< core::Size >
	get_vital_residues() const;

	void
	set_vital_residues( std::set< core::Size > );

	void
	add_vital_residue( core::Size );

	void
	set_const_reference_segment( SmartSegmentCOP ref );

	SmartSegmentCOP
	get_const_reference_segment() const;


	static void
	link_to(SmartSegmentOP n_term, SmartSegmentOP c_term);

	void
	isolate();

	bool
	residue_is_vital( core::Size );

	virtual std::string
	get_name() const { return "SmartSegment"; }

	static void
	become( SmartSegmentOP changing_segment, SmartSegmentCOP src_segment );

private:
	core::Size segment_id_; // this is set to the index in the masterSegmentVector, or 0 if we don't want it in the segmentVector
	SmartSegmentOP n_terminal_neighbor_; //an OP to the segment immediately n-terminal to this one
	SmartSegmentOP c_terminal_neighbor_; //an OP to the segment immediately c-terminal to this one
	utility::vector1<SmartSewingResidueOP> residues_; // the vector actually holding our cut-down residues
	char dssp_code_; // standard S/H/L code
	bool is_a_chimaera_; //exactly what it says. Properly it should be set on construction, but we need to copy without allocating so it is resettable on the fly
	SmartSegmentOP n_terminal_parent_; //like neighbor, but for parents.
	SmartSegmentOP c_terminal_parent_;//like neighbor, but for parents.
	bool is_vital_; //do we need this segment in the assembly? Intended for use in site inclusion and Append runs. Risky to use anywhere else without a cut-down graph.
	bool is_in_Assembly_; // is this segment in an Assembly?
	core::Size length_; // separate from residues_.size() so we can preallocate residue vectors
	//core::Size current_residue_; //a counter so we don't need to allocate one. Value is not guaranteed to reflect anything in particular
	BasisPair basis_pair_ = std::make_pair( Basis( 0, 0), Basis( 0, 0) );
	std::set< core::Size > vital_residues_;
	SmartSegmentCOP const_reference_segment_ = nullptr;
};


} //protocols
} //sewing
} //data_storage


#endif //INCLUDED_protocols_sewing_SmartSegment_hh





