// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/SmartAssembly.hh
/// @brief a SEWING Assembly composed of SmartSegments
/// @author frankdt (frankdt@email.unc.edu)


#ifndef INCLUDED_protocols_sewing_SmartAssembly_hh
#define INCLUDED_protocols_sewing_SmartAssembly_hh

#include <protocols/sewing/data_storage/SmartAssembly.fwd.hh>

#include <protocols/sewing/data_storage/SmartSegment.hh>
#include <protocols/sewing/data_storage/LigandSegment.fwd.hh>
#include <protocols/sewing/data_storage/LigandResidue.fwd.hh>
#include <protocols/sewing/hashing/hasher_data.hh>
#include <protocols/sewing/hashing/BasisMapGenerator.hh>
#include <protocols/sewing/hashing/Hasher.hh>
#include <core/pose/Pose.hh>
// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <numeric/random/random.hh>

namespace protocols {
namespace sewing {
namespace data_storage {

///@brief a SEWING Assembly composed of SmartSegments
class SmartAssembly : public utility::pointer::ReferenceCount {

public:

	SmartAssembly();

	SmartAssembly(hashing::SegmentVectorCOP segment_vector, core::Size window_width=1 );

	SmartAssembly(SmartAssembly const & src);

	virtual ~SmartAssembly();

	SmartAssemblyOP
	clone() const;

	void
	set_starting_segment(SmartSegmentOP start_segment, std::string start_node_vital_segments); // for Append mode

	SmartSegmentOP
	get_n_terminal_segment() const;

	void
	set_n_terminal_segment(SmartSegmentOP new_segment);

	bool
	has_segment( core::Size seg_id );

	SmartSegmentOP
	get_segment( core::Size seg_id );

	SmartSegmentOP
	get_c_terminal_segment() const;

	void
	set_c_terminal_segment(SmartSegmentOP new_segment);

	std::string
	get_start_node_vital_segments();

	void
	set_start_node_vital_segments(std::string start_node_vital_segments);

	bool
	get_last_change_was_n_terminal() const;

	void
	set_last_change_was_n_terminal(bool);

	char
	get_last_change() const;

	SmartSegmentOP
	get_last_chimaera() const;

	BasisPair
	get_last_chimaera_deleted() const;

	void
	set_partner( core::pose::PoseOP partner_pose );

	core::pose::PoseOP
	get_partner() const;

	void
	set_last_change(char);

	void
	set_last_chimaera( SmartSegmentOP );

	void
	set_last_chimaera_deleted( BasisPair );

	bool
	can_delete();

	void
	pdb_segments(  std::map< core::Size, data_storage::SmartSegmentOP >);

	std::map< core::Size, data_storage::SmartSegmentOP > &
	pdb_segments();

	std::map< core::Size, data_storage::SmartSegmentOP > &
	local_segments();

	std::map< core::Size, data_storage::SmartSegmentOP >
	const_pdb_segments() const;

	// bool
	//can_add();

	core::Size
	get_length() const;

	void
	set_length(core::Size);

	core::Size
	get_size() const;

	void
	set_size(core::Size);

	hashing::SegmentVectorCOP
	get_segment_vector() const;

	void
	set_segment_vector(hashing::SegmentVectorCOP);



	bool
	get_modifiable_terminus(char op);

	char
	get_modifiable_terminus();

	void
	set_modifiable_terminus(char modifiable_terminus);

	bool
	is_continuous() const; // false if there is a chain break. Checks between all residues in assembly.

	core::Size
	get_random_segment_id( bool n_terminus ); // n_terminus = true to return an n_terminal seg_id, false for c-terminal seg_id

	std::string
	get_forward_assembly() const;

	std::string
	get_comprehensive_forward_assembly() const;

	std::string
	get_reverse_assembly() const;

	std::string
	get_comprehensive_reverse_assembly() const;

	std::string
	get_comprehensive_assembly_string() const;
	//Additional getters for ligand stuff

	std::map< core::Size, LigandResidueCOP >
	get_partner_ligands() const;

	void
	set_partner_ligands( std::map< core::Size, LigandResidueCOP > );

	std::map< core::Size, LigandResidueOP >
	get_local_ligands();

	std::map< core::Size, LigandResidueOP >
	get_local_ligands() const;

	// these three functions are what might be thought of as "forward" moves in the graph; they should be called directly by the Mover to change the Assembly. All three return pointers to the affected segment.

	bool
	add_segment( bool n_terminus ); //this is a pure random add, and calls the more specified add_segment function

	virtual bool
	add_segment( bool n_terminus, core::Size seg_ID_to_add, core::Size res_ID_1, core::Size res_ID_2 );

	bool
	delete_segment(bool n_terminus);

	bool
	switch_segment(bool n_terminus);

	// these three are the "reverse" moves; they should not normally be called directly, but instead are called during a revert
	/*
	void
	unadd_segment();
	*/
	bool
	undelete_segment();

	void
	unswitch_segment();

	SmartSegmentOP
	recurse_revert_far_n_terminal_parent( SmartSegmentOP current_segment ); // only takes chimaeric segments

	SmartSegmentOP
	recurse_revert_far_c_terminal_parent( SmartSegmentOP current_segment ); // only takes chiameric segments

	void
	unchimerize_ligand_segment( LigandSegmentOP ligseg );


	// and then the methods that actually accept or reject a change.

	void
	revert();

	void
	dump_side_chains();

	utility::vector1< BasisPair >
	get_abbreviated_assembly();

	void
	revert_to_abbreviated_assembly( utility::vector1< BasisPair > );

	void
	pick_random_starting_segment();

	core::pose::Pose
	to_pose(std::string);


	///@brief Transforms the segments in the second basis to the coordinate frame of the first basis. Returns the transformed (second) basis segment.
	bool
	transform_segments( BasisPair basis_pair );

	bool
	chimerize(BasisPair basis_pair, bool n_terminus );

	std::pair< bool, core::Size >
	iterate_over_basis_pairs( SmartSegmentCOP segment_1, SmartSegmentCOP segment_2, bool n_terminus );

	void
	reset_chimaera_contacts( LigandResidueOP ligand, LigandSegmentOP ligand_chimaera, Basis nterm_basis, Basis cterm_basis, bool nterm );

	void
	add_segment_and_neighbors_to_local_segments( SmartSegmentCOP ref_seg, core::Size new_id=0 );

	bool
	sample_ligand();

	void
	unsample_ligand();

	///@brief NOTE This method MUST be called after pdb_segments() is set and before the starting segment is picked for this to work
	void
	load_initial_conformers( data_storage::LigandDescription ligdes );

	core::Size
	get_window_width() const;
	void
	set_window_width( core::Size );
	bool
	get_output_partner() const;
	void
	set_output_partner( bool );
	void
	reconstitute_assembly_from_string(std::string);
	void
	trim_terminal_segments(char,core::Size);

protected:

	void
	set_resID_1( core::Size );
	void
	set_resID_2( core::Size );
	core::Size
	get_resID_1() const;
	core::Size
	get_resID_2() const;
	void
	set_segID_1( core::Size );
	void
	set_segID_2( core::Size );
	core::Size
	get_segID_1() const;
	core::Size
	get_segID_2() const;

private:
	//reference libraries for assembly

	hashing::SegmentVectorCOP segment_vector_;
	//hashing::SegmentVectorOP local_segment_vector_;
	std::map< core::Size, SmartSegmentOP > local_segments_;
	std::map< core::Size, data_storage::SmartSegmentOP > pdb_segments_; //All of these will be in local segments, but not all local segments will be in here
	std::map< core::Size, LigandResidueCOP > const_ligands_;
	std::map< core::Size, LigandResidueCOP > partner_ligands_;
	std::map< core::Size, LigandResidueOP > local_ligands_; //

	//Ligand conformer data
	//Conformer contacts and owners are updated along with the main ligands, so the only thing that should differ among them should be their coordinates
	std::map< core::Size, utility::vector1< LigandResidueOP > > ligand_conformers_; //These will be loaded in (somehow) by the mover and provided to the assembly before generation
	//std::map< core::Size, utility::vector1< LigandResidueCOP > const_ligand_conformers_; //These will be loaded in (somehow) by the mover and provided to the assembly before generation
	std::pair< core::Size, LigandResidueOP > last_sampled_ligand_;

	utility::vector1<core::Size> c_terminal_segments_;
	utility::vector1<core::Size> n_terminal_segments_;

	core::Size window_width_;

	SmartSegmentOP first_segment_;
	SmartSegmentOP n_terminal_segment_;
	SmartSegmentOP c_terminal_segment_;
	//core::Size maximum_segment_count_;

	bool last_change_was_n_terminal_;
	SmartSegmentOP last_chimaera_;
	BasisPair last_chimaera_deleted_;
	SmartSegmentOP last_reverted_parent_chimaera_;
	char last_change_;

	core::Size size_;
	core::Size length_;
	SmartSegmentOP counter_segment_;

	core::Size segID_1_;
	core::Size segID_2_;
	core::Size resID_1_;
	core::Size resID_2_;
	BasisPair current_basis_pair_;
	utility::vector1< std::pair< core::Size, core::Size > > all_basis_pairs_;
	core::pose::PoseOP partner_;
	bool output_partner_;
	char modifiable_terminus_;
	std::string start_node_vital_segments_;
};


} //protocols
} //sewing
} //data_storage


#endif //INCLUDED_protocols_sewing_SmartAssembly_hh





