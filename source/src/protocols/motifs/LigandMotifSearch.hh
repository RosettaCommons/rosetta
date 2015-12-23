// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file LigandMotifSearch.hh
/// @brief Class declaration for ligand motif searching protocol
/// @author mdsmith (mdwsmith@u.washington.edu)

#ifndef INCLUDED_protocols_motifs_LigandMotifSearch_HH
#define INCLUDED_protocols_motifs_LigandMotifSearch_HH

// Unit Headers
#include <protocols/motifs/LigandMotifSearch.fwd.hh>

// Package Headers
#include <protocols/dna/DnaDesignDef.fwd.hh>
#include <protocols/motifs/BuildPosition.fwd.hh>
#include <protocols/motifs/Motif.fwd.hh>
#include <protocols/motifs/MotifLibrary.fwd.hh>

// Project Headers
#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>
#include <set>
#include <string>

#include <core/pack/task/PackerTask.fwd.hh>


namespace protocols {
namespace motifs {

class LigandMotifSearch : public utility::pointer::ReferenceCount
{

public:

	typedef core::Real Real;
	typedef core::Size Size;
	typedef utility::vector1< Size > Sizes;
	typedef core::pose::Pose Pose;
	typedef protocols::dna::DnaDesignDefOP DnaDesignDefOP;
	typedef protocols::dna::DnaDesignDefOPs DnaDesignDefOPs;
	typedef core::pack::task::PackerTask PackerTask;

	// Constructor
	LigandMotifSearch();

	// Destructor
	virtual ~LigandMotifSearch();

	// Copy constructor
	LigandMotifSearch(
		LigandMotifSearch const & src
	);

	LigandMotifSearch &
	operator = (
		LigandMotifSearch const & src
	);


	void run(
		Pose const & pose,
		utility::vector1< Size > & input_BPs
	);

	void run(
		Pose const & pose,
		PackerTask & task
	);

	void run(
		Pose const & pose,
		core::Real & ligand_motif_sphere
	);

	void initialize(
		Pose const & pose
	);

	void initialize(
		Pose const & pose,
		utility::vector1< Size > & input_BPs
	);

	void
	incorporate_motifs(
		Pose const & pose
	);

	core::pack::rotamer_set::Rotamers
	bp_rotamers(
		Size const seqpos
	);

	core::pack::rotamer_set::Rotamers
	get_rotamers();

	bool
	protein_dna_motif();

	void
	position_vector_setup(
		Pose const & pose
	);

	void
	identify_motif_build_positions(
		Pose const & pose,
		Sizes & build_positions
	);

	utility::vector1< Size >
	get_sphere_aa(
		Pose const & pose,
		core::Real cut1
	);

	/*void
	fill_bp_allowed_types(
	Pose const & pose,
	Size const spos,
	std::set< std::string > & allowed_types
	);*/

	void identify_motif_BuildPositions(
		Pose const & pose
	);

	void
	BuildPosition_from_Size(
		Pose const & pose,
		Size const input_BP
	);

	void
	defs2BuildPositions(
		Pose const & pose,
		protocols::dna::DnaDesignDefOPs const & defs
	);

	void
	defs2BuildPositions_findts(
		Pose const & pose,
		protocols::dna::DnaDesignDefOPs const & defs
	);

	utility::vector1< Size >
	map2keyvector(
		std::map< Size, std::set< std::string > > mappositions
	);

	utility::vector1< Size >
	shorten_target_list(
		Pose const & pose,
		Size const bp,
		Sizes & full_tl
	);

	void
	protein_DNA_motif_build_positions_JA(
		Pose const & pose,
		Sizes & build_positions,
		Sizes & target_positions
	);

	void
	protein_DNA_motif_target_positions_JA(
		Pose const & pose,
		Sizes & build_positions,
		Sizes & target_positions,
		Sizes & short_tl
	);

	void
	override_option_input(
		Real const & r1,
		Real const & z1,
		Real const & r2,
		Real const & z2,
		Real const & d1,
		Size const & rlevel
	);

	void
	reset_option_input();

	void
	set_motif_library(
		MotifLibrary & motiflibrary
	);

	// Accessors
	MotifCOPs const & motif_library() const { return motif_library_; }
	Sizes const & dna_positions() const { return dna_positions_; }
	Sizes const & protein_positions() const { return protein_positions_; }
	std::map< Size, std::set< std::string > > const & target_positions() const { return target_positions_; }
	BuildPositionOPs const & build_positionOPs() const { return build_positionOPs_; }
	std::map< std::string, core::conformation::ResidueOPs > const & target_conformers_map() const { return target_conformers_map_; }
	Real const & ztest_cutoff_1() const { return ztest_cutoff_1_;}
	Real const & ztest_cutoff_2() const { return ztest_cutoff_2_;}
	Real const & rmsd_cutoff_1() const { return rmsd_cutoff_1_;}
	Real const & rmsd_cutoff_2() const { return rmsd_cutoff_2_;}
	Real const & dtest_cutoff() const { return dtest_cutoff_;}
	Size const & rot_level() const { return rot_level_;}
	bool const & minimize() const { return minimize_;}

private:
	void init_options();

private:
	MotifCOPs motif_library_;
	Sizes dna_positions_;
	Sizes protein_positions_;
	std::map< Size, std::set< std::string > > target_positions_;
	BuildPositionOPs build_positionOPs_;
	std::map< std::string, core::conformation::ResidueOPs > target_conformers_map_;

	// inputs used for scoring motif, 1s are broader first limit and 2s are a tighter cutoff for final selection
	Real ztest_cutoff_1_;
	Real ztest_cutoff_2_;
	Real rmsd_cutoff_1_;
	Real rmsd_cutoff_2_;
	Real dtest_cutoff_;

	// input rotamer level (4 = all chis explosion, 8 is extra deviations of explosion on all 4 chis)
	Size rot_level_;

	bool minimize_;

	// options set in init_options
	bool bpdata_;
	std::string bpdata_filename_;
	bool output_;
	std::string output_filename_;
	bool data_;
	std::string data_filename_;
	bool quick_and_dirty_;
	bool dump_motifs_;
	bool clear_bprots_;
	Size rots2add_;

};

} // namespace motifs
} // namespace protocols

#endif // INCLUDED_protocols_motifs_LigandMotifSearch
