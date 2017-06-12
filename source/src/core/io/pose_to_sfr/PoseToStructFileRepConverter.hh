// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/pose_to_sfr/PoseToStructFileRepConverter.hh
/// @brief Headers for class to convert a pose to a StructFileRep.
/// @details This conversion is a first step in PDB or mmCIF output.  It could be useful for other
/// input/output, too.
/// @author Vikram K. Mulligan (vmullig@uw.edu), XRW 2016 Team
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com), XRW 2016 Team

#ifndef INCLUDED_core_io_pose_to_sfr_PoseToStructFileRepConverter_hh
#define INCLUDED_core_io_pose_to_sfr_PoseToStructFileRepConverter_hh

// Unit headers
#include <core/io/pose_to_sfr/PoseToStructFileRepConverter.fwd.hh>
#include <core/io/pose_to_sfr/PoseToStructFileRepConverterTests.fwd.hh>  // needed for friendship

// Package headers
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/io/Remarks.fwd.hh>
#include <core/io/ResidueInformation.fwd.hh>

// Project headers
#include <core/types.hh>
#include <core/id/AtomID_Mask.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/PDBInfo.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>
#include <map>


namespace core {
namespace io {
namespace pose_to_sfr {

/// @brief A class to convert a pose to a StructFileRep.
/// @details The StructFileRep is a data structure that serves as an intermediate between
/// Protein Databank file formats (PDB and mmCIF) and Rosetta objects (poses).  This was implemented
/// at the 2016 Chemical XRW (eXtreme Rosetta Workshop), at which this class was popularly referred
/// to as the "pose unbuilder".
class PoseToStructFileRepConverter : public utility::pointer::ReferenceCount
{

	friend class ::PoseToStructFileRepConverterTests; //Needed to allow the unit tests to access private member functions.

public:

	/// @brief Constructor.
	/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
	/// (using the sfr() object).  The options_ object is created and initialized from the options system.
	PoseToStructFileRepConverter();

	/// @brief Constructor with options input.
	/// @details Creates the StructFileRep object, which is subsequently accessible by owning pointer
	/// (using the sfr() object).  The options_ object is duplicated from the input options object.
	PoseToStructFileRepConverter( StructFileRepOptions const &options_in );

	/// @brief Destructor.
	~PoseToStructFileRepConverter() {};

	/// @brief Resets the PoseToStructFileRepConverter object, and reinitializes
	/// it with a fresh StruftFileRep object, returning an owning pointer to the
	/// new object.
	core::io::StructFileRepOP new_sfr();

	/// @brief Fill StructFileRep object using information from given Pose object.
	void init_from_pose(core::pose::Pose const & pose );

	/// @brief Fill StructFileRep object  using information from given Pose object and a set of options.
	void init_from_pose(core::pose::Pose const & pose, StructFileRepOptions const & options);

	/// @brief Fill StructFileRep object using information from given Pose object,
	///  for a specified subset of residues
	void init_from_pose( core::pose::Pose const & pose, utility::vector1< core::Size > const & residue_indices );

	/// @brief Fill StructFileRep object using information from given Pose object,
	///  for a specified subset of residues and atoms
	void init_from_pose( core::pose::Pose const & pose, id::AtomID_Mask const & mask);

	/// @brief Fill StructFileRep object using information from given Pose object,
	///  for a specified subset of residues and atoms
	void init_from_pose( core::pose::Pose const & pose, id::AtomID_Mask const & mask, StructFileRepOptions const & options);

	/// @brief Get additional pose data.
	/// @details This is rewritten from the old core/io/pdb/file_data.cc:write_additional_pdb_data() function.
	/// This was a catch-all for dumping out all sorts of protocol-specific stuff (e.g. membrane info).
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void grab_additional_pose_data( core::pose::Pose const & pose );


	/// @brief Append pdb information to StructFileRep for a single residue.
	void append_residue_to_sfr(
		pose::Pose const & pose,
		core::Size const resnum );


	///  @brief Start atom numbering from given n.  Increments N
	void append_residue_to_sfr(
		pose::Pose const & pose,
		core::Size const resnum,
		core::Size & new_atom_num,
		core::Size const new_tercount );

	/// @brief Get connectivity annotation information from the Pose object and create LinkInformation and
	/// SSBondInformation data as appropriate.
	void get_connectivity_annotation_info( core::pose::Pose const & pose, core::pose::PDBInfoOP const & fresh_pdb_info = nullptr );

	LinkInformation get_link_record( core::pose::Pose const & pose, core::Size ii, core::Size conn, core::pose::PDBInfoOP const & fresh_pdb_info = nullptr );
	SSBondInformation get_ssbond_record( core::pose::Pose const & pose, core::Size ii, core::Size conn, core::pose::PDBInfoOP const & fresh_pdb_info = nullptr );

	/// @brief Get parametric information from the Pose object and add it to the PDB remarks.
	void get_parametric_info( core::io::RemarksOP remarks, core::pose::Pose const & pose );

	/// @brief Grab all the data that makes the pose energies table
	void grab_pose_energies_table( core::pose::Pose const & pose);

	/// @brief Grab all the data that is in pose datacache as string/value pairs
	///  Arbitrary Float data
	///  Arbitrary String data
	void grab_pose_cache_data( core::pose::Pose const & pose);


	/// @brief Debug printing
	friend std::ostream& operator <<(std::ostream &os, StructFileRep const &);


	/// @brief Non-const access to the StructFileRep object.
	///
	StructFileRepOP sfr();

private:

	/// @brief Should this atom be [added]?
	bool add_atom_to_sfr(
		pose::Pose const & pose,
		conformation::Residue const & rsd,
		core::Size atom_index,
		bool use_pdb_info ) const;

	bool use_pdb_info_for_num( pose::Pose const & pose, Size resnum );


	/// @brief Append just residue-based info to StructFileRep
	///  For now, it is only HETNAM Data
	void append_residue_info_to_sfr(
		ResidueInformation const & res_info,
		core::conformation::Residue const & rsd );

	/// @brief Append just atom-based info to StructFileRep.
	///  New atom number will be detected from current SFR or start from 1.
	///  Returns True if successfully added to SFR
	bool append_atom_info_to_sfr(
		core::pose::Pose const & pose,
		ResidueInformation const & res_info,
		core::conformation::Residue const & rsd,
		core::Size const atom_index,
		bool const use_pdb_info);

	/// @brief Append just atom-based info to StructFileRep
	///   Specify the new atom to start from.
	///   Returns True if successfully added to SFR
	bool append_atom_info_to_sfr(
		core::pose::Pose const & pose,
		ResidueInformation const & res_info,
		core::conformation::Residue const & rsd,
		core::Size const atom_index,
		bool const use_pdb_info,
		core::Size const new_atom_num,
		core::Size const new_tercount);

	core::Size get_new_atom_serial_num() const;

private: //PRIVATE FUNCTIONS:

	/// @brief Set whether to write the fold tree, in the
	/// StructFileRepOptions object (options_).
	void set_fold_tree_io( bool const setting);

	/// @brief Get the membrane information from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	/// @param [in] pose The pose.
	/// @param [in] normalize_to_thk Normalized MEM lines, useful for visualizing the boundaries
	/// of the membrane by coupling the NORM and THK coordinates.  Added by Rebecca.
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	grab_membrane_info(
		core::pose::Pose const &pose,
		bool const normalize_to_thk
	);

	/// @brief Get the conect record information from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	/// @param [in] pose The pose, for determining what's bonded to what.
	/// @param [in] atom_index_in_pose The index (in the whole pose) of the atom that we're setting up.
	/// @param [in] res_index The current residue index.
	/// @param [in] atom_index_rsd The index (in the residue) of the atom that we're setting up.
	/// @param [in] dist_cutoff The atom separation, above which a CONECT record is written.
	/// @param [in] write_virtuals Are virtual atoms being written out?
	/// @param [in/out] ai The AtomInformation object that's being set up.  Modified by this function
	/// @author Vikram K. Mulligan (vmullig@uw.edu)
	void
	grab_conect_records_for_atom(
		core::pose::Pose const &pose,
		core::Size const res_index,
		core::Size const atom_index_in_rsd,
		core::io::AtomInformation &ai
	) const;

	/// @brief Get the foldtree from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	void grab_foldtree( core::pose::Pose const &pose, bool const output_foldtree );

	/// @brief Get the pdb parents from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	void grab_pdb_parents( core::pose::Pose const &pose, bool const output_parents );

	/// @brief Get the pdb comments from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	void grab_pdb_comments( core::pose::Pose const &pose, bool const output_comments );

	/// @brief Get the torsion information from the pose and store it in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	/// @details This is an incredibly ugly way to handle the output of these
	/// records.  This abuses the REMARK lines in the PDB format, and should be
	/// refactored.  Keeping as-is for now to preserve legacy behaviour (VKM
	/// 29 January 2016, Chemical XRW).
	void grab_torsion_records( core::pose::Pose const &pose, bool const output_torsions );

	/// @brief Get the pdbinfo labels from the pose and store them in the
	/// StructFileRep for output to pdbs/mmCIF/whatnot.
	void grab_pdbinfo_labels( core::pose::Pose const &pose );

	/// @brief Get the total number of atoms in the SFR.
	///
	core::Size
	total_sfr_atoms( StructFileRep const & sfr) const;

	/// @brief Get the total number of atoms in a chain in the SFR.
	///
	core::Size
	total_sfr_atoms( StructFileRep const & sfr, core::Size const chain_num ) const;

	/// @brief Return the PDB resName, chainID, resSeq, and iCode for the given Rosetta sequence position.
	/// @details Output is res_info.
	void
	get_residue_information(
		core::pose::Pose const & pose,
		core::uint const seqpos,
		bool const use_PDB/*=true*/,
		bool const renumber_chains/*=false*/,
		core::Size const new_tercount,
		ResidueInformation & res_info ) const;

	void
	get_residue_information(
		core::pose::Pose const & pose,
		core::uint const seqpos,
		bool const use_PDB/*=true*/,
		bool const renumber_chains/*=false*/,
		core::Size const new_tercount,
		ResidueInformation & res_info,
		core::pose::PDBInfoOP & fresh_pdb_info
	) const;

	/// @brief fills HELIXInformation and SHEETInformation for SFR
	void generate_secondary_structure_informations(
		core::pose::Pose const & pose
	);

	/// @brief fills one HELIXInformation into SFR
	void generate_HELIXInformation(
		ResidueInformation const & start_info,
		ResidueInformation const & stop_info,
		core::Size const index,
		core::Size const length
	);

	/// @brief fills one SHEETInformation into SFR
	void generate_SHEETInformation(
		ResidueInformation const & start_info,
		ResidueInformation const & stop_info,
		core::Size const index
	);

private: // PRIVATE DATA:

	/// @brief Owning pointer to the StructFileRep object, copied in during
	/// object construction.
	StructFileRepOP sfr_;

	/// @brief The StructFileRepOptions object.
	/// @details Created at PoseToStructFileRepConverter object creation; not read-accessible outside of this object.
	StructFileRepOptions options_;

}; // class PoseToStructFileRep

// Left over from pose energies table.  Remove if possible! JAB
std::string restrict_prec( core::Real inval );

} // namespace pose_to_sfr
} // namespace io
} // namespace core

#endif
