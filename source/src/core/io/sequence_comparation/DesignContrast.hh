// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
// (C) 199x-2008 University of Washington
// (C) 199x-2008 University of California Santa Cruz
// (C) 199x-2008 University of California San Francisco
// (C) 199x-2008 Johns Hopkins University
// (C) 199x-2008 University of North Carolina, Chapel Hill
// (C) 199x-2008 Vanderbilt University

/// @file   core/io/sequence_comparation/DesignContrast.hh
///
/// @brief
/// @author Yi Liu

#ifndef INCLUDED_core_io_sequence_comparation_DesignContrast_hh
#define INCLUDED_core_io_sequence_comparation_DesignContrast_hh

// Unit headers
#include <core/io/sequence_comparation/DesignContrast.fwd.hh>

// mini headers
#include <core/pose/Pose.fwd.hh>

#include <utility/file/FileName.hh>

#include <utility/vector1.hh>


namespace core {
namespace io {
namespace sequence_comparation {

using utility::file::FileName;
using utility::vector1;
/// @brief DesignContrast contains information for comparing the native protein sequence to
/// designed protein sequence. And output the compare resultes to a special formated file which
/// can be used for statistics calculations

class DesignContrast {
public :
	/// @brief default constructor
	DesignContrast (){ }

	// @brief copy constructor
	DesignContrast (DesignContrast const & dc);

	/// @brief default de-constructor
	virtual ~DesignContrast(){
		clear();
	}
	/// @brief Set number of neighbors for all residues in pose
	void setNeighbors(pose::Pose & pose);

	/// @brief Get number of neighbors for all residues in pose
	vector1<int> & getNeighbors();
	vector1<int> const & getNeighbors() const;

	/// @brief Set secondary structure for all residues in pose
	void setSecStruct(pose::Pose & pose);

	/// @brief Get secondary structure for all residues in pose
	vector1<std::string> & getSecStruct();
	vector1<std::string> const & getSecStruct() const;

	/// @brief Get pdb file names from the pdb list files.
	void setNames (); // vector1<std::string> & pdb_file_names ); the pdb_file_names is a private member


	vector1<FileName> & getPdbNames();
	vector1<FileName> const & getPdbNames() const;

	vector1<FileName> & getListNames();
	vector1<FileName> const & getListNames() const;

	void setPdbCodes();

	vector1<std::string> & getPdbCodes();
	vector1<std::string> const & getPdbCodes() const;

	/// @brief this function will output the sequence comparing result between native pose and designed pose
	void output_sqc_file (
		pose::Pose & native_pose,
		pose::Pose & decoy_pose,
		std::string const & single_code,
		std::ofstream & sqc
	);

	/// @brief clear function to clear all datas in this class.
	void clear();

private:
	vector1<FileName> list_file_names_;
	//vector1<std::string> pdb_file_names_;
	vector1<FileName> pdb_file_names_;
	vector1<int> nneighbs_;
	vector1<std::string> secstructs_;
	vector1<std::string> pdb_codes_;
};
} // namespace sequence_comparation
} // namespace io
} // namespace core

#endif //INCLUDED_core_io_sequence_comparation_DesignContrast_HH
