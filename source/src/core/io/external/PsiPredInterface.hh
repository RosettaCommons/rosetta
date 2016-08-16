// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/filters/PsiPredInterface.hh
/// @brief Class that interfaces with psipred for secondary structure prediction
/// @details
/// @author Tom Linsky (tlinsky@uw.edu)

#ifndef INCLUDED_core_io_external_psipredinterface_hh
#define INCLUDED_core_io_external_psipredinterface_hh

// Unit Headers

// Package headers

// Project headers
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <map>
#include <string>

namespace core {
namespace io {
namespace external {

struct PsiPredResult {
	PsiPredResult() :
		nres( 0 ),
		pred_ss( "" )
	{
		psipred_prob.clear();
		psipred2_confidence.clear();
	}
	core::Size nres;
	std::string pred_ss;
	utility::vector1< core::Real > psipred_prob;
	utility::vector1< core::Size > psipred2_confidence;
};
typedef std::map< std::string const, PsiPredResult const > PsiPredResultMap;

/// @brief generates scores for each residue based on psipred confidence and the desired secondary structure
utility::vector1< core::Real >
generate_prob( PsiPredResult const & psipred_result, std::string desired_ss );

/// @brief finds which positions in the blueprint and predicted secondary structure string don't match
/// returns a vector of the residue numbers
utility::vector1< core::Size >
nonmatching_residues( std::string const & blueprint_ss, std::string const & pred_ss );

class PsiPredInterface : public utility::pointer::ReferenceCount {
public:
	PsiPredInterface( std::string const & cmd );

	PsiPredInterface( PsiPredInterface const & rval );

	/// @brief Runs psipred, as specified by "cmd" on the fasta file "fasta_filename", and returns the psipred filename
	/// if cached psipred output has already been saved for the AA sequence of pose, the cached output will be returned.
	PsiPredResult
	run_psipred( core::pose::Pose const & pose, std::string const & blueprint_ss );

#ifdef CXX11  // in its current form PsiPredResultMap is not assignable but compiler tries to generate assigment operator anyway
	PsiPredInterface & operator= ( const PsiPredInterface & ) = delete;
#endif

	//private member functions
private:
	/// @brief converts a pose into a fasta string
	std::string convert_to_fasta( std::string const & pname, std::string const & seq ) const;

	/// @brief dumps fasta sequence for the pose into a file and returns the filename
	/// or exits rosetta if not successful
	std::string create_fasta_file( std::string const & pname, std::string const & seq ) const;

	/// @brief deletes files created by psipred and convert_to_fasta
	void cleanup_after_psipred( std::string const & psipred_filename, std::string const & fasta_filename ) const;

	/// @brief Parses the psipred output and returns the predicted secondary structure and likelihoods of the blueprint secondary structure being present on the pose at each position.
	PsiPredResult
	parse_psipred_output( std::string const & psipred_str,
		std::string const & blueprint_ss,
		std::string const & psipred_horiz_filename ) const;

	/// @brief Parses the psipred pass2 horiz file and returns predicted secondary structure and confidence for each position.
	std::pair< std::string, utility::vector1< core::Size > >
	parse_psipred_horiz_output( std::string const & psipred_horiz_filename ) const;

private:
	/// @brief This class will cache and store Psipred output for each amino acid sequence string. They are stored in this map, which is indexed by sequence string.
	PsiPredResultMap psipred_output_;
	/// @brief The PsiPred command
	std::string cmd_;
};

} //core
} //io
} //external
#endif
