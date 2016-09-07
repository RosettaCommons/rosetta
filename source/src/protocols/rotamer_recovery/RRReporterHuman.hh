// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRReporterHuman.hh
/// @author Matthew O'Meara (mattjomeara@gmail.com)
///
///
///
///The rotamer recovery will be output to the screen. Output looks like:
/// # native_structure_tag1
/// # total = 30
///    resi_idx  nat_bb_bin      pct_bb    nat_rot1    pct_rot1    nat_rot2    pct_rot2    nat_rot3    pct_rot3    nat_rot4    pct_rot4
///           1           E      1.0000           1      1.0000           2      1.0000           1      1.0000         999      0.0000
///           2           B      1.0000           2      1.0000           1      1.0000         999      0.0000         999      0.0000
///    ...
///
/// # native_structure_tag2
/// # total = 30
///    resi_idx  nat_bb_bin      pct_bb    nat_rot1    pct_rot1    nat_rot2    pct_rot2    nat_rot3    pct_rot3    nat_rot4    pct_rot4
///           1           E      1.0000           1      1.0000           2      1.0000           1      1.0000         999      0.0000
///           2           B      1.0000           2      1.0000           1      1.0000         999      0.0000         999      0.0000
///    ...

/// Where the # total is how many proteins compared.
/// resi_idx = residue index
/// nat_bb_bin = dssp naming for bb
/// pct_bb = fraction matching backbone bins
/// nat_rot1 = chi 1
/// pct_rot1 = fraction matching chi bins
/// If 999 appears, that means that the amino acid does not have that chi angle


#ifndef INCLUDED_protocols_rotamer_recovery_RRReporterHuman_HH
#define INCLUDED_protocols_rotamer_recovery_RRReporterHuman_HH

// Unit Headers
#include <protocols/rotamer_recovery/RRReporterHuman.fwd.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pack/dunbrack/DunbrackRotamer.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <map>

//Auto Headers
#include <string>


namespace protocols {
namespace rotamer_recovery {

class PerNativeRRReporterHuman : public utility::pointer::ReferenceCount {

public: // public constructors / destructors

	PerNativeRRReporterHuman();

	PerNativeRRReporterHuman(
		core::pose::Pose const & native_pose
	);

	~PerNativeRRReporterHuman() override;

	PerNativeRRReporterHuman( PerNativeRRReporterHuman const & src );

private: // private helper functions

	/// @brief return one character code for region of Ramachandran plot
	char
	torsion2big_bin(
		core::Real const phi,
		core::Real const psi,
		core::Real const omega);

public: // public interface

	void
	set_native(
		core::pose::Pose const & native_pose);

	void
	report_rotamer_recovery(
		core::pose::Pose const & /*decoy_pose*/,
		core::conformation::Residue const & res,
		core::Real score,
		bool recovered);

	bool
	initialized() const;

	void
	show(
		std::ostream & out,
		core::Size column_width = 12,
		core::Size precision = 4
	) const;

private: // private member data

	bool initialized_;
	core::pose::PoseCOP native_pose_;
	utility::vector1< char > nat_bb_bins_;
	utility::vector1< core::pack::dunbrack::RotVector > nat_rots_;
	utility::vector1< core::pack::dunbrack::RotVector > nat_chis_;

	utility::vector1< utility::vector1< core::Real > > res_scores_;
	utility::vector1< utility::vector1< bool > > res_recovered_;

};


/// @brief The reporting functionality for the rotamer recovery test in
///a human readable format.
class RRReporterHuman : public RRReporter {

public: // constructors destructors

	RRReporterHuman();

	RRReporterHuman( RRReporterHuman const & src );

	~RRReporterHuman() override;

private: // private helper functions


	void
	write_header(
		std::ostream & out
	) const;

public: // public interface

	void
	set_comparer_info(
		std::string const & comparer_name,
		std::string const & comparer_params) override;

	void
	set_protocol_info(
		std::string const & protocol_name,
		std::string const & protocol_params) override;


	void
	reset_recovery() override;


	void
	report_rotamer_recovery(
		core::pose::Pose const & pose1,
		core::pose::Pose const & pose2,
		core::conformation::Residue const & res1,
		core::conformation::Residue const & res2,
		core::Real score,
		bool recovered
	) override;


	core::Real
	recovery_rate() const override;


	void
	show(std::ostream & out ) const override;


	void
	show() const override;

private:

	std::string protocol_name_;
	std::string protocol_params_;
	std::string comparer_name_;
	std::string comparer_params_;
	core::Size column_width_;
	core::Size precision_;
	std::map< std::string, PerNativeRRReporterHuman > per_native_recovery_;
	core::Real residues_considered_;
	core::Size rotamers_recovered_;
	core::Real recovery_score_mean_;
	core::Real recovery_score_m2_;
};


} // namespace rotamer_recovery
} // namespace protocols


#endif // include guard
