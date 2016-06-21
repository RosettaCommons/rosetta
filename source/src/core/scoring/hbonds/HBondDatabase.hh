// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/geometric_solvation/HBondDatabase.cc
/// @brief  Database containing params for HBondEnergy
/// @author John Karanicolas
/// @author Matthew O'Meara

#ifndef INCLUDED_core_scoring_hbonds_HBondDatabase_hh
#define INCLUDED_core_scoring_hbonds_HBondDatabase_hh

// Unit Headers
#include <core/scoring/hbonds/HBondDatabase.fwd.hh>

// Project Headers
#include <core/scoring/hbonds/polynomial.fwd.hh>
#include <core/scoring/hbonds/types.hh>

// Utility Headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/vector1.hh>

// C++
#include <map>

#include <core/scoring/hbonds/FadeInterval.fwd.hh>

#ifdef WIN32
#include <core/scoring/hbonds/polynomial.hh>
#include <core/scoring/hbonds/FadeInterval.hh>
#endif

namespace core {
namespace scoring {
namespace hbonds {


class HBondDatabase : public utility::pointer::ReferenceCount {

private:

	HBondDatabase();
	HBondDatabase(
		//HBondOptionsCOP hbondoptions
		std::string const & hbond_params_database_tag
	);

	HBondDatabase(
		const HBondDatabase& src
	);

public:

	/// @brief only public way to create an HBondDatabase
	static
	HBondDatabaseCOP
	get_database();


	/// @brief only public way to create an HBondDatabase
	static
	HBondDatabaseCOP
	get_database( std::string const & );

	virtual ~HBondDatabase();

	void
	initialize();

	bool
	initialized() const;

	/// @brief read one dimensional polynomial definitions file
	void
	initialize_HBPoly1D();

	/// @brief read table of evaluation types
	void
	initialize_HBEval();

	/// @brief read table of donor bonding strengths
	void
	initialize_don_strength();

	/// @brief read table of acceptor bonding strengths
	void
	initialize_acc_strength();

	/// @brief read table of fade intervals
	void
	initialize_HBFadeInterval();

	/// @brief find polynomial function given name
	FadeIntervalCOP
	HBFadeInterval_from_name( std::string const & name ) const;

	/// @brief find fading function for hbgd_AHdist sort
	FadeIntervalCOP
	AHdist_short_fade_lookup( Size const hb_eval_type ) const;

	/// @brief find fading function for hbgd_AHdist long
	FadeIntervalCOP
	AHdist_long_fade_lookup( Size const hb_eval_type ) const;

	/// @brief find fading function for hbgd_cosBAH
	FadeIntervalCOP
	cosBAH_fade_lookup( Size const hb_eval_type ) const;

	/// @brief find fading function for hbgd_cosAHD
	FadeIntervalCOP
	cosAHD_fade_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial function given name
	Polynomial_1dCOP
	HBPoly1D_from_name( std::string const & name ) const;

	/// @brief find polynomial to hbgd_AHdist dimension
	Polynomial_1dCOP
	AHdist_poly_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial to hbgd_cosBAH dimension when hbgd_AHdist is short
	Polynomial_1dCOP
	cosBAH_short_poly_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial to hbgd_cosBAH dimension when hbgd_AHdist is long
	Polynomial_1dCOP
	cosBAH_long_poly_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial to hbgd_cosAHD dimension when hbgd_AHdist is short
	Polynomial_1dCOP
	cosAHD_short_poly_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial to hbgd_cosAHD dimension when hbgd_AHdist is long
	Polynomial_1dCOP
	cosAHD_long_poly_lookup( Size const hb_eval_type ) const;

	/// @brief find polynomial to hbgd_chi dimension
	Polynomial_1dCOP
	chi_poly_lookup( Size const hb_eval_type ) const;

	/// @brief get the bonding strength of a donor group
	Real
	don_strength( HBDonChemType const don_chem_type ) const;

	/// @brief get the bonding strength of an acceptor group
	Real
	acc_strength( HBAccChemType const ac_chem_type ) const;

	/// @brief find weight type for evaluation type
	HBondWeightType
	weight_type_lookup( Size const hb_eval_type ) const;

	/// @details Signal to use deprecated derivitive calculation in
	///core::scoring::hbonds::hb_energy_deriv_u2(). Once old code has
	///been modified to support the new behavior, remove this
	///option. Since the options are not passe directly to to
	///hb_energy_deriv_u2, access it through the HBondDatabase, rather
	///then messing with the interfaces for the hb_energy_deriv functions.
	inline
	bool
	use_incorrect_deriv() const {
		return false; // APL TEMP return hb_options_->use_incorrect_deriv();
	}

	void
	report_parameter_features_schema_to_db(
		utility::sql_database::sessionOP db_session) const;

private:
	void
	write_hbond_fade_interval_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_polynomial_1d_table_schema(
		utility::sql_database::sessionOP db_session) const;

	void
	write_hbond_evaluation_types_table_schema(
		utility::sql_database::sessionOP db_session) const;

public:
	core::Size
	report_parameter_features(
		utility::sql_database::sessionOP db_session ) const;

private:
	bool initialized_;
	//HBondOptionsCOP hb_options_;
	std::string params_database_tag_; // the tag for the set of files that describe this database

	std::map< const std::string, FadeIntervalCOP > HBFadeInterval_lookup_by_name_;
	utility::vector1< FadeIntervalCOP > HBFadeInterval_lookup_;
	utility::vector1< FadeIntervalCOP > AHdist_short_fade_lookup_;
	utility::vector1< FadeIntervalCOP > AHdist_long_fade_lookup_;
	utility::vector1< FadeIntervalCOP > cosBAH_fade_lookup_;
	utility::vector1< FadeIntervalCOP > cosAHD_fade_lookup_;

	std::map< const std::string, Polynomial_1dCOP > HBPoly1D_lookup_by_name_;
	utility::vector1< Polynomial_1dCOP > HBPoly1D_lookup_;
	utility::vector1< Polynomial_1dCOP > AHdist_poly_lookup_;
	utility::vector1< Polynomial_1dCOP > cosBAH_short_poly_lookup_;
	utility::vector1< Polynomial_1dCOP > cosBAH_long_poly_lookup_;
	utility::vector1< Polynomial_1dCOP > cosAHD_short_poly_lookup_;
	utility::vector1< Polynomial_1dCOP > cosAHD_long_poly_lookup_;
	utility::vector1< Polynomial_1dCOP > chi_poly_lookup_;
	utility::vector1< Real > don_strength_lookup_;
	utility::vector1< Real > acc_strength_lookup_;
	utility::vector1< HBondWeightType > weight_type_lookup_;

	static std::map< const std::string, HBondDatabaseCOP > initialized_databases_;

#ifdef CXX11  // in its current form HBondDatabase is not assignable due to presense of std::map< const std::string, ...> but compiler tries to generate assigment operator anyway
	HBondDatabase & operator= ( const HBondDatabase & ) = delete;
#endif

};

} // hbonds
} // scoring
} // core

#endif
