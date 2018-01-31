// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/NMerRefEnergy.hh
/// @brief  Reference energy method declaration
/// @author Chris King


#ifndef INCLUDED_core_scoring_methods_NMerRefEnergy_hh
#define INCLUDED_core_scoring_methods_NMerRefEnergy_hh

// Unit headers
#include <core/scoring/methods/NMerRefEnergy.fwd.hh>

// Package headers
#include <core/scoring/methods/ContextIndependentOneBodyEnergy.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

// Project headers
#include <core/pose/Pose.fwd.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/DOF_ID.fwd.hh>

#include <utility/vector1.hh>
#include <basic/Tracer.hh>
#include <map>


// Utility headers


namespace core {
namespace scoring {
namespace methods {

class NMerRefEnergy : public ContextIndependentOneBodyEnergy
{
public:
	typedef ContextIndependentOneBodyEnergy parent;

public:


	NMerRefEnergy();


	NMerRefEnergy( utility::vector1< std::map< std::string, core::Real > > const & nmer_ref_energies_in );


	NMerRefEnergy(
		core::Size const nmer_length,
		utility::vector1< std::string > const & fname_vec
	);


	virtual ~NMerRefEnergy();

	virtual
	EnergyMethodOP
	clone() const;

	/////////////////////////////////////////////////////////////////////////////
	// methods for ContextIndependentOneBodyEnergies
	/////////////////////////////////////////////////////////////////////////////

	/// @brief return the reference energy for a given reference list at a nmer frame sequence position
	void
	get_residue_energy_by_table(
		core::pose::Pose const &,
		core::Size const &,
		core::Real &,
		utility::vector1< core::Real > &
	) const;

	virtual
	void
	residue_energy(
		conformation::Residue const & rsd,
		pose::Pose const & pose,
		EnergyMap & emap
	) const;

	bool
	minimize_in_whole_structure_context( pose::Pose const & ) const { return false; }


	virtual
	Real
	eval_dof_derivative(
		id::DOF_ID const & dof_id,
		id::TorsionID const & tor_id,
		pose::Pose const & pose,
		ScoreFunction const & sfxn,
		EnergyMap const & weights
	) const;

	/// @brief DunbrackEnergy is context independent; indicates that no
	/// context graphs are required
	virtual
	void indicate_required_context_graphs( utility::vector1< bool > & ) const;

	/// @brief read reference energy lists from a vector of filenames
	void read_nmer_fname_vector( utility::vector1< std::string > const & );

	//methods to init from outside (for filter construction)
	void nmer_length( Size const nmer_length );
	/// @brief init all options from option flags
	void initialize_from_options();
	/// @brief return number of tables currently loaded in memory
	core::Size n_tables() const;

private:
	utility::vector1< std::map< std::string, core::Real > > nmer_ref_energies_;
	core::Size nmer_length_;
	core::Size nmer_cterm_;

	/// @brief read reference energies from a file containing seq strings and values
	void read_nmer_table( std::string const & ref_fname );
	/// @brief read reference energy lists from a file containing a list of filenames
	void read_nmer_table_list( std::string const & ref_list_fname );
	/// @brief read reference tables from option flags
	void read_nmer_tables_from_options();
	virtual
	core::Size version() const;
};

} // methods
} // scoring
} // core


#endif
