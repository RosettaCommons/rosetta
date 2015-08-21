// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/scoring/rna/RNA_EnergyMethodOptions.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_core_scoring_rna_RNA_EnergyMethodOptions_HH
#define INCLUDED_core_scoring_rna_RNA_EnergyMethodOptions_HH

#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/rna/RNA_EnergyMethodOptions.fwd.hh>
#include <core/types.hh>
#include <string>

namespace core {
namespace scoring {
namespace rna {

class RNA_EnergyMethodOptions: public utility::pointer::ReferenceCount {

public:

	//constructor
	RNA_EnergyMethodOptions();

	//destructor
	~RNA_EnergyMethodOptions();

	void initialize_from_options();

public:

	/// @brief Parameter for adjusting syn-G potential level in RNA_TorsionPotential.
	core::Real syn_G_potential_bonus() const { return syn_G_potential_bonus_; }
	void syn_G_potential_bonus( core::Real setting ){ syn_G_potential_bonus_ = setting; }

	/// @brief Parameter for changing torsion_potential directory name.
	std::string torsion_potential() const { return torsion_potential_; }
	void torsion_potential( std::string setting ){ torsion_potential_ = setting; }

	/// @brief Parameter for changing suiteness_bonus directory name.
	std::string suiteness_bonus() const { return suiteness_bonus_; }
	void suiteness_bonus( std::string setting ){ suiteness_bonus_ = setting; }

	friend
	bool
	operator==( RNA_EnergyMethodOptions const & a, RNA_EnergyMethodOptions const & b );

	friend
	std::ostream &
	operator<< ( std::ostream & out, const RNA_EnergyMethodOptions & options );

	void
	show( std::ostream & out ) const;

private:

	core::Real syn_G_potential_bonus_;
	std::string torsion_potential_;
	std::string suiteness_bonus_;

};

} //rna
} //scoring
} //core

#endif
