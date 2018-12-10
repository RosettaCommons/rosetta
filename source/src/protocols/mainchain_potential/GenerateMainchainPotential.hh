// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/mainchain_potential/GenerateMainchainPotential.hh
/// @brief Headers for a generator for mainchain potentials.  Inputs are a noncanonical residue type with an already-generated
/// sidechain potential; outputs are a potential file suitable for use by the RamaPrePro scoreterm.
/// @author Vikram K. Mulligan (vmulligan@flatironinstitute.org)


#ifndef INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_hh
#define INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_hh

// Project headers
#include <protocols/mainchain_potential/GenerateMainchainPotential.fwd.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialOptions.fwd.hh>
#include <protocols/mainchain_potential/GenerateMainchainPotentialTests.fwd.hh>

// Core headers
#include <core/pose/Pose.fwd.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/mainchain_potential/MainchainScoreTable.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/ScoreType.hh>

// Utility headers
#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

// C++ headers
#include <map>

namespace protocols {
namespace mainchain_potential {

/// @brief A generator for mainchain potentials.  Inputs are a noncanonical residue type with an already-generated sidechain potential;
/// outputs are a potential file suitable for use by the RamaPrePro scoreterm.
class GenerateMainchainPotential : public utility::pointer::ReferenceCount {

	friend class ::GenerateMainchainPotentialTests; //Needed to allow the unit tests for this class access internal data members.

public:

	/// @brief Default constructor.
	/// @details Optionally takes a GenerageMainchainPotentialOptions const-owning pointer.  If null, this object initializes itself
	/// to default values.
	GenerateMainchainPotential( GenerateMainchainPotentialOptionsCOP options = nullptr );

	/// @brief Copy constructor.
	GenerateMainchainPotential(GenerateMainchainPotential const & src);

	/// @brief Destructor.
	virtual ~GenerateMainchainPotential();

	/// @brief Clone function: make a copy of this object and return an owning pointer to the copy.
	GenerateMainchainPotentialOP clone() const;

public:

	/// @brief Entry point into protocol execution.
	void run();

	/// @brief Write the last generated mainchain potential to disk.  The output filename is given in the
	/// options_ object.
	/// @details If we're writing individual scoretables for individual scoreterms, that happens here, too.
	void write_last_generated_to_disk() const;

private: //Functions

	/// @brief Get the scorefunction that we'll be using.
	core::scoring::ScoreFunctionOP generate_sfxn() const;

	/// @brief Generate the one-residue pose, based on options.
	core::pose::PoseOP generate_pose() const;

	/// @brief Can this residue type take protein terminal patches?
	bool protein_patches_can_apply( core::chemical::ResidueTypeCOP restype ) const;

	/// @brief Can this residue type take specified terminal patches?
	bool patches_can_apply( core::chemical::ResidueTypeCOP restype, utility::vector1< core::chemical::VariantType > const & vartypes ) const;

	/// @brief Given a one-residue pose, apply acetylated N-terminus and aminomethylated C-terminus patches.
	/// @details Applies aminomethylated or aminodimethylated C-terminus patches depending on settings in options.
	void apply_protein_patches( core::pose::Pose & pose ) const;

	/// @brief Given a pose, cycle through mainchain dihedrals and generate the mainchain potential.
	void generate_mainchain_potential( core::pose::PoseCOP pose, core::scoring::ScoreFunctionOP sfxn, core::chemical::mainchain_potential::MainchainScoreTableOP newtable, std::map< core::scoring::ScoreType, core::chemical::mainchain_potential::MainchainScoreTableOP >  & last_generated_scoretables_by_scoreterm ) const;

private: //Data

	/// @brief This protocol's options.
	GenerateMainchainPotentialOptionsOP options_;

	/// @brief The last mainchain potential that this object generated.  Nullptr if nothing has been generated.
	core::chemical::mainchain_potential::MainchainScoreTableCOP last_generated_scoretable_;

	/// @brief Mainchain potentials for each scoreterm, separately.
	/// @details Note that these will be kept unnormalized.  Moreover, each scoreterm's weight is treated as though it were 1.0.
	std::map< core::scoring::ScoreType, core::chemical::mainchain_potential::MainchainScoreTableOP > last_generated_scoretables_by_scoreterm_;

};


} //protocols
} //mainchain_potential



#endif //INCLUDED_protocols_mainchain_potential_GenerateMainchainPotential_hh





