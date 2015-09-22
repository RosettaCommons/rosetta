// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rotamers/SingleBasicRotamerLibrary.cc
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

// Unit headers
#include <core/pack/rotamers/SingleBasicRotamerLibrary.hh>
#include <core/pack/rotamers/SingleBasicRotamerLibraryCreator.hh>
#include <core/chemical/rotamers/BasicRotamerLibrarySpecification.hh>

// Package headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TenANeighborGraph.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>
#include <core/pack/dunbrack/ChiSet.fwd.hh>
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/ResidueFactory.hh>

// Utility headers
#include <basic/Tracer.hh>

// C++ headers

namespace core {
namespace pack {
namespace rotamers {

static THREAD_LOCAL basic::Tracer TR("core.pack.rotamers.SingleBasicRotamerLibrary");

///////////////////  Creator Functions ///////////////////////

core::pack::rotamers::SingleResidueRotamerLibraryCOP
SingleBasicRotamerLibraryCreator::create( core::chemical::ResidueType const & ) const {
	return core::pack::rotamers::SingleResidueRotamerLibraryCOP( new SingleBasicRotamerLibrary );
}

std::string
SingleBasicRotamerLibraryCreator::keyname() const {
	return core::chemical::rotamers::BasicRotamerLibrarySpecification::library_name();
}

///////////////////  RotamerLibrary Functions ///////////////////////

SingleBasicRotamerLibrary::SingleBasicRotamerLibrary()
{}

SingleBasicRotamerLibrary::~SingleBasicRotamerLibrary()
{}

using namespace core::pack::dunbrack;

/// @details Single rotamer - no derivative
Real
SingleBasicRotamerLibrary::rotamer_energy_deriv(
	conformation::Residue const & /*rsd*/,
	RotamerLibraryScratchSpace & /*scratch*/
) const
{
	return 0;
}


/// @details Single rotamer - energy is zero
Real
SingleBasicRotamerLibrary::rotamer_energy(
	conformation::Residue const & /*rsd*/,
	RotamerLibraryScratchSpace & /*scratch*/
) const
{
	return 0;
}

/// @details Single rotamer - best is zero
Real
SingleBasicRotamerLibrary::best_rotamer_energy(
	conformation::Residue const & /*rsd*/,
	bool /*curr_rotamer_only*/,
	RotamerLibraryScratchSpace & /*scratch*/
) const
{
	return 0;
}

// @details Not currently implemented -- chokes.
void
SingleBasicRotamerLibrary::assign_random_rotamer_with_bias(
	conformation::Residue const &,// rsd,
	pose::Pose const & /*pose*/,
	RotamerLibraryScratchSpace &,// scratch,
	numeric::random::RandomGenerator &,// RG,
	ChiVector &,// new_chi_angles,
	bool //perturb_from_rotamer_center
) const {
	utility_exit_with_message("SingleBasicRotamerLibrary does not yet implement assign_random_rotamer_with_bias.");
}

void
SingleBasicRotamerLibrary::fill_rotamer_vector(
	pose::Pose const &pose,
	scoring::ScoreFunction const &,
	pack::task::PackerTask const & task,
	graph::GraphCOP,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & /*extra_chi_steps*/,
	bool buried,
	RotamerVector & rotamers //utility::vector1< conformation::ResidueOP >
) const
{
	// Add proton chi rotamers even if there's no rotamer library.
	core::Size nrotadded(0);
	if ( concrete_residue->n_proton_chi() != 0 ) {
		utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
		proton_chi_chisets.push_back( pack::dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );
		for ( Size ii = 1; ii <= concrete_residue->n_proton_chi(); ++ii ) {
			pack::dunbrack::expand_proton_chi(
				task.residue_task( existing_residue.seqpos() ).extrachi_sample_level(
				buried,
				concrete_residue->proton_chi_2_chi( ii ),
				*concrete_residue ),
				concrete_residue,
				ii, proton_chi_chisets);
		}

		rotamers.reserve( rotamers.size() + proton_chi_chisets.size() );
		for ( Size ii = 1; ii <= proton_chi_chisets.size(); ++ii ) {
			conformation::ResidueOP rotamer;
			if ( existing_residue.type().name() == concrete_residue->name() ) {
				rotamer = existing_residue.clone(); //Clone the existing residue if the types match.  (Copies chi angles).
			} else { //Otherwise, build a new residue and align to the existing residue if the types do not match.
				rotamer = core::conformation::ResidueFactory::create_residue( *concrete_residue, existing_residue, pose.conformation(), false ) ;
				//Copy the side-chain chi values if possible:
				for ( core::Size ichi=1, ichimax=rotamer->nchi(); ichi<=ichimax; ++ichi ) {
					if ( ichi > existing_residue.nchi() ) {
						rotamer->set_chi(ichi, 0.0); //Initialize to zero so that values won't be random.
						continue;
					}
					core::conformation::Residue::AtomIndices const & chi_atoms( rotamer->type().chi_atoms( ichi ) );
					core::conformation::Residue::AtomIndices const & existing_chi_atoms( existing_residue.type().chi_atoms( ichi ) );

					if ( existing_residue.atom_name(existing_chi_atoms[1]) == rotamer->atom_name(chi_atoms[1]) &&
							existing_residue.atom_name(existing_chi_atoms[2]) == rotamer->atom_name(chi_atoms[2]) &&
							existing_residue.atom_name(existing_chi_atoms[3]) == rotamer->atom_name(chi_atoms[3]) &&
							existing_residue.atom_name(existing_chi_atoms[4]) == rotamer->atom_name(chi_atoms[4])
							) {
						rotamer->set_chi(ichi, existing_residue.chi(ichi));
					} else {
						rotamer->set_chi(ichi, 0.0); //Initialize to zero in the absence of anything better that we can do.
					}
				}
			}
			for ( Size jj = 1; jj <= concrete_residue->n_proton_chi(); ++jj ) {
				Size jj_protchi = concrete_residue->proton_chi_2_chi( jj );
				rotamer->set_chi(jj_protchi, proton_chi_chisets[ ii ]->chi[ jj_protchi ] );
			}
			rotamers.push_back( rotamer );
			++nrotadded;
		}
	}
	TR.Debug << "Building rotamers for " << concrete_residue->name() << ": Total number suggested: " <<
		nrotadded << std::endl;
}

/// @details Not implemented -- will cause program termination.
/// Is this only used by coarse representations?
void
SingleBasicRotamerLibrary::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	utility_exit_with_message("SingleBasicRotamerLibrary does not currently implement write_to_file()");
}

} // namespace rotamers
} // namespace scoring
} // namespace core
