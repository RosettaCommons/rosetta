// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

	// This loops over a rotamer set and dumps each rotamer as it's own separate pdb
	for( core::Size z(1); z <= rotset->num_rotamers(); ++z ) {
			core::pose::Pose posetest2;
			posetest2.append_residue_by_jump( *(rotset->nonconst_rotamer(z)), 1);
			std::stringstream pose_test_name;
			pose_test_name << "ROTAMERmu_f" << z << ".pdb";
			core::io::pdb::dump_pdb(posetest2, pose_test_name.str() );
		}

// Makes a ResidueOP of whatever type is specified
core::conformation::ResidueOP baseres = core::conformation::ResidueFactory::create_residue( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD )->name_map( protocols::dna::dna_full_name3(motifcop->restype_name2()) ) );

// Displays the mainchain torsions for this residue
std::cout<< "Residue name: " << (new_rots_[i])->name() << std::endl;
utility::vector1<Real> mainchains((new_rots_[i])->mainchain_torsions() );
for( Size v(1); v <= mainchains.size(); ++v) {
	std::cout << "Mainchain torsion #" << v << " is " << mainchains[v] << std::endl;
}
