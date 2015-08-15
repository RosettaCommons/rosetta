// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   rosetta/benchmark/pdb_io.bench.cc
///
/// @brief  Performance benchmark for PDB input and output
/// @author Gordon Lemmon

#include <apps/benchmark/performance/performance_benchmark.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueConnection.hh>

#include <fstream>

#include <utility/vector1.hh>


using namespace core;

bool is_deletable(chemical::ResidueType const& rsd, Size atom_id){
	chemical::AtomIndices const & mainchain_atoms = rsd.mainchain_atoms(); // Can't delete mainchain atoms.
	chemical::AtomIndices::const_iterator const found = find( mainchain_atoms.begin(), mainchain_atoms.end(), atom_id);
	if( found != mainchain_atoms.end() ) return false;

	for(Size i=1; i<= rsd.natoms(); ++i){
		if( rsd.atom_base(i) == atom_id) return false;
		for ( Size j=1; j<= 3; ++j ) {
			chemical::ICoorAtomID const & stub_atom( rsd.icoor(i).stub_atom( j ) );
			Size const atomno = stub_atom.atomno();
			if( atomno == atom_id ) return false;
		}
	}
	//Can't be part of a ResidueConnection's set of stub atoms...
	for( Size i=1; i<= rsd.n_residue_connections(); ++i){
		chemical::AtomICoor const& icoor = rsd.residue_connection(i).icoor();
		for ( Size j = 1; j <= 3; ++j ) {
			if( icoor.stub_atom(j).atomno() == atom_id) return false;
		}
	}
	// Can't be part of a chi...
	for(Size i = 1; i <= rsd.nchi(); ++i){
		chemical::AtomIndices const & chi_atoms = rsd.chi_atoms(i);
		chemical::AtomIndices::const_iterator const found2 = find(chi_atoms.begin(), chi_atoms.end(), atom_id);
		if( found2 != chi_atoms.end()) return false;
	}

	return true;
}

class ResidueTypeBenchmark : public PerformanceBenchmark
{
public:
	chemical::ChemicalManager * cm_;
	chemical::ResidueTypeSetCAP residue_types_;

	ResidueTypeBenchmark(std::string name) : PerformanceBenchmark(name) {};

	virtual void setUp() {
		cm_ = chemical::ChemicalManager::get_instance();
		residue_types_ = cm_->residue_type_set("fa_standard");
	}

	virtual void run(core::Real scaleFactor) {
		chemical::ResidueTypeCOPs::const_iterator const begin = residue_types_->residue_types_DO_NOT_USE().begin();
		chemical::ResidueTypeCOPs::const_iterator const end = residue_types_->residue_types_DO_NOT_USE().end();
		core::Size local_scale_factor = 2 * scaleFactor;
		if( local_scale_factor == 0 ) { local_scale_factor = 1; } // Do at least one repetition, regardless of scaling factor
		for(core::Size i=0; i < local_scale_factor; i++){
			for(chemical::ResidueTypeCOPs::const_iterator iter=begin; iter != end; ++iter) {
				chemical::ResidueTypeCOP res_type = *iter;
				chemical::ResidueTypeOP copy = res_type->clone(); // Tests the copy constructor

				if ( copy->is_RNA() ) continue; // RNA ResidueType can't run finalize after certain atoms are deleted

				for(Size i=copy->natoms(); i >= 1; --i){ // start at the end...
				//while(copy->natoms() > 4){ // residue must have 3 atoms
					//Size const last_index = copy->natoms();
					///TODO fix this so atoms restubbify.
					if( is_deletable(*copy, i)){
						copy->delete_atom( i );
						//if( copy->nchi() == 0 )
						copy->finalize(); // finalize doesn't work if an atom is part of a chi definition
					}
				}
				for(Size i=copy->natoms(); i >= 1; --i){
					copy->delete_atom( i ); // Delete non-deletable atoms, but do not finalize after or else!
				}
			}
		}
	};

	virtual void tearDown() {};

	std::string pdb_string_;
};
