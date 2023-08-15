// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ligand_docking/GALigandDock/MCSAligner.cc
///
/// @brief  Use RDKit maximum common substructure to align a ligand to a reference ligand
/// @author Guangfeng Zhou and Frank DiMaio

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>
#include <core/chemical/AtomRefMapping.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>
#include <core/chemical/MutableResidueType.fwd.hh>
#include <core/chemical/ResidueGraphTypes.hh>

//rdkit
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>


#include <rdkit/GraphMol/ROMol.h>
#include <rdkit/GraphMol/FileParsers/MolWriters.h>
#include <rdkit/GraphMol/FileParsers/FileParsers.h>
#include <rdkit/GraphMol/MonomerInfo.h>

#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/id/AtomID.hh>

#include <core/types.hh>


#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/model_quality/rms.hh>

#include <set>
#include <fstream>
#include <string>
#include <ctime>

#include <protocols/ligand_docking/GALigandDock/LigandConformer.fwd.hh>
#include <protocols/ligand_docking/GALigandDock/GridScorer.hh>
#include <protocols/ligand_docking/GALigandDock/MCSAligner.hh>
#include <protocols/ligand_docking/GALigandDock/util.hh>
#include <core/types.hh>

#include <core/pose/Pose.fwd.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <core/id/AtomID.hh> // AUTO IWYU For AtomID


namespace protocols {
namespace ligand_docking {
namespace ga_ligand_dock {

static basic::Tracer TR( "protocols.ligand_docking.GALigandDock.MCSAligner" );

MCSAligner::MCSAligner(core::pose::Pose const&reference_pose, core::Size reference_ligres_idx, MCSAlignerOptions & options):
	reference_pose_(reference_pose),
	reference_ligres_idx_(reference_ligres_idx)
{
	set_options(options);
}



void
MCSAligner::align_pose(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
	std::map<core::Size, core::Size> const& pair_indices_map, core::Size const& template_idx, core::Size const& ligand_idx){
	utility::vector1< numeric::xyzVector<numeric::Real> > template_coords, ligand_coords;
	core::conformation::Residue const& template_rsd( template_pose.residue(template_idx) );
	core::conformation::Residue const& ligand_rsd( ligand_pose.residue(ligand_idx) );
	core::Size n_points( pair_indices_map.size() );
	utility::vector1< numeric::Real > ww( n_points, 1.0 );
	numeric::xyzMatrix< numeric::Real > R( 0.0 );
	numeric::Real sigma3;
	numeric::xyzVector< core::Real > com_template(0,0,0), com_ligand(0,0,0);

	for ( auto it:pair_indices_map ) {
		template_coords.push_back( template_rsd.xyz(it.second) );
		com_template += template_rsd.xyz(it.second);
		ligand_coords.push_back( ligand_rsd.xyz(it.first) );
		com_ligand += ligand_rsd.xyz(it.first);
	}
	com_template /= n_points;
	com_ligand /= n_points;

	numeric::model_quality::findUU( template_coords, ligand_coords, ww, n_points, R, sigma3 );

	numeric::xyzVector<numeric::Real> new_xyz(0.0);
	for ( core::Size i=1; i<=ligand_rsd.natoms(); i++ ) {
		new_xyz = R*(ligand_rsd.xyz(i)-com_ligand)+com_template;
		ligand_pose.set_xyz( core::id::AtomID( i, ligand_idx ), new_xyz );
	}
	if ( TR.Debug.visible() ) {
		ligand_pose.dump_pdb("pose_after_aligned.pdb");
	}

}

void
MCSAligner::set_torsion_and_align(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
	std::map<core::Size, core::Size> const& pair_indices_map, core::Size const& template_idx, core::Size const& ligand_idx){

	DihedralAtomTuple chi_atoms;
	std::set<DihedralAtomTuple> chi_atomtuple_set_template, chi_atomtuple_set_ligand;
	DihedralAtomTuple2ChiIdxMap atomtuple_chi_idx_map_template, atomtuple_chi_idx_map_ligand;
	core::conformation::Residue const& template_rsd(template_pose.residue(template_idx));
	core::conformation::Residue const& ligand_rsd(ligand_pose.residue(ligand_idx));
	core::Size atm1, atm2, atm3, atm4;
	torsion_in_align_.resize( ligand_rsd.nchi(), false );

	for ( core::Size i=1; i<=ligand_rsd.nchi(); i++ ) {
		atm1 = ligand_rsd.chi_atoms(i)[1];
		atm2 = ligand_rsd.chi_atoms(i)[2];
		atm3 = ligand_rsd.chi_atoms(i)[3];
		atm4 = ligand_rsd.chi_atoms(i)[4];
		if ( pair_indices_map.find(atm1) == pair_indices_map.end() || pair_indices_map.find(atm2) == pair_indices_map.end() ||
				pair_indices_map.find(atm3) == pair_indices_map.end() || pair_indices_map.find(atm4) == pair_indices_map.end() ) {
			if ( TR.Debug.visible() ) {
				TR.Debug << "Cannot find " << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4 << " in ligand chis" << std::endl;
			}
			continue;
		}
		chi_atoms = DihedralAtomTuple(atm1, atm2, atm3, atm4);
		if ( TR.Debug.visible() ) {
			TR.Debug << "Ligand chi atom ids: "  << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4 << std::endl;
		}

		atomtuple_chi_idx_map_ligand[chi_atoms] = i;
	}

	for ( auto it:atomtuple_chi_idx_map_ligand ) {
		atm1 = pair_indices_map.at( it.first.key1() );
		atm2 = pair_indices_map.at( it.first.key2() );
		atm3 = pair_indices_map.at( it.first.key3() );
		atm4 = pair_indices_map.at( it.first.key4() );
		core::Vector const &xyz1( template_rsd.xyz( atm1 ) );
		core::Vector const &xyz2( template_rsd.xyz( atm2 ) );
		core::Vector const &xyz3( template_rsd.xyz( atm3 ) );
		core::Vector const &xyz4( template_rsd.xyz( atm4 ) );

		core::Real angle = numeric::dihedral_radians( xyz1, xyz2, xyz3, xyz4 )*180.0/3.14159216;

		core::Size ichi( atomtuple_chi_idx_map_ligand[it.first] );
		if ( TR.Debug.visible() ) {
			TR.Debug << "Setting ligand chi: " << ichi << ", atom ids: "  << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4
				<< ", value from, to: " << ligand_rsd.chi()[ichi] << ", " << angle
				<< std::endl;
		}

		ligand_pose.set_chi(ichi, ligand_idx, angle);
		torsion_in_align_[ichi] = true;
	}
	if ( TR.Debug.visible() ) ligand_pose.dump_pdb("ligand_pose_after_set_chi.pdb");

	align_pose(template_pose, ligand_pose, pair_indices_map, template_idx, ligand_idx);

}

void
MCSAligner::apply( LigandConformer & lig ){
	using namespace core::chemical;
	using namespace core::chemical::rdkit;

	auto t0 = std::chrono::steady_clock::now();

	RestypeToRDMolOptions &restype2rdmol_options(options_.restype_to_rdmol_options);

	if ( TR.Debug.visible() ) {
		TR.Debug << "reference_pose_ size: " << reference_pose_.size() << " , reference_ligres_idx_ : " << reference_ligres_idx_ <<std::endl;
	}
	MutableResidueType reference_restype(reference_pose_.residue(reference_ligres_idx_).type());
	RestypeToRDMol reference_to_rdmol(reference_restype, restype2rdmol_options);

	if ( lig.ligand_ids().size() != 1 ) {
		utility_exit_with_message("Currently substructure alignment only works for sigle residue ligand.");
	}
	if ( TR.Debug.visible() ) {
		TR.Debug << "liggene ref_pose size: " << lig.get_ref_pose()->size() << " ligand_ids 1st: " << lig.ligand_ids()[1] << std::endl;
	}
	core::conformation::Residue const &ligand = lig.ligand_residue(1);
	MutableResidueType ligand_restype(ligand.type());
	RestypeToRDMol ligand_restype_to_mol(ligand_restype, restype2rdmol_options);
	core::pose::PoseOP ligpose( new core::pose::Pose() );
	lig.to_pose(ligpose);

	::RDKit::RWMOL_SPTR template_rdmol( reference_to_rdmol.Mol() );
	::RDKit::RWMOL_SPTR ligand_rdmol( ligand_restype_to_mol.Mol() );

	for ( auto atom:template_rdmol->atoms() ) {
		::RDKit::AtomPDBResidueInfo const* ari( dynamic_cast< ::RDKit::AtomPDBResidueInfo* >( atom->getMonomerInfo() ));
		if ( TR.Debug.visible() ) {
			if ( ari != nullptr ) {
				TR.Debug << "Template Atom id and name, serial number : " << atom->getIdx() << ", " << ari->getName() << ", "
					<< ari->getSerialNumber() << std::endl;
			} else {
				TR.Debug << "No AtomPDBResidueInfo!" << std::endl;
			}
		}
	}

	for ( auto atom:ligand_rdmol->atoms() ) {
		::RDKit::AtomPDBResidueInfo const* ari( dynamic_cast< ::RDKit::AtomPDBResidueInfo* >( atom->getMonomerInfo() ));
		if ( TR.Debug.visible() ) {
			if ( ari != nullptr ) {
				TR.Debug << "Ligand Atom id and name, serial number : " << atom->getIdx() << ", " << ari->getName() << ", "
					<< ari->getSerialNumber() << std::endl;
			} else {
				TR.Debug << "No AtomPDBResidueInfo!" << std::endl;
			}
		}
	}


	IndexIndexMapping mol_index_mapping;
	mol_index_mapping = find_mapping( template_rdmol, ligand_rdmol );
	mol_index_mapping.show(TR);
	SizePairVec matched_pair_indices;
	std::map<core::Size, core::Size> pair_indices_map;
	for ( IndexIndexMapping::const_iterator iter(mol_index_mapping.begin()), iter_end(mol_index_mapping.end()); iter != iter_end; ++iter ) {
		VD atom_vd(reference_to_rdmol.index_to_vd()[iter->first]);
		core::Size aidx_template(reference_restype.atom_index(atom_vd));
		std::string aname(reference_restype.atom_name(atom_vd));
		if ( TR.Debug.visible() ) {
			TR.Debug << "Template pose atom index and atomname: " << aidx_template << ", " << aname << std::endl;
		}
		atom_vd = ligand_restype_to_mol.index_to_vd()[iter->second];
		core::Size aidx_ligand( ligand_restype.atom_index(atom_vd) );
		aname = ligand_restype.atom_name(atom_vd);
		if ( TR.Debug.visible() ) {
			TR.Debug << "Ligand pose atom index and atomname: " << aidx_ligand << ", " << aname << std::endl;
		}
		matched_pair_indices.emplace_back(aidx_template, aidx_ligand);
		pair_indices_map[aidx_ligand] = aidx_template;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "matched_pair_indices size: " << matched_pair_indices.size() << std::endl;
		TR.Debug << "set_torsion_and_align" << std::endl;
	}

	set_torsion_and_align(reference_pose_, *ligpose, pair_indices_map, reference_ligres_idx_, lig.ligand_ids()[1] );

	auto t1= std::chrono::steady_clock::now();
	TR << "Time for set_torsion_and_align: " << std::chrono::duration<double>(t1-t0).count() << " seconds." << std::endl;

	if ( options_.perturb_rb ) {
		perturb_ligand_rb( *ligpose, lig.ligand_ids(), options_.perturb_rb_translation, options_.perturb_rb_rotation );
	}
	if ( options_.perturb_torsion ) {
		perturb_ligand_torsions( *ligpose, lig.ligand_ids(), torsion_in_align_, options_.perturb_torsion_rotation );
	}
	lig.update_conf(ligpose);
	auto t2= std::chrono::steady_clock::now();
	TR << "Time for perturb_ligand_rb and perturb_ligand_torsions: " << std::chrono::duration<double>(t2-t1).count() << " seconds." << std::endl;
}


}
}
}


