#include <devel/init.hh>

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
#include <core/chemical/rdkit/util.hh>
#include <core/chemical/rdkit/RDMolToRestype.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/ResidueGraphTypes.hh>

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
#include <numeric/model_quality/rms.hh>

#include <set>
#include <fstream>
#include <string>

OPT_1GRP_KEY(String, align, ligandname)

static basic::Tracer TR( "apps.align_substructure" );

typedef utility::vector1<std::pair<core::Size, core::Size>> SizePairVec;
typedef utility::keys::Key4Tuple< core::Size, core::Size, core::Size, core::Size > DihedralAtomTuple;
typedef std::map< utility::keys::Key4Tuple< core::Size, core::Size, core::Size, core::Size >, core::Size> DihedralAtomTuple2ChiIdxMap;

void
align_pose(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
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
		template_coords.push_back( template_rsd.xyz(it.first) );
		com_template += template_rsd.xyz(it.first);
		ligand_coords.push_back( ligand_rsd.xyz(it.second) );
		com_ligand += ligand_rsd.xyz(it.second);
	}
	com_template /= n_points;
	com_ligand /= n_points;

	numeric::model_quality::findUU( template_coords, ligand_coords, ww, n_points, R, sigma3 );

	numeric::xyzVector<numeric::Real> new_xyz(0.0);
	for ( core::Size i=1; i<=ligand_rsd.natoms(); i++ ) {
		new_xyz = R*(ligand_rsd.xyz(i)-com_ligand)+com_template;
		ligand_pose.set_xyz( core::id::AtomID( i, ligand_idx ), new_xyz );
	}
	ligand_pose.dump_pdb("pose_after_aligned.pdb");

}

void
set_torsion_and_align(core::pose::Pose const& template_pose, core::pose::Pose &ligand_pose,
	std::map<core::Size, core::Size> const& pair_indices_map, core::Size const& template_idx, core::Size const& ligand_idx){

	DihedralAtomTuple chi_atoms;
	std::set<DihedralAtomTuple> chi_atomtuple_set_template, chi_atomtuple_set_ligand;
	DihedralAtomTuple2ChiIdxMap atomtuple_chi_idx_map_template, atomtuple_chi_idx_map_ligand;
	core::conformation::Residue const& template_rsd(template_pose.residue(template_idx));
	core::conformation::Residue const& ligand_rsd(ligand_pose.residue(ligand_idx));
	core::Size atm1, atm2, atm3, atm4;

	for ( core::Size i=1; i<=template_rsd.nchi(); i++ ) {
		atm1 = template_rsd.chi_atoms(i)[1];
		atm2 = template_rsd.chi_atoms(i)[2];
		atm3 = template_rsd.chi_atoms(i)[3];
		atm4 = template_rsd.chi_atoms(i)[4];
		if ( pair_indices_map.find(atm1) == pair_indices_map.end() || pair_indices_map.find(atm2) == pair_indices_map.end() ||
				pair_indices_map.find(atm3) == pair_indices_map.end() || pair_indices_map.find(atm4) == pair_indices_map.end() ) {
			TR << "Cannot find " << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4 << " in template chis" << std::endl;
			continue;
		}

		chi_atoms = DihedralAtomTuple(atm1, atm2, atm3, atm4);
		TR << "Template chi atom ids: " << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4 << std::endl;
		atomtuple_chi_idx_map_template[chi_atoms] = i;
	}

	for ( core::Size i=1; i<=ligand_rsd.nchi(); i++ ) {
		atm1 = ligand_rsd.chi_atoms(i)[1];
		atm2 = ligand_rsd.chi_atoms(i)[2];
		atm3 = ligand_rsd.chi_atoms(i)[3];
		atm4 = ligand_rsd.chi_atoms(i)[4];
		chi_atoms = DihedralAtomTuple(atm1, atm2, atm3, atm4);
		TR << "Ligand chi atom ids: "  << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4 << std::endl;
		atomtuple_chi_idx_map_ligand[chi_atoms] = i;
		chi_atoms = DihedralAtomTuple(atm4, atm3, atm2, atm1);
		atomtuple_chi_idx_map_ligand[chi_atoms] = i;
	}

	for ( auto it:atomtuple_chi_idx_map_template ) {
		core::Real chi_value( template_rsd.chi()[it.second] );
		atm1 = pair_indices_map.at( it.first.key1() );
		atm2 = pair_indices_map.at( it.first.key2() );
		atm3 = pair_indices_map.at( it.first.key3() );
		atm4 = pair_indices_map.at( it.first.key4() );
		chi_atoms = DihedralAtomTuple(atm1, atm2, atm3, atm4);
		if ( atomtuple_chi_idx_map_ligand.find(chi_atoms) == atomtuple_chi_idx_map_ligand.end() ) continue;
		core::Size ichi( atomtuple_chi_idx_map_ligand[chi_atoms] );
		TR << "Setting ligand chi: " << ichi << ", atom ids: "  << atm1 << ", " <<  atm2 << ", " << atm3 << ", " << atm4
			<< ", value from, to: " << ligand_rsd.chi()[ichi] << ", " << chi_value
			<< std::endl;
		ligand_pose.set_chi(ichi, ligand_idx, chi_value);
		ligand_pose.dump_pdb("ligand_pose_after_set_chi.pdb");
	}
	align_pose(template_pose, ligand_pose, pair_indices_map, template_idx, ligand_idx);

}

int
main( int argc, char * argv [] ){
	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		using namespace core::chemical::rdkit;

		NEW_OPT( align::ligandname, "the ligand residue name that will be aligned", "LG1" );

		devel::init(argc, argv);
		core::pose::PoseOP template_pose( new core::pose::Pose() );

		if ( option[ in::file::native ].user() ) {
			core::import_pose::pose_from_file( *template_pose, option[ in::file::native ]().name() , core::import_pose::PDB_file);
		} else {
			utility_exit_with_message("Need to specify in::file::native as template.");
		}

		if ( template_pose->size() !=1 ) {
			utility_exit_with_message("Template has to be single residue.");
		}
		core::chemical::IndexIndexMapping molIndexMapping;
		RestypeToRDMolOptions restype2rdmol_options;
		restype2rdmol_options.neutralize = false;
		restype2rdmol_options.keep_hydro = true;
		restype2rdmol_options.sanitize = false;
		restype2rdmol_options.noImplicitHs = true;
		restype2rdmol_options.skipHs = true;
		restype2rdmol_options.aro2double = true;
		MutableResidueType templateResType(template_pose->residue(1).type());
		RestypeToRDMol template_restype_to_mol(templateResType, restype2rdmol_options); //res, neutralize, keep_hydro, sanitize, noImplicitHs, skipHs


		core::conformation::ResidueOP ligand = core::conformation::get_residue_from_name( option[align::ligandname] );
		MutableResidueType ligandResType(ligand->type());
		RestypeToRDMol ligand_restype_to_mol(ligandResType, restype2rdmol_options);
		core::pose::PoseOP ligpose( new core::pose::Pose() );
		ligpose->append_residue_by_jump(*ligand, 1);

		::RDKit::RWMOL_SPTR template_rdmol( template_restype_to_mol.Mol() );
		::RDKit::RWMOL_SPTR ligand_rdmol( ligand_restype_to_mol.Mol() );

		for ( auto atom:template_rdmol->atoms() ) {
			::RDKit::AtomPDBResidueInfo const* ari( dynamic_cast< ::RDKit::AtomPDBResidueInfo* >( atom->getMonomerInfo() ));
			if ( ari != nullptr ) {
				TR << "Template Atom id and name, serial number : " << atom->getIdx() << ", " << ari->getName() << ", "
					<< ari->getSerialNumber() << std::endl;
			} else {
				TR << "No AtomPDBResidueInfo!" << std::endl;
			}
		}

		for ( auto atom:ligand_rdmol->atoms() ) {
			::RDKit::AtomPDBResidueInfo const* ari( dynamic_cast< ::RDKit::AtomPDBResidueInfo* >( atom->getMonomerInfo() ));
			if ( ari != nullptr ) {
				TR << "Ligand Atom id and name, serial number : " << atom->getIdx() << ", " << ari->getName() << ", "
					<< ari->getSerialNumber() << std::endl;
			} else {
				TR << "No AtomPDBResidueInfo!" << std::endl;
			}
		}


		molIndexMapping = find_mapping( template_rdmol, ligand_rdmol );
		molIndexMapping.show(TR);
		SizePairVec matched_pair_indices;
		std::map<core::Size, core::Size> pair_indices_map;
		for ( core::chemical::IndexIndexMapping::const_iterator iter(molIndexMapping.begin()), iter_end(molIndexMapping.end()); iter != iter_end; ++iter ) {
			VD atom_vd(template_restype_to_mol.index_to_vd()[iter->first]);
			core::Size aidx_template(templateResType.atom_index(atom_vd));
			std::string aname(templateResType.atom_name(atom_vd));
			TR << "Template pose atom index and atomname: " << aidx_template << ", " << aname << std::endl;
			atom_vd = ligand_restype_to_mol.index_to_vd()[iter->second];
			core::Size aidx_ligand( ligandResType.atom_index(atom_vd) );
			aname = ligandResType.atom_name(atom_vd);
			TR << "Ligand pose atom index and atomname: " << aidx_ligand << ", " << aname << std::endl;
			matched_pair_indices.emplace_back(aidx_template, aidx_ligand);
			pair_indices_map[aidx_template] = aidx_ligand;
		}

		TR << "matched_pair_indices size: " << matched_pair_indices.size() << std::endl;
		TR << "set_torsion_and_align" << std::endl;
		core::Size template_idx(1), ligand_idx(1);
		set_torsion_and_align(*template_pose, *ligpose, pair_indices_map, template_idx, ligand_idx );

		// ::RDKit::MolToPDBFile(*template_restype_to_mol.Mol(), "template_mol.pdb", 0);

	} catch (utility::excn::Exception const & e) {
		e.display();
		return -1;
	}
	return 0;
}

