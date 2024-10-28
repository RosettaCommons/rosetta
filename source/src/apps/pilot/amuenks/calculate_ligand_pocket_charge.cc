#include <devel/init.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/import_pose.hh>

#include <core/chemical/Atom.hh>
#include <core/chemical/AtomType.hh> 
#include <core/conformation/Residue.hh>
#include <core/conformation/util.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/random/random.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/string_util.hh>

#include <basic/Tracer.hh>
#include <cmath>
#include <limits>

#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/excn/Exceptions.hh>

//for debugging coordinates
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/ResidueFactory.hh>

static basic::Tracer TR("amuenks.charge_in_binding_pocket");

using namespace basic::options;
using namespace basic::options::OptionKeys;

int
main( int argc, char * argv[] ){
    try {
        
        devel::init(argc, argv);

        core::pose::Pose p;
        utility::vector1< utility::file::FileName > fname;
        fname = basic::options::option[in::file::s]();
        core::import_pose::pose_from_file(p, fname[1], core::import_pose::PDB_file);
        core::Size nres = p.size();

	TR << "Ligand is " << p.residue(nres).name() << std::endl;

	core::conformation::Residue const &lig( p.residue( nres ) );
	core::Real distance_to_atom = 0.0;

	numeric::xyzVector< core::Real > lig_com(0,0,0);
	for ( core::Size j=1; j<=lig.natoms(); ++j ) {
		lig_com += lig.xyz(j);
	}
	core::Size nAtms = lig.natoms();

	lig_com /= nAtms;
	TR << lig_com << std::endl;

	core::Size ncharged_residues = 0;
	core::Size ncharged_atoms = 0;
	core::Real absolute_charge_total = 0.0;
	core::Size nearby_atoms = 0;
	for ( core::Size ires = 1; ires < nres; ++ires ) {
		//bool is_charge = false;
		core::conformation::Residue const &resi( p.residue( ires ) );
		core::Real res_distance = 20.0;
		if ( resi.is_protein() ) {
			res_distance = resi.xyz(resi.atom_index("CA")).distance(lig_com);
		}
		if ( resi.is_ligand() ) {
			res_distance = resi.xyz(1).distance(lig_com);
		}
		
		if ( res_distance < 15.0 ) {
			for ( core::Size iatm = resi.first_sidechain_atom(); iatm <= resi.natoms(); ++iatm ) {
				if ( resi.atom_type( iatm ).element() == "N" || resi.atom_type( iatm ).element() == "O" || resi.is_metal() ) {
					bool already_seen = false;
					//TR << resi.name3() << "," << ires << "," << resi.atom_type(iatm).name() << "," << resi.atomic_charge(iatm) << std::endl;
					for ( core::Size latm = 1; latm <= lig.nheavyatoms(); ++latm ) {
						core::Real atom_distance = resi.xyz(iatm).distance(lig.xyz(latm));
						if ( atom_distance <= 5.0 ) {
							//TR << resi.name3() << "," << ires << std::endl;
							if ( !already_seen ) {
								already_seen = true;
								nearby_atoms += 1;
								if ( resi.type().formal_charge(iatm) != 0 ) {//( resi.is_charged() || resi.is_metal() ) && !is_charge ) { 
									ncharged_residues++; 
									//is_charge = true;
									TR << resi.name3() << "," << ires << "," << resi.type().formal_charge(iatm) << "," << resi.atomic_charge(iatm) << std::endl;
								}
								core::Real abs_charge = abs(resi.atomic_charge(iatm));
								absolute_charge_total += abs_charge;
								if ( abs_charge > 0.75 ) { ncharged_atoms++; }
							}
						}
					}
				}		

			}
		}
	}

	TR << "Number of charged residues: " << ncharged_residues << std::endl;
	TR << "Number of charged atoms: " << ncharged_atoms << std::endl;
	TR << "Average absolute atom charge: " << absolute_charge_total/nearby_atoms << std::endl;

    } catch ( utility::excn::Exception const & e ){
        e.display();
        return -1;
    }
    return 0;
}
