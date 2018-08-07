#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ElementSet.hh>
#include <core/chemical/MMAtomTypeSet.hh>
#include <core/chemical/orbitals/OrbitalTypeSet.hh>

#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/residue_io.hh>

#include <core/chemical/ResidueTypeKinWriter.hh>

#include <core/types.hh>

#include <devel/init.hh>

#include <utility/options/FileVectorOption.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <basic/Tracer.hh>

#include <fstream>
#include <string>

using namespace basic::options;

static basic::Tracer TR( "apps.torsion_types" );

int
main( int argc, char * argv [] ){
	try{
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace core::chemical;
		
		devel::init(argc, argv);
		
		utility::options::FileVectorOption & fvec 
			= basic::options::option[ basic::options::OptionKeys::in::file::extra_res_fa ];
			
		core::chemical::ChemicalManager * chem_mang = core::chemical::ChemicalManager::get_instance();
		core::chemical::AtomTypeSetCAP atom_types = chem_mang->atom_type_set("fa_standard");
		core::chemical::ElementSetCAP elements = chem_mang->element_set("default");
		core::chemical::MMAtomTypeSetCAP mm_atom_types = chem_mang->mm_atom_type_set("fa_standard");
		core::chemical::orbitals::OrbitalTypeSetCAP orbital_types = chem_mang->orbital_type_set("fa_standard");
		
		for ( core::Size i = 1, e = fvec.size(); i <= e; ++i ) {
			utility::file::FileName fname = fvec[i];
			std::string filename = fname.name();

			TR << "Processing " << filename << std::endl;
			core::chemical::ResidueTypeOP restype( read_topology_file(
				filenam, atom_types, elements, mm_atom_types, orbital_types ) );
			restype->print_dihedrals();

	} catch (utitlity::excn::Exception const & e) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}