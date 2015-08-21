#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AtomType.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/PackerTask.hh>

//Returns a a rotamer set of all Amino Acids
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////Utility function///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

bool hard_sphere_check( core::pose::Pose & pose, core::conformation::Residue & resi);

std::map< std::string, core::pack::rotamer_set::RotamerSetOP >  build_all_aa_rotamers();

core::pack::rotamer_set::RotamerSetOP get_aa_rotset( std::string & AA);

///Returns a new amino acid of type dictated by threeLetterAA code
core::conformation::Residue new_residue( std::string threeLetterAA );


