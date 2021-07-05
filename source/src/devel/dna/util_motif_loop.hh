#include <core/conformation/Residue.fwd.hh>
#include <core/pack/rotamer_set/RotamerSet.fwd.hh>

#include <core/pose/Pose.fwd.hh> // AUTO IWYU For Pose
#include <map> // AUTO IWYU For map

//Returns a a rotamer set of all Amino Acids
/////////////////////////////////////////////////////////////////////////////////////////////////
//////////Utility function///////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

bool hard_sphere_check( core::pose::Pose & pose, core::conformation::Residue & resi);

std::map< std::string, core::pack::rotamer_set::RotamerSetOP >  build_all_aa_rotamers();

core::pack::rotamer_set::RotamerSetOP get_aa_rotset( std::string & AA);

///Returns a new amino acid of type dictated by threeLetterAA code
core::conformation::Residue new_residue( std::string threeLetterAA );


