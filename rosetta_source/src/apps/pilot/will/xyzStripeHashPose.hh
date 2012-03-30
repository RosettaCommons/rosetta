#ifndef INCLUDED_apps_pilot_will_xyzStripeHashPose_hh
#define INCLUDED_apps_pilot_will_xyzStripeHashPose_hh

#include <apps/pilot/will/xyzStripeHash.hh>

#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>

enum xyzStripeHashPoseMode {
	NBR,
	CB,
	BB,
	BNP,
	HVY,
	ALL
};

class xyzStripeHashPose : public xyzStripeHash<double> {
public:
  xyzStripeHashPose(double radius) : xyzStripeHash<double>(radius) {}
  xyzStripeHashPose(double radius, core::pose::Pose p, xyzStripeHashPoseMode m = BB ) : xyzStripeHash<double>(radius) {
	  init_with_pose(p,BB);
  }
  void init_with_pose(core::pose::Pose const & p, xyzStripeHashPoseMode m = BB) {
	  utility::vector1<double> dummy;
	  init_with_pose(p,dummy,BB);
  }
  void init_with_pose(core::pose::Pose const & p, utility::vector1<double> const & meta_in, xyzStripeHashPoseMode m = BB) {
    int natom = 0;
    for(int ir = 1; ir <= p.n_residue(); ++ir) {
      core::conformation::Residue const & r(p.residue(ir));
      if( NBR==m ) natom++;
      if( CB==m ) if(r.has("CB")) natom++;
      if( BB ==m ) natom += r.has("N")+r.has("CA")+r.has("C")+r.has("O")+r.has("CB");
      if( HVY==m ) natom += r.nheavyatoms();
      if( ALL==m ) natom += r.natoms();
    }
    utility::vector1<numeric::xyzVector<double> > atoms(natom);
    utility::vector1<double>                      meta (natom);
    uint count = 0;
    for(int ir = 1; ir <= p.n_residue(); ++ir) {
      core::conformation::Residue const & r(p.residue(ir));
      if(NBR==m) {
        int ia = r.nbr_atom();
        core::id::AtomID const aid(ia,ir);
        atoms[++count] = p.xyz(aid);
        meta [  count] = r.atom_type(ia).lj_radius();
      } else if(CB==m) {
        if(r.has("CB")){ atoms[++count]=r.xyz("CB"); meta[count]=r.atom_type(r.atom_index("CB")).lj_radius(); }
      } else if(BB==m) {
        if(r.has( "N")){ atoms[++count]=r.xyz( "N"); meta[count]=r.atom_type(r.atom_index( "N")).lj_radius(); }
        if(r.has("CA")){ atoms[++count]=r.xyz("CA"); meta[count]=r.atom_type(r.atom_index("CA")).lj_radius(); }
        if(r.has( "C")){ atoms[++count]=r.xyz( "C"); meta[count]=r.atom_type(r.atom_index( "C")).lj_radius(); }
        if(r.has( "O")){ atoms[++count]=r.xyz( "O"); meta[count]=r.atom_type(r.atom_index( "O")).lj_radius(); }
        if(r.has("CB")){ atoms[++count]=r.xyz("CB"); meta[count]=r.atom_type(r.atom_index("CB")).lj_radius(); }
      } else if(BB==m) {
        if(r.has("CA")){ atoms[++count]=r.xyz("CA"); meta[count]=r.atom_type(r.atom_index("CA")).lj_radius(); }
		if(r.has( "C")){ atoms[++count]=r.xyz( "C"); meta[count]=r.atom_type(r.atom_index( "C")).lj_radius(); }
        if(r.has("CB")){ atoms[++count]=r.xyz("CB"); meta[count]=r.atom_type(r.atom_index("CB")).lj_radius(); }
      } else {
        int natom = (ALL==m) ? r.natoms() : r.nheavyatoms();
        for(int ia = 1; ia <= natom; ++ia) {
          core::id::AtomID const aid(ia,ir);
          atoms[++count] = p.xyz(aid);
          meta [  count] = r.atom_type(ia).lj_radius();
        }
      }
    }
	if(meta_in.size()!=0) {
		init(atoms,meta_in);
	} else {
		init(atoms,meta);
	}
  }

};

#endif
