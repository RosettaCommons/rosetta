#ifdef MAC
#include <OpenCL/cl_platform.h>
#else
#include <CL/cl_platform.h>
#endif

void dump_unique_atoms() {
  core::pose::Pose p;
  core::pose::make_pose_from_sequence(p,"D",core::chemical::FA_STANDARD,false);
  remove_lower_terminus_type_from_pose_residue(p,1);
  remove_upper_terminus_type_from_pose_residue(p,1);
  core::kinematics::Stub s( p.residue(1).xyz("CB"), p.residue(1).xyz("CG"), p.residue(1).xyz("CB"), p.residue(1).xyz("CA") );
  //core::kinematics::Stub s( p.residue(1).xyz("CB"), p.residue(1).xyz("CA"), p.residue(1).xyz("N") );
  vector1<   Vec> seen_x;
  vector1<string> seen_s;
  std::string const aas("ACDEFGHIKLMNPQRSTVWY");
  for(std::string::const_iterator i = aas.begin(); i != aas.end(); ++i) {
    core::pose::make_pose_from_sequence(p,str(*i),core::chemical::FA_STANDARD,false);
    remove_lower_terminus_type_from_pose_residue(p,1);
    remove_upper_terminus_type_from_pose_residue(p,1);
    xform_pose_rev(p,s);
    for(int j = 6; j <= p.residue(1).nchi(); ++j) p.set_chi(j,1,0.0);
    for(int k = 1; k <= p.residue(1).nheavyatoms(); ++k) {
      Vec const x = p.residue(1).xyz(k);
      bool cont(false);      
      vector1<string>::iterator is = seen_s.begin();
      for(vector1<Vec>::const_iterator ix = seen_x.begin(); ix != seen_x.end(); ++ix,++is) {
        if( ix->distance_squared(x) < 0.0001 ) {
          cont = true;
          *is = (*is) + (*i);
        }
      }
      if(cont) continue;
      seen_x.push_back(x);
      seen_s.push_back(str(*i));
    }
  }
  seen_x.clear();
  for(std::string::const_iterator i = aas.begin(); i != aas.end(); ++i) {
    core::pose::make_pose_from_sequence(p,str(*i),core::chemical::FA_STANDARD,false);
    remove_lower_terminus_type_from_pose_residue(p,1);
    remove_upper_terminus_type_from_pose_residue(p,1);
    xform_pose_rev(p,s);
    for(int j = 6; j <= p.residue(1).nchi(); ++j) p.set_chi(j,1,0.0);
    p.dump_pdb("out/"+str(*i)+".pdb");
    for(int k = 1; k <= p.residue(1).nheavyatoms(); ++k) {
      Vec const x = p.residue(1).xyz(k);
      bool cont(false);      
      for(vector1<Vec>::const_iterator ix = seen_x.begin(); ix != seen_x.end(); ++ix) {
        if( ix->distance_squared(x) < 0.0001 ) {
          cont = true;
        }
      }
      if(cont) continue;
      string aname (p.residue(1).atom_name(k));
      aname.erase(remove_if(aname.begin(),aname.end(),isspace),aname.end());
      cout << "ATOM " << LJ(23,aname+"_"+seen_s[seen_x.size()+1]) << " " << F(10,7,x.x()) << " " << F(10,7,x.y()) << " " << F(10,7,x.z()) << std::endl;
      seen_x.push_back(x);
    }
    
  }
}


// point in POI octree x,y,z 8/8/5/11

/*
  ILE CG1 -> CG
  THR OG1 -> OG CG2 -> CG
  VAL CG1 <-> CG2 so that CG & CG2 are same as others
  LEU CD1 -> CD
  ILE CD1 -> CD

  GLY N CA C O
  ALA N CA C O  CB
  CYS N CA C O  CB                            SG
  ASN N CA C O -CB- CG      OD1               ND2
  ASP N CA C O -CB- CG      OD1               OD2
  GLN N CA C O -CB- CG  CD-                                             OE1          NE2
  GLU N CA C O -CB- CG  CD-                                             OE1          OE2
  VAL N CA C O -CB- CG               CG2
  LEU N CA C O -CB- CG  CD                    CD2
  ILE N CA C O -CB- CG  CD           CG2
  THR N CA C O -CB-           OG     CG2
  SER N CA C O -CB-           OG
  ARG N CA C O -CB- CG  CD                                                     NE  CZ NH1 NH2
  LYS N CA C O -CB- CG  CD                                                     CE  NZ
  MET N CA C O -CB- CG                        SD                               CE
  HIS N CA C O -CB- CG                        ND1 CD2 CE1 NE2
  TRP N CA C O -CB- CG                        CD1 CD2 NE1 CE2 CE3 CZ2 CZ3 CH2
  PHE N CA C O -CB- CG  CD1 CE1 CZ CD2 CE2
  TYR N CA C O -CB- CG  CD1 CE1 CZ CD2 CE2    OH

  CG
  CD
  OD1
  OG
  CG2
  CD1
  CE1
  CZ
  CD2
  CE2
  SG
  ND2
  OD2
  CD2_L
  SD
  ND1_H
  CD2_H
  CE1_N
  NE2_H
  OH
  CD1_W
  CD2_W
  NE1_W
  CE2_W
  CE3_W
  CZ2_W
  CZ3_W
  CH2_W
*/
