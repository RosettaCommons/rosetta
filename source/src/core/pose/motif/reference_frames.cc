// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/motif/reference_frames.hh
/// @brief  compute motif-related reference frames for residues, chemical groups, etc
/// @author Will Sheffler

// Unit headers

#include <core/pose/motif/reference_frames.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>
#include <utility/tools/make_vector1.hh>

namespace core {
namespace pose {
namespace motif {

using numeric::xyzVector;
using numeric::xyzTransform;
using namespace core::chemical;
using core::Size;
using core::Real;
using std::string;
using Xform = xyzTransform<Real>;
using Vec = xyzVector<Real>;
using PoseCR = const Pose &;
using SizeCR = const Size &;
using utility::tools::make_vector1;
using core::id::AtomID;
using std::cout;
using std::endl;
using AIDs = utility::vector1<core::id::AtomID>;


static
AIDs get_AIDs(
	Pose const & pose,
	core::Size const & ir,
	utility::vector1<std::string> const & anames
){
	AIDs aids;
	for ( std::string const & aname : anames ) {
		if ( pose.residue_type(ir).has(aname) ) {
			aids.push_back( AtomID(pose.residue_type(ir).atom_index(aname),ir) );
		} else {
			std::cout << "WARNING! get_sidechain_reference_frame_downstream_atoms: res "<<ir<<" "<< pose.residue_type(ir).name() << " has no atom named '" << aname<<"'!" << std::endl;
			aids.push_back( core::id::AtomID::BOGUS_ATOM_ID() );
			utility_exit_with_message("atom name mismatch somewhere!!!");
		}
	}
	return aids;
}

Vec a(Pose const & p, numeric::Size const & i, string const & s){
	if ( ! p.residue(i).has(s) ) utility_exit_with_message("no "+s+" atom in residue "+p.residue(i).name()+"!");
	return p.residue(i).xyz(s);
}

Xform get_backbone_reference_frame(
	xyzVector<Real> const & N,
	xyzVector<Real> const & CA,
	xyzVector<Real> const & C
){
	Vec CEN(-0.865810,-1.764143,1.524857); // average 'CEN' icoor
	CEN = Xform().from_four_points( CA, N, CA, C ) * CEN;
	Vec const DIR1 = C-N;
	Vec const CEN2 = (C+N)/2;
	return Xform().from_four_points( CEN, CEN2, CA, CA+DIR1 );
}
Xform get_backbone_reference_frame(Pose const & p, numeric::Size const & i){
	runtime_assert( i >= 1 );
	runtime_assert( i <= p.size() );
	if ( !p.residue(i).is_protein() ) return Xform::BAD_XFORM();
	return get_backbone_reference_frame( a(p,i,"N"), a(p,i,"CA"), a(p,i,"C") );
}


Xform get_nterminal_peptide_bond_reference_frame(
	xyzVector<Real> const & H,
	xyzVector<Real> const & N,
	xyzVector<Real> const & CA
){
	return Xform().from_four_points( H, H , N, CA );
}
Xform get_nterminal_peptide_bond_reference_frame(Pose const & p, numeric::Size const & i){
	if ( !p.residue(i).is_protein() ) return Xform::BAD_XFORM();
	if ( p.residue(i).is_lower_terminus() ) return Xform::BAD_XFORM();
	// if( 1==i ) return Xform::BAD_XFORM();
	if ( p.residue(i).aa()==aa_pro ) return Xform::BAD_XFORM();
	return get_nterminal_peptide_bond_reference_frame( a(p,i, "H" ), a(p,i,"N"), a(p,i,"CA") );
}

Xform get_cterminal_peptide_bond_reference_frame(
	xyzVector<Real> const & O,
	xyzVector<Real> const & C,
	xyzVector<Real> const & CA
){
	return Xform().from_four_points( O, O, C, CA );
}
Xform get_cterminal_peptide_bond_reference_frame(Pose const & p, numeric::Size const & i){
	if ( !p.residue(i).is_protein() ) return Xform::BAD_XFORM();
	if ( p.residue(i).is_upper_terminus() ) return Xform::BAD_XFORM();
	// if( i==p.size() ) return Xform::BAD_XFORM();
	return get_cterminal_peptide_bond_reference_frame( a(p,i,"O"), a(p,i,"C"), a(p,i,"CA") );
}


AIDs
get_nterminal_peptide_bond_reference_frame_atomids(Pose const & p, numeric::Size const & i, bool extra){
	if ( !p.residue(i).is_protein()       ) return AIDs();
	if ( p.residue(i).is_lower_terminus() ) return AIDs();
	if ( p.residue(i).aa()==aa_pro ) return AIDs();
	AIDs aids = get_AIDs(p,i,make_vector1<string>("H","N","CA"));
	if ( extra && i > 1 ) aids.push_back( AtomID(p.residue(i-1).atom_index("C"),i-1) );
	return aids;
}
AIDs
get_cterminal_peptide_bond_reference_frame_atomids(Pose const & p, numeric::Size const & i, bool extra){
	if ( !p.residue(i).is_protein()       ) return AIDs();
	if ( p.residue(i).is_upper_terminus() ) return AIDs();
	AIDs aids = get_AIDs(p,i,make_vector1<string>("CA","C","O"));
	if ( extra && i < p.size() ) aids.push_back( AtomID(p.residue(i+1).atom_index("N"),i+1) );
	return aids;
}


Xform get_sidechain_reference_frame(Pose const & p, numeric::Size const & i){
	runtime_assert( i >= 1 );
	runtime_assert( i <= p.size() );
	switch(p.residue(i).aa()){
	case aa_ala : return get_frame_ala(p,i);
	case aa_cys : return get_frame_cys(p,i);
	case aa_asp : return get_frame_asp(p,i);
	case aa_glu : return get_frame_glu(p,i);
	case aa_phe : return get_frame_phe(p,i);
	case aa_gly : return get_frame_gly(p,i);
	case aa_his : return get_frame_his(p,i);
	case aa_ile : return get_frame_ile(p,i);
	case aa_lys : return get_frame_lys(p,i);
	case aa_leu : return get_frame_leu(p,i);
	case aa_met : return get_frame_met(p,i);
	case aa_asn : return get_frame_asn(p,i);
	case aa_pro : return get_frame_pro(p,i);
	case aa_gln : return get_frame_gln(p,i);
	case aa_arg : return get_frame_arg(p,i);
	case aa_ser : return get_frame_ser(p,i);
	case aa_thr : return get_frame_thr(p,i);
	case aa_val : return get_frame_val(p,i);
	case aa_trp : return get_frame_trp(p,i);
	case aa_tyr : return get_frame_tyr(p,i);
	default :
		return Xform::BAD_XFORM();
		//utility_exit_with_message("don't know what to do with unknown AA");
	}
	return Xform::BAD_XFORM();
}

Xform get_frame_ala(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_ala); return Xform().from_four_points(a(p,i,"CB" ),a(p,i,"CB" ),a(p,i,"CA" ),a(p,i,"N"  )); }
Xform get_frame_cys(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_cys); return Xform().from_four_points(a(p,i,"SG" ),a(p,i,"SG" ),a(p,i,"CB" ),a(p,i,"CA" )); }
Xform get_frame_gly(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_gly); return Xform().from_four_points(a(p,i,"CA" ),a(p,i,"N"  ),a(p,i,"CA" ),a(p,i,"C"  )); }
Xform get_frame_ile(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_ile); return Xform().from_four_points(a(p,i,"CG1"),a(p,i,"CD1"),a(p,i,"CG1"),a(p,i,"CB" )); }
Xform get_frame_lys(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_lys); return Xform().from_four_points(a(p,i,"NZ" ),a(p,i,"NZ" ),a(p,i,"CE" ),a(p,i,"CD" )); }
Xform get_frame_met(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_met); return Xform().from_four_points(a(p,i,"SD" ),a(p,i,"CE" ),a(p,i,"SD" ),a(p,i,"CG" )); }
Xform get_frame_asn(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_asn); return Xform().from_four_points(a(p,i,"CG" ),a(p,i,"OD1"),a(p,i,"CG" ),a(p,i,"ND2")); }
Xform get_frame_pro(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_pro); return Xform().from_four_points(a(p,i,"CG" ),a(p,i,"CD" ),a(p,i,"CG" ),a(p,i,"CB" )); }
Xform get_frame_gln(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_gln); return Xform().from_four_points(a(p,i,"CD" ),a(p,i,"OE1"),a(p,i,"CD" ),a(p,i,"NE2")); }
Xform get_frame_arg(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_arg); return Xform().from_four_points(a(p,i,"CZ" ),a(p,i,"NH1"),a(p,i,"CZ" ),a(p,i,"NH2")); }
Xform get_frame_ser(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_ser); return Xform().from_four_points(a(p,i,"OG" ),a(p,i,"OG" ),a(p,i,"CB" ),a(p,i,"CA" )); }
Xform get_frame_thr(PoseCR p, SizeCR i){ runtime_assert( p.residue(i).aa() == aa_thr); return Xform().from_four_points(a(p,i,"OG1"),a(p,i,"OG1"),a(p,i,"CB" ),a(p,i,"CG2")); }

Xform get_frame_vee(Vec const & a, Vec const & b, Vec const & c){
	Vec const CEN = (a+b)/2.0;
	Vec const DIR = CEN-c;
	Vec const DR2 = a-b;
	return Xform().from_four_points(CEN,CEN+DIR,CEN,CEN-DR2);
}

Xform get_frame_val(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_val);
	return get_frame_vee( a(p,i,"CG1"), a(p,i,"CG2"), a(p,i,"CB"));
}
Xform get_frame_leu(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_leu);
	return get_frame_vee( a(p,i,"CD1"), a(p,i,"CD2"), a(p,i,"CG"));
}
Xform get_frame_glu(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_glu);
	return get_frame_vee( a(p,i,"OE1"), a(p,i,"OE2"), a(p,i,"CD"));
}
Xform get_frame_asp(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_asp);
	return get_frame_vee( a(p,i,"OD2"), a(p,i,"OD1"), a(p,i,"CG"));
}

Xform get_frame_his(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_his);
	Vec CEN((a(p,i,"CG" )+a(p,i,"ND1")+a(p,i,"CE1")+a(p,i,"CD2")+a(p,i,"NE2"))/5.0);
	return Xform().from_four_points(CEN,a(p,i,"NE2"),CEN,a(p,i,"ND1"));
}
Xform get_frame_phe(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_phe);
	Vec CEN((a(p,i,"CG")+a(p,i,"CD1")+a(p,i,"CE1")+a(p,i,"CD2")+a(p,i,"CE2")+a(p,i,"CZ" ))/6.0);
	Vec OUT((a(p,i,"CG")-CEN).cross(a(p,i,"CD1")-CEN));
	return Xform().from_four_points( CEN, a(p,i,"CG"), CEN, CEN-OUT );
}
Xform get_frame_tyr(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_tyr);
	Vec CEN((a(p,i,"CG" )+a(p,i,"CD1")+a(p,i,"CE1")+a(p,i,"CD2")+a(p,i,"CE2")+a(p,i,"CZ" ))/6.0);
	Vec OUT((a(p,i,"CG" )-CEN).cross(a(p,i,"CD1")-CEN));
	return Xform().from_four_points( CEN, a(p,i,"CG"), CEN, CEN-OUT );
}
Xform get_frame_trp(PoseCR p, SizeCR i){
	runtime_assert( p.residue(i).aa() == aa_trp);
	Vec CEN((a(p,i,"CE2")+a(p,i,"CD2")+a(p,i,"CE3")+a(p,i,"CZ3")+a(p,i,"CH2")+a(p,i,"CZ2"))/6.0);
	Vec OUT((a(p,i,"CE2")-CEN).cross(a(p,i,"CD2")-CEN));
	return Xform().from_four_points( CEN, CEN+OUT, CEN,a(p,i,"CD2") );
}


AIDs
get_backbone_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir
){
	runtime_assert( ir >= 1 );
	runtime_assert( ir <= pose.size() );
	AIDs aids;
	aids.push_back(AtomID(1,ir));
	aids.push_back(AtomID(2,ir));
	aids.push_back(AtomID(3,ir));
	if ( pose.residue(ir).aa()!=aa_gly ) aids.push_back(AtomID(5,ir));
	else                              aids.push_back(id::AtomID::BOGUS_ATOM_ID());
	return aids;
}

AIDs
get_backbone_reference_frame_atomids_with_downstream(
	Pose const & pose,
	core::Size const & ir
){
	runtime_assert( ir >= 1 );
	runtime_assert( ir <= pose.size() );
	AIDs aids;
	aids.push_back(AtomID(1,ir));
	aids.push_back(AtomID(2,ir));
	aids.push_back(AtomID(3,ir));
	for ( Size ia = 5; ia <= pose.residue_type(ir).natoms(); ++ia ) {
		// if( pose.residue_type(ir).atom_is_polar_hydrogen(ia) ){
		//  for(Size j = 1; j <= pose.residue_type(ir).natoms(); ++j){
		//   cout << pose.residue(ir).name() << " " << pose.residue(ir).atom_name(j) << ir << " " << j << " " << pose.residue_type(ir).atom_is_polar_hydrogen(j) << endl;
		//  }
		//  cout << endl;
		// }
		// std::cout << "MASK CHECK: "
		// << pose.residue_type(ir).atom_name(ia) << " "
		// << pose.residue_type(ir).atom_is_backbone(ia) << " "
		// << pose.residue_type(ir).atom_is_hydrogen(ia) << " "
		// << pose.residue_type(ir).atom_is_polar_hydrogen(ia) << " "
		// << std::endl;
		if ( pose.residue_type(ir).atom_is_backbone(ia) ) continue;
		if ( pose.residue_type(ir).atom_is_hydrogen(ia) && !pose.residue_type(ir).atom_is_polar_hydrogen(ia) ) continue;
		aids.push_back(AtomID(ia,ir));
	}
	return aids;
}

AIDs
get_sidechain_reference_frame_atomids(
	Pose const & pose,
	core::Size const & ir
){
	runtime_assert( ir >= 1 );
	runtime_assert( ir <= pose.size() );
	switch(pose.residue(ir).aa()){
	case aa_ala : return get_atoms_ala(pose,ir);
	case aa_cys : return get_atoms_cys(pose,ir);
	case aa_asp : return get_atoms_asp(pose,ir);
	case aa_glu : return get_atoms_glu(pose,ir);
	case aa_phe : return get_atoms_phe(pose,ir);
	case aa_gly : return get_atoms_gly(pose,ir);
	case aa_his : return get_atoms_his(pose,ir);
	case aa_ile : return get_atoms_ile(pose,ir);
	case aa_lys : return get_atoms_lys(pose,ir);
	case aa_leu : return get_atoms_leu(pose,ir);
	case aa_met : return get_atoms_met(pose,ir);
	case aa_asn : return get_atoms_asn(pose,ir);
	case aa_pro : return get_atoms_pro(pose,ir);
	case aa_gln : return get_atoms_gln(pose,ir);
	case aa_arg : return get_atoms_arg(pose,ir);
	case aa_ser : return get_atoms_ser(pose,ir);
	case aa_thr : return get_atoms_thr(pose,ir);
	case aa_val : return get_atoms_val(pose,ir);
	case aa_trp : return get_atoms_trp(pose,ir);
	case aa_tyr : return get_atoms_tyr(pose,ir);
	default : utility_exit_with_message("don't know what to do with unknown AA");
	}
	utility_exit_with_message("don't know what to do with unknown AA");
}

AIDs
get_sidechain_reference_frame_atomids_with_downstream(
	Pose const & pose,
	core::Size const & ir
){
	runtime_assert( ir >= 1 );
	runtime_assert( ir <= pose.size() );
	switch(pose.residue(ir).aa()){
	case aa_ala : return get_atoms_ala_downstream(pose,ir);
	case aa_cys : return get_atoms_cys_downstream(pose,ir);
	case aa_asp : return get_atoms_asp_downstream(pose,ir);
	case aa_glu : return get_atoms_glu_downstream(pose,ir);
	case aa_phe : return get_atoms_phe_downstream(pose,ir);
	case aa_gly : return get_atoms_gly_downstream(pose,ir);
	case aa_his : return get_atoms_his_downstream(pose,ir);
	case aa_ile : return get_atoms_ile_downstream(pose,ir);
	case aa_lys : return get_atoms_lys_downstream(pose,ir);
	case aa_leu : return get_atoms_leu_downstream(pose,ir);
	case aa_met : return get_atoms_met_downstream(pose,ir);
	case aa_asn : return get_atoms_asn_downstream(pose,ir);
	case aa_pro : return get_atoms_pro_downstream(pose,ir);
	case aa_gln : return get_atoms_gln_downstream(pose,ir);
	case aa_arg : return get_atoms_arg_downstream(pose,ir);
	case aa_ser : return get_atoms_ser_downstream(pose,ir);
	case aa_thr : return get_atoms_thr_downstream(pose,ir);
	case aa_val : return get_atoms_val_downstream(pose,ir);
	case aa_trp : return get_atoms_trp_downstream(pose,ir);
	case aa_tyr : return get_atoms_tyr_downstream(pose,ir);
	default : utility_exit_with_message("don't know what to do with unknown AA");
	}
	utility_exit_with_message("don't know what to do with unknown AA");
}

AIDs get_atoms_ala(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CB" ,"CA" ,"N","C"                                      )); }
AIDs get_atoms_arg(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("NH1","CZ" ,"NH2","NE","CD"                              )); }
AIDs get_atoms_asn(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OD1","CG" ,"ND2","CB"                                   )); }
AIDs get_atoms_asp(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CG","CB"                                    )); }
AIDs get_atoms_cys(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("SG" ,"CB" ,"CA"                                         )); }
AIDs get_atoms_gln(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OE1","CD" ,"NE2","CG"                                   )); }
AIDs get_atoms_glu(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD","CG"                                    )); }
AIDs get_atoms_gly(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("N"  ,"CA" ,"C"                                          )); }
AIDs get_atoms_his(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD2","NE2","CE1","ND1","CG","CB"                        )); }
AIDs get_atoms_ile(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD1","CG1","CB"                                         )); }
AIDs get_atoms_leu(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CG","CB"                                    )); }
AIDs get_atoms_lys(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("NZ" ,"CE" ,"CD"                                         )); }
AIDs get_atoms_met(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CE" ,"SD" ,"CG"                                         )); }
AIDs get_atoms_phe(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CZ","CG","CB"                  )); }
AIDs get_atoms_pro(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD" ,"CG" ,"CB","CA","N","C"                            )); }
AIDs get_atoms_ser(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OG" ,"CB" ,"CA"                                         )); }
AIDs get_atoms_thr(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OG1","CB" ,"CG2","CA"                                   )); }
AIDs get_atoms_trp(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD2","CE3","CZ3","CH2","CZ2","CE2","NE1","CD1","CG","CB")); }
AIDs get_atoms_tyr(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OH","CZ","CG","CB"             )); }
AIDs get_atoms_val(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CB","CA"                                    )); }

AIDs get_atoms_ala_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CB" ,"CA" ,"N","C"                                                       )); }
AIDs get_atoms_arg_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("NH1","CZ" ,"NH2","NE","CD"              ,"1HH1","2HH1","1HH2","2HH2","HE")); }
AIDs get_atoms_asn_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OD1","CG" ,"ND2","CB"                                      ,"1HD2","2HD2")); }
AIDs get_atoms_asp_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OD1","OD2","CG","CB"                                                     )); }
AIDs get_atoms_cys_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("SG" ,"CB" ,"CA"                                                          )); }
AIDs get_atoms_gln_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OE1","CD" ,"NE2","CG"                                      ,"1HE2","2HE2")); }
AIDs get_atoms_glu_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OE1","OE2","CD","CG"                                                     )); }
AIDs get_atoms_gly_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("N"  ,"CA" ,"C"                                                           )); }
AIDs get_atoms_his_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD2","NE2","CE1","ND1","CG","CB"                                   ,"HE2")); }
AIDs get_atoms_ile_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD1","CG1","CB"                                                          )); }
AIDs get_atoms_leu_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD1","CD2","CG","CB"                                                     )); }
AIDs get_atoms_lys_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("NZ" ,"CE" ,"CD"                                        ,"1HZ","2HZ","3HZ")); }
AIDs get_atoms_met_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CE" ,"SD" ,"CG"                                                          )); }
AIDs get_atoms_phe_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CZ","CE1" ,"CE2","CD1","CD2","CG","CB"                                   )); }
AIDs get_atoms_pro_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD" ,"CG" ,"CB","CA","N","C"                                             )); }
AIDs get_atoms_ser_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OG" ,"CB" ,"CA"                                                     ,"HG")); }
AIDs get_atoms_thr_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OG1","CB" ,"CG2","CA"                                              ,"HG1")); }
AIDs get_atoms_trp_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CD2","CE3","CZ3","CH2","CZ2","CE2","NE1","CD1","CG","CB"           ,"HE1")); }
AIDs get_atoms_tyr_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("OH","CZ","CE1" ,"CE2","CD1","CD2","CG","CB"                         ,"HH")); }
AIDs get_atoms_val_downstream(PoseCR p, SizeCR i){ return get_AIDs(p,i,make_vector1<std::string>("CG1","CG2","CB","CA"                                                     )); }


} // motif
} // pose
} // core
