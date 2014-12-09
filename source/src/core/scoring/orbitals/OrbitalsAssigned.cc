// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

///@author Steven Combs

#include <core/scoring/orbitals/OrbitalsLookup.hh>
#include <core/scoring/orbitals/OrbitalsAssigned.hh>
#include <core/scoring/orbitals/OrbitalsScore.hh>
#include <core/scoring/orbitals/OrbitalsScoreCreator.hh>
#include <core/scoring/orbitals/OrbitalsAssigned.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <numeric/constants.hh>
#include <core/chemical/AtomType.hh>

//option headers
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


namespace core{
namespace scoring{
namespace orbitals{




utility::vector1< numeric::xyzVector<core::Real> > OrbitalsAssigned::get_lp_xyz(
		core::conformation::Residue const & residue
)
{
	core::Real dist = 0.7;
	utility::vector1< numeric::xyzVector<core::Real> > lp_xyz;
	utility::vector1< numeric::xyzVector<core::Real> > lp_holder_xyz;
	numeric::xyzVector<core::Real> vector_d;
	numeric::xyzVector<core::Real> vector_f;


	if(residue.is_aromatic() ){
		vector_d = residue.atom("CG").xyz() - residue.atom("CD1").xyz();
		vector_f = residue.atom("CG").xyz() - residue.atom("CD2").xyz();
		lp_holder_xyz = cp_function("aroC", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
		//utility::vector1< numeric::xyzVector<core::Real> > lp_holder_xyz2 = aromatic_ring_center(vector_d, vector_f, residue, dist);
		//lp_xyz.insert(lp_xyz.end(), lp_holder_xyz2.begin(), lp_holder_xyz2.end());
		if(residue.name3() == "TRP"){
				lp_holder_xyz = cp_function("Ntrp", vector_d, vector_f, residue, dist);
				lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());


			}
		if(residue.name3() == "TYR"){
				lp_holder_xyz = (CoordinatesTetrahedral("OH", "CZ", "HH", dist, residue));
				lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
			}
	}
	else if(residue.name3() == "GLU"){
		vector_d = residue.atom("CD").xyz() - residue.atom("OE1").xyz();
		vector_f = residue.atom("CD").xyz() - residue.atom("OE2").xyz();

		lp_holder_xyz = cp_function("COO", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = cp_function("OOC", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = CoordinatesDihedral("OE1", "CD", "OE2", dist,residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = CoordinatesDihedral("OE2", "CD", "OE1", dist,residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

	}
	else if(residue.name3() == "ASP"){
		vector_d = residue.atom("CG").xyz() - residue.atom("OD1").xyz();
		vector_f = residue.atom("CG").xyz() - residue.atom("OD2").xyz();

		lp_holder_xyz = cp_function("COO", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = cp_function("OOC", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = CoordinatesDihedral("OD1", "CG", "OD2", dist, residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = CoordinatesDihedral("OD2", "CG", "OD1", dist, residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
	}
	else if(residue.name3() == "ASN"){
		vector_d = residue.atom("CG").xyz() - residue.atom("OD1").xyz();
		vector_f = residue.atom("CG").xyz() - residue.atom("ND2").xyz();
		lp_holder_xyz = cp_function("CNH2", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = cp_function("ONH2", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());


		lp_holder_xyz = CoordinatesDihedral("OD1", "CG", "ND2", dist, residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		//vector_d = residue.atom("ND2").xyz() - residue.atom("CG").xyz();
		//vector_f = residue.atom("ND2").xyz() - residue.atom("OD1").xyz();
		lp_holder_xyz = cp_function("NH2O",  vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());



	}
	else if(residue.name3() == "GLN"){
		vector_d = residue.atom("CD").xyz() - residue.atom("OE1").xyz();
		vector_f = residue.atom("CD").xyz() - residue.atom("NE2").xyz();

		lp_holder_xyz = cp_function("CNH2", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = cp_function("ONH2", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = CoordinatesDihedral("OE1", "CD", "NE2", dist, residue);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		//vector_d = residue.atom("NE2").xyz() - residue.atom("CD").xyz();
		//vector_f = residue.atom("NE2").xyz() - residue.atom("OE1").xyz();
		lp_holder_xyz = cp_function("NH2O", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

	}
	else if(residue.name3() == "HIS"){

		vector_d = residue.atom("ND1").xyz() - residue.atom("CE1").xyz();
		vector_f = residue.atom("ND1").xyz() - residue.atom("CG").xyz();


		//vector_d = residue.atom("CG").xyz() - residue.atom("CE1").xyz();
		//vector_f = residue.atom("CG").xyz() - residue.atom("CD2").xyz();
		lp_holder_xyz = cp_function("aroC", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
		//utility::vector1< numeric::xyzVector<core::Real> > lp_holder_xyz2 = aromatic_ring_center(vector_d, vector_f, residue, dist);
		//lp_xyz.insert(lp_xyz.end(), lp_holder_xyz2.begin(), lp_holder_xyz2.end());


		numeric::xyzVector<core::Real> vector_df_norm = vector_d.normalized()+vector_f.normalized();

		lp_holder_xyz = cp_function("Nhis", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_holder_xyz = cp_function("Ntrp", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());

		lp_xyz.push_back((vector_df_norm.normalized()*dist)+residue.atom("ND1").xyz());

		numeric::xyzVector<core::Real> vector_a = residue.atom("NE2").xyz() - residue.atom("CD2").xyz();
		numeric::xyzVector<core::Real> vector_b = residue.atom("NE2").xyz() - residue.atom("CE1").xyz();
		numeric::xyzVector<core::Real> vector_ab_norm = vector_a.normalized()+vector_b.normalized();

		lp_xyz.push_back((vector_ab_norm.normalized()*dist)+residue.atom("NE2").xyz());

	}



	else if(residue.name3() == "ARG"){
		vector_d = residue.atom("NH1").xyz() - residue.atom("CZ").xyz();
		vector_f = residue.atom("NH1").xyz() - residue.atom("NH2").xyz();
		lp_holder_xyz = cp_function("Narg", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
		lp_holder_xyz = cp_function("aroC", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());


	}
	else if(residue.name3() == "PRO"){
		vector_d = residue.atom("N").xyz() - residue.atom("CA").xyz();
		vector_f = residue.atom("N").xyz() - residue.atom("CD").xyz();
		lp_holder_xyz = cp_function("Npro", vector_d, vector_f, residue, dist);
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());


	}

	else if(residue.name3() == "SER"){
		lp_holder_xyz = (CoordinatesTetrahedral("OG", "CB", "HG", dist, residue));
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
	}
	else if(residue.name3() == "THR"){
		lp_holder_xyz = (CoordinatesTetrahedral("OG1", "CB", "HG1", dist, residue));
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
	}


	else if(residue.name3() == "CYS"){
		lp_holder_xyz = (CoordinatesTetrahedral("SG", "CB", "HG", dist, residue));
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
	}
	else if(residue.name3() == "MET"){
		lp_holder_xyz = (CoordinatesTetrahedral("SD", "CG", "CE", dist, residue));
		lp_xyz.insert(lp_xyz.end(), lp_holder_xyz.begin(), lp_holder_xyz.end());
	}

		std::cout << residue.name3() << residue.seqpos() << std::endl;
		for(core::Real i = 1; i <= lp_xyz.size(); ++i){
			printf("ATOM      1  H   FLR D %3i      %3.3f  %3.3f  %3.3f \n", residue.seqpos(), lp_xyz[i].x(), lp_xyz[i].y(), lp_xyz[i].z() );
		}


	return lp_xyz;
}

//Taken from BCL!!!
utility::vector1< numeric::xyzVector<core::Real> >  OrbitalsAssigned::CoordinatesDihedral(
		std::string atom1,
		std::string atom2,
		std::string atom3,
		core::Real distance_xa,
		core::conformation::Residue const & residue
){
	utility::vector1< numeric::xyzVector<core::Real> > lp_xyz;

	numeric::xyzVector<core::Real> vector_a = residue.atom(atom1).xyz();
	numeric::xyzVector<core::Real> vector_b = residue.atom(atom2).xyz();
	numeric::xyzVector<core::Real> vector_c = residue.atom(atom3).xyz();

	//core::Real distance_xa = 01.0;
	core::Real angle_xab = numeric::constants::r::pi_2_over_3; //120 degrees
	//for one point it should be 180 for another it should be 0
	//core::Real dihedral_xabc = numeric::constants::r::pi;

	numeric::xyzVector<core::Real> a( (vector_a - vector_b).normalized());
	numeric::xyzVector<core::Real> b((vector_b - vector_c).normalized());
	numeric::xyzVector<core::Real> c(cross_product(a,b).normalized());
	numeric::xyzVector<core::Real> d(cross_product(a,c).normalized());

	core::Real dihedral_xabc1 = numeric::constants::r::pi;
	core::Real dihedral_xabc2 = 0;

	numeric::xyzVector<core::Real> v1 = a * std::cos(numeric::constants::r::pi - angle_xab);
	numeric::xyzVector<core::Real> v2 = c * std::sin(numeric::constants::r::pi - angle_xab);
	numeric::xyzVector<core::Real> v3 = d * std::sin(numeric::constants::r::pi - angle_xab);

	numeric::xyzVector<core::Real> x1
	(
			(
					v1 - v2 * std::sin(dihedral_xabc1) +	v3 * std::cos(dihedral_xabc1)
			) * distance_xa
	);

	numeric::xyzVector<core::Real> x2
	(
			(
					v1 - v2 * std::sin(dihedral_xabc2) +	v3 * std::cos(dihedral_xabc2)
			) * distance_xa
	);


	lp_xyz.push_back(vector_a+x1);
	lp_xyz.push_back(vector_a+x2);
	return lp_xyz;
}

//Taken from the BCL!!
utility::vector1< numeric::xyzVector<core::Real> >  OrbitalsAssigned::CoordinatesTetrahedral(
		std::string atom1,
		std::string atom2,
		std::string atom3,
		core::Real distance,
		core::conformation::Residue const & residue

){
	utility::vector1< numeric::xyzVector<core::Real> > lp_xyz;
	numeric::xyzVector<core::Real> foot_point;
	numeric::xyzVector<core::Real> offset;

	numeric::xyzVector<core::Real> vector_a = residue.atom(atom1).xyz();
	numeric::xyzVector<core::Real> vector_b = residue.atom(atom2).xyz();
	numeric::xyzVector<core::Real> vector_c = residue.atom(atom3).xyz();

	core::Real distance_xa = (distance * std::cos( 54.75 / 180 * numeric::constants::r::pi));

	numeric::xyzVector<core::Real> x(  (((vector_a - vector_b).normalized() + (vector_a - vector_c).normalized() ).normalized() * distance_xa));

	foot_point = (vector_a + x);

	core::Real distance_x_a = (distance * std::sin( 54.75 / 180 * numeric::constants::r::pi));

	offset = (distance_x_a * cross_product(vector_a-vector_b, vector_a-vector_c).normalized());

	lp_xyz.push_back(foot_point + offset);
	lp_xyz.push_back(foot_point - offset);

	return lp_xyz;


}





utility::vector1< numeric::xyzVector<core::Real> > OrbitalsAssigned::aromatic_ring_center(
		numeric::xyzVector<core::Real> vector_d,
		numeric::xyzVector<core::Real> vector_f,
		core::conformation::Residue const & residue,
		core::Real dist

){
	utility::vector1< numeric::xyzVector<core::Real> > lp_xyz;
	numeric::xyzVector<core::Real> lp_holder_xyz_right = cross_product(vector_d, vector_f);
	numeric::xyzVector<core::Real> lp_holder_xyz_left = cross_product(-vector_d, vector_f);
	numeric::xyzVector<core::Real> cg_right;
	numeric::xyzVector<core::Real> cz_right;
	numeric::xyzVector<core::Real> cg_left;
	numeric::xyzVector<core::Real> cz_left;
	numeric::xyzVector<core::Real> ce2_r;
	numeric::xyzVector<core::Real> cd2_r;
	numeric::xyzVector<core::Real> cz3_r;
	numeric::xyzVector<core::Real> cd1_r;
	numeric::xyzVector<core::Real> ce2_l;
	numeric::xyzVector<core::Real> cd2_l;
	numeric::xyzVector<core::Real> cz3_l;
	numeric::xyzVector<core::Real> cd1_l;
	for(core::Real i=1; i <= residue.atoms().size(); ++i){
		if(residue.atom_name(i) == " CZ " && residue.name3() == "TYR" || residue.atom_name(i) == " CG " && residue.name3() == "PHE"){
			cz_right = residue.atom(i).xyz();
			cz_left = residue.atom(i).xyz();
		}
		if(
				residue.atom_name(i) == " CG " && residue.name3() == "TYR" || residue.atom_name(i) == " CG " && residue.name3() == "PHE"
		){
			cg_right = residue.atom(i).xyz();
			cg_left = residue.atom(i).xyz();

		}
		if(residue.name3() == "TRP"){
			if(residue.atom_name(i) == " CE2"){
				ce2_r = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CD2"){
				cd2_r = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CZ3"){
				cz3_r = residue.atom(i).xyz();

			}
			if(residue.atom_name(i) == " CD1"){
				cd1_r = residue.atom(i).xyz();

			}
			if(residue.atom_name(i) == " CE2"){
				ce2_l = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CD2"){
				cd2_l = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CZ3"){
				cz3_l = residue.atom(i).xyz();

			}
			if(residue.atom_name(i) == " CD1"){
				cd1_l = residue.atom(i).xyz();
			}
		}
		if(residue.name3() == "HIS"){
			if(residue.atom_name(i) == " CD2"){
				cd2_r = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CG "){
				cg_right = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CE1"){
				ce2_r = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CD2"){
				cd2_l = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CG "){
				cg_left = residue.atom(i).xyz();
			}
			if(residue.atom_name(i) == " CE1"){
				ce2_l = residue.atom(i).xyz();
			}
		}
	}

		if(residue.name3() == "TYR" || residue.name3() == "PHE"){
			lp_xyz.push_back( ((cg_right+cz_right)/2)  + (lp_holder_xyz_right.normalized() * dist));
			lp_xyz.push_back( ((cg_left+cz_left)/2) + (lp_holder_xyz_left.normalized() * dist));
		}
		if(residue.name3() == "TRP"){
			lp_xyz.push_back( (ce2_r+cz3_r)/2 + (lp_holder_xyz_right.normalized() * dist));
			lp_xyz.push_back(  ((ce2_r+cd2_r)/2 + cd1_r)/2  + (lp_holder_xyz_right.normalized() * dist));
			lp_xyz.push_back( (ce2_l+cz3_l)/2  + (lp_holder_xyz_left.normalized() * dist));
			lp_xyz.push_back(  ((ce2_l+cd2_l)/2 + cd1_l)/2  + (lp_holder_xyz_left.normalized() * dist));
		}
		if(residue.name3() == "HIS"){
			lp_xyz.push_back( ((cd2_r+cg_right)/2 + ce2_r)/2  + (lp_holder_xyz_right.normalized() * dist));
			lp_xyz.push_back( ((cd2_l+cg_left)/2 + ce2_l)/2 + (lp_holder_xyz_left.normalized() * dist)  );
		}
		return lp_xyz;
}







utility::vector1< numeric::xyzVector<core::Real> > OrbitalsAssigned::cp_function(
		std::string atomtype,
		numeric::xyzVector<core::Real> vector_d,
		numeric::xyzVector<core::Real> vector_f,
		core::conformation::Residue const & residue,
		core::Real dist

){
	utility::vector1< numeric::xyzVector<core::Real> > lp_xyz;
	numeric::xyzVector<core::Real> lp_holder_xyz_right = cross_product(vector_d, vector_f);
	numeric::xyzVector<core::Real> lp_holder_xyz_left = cross_product(-vector_d, vector_f);
	for(core::Real i=1; i <= residue.atoms().size(); ++i)
	{
		if( residue.atom_type(i).name() == atomtype){
			lp_xyz.push_back( (lp_holder_xyz_right.normalized() * dist)+residue.atom(i).xyz());
			lp_xyz.push_back( (lp_holder_xyz_left.normalized() * dist)+residue.atom(i).xyz());

		}
	}
	return lp_xyz;
}

utility::vector1< std::pair< numeric::xyzVector< core::Real >,  std::string > > OrbitalsAssigned::get_hydrogens(
		core::conformation::Residue const & resid1
){
	utility::vector1< std::pair< numeric::xyzVector< core::Real >,  std::string > > resid1_H_xyz;
	for(core::Real i=1; i <= resid1.atoms().size(); ++i)
	{
		if( resid1.atom_type(i).name() == "Haro" ||  resid1.atom_type(i).name() == "Hpol"){
			std::pair< numeric::xyzVector< core::Real >, std::string > new_pair(resid1.atom(i).xyz(), resid1.atom_type(i).name());
			resid1_H_xyz.push_back( new_pair );

		}
	}
	return resid1_H_xyz;
}

}
}
}
