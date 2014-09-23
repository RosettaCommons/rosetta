// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/chiral/ChiralMover.cc
/// @brief ChiralMover methods implemented
/// @author Kevin Drew, kdrew@nyu.edu

// Unit Headers
#include <protocols/simple_moves/chiral/ChiralMover.hh>
// Package Headers

// Project Headers
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/Patch.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Residue.functions.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/util.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>
// Utility Headers
#include <numeric/xyz.functions.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
//#include <core/types.hh>

// C++ Headers

using basic::T;
using basic::Error;
using basic::Warning;

static thread_local basic::Tracer TR( "protocols.simple_moves.chiral.ChiralMover" );

using namespace core;
using namespace conformation;
using namespace chemical;
using namespace core::id;

namespace protocols {
namespace simple_moves {
namespace chiral {


std::pair< std::string, std::string > L2DChiralData[] = { 
			std::make_pair("A04","D04"),
			std::make_pair("A05","D05"),
			std::make_pair("A06","D06"),
			std::make_pair("A07","D07"),
			std::make_pair("A12","D12"),
			std::make_pair("A20","D20"),
			std::make_pair("A24","D24"),
			std::make_pair("A30","D30"),
			std::make_pair("A31","D31"),
			std::make_pair("A33","D33"),
			std::make_pair("A34","D34"),
			std::make_pair("A43","D43"),
			std::make_pair("A44","D44"),
			std::make_pair("A45","D45"),
			std::make_pair("A48","D48"),
			std::make_pair("A68","D68"),
			std::make_pair("A69","D69"),
			std::make_pair("A78","D78"),
			std::make_pair("A80","D80"),
			std::make_pair("A82","D82"),
			std::make_pair("A83","D83"),
			std::make_pair("A84","D84"),
			std::make_pair("A91","D91"),
			std::make_pair("A92","D92"),
			std::make_pair("A94","D94"),
			std::make_pair("A96","D96"),
			std::make_pair("A97","D97"),
			std::make_pair("A98","D98"),
			std::make_pair("A99","D99"),
			std::make_pair("ABA","DABA"),
			std::make_pair("ALA","DALA"),
			std::make_pair("APA","DAPA"),
			std::make_pair("ARG","DARG"),
			std::make_pair("ASN","DASN"),
			std::make_pair("ASP","DASP"),
			std::make_pair("B00","E00"),
			std::make_pair("B01","E01"),
			std::make_pair("B02","E02"),
			std::make_pair("B03","E03"),
			std::make_pair("B04","E04"),
			std::make_pair("B05","E05"),
			std::make_pair("B06","E06"),
			std::make_pair("B07","E07"),
			std::make_pair("B12","E12"),
			std::make_pair("B19","E19"),
			std::make_pair("B21","E21"),
			std::make_pair("B27","E27"),
			std::make_pair("B28","E28"),
			std::make_pair("B30","E30"),
			std::make_pair("B31","E31"),
			std::make_pair("B35","E35"),
			std::make_pair("B36","E36"),
			std::make_pair("B38","E38"),
			std::make_pair("B40","E40"),
			std::make_pair("B44","E44"),
			std::make_pair("B47","E47"),
			std::make_pair("B48","E48"),
			std::make_pair("B49","E49"),
			std::make_pair("B50","E50"),
			std::make_pair("B51","E51"),
			std::make_pair("B53","E53"),
			std::make_pair("B54","E54"),
			std::make_pair("B56","E56"),
			std::make_pair("B57","E57"),
			std::make_pair("B58","E58"),
			std::make_pair("B59","E59"),
			std::make_pair("B60","E60"),
			std::make_pair("B61","E61"),
			std::make_pair("B62","E62"),
			std::make_pair("B63","E63"),
			std::make_pair("B67","E67"),
			std::make_pair("B74","E74"),
			std::make_pair("B92","E92"),
			std::make_pair("B93","E93"),
			std::make_pair("B94","E94"),
			std::make_pair("B95","E95"),
			std::make_pair("B96","E96"),
			std::make_pair("B97","E97"),
			std::make_pair("B99","E99"),
			std::make_pair("C00","F00"),
			std::make_pair("C01","F01"),
			std::make_pair("C02","F02"),
			std::make_pair("C03","F03"),
			std::make_pair("C04","F04"),
			std::make_pair("C05","F05"),
			std::make_pair("C11","F11"),
			std::make_pair("C12","F12"),
			std::make_pair("C15","F15"),
			std::make_pair("C16","F16"),
			std::make_pair("C20","F20"),
			std::make_pair("C26","F26"),
			std::make_pair("C27","F27"),
			std::make_pair("C30","F30"),
			std::make_pair("C36","F36"),
			std::make_pair("C40","F40"),
			std::make_pair("C41","F41"),
			std::make_pair("C42","F42"),
			std::make_pair("C43","F43"),
			std::make_pair("C53","F53"),
			std::make_pair("C54","F54"),
			std::make_pair("C55","F55"),
			std::make_pair("C60","F60"),
			std::make_pair("C61","F61"),
			std::make_pair("C80","F80"),
			std::make_pair("C81","F81"),
			std::make_pair("C83","F83"),
			std::make_pair("C84","F84"),
			std::make_pair("C85","F85"),
			std::make_pair("C86","F86"),
			std::make_pair("C87","F87"),
			std::make_pair("C88","F88"),
			std::make_pair("C89","F89"),
			std::make_pair("C90","F90"),
			std::make_pair("C91","F91"),
			std::make_pair("C92","F92"),
			std::make_pair("C93","F93"),
			std::make_pair("C94","F94"),
			std::make_pair("CYD","DCYD"),
			std::make_pair("CYS","DCYS"),
			std::make_pair("GLN","DGLN"),
			std::make_pair("GLU","DGLU"),
			std::make_pair("HIS","DHIS"),
			std::make_pair("HIS_D","DHIS_D"),
			std::make_pair("HLU","DHLU"),
			std::make_pair("HPR","DHPR"),
			std::make_pair("HTY","DHTY"),
			std::make_pair("ILE","DILE"),
			std::make_pair("LEU","DLEU"),
			std::make_pair("LYS","DLYS"),
			std::make_pair("MAL","DMAL"),
			std::make_pair("MET","DMET"),
			std::make_pair("MPA","DMPA"),
			std::make_pair("MTP","DMTP"),
			std::make_pair("NLU","DNLU"),
			std::make_pair("NVL","DNVL"),
			std::make_pair("PHE","DPHE"),
			std::make_pair("PRO","DPRO"),
			std::make_pair("SER","DSER"),
			std::make_pair("THR","DTHR"),
			std::make_pair("TRP","DTRP"),
			std::make_pair("TYR","DTYR"),
			std::make_pair("VAL","DVAL")
	};

//kdrew: convert data to map
std::map<std::string, std::string> L2DChiralMap(L2DChiralData, L2DChiralData + sizeof L2DChiralData / sizeof L2DChiralData [0]);


std::pair< std::string, std::string > D2LChiralData[] = { 
			std::make_pair("D04","A04"),
			std::make_pair("D05","A05"),
			std::make_pair("D06","A06"),
			std::make_pair("D07","A07"),
			std::make_pair("D12","A12"),
			std::make_pair("D20","A20"),
			std::make_pair("D24","A24"),
			std::make_pair("D30","A30"),
			std::make_pair("D31","A31"),
			std::make_pair("D33","A33"),
			std::make_pair("D34","A34"),
			std::make_pair("D43","A43"),
			std::make_pair("D44","A44"),
			std::make_pair("D45","A45"),
			std::make_pair("D48","A48"),
			std::make_pair("D68","A68"),
			std::make_pair("D69","A69"),
			std::make_pair("D78","A78"),
			std::make_pair("D80","A80"),
			std::make_pair("D82","A82"),
			std::make_pair("D83","A83"),
			std::make_pair("D84","A84"),
			std::make_pair("D91","A91"),
			std::make_pair("D92","A92"),
			std::make_pair("D94","A94"),
			std::make_pair("D96","A96"),
			std::make_pair("D97","A97"),
			std::make_pair("D98","A98"),
			std::make_pair("D99","A99"),
			std::make_pair("DABA","ABA"),
			std::make_pair("DALA","ALA"),
			std::make_pair("DAPA","APA"),
			std::make_pair("DARG","ARG"),
			std::make_pair("DASN","ASN"),
			std::make_pair("DASP","ASP"),
			std::make_pair("E00","B00"),
			std::make_pair("E01","B01"),
			std::make_pair("E02","B02"),
			std::make_pair("E03","B03"),
			std::make_pair("E04","B04"),
			std::make_pair("E05","B05"),
			std::make_pair("E06","B06"),
			std::make_pair("E07","B07"),
			std::make_pair("E12","B12"),
			std::make_pair("E19","B19"),
			std::make_pair("E21","B21"),
			std::make_pair("E27","B27"),
			std::make_pair("E28","B28"),
			std::make_pair("E30","B30"),
			std::make_pair("E31","B31"),
			std::make_pair("E35","B35"),
			std::make_pair("E36","B36"),
			std::make_pair("E38","B38"),
			std::make_pair("E40","B40"),
			std::make_pair("E44","B44"),
			std::make_pair("E47","B47"),
			std::make_pair("E48","B48"),
			std::make_pair("E49","B49"),
			std::make_pair("E50","B50"),
			std::make_pair("E51","B51"),
			std::make_pair("E53","B53"),
			std::make_pair("E54","B54"),
			std::make_pair("E56","B56"),
			std::make_pair("E57","B57"),
			std::make_pair("E58","B58"),
			std::make_pair("E59","B59"),
			std::make_pair("E60","B60"),
			std::make_pair("E61","B61"),
			std::make_pair("E62","B62"),
			std::make_pair("E63","B63"),
			std::make_pair("E67","B67"),
			std::make_pair("E74","B74"),
			std::make_pair("E92","B92"),
			std::make_pair("E93","B93"),
			std::make_pair("E94","B94"),
			std::make_pair("E95","B95"),
			std::make_pair("E96","B96"),
			std::make_pair("E97","B97"),
			std::make_pair("E99","B99"),
			std::make_pair("F00","C00"),
			std::make_pair("F01","C01"),
			std::make_pair("F02","C02"),
			std::make_pair("F03","C03"),
			std::make_pair("F04","C04"),
			std::make_pair("F05","C05"),
			std::make_pair("F11","C11"),
			std::make_pair("F12","C12"),
			std::make_pair("F15","C15"),
			std::make_pair("F16","C16"),
			std::make_pair("F20","C20"),
			std::make_pair("F26","C26"),
			std::make_pair("F27","C27"),
			std::make_pair("F30","C30"),
			std::make_pair("F36","C36"),
			std::make_pair("F40","C40"),
			std::make_pair("F41","C41"),
			std::make_pair("F42","C42"),
			std::make_pair("F43","C43"),
			std::make_pair("F53","C53"),
			std::make_pair("F54","C54"),
			std::make_pair("F55","C55"),
			std::make_pair("F60","C60"),
			std::make_pair("F61","C61"),
			std::make_pair("F80","C80"),
			std::make_pair("F81","C81"),
			std::make_pair("F83","C83"),
			std::make_pair("F84","C84"),
			std::make_pair("F85","C85"),
			std::make_pair("F86","C86"),
			std::make_pair("F87","C87"),
			std::make_pair("F88","C88"),
			std::make_pair("F89","C89"),
			std::make_pair("F90","C90"),
			std::make_pair("F91","C91"),
			std::make_pair("F92","C92"),
			std::make_pair("F93","C93"),
			std::make_pair("F94","C94"),
			std::make_pair("DCYD","CYD"),
			std::make_pair("DCYS","CYS"),
			std::make_pair("DGLN","GLN"),
			std::make_pair("DGLU","GLU"),
			std::make_pair("DHIS","HIS"),
			std::make_pair("DHIS_D","HIS_D"),
			std::make_pair("DHLU","HLU"),
			std::make_pair("DHPR","HPR"),
			std::make_pair("DHTY","HTY"),
			std::make_pair("DILE","ILE"),
			std::make_pair("DLEU","LEU"),
			std::make_pair("DLYS","LYS"),
			std::make_pair("DMAL","MAL"),
			std::make_pair("DMET","MET"),
			std::make_pair("DMPA","MPA"),
			std::make_pair("DMTP","MTP"),
			std::make_pair("DNLU","NLU"),
			std::make_pair("DNVL","NVL"),
			std::make_pair("DPHE","PHE"),
			std::make_pair("DPRO","PRO"),
			std::make_pair("DSER","SER"),
			std::make_pair("DTHR","THR"),
			std::make_pair("DTRP","TRP"),
			std::make_pair("DTYR","TYR"),
			std::make_pair("DVAL","VAL")
	};

//kdrew: convert data to map
std::map<std::string, std::string> D2LChiralMap(D2LChiralData, D2LChiralData + sizeof D2LChiralData / sizeof D2LChiralData [0]);



bool is_d_chiral( core::chemical::ResidueType restype )
{
	std::string const base_name( residue_type_base_name( restype ) );
	std::map< std::string, std::string >::iterator it2 = D2LChiralMap.find( base_name ); 
	TR << "residue base_name: " << base_name << " is_d_chiral: " << (it2 != D2LChiralMap.end() ) << std::endl;
	return ( it2 != D2LChiralMap.end() );
}
bool is_l_chiral( core::chemical::ResidueType restype )
{
	std::string const base_name( residue_type_base_name( restype ) );
	std::map< std::string, std::string >::iterator it2 = L2DChiralMap.find( base_name ); 
	TR << "residue base_name: " << base_name << " is_l_chiral: " << (it2 != L2DChiralMap.end() ) << std::endl;
	return ( it2 != L2DChiralMap.end() );
}


ResidueType const & get_chiral_residue_type( ResidueType const & rt, Chirality chirality )
{
	//kdrew: first letters of a residuetype name (before '_p') are the letter code for the aa and is what is stored in the map
	//std::string base_name;
	//std::string patch_name;
	
	//kdrew: use for finding base_name instead of substr
	std::string const base_name( residue_type_base_name( rt ) );
	std::string const patch_name( residue_type_all_patches_name( rt ) );
	TR << "base_name: " << base_name << " patch_name: " << patch_name << std::endl;
	
	//kdrew: is residuetype patched?
	//Size base_end_pos = rt.name().find("_p");	
	//TR << "base_end_pos: " << base_end_pos << std::endl;
	//if( base_end_pos != std::string::npos )
	//{
	//	base_name = rt.name().substr( 0, base_end_pos );
	//	//kdrew: the remaining porition of the string is the patch, reapply this porition to the end of the mapped string
	//	patch_name = rt.name().substr( base_end_pos, rt.name().size() );
	//}
	////kdrew: if not patched just use the residuetype name
	//else
	//{
	//	base_name = rt.name();
	//	patch_name = "";
	//}

	TR << "restype: " << rt.name() << " " << rt.aa() << std::endl;
	std::string chiral_name;
	
	std::map< std::string, std::string >::iterator it = L2DChiralMap.find( base_name ); 
	if ( it != L2DChiralMap.end() && chirality != L_CHIRALITY )
	{
		chiral_name = it->second;
		TR << "chiral_name: " << chiral_name << std::endl;
		chiral_name.append( patch_name );
		TR << "chiral_name(patched): " << chiral_name << std::endl;
		ResidueTypeSetCOP fa_standard(ChemicalManager::get_instance()->residue_type_set(FA_STANDARD));
		ResidueType const & d_rsd_type( fa_standard->name_map( chiral_name ) );
		return d_rsd_type;
	}
	else
	{
		TR << " not found in L2D map, checking D2L" <<  std::endl;
		std::map< std::string, std::string >::iterator it2 = D2LChiralMap.find( base_name ); 
		if ( it2 != D2LChiralMap.end() && chirality != D_CHIRALITY )
		{
			chiral_name = it2->second;
			TR << "chiral_name: " << chiral_name << std::endl;
			chiral_name.append( patch_name );
			TR << "chiral_name(patched): " << chiral_name << std::endl;
			ResidueTypeSetCOP fa_standard(ChemicalManager::get_instance()->residue_type_set(FA_STANDARD));
			ResidueType const & d_rsd_type( fa_standard->name_map( chiral_name ) );
			return d_rsd_type;
		}
		else
		{
			TR << " not found in D2L" <<  std::endl;
			TR << " possibly achiral (ex GLY) or not listed in map" <<  std::endl;
			return rt;
		}
	}

}

/*
kdrew: the apply function changes a single residue's chirality
pose: pose to make change to
chiral_seq_pos: position of residue to change chirality
*/
void ChiralMover::apply( core::pose::Pose & pose )
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	//kdrew: assert for validity of parameters
	//runtime_assert ( chiral_seq_pos_ != 1 );

	Real phi_angle = pose.phi( chiral_seq_pos_ );
	Real psi_angle = pose.psi( chiral_seq_pos_ );


	TR << "phi_angle: " << phi_angle << " psi_angle: " << psi_angle << std::endl;

	//kdrew: get residue type
	ResidueType rtype = pose.residue_type( chiral_seq_pos_ );
	TR << "Current residue type: " << rtype.name()  << std::endl;
	TR << "Current residue type lower terminus: " << rtype.is_lower_terminus()  << std::endl;
	TR << "Current residue type upper terminus: " << rtype.is_upper_terminus()  << std::endl;

	AtomID const atom1( rtype.atom_index( "N" ), chiral_seq_pos_ );
	AtomID const atom2( rtype.atom_index( "CA" ), chiral_seq_pos_ );
	AtomID const atom3( rtype.atom_index( "C" ), chiral_seq_pos_ );
	AtomID const atom4( rtype.atom_index( "O" ), chiral_seq_pos_ );

	Real pseudo_psi = 0.0;
	if( rtype.is_upper_terminus() )
	{
    	pseudo_psi = pose.conformation().torsion_angle( atom1, atom2, atom3, atom4);
		TR << "pseudo_psi: " << pseudo_psi << std::endl;
	}
	//kdrew: get chiral residue type
	ResidueType const & chiral_rtype = get_chiral_residue_type( rtype, chirality_ );
	TR << "Flipped residue type: " << chiral_rtype.name()  << " " << chiral_rtype.aa() << std::endl;
	if( chiral_rtype.name() == rtype.name() )
	{
		TR << " not making chiral change" << std::endl;
		return;
	}
	conformation::Residue res = pose.residue(chiral_seq_pos_);
	//kdrew: mutate to chiral residue type
	conformation::Residue replace_res ( chiral_rtype, true );

	if(orient_functional_group_)
	{
		AtomIndices rtype_sidechain_atoms = rtype.all_sc_atoms();

		AtomIndices chiral_neighbor_ids = chiral_rtype.bonded_neighbor(chiral_rtype.first_sidechain_atom());
		AtomIndices neighbor_ids = rtype.bonded_neighbor(rtype.first_sidechain_atom());

		utility::vector1< std::pair< std::string, std::string > > atom_pairs;
		
		if( !rtype.atom_type(rtype.first_sidechain_atom()).is_hydrogen() && !chiral_rtype.atom_type(chiral_rtype.first_sidechain_atom()).is_hydrogen()  )
		{
			TR << rtype.name() << " " << rtype.atom_name(rtype.first_sidechain_atom()) << std::endl;
			atom_pairs.push_back( std::make_pair( rtype.atom_name(rtype.first_sidechain_atom()), chiral_rtype.atom_name(chiral_rtype.first_sidechain_atom() ) ));
		}
		else
		{
        	TR << "First sidechain atoms are not heavy atoms " << std::endl;
		}

		//kdrew: loop through all neighboring atoms to first sidechain atom looking for the first two heavy atoms to use for alignment
		for (core::Size j = 1; j <= neighbor_ids.size() && atom_pairs.size() < 3; ++j) 
		{
			//kdrew: if heavy atom
			if( !rtype.atom_type(neighbor_ids[j]).is_hydrogen()) 
			{
				TR << rtype.name() << " " << rtype.atom_name(neighbor_ids[j]) << std::endl;
				//TR << chiral_rtype.name() << " " << chiral_rtype.atom_name(chiral_neighbor_ids[j]) << std::endl;
				
				//kdrew: assumes atom names are the same in both L and D residue types
				atom_pairs.push_back( std::make_pair( rtype.atom_name( neighbor_ids[j] ), rtype.atom_name( neighbor_ids[j] )));
			}
		}
			//if( !chiral_rtype.atom_type(chiral_neighbor_ids[j]).is_hydrogen() )

		//kdrew: if there are not enough side chain atoms (ex. ALA) choose an atom from the backbone (ex. N) for alignment
		if ( atom_pairs.size() < 3 )
		{
			for (core::Size j = 1; j <= neighbor_ids.size() && atom_pairs.size() < 3; ++j) 
			{
				TR << "looping through neighbors: " << rtype.name() << " " << rtype.atom_name(neighbor_ids[j]) << std::endl;

				//kdrew: if heavy atom
				if( !rtype.atom_type(neighbor_ids[j]).is_hydrogen()) 
				{
					AtomIndices neighbor_neighbor_ids = rtype.bonded_neighbor( neighbor_ids[j]  );
					for (core::Size jj = 1; jj <= neighbor_neighbor_ids.size() && atom_pairs.size() < 3; ++jj) 
					{
						if( !rtype.atom_type(neighbor_neighbor_ids[jj]).is_hydrogen()) 
						{
								TR << "looping through neighbor neighbors: " << rtype.name() << " " << rtype.atom_name(neighbor_neighbor_ids[jj]) << std::endl;
								atom_pairs.push_back( std::make_pair( rtype.atom_name( neighbor_neighbor_ids[jj] ), rtype.atom_name( neighbor_neighbor_ids[jj] )));
						}
					}

				}
			}

		}

		//kdrew: assert that residue type has three side chain atoms
		runtime_assert_msg ( atom_pairs.size() == 3 , "not enough heavy atoms to align residues" );

		core::conformation::idealize_hydrogens( replace_res, pose.conformation() );
		replace_res.orient_onto_residue( res, atom_pairs );

		//kdrew: for all sidechain atoms in residue_type, copy xyz
		for (core::Size j = 1; j <= rtype_sidechain_atoms.size(); ++j) 
		{
			replace_res.atom(rtype_sidechain_atoms[j]).xyz( res.atom(rtype_sidechain_atoms[j]).xyz() );
		}
		
		pose.replace_residue( chiral_seq_pos_ , replace_res, false );
	}
	else
	{

		pose.replace_residue( chiral_seq_pos_ , replace_res, true );
		//pose.dump_pdb( "rosetta_out_chiral_preidealized.pdb" );

		//kdrew: idealize alpha carbon hydrogen
		core::conformation::ResidueOP iires = pose.residue( chiral_seq_pos_ ).clone();
		core::conformation::idealize_hydrogens( *iires, pose.conformation() );
		pose.replace_residue( chiral_seq_pos_, *iires, false );
	}

	//pose.dump_pdb( "rosetta_out_chiral_postidealized.pdb" );
	pose.set_phi( chiral_seq_pos_, (-1.0 * phi_angle ) );
	//pose.dump_pdb( "rosetta_out_chiral_moved_phi.pdb" );
	pose.set_psi( chiral_seq_pos_, (-1.0 * psi_angle ) );
	//pose.dump_pdb( "rosetta_out_chiral_moved_phi_psi.pdb" );
	
	 AtomID const atom1c( chiral_rtype.atom_index( "N" ), chiral_seq_pos_ );
	 AtomID const atom2c( chiral_rtype.atom_index( "CA" ), chiral_seq_pos_ );
	 AtomID const atom3c( chiral_rtype.atom_index( "C" ), chiral_seq_pos_ );
	 AtomID const atom4c( chiral_rtype.atom_index( "O" ), chiral_seq_pos_ );

	if( chiral_rtype.is_upper_terminus() )
	{
    	Real chiral_pseudo_psi = pose.conformation().torsion_angle( atom1c, atom2c, atom3c, atom4c);
		TR << "chiral pseudo_psi: " << chiral_pseudo_psi << std::endl;
    	pose.conformation().set_torsion_angle( atom1c, atom2c, atom3c, atom4c, -1.0*pseudo_psi);
	}

	TR << "chiral phi_angle: " << pose.phi( chiral_seq_pos_ ) << " chiral psi_angle: " << pose.psi( chiral_seq_pos_ ) << std::endl;

	TR<< "exiting apply" << std::endl;
}

std::string
ChiralMover::get_name() const {
	return "ChiralMover";
}

///@brief
ChiralMover::ChiralMover( 
		core::Size chiral_seq_position 
	): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( FLIP_CHIRALITY ), orient_functional_group_(false)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover( 
		core::Size chiral_seq_position,
		Chirality chirality 
	): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( chirality ), orient_functional_group_(false)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover( 
		core::Size chiral_seq_position,
		bool orient_functional_group 
	): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( FLIP_CHIRALITY ), orient_functional_group_(orient_functional_group)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::ChiralMover( 
		core::Size chiral_seq_position,
		Chirality chirality,
		bool orient_functional_group 
	): Mover(), chiral_seq_pos_( chiral_seq_position ), chirality_( chirality ), orient_functional_group_(orient_functional_group)
{
	Mover::type( "ChiralMover" );
}

ChiralMover::~ChiralMover(){}

}//chiral
}//simple_moves
}//protocols

