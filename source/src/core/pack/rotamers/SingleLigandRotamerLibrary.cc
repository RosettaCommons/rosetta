// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/rotamers/SingleLigandRotamerLibrary.cc
///
/// @brief
/// @author Ian W. Davis
#include <utility/fixedsizearray1.hh>
// Unit headers
#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>

// Package headers
#include <core/pack/dunbrack/ChiSet.hh>
#include <core/pack/dunbrack/DunbrackRotamer.hh>
#include <core/pack/dunbrack/RotamerLibraryScratchSpace.hh>

// Project headers
//#include <core/chemical/automorphism.hh>
//#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Residue.functions.hh>

#include <core/io/pdb/pdb_writer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Utility headers
//#include <numeric/xyz.functions.hh>
//#include <numeric/xyzMatrix.hh>
//#include <numeric/model_quality/rms.hh>
//#include <ObjexxFCL/FArray1D.hh>
//#include <ObjexxFCL/FArray2D.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>


// C++ headers
#include <cstdlib>
#include <fstream>
#include <string>
#include <set>

#include <utility/vector1.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>


namespace core {
namespace pack {
namespace rotamers {


static THREAD_LOCAL basic::Tracer TR( "core.pack.rotamers.SingleLigandRotamerLibrary" );


// helper for debugging
void dump_library(std::string filename, RotamerVector const & rotamers)
{
	std::ofstream out( filename.c_str() );
	for ( Size i = 1; i <= rotamers.size(); ++i ) {
		out << "MODEL \n";
		Size atomno = 1;
		core::io::pdb::dump_pdb_residue(*rotamers[i], atomno, out);
		out << "ENDMDL\n";
	}
	out.close();
}


SingleLigandRotamerLibrary::SingleLigandRotamerLibrary():
	SingleResidueRotamerLibrary(),
	atom_positions_(),
	ref_energy_(0.0)
	// rigid_frags_(),
	// automorphs_(),
	// frag_automorphs_(),
	// total_superpos_(0)
{}

SingleLigandRotamerLibrary::~SingleLigandRotamerLibrary()
{}

/// @details Reads conformers from PDB-format file.
/// Chain ID, residue name and number, etc are all ignored -- must have TER records.
void
SingleLigandRotamerLibrary::init_from_file(
	std::string const & filename,
	chemical::ResidueType const & restype
)
{
	//std::cout << "Loading from: " << filename << "\n";
	utility::io::izstream data( filename.c_str() );
	if ( !data.good() ) {
		utility_exit_with_message( "Unable to open file: " + filename + '\n' );
	}
	atom_positions_.clear();

	NamePosMap name_map;
	bool found_ref_energy = false;
	std::string line;
	// This code is not currently as smart as the general-purpose PDB reader.
	// Any atoms in the residue that don't have coordinate entries will be
	// left with their default values, leading to really weird bugs.
	// We can do a limited building from ideal coordinates, for hydrogens and virtual atoms.
	while ( std::getline( (std::istream&)data, line) ) {
		if ( utility::startswith(line, "ATOM  ") || utility::startswith(line, "HETATM") ) {
			if ( line.length() < 54 ) {
				TR << "ATOM/HETATM line too short in PDB-format rotamer file!" << std::endl;
				continue; // to next line
			}
			std::string atom_name = line.substr(12,4);
			core::Real x, y, z;
			x = std::atof( line.substr(30,8).c_str() );
			y = std::atof( line.substr(38,8).c_str() );
			z = std::atof( line.substr(46,8).c_str() );
			//std::cout << x << " " << y << " " << z << "\n";

			if ( name_map.count( atom_name ) != 0 ) {
				//TODO: cache so we only print this once.
				TR.Warning << "ATOM name " << atom_name << " found more than once in rotamer - using later position. " << std::endl;
			}
			name_map[ atom_name ] = core::Vector( x, y, z );

		} else if ( utility::startswith(line, "REF_EN") ) {
			if ( found_ref_energy ) {
				TR.Error << "Reference energy specified more than once in PDB-format rotamer file!" << std::endl;
			}
			found_ref_energy = true;
			ref_energy_ = std::atof( line.substr(6).c_str() );
			TR << "Reference energy for " << restype.name() << " is " << ref_energy_ << std::endl;
		} else { // e.g. TER lines
			if ( name_map.size() ) {
				atom_positions_.push_back( name_map );
			}
			name_map.clear();
		}
	}
	// Catch the last entry if we don't end with a TER
	if ( name_map.size() ) {
		atom_positions_.push_back( name_map );
	}

	TR << "Read in " << atom_positions_.size() << " rotamers from " << filename << " !" << std::endl;
	data.close();

	// Breaking the ligand into rigid fragments that would supply (putative) pharamacophores
	// to superimpose on was a nice idea, but it breaks the packer assumption that nbr_atom doesn't move.

	//find_fragments(restype);
	//list_automorphisms(restype);
	//unique_auto_for_frags();
}

/// @brief Reads conformers from a vector of name:coordinate maps
void
SingleLigandRotamerLibrary::init_from_vector(
	utility::vector1< NamePosMap > const & coordinates
) {
	atom_positions_.append( coordinates );
}

/// @details Not currently implemented -- returns 0.
Real
SingleLigandRotamerLibrary::rotamer_energy_deriv(
	conformation::Residue const & rsd,
	dunbrack::RotamerLibraryScratchSpace & scratch
) const
{
	dunbrack::Real5 & dE_dbb( scratch.dE_dbb() );
	dunbrack::Real4 & dE_dchi( scratch.dE_dchi() );
	std::fill( dE_dbb.begin(), dE_dbb.end(), 0 );
	std::fill( dE_dchi.begin(), dE_dchi.end(), 0 );
	return rotamer_energy(rsd, scratch);
}


/// @details Not currently implemented -- returns 0.
/// For most things, it's actually quite hard to correctly rank the conformers by energy.
/// Thus, I've chosen not to introduce noise by attempting to rank them and doing it badly.
Real
SingleLigandRotamerLibrary::rotamer_energy(
	conformation::Residue const & /*rsd*/,
	dunbrack::RotamerLibraryScratchSpace & /*scratch*/
) const
{
	return ref_energy_; //0.0;
}


/// @details Not currently implemented -- returns 0.
/// For most things, it's actually quite hard to correctly rank the conformers by energy.
/// Thus, I've chosen not to introduce noise by attempting to rank them and doing it badly.
Real
SingleLigandRotamerLibrary::best_rotamer_energy(
	conformation::Residue const & /*rsd*/,
	bool /*curr_rotamer_only*/,
	dunbrack::RotamerLibraryScratchSpace & /*scratch*/
) const
{
	return ref_energy_; //0.0;
}


/// Helper function for superposition
//void SingleLigandRotamerLibrary::superimpose(
// conformation::Residue const & existing,
// conformation::Residue & conformer,
// Fragment const & frag,
// Automorphism const & morph
//) const
//{
// //utility_exit_with_message("not implemented yet");
// using ObjexxFCL::FArray1D;
// using ObjexxFCL::FArray2D;
// using namespace numeric;
//
// Size const natoms = frag.size();
// FArray2D< Real > xx( 3, natoms, 0. );
// FArray2D< Real > yy( 3, natoms, 0. );
// FArray1D< Real > ww( natoms, 1.0 ); // uniform weighting
// FArray2D< Real > uu( 3, 3, 0.0 );
// Real ctx(0); // not really used
//
// Vector e_ctr(0.), c_ctr(0.);
// for(Size i = 1; i <= natoms; ++i) {
//  Size const e_atom = frag[i];
//  Size const c_atom = morph[e_atom];
//  for(Size j = 1; j <= 3; ++j) {
//   xx(j, i) =  existing.xyz(e_atom)(j);
//   yy(j, i) = conformer.xyz(c_atom)(j);
//  }
//  e_ctr +=  existing.xyz(e_atom);
//  c_ctr += conformer.xyz(c_atom);
// }
//debug_assert(natoms > 0);
// e_ctr /= natoms;
// c_ctr /= natoms;
//
// // This is not actually very accurate, in my experience so far!
// numeric::model_quality::findUU( xx, yy, ww, natoms, uu, ctx );
//
// typedef xyzMatrix< Real > Rotation;
// Rotation rot(Rotation::rows( uu(1,1), uu(1,2), uu(1,3), uu(2,1), uu(2,2), uu(2,3), uu(3,1), uu(3,2), uu(3,3) ));
// for(Size i = 1, i_end = existing.natoms(); i <= i_end; ++i) {
//  conformer.set_xyz( i, (rot * (conformer.xyz(i) - c_ctr)) + e_ctr );
// }
//}


/// @brief Helper function, combines existing's metadata with conformer's conformation.
conformation::ResidueOP
dup_residue(
	conformation::Residue const & existing,
	conformation::Residue const & conformer
)
{
	// This is bad:  fields like seqpos, chain, etc. don't match existing residue!
	//conformation::ResidueOP newrsd = rotamers_[i]->clone();

	// Could start by cloning either one, but I think people are more likely to introduce
	// new metadata than new conformational data, so I'll let clone() copy the metadata.
	//conformation::ResidueOP newrsd = existing.clone();
	//newrsd->atoms() = conformer.atoms();
	//newrsd->chi() = conformer.chi();
	//newrsd->mainchain_torsions() = conformer.mainchain_torsions();
	//newrsd->actcoord() = conformer.actcoord();

	// The above is also bad:  existing may not be the same residue type as conformer!
	conformation::ResidueOP newrsd = conformer.clone();
	newrsd->chain( existing.chain() );
	newrsd->seqpos( existing.seqpos() );
	newrsd->copy_residue_connections_from( existing ); // this is probably not good enough if residue types diverge more than protonation state...

	return newrsd;
}


void
SingleLigandRotamerLibrary::fill_rotamer_vector(
	pose::Pose const & pose,
	scoring::ScoreFunction const &,
	pack::task::PackerTask const & task,
	graph::GraphCOP,
	chemical::ResidueTypeCOP concrete_residue,
	conformation::Residue const& existing_residue,
	utility::vector1< utility::vector1< Real > > const & /*extra_chi_steps*/,
	bool buried,
	RotamerVector & rotamers //utility::vector1< conformation::ResidueOP >
) const
{
	//std::cout << "SingleLigandRotamerLibrary :: fill_rotamer_vector() being called...\n";
	//rotamers.clear(); // am I supposed to do this?  No: might contain rotamers for other residue types already!
	int const start_size = rotamers.size();

	// Build the base set of rotamer from PDB data

	RotamerVector base_rotamers;
	build_base_rotamers( *concrete_residue, base_rotamers );

	// Expand base rotamers with extra proton chi sampling

	bool expand_proton_chi = ( concrete_residue->n_proton_chi() != 0 );
	Size const max_total_rotamers = 21654; // = 401 rotamers * 54 hydroxyl variations = 401 * (2 * 3^3)
	RotamerVector new_rotamers;

	if ( basic::options::option[ basic::options::OptionKeys::packing::ignore_ligand_chi]() == true ) {
		expand_proton_chi = false;
	}

	// Logic for creating proton_chi rotamers copied from APL (RotamerSet_.cc)
	utility::vector1< pack::dunbrack::ChiSetOP > proton_chi_chisets;
	if ( expand_proton_chi ) {
		proton_chi_chisets.push_back( dunbrack::ChiSetOP( new pack::dunbrack::ChiSet( concrete_residue->nchi() ) ) );
		for ( Size ii = 1; ii <= concrete_residue->n_proton_chi(); ++ii ) {
			pack::dunbrack::expand_proton_chi(
				task.residue_task( existing_residue.seqpos() ).extrachi_sample_level(
				buried,
				concrete_residue->proton_chi_2_chi( ii ),
				*concrete_residue ),
				concrete_residue,
				ii, proton_chi_chisets);
			// In a pathological case, I've seen 30 rotamers * 20,000 proton chi variations = 600,000 rotamers = out of memory
			// That's 9, count 'em 9, hydroxyls in the ligand for PDB 1u33.
			// Wait, wait -- I can do better -- 19 hydroxyls in PDB 1xd1.
			if ( base_rotamers.size()*proton_chi_chisets.size() > max_total_rotamers ) {
				TR.Warning << "Aborting proton_chi expansion for " << concrete_residue->name() << " because we would exceed " << max_total_rotamers << " rotamers!" << std::endl;
				proton_chi_chisets.resize( max_total_rotamers / base_rotamers.size() );
				break;
			}
		}
		new_rotamers.reserve( base_rotamers.size()*proton_chi_chisets.size() );
	} else {
		new_rotamers.reserve( base_rotamers.size() );
	}


	// Fill new_rotamers with new Residues, including proton_chi expansions
	for ( Size i = 1; i <= base_rotamers.size(); ++i ) {
		debug_assert( concrete_residue->name() == base_rotamers[i]->name() );
		if ( concrete_residue->in_residue_type_set() ) {
			debug_assert( concrete_residue->residue_type_set()->name() == base_rotamers[i]->residue_type_set()->name() ); // fa_standard / centroid
		}
		if ( expand_proton_chi ) {
			for ( Size ii = 1; ii <= proton_chi_chisets.size(); ++ii ) {
				conformation::ResidueOP newrsd = dup_residue( existing_residue, *base_rotamers[i] );
				new_rotamers.push_back( newrsd );
				for ( Size jj = 1; jj <= concrete_residue->n_proton_chi(); ++jj ) {
					newrsd->set_chi(
						concrete_residue->proton_chi_2_chi( jj ),
						proton_chi_chisets[ ii ]->chi[ jj ] );
				}
			}
		} else {
			// This is bad:  fields like seqpos, chain, etc. don't match existing residue!
			//conformation::ResidueOP newrsd = base_rotamers[i]->clone();
			conformation::ResidueOP newrsd = dup_residue( existing_residue, *base_rotamers[i] );
			new_rotamers.push_back( newrsd );
		}
	}

	// Superimpose
	// This was a nice idea, superimposing on the "pharmacophores" of the ligand,
	// and just trying them all.  But no doubt it violates the packer's expectation
	// that the residue's NBR_ATOM does not move during repacking!
	// So, this code is not sound during scoring -- DO NOT USE.
	//bool const do_multiple_superpos = false;
	//if( do_multiple_superpos ) {
	// using utility::vector1;
	// rotamers.reserve( rotamers.size() + new_rotamers.size()*total_superpos_ );
	// // For all rigid fragments...
	// for(Size i = 1, i_end = rigid_frags_.size(); i <= i_end; ++i) {
	//  Fragment const & frag = rigid_frags_[i];
	//  // And for all of their automorphisms...
	//  vector1< Automorphism * > const & morphs = frag_automorphs_[i];
	//  for(Size j = 1, j_end = morphs.size(); j <= j_end; ++j) {
	//   Automorphism const & morph = *(morphs[j]);
	//   // Try superimposing each conformer using that grouping of atoms!
	//   for(Size k = 1, k_end = new_rotamers.size(); k <= k_end; ++k) {
	//    conformation::ResidueOP newrsd = dup_residue( existing_residue, *new_rotamers[k] );
	//    superimpose( existing_residue, *newrsd, frag, morph );
	//    rotamers.push_back(newrsd);
	//   }
	//  }
	// }
	//} else {
	rotamers.reserve( rotamers.size() + new_rotamers.size() );
	for ( Size k = 1, k_end = new_rotamers.size(); k <= k_end; ++k ) {
		conformation::ResidueOP newrsd = new_rotamers[k];
		// Superimposes on nbr_atom and 2 of its neighbors
		newrsd->place( existing_residue, pose.conformation() );
		rotamers.push_back(newrsd);
	}
	//}

	int const end_size = rotamers.size();
	TR << "Added " << end_size - start_size << " rotamers for " << concrete_residue->name() << std::endl;

	// Debugging:
	//dump_library(concrete_residue->name()+".rotlib.pdb", rotamers);
}

/// @brief Build a set of rotamers for the given ResidueType
void
SingleLigandRotamerLibrary::build_base_rotamers( chemical::ResidueType const & restype, RotamerVector & base_rotamers ) const {

	// Why do we do this on the fly, instead of when loading?
	// Because the same SingleLigandRotamerLibrary might be applied to more than one ResidueType,
	// which means the Residues being built would be different.
	// (Storing Residue objects also means that references to orphan ResidueTypes hang around after they're destroyed.)

	base_rotamers.reserve( base_rotamers.size() + atom_positions_.size() ); // Adding rotamers

	// Any atoms in the residue that don't have coordinate entries will be
	// left with their default values, leading to really weird bugs.
	// We can do a limited building from ideal coordinates, for hydrogens and virtual atoms.
	std::set< std::string > skipped_atom_names;
	utility::vector1< bool > missed(restype.natoms(),false); // Don't reset - only notify once, instead of for each library entry
	for ( core::Size resn(1); resn <= atom_positions_.size(); ++resn ) {
		NamePosMap const & name_map( atom_positions_[resn] );
		conformation::ResidueOP rsd = conformation::ResidueFactory::create_residue( restype );
		core::Size set_xyzs = 0;
		utility::vector1< bool > missing(rsd->natoms(),true);
		for ( NamePosMap::const_iterator iter( name_map.begin() ), iter_end( name_map.end() ); iter != iter_end; ++iter ) {
			std::string const & atom_name( iter->first );
			core::Vector const & pos( iter->second );
			if ( rsd->has( atom_name ) ) {
				rsd->set_xyz( atom_name, pos );
				missing[ rsd->atom_index(atom_name) ] = false;
				set_xyzs += 1;
			} else if ( skipped_atom_names.count(atom_name) == 0 ) {
				TR.Warning << "Skipping unrecognized atom '" << atom_name << "' in library for " << restype.name() << std::endl;
				skipped_atom_names.insert(atom_name);
			}
		}
		if ( set_xyzs < rsd->natoms() ) { fill_missing_atoms( missing, rsd, missed ); }
		// Torsion angles are not automatically calculated from the coordinates.
		// We should manually assign chi() and mainchain_torsions() for each conformer.
		conformation::set_chi_according_to_coordinates( *rsd );
		base_rotamers.push_back( rsd );
	}

}

/// @brief Fills in missing hydrogens/virtual atoms from library load
void
SingleLigandRotamerLibrary::fill_missing_atoms( utility::vector1< bool > missing, conformation::ResidueOP rsd, utility::vector1< bool > & missed ) const
{
	debug_assert( rsd );
	debug_assert( missing.size() == rsd->natoms() );
	//Unlike Residue::fill_missing_atoms(), only do a single pass -
	// The residue should be constructed so that any atoms which would be typically missing
	// (i.e. hydrogens and virtual atoms) are either built from present atoms, or can be built
	// from "missing" atoms which are built earlier
	for ( Size i=1; i<= rsd->natoms(); ++i ) {
		if ( missing[i] ) {
			if ( ! rsd->atom_is_hydrogen( i ) && ! rsd->is_virtual( i ) ) {
				utility_exit_with_message("Non-virtual heavy atom "+rsd->atom_name(i)+" is missing in rotamer library for residue "+rsd->name()+"!");
			}

			chemical::AtomICoor const & ic( rsd->icoor(i) );
			// check to see if any of our stub atoms are missing:
			for ( Size j=1; j<= 3; ++j ) {
				Size stubno( ic.stub_atom(j).atomno() );
				if ( missing[ stubno ] ) {
					TR.Error << "[ ERROR ] Missing atom " << stubno << " (" << rsd->atom_name(stubno) << ") when trying to place atom " <<
						i << " (" << rsd->atom_name(i) << ") in " << rsd->name() << std::endl;
					utility_exit_with_message("Cannot build missing atoms in ligand rotamer library");
				}
			}

			// no stub atoms missing: build our ideal coordinates
			missing[i] = false; // In case we're building later residues off of this one.
			rsd->set_xyz( i, ic.build( *rsd ) ); // We just checked that all stub atoms exist in this residue, so we should be safe with the build call.
			if ( ! missed[ i ] ) {
				missed[ i ] = true;
				TR << "Atom " << rsd->atom_name(i) << " from residue " << rsd->name() << " not found in PDB_ROTAMERS library, creating based on idealized geometry." << std::endl;
			}
		}
	}
}


/// @details Not implemented -- will cause program termination.
/// Is this only used by coarse representations?
void
SingleLigandRotamerLibrary::write_to_file( utility::io::ozstream & /*out*/ ) const
{
	utility_exit_with_message("Ligand residue rotamers can't be written to file!");
}


// Helper function
//bool bond_is_rotatable(chemical::ResidueTypeCOP restype, core::Size a1, core::Size a2)
//{
// for(core::Size i = 1, i_end = restype->nchi(); i <= i_end; ++i) {
//  chemical::AtomIndices chi = restype->chi_atoms(i);
// debug_assert( chi.size() == 4 );
//  if( (chi[2] == a1 && chi[3] == a2) || (chi[2] == a2 && chi[3] == a1) ) return true;
// }
// return false;
//}


//void SingleLigandRotamerLibrary::find_fragments(chemical::ResidueTypeCOP restype)
//{
// // Fragments are delimited by rotatable bonds, but the atom on the far end
// // of the bond is still part of the fragment.
// // Thus, fragments may overlap slightly.
// using core::Size;
// using utility::vector1;
// using namespace core::chemical;
//
// Size const natoms = restype->nheavyatoms();
// vector1<bool> in_frag_core(natoms, false); // atoms on far side of rot bond aren't in "core"
// for(Size root = 1; root <= natoms; ++root) {
//  if( in_frag_core[root] ) continue; // already got this fragment!
//  Fragment the_frag;
//  vector1<bool> visited(natoms, false);
//  AtomIndices to_visit;
//  to_visit.push_back(root); // will only hold atoms in fragment core
//  while( !to_visit.empty() ) {
//   Size const curr = to_visit.back();
//   to_visit.pop_back();
//   if( visited[curr] ) continue;
//   visited[curr] = true;
//   the_frag.push_back(curr);
//   in_frag_core[curr] = true;
//   AtomIndices const & nbrs = restype->nbrs(curr);
//   for(Size i = 1; i <= nbrs.size(); ++i) {
//    Size const nbr = nbrs[i];
//    if( nbr > natoms ) continue; // e.g. hydrogens!
//    if( visited[nbr] ) continue;
//    if( bond_is_rotatable(restype, curr, nbr) ) {
//     // don't "recurse", handle everything here
//     visited[nbr] = true;
//     the_frag.push_back(nbr);
//    } else { // non-rotatable bond, part of fragment core
//     // "recurse"; flags will be set in main loop
//     to_visit.push_back(nbr);
//    }
//   }
//  }
//  // Only keep frags with 3+ atoms -- can't superimpose on less.
//  if( the_frag.size() >= 3 ) rigid_frags_.push_back(the_frag);
// }
//debug_assert( rigid_frags_.size() <= restype->nchi()+1 );
//
// for(Size i = 1; i <= rigid_frags_.size(); ++i) {
//  TR << "Fragment " << i << ":";
//  Fragment const & frag = rigid_frags_[i];
//  for(Size j = 1; j <= frag.size(); ++j) TR << " " << restype->atom_name(frag[j]);
//  TR << std::endl;
// }
//}


//void SingleLigandRotamerLibrary::list_automorphisms(chemical::ResidueTypeCOP restype)
//{
// using namespace core::chemical;
// AutomorphismIterator ai( restype, false /*don't include H*/ );
// while(true) {
//  Automorphism a = ai.next();
//  if( a.empty() ) break;
//  automorphs_.push_back( a );
// }
//debug_assert( automorphs_.size() > 0 );
// TR << "Ligand has " << automorphs_.size() << " automorphisms" << std::endl;
//}


//void SingleLigandRotamerLibrary::unique_auto_for_frags()
//{
// using utility::vector1;
// total_superpos_ = 0;
// frag_automorphs_.resize( rigid_frags_.size() ); // one entry per fragment
// // For each fragment...
// for(Size i = 1, i_end = rigid_frags_.size(); i <= i_end; ++i) {
//  Fragment const & frag = rigid_frags_[i];
//  vector1< Automorphism * > const & frag_morphs = frag_automorphs_[i];
//  // Look at each whole-molecule automorphism...
//  for(Size j = 1, j_end = automorphs_.size(); j <= j_end; ++j) {
//   Automorphism /*const*/ & new_morph = automorphs_[j];
//   bool already_have_it = false;
//   // And add it if it's not equivalent to some other one we already added.
//   for(Size k = 1, k_end = frag_morphs.size(); k <= k_end && !already_have_it; ++k) {
//    Automorphism const & old_morph = *(frag_morphs[k]);
//   debug_assert( new_morph.size() == old_morph.size() );
//    bool are_same = true;
//    // Two automorphisms are the same from a fragment's point of view
//    // if they contain the same mapping for all fragment atom positions.
//    for(Size l = 1, l_end = frag.size(); l <= l_end; ++l) {
//     Size const atom = frag[l];
//     if( new_morph[atom] != old_morph[atom] ) {
//      are_same = false;
//      break;
//     }
//    } // end compare two automorphisms
//    if( are_same ) already_have_it = true;
//   } // end search over existing automorphisms for fragment
//   if( !already_have_it ) {
//    frag_automorphs_[i].push_back( &new_morph );
//    total_superpos_ += 1;
//   }
//  } // end search over whole-molecule automorphisms
// } // end search over fragments
//debug_assert( total_superpos_ > 0 );
// TR << total_superpos_ << " unique possible rotamer-substitution superpositions" << std::endl;
//}


} // namespace rotamers
} // namespace scoring
} // namespace core
