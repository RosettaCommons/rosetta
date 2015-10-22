// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/atomtype_support.cc
/// @author Rocco Moretti (rmorettiase@gmail.com)

#include <core/chemical/atomtype_support.hh>
#include <core/chemical/AtomType.hh>

#include <core/chemical/ResidueType.hh>

namespace core {
namespace chemical {

/// @brief An atom is aromatic if it has any aromatic bonds to a non-virtual atom.
/// TODO: We need better aromatic ring detection.
bool retype_is_aromatic(VD const & atom, ResidueGraph const & graph) {
	OutEdgeIter bonds, bonds_end;
	for ( boost::tie(bonds, bonds_end) = boost::out_edges(atom,graph); bonds != bonds_end; ++bonds ) {
		if ( graph[ *bonds ].bond_name() == AromaticBond ) {
			VD const & tvd( boost::target( *bonds, graph) );
			Atom const & t( graph[tvd] );
			debug_assert( t.element_type() );
			if ( t.element_type()->get_chemical_name() != "Virt" ) {
				return true;
			}
		}
	}
	return false;
}

/// @brief Reassign Rosetta atom types based on the current fullatom heuristics.
///
/// If preserve is true, only retype those atoms which have an atom_type_index of zero.
/// @details The logic here comes from molfile_to_params.py
/// Which is itself based on Rosetta++ ligand_ns.cc set_rosetta_atom_types(),
/// and has been validated against the Meiler and Baker 2006 cross docking test set
/// assignments.
///
/// I'm not saying the logic is good, but it's the logic we're using.
///
/// This function assumes that:
///   * All bonds and atoms exist.
///   * Bond types (bond_name) are correctly set
///   * The appropriate element objects have been set in Atoms.
void
rosetta_retype_fullatom(ResidueType & restype, bool preserve/*=false*/) {
	using namespace core::chemical;

	// For each atom, analyze bonding pattern to determine type
	ResidueGraph const & graph( restype.graph() );
	VDs deferredHs; // Hydrogens which could be typed as Haro, depending on how what they're attached to gets typed.
	VIter itr, itr_end;
	for ( boost::tie(itr, itr_end) = vertices(graph); itr != itr_end; ++itr ) {
		VD vd ( *itr );
		Atom const & a( graph[*itr] );
		if ( preserve && a.atom_type_index() != 0 ) {
			continue;
		}
		ElementCOP element_type( a.element_type() );
		debug_assert( element_type );
		element::Elements element( element_type->element() );
		// H, C, O, N have complicated rules.
		// Everything else maps to a single atom type.
		if ( element_type->get_chemical_name() == "Virt" ) {
			restype.set_atom_type(vd, "VIRT");
		} else if ( element_type->get_chemical_name() == "Repl" ) {
			restype.set_atom_type(vd, "REPL");
		} else if ( element_type->get_chemical_name() == "Suck" ) {
			restype.set_atom_type(vd, "SUCK");
		} else if ( element == element::H ) {
			OutEdgeIter bonds, bonds_end;
			core::Size num_aro_C(0), num_C(0), num_NOS(0);
			for ( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph) );
				Atom const & t( graph[tvd] );
				debug_assert( t.element_type() );
				element::Elements t_element( t.element_type()->element() );
				if ( t_element == element::N || t_element == element::O || t_element == element::S ) { ++num_NOS; }
				// Instead of also counting number of C's typed as aroC's here (which may depend on atom iteration ordering)
				// We annotate the hydrogens as best we can, and notate ones which may switch later.
				if ( t_element == element::C ) {
					++num_C;
					if ( retype_is_aromatic(tvd,graph) ) {
						++num_aro_C;
					}
				}
			}
			if ( num_NOS >=1 ) {
				restype.set_atom_type(vd, "Hpol");
			} else if ( num_aro_C >= 1 ) {
				restype.set_atom_type(vd, "Haro");
			} else {
				restype.set_atom_type(vd, "Hapo");
				if ( num_C >= 1 ) {
					deferredHs.push_back(vd); // Could possibly be Haro - test later.
				}
			}
		} else if ( element == element::C ) {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), num_dbl_nonO(0), num_aro_nonO(0), num_aro_N(0);
			for ( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph) );
				Atom const & t( graph[tvd] );
				debug_assert( t.element_type() );
				if ( t.element_type()->get_chemical_name() == "Virt" ) { continue; }
				element::Elements t_element( t.element_type()->element() );
				switch( graph[*bonds].bond_name() ) {
				case SingleBond :
					if ( t_element == element::H ) { ++num_H; }
					break;
				case DoubleBond :
					saturated = false;
					if ( t_element != element::O ) { ++num_dbl_nonO; }
					break;
				case TripleBond :
					saturated = false;
					break;
				case AromaticBond :
					saturated = false;
					if ( t_element != element::O ) { ++num_aro_nonO; }
					if ( t_element == element::N ) { ++num_aro_N; } // really if, not else if
					break;
				default :
					break;
				}
			}
			if ( saturated ) {
				if ( num_H >= 3 ) {
					restype.set_atom_type(vd, "CH3 ");
				} else if ( num_H == 2 ) {
					restype.set_atom_type(vd, "CH2 ");
				} else {
					restype.set_atom_type(vd, "CH1 ");
				}
			} else { // unsaturated
				if ( num_aro_nonO >= 2 ) {
					restype.set_atom_type(vd, "aroC");
				} else if ( num_dbl_nonO >= 1 ) {
					restype.set_atom_type(vd, "aroC");
				} else if ( num_aro_N >= 1 ) {
					restype.set_atom_type(vd, "CNH2");
				} else {
					restype.set_atom_type(vd, "COO ");
				}
			}
		} else if ( element == element::N ) {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), heavy_nbrs(0);
			for ( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph) );
				Atom const & t( graph[tvd] );
				debug_assert( t.element_type() );
				if ( t.element_type()->get_chemical_name() == "Virt" ) { continue; }
				element::Elements t_element( t.element_type()->element() );
				if ( t_element == element::H ) { ++num_H; }
				else { ++heavy_nbrs; } // We've already ignored all the virtual atoms.
				if ( graph[*bonds].bond_name() != SingleBond ) { saturated = false; }
			}

			if ( num_H >= 3 ) {
				restype.set_atom_type(vd, "Nlys"); // carries a VERY high desolvation penalty
			} else if ( num_H == 2 ) {
				// Not totally sure about this one, may want Ntrp instead if more than one heavy neighbor:
				restype.set_atom_type(vd, "NH2O"); // Narg would also be a possibility, but they're fairly similar
			} else if ( num_H == 1 ) {
				if ( heavy_nbrs <= 2 ) {
					restype.set_atom_type(vd, "Ntrp"); // should always be 2 neighbors, not less
				} else {
					restype.set_atom_type(vd, "Ntrp"); // Npro? protonated tertiary amine
				} // I know they're the same -- I'm just copying molfile_to_params, which splits the case.
			} else {
				if ( heavy_nbrs <= 2 ) {
					restype.set_atom_type(vd, "Nhis");
				} else if ( heavy_nbrs == 3 ) {
					if ( saturated ) {
						restype.set_atom_type(vd, "Nhis"); // deprotonated tertiary amine; need an sp3 hybrid H-bond acceptor type...
					} else { // This also catches nitro groups -- is that what we want here?
						restype.set_atom_type(vd, "Npro"); // X=[N+](X)X, including nitro groups
					}
				} else {
					restype.set_atom_type(vd, "Npro"); // quaternary amine
				}
			}
		} else if ( element == element::O ) {
			OutEdgeIter bonds, bonds_end;
			bool saturated(true);
			core::Size num_H(0), num_bonds(0), bonded_to_N(0), bonded_to_C_to_N(0), unsat_nbrs(0);
			for ( boost::tie(bonds, bonds_end) = boost::out_edges(*itr,graph); bonds != bonds_end; ++bonds ) {
				VD const & tvd( boost::target( *bonds, graph) );
				Atom const & t( graph[tvd] );
				debug_assert( t.element_type() );
				if ( t.element_type()->get_chemical_name() == "Virt" ) { continue; }
				element::Elements t_element( t.element_type()->element() );

				++num_bonds; // Bonds to non-virtual atoms.
				if ( graph[*bonds].bond_name() != SingleBond ) { saturated = false; }
				if ( t_element == element::H ) { ++num_H; }
				else if ( t_element == element::N ) { ++bonded_to_N; }
				OutEdgeIter bonds2, bonds_end2; // second degree bonds.
				bool sat_neighbor = true;
				for ( boost::tie(bonds2, bonds_end2) = boost::out_edges(tvd,graph); bonds2 != bonds_end2; ++bonds2 ) {
					// Ignore the bond back to the atom we're typing.
					//if( boost::target( *bonds2, graph) == *itr ) { continue; }
					VD const & tvd2( boost::target( *bonds2, graph) );
					Atom const & t2( graph[tvd2] );
					debug_assert( t2.element_type() );
					if ( t2.element_type()->get_chemical_name() == "Virt" ) { continue; }
					element::Elements t2_element( t2.element_type()->element() );

					if ( t_element == element::C && t2_element == element::N ) { ++bonded_to_C_to_N; }
					if ( graph[*bonds2].bond_name() != SingleBond ) { sat_neighbor = false; }
				}
				if ( ! sat_neighbor ) { ++unsat_nbrs; }
			}
			if ( saturated ) {
				if ( num_bonds < 2 ) {
					restype.set_atom_type(vd, "OOC "); // catches C(=O)[O-] (Kekule form) -- new rule by IWD
				} else {
					core::Size ring_size( restype.smallest_ring_size( *itr ) );
					if ( num_H > 0 ) {
						restype.set_atom_type(vd, "OH  "); // catches C(=O)OH (Kekule form)
					} else if ( ring_size < 5 ) {
						restype.set_atom_type(vd, "OH  "); // small, strained rings leave the O more exposed? (IWD, see 1p8d)
					} else if ( ring_size < 999999 && unsat_nbrs > 0 ) {
						restype.set_atom_type(vd, "Oaro"); // catches aromatic O in furan-like rings, though I rarely see these H-bond (IWD)
					} else {
						restype.set_atom_type(vd, "OH  "); // catches ethers, ROR (IWD, see comment)
						// The lone pairs on ethers are capable of H-bonding in the same way that alcohols are.
						// While alkyl ethers are quite non-polar, many others seem to make Hbonds,
						// such as those attached to phosphates (R-O-PO3), methyls (R-O-CH3), and aromatic rings (R-O-Ph).
						// It is unclear from the literature how strong these are, and is probably very situation dependent.
					}
				}
			} else if ( num_H > 0 ) {
				restype.set_atom_type(vd, "OH  "); // catches c(o)oH (aromatic bonds to both O)
			} else if ( bonded_to_N ) {
				restype.set_atom_type(vd, "ONH2");
			} else if ( bonded_to_C_to_N ) { // This is a non-standard rule introduced by IWD, agreed to by KWK:
				restype.set_atom_type(vd, "ONH2");
			} else {
				restype.set_atom_type(vd, "OOC ");
			}
		} else if ( element::S  == element ) {
			restype.set_atom_type(vd, "S   ");
		} else if ( element::P  == element ) {
			restype.set_atom_type(vd, "Phos");
		} else if ( element::F  == element ) {
			restype.set_atom_type(vd, "F   ");
		} else if ( element::Cl == element ) {
			restype.set_atom_type(vd, "Cl  ");
		} else if ( element::Br == element ) {
			restype.set_atom_type(vd, "Br  ");
		} else if ( element::I  == element ) {
			restype.set_atom_type(vd, "I   ");
		} else if ( element::Na == element ) {
			restype.set_atom_type(vd, "Na1p");
		} else if ( element::K  == element ) {
			restype.set_atom_type(vd, "K1p ");
		} else if ( element::Mg == element ) {
			restype.set_atom_type(vd, "Mg2p");
		} else if ( element::Fe == element ) {
			restype.set_atom_type(vd, "Fe3p");
		} else if ( element::Ca == element ) {
			restype.set_atom_type(vd, "Ca2p");
		} else if ( element::Zn == element ) {
			restype.set_atom_type(vd, "Zn2p");
		} else {
			utility_exit_with_message("Cannot type atom with element '"+element::name_from_elements(element)+"'");
		}
	} // For vertices in graph

	// Hydrogens attached to aroCs == Haro
	// We only look at the deferred atoms, which will be hydrogens attached to carbons, already typed as Hapo.
	// As it's only the deferred atoms, we don't need to re-check the preserve state.
	for ( VDs::const_iterator aroit(deferredHs.begin()), aroend(deferredHs.end()); aroit != aroend; ++aroit ) {
		OutEdgeIter bonds, bonds_end;
		for ( boost::tie(bonds, bonds_end) = boost::out_edges(*aroit,graph); bonds != bonds_end; ++bonds ) {
			VD const & tvd( boost::target( *bonds, graph) );
			Atom const & t( graph[tvd] );
			if ( t.atom_type_index() && restype.atom_type_set()[ t.atom_type_index() ].atom_type_name() == "aroC" ) {
				restype.set_atom_type( *aroit, "Haro" );
				break;
			}
		}
	}

#ifndef NDEBUG
	// Make sure all the types have been set.
	for ( boost::tie(itr, itr_end) = vertices(graph); itr != itr_end; ++itr ) {
		Atom const & a( graph[*itr] );
		debug_assert( a.atom_type_index() != 0 );
	}
#endif

}

} // chemical
} // core
