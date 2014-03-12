// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
/// @file ResidueGraphTypes.hh
///
/// @brief
/// Graph structure for ResidueType
///
/// @details
/// This is the typedefs and filtered graphs for the graph implementation of ResidueType based on boost graphs.
/// Filtered graphs are graph structures that have been filtered based on a certain criteria. For example, the
/// Acceptor atom graph has been filtered so that every node and edge in the graph is associated with an acceptor
/// atom. The properties of the filtered graphs can be determined by any criteria. Currently, atom types are used
/// as the metric to filter the graphs. This does not have to be the case. Graphs can be filtered based on the
/// atoms, orbitals, etc etc. It is up to your immagination. The unit tests for these show examples of how to use
/// the filtered graphs.
///
/// Each filter graph has an operator that is used to determine if a node should be in the graph. An iterator through
/// each node and edge of the graph is available. Specifically, if you want to iterate through the graph nodes, you would
/// use this method: for(HeavyAtomVIterPair vp = boost::vertices(heavy_atom_graph); vp.first != vp.second; ++vp.first){}
///
/// @author Steven Combs
////////////////////////////////////////////////////////////////////////
// Unit headers
#include <core/chemical/ResidueGraphTypes.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/chemical/AtomType.hh>

// Package headers

namespace core {
    namespace chemical {
        
        /////////////////////////////////////////////////////////////
        ////////// PREDICATES for FILTERED GRAPHS ///////////////////
        ////////////////////////////////////////////////////////////
        


        bool HeavyAtomFilter::operator()(VD const vd) const{
            return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_heavyatom();
        }
        
        bool AcceptorAtomFilter::operator()(VD const vd) const{
            return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_acceptor();
        }

        bool HeavyAtomWithPolarHydrogensFilter::operator()(VD const vd) const{
            
            for(OutEdgeIterPair ep = boost::out_edges(vd, *graph_); ep.first != ep.second; ++ep.first){
            	OutEdgeIter e_iter= ep.first;
            	ED ed = *e_iter;
            	VD target = boost::target(ed, *graph_);
                Atom const& a =  (*graph_)[target];
                AtomType const& at = (*atom_types_)[ a.atom_type_index() ];
                if( at.is_polar_hydrogen() ) return true;
            }
            return false;
        }

        bool HeavyAtomWithHydrogensFilter::operator()(VD const vd) const{
            
            for(OutEdgeIterPair ep = boost::out_edges(vd, *graph_); ep.first != ep.second; ++ep.first){
            	OutEdgeIter e_iter= ep.first;
                ED ed = *e_iter;
                VD target = boost::target(ed, *graph_);
                Atom const& a =  (*graph_)[target];
                AtomType const& at = (*atom_types_)[ a.atom_type_index() ];
                if( at.is_hydrogen() ) return true;
            }
            return false;
        }
        

        bool HydrogenAtomFilter::operator()(VD const vd) const{
            return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_hydrogen();
        }
        
        bool AromaticAtomFilter::operator()(VD const vd) const{
            return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_aromatic();
        }
        bool PolarHydrogenFilter::operator()(VD const vd) const{
            return (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_polar_hydrogen();
        }

        bool APolarHydrogenFilter::operator()(VD const vd) const{
            return  (*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_hydrogen() && !(*atom_types_)[ (*graph_)[vd].atom_type_index() ].is_polar_hydrogen();
        }

        
        
    }
}
///////////////////////////////////////////////////////////////


