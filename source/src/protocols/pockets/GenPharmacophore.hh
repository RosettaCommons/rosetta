// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief A class to generate pharmacophore from RNA binding proteins
/// @author Ragul Gowthaman (ragul@ku.edu)
/// @author Yan Xia (seanxiay@ku.edu)

#ifndef INCLUDED_protocols_pockets_GenPharmacophore_hh
#define INCLUDED_protocols_pockets_GenPharmacophore_hh

#include <sstream>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/pose/Pose.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cmath>

using namespace std;

namespace protocols {
    namespace pockets {
        class SmallMol {
            private:
                string molName;
                string pdbContent;
                vector< vector<core::Real> > coordinates;
                core::Real rmsd;
                vector<core::Real> cen;
                int s; // size of the cluster
                SmallMol *parent; // parent of the cluster
                bool visited;

                vector<string> tokenize_connected(vector<string> &fields);

            public:
                SmallMol() : molName("no_name"), pdbContent(""), cen(3), s(1), parent(this), visited(0) { }
                SmallMol(const SmallMol &other);
                ~SmallMol();
                void add_atom(string line);
                void update_center();
                void set_name(string name) { molName = name; }
                string get_name() { return molName; }
                static core::Real calRMSD(SmallMol &mol1, SmallMol &mol2);
                int numberOfAtoms() const { return (int) coordinates.size(); }
                vector< vector<core::Real> > const & get_coordinates() { return coordinates; }
                void printCoordinates() const;
                void printContent() const;
                string getContent() const;
                core::Real get_rmsd() { return rmsd; }
                core::Real get_center(int c);
                core::Real cal_distance(SmallMol *other);
                core::Real cal_min_dist(SmallMol *other);
                bool operator < (const SmallMol &other) const { return rmsd < other.rmsd; }
                SmallMol *findRoot();
                bool connected(SmallMol *m);
                void connect(SmallMol *m);
                int get_size() { return s; }
                void set_size(int t) { s = t; }
                SmallMol *get_parent() { return parent; }
                void set_parent(SmallMol *p) { parent = p; }
                SmallMol *get_root() { return findRoot(); }
                bool get_visited() { return visited; }
                void set_visited(bool v) { visited = v; }
        };

        class UnionEdge {
            private:
                SmallMol *a;
                SmallMol *b;
                core::Real dist;
                // int s;

            public:
                UnionEdge(SmallMol *x, SmallMol *y) : a(x), b(y) { dist = a->cal_min_dist(b); }
                SmallMol *get_a() { return a; }
                SmallMol *get_b() { return b; }
                core::Real get_dist() const { return dist; }
                bool operator< (UnionEdge const &e) const { return dist < e.get_dist(); }
        };

        class GenPharmacophore : public utility::pointer::ReferenceCount {
            public:

                bool is_buried_ring(core::conformation::Residue const & rsd, core::Real const & ring_sasa, core::Real const & sasa_cutoff);
                core::Real get_RNAring_sasa( core::conformation::Residue const & rsd, int const & rsdno, core::id::AtomID_Map<core::Real> const & pose_atom_sasa );
                void get_ideal_hydrogenBond_atoms(core::pose::Pose const & protein_pose);
                void cluster_KeyFeatures(std::string const & input_filename, std::string const & output_filename) const;
                std::string extract_rna_rings_from_protein_rna_complex(core::pose::Pose const & protein_pose, core::pose::Pose const & rna_pose);
                std::string extract_Hbond_atoms_from_protein_rna_complex(core::pose::Pose const & protein_pose, core::pose::Pose const & rna_pose);
                std::string make_compatible_with_ROCS_custom_ForceField(std::string & input_filename);
                void print_string_to_PDBfile(std::string const & input_filename, std::string const & output_filename) const;
                std::string ReplaceString(std::string subject, const std::string& search, const std::string& replace);
        };

    } // pockets
} // protocols

#endif
