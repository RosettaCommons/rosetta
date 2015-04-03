// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/Mover.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_MOVER_HH
#define DEVEL_INVKINLIGLOOPDESIGN_MOVER_HH

#ifdef WIN32
	#define _USE_MATH_DEFINES
#endif
#include <cmath>


#include <core/pack/dunbrack/RotamerLibraryScratchSpace.fwd.hh>

#include <iosfwd>
#include <cmath>
#include <devel/inv_kin_lig_loop_design/Fragment.hh>

#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>


namespace devel {

  namespace inv_kin_lig_loop_design {

    using namespace std;
    using core::conformation::Residue;
    using platform::Size;

    struct Mover {

      Mover(core::pose::Pose* pose_);

		~Mover();

      //void setPhiPsi(Peptide* p, const double phi, const double psi ) const;
      //void setPhiPsi(const int seqpos, const double phi, const double psi ) const;

      void clear();

      // Fragments

			//       void setFragment(const Residue* r, const int frag_len, const Fragment::Entry& e) const;
			//       void setFragmentFile(const Residue* r, const int frag_len, const Fragment::File* frag_file) const;
			//       const Fragment::File* getFragmentFile(const Residue* r, const int len) const;
			//       void clearFragmentFiles() const;
			//       void randomFragments(const char chain_id, const int from, const int to, const int frag_len) const;
			//       void insertFragment (const char chain_id, const int from, const int to, const int frag_len) const;
			//       void randomFragments(const Residue* a, const Residue* b, const int frag_len) const;
			//       void insertFragment (const Residue* a, const Residue* b, const int frag_len) const;


      void setPhiPsi(const int seqpos, const double phi, const double psi ) const;
      void setFragment(const int seqpos, const int frag_len, const Fragment::Entry& e) const;
      void insertFragment(const int seqpos_a, const int seqpos_b, const int frag_len ) const;
      void randomFragments(const int seqpos_a, const int seqpos_b, const int frag_len ) const;
      const Fragment::File* getFragmentFile(int seqpos, const int len) const;

			//     // Small moves
      void smallMoves(const int seqpos_a, const int seqpos_b, const double sigma = ( M_PI / 180.0 ) ) const ;
			//     void smallMoves(const Residue* a, const Residue* b, const double sigma = DEG2RAD*1.0) const;

			//     // Rotamers
			//     const vector<Kin::Rotamer>& getRotamers(const char chain_id, const int res_num ) const;
      void randomRotamer(const int seqpos) const;
			//     void randomRotamer(const Residue* r) const;

      void setExtended(const Residue* a, const Residue* b) const;
      void setHelix(const Residue* a, const Residue* b) const;

      // Chainbreak
      void setChainbreakReq(const Residue* a, const Residue* b); // if there's a chainbreak in a,b the requested chainbreak gets set to the actual value

      //   // Ccd move
      //   private:
      //     typedef map<const Kin::Dihedral*,MoverInfo::CcdTarget> m_ccd_targets_t;
      //     mutable m_ccd_targets_t mCcdTargets;
      //     mutable Tube tube_ccd;
      //     mutable std::size_t num_res_ccd;
      //     mutable Polypeptide* pp_ccd;
      //     mutable Peptide* p_ccd;
      //     mutable Peptide* q_ccd;

      //     void addCcdPoints(vector<pair<Point3,Point3> >& v, const vector<const Residue*>& vRes ) const;

      //   public:

      //     void setCcdNumRes(const std::size_t num_res); // sets the number of residues on either side of the chainbreak to consider when performing the ccd

      //     //void setCcdTargets(const Residue* from, const Residue* to, const Residue* target) const;
      //     void clearCcdTargets() const;
      //     void setCcdTargets(const Abng::Graph& g);

      //     void ccdMove(const Kin::Dihedral* dih) const;
      //     void ccdMove(const Residue* r) const;
      //     void ccdMove(const Residue* a, const Residue* b) const;
      //     void ccdMove(const char chain_id, const int res_num) const;
      //     void ccdMove(const char chain_id, const int from, const int to) const;

			//     void ccdClose(const Residue* a, const Residue* b, Abng::Graph& g) const;

			// //     void addDipeptide_ccd(const Residue* r);

    private:
      //mutable map<const Residue*, vector<Kin::Rotamer> > mvRot;
      mutable map<pair<int,int>,const Fragment::File*> mFragmentFiles; // will always map to a librarian supplied fragfile


		core::pose::Pose* pose; // unacceptable
		core::pack::dunbrack::RotamerLibraryScratchSpaceOP scratch_;

		}; // class Mover

  } // LoopDesign

} // devel


#endif // UTIL_MOVER_HH
