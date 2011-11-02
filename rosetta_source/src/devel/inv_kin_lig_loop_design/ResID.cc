// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   devel/InvKinLigLoopDesign/ResID.cc
///
/// @brief
/// @author

// Unit headers
#include <devel/inv_kin_lig_loop_design/ResID.hh>

// Project headers
#include <core/pose/PDBInfo.hh>

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


namespace devel {

  namespace inv_kin_lig_loop_design {

    // ===============================================
    // ==================== ResID ====================
    // ===============================================

    const char* chain_chars = " ABCDEFGHIJKLMNOPQRSTUVWXYZ01234567890"; // this is kind of dumb

    ResID::ResID() : chain_id('_'), res_num(0) {}

    ResID::ResID( char chain_id, int res_num ) : chain_id(chain_id), res_num(res_num) {
      if( this->chain_id == ' ' ) {
				this->chain_id = '_';
      }
    }

    ResID::ResID( core::conformation::Residue const& r ) : chain_id( chain_chars[ r.chain() ] ), res_num(r.seqpos()) {
      // !!! don't assume that this will match the resid obtained directly from the pdb
      if( chain_id == ' ' ) {
				chain_id = '_';
      }
    }

    bool ResID::operator==(ResID const& other ) const {
      return chain_id == other.chain_id && res_num == other.res_num;
    }

    bool ResID::operator<(ResID const& other ) const {
      return chain_id < other.chain_id || (chain_id == other.chain_id && res_num < other.res_num );
    }

    //   bool operator<(ResID const& a, ResID const& b ) {
    //     return a::operator<( b );
    //   }

    istream& operator>>(istream& in, ResID& res_id ) {

      std::string s;
      in >> s;

      size_t p0 = s.find(':');
      //size_t p1 = s.find(':',p0+1);
      assert( p0 == 1 );
      assert( s.size() >= 2 );


      res_id.chain_id = s[0]; // NB: this will never be " "
      //res_id.res_num  = atoi( s.substr(p0+1,p1).c_str() );
      res_id.res_num  = atoi( s.substr(p0+1,s.size()).c_str() );
      //     if( p1 < s.size() ) {
      //       name = AtomType( s.substr(p1) );
      //     }

      return in;
    }

    ostream& operator<<(ostream& out, ResID const& res_id ) {
      return out << res_id.chain_id << ":" << res_id.res_num;
    }

    // ====================================================
    // ==================== get_resids ====================
    // ====================================================

    resids_type ResID::get_resids( core::pose::Pose& pose ) {

      resids_type rval;

      size_t n = pose.n_residue();
      assert( n == pose.pdb_info()->nres() );

      for( size_t i = 1; i <= n; ++i ) {
				ResID const res_id( pose.pdb_info()->chain(i), pose.pdb_info()->number(i));
				rval[ res_id ] = const_cast<core::conformation::Residue*>(&pose.residue(i));
				//cout << res_id << " => " << &pose.residue(i) << endl;
      }

      return rval;
    }

  } // LoopDesign

} // devel
