// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/TMscore.hh
/// @brief  reimplementation of TMscore sequence-based structure superposition
/// @author Hahnbeom Park

#ifndef INCLUDED_core_scoring_TMscore_HH
#define INCLUDED_core_scoring_TMscore_HH

// Project headers
#include <core/types.hh>
#include <core/conformation/Residue.fwd.hh>
#include <core/pose/Pose.fwd.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>

namespace core {
namespace scoring {

////////////
class TMscoreStore
{
public:

	TMscoreStore()
	{
		clear();

		TM_max = 0.0;
		TM10_max = 0.0;
		maxsub_max = 0.0;
		n_GDT05_max = 0;
		n_GDT1_max = 0;
		n_GDT2_max = 0;
		n_GDT4_max = 0;
		n_GDT8_max = 0;
		GDTTS = 0.0;
		GDTHA = 0.0;
		maxsub = 0.0;
		TM = 0.0;
		TM10 = 0.0;
	}
	~TMscoreStore(){}

	void
	set_nseq(core::Size const nseq ){ nseq_ = (core::Real)(nseq); }

	inline
	void clear(){
		n_GDT05 = 0;                 ///for GDT-score, # of dis<0.5;
		n_GDT1 = 0;                  ///for GDT-score, # of dis<1;
		n_GDT2 = 0;                  ///for GDT-score, # of dis<2;
		n_GDT4 = 0;                  ///for GDT-score, # of dis<4;
		n_GDT8 = 0;                  ///for GDT-score, # of dis<8;
		maxsub_sum = 0;
		TM_sum = 0;
		TM_sum10 = 0;
	}

	inline
	void add_residue_dis( core::Real const d0,
		core::Real const dis )
	{
		if ( dis<=8 ) n_GDT8++;
		if ( dis<=4 ) n_GDT4++;
		if ( dis<=2 ) n_GDT2++;
		if ( dis<=1 ) n_GDT1++;
		if ( dis<=0.5 ) n_GDT05++;

		/// for MAXsub-score:
		if ( dis<3.5 ) maxsub_sum += 1/(1+(dis/3.5)*(dis/3.5));

		/// for TM-score:
		TM_sum += 1/(1+(dis/d0)*(dis/d0));

		/// for TM-score10:
		if ( dis < 10.0 ) TM_sum10 += 1/(1+(dis/d0)*(dis/d0));
	}

	inline
	void update()
	{
		// update score if get better
		if ( TM_max < TM ) TM = TM_max;
		if ( TM10_max< TM10 ) TM10_max=TM10;
		if ( maxsub_max<maxsub ) maxsub_max=maxsub;
		if ( n_GDT05_max < n_GDT05 ) n_GDT05_max=n_GDT05;
		if ( n_GDT1_max < n_GDT1 ) n_GDT1_max=n_GDT1;
		if ( n_GDT2_max < n_GDT2 ) n_GDT2_max=n_GDT2;
		if ( n_GDT4_max < n_GDT4 ) n_GDT4_max=n_GDT4;
		if ( n_GDT8_max < n_GDT8 ) n_GDT8_max=n_GDT8;
	}


	// Multiply by 100.0 not to lose precision when printed out as silent!
	inline
	void apply()
	{
		maxsub = 100.0* maxsub_sum/nseq_; ///MAXsub-score;
		TM =     100.0* TM_sum/nseq_; ///TM-score;
		TM10 =   100.0* TM_sum10/nseq_; ///TM-score10;
		GDTTS =  100.0*(n_GDT1_max +n_GDT2_max+n_GDT4_max+n_GDT8_max)/(4.0*nseq_);
		GDTHA =  100.0*(n_GDT05_max+n_GDT1_max+n_GDT2_max+n_GDT4_max)/(4.0*nseq_);
	}

public:

	core::Real n_cut;
	core::Real nseq_;
	core::Real nseq_ref_;

	core::Real TM, TM10;
	core::Real GDTTS;
	core::Real GDTHA;
	core::Real maxsub;

	core::Real TM_sum;
	core::Real TM_sum10;
	core::Real maxsub_sum;
	core::Size n_GDT05, n_GDT1, n_GDT2, n_GDT4, n_GDT8;

	core::Real TM_max;
	core::Real TM10_max;
	core::Real maxsub_max;
	core::Size n_GDT05_max;        ///number of residues<0.5;
	core::Size n_GDT1_max;         ///number of residues<1;
	core::Size n_GDT2_max;         ///number of residues<2;
	core::Size n_GDT4_max;         ///number of residues<4;
	core::Size n_GDT8_max;         ///number of residues<8;

}; // end class TMscoreStore

class TMscore
{

public:
	TMscore( ObjexxFCL::FArray2D< core::Real > const &p1 );
	~TMscore();

	void
	apply( ObjexxFCL::FArray2D< core::Real > const &p2 );

	void
	apply( ObjexxFCL::FArray2D< core::Real > const &p2,
		utility::vector0< core::Vector > &u,
		core::Vector &t,
		bool const get_ut );

	// Accessor
	core::Real get_TMscore() const { return score_.TM_max; }
	core::Real get_GDTTS() const { return score_.GDTTS; }
	core::Real get_GDTHA() const { return score_.GDTHA; }

private:

	core::Size iA( core::Size const i ) const { return iA_[i]; }
	core::Size iB( core::Size const i ) const { return iB_[i]; }

	void
	set_default();

	void
	convert_FArray2D_to_vector0( ObjexxFCL::FArray2D< core::Real > const & p1,
		utility::vector0< core::Vector > &xyz
	);

	void
	get_ali_params();

	utility::vector0< core::Size >
	score_fun( core::Real const d,
		utility::vector0< core::Vector > const & vt,
		TMscoreStore &score ) const;


	void
	extend( core::Real const &d,
		TMscoreStore &score,
		utility::vector0< core::Vector > &u,
		core::Vector &t,
		utility::vector0< core::Size > &i_ali,
		utility::vector0< core::Size > &k_ali0
	) const;

	utility::vector0< core::Vector >
	get_transrot_ref( utility::vector0< core::Size > const &k_ali,
		utility::vector0< core::Vector > &u,
		core::Vector &t
	) const;

private:

	TMscoreStore score_;
	core::Size nseq_;
	core::Real d0_;
	core::Real d0_search_;
	utility::vector1< core::Size > L_ini_; // Only the vector1 being used!
	utility::vector0< core::Size > iA_, iB_; // Residue indexer for ref & pose
	utility::vector0< core::Vector > xyza_, xyzb_;


};

}
}

#endif
