// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief

#include <protocols/viewer/viewers.hh>

#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>

#include <core/chemical/AA.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/optimization/Multifunc.hh>
#include <core/optimization/Minimizer.hh>

#include <basic/options/util.hh>
#include <devel/init.hh>

#include <boost/format.hpp>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/random/random.hh>
#include <numeric/constants.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/string.functions.hh>
#include <core/kinematics/Jump.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <queue>


OPT_1GRP_KEY(Integer, pseudoca, ncas)
OPT_1GRP_KEY(Integer, pseudoca, ncycles)
OPT_1GRP_KEY(Integer, pseudoca, beamwidth)


struct pseudoCA {
	pseudoCA( core::Vector X_in, core::Real score_in ) { X=X_in; score=score_in; }
	core::Vector X;
	core::Real score;
};

struct pseudoCAbond {
	pseudoCAbond( core::Size i_in, core::Size j_in, core::Real score_in ) {
		i=i_in;
		j=j_in;
		score=score_in;
	}
	core::Size i, j;
	core::Real score;
};

// indexing
// traceB:  N (N-1) ... 2 1
// traceF:  N (N+1) ... (N+M-1) (N+M)
class pseudoTrace {
private:
	utility::vector1< core::Size > traceF_;
	utility::vector1< core::Size > traceB_;
	core::Real score_;

	core::Real DENS_WT, DIST_WT, ANGLE_WT, OVERLAP_WT;

public:
	pseudoTrace(core::Size seed, utility::vector1<pseudoCA> const &calist) {
		traceF_.push_back(seed);
		traceB_.push_back(seed);

		DENS_WT = 5.0;
		DIST_WT = 0.1;
		ANGLE_WT = 1.0;
		OVERLAP_WT = 10.0;

		score_ = DENS_WT*calist[seed].score;
	}

	void add_to_C(core::Size ext, utility::vector1<pseudoCA> const &calist) {
		traceF_.push_back(ext);
		score_ += DENS_WT*calist[ext].score;

		// distance
		core::Real dist = (calist[get(size())].X-calist[get(size()-1)].X).length() - 3.8;
		score_ += DIST_WT*(dist*dist);

		// angle
		if (size()>=3) {
				core::Real angle = numeric::angle_radians(calist[get(size())].X, calist[get(size()-1)].X, calist[get(size()-2)].X);
				score_ += ANGLE_WT*(angle-1.57079633)*(angle-1.57079633);
		}

		// clash
		for (int i=1; i<=(int)size()-3; ++i) {
			core::Real overlap = std::max( 3.0 - (calist[get(1)].X-calist[get(i)].X).length() , 0.0);
			score_ += OVERLAP_WT*(overlap*overlap);
		}
	}

	void add_to_N(core::Size ext, utility::vector1<pseudoCA> const &calist) {
		traceB_.push_back(ext);
		score_ += DENS_WT*calist[ext].score;

		// distance
		core::Real dist = (calist[get(1)].X-calist[get(2)].X).length() - 3.8;
		score_ += DIST_WT*(dist*dist);

		// angle
		if (size()>=3) {
				core::Real angle = numeric::angle_radians(calist[get(1)].X, calist[get(2)].X, calist[get(3)].X);
				score_ += ANGLE_WT*(angle-1.57079633)*(angle-1.57079633);
		}

		// clash
		for (core::Size i=4; i<=size(); ++i) {
			core::Real overlap = std::max( 3.0 - (calist[get(1)].X-calist[get(i)].X).length() , 0.0);
			score_ += OVERLAP_WT*(overlap*overlap);
		}
	}

	core::Size get_C() { return traceF_[traceF_.size()]; }
	core::Size get_N() { return traceB_[traceB_.size()]; }
	core::Size size() const { return traceF_.size() + traceB_.size() - 1; }
	core::Size get(core::Size i) {
		if (i<=traceB_.size())
			return traceB_[traceB_.size()+1-i];
		else
			return traceF_[i-traceB_.size()+1];
	}

	core::Real score() const  {
		return score_;
	}

	bool is_subset(const pseudoTrace& other) const {
		if (size() > other.size()) return false;

		// ensure every point I have is also in other
		for (core::Size i=1; i<=traceF_.size(); ++i) {
			if (!other.contains(traceF_[i])) return false;
		}

		return true;
	}

	bool contains(core::Size idx) const {
		if (std::find( traceF_.begin(), traceF_.end(), idx ) != traceF_.end() ) return true;
		if (std::find( traceB_.begin(), traceB_.end(), idx ) != traceB_.end() ) return true;
		return false;
	}
};

////////////////////////////////////

class pseudoTraceRecords {
public:
	pseudoTraceRecords(core::Size N) {
		N_ = N;
		data_.reserve( N_+1 );
		high_idx_ = 0;
	}

	void
	addTrace( pseudoTrace const &newtrace ) {
		if (size() == N_ && newtrace.score() > data_[high_idx_].score() )
			return;

		// equality check
		for (core::Size i=1; i<=data_.size(); ++i) {
			if (newtrace.is_subset( data_[i] ))
				return;
			if (data_[i].is_subset( newtrace )) {
				data_[i] = newtrace;
				find_highest_score();
			}
		}

		if (size() < N_) {
			data_.push_back( newtrace );
			find_highest_score();
		} else {
			data_[high_idx_] = newtrace;
			find_highest_score();
		}
	}

	pseudoTrace get( core::Size i ) { return data_[i];	}

	core::Size size() { return data_.size(); }

	void clear() {
		high_idx_ = 0;
		data_.clear();
	}

	core::Real maxscore( ) { return data_[high_idx_].score(); }

private:
	void
	find_highest_score( ) {
		high_idx_ = 1;
		for (core::Size i=2; i<=data_.size(); ++i) {
			if (data_[i].score() > data_[high_idx_].score())
				high_idx_ = i;
		}
	}

	core::Size N_;
	utility::vector1<pseudoTrace> data_;

	// highest-scoring element + location
	core::Size high_idx_;
};


bool operator<(const pseudoCA& lhs, const pseudoCA& rhs) {
  return lhs.score < rhs.score;
}

bool operator<(const pseudoCAbond& lhs, const pseudoCAbond& rhs) {
  return lhs.score < rhs.score;
}

bool operator<(const pseudoTrace& lhs, const pseudoTrace& rhs) {
  return lhs.score() < rhs.score();
}


void
run_trace(utility::vector1<pseudoCA> selectPseudos, utility::vector1<pseudoTrace> & pseudoTraces) {
	core::Size beamwidth = basic::options::option[ basic::options::OptionKeys::pseudoca::beamwidth ]();

	pseudoTraces.clear();

	pseudoTraceRecords Q(beamwidth), Qnew(beamwidth);

	// start with single-point traces
	for (core::Size i=1; i<=selectPseudos.size(); ++i) {
		Q.addTrace( pseudoTrace(i,selectPseudos) );
	}

	core::Real best_score = Q.maxscore();
	std::cerr << "cycle 0 best score = " << best_score << std::endl;
	core::Real last_score = 99999;
	core::Size cyc=1;
	while (last_score != best_score) {
		last_score = best_score;

		for (core::Size i=1; i<=Q.size(); ++i) {
			pseudoTrace Q_i = Q.get(i);

			Qnew.addTrace( Q_i );

			// try to grow
			for (core::Size i=1; i<=selectPseudos.size(); ++i) {
				core::Real distN2 = (selectPseudos[i].X - selectPseudos[Q_i.get_N()].X).length_squared();
				if (distN2 <= 4.3*4.3 && distN2 >= 3.3*3.3) {
					pseudoTrace Q_i_N = Q_i;
					Q_i_N.add_to_N(i, selectPseudos);
					Qnew.addTrace( Q_i );
				}

				core::Real distC2 = (selectPseudos[i].X - selectPseudos[Q_i.get_C()].X).length_squared();
				if (distC2 <= 4.3*4.3 && distC2 >= 3.3*3.3) {
					Q_i.add_to_C(i, selectPseudos);
					Qnew.addTrace( Q_i );
				}
			}
		}

		Q = Qnew;
		best_score = Q.maxscore();
		Qnew.clear();

		std::cerr << "cycle " << cyc++ << " best score = " << best_score << std::endl;
	}

	for (core::Size i=1; i<=Q.size(); ++i) {
		pseudoTraces.push_back(Q.get(i));
	}
	std::sort( pseudoTraces.begin(), pseudoTraces.end() );
}


/// Bfactor multifunc
class PseudoCAMultifunc : public core::optimization::Multifunc {
public:
	PseudoCAMultifunc(
		utility::vector1< core::Vector > const & /*pseudo_cas*/,
		utility::vector1< std::pair< core::Size, core::Size > > /*bonded_cas*/
	) {

	}

	virtual ~PseudoCAMultifunc() {}

	virtual	core::Real
	operator ()( core::optimization::Multivec const & vars ) const;

	virtual	void
	dfunc( core::optimization::Multivec const & vars, core::optimization::Multivec & dE_dvars ) const;

	virtual	void
	dump( core::optimization::Multivec const & vars1, core::optimization::Multivec const & vars2 ) const;

private:
};

///////////////////////////////////////////////////////////////////////////////

void
traceCAs() {
	core::Size ncas = basic::options::option[ basic::options::OptionKeys::pseudoca::ncas ]();

	runtime_assert(ncas != 0);

	// initial placement
	core::scoring::electron_density::ElectronDensity &edm = core::scoring::electron_density::getDensityMap();

	// sample at 0.333A grid saving best #ca density values
	core::Real xlen = (edm.get_f2c()*core::Vector(1.0,0,0)).length();
	core::Real ylen = (edm.get_f2c()*core::Vector(0,1.0,0)).length();
	core::Real zlen = (edm.get_f2c()*core::Vector(0,0,1.0)).length();

	numeric::xyzVector< core::Size > nsteps(
		(core::Size)std::floor( 0.5 + 3*xlen ),
		(core::Size)std::floor( 0.5 + 3*ylen ),
		(core::Size)std::floor( 0.5 + 3*zlen )
	);


	std::priority_queue<pseudoCA> Q;

	for (core::Size nx = 0; nx<nsteps[0]; ++nx)
	for (core::Size ny = 0; ny<nsteps[1]; ++ny)
	for (core::Size nz = 0; nz<nsteps[2]; ++nz) {
		core::Vector fX (
			(core::Real)nx/(core::Real)nsteps[0],
			(core::Real)ny/(core::Real)nsteps[1],
			(core::Real)nz/(core::Real)nsteps[2]
		);
		core::Vector cX = edm.get_f2c()*fX;
		core::Real score = -edm.matchPointFast( cX );
		if (Q.size() < 20*ncas) {  // may be too permissive (?)
			Q.push( pseudoCA(cX,score) );
		} else {
			if (score < Q.top().score) {
				Q.pop();
				Q.push( pseudoCA(cX,score) );
			}
		}
	}

	// make a bit sparser
 	core::Real cutoffdist = 1;
 	utility::vector1<pseudoCA> topPseudos, selectPseudos;
 	while (!Q.empty()) {
 		topPseudos.push_back(Q.top());
 		Q.pop();
 	}
 	std::sort(topPseudos.begin(), topPseudos.end());
 	for (core::Size i=1; i<=topPseudos.size(); ++i) {
 		bool tooclose=false;
 		for (core::Size j=1; j<=selectPseudos.size() && !tooclose; ++j) {
 			if ( (topPseudos[i].X - selectPseudos[j].X).length_squared() < cutoffdist*cutoffdist) tooclose=true;
 		}
 		if (!tooclose) {
 			selectPseudos.push_back( topPseudos[i] );
 			if (selectPseudos.size() >= ncas) break;
 		}
	}

	// dump as pdb
	std::ofstream oss( "ca_trace.pdb" );
	for (core::Size i=1; i<=selectPseudos.size(); ++i) {
		oss << boost::format("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n")
			% "ATOM" % i % " CA " % " " % "ALA" % "A" % i % " " % selectPseudos[i].X[0] % selectPseudos[i].X[1] % selectPseudos[i].X[2] % 1.0 % 1.0 % " C" % "  ";
	}

	// detect CA connectivity
	utility::vector1<  pseudoTrace > pseudoTraces;
	run_trace(selectPseudos, pseudoTraces);

	// minimize cycles
	// TO DO

	for (core::Size i=1; i<=std::min((core::Size)100,pseudoTraces.size()); ++i) {
		std::string name = "ca_trace_"+utility::to_string(i)+".pdb";
		std::ofstream oss_i( name.c_str() );
		for (core::Size j=1; j<=pseudoTraces[i].size(); ++j) {
				oss_i << boost::format("%-6s%5d %4s%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n")
					% "ATOM" % j % " CA " % " " % "ALA" % "A" % j % " "
					% selectPseudos[pseudoTraces[i].get(j)].X[0] % selectPseudos[pseudoTraces[i].get(j)].X[1] % selectPseudos[pseudoTraces[i].get(j)].X[2]
					% 1.0 % 1.0 % " C" % "  ";
		}
	}
}


///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char * argv [] )
{
	try {
	// options, random initialization
	NEW_OPT(pseudoca::ncas, "number of cas", 100);
	NEW_OPT(pseudoca::ncycles, "number of cycles", 10);
	NEW_OPT(pseudoca::beamwidth, "number of cycles", 100);

	devel::init( argc, argv );
	traceCAs();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
