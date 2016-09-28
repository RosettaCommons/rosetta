// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/Burial_v2Energy.cc
/// @details I decided to make this a whole body method because I'm trying to calculate this on a non-canonical aa that is much bigger then 12 angstroms
/// @author TJ Brunette

#include <core/scoring/methods/Burial_v2Energy.hh>
#include <core/scoring/methods/Burial_v2EnergyCreator.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/Atom.hh>
#include <core/id/AtomID.hh>

#include <core/pose/metrics/CalculatorFactory.hh>
#include <core/pose/Pose.hh>

#include <basic/MetricValue.hh>
#include <basic/prof.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/EnergyMap.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <utility/io/izstream.hh>
#include <utility/file/FileName.hh>
#include <utility/string_util.hh>

#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.Burial_v2Energy" );

namespace core {
namespace scoring {
namespace methods {

/// @details This must return a fresh instance of the Burial_v2Energy class,
/// never an instance already in use
methods::EnergyMethodOP
Burial_v2EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new Burial_v2Energy );
}

ScoreTypes
Burial_v2EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back(burial_v2);
	return sts;
}



/// clone
EnergyMethodOP Burial_v2Energy::clone() const {
	return EnergyMethodOP( new Burial_v2Energy );
}

core::Real Burial_v2Energy::using_atom_distance(core::pose::Pose const & pose) const{
	Real score = 0;
	Real dist_cutoff = 8;
	Size residue_gap = 6; //number of residues to not measure. This is so I don't count contacts in the same chain as burial
	for(Size resid=1; resid<=residue_ids_.size(); ++resid){
		//get the number of cb within 8 angstroms of an atom in the side-chain
		conformation::Residue heme_rsd = pose.residue( residue_ids_[resid] );
		//std::cout << "examining resid" << residue_ids_[resid] << std::endl;
		for(Size ii=1; ii<=pose.total_residue(); ++ii){
			bool found = false;
			if(std::abs((int)residue_ids_[resid]-(int)ii)>(int)residue_gap){//don't measure any residue within 6 of the selected residue. Looking for long range measurements
				conformation::Residue protein_rsd = pose.residue(ii);
				numeric::xyzVector<core::Real> protein_xyz = pose.residue(ii).nbr_atom_xyz();
				for ( Size atomno = 1; atomno <= pose.residue_type(residue_ids_[resid]).natoms() && found==false; ++atomno ) {//ncaa atom
					//get the number of cb within 8 angstroms of an atom in the side-chain
					numeric::xyzVector<core::Real> heme_xyz = heme_rsd.xyz(atomno);
					Real dist = heme_xyz.distance(protein_xyz);
					if(dist < dist_cutoff){
						//std::cout << "non-heme:" << ii <<"," << dist << std::endl;
						found = true;
					}
				}
			}
			if(found)
				score--;
		}
	}
	return(score);
}

core::Real Burial_v2Energy::using_totalSasa(core::pose::Pose const & pose) const{
	using namespace core::pose::metrics;
	using basic::MetricValue;
	using namespace core;
	Real score=0;
	basic::MetricValue< utility::vector1< core::Real > > residue_sasa;
	pose.metric( "sasa", "residue_sasa", residue_sasa ); //this does not represent a new calculation since pose metric calculators are smart enough to figure it out
	runtime_assert( pose.total_residue() == (residue_sasa.value()).size() );
	for(Size resid=1; resid<=residue_ids_.size(); ++resid){
		core::Real this_sasa = residue_sasa.value()[residue_ids_[resid]];
		//std::cout << "this_sasa" << this_sasa << std::endl;
		score+=this_sasa;
	}
	//std::cout << "fa score" << score << std::endl;
	return(score);

}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void Burial_v2Energy::finalize_total_energy(
	core::pose::Pose & pose,
	ScoreFunction const & ,
	EnergyMap & totals
) const {
	if(pose.is_fullatom())
		totals[burial_v2] =using_totalSasa(pose)/10; //For HEM the sasa score was about 10x higher.
	else
		totals[burial_v2] =using_atom_distance(pose);
}


core::Size Burial_v2Energy::version() const
{
	return 1;
}


void Burial_v2Energy::init_from_file() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	std::string const & burial_fn( option[ in::file::burial ]()[1] );
	utility::io::izstream input(burial_fn);
	if ( !input.good() ) {
		std::string const msg( "Error opening file: " + burial_fn );
		utility_exit_with_message( msg );
	}
	using std::string;
	string line;
	getline(input,line); // header
	residue_ids_ = utility::string_split(line,',',core::Size());
	// for(Size ii=1; ii<=residue_ids_.size(); ++ii)
	// 	std::cout << "residue_ids_[1]" << residue_ids_[ii] << std::endl;
}

} // methods
} // scoring
} // core
