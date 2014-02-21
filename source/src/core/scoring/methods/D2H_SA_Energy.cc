// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/D2H_SA_Energy.cc
/// @brief  sa energy function definition.
/// @author Björn Wallner


// Unit headers
#include <core/scoring/methods/D2H_SA_Energy.hh>
#include <core/scoring/methods/D2H_SA_EnergyCreator.hh>

// Package headers
#include <core/scoring/sasa.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/EnvPairPotential.hh>
#include <basic/datacache/CacheableData.hh> //TMP HACK
#include <core/scoring/EnergiesCacheableDataType.hh>
#include <core/pose/datacache/CacheableDataType.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <numeric/statistics.functions.hh>
#include <numeric/util.hh>
//#include <math.h>

namespace core {
namespace scoring {
namespace methods {

static basic::Tracer TR("core.scoring.methods.D2H_SA_Energy.cc");

/// @details This must return a fresh instance of the D2H_SA_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
D2H_SA_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return new D2H_SA_Energy;
}

ScoreTypes
D2H_SA_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( d2h_sa );
	return sts;
}

/// c-tor
D2H_SA_Energy::D2H_SA_Energy() :
	parent( new D2H_SA_EnergyCreator )
{
	HDX_data_defined_=false;
	chain_=0;
	rsa_index_start_=0;
	rsa_index_end_=0;
	reweight_=1;
	
	Size total_length(0);
	std::string line;
	if(basic::options::option[ basic::options::OptionKeys::in::file::d2h_sa_reweight ].user()) {
		reweight_=basic::options::option[ basic::options::OptionKeys::in::file::d2h_sa_reweight]();
	}
	if(basic::options::option[ basic::options::OptionKeys::in::file::HDX ].user()) {
		TR << "Reading Hydrogen Exchange data\n";
		std::string HDX_datafile(basic::options::option[ basic::options::OptionKeys::in::file::HDX ]());
		

		utility::io::izstream stream (HDX_datafile);
		getline(stream,line);
		TR << line << std::endl;
		getline(stream,line);
		{
			std::istringstream l(line);
			TR << line << std::endl;
			l>> total_length >> chain_ >> rsa_index_start_ >> rsa_index_end_;
		}
		//Reading to local variables to allow for missing data, which is not implemented ... :-)
		//utility::vector1 < core::Size > positions,group; //(total_length,0);
		//		utility::vector1 < core::Size > group(total_length,0);
		//		utility::vector1 < core::Real > data; //(total_length,0.0);
		std::map<core::Size, core::Real> data_map;
		std::map<core::Size, core::Size> group_map;
		

		//	Size pos;
		//Real tmp;
		for(Size i=1;i<=total_length;i++)
			{
				getline(stream,line);
				{
				std::istringstream l(line);
				Size pos;
				Real tmp;
				Size g;
				l >> pos >> tmp >> g; 
				data_map[pos]=tmp;
				group_map[pos]=g;

				//positions.push_back(pos);
				//data.push_back(tmp);
				//group.push_back(g);



			}
			TR << line << std::endl ;
		}


		

		if(total_length==0)
			{
				utility_exit_with_message("bad format for d2h file total_length=0");
			}
		//	std::cout << sd << " " << mean << " "<< sd2 << " " << "\n";
		utility::vector1 < core::Size > keys;
		for (std::map<Size, Real>::iterator i = data_map.begin(); i != data_map.end(); ++i)
    {
        keys.push_back(i->first);
    }
    sort(keys.begin(), keys.end()); 
		for(Size i=1;i<=keys.size();i++)
			{
				data_.push_back(data_map[keys[i]]);
				position_.push_back(keys[i]);
				group_.push_back(group_map[keys[i]]);
			}
		
//all_data_.resize(positions[total_length],-1);
//all_group_.resize(positions[total_length],-1);
//chain_=5;
//for(Size i=1;i<=total_length;i++)
//	{
//		all_data_[positions[i]]=data[i]; //-mean_;
//		all_group_[positions[i]]=group[i]; //-mean_;
//	}
//for(Size i=1;i<=all_data_.size();i++)
//	{
//		TR.Debug << i << " " << all_data_[i] << std::endl;
//		if(all_data_[i]!=-1) {
//			data2_.push_back(all_data_[i]);
//			position2_.push_back(i);
//			group2_.push_back(all_group_[i]);
//		}
//
//	}
//
//		for(Size i=1;i<=data_.size();i++)
//			{
//				TR.Debug <<  data_[i] << " " << position_[i] << " - " << group_[i] << std::endl;
//			}


		Real mean = numeric::statistics::mean(data_.begin(),data_.end(),0.0);
		Real sd = numeric::statistics::std_dev_with_provided_mean(data_.begin(),data_.end(),mean);
		mean_=mean;
		sd_=sd;
		HDX_data_defined_=true;
		stream.close();
		stream.clear();
	}

}

/// clone
EnergyMethodOP
D2H_SA_Energy::clone() const
{
	return new D2H_SA_Energy;
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////

void
D2H_SA_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const &,
	EnergyMap & totals
) const {
	using namespace conformation;

	if(!HDX_data_defined_) {
		totals [ d2h_sa ] = 0;
		return;
	}
	bool fullatom=pose.is_fullatom();


	Size len(data_.size());
	Size nres( pose.total_residue() );
	//Size max_res( pose.total_residue() );  // unused
	const Real probe_radius(1.4); //default water probe
	//std::cout << "Length " << len << "\n";
	Size chain_for_rsa(1);
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		using namespace core::conformation::symmetry;
		SymmetricConformation const & symm_conf(dynamic_cast< SymmetricConformation const & > ( pose.conformation() ) );
		SymmetryInfoCOP symm_info( symm_conf.Symmetry_Info() );
		//max_res = symm_info->num_independent_residues();  // unused
		chain_for_rsa=pose.residue(symm_info->bb_follows(1)).chain();
		//for(Size i=1;i<=nres;i++) {
		//	TR.Debug  << "BB_FOLLOWS " << i << " " << symm_info->bb_follows(i) << std::endl;
		//}
	}

	//TR.Debug << "max_res: " << max_res << " nres: " << nres << " chain_for_rsa: " << chain_for_rsa << std::endl;

	id::AtomID_Map< Real > atom_sasa;
	utility::vector1< Real > rsd_sasa(nres,0.0);
	utility::vector1< Real > rsd_sasa_raw(len,0.0);
	utility::vector1< Real > debug;
	utility::vector1< Real > rsd_sasa_grouped(len,0.0);
	utility::vector1< Real > rsd_sasa_grouped_median(len,0.0);
	utility::vector1< Real > tmp;
	utility::vector1< Real > group_mean;
	utility::vector1< Real > group_median;
	Size rsa_index=0;
	Size rsa_index_start=0;
	Size rsa_index_end=0;

	pose::PDBInfoCOP pdb_info = pose.pdb_info();
	if(rsa_index_start_ == 0) { 
		for(Size i=1;i<=nres;i++) {
			conformation::Residue const & rsd( pose.residue(i) );
			if(rsd.chain() == chain_for_rsa) {
				rsa_index_end=i;
				if(rsa_index==0) {
					rsa_index=i;
					rsa_index_start=i;
				}
			}
		}
		//set_rsa_range(rsa_index_start,rsa_index_end);
	} else {
		rsa_index_end=rsa_index_end_;
		rsa_index=rsa_index_start_;
	}

	TR.Debug << "Look for chain: " << chain_for_rsa << " found it between " << rsa_index << " and " << rsa_index_end << std::endl;
	Size counter2(0);

	if(fullatom) {
		//TR.Debug << "fullatom" << std::endl; 
		calc_per_atom_sasa( pose, atom_sasa, rsd_sasa, probe_radius );
	} else {
		//TR.Debug << "centroid" << std::endl;
		utility::vector1< Size > cen8(nres,0);
		utility::vector1< Size > cen8_2(nres,0);

		for ( Size i = rsa_index_start; i <= rsa_index_end;++i ) {
			Vector const v1( pose.residue(i).nbr_atom_xyz() );
			for ( Size j = 1; j <= nres; ++j ) {
				Size aa_i=pose.residue(i).nbr_atom();
				Size aa_j=pose.residue(j).nbr_atom();
				if(!pose.residue(i).atom_type(aa_i).is_virtual() &&
					 !pose.residue(j).atom_type(aa_j).is_virtual() && 
					 !(j>=rsa_index_start && j<=rsa_index_end && i>=j)) {
					counter2++;
					Vector const v2( pose.residue(j).nbr_atom_xyz() );
					Real dist(v1.distance_squared( v2 ));
					if ( dist < 64 ) {

						cen8[i]++;
						cen8[j]++;
					}
				}
			}
		}
		for ( Size i = rsa_index_start; i <= rsa_index_end;++i ) {
			rsd_sasa[i]=-9.8017*cen8[i]+137.1306;
			if(rsd_sasa[i] < 0 || rsd_sasa[i]>137) {
				rsd_sasa[i]=0;
			}
		}
	}

	rsd_sasa_raw[1]=rsd_sasa[rsa_index_start+position_[1]-1];
	tmp.push_back(rsd_sasa[rsa_index_start+position_[1]-1]);
	utility::vector1< Size > group_index(len,0);
	for(Size i=2;i<=len;i++) {
		group_index[i-1]=group_mean.size()+1;
		if(group_[i-1]  != group_[i]) {
			group_mean.push_back(numeric::statistics::mean(tmp.begin(),tmp.end(),0.0)); 
			group_median.push_back(numeric::median(tmp));
			tmp.clear();
		}
		
		tmp.push_back(rsd_sasa[rsa_index_start+position_[i]-1]);
		rsd_sasa_raw[i]=rsd_sasa[rsa_index_start+position_[i]-1];
	}
	group_index[len]=group_mean.size()+1;
	group_mean.push_back(numeric::statistics::mean(tmp.begin(),tmp.end(),0.0));
	group_median.push_back(numeric::median(tmp));

	for(Size i=1;i<=len;i++) {
		rsd_sasa_grouped[i]=group_mean[group_index[i]];
		rsd_sasa_grouped_median[i]=group_median[group_index[i]];
		TR.Debug << "MEAN: " << i << " " << position_[i] << " " << rsa_index_start+position_[i]-1 << " " << rsd_sasa[rsa_index_start+position_[i]-1] << " " << rsd_sasa_raw[i] << " " << rsd_sasa_grouped[i] << " " << rsd_sasa_grouped_median[i] << " " << data_[i] << " " << group_[i] << " " << group_index[i] << " " << len << std::endl;
	}

	Real m = numeric::statistics::mean(rsd_sasa_grouped.begin(),rsd_sasa_grouped.end(),0.0);
	Real sd = numeric::statistics::std_dev_with_provided_mean(rsd_sasa_grouped.begin(),rsd_sasa_grouped.end(),m);
	Real corr=numeric::statistics::corrcoef_with_provided_mean_and_std_dev(data_,mean_,sd_,rsd_sasa_grouped,m,sd);

	Real z = -50 * log((1+corr)/(1-corr)); //Fisher transform, slightly better behavior for significant correlations times -100.
	totals [ d2h_sa ] = reweight_*z; 
} // finalize_total_energy


core::Size
D2H_SA_Energy::version() const
{
	return 1; // Initial versioning
}

} // methods
} // scoring
} // core
