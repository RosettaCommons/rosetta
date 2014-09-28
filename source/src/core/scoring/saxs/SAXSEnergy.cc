// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/SAXSEnergy.cc
/// @brief  "Energy" based on a similarity of theoretical SAXS spectrum computed for a pose and the experimental data
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

// Unit headers
#include <core/scoring/saxs/FormFactorManager.hh>
#include <core/scoring/saxs/FormFactor.hh>
#include <core/scoring/saxs/SinXOverX.hh>
#include <core/scoring/saxs/SAXSEnergy.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorFA.hh>
#include <core/scoring/saxs/SAXSEnergyCreatorCEN.hh>
#include <core/scoring/saxs/SAXSEnergyCreator.hh>
#include <core/scoring/saxs/DistanceHistogram.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/io/pdb/file_data.hh>

#include <numeric/interpolation/spline/Interpolator.hh>
#include <numeric/interpolation/spline/SplineGenerator.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.hh>


// Utility headers
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>

// C++
#include <iomanip>
#include <string>

#include <utility/vector1.hh>

//Auto Headers
//#include <core/import_pose/import_pose.hh>

namespace core {
namespace scoring {
namespace saxs {

static thread_local basic::Tracer trSAXSEnergy( "core.scoring.saxs.SAXSEnergy" );

std::string SAXSEnergy::fa_cfg_file_("ff-rosetta-fa.cfg");
std::string SAXSEnergy::cen_cfg_file_("ff-rosetta-cen.cfg");

ScoreTypes SAXSEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( saxs_score );
	return sts;
}

methods::EnergyMethodOP SAXSEnergyCreator::create_energy_method( methods::EnergyMethodOptions const &) const {

	methods::EnergyMethodCreatorOP creator( new SAXSEnergyCreatorCEN );
	return methods::EnergyMethodOP( new SAXSEnergy(SAXSEnergy::cen_cfg_file_,
	    chemical::ChemicalManager::get_instance()->residue_type_set(core::chemical::CENTROID),
	    saxs_cen_score, creator ) );
}

/// c-tor
SAXSEnergy::SAXSEnergy(std::string & config_file,core::chemical::ResidueTypeSetCAP rsd_set_ap,
		ScoreType the_variant, methods::EnergyMethodCreatorOP creator) : WholeStructureEnergy( creator ) {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;
    
    core::chemical::ResidueTypeSetCOP rsd_set( rsd_set_ap );

    saxs_score_variant_ = the_variant;
    trSAXSEnergy.Warning << "SAXS score setup : " << the_variant << std::endl;
    init_ff(config_file);
    set_up_q();
    the_config_file_ = config_file;

    score_bias_ = option[score::saxs::min_score]();
    score_bias_ = std::pow(10,score_bias_);

    if(option[score::saxs::ref_spectrum].user()) {
	std::string file = option[score::saxs::ref_spectrum]();
	trSAXSEnergy << "Reading reference spectrum: " << file << std::endl;
	read_intensities(file,reference_intensities_);
    } else {
	if ( basic::options::option[in::file::native].user() ) {
	    core::pose::Pose reference_pose;
	    trSAXSEnergy << "Using "<<basic::options::option[in::file::native]()<<" as a reference for SAXS energy"<<std::endl;
	    //core::import_pose::pose_from_pdb(reference_pose, *rsd_set,basic::options::option[in::file::native]());
		 core::io::pdb::build_pose_from_pdb_as_is(reference_pose, *rsd_set, basic::options::option[in::file::native]());
	    compute_intensities(reference_pose,reference_intensities_);
	    trSAXSEnergy << "Calculated reference spectrum from a native: "
		    << basic::options::option[in::file::native]() << std::endl;
	}
    }
}

/// c-tor
SAXSEnergy::SAXSEnergy(const std::string & config_file,
		const utility::vector1<Real> & source_q,
		const utility::vector1<Real> & reference_spectrum,
		ScoreType score_variant,
		methods::EnergyMethodCreatorOP the_creator) :
		WholeStructureEnergy( the_creator ) {

    saxs_score_variant_ = score_variant;
    the_config_file_ = config_file;
    init_ff(config_file);
    set_up_q(source_q);
    fit_intensities(source_q,reference_spectrum,reference_intensities_);
}

void SAXSEnergy::set_up_q() {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    Real s_min_ = basic::options::option[score::saxs::q_min]();
    Real s_max_ = basic::options::option[score::saxs::q_max]();
    Real s_step_ = basic::options::option[score::saxs::q_step]();
    q_.clear();
    for(Real s=s_min_;s<=s_max_;s+=s_step_)
    	    q_.push_back(s);
    pose_intensities_.resize(q_.size());
    reference_intensities_.resize( q_.size() );

    trSAXSEnergy.Warning << "SAXS score q-set : " << q_.size() << " points from a file given at cmdline"<< std::endl;
}

void SAXSEnergy::set_up_q(const utility::vector1<Real> & source_q) {

    q_.resize( source_q.size() );
    for(Size i=1;i<=q_.size();i++)
      q_[i] = source_q[i];
    pose_intensities_.resize( source_q.size() );
    reference_intensities_.resize( source_q.size() );
    trSAXSEnergy.Warning << "SAXS score q-set : " << q_.size() << " points as a deep copy"<< std::endl;
}

void SAXSEnergy::init_ff(const std::string & config_file) {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    ff_manager_ = FormFactorManager::get_manager();
    if ( basic::options::option[score::saxs::custom_ff].user() ) {
        trSAXSEnergy << "Loading custom FF from "<<basic::options::option[score::saxs::custom_ff]()<<std::endl;
        ff_manager_->load_ff( basic::options::option[score::saxs::custom_ff]() );
    } else {
	ff_manager_->load_ff_from_db(basic::database::full_name( "scoring/score_functions/saxs/" + config_file ));
    }

    if_hydrogens_=false;
}


Real SAXSEnergy::total_energy(const pose::Pose & pose) const {

    if( q_.size() != pose_intensities_.size() )
	  pose_intensities_.resize( q_.size() );
    compute_intensities(pose,pose_intensities_);

    return compute_chi(pose_intensities_,reference_intensities_);
}

void SAXSEnergy::finalize_total_energy(pose::Pose & pose, ScoreFunction const &,
	EnergyMap & totals ) const {

//	PROF_START( basic::SAXS );
	if(q_.size() != pose_intensities_.size())
	  pose_intensities_.resize( q_.size() );
	compute_intensities(pose,pose_intensities_);
//	totals[ saxs_score_variant_ ] = compute_L1(pose_intensities,reference_intensities_);
	totals[ saxs_score_variant_ ] = compute_chi(pose_intensities_,reference_intensities_);
//	PROF_STOP( basic::SAXS );

}

void SAXSEnergy::read_spectrum(std::string & file_name,utility::vector1<Real> & q,utility::vector1<Real> & I) const {

	Real tX,tY;
	utility::io::izstream input(file_name.c_str());
	if ( !input.good() )
		utility_exit_with_message( "Unable to open reference spectrum file: " + file_name );
	std::string line;
        while( getline( input, line ) ) {
	    if ( line.substr(0,1) == "#" ) continue;
	    if ( line.length() < 4 ) continue;
	    std::istringstream line_stream( line );
	    line_stream >> tX >> tY;
	    q.push_back( tX );
	    I.push_back( tY );
	    trSAXSEnergy.Trace << "Reference SAXS data: "<<tX<<" "<<tY<<std::endl;
	}
}

void SAXSEnergy::fit_intensities(const utility::vector1<Real> & q,const utility::vector1<Real> & I,
					utility::vector1<Real> & result) const {

	Real minX,minY;
	Real maxX,maxY;

	minX = maxX = q[1];
	maxY = minY = I[1];
	for(Size i=2;i<=q.size();i++) {
		if(q[i] < minX) minX = q[i];
		if(q[i] > maxX) maxX = q[i];
		if(I[i] < minY) minY = I[i];
		if(I[i] > maxY) maxY = I[i];
	}

	Real delta = 0.1;
	numeric::interpolation::spline::SplineGenerator gen( minX-delta, minY-delta, 0, maxX+delta, maxY+delta, 0 );
	for (Size i = 1; i <= q.size(); ++i) {
		gen.add_known_value( q[i],I[i] );
	}
	utility::pointer::shared_ptr< numeric::interpolation::spline::Interpolator > spline_interpolator = gen.get_interpolator();

	if( result.size() != q_.size() ) {
	    result.clear();
	    result.resize( q_.size() );
	}
        Real dy = 0.0;
	for(Size i_s=1;i_s<=q_.size();++i_s) {
	    Real r;
	    spline_interpolator->interpolate(q_[i_s], r, dy);
	    result[i_s] = r;
	}
}

void SAXSEnergy::read_intensities(std::string & file_name, utility::vector1<Real> & result) const {

	utility::vector1<Real> x;
	utility::vector1<Real> y;
	read_spectrum(file_name,x,y);
	fit_intensities(x,y,result);
}


void SAXSEnergy::rehash_form_factors(const core::pose::Pose & pose) const {

    // ---------- Tabulate form factors
    ff_manager_->tabulate(q_);
    std::set<FormFactorOP> ff_set;
    // ---------- Find the unique set of atom types
    for ( Size i = 1; i <= pose.total_residue(); ++i ) {
	core::conformation::Residue resi = pose.residue(i);
        for ( Size m = 1; m <= resi.natoms(); ++m ) {
    	    trSAXSEnergy.Trace << "rehash: trying "<<resi.atom_type(m).name();
    	    if((! if_hydrogens_)&&( resi.atom_type(m).is_hydrogen())) {
	        trSAXSEnergy.Trace << " rejected hydrogen"<<std::endl;
	    	continue;
	    }
	    if( ! ff_manager_->is_known_atom( resi.atom_type(m).name() ) ) {
		trSAXSEnergy.Trace << " rejected unknown"<<std::endl;
		continue;
	    }
	    FormFactorOP fi =  ff_manager_->get_ff(resi.atom_type(m).name());
	    ff_set.insert( fi );
	}
    }

    // ---------- Repack the set of unique FFs to a vactor and hash their indexes in with a map
    std::set<FormFactorOP>::iterator it;
    ff_ops_.clear();
    ff_map_.clear();
    Size i = 1;
    for (it=ff_set.begin(); it!=ff_set.end(); it++) {
      ff_map_.insert( std::pair<FormFactorOP,Size>(*it,i) );
      ff_ops_.push_back( *it );
      i++;
    }

    // ---------- Create a matrix of distance histograms
    for(Size i=1;i<=ff_ops_.size();i++) {
      utility::vector1<DistanceHistogramOP> row;
      dhist_.push_back( row );
      for(Size j=1;j<=ff_ops_.size();j++) {
        dhist_[i].push_back( DistanceHistogramOP( new DistanceHistogram() ) );
      }
    }

    // ---------- Prepare atoms, residues and assign atom types
    for ( Size i = 1; i <= pose.total_residue(); ++i ) {
	core::conformation::Residue resi = pose.residue(i);
        for ( Size m = 1; m <= resi.natoms(); ++m ) {
    	    if((! if_hydrogens_)&&( resi.atom_type(m).is_hydrogen()))
	    	continue;
	    if( ! ff_manager_->is_known_atom( resi.atom_type(m).name() ) )
		continue;
	    r_ids_.push_back( i );
	    a_ids_.push_back( m );
	    FormFactorOP fi =  ff_manager_->get_ff(resi.atom_type(m).name());
	    atom_ff_types_.push_back( ff_map_[fi] );
	}
    }
    trSAXSEnergy.Debug << "Found "<<atom_ff_types_.size()<<" atoms suitable for SAXS computations"<<std::endl;
    compute_distance_histogram(pose);
    zero_ = compute_zero_intensity();
    trSAXSEnergy.Debug << "Zero-intensity is: "<<zero_<<std::endl;
}

void SAXSEnergy::compute_distance_histogram(const core::pose::Pose & pose) const {


	for(Size i=1;i<=ff_ops_.size();i++)
	  for(Size j=1;j<=ff_ops_.size();j++)
	    dhist_[i][j]->zeros();

	/*********** Compute distance histogram ************/
    	for ( Size i = 2; i <= r_ids_.size(); ++i ) {
		Size res_i = r_ids_[i];
		Size atm_i = a_ids_[i];
		Size ti = atom_ff_types_[i];
		numeric::xyzVector<Real> ai = pose.residue(res_i).xyz(atm_i);
    		for ( Size j = 1; j < i; ++j ) {
		    Size res_j = r_ids_[j];
		    Size atm_j = a_ids_[j];
		    Size tj = atom_ff_types_[j];
		    Real d = ai.distance( pose.residue(res_j).xyz(atm_j) );
		    if( tj <= ti )
		      dhist_[ ti ][ tj  ]->insert(d);
		    else
		      dhist_[ tj ][ ti  ]->insert(d);
		}
        }

}

void SAXSEnergy::compute_intensities(const core::pose::Pose & pose,utility::vector1<Real> & result) const {


	SinXOverX *sin_x_by_x_ = SinXOverX::get_instance();

	if(ff_ops_.size() == 0) rehash_form_factors(pose);

	compute_distance_histogram(pose);
	zero_ = compute_zero_intensity();

	for(Size i_s=1;i_s<=q_.size();++i_s) {

    	    Real sum(0);
	    Real val_s = q_[i_s];
    	    for ( Size i = 1; i <= ff_ops_.size(); ++i ) {
    		Real fi = ff_ops_[i]->get(i_s);
        	for ( Size j = 1; j <= i; ++j ) {
        	    Real fij = ff_ops_[j]->get(i_s) * fi;
        	    DistanceHistogramOP dh = dhist_[i][j];
		    for(Size i_r=1;i_r<=  dh->last_nonempty_bin();++i_r) {  // ---------- Go through histogram bins
			Size dn = dh->get( i_r );
			if( dn > 0 )
	    		  sum += fij * dn * sin_x_by_x_->evaluate( dh->distance( i_r ) * val_s );
	    	    }
    	        }
	    }

// 	    trSAXSEnergy.Trace <<std::setw(8)<<std::setprecision(4)<<val_s<<" "
//		<< std::setw(8)<<std::setprecision(4)<<sum<<std::endl;
            result[i_s] = sum;
        }
//	trSAXSEnergy.Debug << "SAXS score: "<<ff_ops_.size()<<" atoms suitable for SAXS computations"<<std::endl;
}

Real SAXSEnergy::compute_zero_intensity() const {

	Real sum(0);
	for(Size i=1;i<=ff_ops_.size();i++)
	  for(Size j=1;j<=ff_ops_.size();j++) {
	    sum += dhist_[i][j]->total() * ff_ops_[i]->get(1) * ff_ops_[j]->get(1);
	  }
        return 2.0 * sum;
}


Real SAXSEnergy::compute_L1(utility::vector1<Real> const & saxs_scored, utility::vector1<Real> const & saxs_reference) const {

	Real chi = -1.0;
        Real sum = 0;
        assert( saxs_scored.size() == saxs_reference.size() );
        for(Size i=1;i<=saxs_scored.size();++i)
    	    sum += saxs_reference[i] / saxs_scored[i];
        Real lambda = sum / ((Real) saxs_scored.size());

	trSAXSEnergy.Trace << "\nComputing SAXS energy:\n";
        trSAXSEnergy.Trace << "     q     pose  reference\n";
	for(Size i=1;i<=saxs_reference.size();++i) {
	    Real tmp = saxs_scored[i] * lambda - saxs_reference[i];
	    trSAXSEnergy.Trace << i<<" " << saxs_scored[i] * lambda <<" "<< saxs_reference[i] << "\n";
	    if( fabs(tmp) > chi )
        	chi = fabs(tmp);
	}
	trSAXSEnergy.Trace <<std::endl;

	trSAXSEnergy.Debug << "\nSAXS energy: " <<  chi << std::endl;
	return chi;
}

Real SAXSEnergy::compute_chi(utility::vector1<Real> const & saxs_scored, utility::vector1<Real> const & saxs_reference) const {

	Real chi = 0;
	//Real sum = 0;
        assert( saxs_scored.size() == saxs_reference.size() );

	trSAXSEnergy.Trace << "\nComputing SAXS energy:\n";
	trSAXSEnergy.Trace << "norm(0): "<<zero_<<"\n";
        trSAXSEnergy.Trace << "     q     pose  pose*lambda reference\n";
	for(Size i=1;i<=saxs_reference.size();++i) {
	    Real tmp = saxs_scored[i] - saxs_reference[i];
	    tmp /= zero_;
	    trSAXSEnergy.Trace << i<<" " << saxs_scored[i]/zero_ <<" "<< saxs_reference[i]/ zero_
		<< " " << (saxs_scored[i] / zero_ - saxs_reference[i]/ zero_) << "\n";
            chi += tmp*tmp / saxs_reference[i] * zero_;
	}
	trSAXSEnergy.Trace <<std::endl;
	chi /= ((Real) q_.size()) ;

//	Real energy = chi;
	Real energy = log10(chi + score_bias_);
	trSAXSEnergy.Debug << "\nSAXS chi2, energy, n_atoms: " <<  chi << " " << energy << " " << atom_ff_types_.size() << std::endl;
	return energy;
}


Real SAXSEnergy::compute_chi_with_fit(utility::vector1<Real> const & saxs_scored, utility::vector1<Real> const & saxs_reference) const {

	Real chi = 0;
        Real sum = 0;
        assert( saxs_scored.size() == saxs_reference.size() );
        for(Size i=1;i<=saxs_scored.size();++i)
    	    sum += saxs_reference[i] / saxs_scored[i];
        Real lambda = sum / ((Real) saxs_scored.size());

	trSAXSEnergy.Trace << "\nComputing SAXS energy:\n";
	trSAXSEnergy.Trace << "lambda:  "<<lambda<<"\n";
	trSAXSEnergy.Trace << "norm(0): "<<zero_<<"\n";
        trSAXSEnergy.Trace << "     q     pose  pose*lambda reference\n";
	for(Size i=1;i<=saxs_reference.size();++i) {
	    Real tmp = saxs_scored[i] * lambda / zero_ - saxs_reference[i] / zero_;
	    trSAXSEnergy.Trace << i<<" " << saxs_scored[i]<<" "<<saxs_scored[i] * lambda/ zero_ <<" "<< saxs_reference[i]/ zero_
		<< " " << (saxs_scored[i] * lambda / zero_ - saxs_reference[i]/ zero_) << "\n";
            chi += tmp*tmp;
	}
	trSAXSEnergy.Trace <<std::endl;
	chi /= ((Real) q_.size());

//	Real energy = chi;
	Real energy = log10(chi + score_bias_);
	trSAXSEnergy.Debug << "\nSAXS chi2, energy, n_atoms: " <<  chi << " " << energy << " " << atom_ff_types_.size() << std::endl;
	return energy;
}

core::Size
SAXSEnergy::version() const
{
	return 1; // Initial versioning
}


} // saxs
} // scoring
} // core

