// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/saxs/PDDFEnergy.cc
/// @brief  "Energy" based on a similarity of theoretical PDDF (pairwise distance distribution function)
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


// Unit headers
#include <core/scoring/saxs/FormFactorManager.hh>
#include <core/scoring/saxs/FormFactor.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergy.hh>
#include <protocols/scoring/methods/saxs/PDDFEnergyCreator.hh>

#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/WholeStructureEnergy.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/id/NamedAtomID.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>

#include <numeric/interpolation/spline/Interpolator.hh>
// AUTO-REMOVED #include <numeric/interpolation/spline/SplineGenerator.hh>

// Options
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

#include <core/chemical/ChemicalManager.hh>


// Utility headers
#include <basic/prof.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
#include <basic/Tracer.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

// C++
#include <iomanip>

#include <core/import_pose/import_pose.hh>
#include <core/kinematics/Jump.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace saxs {

static thread_local basic::Tracer trPDDFEnergy( "protocols.scoring.methods.saxs.PDDFEnergy" );

core::scoring::ScoreTypes PDDFEnergyCreator::score_types_for_method() const {
	core::scoring::ScoreTypes sts;
	sts.push_back( core::scoring::pddf_score );
	return sts;
}

core::scoring::methods::EnergyMethodOP PDDFEnergyCreator::create_energy_method( core::scoring::methods::EnergyMethodOptions const &) const {

	return core::scoring::methods::EnergyMethodOP( new PDDFEnergy() );
}

/// c-tors
PDDFEnergy::PDDFEnergy(utility::vector1<core::Real> const & reference_d,utility::vector1<core::Real> const & reference_pddf) : WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new PDDFEnergyCreator ) ) {

	if_fit_area_ = false;
	d_.clear();
	reference_pddf_.clear();
	for(core::Size i=1;i<=reference_pddf.size();i++) {
	    d_.push_back( reference_d[i] );
	    reference_pddf_.push_back(reference_pddf[i]);
	}
	pose_pddf_.resize( d_.size(), 0 );
	bin_size_ = (d_[ d_.size()  ] - d_[1] ) / (d_.size() - 1);
        min_bin_ = (core::Size)(d_[1] / bin_size_);
	max_bin_ = (core::Size)(d_[d_.size()] / bin_size_);
	if(max_bin_ > reference_pddf_.size())
	    max_bin_ = reference_pddf_.size();
}

PDDFEnergy::PDDFEnergy() : WholeStructureEnergy( core::scoring::methods::EnergyMethodCreatorOP( new PDDFEnergyCreator ) ) {

    using namespace basic::options;
    using namespace basic::options::OptionKeys;

    norm_ = 1000.0;
    if_hydrogens_ = false;
    if_fit_area_ = false;

/*
    ff_manager_ = FormFactorManager::get_manager();
    if ( basic::options::option[score::saxs::custom_ff].user() ) {
        trPDDFEnergy << "Loading custom FF from "<<basic::options::option[score::saxs::custom_ff]()<<std::endl;
        ff_manager_->load_ff( basic::options::option[score::saxs::custom_ff]() );
    } else {
	ff_manager_->load_ff(basic::database::full_name( "scoring/score_functions/saxs/ff-rosetta.cfg" ));
    }

*/
    if( basic::options::option[score::saxs::fit_pddf_area]() )
	if_fit_area_ = true;

    core::Real min_d = basic::options::option[score::saxs::d_min]();
    core::Real max_d = basic::options::option[score::saxs::d_max]();
    core::Real bin_d = basic::options::option[score::saxs::d_step]();

    if(option[score::saxs::ref_pddf].user()) {
	std::string file = option[score::saxs::ref_pddf]();

	trPDDFEnergy << "Reading reference PDDF: " << file << std::endl;
	read_pddf(file);
    min_bin_ = (core::Size)(min_d / bin_size_);   //These should be static casts not c style casts
	max_bin_ = (core::Size)(max_d / bin_size_);
	if(max_bin_ > reference_pddf_.size())
	    max_bin_ = reference_pddf_.size();
    } else {
	bin_size_ = bin_d;
	if(option[in::file::native].user()) {

	    core::pose::Pose reference_pose;
	    trPDDFEnergy << "Using "<<basic::options::option[in::file::native]()<<" as a reference"<<std::endl;
    	    core::chemical::ResidueTypeSetCOP rsd_set;
            rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
	                      basic::options::option[ in::file::residue_type_set ]());

	    core::import_pose::pose_from_pdb(reference_pose, *rsd_set,
		basic::options::option[in::file::native]());

	    core::util::switch_to_residue_type_set( reference_pose, core::chemical::CENTROID );

	    create_pddf( reference_pose,min_d,bin_d,max_d );
    	    min_bin_ = (core::Size)(min_d / bin_size_);  //These should be static casts not c style casts
	    max_bin_ = (core::Size)(max_d / bin_size_);
	} else {
	    utility_exit_with_message( "[ERROR] Unable to load reference PDDF: " );
	}
    }
}


void PDDFEnergy::finalize_total_energy(core::pose::Pose & pose, core::scoring::ScoreFunction const &,
	core::scoring::EnergyMap & totals ) const {

	PROF_START( basic::SAXS );
	compute_pddf_without_ff(pose);
//	totals[ pddf_score ] = compute_L1(pose_intensities,reference_intensities_);
	totals[ core::scoring::pddf_score ] = compute_chi(pose_pddf_,reference_pddf_);
	PROF_STOP( basic::SAXS );
}

core::Real PDDFEnergy::evaluate_pddf_energy(const core::pose::Pose & pose) const {

	compute_pddf_without_ff(pose);
	return compute_chi(pose_pddf_,reference_pddf_);
}


void PDDFEnergy::create_pddf(core::pose::Pose & ref_pose,core::Real d_min,core::Real d_step,core::Real d_max) {

	d_.clear();
	bin_size_ = d_step;
	for(core::Real d=d_min;d<=d_max;d+=d_step)
	    d_.push_back(d);
	pose_pddf_.resize(d_.size());
	reference_pddf_.resize(d_.size(),0.0);

	compute_pddf_without_ff(ref_pose);

	for(core::Size i=1;i<=pose_pddf_.size();i++)
	    reference_pddf_.push_back(pose_pddf_[i]);
}

void PDDFEnergy::read_pddf(std::string file_name) {

	d_.clear();
	reference_pddf_.clear();
	core::Real tX,tY;
	utility::io::izstream input(file_name.c_str());
	if ( !input.good() )
		utility_exit_with_message( "Unable to open reference spectrum file: " + file_name );
	std::string line;
    core::Real norm = 0.0;
    while( getline( input, line ) ) {
	    if ( line.substr(0,1) == "#" ) continue;
	    if ( line.length() < 4 ) continue;
	    std::istringstream line_stream( line );
	    line_stream >> tX >> tY;
	    d_.push_back( tX );
	    reference_pddf_.push_back( tY );
	    norm += tY;
	    trPDDFEnergy.Trace << "Reference PDDF data: "<<tX<<" "<<tY<<std::endl;
	}

	pose_pddf_.resize( d_.size(), 0 );
	bin_size_ = (d_[ d_.size()  ] - d_[1] ) / (d_.size() - 1);

//	  trPDDFEnergy.Debug << "Changing reference PDDF normalisation from: "<<norm<<" to "<<norm_<<std::endl;
//	  norm = norm_ / norm;
//	for(Size i=1;i<=reference_pddf_.size();i++)
//	    reference_pddf_[i] *= norm;
}

utility::vector1<core::Real> & PDDFEnergy::compute_pddf_without_ff(const core::pose::Pose & pose) const {

        a_ids_.clear();
        r_ids_.clear();
	for(core::Size i=1;i<= pose_pddf_.size();i++) {
	    pose_pddf_[i] = 0.0;
	}
	/*********** Rehashing atoms, preparig data structure ************/
        for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
	        core::conformation::Residue resi = pose.residue(i);
        	for ( core::Size m = 1; m <= resi.natoms(); ++m ) {
//        	    trPDDFEnergy.Trace << "rehash: trying "<<resi.atom_type(m).name();
    			if((! if_hydrogens_)&&( resi.atom_type(m).is_hydrogen())) {
//		        trPDDFEnergy.Trace << " rejected hydrogen"<<std::endl;
	    		continue;
	    	    }
		    a_ids_.push_back( m );
		    r_ids_.push_back( i );
		}
	}
	/*********** Compute distance matrix ************/
	trPDDFEnergy.Debug << "computing PDDF based on "<<r_ids_.size()<<" atoms"<<std::endl;
    	for ( core::Size i = 2; i <= r_ids_.size(); ++i ) {
		core::Size res_i = r_ids_[i];
		core::Size atm_i = a_ids_[i];
    		for ( core::Size j = 1; j < i; ++j ) {
		    core::Size res_j = r_ids_[j];
		    core::Size atm_j = a_ids_[j];

            	    core::Real d = pose.residue(res_i).xyz(atm_i).distance( pose.residue(res_j).xyz(atm_j) );
            	    core::Size bin = ((core::Size)(d/bin_size_))+1;
            	    if(bin>reference_pddf_.size())
            	      bin = reference_pddf_.size();
        	    pose_pddf_[ bin ] ++;
		}
    }

        return pose_pddf_;
}


utility::vector1<core::Real> & PDDFEnergy::compute_pddf(const core::pose::Pose & pose) const {

	/*********** Rehashing atoms, preparig data structure ************/
        a_ids_.clear();
        r_ids_.clear();
	if(factors_.size() == 0) {
            for ( core::Size i = 1; i <= pose.total_residue(); ++i ) {
	        core::conformation::Residue resi = pose.residue(i);
        	for ( core::Size m = 1; m <= resi.natoms(); ++m ) {
        	    trPDDFEnergy.Trace << "rehash: trying "<<resi.atom_type(m).name();
    			if((! if_hydrogens_)&&( resi.atom_type(m).is_hydrogen())) {
		        trPDDFEnergy.Trace << " rejected hydrogen"<<std::endl;
	    		continue;
	    	    }
		    if( ! ff_manager_->is_known_atom( resi.atom_type(m).name() ) ) {
			trPDDFEnergy.Trace << " rejected unknown"<<std::endl;
			continue;
		    }
		    core::scoring::saxs::FormFactorOP fi =  ff_manager_->get_ff(resi.atom_type(m).name());
		    utility::vector1<core::Real> row;
		    for(core::Size i_s=1;i_s<=d_.size();++i_s) {
			row.push_back(fi->ff(d_[i_s]));
		    }
		    factors_.push_back( row );
		    a_ids_.push_back( m );
		    r_ids_.push_back( i );
		    if( ( resi.atom_type(m).name()[0] == 'C' ) && ( resi.atom_type(m).name()[1] == 'E' ) && ( resi.atom_type(m).name()[2] == 'N' ) ) {
			is_glob_.push_back(true);
			trPDDFEnergy.Trace << " .. is glob ..";
		    }
		    else {
			is_glob_.push_back(false);
			trPDDFEnergy.Trace << " .. is atom ..";
		    }
		    trPDDFEnergy.Trace << " OK."<<std::endl;
		}
	    }
	    trPDDFEnergy.Trace << "Found "<<factors_.size()<<" atoms suitable for PDDF computations"<<std::endl;

	    dmatrix_.resize(factors_.size());
	    for(core::Size i=1;i<=dmatrix_.size();i++)
		dmatrix_[i].resize(dmatrix_.size());
	}

	/*********** Compute distance matrix ************/
    	for ( core::Size i = 2; i <= factors_.size(); ++i ) {
		core::Size res_i = r_ids_[i];
		core::Size atm_i = a_ids_[i];
    		for ( core::Size j = 1; j < i; ++j ) {
		    core::Size res_j = r_ids_[j];
		    core::Size atm_j = a_ids_[j];

            	    core::Real d = pose.residue(res_i).xyz(atm_i).distance( pose.residue(res_j).xyz(atm_j) );
        	    pose_pddf_[ ((core::Size)(d/bin_size_)) ] ++;
		}
        }
        return pose_pddf_;
}

core::Real PDDFEnergy::compute_L1(utility::vector1<core::Real> const & pddf_scored, utility::vector1<core::Real> const & pddf_reference) const {

	core::Real chi = -1.0;
        core::Real sum = 0;
        assert( pddf_scored.size() == pddf_reference.size() );
        for(core::Size i=1;i<=pddf_scored.size();++i)
    	    sum += pddf_reference[i] / pddf_scored[i];
        core::Real lambda = sum / ((core::Real) pddf_scored.size());

	trPDDFEnergy.Trace << "Computing PDDF energy:\n";
        trPDDFEnergy.Trace << "     q     pose  reference\n";
	for(core::Size i=1;i<=pddf_reference.size();++i) {
	    core::Real tmp = pddf_scored[i] * lambda - pddf_reference[i];
	    if( fabs(tmp) > chi )
        	chi = fabs(tmp);
	}

	trPDDFEnergy.Debug << "\nPDDF energy: " <<  chi << std::endl;
	return chi;
}


core::Real PDDFEnergy::compute_chi(utility::vector1<core::Real> const & pddf_scored, utility::vector1<core::Real> const & pddf_reference) const {

	core::Real chi = 0;
        assert( pddf_scored.size() == pddf_reference.size() );
	core::Real lambda = 1.0;
        if ( if_fit_area_ ) {
    	    core::Real sum_scored = 0;
    	    core::Real sum_ref = 0;
    	    for(core::Size i=min_bin_;i<=max_bin_;++i) {
    		sum_scored += pddf_scored[i];
		sum_ref += pddf_reference[i];
	    }
    	    lambda = sum_ref / sum_scored;
        }

	trPDDFEnergy.Trace << "Computing PDDF energy:\n";
        trPDDFEnergy.Trace << "     r_i   pose  reference  avg   delta-chi  sum-chi\n";
        core::Real nonzero = 0.0;
	for(core::Size i=min_bin_;i<=max_bin_;++i) {
	    if( (pddf_scored[i]<0.00001) && ( pddf_reference[i]<0.000001) )
		continue;
	    nonzero++;
	    core::Real avg = (pddf_scored[i] * lambda + pddf_reference[i]) / 2.0;
	    core::Real tmp = pddf_scored[i] * lambda - pddf_reference[i];
            chi += tmp*tmp/avg;
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<<d_[i]<<" ";
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<<pddf_scored[i]*lambda<<" ";
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<<pddf_reference[i]<<" ";
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<< avg<<" ";
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<< (tmp*tmp)/avg<<" ";
	    trPDDFEnergy.Trace << std::setw(9)<<std::setprecision(4)<<chi<<"\n";
	}
	chi /= nonzero;

	trPDDFEnergy.Debug << "\nPDDF energy: " <<  chi << std::endl;
	return chi;
}
core::Size
PDDFEnergy::version() const
{
	return 1; // Initial versioning
}


} // saxs
} // scoring
} // methods
} // protocols
