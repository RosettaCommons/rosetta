// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/EnvPairPotential.cc
/// @brief  Membrane Potential
/// @author Bjorn Wallner


// Unit headers
#include <core/scoring/ProQPotential.hh>
#include <core/scoring/ProQPotential.fwd.hh>
//#include <core/scoring/MembraneTopology.hh>

// Package headers

//#include <core/scoring/EnvPairPotential.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <basic/database/open.hh>
#include <core/pose/Pose.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>


//symmetry
#include <core/pose/symmetry/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>

//options
#include <basic/options/option.hh>
#include <basic/options/keys/ProQ.OptionKeys.gen.hh>
#include <basic/options/keys/membrane.OptionKeys.gen.hh>

// Utility headers
#include <utility/io/izstream.hh>
#include <core/types.hh>


#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>
#include <numeric/conversions.hh>

#include <basic/Tracer.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <string>

namespace core {
namespace scoring {
static THREAD_LOCAL basic::Tracer TR( "core.scoring.ProQPotential" );

ProQPotential::ProQPotential()
{
	//  load the data
	// Size const num_models( 5 );
	//Size const num_features( 260 );
	Size const num_head( 11 );
	//Size const num_features_=260;
	//Size const num_features_proq2_=174;


	//Size num_models_=5;

	cross_val_=false;


	if ( basic::options::option[ basic::options::OptionKeys::ProQ::svmmodel].user() ) {
		svm_model_=basic::options::option[basic::options::OptionKeys::ProQ::svmmodel]();
		cross_val_=true;
	}


	std::string tag,line;
	//std::string aa;
	Size pos;
	// chemical::AA aa;
	//
	linear_weights_.dimension(num_models_,num_features_proqm_)=0;
	b_.dimension(num_models_);
	utility::io::izstream stream;
	for ( Size i=1; i<= num_models_; ++i ) {
		std::ostringstream i_stream;
		i_stream << i;;
		//std::cout << file << "\n";
		//utility::io::izstream stream;
		std::string file("scoring/ProQ/ProQM_model." + i_stream.str() + ".linear_retrained_for_rosetta");
		//std::cout << "file: " << file << std::endl;
		//std::string file("scoring/ProQ/ProQM_model." + i_stream.str() + ".linear");
		basic::database::open( stream, file);
		for ( Size j=1; j<=num_features_proqm_+num_head; ++j ) {
			getline( stream, line );
			std::istringstream l(line);
			//std::cout << tag << " " << aa << std::endl;
			if ( j==1 ) {
				l >> tag;
				if ( l.fail() || tag != "SVM-light"  ) utility_exit_with_message("bad format for " + file);
			}
			if ( j==num_head ) {
				l >> b_(i);
			}
			if ( j>num_head ) {
				l >> tag;
				l >> pos;
				l >> linear_weights_(i,pos);
				//TR.Debug << j << " " << pos << " " << linear_weights_(i,pos) << " " << b_(i) << std::endl;
			}
		}
		//std::cout << "Before outputing features..." << std::endl;
		//for(Size j=1;j<=num_features_;++j) {
		// std::cout << j << " " << linear_weights_(i,j) << " " << b_(i) << std::endl;
		//}
	}

	linear_weights_proq2_.dimension(num_models_,num_features_proq2_);
	b_proq2_.dimension(num_models_);
	//utility::io::izstream stream;
	for ( Size i=1; i<= num_models_; ++i ) {
		std::ostringstream i_stream;
		i_stream << i;;
		//std::cout << file << "\n";
		//utility::io::izstream stream;
		std::string file("scoring/ProQ/ProQ2_model." + i_stream.str() + ".linear_retrained_for_rosetta"); //_retrained_for_rosetta");
		basic::database::open( stream, file);
		for ( Size j=1; j<=num_features_proq2_+num_head; ++j ) {
			getline( stream, line );
			std::istringstream l(line);
			//std::cout << tag << " " << aa << std::endl;
			if ( j==1 ) {
				l >> tag;
				if ( l.fail() || tag != "SVM-light"  ) utility_exit_with_message("bad format for " + file);
			}
			if ( j==num_head ) {
				l >> b_proq2_(i);
			}
			if ( j>num_head ) {
				l >> tag;
				l >> pos;
				l >> linear_weights_proq2_(i,pos);
			}
		}
		/* for(Size j=1;j<=num_features_proq2_;++j) {
		TR.Debug << j << " " << linear_weights_proq2_(i,j) << " " << b_proq2_(i) << std::endl;
		}
		*/
	}
}


void
ProQPotential::setup_for_scoring( pose::Pose & ) const
{}

void
ProQPotential::score(pose::Pose & pose,
	ObjexxFCL::FArray2D< Real > & feature_vector,
	ObjexxFCL::FArray1D< Real > & score,
	bool ProQ2) const
{

	Size nres=pose.total_residue();
	Size n=1; //if position in sequence.
	if ( !ProQ2 ) {
		for ( Size l=1; l<=nres; ++l ) {
			if ( pose.residue(l).type().aa() == core::chemical::aa_vrt ) continue;
			//if(pose.residue(l).is_virtual_residue()) continue;
			score(n)=0;
			for ( Size i=1; i<=num_models_; ++i ) {
				if ( cross_val_ && i!=svm_model_ ) continue;
				Real pred(-b_(i));
				for ( Size k=1; k<=num_features_proqm_; ++k ) {
					//Real factor=feature_vector(l,k)*linear_weights_(i,k);
					//std::cout << "SVMcalc: " << l << " " << i << " " << k << " " << feature_vector(l,k) << " " << linear_weights_(i,k) << " " << pred << std::endl;
					pred+=feature_vector(l,k)*linear_weights_(i,k);
				}
				score(n)+=pred;
			}
			if ( !cross_val_ ) {
				score(n)/=num_models_;
			}
			n++;
		}
	} else { //ProQ2
		for ( Size l=1; l<=nres; ++l ) {
			//if(pose.residue(l).is_virtual_residue()) continue; //This was not true for some VRT residue !
			//std::cout << l << " " << pose.residue(l).type().aa() << std::endl;
			if ( pose.residue(l).type().aa() == core::chemical::aa_vrt ) continue;
			//std::cout << "PASSED: " << l << std::endl;
			//if(pose.residue(l).name().compare("VRT")==0) continue;
			//std::cout << pose.residue(l) << "dimension " << dim1 << " " << dim2 << "l: " << l << " n: " << l << " virtual: " << pose.residue(l).is_virtual_residue() << " virt_from_name: " << pose.residue(l).name().compare("VRT") << " res: " << pose.residue(l).name3() << " " << pose.residue(l).name() << std::endl;
			score(n)=0;
			for ( Size i=1; i<=num_models_; ++i ) {
				if ( cross_val_ && i!=svm_model_ ) continue;
				Real pred(-b_proq2_(i));
				for ( Size k=1; k<=num_features_proq2_; ++k ) {
					//Real factor=feature_vector(l,k)*linear_weights_proq2_(i,k);
					//std::cout << "SVMcalc: " << l << " " << i << " " << k << " " << feature_vector(l,k) << " " << linear_weights_(i,k) << " " << pred << std::endl;
					pred+=feature_vector(l,k)*linear_weights_proq2_(i,k);
				}
				score(n)+=pred;
			}
			if ( !cross_val_ ) {
				score(n)/=num_models_;
			}
			n++;
		}

	}

	//std::cout << "IN SCORE: " << b_.size() << " " << feature_vector.size() << std::endl;
	//score/=num_models_; //#0.750107;
}

Size
ProQPotential::num_features_proqm() const {
	return num_features_proqm_;
}

Size
ProQPotential::num_features_proq2() const {
	return num_features_proq2_;
}

}
}

