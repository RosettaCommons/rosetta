// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/ProQ_Energy.cc
/// @brief  ProQ energy function definition.
/// @author Björn Wallner


// Unit headers
#include <core/scoring/methods/ProQ_Energy.hh>
#include <core/scoring/methods/ProQ_EnergyCreator.hh>

// Package headers
#include <core/scoring/sasa.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/TwelveANeighborGraph.hh>
#include <core/scoring/ContextGraphTypes.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/datacache/CacheableString.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/AtomType.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/ProQ.OptionKeys.gen.hh>
#include <basic/Tracer.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/statistics/functions.hh>
#include <numeric/util.hh>
//#include <math.h>

//#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>


namespace core {
namespace scoring {
namespace methods {

static THREAD_LOCAL basic::Tracer TR( "core.scoring.methods.ProQ_Energy.cc" );

/// @details This must return a fresh instance of the ProQ_Energy class,
/// never an instance already in use
methods::EnergyMethodOP
ProQ_EnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new ProQ_Energy );
}

ScoreTypes
ProQ_EnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( ProQ ); // water (ProQ2)
	sts.push_back( ProQM ); // membrane (ProQM)
	// The only reason to seperate the two terms is to avoid mixup of the two.
	// One of them will always be zero, and it will be clear just from the scorefile which version was run.
	return sts;
}

/// c-tor
ProQ_Energy::ProQ_Energy() :
	parent( methods::EnergyMethodCreatorOP( new ProQ_EnergyCreator ) ),
	potential_( ScoringManager::get_instance()->get_ProQPotential() )
{
	initialize();
}

ProQ_Energy::ProQ_Energy( ProQ_Energy const & src) :
	parent( src ),
	potential_( src.potential_ )
{
	nres_=src.nres_;
	entropy_=src.entropy_;
	topology_=src.topology_;
	prob_profile_=src.prob_profile_;
	scaled_logodds_profile_=src.scaled_logodds_profile_;
	ss_pred_=src.ss_pred_;
	ss1_=src.ss1_;
	z_pred_=src.z_pred_;
	rsa_pred_=src.rsa_pred_;
	rsa_class_pred_=src.rsa_class_pred_;
	all_inputs_ProQM_=src.all_inputs_ProQM_;
	all_inputs_ProQ2_=src.all_inputs_ProQ2_;
}


/// clone
EnergyMethodOP
ProQ_Energy::clone() const
{
	//std::cout << "Cloning ProQ... " << nres_ << "\n";
	return EnergyMethodOP( new ProQ_Energy ( *this ) );
}

/////////////////////////////////////////////////////////////////////////////
// scoring
/////////////////////////////////////////////////////////////////////////////
void
ProQ_Energy::setup_for_scoring( pose::Pose & pose, ScoreFunction const & ) const
{
	//std::cout << "SecondaryStructureEnergy::setup_for_scoring" << std::endl;
	pose.update_residue_neighbors();
	potential_.setup_for_scoring( pose );
}

void
ProQ_Energy::finalize_total_energy(
	pose::Pose & pose,
	ScoreFunction const & scorefxn,
	EnergyMap & totals
) const {
	using namespace conformation;
	Size nres2=pose.size();
	Size nres=nres_;

	TR.Debug << "NRES: " << nres_ <<  " " << nres2 << " " << entropy_.size() << std::endl;
	const Size num_features_ProQM(potential_.num_features_proqm());
	const Size num_features_ProQ2(potential_.num_features_proq2());
	totals[ ProQ ] = 0;
	totals[ ProQM ] = 0;
	Real normalizing_factor=basic::options::option[ basic::options::OptionKeys::ProQ::normalize ]();

	if ( scorefxn.get_weight(ProQM)!=0 && all_inputs_ProQM_ ) {
		TR.Debug << "ProQM: " << std::endl;
		ObjexxFCL::FArray2D< Real > feature_vector(nres,num_features_ProQM);
		ObjexxFCL::FArray1D< Real > proq(nres);
		feature_vector=0;
		TR.Debug << "ProQM features: " << std::endl;
		calculate_feature_vector(pose,feature_vector);
		TR.Debug << "ProQM potential: " << std::endl;
		potential_.score( pose,feature_vector, proq, false);
		Real proq_mean=0;
		for ( Size i=1; i<=nres; ++i ) {
			TR.Debug << "ProQM: " << i << " " << proq(i) << std::endl;
			proq_mean+=proq(i);
		}
		totals[ ProQM ] = proq_mean/normalizing_factor;
		if ( basic::options::option[ basic::options::OptionKeys::ProQ::output_local_prediction ]() ) {
			output_local_prediction(pose,proq,"ProQM");
		}
	}
	if ( scorefxn.get_weight(ProQ) !=0 && all_inputs_ProQ2_ ) {
		ObjexxFCL::FArray2D< Real > feature_vector(nres,num_features_ProQ2);
		ObjexxFCL::FArray1D< Real > proq(nres);
		feature_vector=0;
		TR.Debug << "ProQ features: " << std::endl;
		calculate_feature_vector_proq2(pose,feature_vector);
		//std::cout << "Done feature vector" << std::endl;
		TR.Debug << "ProQ potential: " << std::endl;
		potential_.score( pose,feature_vector, proq,true);
		Real proq_mean=0;
		for ( Size i=1; i<=nres; ++i ) {
			TR.Debug << "ProQ: " << i << " " << proq(i) << std::endl;
			//std::cout << "ProQ: " << i << " " << proq(i) << std::endl;
			proq_mean+=proq(i);
		}
		totals[ ProQ ] = proq_mean/normalizing_factor;
		if ( basic::options::option[ basic::options::OptionKeys::ProQ::output_local_prediction ]() ) {
			output_local_prediction(pose,proq,"ProQ2");
		}
	}
} // finalize_total_energy


/// @brief ProQ_Energy distance cutoff
Distance
ProQ_Energy::atomic_interaction_cutoff() const
{
	return 6.0; /// now subtracted off 6.0 from cutoffs in centroid params files
}

core::Size
ProQ_Energy::version() const
{
	return 1; // Initial versioning
}

void
ProQ_Energy::initialize() {
	std::cout << "Reading all sequence specific stuff..." << std::endl;

	std::string basename("ProQM/idealized_test.orig/1e12A_0001.pdb");
	if ( basic::options::option[ basic::options::OptionKeys::ProQ::basename ].user() ) {
		basename=basic::options::option[ basic::options::OptionKeys::ProQ::basename]();
	} else {
		utility_exit_with_message( "Need ProQ::basename to read in the seq. specific stuff...");
	}

	std::string ss2file(basename + ".ss2");
	std::string profile(basename + ".psi");
	std::string mtxfile(basename + ".mtx");
	Size tmp=read_profiles_and_entropy(profile,mtxfile); //this will also initialize nres_ needs to called first, stupied... yes..
	TR.Debug << "AFTER MTX READ: " << nres_ << " " << tmp << " " << entropy_.size() << std::endl;
	read_ss2(ss2file);


	if ( basic::options::option[ basic::options::OptionKeys::ProQ::membrane ]() ) {
		std::string spanfile(basename + ".span");
		// std::string chkfile(basename + ".chk");
		std::string zpredfile(basename + ".zpred");
		std::string mpSAfile(basename + ".mpSA");
		topology_.initialize(spanfile);

		read_zpred(zpredfile);
		read_mpSA(mpSAfile);
		all_inputs_ProQM_=true;
	} else {
		std::string accfile(basename + ".acc");
		read_acc(accfile);
		all_inputs_ProQ2_=true;
	}
	// utility::vector1< Real > row(profile_.prof_row(1));
	//for(Size i=1;i<=20;++i) {
	// std::cout << "check: " << i << " " << row[i] << std::endl;
	//}
}

void
ProQ_Energy::output_local_prediction(pose::Pose & pose,
	ObjexxFCL::FArray1D< Real > & proq,
	std::string name) const {
	std::string outfile(basic::options::option[ basic::options::OptionKeys::ProQ::prefix]());
	using namespace core::pose::datacache;
	std::string tag( "empty_tag" );
	if ( pose.data().has( CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
		tag = static_cast< basic::datacache::CacheableString const & >
			( pose.data().get( CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
	}
	std::string gz("");
	if ( basic::options::option[ basic::options::OptionKeys::ProQ::use_gzip ]() ) {
		gz=".gz";
	}
	outfile+=name+"."+tag+gz;
	std::cout << "Outfile: " << outfile << std::endl;
	utility::io::ozstream output(outfile);
	if ( !output.good() ) {
		std::cerr << "Could not make " + outfile << std::endl;
		exit(1);
	}
	for ( Size i=1; i<=nres_; ++i ) {
		output << i << " " << proq(i) << std::endl;
	}
	output.close();
}

void
ProQ_Energy::read_ss2(std::string ss2file) {
	utility::io::izstream stream;

	std::string line;

	std::string buff_str;
	std::cout << "Reading ss2 from " << ss2file << std::endl;
	ss_pred_.dimension(nres_,3);
	ss1_.dimension(nres_);
	stream.open(ss2file);
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + ss2file);
	}

	if ( stream.get() == '#' ) {
		//std::cout << "ss2 format with header: "<< std::endl;
		getline(stream,line);
		//std::cout << "#" << line << std::endl;
		getline(stream,line);
		//std::cout << line << std::endl;
	}
	//else {
	// //std::cout << "ss2 format without header"<< std::endl;
	//}
	for ( Size i=1; i<=nres_; ++i ) {
		getline(stream,line);
		std::istringstream l(line);
		l >>  buff_str >> buff_str >> ss1_(i);
		l >> ss_pred_(i,1) >> ss_pred_(i,2) >> ss_pred_(i,3); // 1.coil,2.helix,3.sheet
		if ( ss1_(i) == 'C' ) {
			ss1_(i)='L'; // This assignment does not always correspond to the highest probability, e.g. HLH, will be HHH
		}
		/*if(ss_pred_(i,2) > ss_pred_(i,1)
		&& ss_pred_(i,2) > ss_pred_(i,3))
		ss='H';
		if(ss_pred_(i,3) > ss_pred_(i,1)
		&& ss_pred_(i,3) > ss_pred_(i,2))
		ss='E';
		std::cout << i << " " << ss_temp << " " << ss << " " << ss_pred_(i,1) << " " <<  ss_pred_(i,2) << " " << ss_pred_(i,3) << std::endl;

		ss1_(i)=ss_temp;
		*/
	}
	//std::exit(1);
	//for(Size i=1;i<=nres_;++i) {
	// std::cout << ss_pred_(i,1) << " " <<  ss_pred_(i,2) << " " << ss_pred_(i,3) << std::endl;
	//}
}
void
ProQ_Energy::read_acc(std::string accfile) {
	std::cout << "Read acc " << accfile << std::endl;
	utility::io::izstream stream;
	std::string line;
	rsa_class_pred_.dimension(nres_);
	stream.open(accfile);
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + accfile);
	}
	/*for(Size i=1;i<=5;++i) { //skip first 5 lines
	getline(stream,line);
	}*/
	getline(stream,line);
	std::istringstream l(line);
	for ( Size i=1; i<=nres_; ++i ) {
		l >>  rsa_class_pred_(i);
		runtime_assert(rsa_class_pred_(i) == 'e' ||
			rsa_class_pred_(i) == 'b');
	}
	/*
	for(Size i=1;i<=nres_;++i) {
	std::cout << "ACC: " << i << " " << rsa_class_pred_(i) << std::endl;
	}*/
}
void
ProQ_Energy::read_mpSA(std::string mpSAfile) {
	std::cout << "Read mpSA " << mpSAfile << std::endl;
	utility::io::izstream stream;
	std::string line;
	std::string buff_str;
	rsa_pred_.dimension(nres_);
	rsa_class_pred_.dimension(nres_);
	stream.open(mpSAfile);
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + mpSAfile);
	}
	/*for(Size i=1;i<=5;++i) { //skip first 5 lines
	getline(stream,line);
	}*/
	for ( Size i=1; i<=nres_; ++i ) {
		getline(stream,line);
		std::istringstream l(line);
		//l >>  buff_str >> rsa_class_(i) >> rsa_pred_(i);
		l >>  buff_str >> buff_str >> rsa_pred_(i);
		if ( rsa_pred_(i) > 25 ) {
			rsa_class_pred_(i)='e';
		} else {
			rsa_class_pred_(i)='b';
		}
	}
	/*for(Size i=1;i<=nres_;++i) {
	std::cout << i << " " << rsa_pred_(i) << std::endl;
	}*/
}

void
ProQ_Energy::read_zpred(std::string zpredfile) {
	std::cout << "Reading zpred from " << zpredfile << std::endl;
	utility::io::izstream stream;
	std::string line;
	z_pred_.dimension(nres_);
	stream.open(zpredfile);
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + zpredfile);
	}
	for ( Size i=1; i<=nres_; ++i ) {
		getline(stream,line);
		std::istringstream l(line);
		l >>  z_pred_(i);
	}
	//for(Size i=1;i<=nres_;++i) {
	// std::cout << z_pred_(i) << std::endl;
	//}
}


Size
ProQ_Energy::read_profiles_and_entropy(std::string profile,
	std::string mtxfile) {

	std::cout << "Reading profile and entropy from " << profile << std::endl;
	utility::io::izstream stream;
	std::string line;
	int buff_int;
	std::string buff_str;

	//initialize
	Size nres;
	stream.open(mtxfile);
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + mtxfile);
	}
	getline(stream,line);
	std::istringstream l(line);
	l >> nres_;
	nres=nres_;

	TR.Debug << "NRES from  " << mtxfile << " : " << nres_ <<std::endl;
	//initialize
	prob_profile_.dimension(nres_,20);
	scaled_logodds_profile_.dimension(nres_,20);
	entropy_.dimension(nres_);

	//Skip next 13 lines
	for ( Size i=1; i<=13; ++i ) {
		getline(stream,line);
	}
	for ( Size i=1; i<=nres_; ++i ) {
		getline(stream,line);
		std::istringstream l(line);
		// std::cout << line << std::endl;
		l >> buff_int;
		l >> buff_int;

		scaled_logodds_profile_(i,1)=(float)buff_int/100;
		//std::cout << i << " " <<scaled_logodds_profile_(i,1) <<  " " ;
		l >> buff_int;
		for ( Size j=2; j<=19; ++j ) {
			l >> buff_int;
			scaled_logodds_profile_(i,j)=(float)buff_int/100;
			//std::cout << scaled_logodds_profile_(i,j) <<  " " ;
		}
		l >> buff_int;
		l >> buff_int;
		scaled_logodds_profile_(i,20)=(float)buff_int/100;
		// std::cout << scaled_logodds_profile_(i,20) <<  std::endl;
	}
	stream.close();
	for ( Size i=1; i<=nres_; ++i ) {
		//  std::cout << "LOGODDS "  << i << " ";
		for ( Size j=1; j<=20; ++j ) {
			// std::cout << scaled_logodds_profile_(i,j) << " -> " ;
			scaled_logodds_profile_(i,j)=1/(1+exp((-1)*scaled_logodds_profile_(i,j)));
			// std::cout << scaled_logodds_profile_(i,j) << std::endl;
		}
	}
	stream.close();
	stream.open(profile); //entropy and prob_profile
	if ( stream.fail() ) {
		utility_exit_with_message( "ERROR: Unable to open: " + profile);
	}
	getline(stream,line);
	getline(stream,line);
	getline(stream,line);

	for ( Size i=1; i<=nres_; ++i ) {
		getline(stream,line);
		std::istringstream l(line);
		for ( Size j=1; j<=22; ++j ) { //low precision log_odds profile skipped
			l >> buff_str;
		}
		//std::cout << line << std::endl;
		for ( Size j=1; j<=20; ++j ) { //prob_matrix
			l >> buff_int;
			prob_profile_(i,j)=(float)buff_int/100;
			//std::cout << prob_profile_(i,j) << " ";
		}
		l >> entropy_(i);
		// std::cout << line << std::endl;
		// std::cout << entropy_(i) << std::endl;
	}
	stream.close();
	return nres;
}
void
ProQ_Energy::calculate_feature_vector_proq2(pose::Pose & pose,
	ObjexxFCL::FArray2D< Real > & feature_vector) const{

	//std::cout << "Calculating Feature vector ProQ2" << std::endl;
	const int nres=(int)nres_;
	//Calculate secondary structure
	ObjexxFCL::FArray2D< Real > feature_vector_orig(nres,174);


	dssp::Dssp ss(pose);
	ss.dssp_reduced();

	utility::vector1<Real> rsd_sasa_rel(nres,0.0);
	calc_per_atom_sasa_sc( pose, rsd_sasa_rel, true /*normalize*/ );
	//for(int i=1;i<=nres_;++i) {
	// std::cout << "RSA: " << rsd_sasa_rel[i] << std::endl;
	//}

	//ATOM 91
	atom_feature(pose,feature_vector,1,23);
	//RES 21
	res_feature(pose,feature_vector,92,23);
	//SURF 24 4x6
	surf_feature(pose,rsd_sasa_rel,feature_vector,113,23);
	gss_sc_feature(pose,ss,feature_vector,154);
	grsa_sc_feature(pose,rsd_sasa_rel,feature_vector,172);

	for ( Size i=1; i<=nres_; ++i ) {

		//stride 33 3x11 (dssp) binary [helix -5, sheet -5, coil -5] .. [helix 5, sheet 5, coil 5]
		stride_feature(pose,ss,feature_vector,i,137,5); // win 5
		//ss 1 prob for the specific ss.
		ss_feature(pose,ss,feature_vector,i,152,1); //win 1
		// ss_sc 1 Q3 over a 21 residue window
		ss_sc_feature(pose,ss,feature_vector,i,153,21); //win 21, same as default ProQM

		entropy_feature(pose,feature_vector,i,155,3); //win 3
		// profile 3x20, -1..1
		//profile_feature(pose,feature_vector,i,196);
		// rsa 1
		rsa_feature(pose,feature_vector,rsd_sasa_rel,i,158,13);
		//prsa_feature(pose,feature_vector,i,257);
		// rsa_sc 1
		rsa_sc_feature(pose,feature_vector,rsd_sasa_rel,i,171,21);
		// terminiN 1 + // terminiC 1
		termini_feature(pose,feature_vector,i,173,5);
	}

	if ( basic::options::option[ basic::options::OptionKeys::ProQ::output_feature_vector ]() ) {
		using namespace core::pose::datacache;
		std::string tag( "empty_tag" );
		if ( pose.data().has( CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
			tag =
				static_cast< basic::datacache::CacheableString const & >
				( pose.data().get( CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
		}
		//std::cout << "SVM START " << tag << std::endl;
		for ( Size i=1; i<=nres_; ++i ) {
			std::cout << "SVM " << tag << " ";
			//std::cout << "SVM 0 ";
			for ( Size j=1; j<=174; ++j ) {
				std::cout << j << ":" << feature_vector(i,j) << " ";
			}
			std::cout << std::endl;
		}
	}
}

void
ProQ_Energy::calculate_feature_vector(pose::Pose & pose,
	ObjexxFCL::FArray2D< Real > & feature_vector) const{

	//std::cout << "Calculating Feature vector for ProQM" << std::endl;
	//const int nres=(int)pose.size();
	const int nres=(int)nres_;
	//Calculate secondary structure
	ObjexxFCL::FArray2D< Real > feature_vector_orig(nres,260);
	dssp::Dssp ss(pose);
	ss.dssp_reduced();

	utility::vector1<Real> rsd_sasa_rel(nres,0.0);
	calc_per_atom_sasa_sc( pose, rsd_sasa_rel, true /*normalize*/ );
	bool use_naccess=false; //just to verify with original method, since the sasa routine is not exactly like naccess.
	if ( use_naccess ) {
		utility::io::izstream stream;
		std::string line;
		std::string buff_str;
		stream.open("ProQM/idealized_test.orig/1e12A_0001.pdb.rsa");
		if ( stream.fail() ) {
			utility_exit_with_message( "ERROR: Unable to open: ProQM/idealized_test.orig/1e12A_0001.pdb.rsa");
		}
		getline(stream,line);
		getline(stream,line);
		getline(stream,line);
		getline(stream,line);

		for ( int i=1; i<=nres; ++i ) {
			getline(stream,line);
			std::istringstream l(line);
			l >> buff_str >> buff_str >> buff_str ;
			l >> buff_str >> buff_str ;
			l >> buff_str >> rsd_sasa_rel[i];
			//  rsd_sasa_rel[i]=rsa_sasa_rel[i];
			//std::cout << line << std::endl << rsa_sasa_all[i] << " " << rsa_sasa_all[i] << " " << rsa_sasa[i] << std::endl;
		}
	}

	ObjexxFCL::FArray1D< Real > Z;
	calculateZ(pose,Z);


	// for(Size i=1;i<=pose.size();++i) {
	// std::cout << "DSSP " << i << " " << ss.get_dssp_secstruct(i) << std::endl;
	//}

	//ATOM 91
	atom_feature(pose,feature_vector,1);
	//RES 21
	res_feature(pose,feature_vector,92);

	//SURF 24 4x6
	surf_feature(pose,rsd_sasa_rel,feature_vector,113);

	for ( Size i=1; i<=nres_; ++i ) {


		//stride 33 3x11 (dssp) binary [helix -5, sheet -5, coil -5] .. [helix 5, sheet 5, coil 5]
		stride_feature(pose,ss,feature_vector,i,137);
		//ss 1 prob for the specific ss.
		ss_feature(pose,ss,feature_vector,i,170);
		// ss_sc 1 Q3 over a 21 residue window
		ss_sc_feature(pose,ss,feature_vector,i,171);
		// 11 topology -5..5 /M=1 nonM=0
		topology_feature(pose,feature_vector,i,172);
		// zpred normalized (z-5)/20, 1
		zpred_feature(pose,feature_vector,i,183);
		// z normalized (z-5)/20
		z_feature(pose,Z,feature_vector,i,184);
		// entropy 11, -5..5
		entropy_feature(pose,feature_vector,i,185);
		// profile 3x20, -1..1
		profile_feature(pose,feature_vector,i,196);
		// rsa 1
		rsa_feature(pose,feature_vector,rsd_sasa_rel,i,256);
		// prsa 1
		prsa_feature(pose,feature_vector,i,257);
		// rsa_sc 1
		rsa_sc_feature(pose,feature_vector,rsd_sasa_rel,i,258);
		// terminiN 1 + // terminiC 1
		termini_feature(pose,feature_vector,i,259);
	}

	if ( basic::options::option[ basic::options::OptionKeys::ProQ::output_feature_vector ]() ) {
		using namespace core::pose::datacache;
		std::string tag( "empty_tag" );
		if ( pose.data().has( CacheableDataType::JOBDIST_OUTPUT_TAG ) ) {
			tag =
				static_cast< basic::datacache::CacheableString const & >
				( pose.data().get( CacheableDataType::JOBDIST_OUTPUT_TAG ) ).str();
		}
		//std::cout << "SVM START " << tag << std::endl;
		for ( Size i=1; i<=nres_; ++i ) {
			std::cout << "SVM " << tag << " ";
			for ( Size j=1; j<=260; ++j ) {
				std::cout << j << ":" << feature_vector(i,j) << " ";
			}
			std::cout << std::endl;
		}
	}
}


void ProQ_Energy::atom_feature(pose::Pose & pose,ObjexxFCL::FArray2D< Real > & vec,Size index,int windowsize) const{
	//const int nres=(int)pose.size();
	const int nres=(int)nres_;

	//const int windowsize=21;
	ObjexxFCL::FArray3D< Size >atomcontacts(nres,13,13);
	atomcontacts=0;
	// ObjexxFCL::FArray2D< Real >res_contact_dist(nres,nres);
	for ( int i = 1; i <= nres; i=i+1 ) {
		conformation::Residue const & rsd_i( pose.residue(i) );
		const Energies & energies( pose.energies() );
		const TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph() );
		for ( utility::graph::Graph::EdgeListConstIter
				ir  = graph.get_node(i)->const_edge_list_begin(),
				ire = graph.get_node(i)->const_edge_list_end();
				ir != ire; ++ir ) {
			int const j( (*ir)->get_other_ind( i ) );
			//    for(int j=i+1;j<=nres;++j) {
			if ( std::abs(i-j)<=1 ) continue;

			conformation::Residue const & rsd_j( pose.residue(j) );
			for ( Size ii=1; ii<=rsd_i.nheavyatoms(); ++ii ) {
				core::conformation::Atom const & atom_i = rsd_i.atom(ii);
				int atype1=atom13(rsd_i,ii);
				//std::cout << "ATOMTYPES " << rsd_i.name3() << rsd_i.atom_name(ii) << " " << atype1 << std::endl;
				if ( atype1==0 ) {
					// std::cout << "Skipping: " << rsd_i.atom_name(ii) << std::endl;
					continue;

				}

				for ( Size jj=1; jj<=rsd_j.nheavyatoms(); ++jj ) {
					core::conformation::Atom const & atom_j = rsd_j.atom(jj);
					int atype2=atom13(rsd_j,jj);
					//std::cout << "ATOMTYPES " << rsd_j.name3() << rsd_j.atom_name(jj) << " " << atype2 << std::endl;
					if ( atype2==0 ) {
						continue;
					}
					//std::cout << i << " " << j << " jj: " << jj << " ii:" << ii << std::endl;
					Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());
					if ( sqdist >= 16 ) continue;

					/*if(i<=11 || j<=11) {
					if(atype1==1 && atype2==3) {
					printf("Contact between %d %d - %d %d dist: %f\n",i,atype1,j,atype2,sqrt(sqdist));
					}
					}*/
					// std::cout <<  rsd_i.name3() << " " << rsd_i.atom_name(ii) << " " << atype1 << " - " << rsd_j.name3() << " " << rsd_j.atom_name(jj) << " " << atype2 << " " << sqrt(sqdist) << std::endl;
					atomcontacts(i,atype1,atype2)++;
					//atomcontacts(j,atype2,atype1)++;
				}
			}
		}
	}
	int offset_bug=0;
	if ( basic::options::option[ basic::options::OptionKeys::ProQ::prof_bug ]() ) {
		offset_bug=1;
	}
	ObjexxFCL::FArray2D< Real >atom(13,13);
	ObjexxFCL::FArray1D< Real >atomfrac(91);
	for ( int i=1; i<=nres; ++i ) {  //central residue
		atom=0;
		atomfrac=0;
		Real atomcount=0;
		for ( int j=i-(windowsize-1)/2; j<=i+(windowsize-1)/2; ++j ) { // each residue in window /retrain with a tertiary window <12Å
			if ( j>=1+offset_bug && j<=nres ) {
				for ( int k=1; k<=13; ++k ) {
					for ( int l=1; l<=13; l++ ) {
						atom(k,l)+=atomcontacts(j-offset_bug,k,l);
						atomcount+=atomcontacts(j-offset_bug,k,l);
					}
				}

			}
		}
		//normalize
		//  std::cout << "ATOMCOUNT: " << atomcount << std::endl;
		if ( atomcount!=0 ) {
			int n=1;
			for ( int k=1; k<=13; ++k ) {
				atomfrac(n)=atom(k,k)/atomcount;
				//   if(n==56) {
				//    std::cout << "ATOM56: " << k << "," <<  k << std::endl;
				//   }
				++n;
				for ( int l=k+1; l<=13; l++ ) {
					atomfrac(n)=(atom(k,l)+atom(l,k))/atomcount;
					//   if(n==56) {
					//    std::cout << "ATOM56: " << k << "," <<  l << std::endl;
					//   }
					++n;
				}
			}

		}
		// std::cout << "POS " << atomcount << " " << i << ": ";
		for ( int j=1; j<=91; ++j ) {
			//int pos=index+j-1;
			//   std::cout << index+j-1 << ":" << atomfrac(j) << " " ;
			vec(i,index+j-1)=atomfrac(j);
		}
		//  std::cout << std::endl;
	}
}

void ProQ_Energy::res_feature(pose::Pose & pose, ObjexxFCL::FArray2D< Real > & vec,Size index, int windowsize) const{
	const int nres=(int)nres_;
	ObjexxFCL::FArray2D< bool >res_contact_map(nres,nres);
	ObjexxFCL::FArray2D< Real >res_contact_dist(nres,nres);
	res_contact_map=false;

	for ( int i = 1; i <= nres; i=i+1 ) {
		//i=5;
		// get the appropriate residue from the pose.
		//conformation::Residue const & rsd_i( pose.residue(i) );
		//  std::cout << "RES: " << res6(rsd_i) << " " << rsd_i << std::endl;
		//Size const atomindex_i = rsd.atom_index("CA"); //representative_atom_name( rsd.aa() ));
		//core::conformation::Atom const & atom_i = rsd.atom(atomindex_i);
		const Energies & energies( pose.energies() );
		const TwelveANeighborGraph & graph ( energies.twelveA_neighbor_graph() );
		// iterate across neighbors within 12 angstroms
		//  for ( int j = i+1; j <= nres; j=j+1 ) {
		//j=53;
		for ( utility::graph::Graph::EdgeListConstIter
				ir  = graph.get_node(i)->const_edge_list_begin(),
				ire = graph.get_node(i)->const_edge_list_end();
				ir != ire; ++ir ) {
			int const j( (*ir)->get_other_ind( i ) );
			if ( std::abs(i-j)<=5 ) continue;

			Real dist=crd(pose,i,j);
			if ( dist<36 ) { // 6Å
				//   conformation::Residue const & rsd_j( pose.residue(j) );
				//Size atomindex_j( rsd_j.type().nbr_atom() );
				//core::conformation::Atom const & atom_j = rsd_j.atom(atomindex_j);
				//Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());
				res_contact_map(i,j)=true;
				res_contact_map(j,i)=true;
				res_contact_dist(i,j)=sqrt(dist);
				res_contact_dist(j,i)=sqrt(dist);
			}
		}
	}

	ObjexxFCL::FArray2D< Real >res_contact_prof(6,6);
	ObjexxFCL::FArray1D< Real >resfrac(21);
	for ( int i=1; i<=nres; ++i ) { //central residue
		// std::cout << i << " " << pose.residue(i).name3() << " " << res6(pose.residue(i).name1()) << " " << res6(pose.residue(i)) << ": " << std::endl;;
		res_contact_prof=0;
		resfrac=0;
		Size total_res_contacts_count(0);
		Real total_res_contacts(0);
		// std::cout << "Done init... " << std::endl;

		for ( int j=i-(windowsize-1)/2; j<=i+(windowsize-1)/2; ++j ) { // each residue in window /retrain with a tertiary window <12Å
			if ( j>=1 && j<=nres ) {
				//    std::cout << j << " ";
				for ( int k=1; k<=nres; ++k ) {
					if ( res_contact_map(j,k) ) {
						//std::cout << " Contact between " << j << " and " << k << " dist: " << res_contact_dist(j,k) << ",restype1: (" << pose.residue(j).name1()<< ")" << res_j << ", restype2: " << " (" << pose.residue(k).name1() << ") " << res_k << std::endl;
						total_res_contacts_count++;
						for ( int l=1; l<=20; l++ ) {
							int res1=res6(profile_index_to_aa(l));

							for ( int n=1; n<=20; n++ ) {
								int res2=res6(profile_index_to_aa(n));
								if ( basic::options::option[ basic::options::OptionKeys::ProQ::prof_bug ]() ) {
									res2=res6(pose.residue(k));
								}
								Real joined_prob=prob_profile_(j,l)*prob_profile_(k,n);
								//  std::cout << "PROF: " << i << " " << j << " " << aa1 << " " << aa2 << " " << l << " " << n << " " << res1 << " " << res2 << " " <<  prob_profile_(j,l) << " " << prob_profile_(k,n) << std::endl;
								res_contact_prof(res1,res2)+=joined_prob;
								total_res_contacts+=joined_prob;
								//         std::cout << aa1 << " " << res1 << " " << aa2 << " " << res2 << std::endl;
							}
						}
					}
				}
			}
		}
		if ( total_res_contacts!=0 ) {
			int n=1;
			for ( int j = 1; j <= 6; ++j ) {
				resfrac(n)=res_contact_prof(j,j)/total_res_contacts;
				//if(n==17 || n==14 || n==19) {
				// std::cout << "RES: " << n << " " <<j << "," << j << std::endl;

				//}
				// std::cout << i << " " << j << " " << j << " " << res_contact_prof(j,j) << " " << total_res_contacts <<  " " << resfrac(n) << std::endl;
				++n;
				for ( int k =j+1; k <= 6; ++k ) {
					resfrac(n)=(res_contact_prof(k,j)+res_contact_prof(j,k))/total_res_contacts;
					//if(n==17 || n==14 || n==19) {
					// std::cout << "RES: " << n << " " << j << "," << k << std::endl;
					//}
					//std::cout << i << " " << j << " " << k << " " << res_contact_prof(k,j)+res_contact_prof(j,k) << " " << total_res_contacts <<  " " << resfrac(n) << std::endl;
					++n;
				}
			}
		}
		for ( int j = 1; j <= 21; ++j ) {
			//std::cout << pos << ":" << resfrac(j) << " " ;
			vec(i,index+j-1)=resfrac(j);
		}
		//std::cout << std::endl;
	}

	// for(int i=-(windowsize-1)/2;i<nres-(windowsize-1)/2;++i) {
	//  for ( int j = i; j < i+windowsize; ++j ) {
	//   std::cout << j << " ";
	//   //   if(res_contact_map(i,j))
	//  }
	//  std::cout << std::endl;
	// }
	//std::exit(1);
	//std::cout << "TMP " << tmp << std::endl;


}
int ProQ_Energy::atom13(conformation::Residue const & rsd,Size const atom_i) const {
	//int type=atom13_0(rsd,atom_i)+1;
	//std::cout << "Type: " << type << std::endl;
	//return type;
	return atom13_0(rsd,atom_i)+1;
}

int ProQ_Energy::atom13_0(conformation::Residue const & rsd,Size const atom_i) const {

	/* Types  - Return number/
	*
	* C         - 0
	* N         - 1
	* O         - 2
	* CA        - 3
	* CH3       - 4
	* CH/CH2    - 5
	* C(OO-)    - 6
	* NH        - 7
	* NH2       - 8
	* (C)OO-    - 9
	* =O        - 10
	* OH        - 11
	* S         - 12
	* OXT       - 13
	*/
	//printf("%s %s\n",name,res);
	//if(rsd.name().compare("VRT")==0) {
	// std::cout << rsd << atom_i << " " << rsd.atom_type(atom_i).is_virtual() << std::endl;
	//}
	if ( rsd.aa()==core::chemical::aa_vrt || rsd.atom_type(atom_i).is_virtual() ) {
		return -1;
	}

	std::string const & name=rsd.atom_name(atom_i);
	std::string const & res=rsd.name3();

	if ( name==" C  " ) {
		return 0;
	}
	if ( name==" N  " ) {
		return 1;
	}
	if ( name==" O  " ) {
		return 2;
	}
	if ( name==" CA " ) {
		return 3;
	}
	if (
			(res=="ILE" && (name==" CD1" || name==" CG2")) ||
			(res=="LEU" && (name==" CD1" || name==" CD2")) ||
			(res=="MET" && name==" CE ") ||
			(res=="THR" && name==" CG2") ||
			(res=="ALA" && name==" CB ") ||
			(res=="VAL" && (name==" CG1" || name==" CG2")) ) {
		return 4;
	}
	if ( (res=="ASP" && name==" CG ") ||
			(res=="GLU" && name==" CG ") ) {
		return 6;
	}
	if ( name==" NE " || name==" ND1" || name==" NE1" ) {
		return 7;
	}
	if ( name==" NH1" || name==" NH2" || name==" ND2"
			|| name==" NE2" || name==" NZ " ) {
		return 8;
	}
	if ( (res=="ASP" && (name==" OD2" || name==" OD1")) ||
			(res=="GLU" && (name==" OE1" || name==" OE2")) ) {
		return 9;
	}
	if ( (res=="ASN" && name==" OD1") ||
			(res=="GLN" && name==" OE1") ) {
		return 10;
	}
	if ( name==" OG " || name==" OH " || name==" OG1" ) {
		return 11;
	}
	if ( name==" SG " || name==" SD " ) {
		return 12;
	}
	// if(name=="OXT")
	// return 13;
	// std::cout << "ATOM13: '" << res << "' '" << name << "'" << std::endl;
	if ( name.compare(1,1,"C")==0 ) {
		//std::cout << "ATOM13 res: " << res << " atom: " <<  name << std::endl;
		return 5;
	} else {
		//std::cout << "ATOM13_NOT_SET: '" << res << "' '" << name << "'" << std::endl;
		//std::cout << rsd.atom_type(atom_i) ;;
	}
	return -1;
}

int ProQ_Energy::res6(conformation::Residue const & rsd) const {
	const char aa(rsd.name1());
	return res6(aa);
}

int ProQ_Energy::res6(char const aa) const {

	switch(aa) {
	case 'R' : return 1;
	case 'K' : return 1;
	case 'D' : return 2;
	case 'E' : return 2;
	case 'H' : return 3;
	case 'F' : return 3;
	case 'W' : return 3;
	case 'Y' : return 3;
	case 'N' : return 4;
	case 'Q' : return 4;
	case 'S' : return 4;
	case 'T' : return 4;
	case 'A' : return 5;
	case 'I' : return 5;
	case 'L' : return 5;
	case 'M' : return 5;
	case 'V' : return 5;
	case 'C' : return 5;
	case 'G' : return 6;
	case 'P' : return 6;

	}
	return 0;

}
char ProQ_Energy::profile_index_to_aa(int i) const {

	static char aa[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};
	debug_assert(i>=1 && i<=20);
	return aa[i-1];
}

//Closest residue distance as defined in ProQres.c
Real ProQ_Energy::crd(pose::Pose & pose,Size i,Size j) const {
	conformation::Residue const & rsd_i( pose.residue(i) );
	conformation::Residue const & rsd_j( pose.residue(j) );
	Real mindist(99999999);
	//ignore backbone except CA
	for ( Size ii=1; ii<=rsd_i.nheavyatoms(); ++ii ) {
		//std::cout << "ii: " << ii << " nheavy:" << rsd_i.nheavyatoms() << " " << "name: " << rsd_i.atom_type(ii).name() << " " << rsd_i.atom_is_backbone(ii) << " " << rsd_i.atom_type(ii).name().rfind("bb") << " " << rsd_i.atom_type(ii).name().find("CA") << std::endl;
		if ( rsd_i.atom_is_backbone(ii) && rsd_i.atom_name(ii).compare(1,2,"CA") != 0 ) continue;

		core::conformation::Atom const & atom_i = rsd_i.atom(ii);

		//std::cout << "ii: " << ii << " " << atom_i << ", name: " << rsd_i.atom_type(ii).name() << " " << rsd_i.atom_type(ii).atom_type_name() << std::endl;
		//std::cout << rsd_i.atom_type(ii);

		for ( Size jj=1; jj<=rsd_j.nheavyatoms(); ++jj ) {
			if ( rsd_j.atom_is_backbone(jj) && rsd_j.atom_name(jj).compare(1,2,"CA")!=0 ) continue;

			core::conformation::Atom const & atom_j = rsd_j.atom(jj);
			Real sqdist = atom_i.xyz().distance_squared(atom_j.xyz());

			//std::cout << "jj: " << jj << " " << atom_j << ", name: " << rsd_i.atom_type(jj).name() << " " << rsd_i.atom_type(jj).atom_type_name() << std::endl;
			//std::cout << i << ":" << ii << " " << rsd_i.atom_name(ii) << " " << j << ":" << jj << " " << rsd_j.atom_name(jj) << " " << sqrt(sqdist) << std::endl;
			if ( sqdist<mindist ) {
				mindist=sqdist;
			}
		}
	}

	return mindist;
}

void ProQ_Energy::sum_profile(Size j,ObjexxFCL::FArray1D< Real > & vec) const {
	for ( int k=1; k<=20; k++ ) {
		int res=res6(profile_index_to_aa(k));
		vec(res)+=prob_profile_(j,k);
	}
}

void ProQ_Energy::surf_feature(pose::Pose &, utility::vector1< Real > & rsd_sasa_rel, ObjexxFCL::FArray2D< Real > & vec,
	Size index, int windowsize) const
{
	const int nres=(int)nres_;

	ObjexxFCL::FArray1D< Real > surface25_prof(6);
	ObjexxFCL::FArray1D< Real > surface50_prof(6);
	ObjexxFCL::FArray1D< Real > surface75_prof(6);
	ObjexxFCL::FArray1D< Real > surface100_prof(6);
	int start=index;
	//Real surface_total_prof=0;
	for ( int i=1; i<=nres; ++i ) {
		int surface_total=0;
		index=start;
		surface25_prof=0;
		surface50_prof=0;
		surface75_prof=0;
		surface100_prof=0;
		for ( int j=i-(windowsize-1)/2; j<=i+(windowsize-1)/2; ++j ) {
			if ( j>=1 && j<=nres ) {
				surface_total++;
				if ( rsd_sasa_rel[j]<=25 ) {
					sum_profile(j,surface25_prof);
				} else if ( rsd_sasa_rel[j]>25 && rsd_sasa_rel[j]<=50 ) {
					sum_profile(j,surface50_prof);
				} else if ( rsd_sasa_rel[j]>50 && rsd_sasa_rel[j]<=75 ) {
					sum_profile(j,surface75_prof);
				} else if ( rsd_sasa_rel[j]>75 ) {
					sum_profile(j,surface100_prof);
				}
			}
		}

		//normalize
		for ( int j=1; j<=6; ++j ) {
			surface25_prof(j)/=surface_total;
			surface50_prof(j)/=surface_total;
			surface75_prof(j)/=surface_total;
			surface100_prof(j)/=surface_total;
		}
		for ( int j=1; j<=6; ++j,++index ) {
			vec(i,index)=surface25_prof(j);
		}
		for ( int j=1; j<=6; ++j,++index ) {
			vec(i,index)=surface50_prof(j);
		}
		for ( int j=1; j<=6; ++j,++index ) {
			vec(i,index)=surface75_prof(j);
		}
		for ( int j=1; j<=6; ++j,++index ) {
			vec(i,index)=surface100_prof(j);
		}
		/*std::cout << i << " SURF: ";
		for(int j=1,index=start;j<=24;++j,++index) {
		std::cout << index << ":" << vec(i,index) << " ";;

		}
		std::cout << std::endl;
		*/

	}
	//std::exit(1);
}

void ProQ_Energy::stride_feature(pose::Pose &, dssp::Dssp & ss, ObjexxFCL::FArray2D< Real > & vec, int pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;

	Real h(0);
	Real e(0);
	Real c(0);
	//std::cout << pos << ": ";
	for ( int j=pos-(windowsize-1)/2; j<=pos+(windowsize-1)/2; ++j ) {
		h=0;
		e=0;
		c=0;
		if ( j>=1 && j<=nres ) {
			if ( ss.get_dssp_secstruct(j) == 'H' ) {
				h=1;
			} else if ( ss.get_dssp_secstruct(j) == 'E' ) {
				e=1;
			} else if ( ss.get_dssp_secstruct(j) == 'L' ) {
				c=1;
			}

		}
		vec(pos,index++)=h;
		vec(pos,index++)=e;
		vec(pos,index++)=c;
		/*std::cout << "(" << j << ":" << ss.get_dssp_secstruct(j) << ":" << h << "," << e << "," << c << "): ";
		std::cout << index-3 << ":" << vec(pos,index-3) << " "
		<< index-2 << ":" << vec(pos,index-2) << " "
		<< index-1 << ":" << vec(pos,index-1) << " "; */
	}
	//std::cout << std::endl;
}

void ProQ_Energy::ss_feature(pose::Pose &, dssp::Dssp & ss, ObjexxFCL::FArray2D< Real > & vec, int pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;
	int ss_type(0);
	Real ss_prob(0);
	//std::cout << pos << ": ";
	for ( int j=pos-(windowsize-1)/2; j<=pos+(windowsize-1)/2; ++j ) {
		//1.coil,2.helix,3.sheet
		if ( j>=1 && j<=nres ) {
			if ( ss.get_dssp_secstruct(j) == 'H' ) {
				ss_type=2;
			} else if ( ss.get_dssp_secstruct(j) == 'E' ) {
				ss_type=3;
			} else if ( ss.get_dssp_secstruct(j) == 'L' ) {
				ss_type=1;
			}
			ss_prob=ss_pred_(j,ss_type);
		}
		//std::cout << "(" << j << ":" << ss.get_dssp_secstruct(j) << ":" << ss1_(j)<< ":" << ss_pred_(j,1) << "," << ss_pred_(j,2) << "," << ss_pred_(j,3) << "): ";
		vec(pos,index++)=ss_prob;
		//std::cout << index-1 << ":" << vec(pos,index-1) << " ";
	}
	//std::cout << std::endl;
}

void ProQ_Energy::gss_sc_feature(pose::Pose &, dssp::Dssp & ss, ObjexxFCL::FArray2D< Real > & vec, Size index) const
{
	const int nres=(int)nres_;

	int match(0);
	for ( int i=1; i<=nres; ++i ) {
		if ( ss.get_dssp_secstruct(i) == ss1_(i) ) {
			++match;
		}
	}
	Real Q3=(Real)match/nres;
	for ( int i=1; i<=nres; ++i ) {
		vec(i,index)=Q3;
	}
}

void ProQ_Energy::grsa_sc_feature(pose::Pose &, utility::vector1< Real> & rsd_sasa_rel,
	ObjexxFCL::FArray2D< Real > & vec, Size index) const
{
	const int nres=(int)nres_;

	int match(0);
	for ( int i=1; i<=nres; ++i ) {
		//std::cout << "GRSA: " << rsd_sasa_rel[i] << " " << rsa_class_pred_(i) ;
		if ( (rsd_sasa_rel[i] < 25 && rsa_class_pred_(i) == 'b') ||
				(rsd_sasa_rel[i] > 25 && rsa_class_pred_(i) == 'e') ) {
			match++;
			//std::cout << " MATCH ";
		}
		//std::cout << std::endl;
	}

	Real Q3=(Real)match/nres;
	for ( int i=1; i<=nres; ++i ) {
		vec(i,index)=Q3;
	}
}

void ProQ_Energy::ss_sc_feature(pose::Pose &, dssp::Dssp & ss, ObjexxFCL::FArray2D< Real > & vec, int const pos,
	Size index, int windowsize) const
{
	const int nres=(int)nres_;

	int count(0);
	int correct(0);
	for ( int j=pos-(windowsize-1)/2; j<=pos+(windowsize-1)/2; ++j ) {
		if ( j>=1 && j<=nres ) {
			if ( ss.get_dssp_secstruct(j) == ss1_(j) ) {
				++correct;
			}
			++count;
		}
	}
	vec(pos,index)=(Real)correct/count;
	//std::cout << "Q3: " << pos << ":" << vec(pos,index) << std::endl;;
}

void ProQ_Energy::topology_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index) const
{
	const int win=11;
	const int nres=(int)nres_;

	// std::cout << "TOPOLOGY " << pos << " " << topology_.tmregion(pos) << " : " ;
	for ( int i=(int)pos-(win-1)/2; i<=(int)pos+(win-1)/2; ++i,++index ) {
		//  std::cout << pos ;
		vec(pos,index)=0;
		if ( i>=1 && i <= nres ) {
			vec(pos,index)=topology_.tmregion(i);
		}
		//std::cout << " " << i << ":" << index << "(" << vec(pos,index) << ")";
	}

	//std::cout << std::endl;
}

void ProQ_Energy::zpred_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index) const
{
	const int nres=(int)nres_;
	const int win=1;
	for ( int i=(int)pos-(win-1)/2; i<=(int)pos+(win-1)/2; ++i,++index ) {
		Real z_pred(0);
		if ( i<1 ) {
			z_pred=z_pred_(1);
		} else if ( i>nres ) {
			z_pred=z_pred_(nres);
		} else {
			z_pred=z_pred_(i);
		}
		vec(pos,index)=(z_pred-5)/20;  //min zpred is 5 => 0 in membrane > 1 outside 20Å.
		//std::cout << "ZPRED " << pos << ": " << vec(pos,index) << std::endl;
	}
}

void ProQ_Energy::calculateZ(pose::Pose & pose,ObjexxFCL::FArray1D< Real > & Z) const {

	if ( !basic::options::option[ basic::options::OptionKeys::ProQ::membrane ]() ) {
		return;
	}

	//Determine the membrane center and normal
	Vector normal(0);
	Vector center(0);
	Vector inside(0);
	Vector outside(0);
	for ( Size i=1; i<=topology_.tmhelix(); ++i ) {
		Vector const & start( pose.residue( topology_.span_begin(i) ).atom( 2 ).xyz());
		Vector const & end( pose.residue( topology_.span_end(i) ).atom( 2 ).xyz());
		// all odd helices goes from outside in (from c++)
		if ( topology_.helix_id(i) % 2 == 0 ) {
			inside+=start;
			outside+=end;
		} else {
			outside+=start;
			inside+=end;
		}
	}
	normal=outside-inside;
	normal.normalize();
	center=(Real)0.5*(outside+inside)/topology_.tmhelix();
	//std::cout << "NORMAL: " << normal.x() << ","  << normal.y() << "," << normal.z() << std::endl;
	//std::cout << "CENTER: " << center.x() << ","  << center.y() << "," << center.z() << std::endl;

	//project CA on the membrane normal and subtract center
	Z.dimension(nres_);
	for ( Size i=1; i<=nres_; ++i ) {
		Vector const & coord( pose.residue( i ).atom( 2 ).xyz());
		Z(i)=fabs(normal.dot(coord-center));
		//std::cout << i << " " << Z(i) << std::endl;
	}
}

void ProQ_Energy::z_feature(pose::Pose &, ObjexxFCL::FArray1D< Real > & Z, ObjexxFCL::FArray2D< Real > & vec,
	int const pos, Size index) const
{
	const int nres=(int)nres_;

	const int win=1;
	for ( int i=(int)pos-(win-1)/2; i<=(int)pos+(win-1)/2; ++i,++index ) {
		Real z(0);
		if ( i<1 ) {
			z=Z(1);
		} else if ( i>nres ) {
			z=Z(nres);
		} else {
			z=Z(i);
		}

		if ( z<5 ) {
			z=5;
		} else if ( z>25 ) {
			z=25;
		}

		vec(pos,index)=(z-5)/20;  //min zpred is 5 => 0 in membrane > 1 outside 20Å.
		//std::cout << "Z " << pos << ": " << vec(pos,index) << std::endl;
	}
}

void ProQ_Energy::entropy_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;

	for ( int i=(int)pos-(windowsize-1)/2; i<=(int)pos+(windowsize-1)/2; ++i,++index ) {
		vec(pos,index)=0;
		if ( i>=1 && i<=nres ) {
			vec(pos,index)=entropy_(i);
		}
		//std::cout << "Entropy " << pos << "," << index << " " << i << " " << " : " << vec(pos,index) << std::endl;
	}
}

void ProQ_Energy::profile_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;

	for ( int i=(int)pos-(windowsize-1)/2; i<=(int)pos+(windowsize-1)/2; ++i ) {
		if ( i>=1 && i<=nres ) {
			for ( int j=1; j<=20; ++j,++index ) {
				vec(pos,index)=scaled_logodds_profile_(i,j);
				//std::cout << "profile " << pos << "," << index << " " << i << " : " << vec(pos,index) << std::endl;
			}
		} else {
			for ( int j=1; j<=20; ++j,++index ) {
				vec(pos,index)=0;
				//std::cout << "profile " << pos << "," << index << " " << i << " : " << vec(pos,index) << std::endl;
			}
		}
	}
}

void ProQ_Energy::rsa_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, utility::vector1< Real> & rsd_sasa_rel,
	int const pos, Size index, int windowsize) const
{
	const int nres=(int)nres_;
	for ( int i=(int)pos-(windowsize-1)/2; i<=(int)pos+(windowsize-1)/2; ++i,++index ) {
		vec(pos,index)=0;
		if ( i>=1 && i<=nres ) {
			vec(pos,index)=rsd_sasa_rel[i]/100;
		}
		//std::cout << "RSA_FEAT " << pos << "," << index << " " << i << " : " << vec(pos,index) << std::endl;
	}
}

void ProQ_Energy::prsa_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;

	for ( int i=(int)pos-(windowsize-1)/2; i<=(int)pos+(windowsize-1)/2; ++i,++index ) {
		vec(pos,index)=0;
		if ( i>=1 && i<=nres ) {
			vec(pos,index)=rsa_pred_(i)/100;
		}
		//std::cout << "PRSA_FEAT " << pos << "," << index << " " << i << " : " << vec(pos,index) << std::endl;
	}
}

void ProQ_Energy::rsa_sc_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec,
	utility::vector1< Real> & rsd_sasa_rel, int const pos, Size index, int windowsize) const
{
	const int nres=(int)nres_;

	int count(0);
	int agree(0);
	vec(pos,index)=0;
	for ( int i=(int)pos-(windowsize-1)/2; i<=(int)pos+(windowsize-1)/2; ++i ) {
		if ( i>=1 && i<=nres ) {
			if ( (rsd_sasa_rel[i] < 25 && rsa_class_pred_(i) == 'b') ||
					(rsd_sasa_rel[i] > 25 && rsa_class_pred_(i) == 'e') ) {
				++agree;
			}
			++count;
		}
	}

	vec(pos,index)=(Real)agree/count;
	//std::cout << "RSA_SC " << pos << "," << index << " " << " : " << vec(pos,index) << " " << agree << " " << count << std::endl;
}

void ProQ_Energy::termini_feature(pose::Pose &, ObjexxFCL::FArray2D< Real > & vec, int const pos, Size index,
	int windowsize) const
{
	const int nres=(int)nres_;

	Real distN=(Real)(pos-1)/windowsize;
	Real distC=(Real)(nres-pos)/windowsize;
	if ( distN>1 ) {
		distN=1;
	}
	if ( distC>1 ) {
		distC=1;
	}
	vec(pos,index)=distN;
	vec(pos,index+1)=distC;
	//std::cout << "TERMINI " << pos << ": " << vec(pos,index) << " " << vec(pos,index+1) << std::endl;
}

} // methods
} // scoring
} // core

