// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

#include <core/types.hh>
#include <devel/init.hh>

#include <basic/options/option.hh>

#include <core/pose/Pose.hh>
#include <core/io/pdb/pose_io.hh>

#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>

#include <core/scoring/packstat/types.hh>
#include <core/scoring/packstat/PackingScore.hh>
#include <core/scoring/packstat/compute_sasa.hh>
#include <core/scoring/packstat/LeeRichards.hh>
#include <core/scoring/packstat/AtomRadiusMap.hh>
#include <core/scoring/packstat/SimplePDB.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <core/kinematics/MoveMap.hh>

#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/xyz.io.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/random/random.hh>

#include <core/scoring/rms_util.hh>

#include <time.h>

#include <pthread.h>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packstat.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>





using core::Real;
using core::Size;
using utility::vector1;
using numeric::xyzVector;

vector1<Real> uniforms;

void init_myuniform() {
	uniforms.clear();
	uniforms.push_back( 0.0331743322312832 ); uniforms.push_back( 0.212354763410985 );
	uniforms.push_back( 0.465039731701836 ); uniforms.push_back( 0.820828061783686 );
	uniforms.push_back( 0.568490081233904 ); uniforms.push_back( 0.761286014690995 );
	uniforms.push_back( 0.912834479706362 ); uniforms.push_back( 0.491649915929884 );
	uniforms.push_back( 0.279530025785789 ); uniforms.push_back( 0.326601809356362 );
	uniforms.push_back( 0.0305851215962321 ); uniforms.push_back( 0.892525085946545 );
	uniforms.push_back( 0.347199349664152 ); uniforms.push_back( 0.0469867852516472 );
	uniforms.push_back( 0.466409785673022 ); uniforms.push_back( 0.932394379982725 );
	uniforms.push_back( 0.661947170970961 ); uniforms.push_back( 0.922915467992425 );
	uniforms.push_back( 0.603194614639506 ); uniforms.push_back( 0.52322421502322 );
	uniforms.push_back( 0.87211511656642 ); uniforms.push_back( 0.749204269144684 );
	uniforms.push_back( 0.0501798943150789 ); uniforms.push_back( 0.647677037632093 );
	uniforms.push_back( 0.212808817625046 ); uniforms.push_back( 0.43113399296999 );
	uniforms.push_back( 0.354633759474382 ); uniforms.push_back( 0.258986184839159 );
	uniforms.push_back( 0.982468512142077 ); uniforms.push_back( 0.159200440626591 );
	uniforms.push_back( 0.650892873760313 ); uniforms.push_back( 0.152489541098475 );
	uniforms.push_back( 0.653112890431657 ); uniforms.push_back( 0.570019864244387 );
	uniforms.push_back( 0.0838310152757913 ); uniforms.push_back( 0.910201397491619 );
	uniforms.push_back( 0.218164353631437 ); uniforms.push_back( 0.0525630488991737 );
	uniforms.push_back( 0.0118643417954445 ); uniforms.push_back( 0.413506688550115 );
	uniforms.push_back( 0.932123879203573 ); uniforms.push_back( 0.669445681385696 );
	uniforms.push_back( 0.169370586052537 ); uniforms.push_back( 0.296843440970406 );
	uniforms.push_back( 0.71213926654309 ); uniforms.push_back( 0.93877394683659 );
	uniforms.push_back( 0.225173749262467 ); uniforms.push_back( 0.289904457982630 );
	uniforms.push_back( 0.372512749163434 ); uniforms.push_back( 0.846214164746925 );
	uniforms.push_back( 0.978809693828225 ); uniforms.push_back( 0.956038564676419 );
	uniforms.push_back( 0.312367196427658 ); uniforms.push_back( 0.806694239610806 );
	uniforms.push_back( 0.299944201717153 ); uniforms.push_back( 0.119206607108936 );
	uniforms.push_back( 0.484675208339468 ); uniforms.push_back( 0.694195500342175 );
	uniforms.push_back( 0.258086859248579 ); uniforms.push_back( 0.974756510462612 );
	uniforms.push_back( 0.174823179841042 ); uniforms.push_back( 0.0478785356972367 );
	uniforms.push_back( 0.842924054479226 ); uniforms.push_back( 0.845733262598515 );
	uniforms.push_back( 0.683667376404628 ); uniforms.push_back( 0.905862155603245 );
	uniforms.push_back( 0.242450795834884 ); uniforms.push_back( 0.76471592951566 );
	uniforms.push_back( 0.544274961343035 ); uniforms.push_back( 0.3281284798868 );
	uniforms.push_back( 0.304273429326713 ); uniforms.push_back( 0.483715623384342 );
	uniforms.push_back( 0.156516925431788 ); uniforms.push_back( 0.509351484943181 );
	uniforms.push_back( 0.146319671766832 ); uniforms.push_back( 0.0990302439313382 );
	uniforms.push_back( 0.0902824332006276 ); uniforms.push_back( 0.854539727093652 );
	uniforms.push_back( 0.441274259705096 ); uniforms.push_back( 0.71080990252085 );
	uniforms.push_back( 0.873236656654626 ); uniforms.push_back( 0.761584523133934 );
	uniforms.push_back( 0.903422515373677 ); uniforms.push_back( 0.0540094620082527 );
	uniforms.push_back( 0.330672544427216 ); uniforms.push_back( 0.478437855839729 );
	uniforms.push_back( 0.46681995340623 ); uniforms.push_back( 0.34077229979448 );
	uniforms.push_back( 0.968053152784705 ); uniforms.push_back( 0.512479383731261 );
	uniforms.push_back( 0.522685769014060 ); uniforms.push_back( 0.860231585334986 );
	uniforms.push_back( 0.131522696232423 ); uniforms.push_back( 0.861790018621832 );
	uniforms.push_back( 0.518712243065238 ); uniforms.push_back( 0.994612438371405 );
	uniforms.push_back( 0.46870283363387 ); uniforms.push_back( 0.965584655292332 );
	uniforms.push_back( 0.571425111265853 ); uniforms.push_back( 0.450497290119529 );
	uniforms.push_back( 0.596647670958191 ); uniforms.push_back( 0.435171453049406 );
	uniforms.push_back( 0.462421316420659 ); uniforms.push_back( 0.0265133723150939 );
	uniforms.push_back( 0.06018988462165 ); uniforms.push_back( 0.136005905224010 );
	uniforms.push_back( 0.456596181029454 ); uniforms.push_back( 0.636125463526696 );
	uniforms.push_back( 0.347120713442564 ); uniforms.push_back( 0.988651815569028 );
	uniforms.push_back( 0.317874954780564 ); uniforms.push_back( 0.767821838380769 );
	uniforms.push_back( 0.106660933699459 ); uniforms.push_back( 0.900170157896355 );
	uniforms.push_back( 0.674955864902586 ); uniforms.push_back( 0.7744529126212 );
	uniforms.push_back( 0.0682156919501722 ); uniforms.push_back( 0.16353478631936 );
	uniforms.push_back( 0.47398390294984 ); uniforms.push_back( 0.0750290653668344 );
	uniforms.push_back( 0.546150231733918 ); uniforms.push_back( 0.103311127051711 );
	uniforms.push_back( 0.651617111172527 ); uniforms.push_back( 0.295266589149833 );
	uniforms.push_back( 0.36726274411194 ); uniforms.push_back( 0.302034016698599 );
	uniforms.push_back( 0.806727706454694 ); uniforms.push_back( 0.8731925918255 );
	uniforms.push_back( 0.658777604810894 ); uniforms.push_back( 0.858751825988293 );
	uniforms.push_back( 0.987147206207737 ); uniforms.push_back( 0.765560046071187 );
	uniforms.push_back( 0.107388277072459 ); uniforms.push_back( 0.476451081689447 );
	uniforms.push_back( 0.453372968826443 ); uniforms.push_back( 0.477047462016344 );
	uniforms.push_back( 0.511695447377861 ); uniforms.push_back( 0.741378476610407 );
	uniforms.push_back( 0.268000422744080 ); uniforms.push_back( 0.67332601454109 );
	uniforms.push_back( 0.550826609833166 ); uniforms.push_back( 0.730828423053026 );
	uniforms.push_back( 0.158377447165549 ); uniforms.push_back( 0.969179981388152 );
	uniforms.push_back( 0.0342681000474840 ); uniforms.push_back( 0.130203500390053 );
	uniforms.push_back( 0.96198905352503 ); uniforms.push_back( 0.151295251445845 );
	uniforms.push_back( 0.341955640818924 ); uniforms.push_back( 0.303998979507014 );
	uniforms.push_back( 0.126801372738555 ); uniforms.push_back( 0.470698264893144 );
	uniforms.push_back( 0.920687223551795 ); uniforms.push_back( 0.359546367544681 );
	uniforms.push_back( 0.868214765097946 ); uniforms.push_back( 0.000346498331055045 );
	uniforms.push_back( 0.32511946791783 ); uniforms.push_back( 0.404821055941284 );
	uniforms.push_back( 0.725521486485377 ); uniforms.push_back( 0.54769479460083 );
	uniforms.push_back( 0.861707922304049 ); uniforms.push_back( 0.097649555420503 );
	uniforms.push_back( 0.87782586272806 ); uniforms.push_back( 0.0868273449596018 );
	uniforms.push_back( 0.99655256420374 ); uniforms.push_back( 0.0535810000728816 );
	uniforms.push_back( 0.606473366729915 ); uniforms.push_back( 0.934530368074775 );
	uniforms.push_back( 0.654625355731696 ); uniforms.push_back( 0.596238830825314 );
	uniforms.push_back( 0.134928556857631 ); uniforms.push_back( 0.762741550104693 );
	uniforms.push_back( 0.105536070885137 ); uniforms.push_back( 0.359121092362329 );
	uniforms.push_back( 0.689746178453788 ); uniforms.push_back( 0.665619196835905 );
	uniforms.push_back( 0.0581016337964684 ); uniforms.push_back( 0.814952733926475 );
	uniforms.push_back( 0.329314556671306 ); uniforms.push_back( 0.0406567521858960 );
	uniforms.push_back( 0.0337584835942835 ); uniforms.push_back( 0.440054956125095 );
	uniforms.push_back( 0.628092234255746 ); uniforms.push_back( 0.867790650576353 );
	uniforms.push_back( 0.903888480970636 ); uniforms.push_back( 0.681713405996561 );
	uniforms.push_back( 0.398803260875866 ); uniforms.push_back( 0.223954246379435 );
	uniforms.push_back( 0.994349050102755 ); uniforms.push_back( 0.361061308067292 );
	uniforms.push_back( 0.841436631744727 ); uniforms.push_back( 0.58085672208108 );
	uniforms.push_back( 0.781777318101376 ); uniforms.push_back( 0.189785392023623 );
	uniforms.push_back( 0.814942180877551 ); uniforms.push_back( 0.695940729230642 );
	uniforms.push_back( 0.81591795058921 ); uniforms.push_back( 0.670843106228858 );
	uniforms.push_back( 0.574883400462568 ); uniforms.push_back( 0.933504132786766 );
}

Real myuniform() {
	// return numeric::random::uniform();
	Real back = uniforms.back();
	uniforms.pop_back();
	return back;
}

numeric::xyzVector<Real> randvec()
{
	numeric::xyzVector<Real> xyz(myuniform(),myuniform(),myuniform());
	while( xyz.length() > 1.0 )	{
		xyz = numeric::xyzVector<Real>(myuniform(),myuniform(),myuniform());
	}
	xyz.normalize();
	return xyz;
}

numeric::xyzMatrix<Real> rand_rot()
{
	using numeric::random::uniform;
	numeric::xyzVector<core::Real> axis(uniform(),uniform(),uniform());
	while( axis.length() > 1.0 ) axis = numeric::xyzVector<core::Real>(uniform(),uniform(),uniform());
	core::Real mag = uniform() * 2.0 * numeric::NumericTraits<Real>::pi();
	return numeric::rotation_matrix<core::Real>( axis, mag );
}

Real get_dots_area( core::scoring::packstat::Spheres spheres, Real probe_radius, Size N = 10, bool csa=false )
{
	Real area = 0.0;
	for( Size i = 1; i <= N; ++i ) {
		numeric::xyzMatrix<Real> rot = rand_rot();
		for( core::scoring::packstat::SphereIter s = spheres.begin(); s != spheres.end(); ++s ){
			numeric::xyzVector<Real> tmp = s->xyz;
			tmp = rot * tmp;
			s->xyz = tmp;
		}
		Real d = compute_sasa_generic( spheres, probe_radius, csa );
		area += d / (Real)N;
	}
	return area;
}

void test_known_value()
{
	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace packstat;

	PosePackDataOP pd = new PosePackData;

	pd->spheres.clear();
	for( int i = 0; i < 10; ++i )
		for( int j = 0; j < 10; ++j )
			for( int k = 0; k < 10; ++k )
				pd->spheres.push_back( Sphere( numeric::xyzVector<Real>( i*2.0, j*4.0 , k*4.0 ), 2.0 ) );

	for( Real pr = 0.0; pr < 3.01; pr += 0.1 ) {
		bool CSA = true;
		std::cerr << pr << " " << get_dots_area(pd->spheres,pr,1,CSA) << " " <<
			compute_surface_area_leerichards(pd,0.1,pr,CSA) << std::endl;
	}

}

void test_vs_dots( std::string fname )
{
	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace packstat;

  SimplePDB pdb;
  utility::io::izstream in(fname.c_str());
  in >> pdb;
	PosePackDataOP pd( pdb.get_pose_pack_data() );

	Reals slices;
	slices.push_back(1.0);	slices.push_back(0.5);	slices.push_back(0.2);
	slices.push_back(0.1);	slices.push_back(0.05);//	slices.push_back(0.02);	slices.push_back(0.01);
	for( Size i = 1; i <= slices.size(); ++i ) {
		for( Size j = 1; j <= 50; ++j ) {
			time_t t1 = clock();
			Real ps  = compute_packing_score_leerichards(pd,slices[i],randvec());
			t1 = clock() - t1;
			time_t t2 = clock();
			Real ps2 = compute_packing_score( *pd, i-1 );
			t2 = clock() - t2;
			std::cerr << "LeeRichards packing score: " << fname << " " << pd->spheres.size() << " " << pd->centers.size() << " "
								<< slices[i] << " " << ps  << " " << Real(t1)/1000000.0
								<< i-1       << " " << ps2 << " " << Real(t2)/1000000.0 << std::endl;
		}
	}

}


std::map<std::string,utility::io::ozstream*> outs;

struct WorkDat {
	core::pose::Pose *pose;
	core::scoring::packstat::MultiProbePoseAccumulator *accum;
	core::Real spacing,alpha;
	core::Size pri;
};

void *do_work(void *dat) {
	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace packstat;

	WorkDat *d = (WorkDat*)dat;
	Pose & pose(*(d->pose));
	MultiProbePoseAccumulator *accum = d->accum;
	Real spacing = d->spacing;
	Real alpha = d->alpha;
	Size pri = d->pri;

	accum->set_pr_idx(pri);
	time_t t = clock();
	LeeRichards *lr = new LeeRichards(pose,accum,spacing,alpha,true,randvec());
	while( accum->failed() ) {
		std::cerr << "LR failed trying again!" << std::endl;
		delete lr;
		accum->set_pr_idx(pri);
		lr = new LeeRichards(pose,accum,spacing,alpha,true,randvec());
	}
	delete lr;
	std::cout << "pr " << pri << " " << alpha << " " << sqrt(alpha*alpha+1.85*1.85)-1.85 << " " << Real(clock()-t)/1000000.0 << std::endl;
	// return NULL;
	delete d;
	pthread_exit(NULL);
}


void test( std::string fname )
{

	using namespace core;
	using namespace pose;
	using namespace scoring;
	using namespace packstat;
	using namespace basic::options;

	// Pose pose;
	// core::import_pose::pose_from_pdb(pose,fname);
	// PosePackDataOP pd = new PosePackData( pose_to_pack_data(pose) );

	if(1) {

		init_myuniform();

		std::string pdbid = fname.substr(fname.rfind('/')+1).substr(3,4);
		std::string outname = ("out/"+pdbid.substr(1,2)+"/"+pdbid);
		std::cerr << "checking for " << outname << std::endl;
		if( utility::file::file_exists(outname) ) {
			std::cerr << "file exists! skipping " << fname << std::endl;
			return;
		};

		Pose pose;
		std::cerr << fname << std::endl;
		core::import_pose::pose_from_pdb(pose,fname);
		// std::cerr << "test_leerichards: PDBInfo: " << pose.pdb_info()() << std::endl;
		// std::cerr << "test_leerichards: pose has " << pose.pdb_info()->get_unrecognized_atoms().size() << " unrecognized atoms" << std::endl;
		MultiProbePoseAccumulator *accum = new MultiProbePoseAccumulator(pose,fname);
		AccumulatorOP accumOP(accum);
		Real spacing = 0.05;

		Reals prs;
		prs.push_back(0.0); prs.push_back(0.1);	prs.push_back(0.2);	prs.push_back(0.3);	prs.push_back(0.4);
		prs.push_back(0.5);	prs.push_back(0.6);	prs.push_back(0.7);	prs.push_back(0.8);	prs.push_back(0.9);
		prs.push_back(1.0);	prs.push_back(1.1);	prs.push_back(1.2);	prs.push_back(1.4);	prs.push_back(1.6);
		prs.push_back(1.8);	prs.push_back(2.1);	prs.push_back(2.4);	prs.push_back(2.7);	prs.push_back(3.0);


		Reals alphas; // produce the probe growth schedule commented above if r0=1.85 ~ average radius
		// alphas.push_back(0.0000000); alphas.push_back(0.1001000); alphas.push_back(0.2001776); alphas.push_back(0.3003155);
		// alphas.push_back(0.4005605); alphas.push_back(0.5009956); alphas.push_back(0.6017685); alphas.push_back(0.7031414);
		// alphas.push_back(0.8055801); alphas.push_back(0.9099120); alphas.push_back(1.0176068); alphas.push_back(1.1312753);
		// alphas.push_back(1.2555547); alphas.push_back(1.3986827); alphas.push_back(1.5752916); alphas.push_back(1.8113733);
		// alphas.push_back(2.1530971); alphas.push_back(2.6824750); alphas.push_back(3.5451855); alphas.push_back(5.0000000);
		alphas.push_back(0.0000000); alphas.push_back(0.6167577); alphas.push_back(0.8835883); alphas.push_back(1.0960643);
		alphas.push_back(1.2816093); alphas.push_back(1.4507516); alphas.push_back(1.6089340); alphas.push_back(1.7595542);
		alphas.push_back(1.9051524); alphas.push_back(2.0480757); alphas.push_back(2.1910429); alphas.push_back(2.3378414);
		alphas.push_back(2.4943877); alphas.push_back(2.6704755); alphas.push_back(2.8827284); alphas.push_back(3.1596130);
		alphas.push_back(3.5499700); alphas.push_back(4.1377324); alphas.push_back(5.0680891); alphas.push_back(6.5954530);

		time_t ttot = clock();
		int num_threads = option[ OptionKeys::packstat::threads ]();
		if( num_threads ) {
			pthread_t threads[20];
		  	pthread_attr_t attr;
		  	pthread_attr_init(&attr);
		  	pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			for( Size pri = 1; pri <= prs.size(); ++pri ) {
				Real alpha = alphas[pri];
				if(pri > 10) spacing *= 1.1;
				WorkDat *d = new WorkDat;
				d->pose = &pose;
				d->accum = accum;
				d->spacing = spacing;
				d->alpha = alpha;
				d->pri = pri;
				// do_work(d);
				pthread_create(&threads[pri-1], &attr, do_work, (void *) d);
			}
			for(int i=0;i<20;++i) {
				pthread_join(threads[i],NULL);
			}
			pthread_attr_destroy(&attr);
			// pthread_exit(NULL);
		} else {
			for( Size pri = 1; pri <= prs.size(); ++pri ) {
				Real alpha = alphas[pri];
				accum->set_pr_idx(pri);
				time_t t = clock();
				// if( pri > 12 ) adj_spacing *= 2; if( pri > 15 ) adj_spacing *= 2; if( pri > 18 ) adj_spacing *= 2;
				if(pri > 10) spacing *= 1.1;
				LeeRichards *lr = new LeeRichards(pose,accum,spacing,alpha,true,randvec());
				while( accum->failed() ) {
					std::cerr << "LR failed trying again!" << std::endl;
					delete lr;
					accum->set_pr_idx(pri);
					lr = new LeeRichards(pose,accum,spacing,alpha,true,randvec());
				}
				delete lr;
				std::cout << "pr " << pri << " " << alpha << " " << sqrt(alpha*alpha+1.85*1.85)-1.85 << " " << Real(clock()-t)/1000000.0 << std::endl;
			}
		}
		// system(("mkdir -p out/"+pdbid.substr(1,2)).c_str());
		// utility::io::ozstream out( outname.c_str() );
		// accum->show( out );
		std::cerr << pdbid << " took " << Real(clock()-ttot)/Real(CLOCKS_PER_SEC) << "sec" << std::endl;
		std::string pref = basic::options::option[ basic::options::OptionKeys::out::prefix ]();
		accum->save_training_data( pref, outs );
		// out.close();
		std::cerr << "DONE with test" << std::endl;
		return;
	}

	if(1) {
	  SimplePDB pdb;
	  utility::io::izstream in(fname.c_str());
	  in >> pdb;
		PosePackDataOP pd = pdb.get_pose_pack_data();
		Real slice = 0.2;
		// bool CSA = true;
		// Real buried_area = 0.0;

		std::cerr << "packscore lr " << compute_packing_score_leerichards(pd,slice) << std::endl;
		return;
	}

	if( 1 ) {
		for( Real pr = 0.0; pr <= 3.01; pr += 0.2 ) {
		  SimplePDB pdb;
		  utility::io::izstream in(fname.c_str());
		  in >> pdb;
			PosePackDataOP pd = pdb.get_pose_pack_data();
			Real slice = 0.2;
			Real probe = pr;
			bool CSA = true;
			Real buried_area = 0.0;

			std::cerr << "packscore lr " << compute_packing_score_leerichards(pd,slice) << std::endl;

			Real area  = compute_surface_area_leerichards(buried_area,pd,slice,probe,CSA);
			// XYZs deriv = compute_surface_area_leerichards_deriv(pd,slice,probe,CSA);
			// XYZs check =   check_surface_area_leerichards_deriv(pd,slice,probe,CSA,9999);

			std::cerr << "PR " << pr << " " << area << " " << buried_area << std::endl;
			// for( Size i = 1; i <= check.size(); ++i ) {
			// 	std::cerr << i << " " << deriv[i].x() << " " << deriv[i].y() << " " << deriv[i].z()
			// 	           << "     " << check[i].x() << " " << check[i].y() << " " << check[i].z() << std::endl;
			// }
		}
		return;
	}

	PosePackDataOP pd = new PosePackData;//pdb.get_pose_pack_data();

	Real ST = 0.001, SP = 3.999;
	for( Real sep = ST; sep <= SP; sep += (SP-ST)/9.0 ) {

		pd->spheres.clear();
		Real SCALE = 1.0;
		pd->spheres.push_back( Sphere( XYZ( -sep/2.0, 0.0, 0.0 )*SCALE, 10.0/10.0*SCALE ) );
		pd->spheres.push_back( Sphere( XYZ(  sep/2.0, 0.0, 0.0 )*SCALE, 10.0/10.0*SCALE ) );
		// pd->spheres.push_back( Sphere( XYZ( 0, -sep/2.0, 0.0 )*SCALE, 1.0*SCALE ) );
		// pd->spheres.push_back( Sphere( XYZ( 0,  sep/2.0, 0.0 )*SCALE, 1.0*SCALE ) );
		// pd->spheres.push_back( Sphere( XYZ( 0.0, 0.0, -sep/2.0 )*SCALE, 1.0*SCALE ) );
		// pd->spheres.push_back( Sphere( XYZ( 0.0, 0.0,  sep/2.0 )*SCALE, 1.0*SCALE ) );

		Real slice = 0.001234567789;
		Real probe = 0.0;
		bool CSA = false;
		XYZ plane(0,0,1);

		Real area  = compute_surface_area_leerichards      (pd,slice,probe,CSA);
		XYZs deriv = compute_surface_area_leerichards_deriv(pd,slice,probe,CSA);
		XYZs check =   check_surface_area_leerichards_deriv(pd,slice,probe,CSA);

		// std::cerr << "area " << area << std::endl;
		// for( Size i = 1; i <= deriv.size(); i++ ) {
		// 	std::cerr << "deriv " << i << " " << deriv[i] << std::endl;
		// 	std::cerr << "check " << i << " " << check[i] << std::endl;
		// }

		std::cerr << "area " << area << std::endl;
		std::cerr << sep << " " << deriv[1].x() << " " << deriv[1].y() << " " << deriv[1].z()
		             << "     " << check[1].x() << " " << check[1].y() << " " << check[1].z() << std::endl;
		std::cerr << sep << " " << deriv[2].x() << " " << deriv[2].y() << " " << deriv[2].z()
		             << "     " << check[2].x() << " " << check[2].y() << " " << check[2].z() << std::endl;
	}

	// std::cerr << fname << " " << compute_packing_score_leerichards(pd,0.33,randvec()) << std::endl;


}

int
main (int argc, char *argv[])
{

	try {



	devel::init( argc, argv );

  using namespace basic::options;
  using namespace basic::options::OptionKeys;
  using namespace utility;

  // test_io();

	// test_sasa_dots();

	// core::pose::Pose native_pose;
	// core::import_pose::pose_from_pdb( native_pose, core::options::option[ core::options::OptionKeys::in::file::native ]() );

	if( option[ in::file::s ].user() ) {
  	vector1<file::FileName> files( option[ in::file::s ]() );
  	for( size_t i = 1; i <= files.size(); ++i ) {
    	test( files[i] );
  	}
	} else if( option[ in::file::l ].user() ) {
  		vector1<file::FileName> files( option[ in::file::l ]() );
  		for( size_t i = 1; i <= files.size(); ++i ) {
			utility::io::izstream list( files[i] );
			std::string fname;
			while( list >> fname ) {
				// std::cerr << "'" << fname << "'" << std::endl;
    		test( fname );
			}
  		}
	}

	for( std::map<std::string,utility::io::ozstream*>::iterator i = outs.begin(); i != outs.end(); ++i ) {
		i->second->close();
		delete i->second;
	}

	pthread_exit(NULL);
	return 0;



	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}
