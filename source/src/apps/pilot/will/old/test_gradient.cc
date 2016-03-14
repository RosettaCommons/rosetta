// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.


// option key includes

#include <basic/options/keys/holes.OptionKeys.gen.hh>


vector1<Real>
get_area( vector1<Ball> const & balls ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	vector1<Real> area( balls.size() );
	redi::pstream proc( string( basic::options::option[ holes::dalphaball ]() ) );

	proc << "NPOINTS" << std::endl << balls.size() << std::endl << "COORDS" << std::endl;
	for( vector1<Ball>::const_iterator i = balls.begin(); i != balls.end(); ++i ) {
		proc << i->x() << " " << i->y() << " " << i->z() << " " << i->r() << " " << std::endl;
	}
	proc << "WEIGHTS" << std::endl;
	for( Size i = 1; i <= balls.size(); i++ ) {
		// for( Size j = 1; j <= 20; j++ ) {
		// 	proc << params.param("TRP")[24*(i-1)+4+j] << " ";
		// }
		// proc << std::endl;
		proc << "0.05509491 -0.38881221  1.35370724  0.68471421 -0.62791115  1.86401631  0.58868717 -1.17232821 -0.31347622  0.51271984  1.84685186 -0.16408191 -0.03606228 -0.63129269 -0.72709036  0.34628173  0.15434081  0.57971898  1.20265581 -0.74419286" << std::endl;
	}
	proc << "END" << std::endl << redi::peof;

	// int index,ialpha;
	Real val;
	for( Size a = 1; a <= 20; a++ ) {
		for( vector1<Real>::iterator i = area.begin(); i != area.end(); ++i ) {
			proc /*>> ialpha >> index*/ >> val;
			(*i) += val;
			// std::cerr << val << std::endl;
		}
	}

	return area;
}

Real get_area_tot( vector1<Ball> const & balls ) {
	vector1<Real> area = get_area(balls);
	Real sum = 0;
	for( vector1<Real>::const_iterator i = area.begin(); i != area.end(); ++i ) sum += *i;
	return sum;
}

vector1<xyzVector<Real> >
get_deriv( vector1<Ball> const & balls ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	vector1< xyzVector<Real> > deriv( balls.size(), xyzVector<Real>(0,0,0) );
	redi::pstream proc( string(option[ holes::dalphaball ]()) + " DERIV" );

	proc << "NPOINTS" << std::endl << balls.size() << std::endl << "COORDS" << std::endl;
	for( vector1<Ball>::const_iterator i = balls.begin(); i != balls.end(); ++i ) {
		proc << i->x() << " " << i->y() << " " << i->z() << " " << i->r() << " " << std::endl;
	}
	proc << "WEIGHTS" << std::endl;
	for( Size i = 1; i <= balls.size(); i++ ) {
		// for( Size j = 1; j <= 20; j++ ) {
		// 	proc << params.param("TRP")[24*(i-1)+4+j] << " ";
		// }
		// proc << std::endl;
		proc << "0.05509491 -0.38881221  1.35370724  0.68471421 -0.62791115  1.86401631  0.58868717 -1.17232821 -0.31347622  0.51271984  1.84685186 -0.16408191 -0.03606228 -0.63129269 -0.72709036  0.34628173  0.15434081  0.57971898  1.20265581 -0.74419286" << std::endl;
	}
	proc << "END" << std::endl << redi::peof;

	// int index,ialpha;
	Real /*val,*/dx,dy,dz;
	for( Size a = 1; a <= 20; a++ ) {
		for( vector1< xyzVector<Real> >::iterator i = deriv.begin(); i != deriv.end(); ++i ) {
				// char buf[999];
				// proc.getline(buf,999);
				// std::istringstream iss(buf);
				// std::cerr << "DAlphaBall: " << iss.str() << std::endl;
				// iss >> ialpha >> index >> val >> dx >> dy >> dz;
				// std::cerr << "Read in: " << dx << " " << dy << " " << dz << std::endl;
				// i->x() += dx; i->y() += dy; i->z() += dz;
				// std::cerr << "Read in: " << i->x() << " " << i->y() << " " << i->z() << std::endl;
				proc /*>> ialpha >> index >> val*/ >> dx >> dy >> dz;
				i->x() += dx; i->y() += dy; i->z() += dz;
		}
	}

	return deriv;
}

void
test_deriv( Real d, Real & max_mag_diff, Real & max_frac_diff ) {
	using numeric::random::gaussian;

	vector1<Ball> balls;
	for( Size i = 1; i <= 10; i++ ) {
		balls.push_back( Ball(gaussian(),gaussian(),gaussian(),1.0) );
	}
	// balls.push_back( Ball(gaussian(),gaussian(),gaussian(),1.0) );
	// balls.push_back( Ball(-0.196234188151175,-1.28970865616930 , 0.31406794256922  , 1.0 ) );
	// balls.push_back( Ball(-1.64685693841249 , 0.296228386899733, 0.0756031294018824, 1.0 ) );
	// balls.push_back( Ball( 0.939475667326099,-0.252929547228483,-0.59227524517714  , 1.0 ) );
	// balls.push_back( Ball(-1.26848083466054 ,-0.71866424678762 ,-0.0297940702068403, 1.0 ) );

	// vector1<Real> area = get_area(balls);
	// for( Size i = 1; i <= 10; i++ ) {
	// 	std::cout << "area " << i << " " << area[i] << std::endl;
	// }

	vector1<xyzVector<Real> > deriv = get_deriv(balls);
	// for( Size i = 1; i <= 10; i++ ) {
	// 	std::cout << "area " << i << " " << deriv[i].x() << " " << deriv[i].y() << " " << deriv[i].z() << std::endl;
	// }

	Real area0 = get_area_tot(balls);
	std::cout << area0 << " ";
	vector1<xyzVector<Real> > nderiv(balls.size());
	for( Size i = 1; i <= balls.size(); i++ ) {
		balls[i].x() += d; nderiv[i].x() = (get_area_tot(balls)-area0) / d; balls[i].x() -= d;
		balls[i].y() += d; nderiv[i].y() = (get_area_tot(balls)-area0) / d; balls[i].y() -= d;
		balls[i].z() += d; nderiv[i].z() = (get_area_tot(balls)-area0) / d; balls[i].z() -= d;
	}

	max_mag_diff = 0.0;
	max_frac_diff = 0.0;
	for( Size i = 1; i <= balls.size(); i++ ) {
		// std::cout << "d_area " << i << " x " << deriv[i].x() << " " << nderiv[i].x() << std::endl;
		// std::cout << "d_area " << i << " y " << deriv[i].y() << " " << nderiv[i].y() << std::endl;
		// std::cout << "d_area " << i << " z " << deriv[i].z() << " " << nderiv[i].z() << std::endl;
		Real  magdiff = (nderiv[i]-deriv[i]).length();
		Real fracdiff = magdiff / ( nderiv[i].length()+deriv[i].length() );
		if(  magdiff > max_mag_diff  ) max_mag_diff = magdiff;
		if( fracdiff > max_frac_diff ) max_frac_diff = fracdiff;
	}

	// std::cout << "max numeric / analytical diff " << max_mag_diff << std::endl;

}

numeric::xyzVector<Real>
get_numeric_deriv( id::AtomID aid, pose::Pose & pose, PoseBalls & pb, HolesParams & params, Real score0 ) {
	pb.ball(aid).xyz() += xyzVector<Real>(0.0001,0,0);
	Real ndx = (get_holes_score(pose,pb,params) - score0) / 0.0001;
	pb.ball(aid).xyz() -= xyzVector<Real>(0.0001,0,0);
	pb.ball(aid).xyz() += xyzVector<Real>(0,0.0001,0);
	Real ndy = (get_holes_score(pose,pb,params) - score0) / 0.0001;
	pb.ball(aid).xyz() -= xyzVector<Real>(0,0.0001,0);
	pb.ball(aid).xyz() += xyzVector<Real>(0,0,0.0001);
	Real ndz = (get_holes_score(pose,pb,params) - score0) / 0.0001;
	pb.ball(aid).xyz() -= xyzVector<Real>(0,0,0.0001);
	return xyzVector<Real>(ndx,ndy,ndz);
}

void
test_deriv_pose( pose::Pose & pose ) {

	HolesParams params;

	PoseBalls pb(pose);

	core::id::AtomID_Map< numeric::xyzVector<core::Real> > derivs;
	core::id::initialize_heavy_only( derivs, pose, xyzVector<Real>(0,0,0) );

	Real score0 = get_holes_score(pose,pb,params);
	get_holes_deriv(pose,pb,params,derivs);

	std::cout << "score: " << score0 << " " << score0 / (pose.total_residue()-2) << std::endl;

	for( int r = 2; r <= (int)pose.total_residue()-1; r++ ) {
		for( int i = 1; i <= 4; i++ ) {
			id::AtomID aid(i,r);
			xyzVector<Real> nderiv = get_numeric_deriv( aid, pose, pb, params, score0 );
			Real fracdiff = (derivs[aid]-nderiv).length() / ( nderiv.length() + derivs[aid].length() );
			std::cerr << " deriv " << r << " " << i << " " << fracdiff << " "
						 << derivs[aid].x() << " " << nderiv.x() << " "
						 << derivs[aid].y() << " " << nderiv.y() << " "
						 << derivs[aid].z() << " " << nderiv.z() << " "
						 << std::endl;
		}
	}
}

core::Real
balls_rms( PoseBalls & a, PoseBalls & b ) {
	ObjexxFCL::FArray2D<Real> fa(3,a.nballs());
	ObjexxFCL::FArray2D<Real> fb(3,a.nballs());
	for( Size i = 1; i<= a.nballs(); ++i ) {
		fa(1,i) = a.ball(i).x();
		fa(2,i) = a.ball(i).y();
		fa(3,i) = a.ball(i).z();
		fb(1,i) = b.ball(i).x();
		fb(2,i) = b.ball(i).y();
		fb(3,i) = b.ball(i).z();
	}
	return numeric::model_quality::rms_wrapper(a.nballs(),fa,fb);
}

void
balls_pdb( PoseBalls & pb, std::string fname ) {
	std::ofstream out(fname.c_str());
	for( Size i = 1; i <= pb.nballs(); i++ ) {
    	out << "HETATM" + I( 5, ( min( 9999, (int)i) ) ) + " "+string("ATOM")+" "+string("RES")+" Z"
		+ I( 4, pb.res_num(i) ) + "    "
		+ F( 8, 3, pb.ball(i).x() ) + F( 8, 3, pb.ball(i).y() ) + F( 8, 3, pb.ball(i).z() )
		+ F( 6, 2, 1.0 ) + ' ' + F( 5, 2, 1.0 );
		out << std::endl;
	}
	out.close();
}

void
test_gradient() {
	using namespace pose;
	using namespace io::pdb;
	using namespace core::id;

	core::scoring::packing::HolesParams params;

	Pose bad,good;
	core::import_pose::pose_from_file(bad ,"hivp_low_1rv7.pdb.clean.align", core::import_pose::PDB_file);
	core::import_pose::pose_from_file(good,"hivp_high_1tw7.pdb.clean.align", core::import_pose::PDB_file);

	PoseBalls bbad(bad);
	PoseBalls bgood(good);
	PoseBalls bbad0(bad);
	PoseBalls bgood0(good);

	std::cerr << "rms " << balls_rms(bbad,bgood) << " nres " << good.total_residue() << " " << bad.total_residue() << std::endl;

	AtomID_Map<xyzVector<Real> > dbad,dgood;
	initialize_heavy_only(dbad ,bad ,xyzVector<Real>(0,0,0));
	initialize_heavy_only(dgood,good,xyzVector<Real>(0,0,0));

	Real bscore=0.0,gscore=0.0;
	for( int iter = 1; iter <= 10000; iter++ ) {
		Real SCALE = min(1/100.0,0.2/Real(iter));
		Real bscore = get_holes_deriv( bad , bbad , params, dbad );
		Real gscore = get_holes_deriv( good, bgood, params, dgood );
		xyzVector<Real> G;
		for( Size b = 1; b <= bbad.nballs(); b++ ) {
			G = dbad[bbad.index_to_id(b)] * SCALE;
			if( G.length() > 0.1 ) G /= (G.length()*10.0);
			bbad.ball(b).xyz() -= G;
			bbad.ball(b).xyz() += (bbad0.ball(b).xyz() - bbad.ball(b).xyz()) *SCALE*10.0;
			G = dgood[bgood.index_to_id(b)] * SCALE;
			if( G.length() > 0.1 ) G /= (G.length()*10.0);
			bgood.ball(b).xyz() -= G;
			bgood.ball(b).xyz() += (bgood0.ball(b).xyz() - bgood.ball(b).xyz()) *SCALE*10.0;
			// bgood.ball(b).xyz() -= dgood[bgood.index_to_id(b)]/100.0;
		}
		std::cerr << "rnd " << iter << " score " << bscore << " " << gscore
		          << " rms " << balls_rms(bbad0,bbad) << " " << balls_rms(bbad,bgood) << " " << balls_rms(bgood,bgood0) << std::endl;
		// balls_pdb(bbad,"bad"+string_of(iter)+".pdb");
	}
	bscore = get_holes_score( bad, bbad, params );
	std::cerr << "final score " << bscore << " rms " << balls_rms(bgood,bbad) << std::endl;

}
