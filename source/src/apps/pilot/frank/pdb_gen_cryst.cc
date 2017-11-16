/// @file
/// @brief

#include <apps/pilot/frank/spacegroup.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/CrystInfo.hh>
#include <core/id/AtomID.hh>
#include <core/id/AtomID_Map.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/pdb/pdb_writer.hh>
#include <core/types.hh>
#include <core/pack/task/ResfileReader.hh>

#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/moves/MoverContainer.hh>

#include <utility/vector1.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/mpistream.hh>
#include <utility/string_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/file_sys_util.hh>

#include <numeric/random/random.hh>
#include <numeric/NumericTraits.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/constants.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>

#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>

#include <devel/init.hh>

#include <utility/string_constants.hh>
#include <utility/excn/Exceptions.hh>

#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/option_macros.hh>
#include <basic/basic.hh>
#include <basic/database/open.hh>

// option includes
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <fstream>
#include <iostream>
#include <math.h>
#include <cmath>

#include <sstream>
#include <string>
#include <queue>
#include <cstdarg>

using namespace basic;
using namespace core;
using namespace core::pose;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;

#define DEG2RAD 0.0174532925199433
#define RAD2DEG 57.295779513082323

////////////////////////////////////////////////
// helper functions
inline int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
inline Real pos_mod(Real x,Real y) {
	Real r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline int min_mod(int x,int y) {
	int r=x%y; if ( r<-y/2 ) { r+=y; } if ( r>=y/2 ) { r-=y; }
	return r;
}
inline double min_mod(double x,double y) {
	double r=std::fmod(x,y); if ( r<-0.5*y ) { r+=y; } if ( r>=0.5*y ) { r-=y; }
	return r;
}

static basic::Tracer TR( "cryst.gen" );

//missing string printf
//this is safe and convenient but not exactly efficient
inline std::string formatstr(const char* fmt, ...){
	int size = 512;
	char* buffer = 0;
	buffer = new char[size];
	va_list vl;
	va_start(vl, fmt);
	int nsize = vsnprintf(buffer, size, fmt, vl);
	if ( size<=nsize ) { //fail delete buffer and try again
		delete[] buffer;
		buffer = 0;
		buffer = new char[nsize+1]; //+1 for /0
		nsize = vsnprintf(buffer, size, fmt, vl);
	}
	std::string ret(buffer);
	va_end(vl);
	delete[] buffer;
	return ret;
}


struct pdbline {
	std::string line;
	numeric::xyzVector<core::Real> X;

	pdbline( std::string linein ){
		line = linein;
		X[0] = atof(line.substr(30,8).c_str());
		X[1] = atof(line.substr(38,8).c_str());
		X[2] = atof(line.substr(46,8).c_str());
	}

	std::string
	getline(char chain) {
		line.replace( 30,8, formatstr( "%8.3f", X[0] ));
		line.replace( 38,8, formatstr( "%8.3f", X[1] ));
		line.replace( 46,8, formatstr( "%8.3f", X[2] ));
		line[21]=chain;
		return line;
	}
};

// quick and dirty PDB read where we only care about heavyatom locations and atom ids
void
readPDBcoords(std::string filename, Spacegroup &sg, utility::vector1<pdbline> &pdblines, utility::vector1<Vector> &CAs) {
	std::ifstream inpdb(filename.c_str());
	std::string buf;

	while ( std::getline(inpdb, buf ) ) {
		if ( buf.substr(0,6)=="CRYST1" ) {
			sg.set_spacegroup( buf.substr(55,10) );
			sg.set_parameters(
				atof(buf.substr(7,8).c_str()), atof(buf.substr(16,8).c_str()), atof(buf.substr(25,8).c_str()),
				atof(buf.substr(34,6).c_str()), atof(buf.substr(41,6).c_str()), atof(buf.substr(47,6).c_str()) );
			continue;
		}

		if ( buf.substr(0,6)!="ATOM  " && buf.substr(0,6)!="HETATM" ) continue;
		if ( buf.substr(16,1)!=" " && buf.substr(16,1)!="A" ) continue;
		pdblines.push_back( pdbline(buf) );

		if ( buf.substr(12,4) == " CA " ) {
			CAs.push_back( pdblines[ pdblines.size() ].X );
		}
	}
}

void
apply_transform( numeric::xyzMatrix<Real> R, numeric::xyzVector<Real> T, utility::vector1<pdbline> &pdblines, utility::vector1<Vector> &CAs) {
	for ( int i=1; i<=(int)pdblines.size(); ++i ) {
		pdblines[i].X = R*pdblines[i].X + T;
	}
	for ( int i=1; i<=(int)CAs.size(); ++i ) {
		CAs[i] = R*CAs[i] + T;
	}
}

static std::string chains = utility::ALPHANUMERICS + "!@#$%^&*() ";

int
main( int argc, char * argv [] ) {
	try {
		Real contact_dist=12.0;
		std::string pdbfile;

		if ( argc<=1 ) {
			std::cerr << "USAGE: " << argv[0] << " pdbfile [clashdist]" << std::endl;
			exit (1);
		} else {
			pdbfile = argv[1];
			if ( argc>=3 ) contact_dist = atof(argv[2]);
		}

		// load pose using minimal reader
		utility::vector1<pdbline> pdblines;
		utility::vector1<Vector> CAs;
		Spacegroup sg;
		readPDBcoords( pdbfile, sg, pdblines, CAs );

		Vector com(0,0,0);
		Size nres = CAs.size();

		if ( nres == 0 ) {
			std::cerr << "Error reading PDB information from " << pdbfile << std::endl;
			exit(1);
		}
		for ( Size i=1; i<= nres; ++i ) com += CAs[i];
		com /= nres;

		Real mindis2(1e6);
		for ( Size i=1; i<= nres; ++i ) {
			Real const dis2( com.distance_squared(  CAs[i] ) );
			if ( dis2 < mindis2 ) {
				mindis2 = dis2;
			}
		}
		Size nsymm = sg.nsymmops();
		Size bestxform=0;
		Vector bestoffset(0,0,0);
		mindis2=1e6;
		com = sg.c2f()*com;
		for ( Size i=1; i<=nsymm; ++i ) {
			Vector foffset = sg.symmop(i).get_rotation()*com + sg.symmop(i).get_translation(), rfoffset;
			rfoffset[0] = min_mod( foffset[0], 1.0 );
			rfoffset[1] = min_mod( foffset[1], 1.0 );
			rfoffset[2] = min_mod( foffset[2], 1.0 );
			Real dist = (sg.f2c()*rfoffset).length_squared();
			if ( dist<mindis2 ) {
				mindis2=dist;
				bestxform=i;
				bestoffset = foffset - rfoffset;
			}
		}
		numeric::xyzMatrix<Real> Rx = sg.f2c()*sg.symmop(bestxform).get_rotation()*sg.c2f();
		numeric::xyzVector<Real> Tx = sg.f2c()*(sg.symmop(bestxform).get_translation() - bestoffset);
		apply_transform( Rx,Tx, pdblines, CAs );

		// write main chain
		int nchains = chains.length();
		int currchain = 0;
		for ( int i=1; i<=(int)pdblines.size(); ++i ) {
			std::cout <<pdblines[i].getline(chains[currchain])<<std::endl;
		}
		currchain = (currchain+1)%nchains;

		//// find contacts
		Real radius = 0;
		for ( Size i=1; i<= nres; ++i ) {
			radius = std::max( (CAs[i]).length_squared() , radius );
		}
		radius = sqrt(radius);

		for ( int s=1; s<=(int)sg.nsymmops(); ++s ) {
			numeric::xyzMatrix<Real> R_i = sg.symmop(s).get_rotation();

			for ( int a=-1; a<=1; ++a ) {
				for ( int b=-1; b<=1; ++b ) {
					for ( int c=-1; c<=1; ++c ) {
						if ( s==1 && a==0 && b==0 && c==0 ) continue;

						numeric::xyzVector<Real> T_i = sg.symmop(s).get_translation() + numeric::xyzVector<Real>(a,b,c);

						// pass 1 check vrt-vrt dist to throw out very distant things
						Real disVRT = (sg.f2c()*T_i).length();
						if ( disVRT>contact_dist+2*radius ) continue;

						// pass 2 check ca-ca dists
						numeric::xyzMatrix<core::Real> R_i_realspace =  sg.f2c()*R_i*sg.c2f();
						bool contact=false;
						for ( Size j=1; j<= nres && !contact; ++j ) {
							Vector Xi = R_i_realspace*CAs[j] + sg.f2c()*T_i;
							for ( Size k=1; k<= nres && !contact; ++k ) {
								if ( (Xi-CAs[k]).length_squared() < contact_dist*contact_dist ) {
									contact=true;
								}
							}
						}

						if ( contact ) {
							utility::vector1<pdbline> pdblines_i = pdblines;
							utility::vector1<Vector> CAs_i; // dummy
							apply_transform( R_i_realspace, sg.f2c()*T_i, pdblines_i, CAs_i );

							for ( int i=1; i<=(int)pdblines_i.size(); ++i ) {
								std::cout <<pdblines_i[i].getline(chains[currchain])<<std::endl;
							}
							currchain = (currchain+1)%nchains;
						}
					}
				}
			}
		}


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

