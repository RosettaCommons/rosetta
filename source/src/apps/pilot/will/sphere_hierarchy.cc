#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
#include <basic/database/open.hh>
#include <core/conformation/symmetry/util.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/ScoreFileSilentStruct.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/sasa.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/packing/compute_holes_score.hh>
#include <numeric/conversions.hh>
#include <numeric/model_quality/rms.hh>
#include <numeric/random/random.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyz.io.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/FArray3D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/io/izstream.hh>
#include <utility/io/ozstream.hh>

#include <apps/pilot/will/will_util.ihh>

typedef numeric::xyzVector<Real> Vec;
typedef numeric::xyzMatrix<Real> Matf;
using core::id::AtomID;
using basic::options::option;
using core::pose::Pose;
using core::Real;
using core::scoring::ScoreFunctionOP;
using core::Size;
using numeric::max;
using numeric::min;
using numeric::random::gaussian;
using numeric::random::uniform;
using numeric::rotation_matrix_degrees;
using numeric::conversions::radians;
using numeric::conversions::degrees;
using ObjexxFCL::format::F;
using ObjexxFCL::format::I;
using ObjexxFCL::string_of;
using std::cerr;
using std::cout;
using std::string;
using utility::io::izstream;
using utility::io::ozstream;
using utility::vector1;
using std::endl;
using core::import_pose::pose_from_pdb;
using core::kinematics::Stub;
using core::conformation::ResidueOP;

// max ushort 65535
// 12  32  122  482  1922   7682  30722
// 12  44  166  648  2570  10252  40974


vector1<Size> getclosest(Vec const & v, vector1<Vec> const & sph) {
	Size NTOP = 7;
	vector1<Real> mnv(NTOP,9e9);
	vector1<Size> mni(NTOP,0);
	for(Size j = 1; j <= (Size)sph.size(); ++j) {		
		Real d = v.distance_squared(sph[j]);
		for(Size k = 1; k <= NTOP; k++) {
			if(mnv[k] > d) {
				for(Size l = NTOP; l > k; --l) {
					mnv[l] = mnv[l-1];
					mni[l] = mni[l-1];						
				}
				mnv[k] = d;
				mni[k] = j;
				break;
			}
		}
	}
	if( v.distance(sph[mni[mni.size()]]) > 1.5*v.distance(sph[mni[mni.size()-1]]) ) {
		mni[mni.size()] = 65535;
	}
	return mni;
}

int levelof(int i) {
	if(i <=    12) return 1;
	if(i <=    44) return 2;
	if(i <=   166) return 3;
	if(i <=   648) return 4;
	if(i <=  2570) return 5;
	if(i <= 10252) return 6;
	if(i <= 40974) return 7;
	return 65535;
}

int main (int argc, char const *argv[])
{

	vector1<vector1<Vec> > coord(7);
	for(int i = 1; i <= 7; ++i){
		vector1<int> sizes(7,12); sizes[2]=32; sizes[3]=122; sizes[4]=482; sizes[5]=1922; sizes[6]=7682; sizes[7]=30722;
		coord[i].resize(sizes[i]);
		izstream is("/Users/sheffler/Dropbox/project/sphere_hierarchy/att/sphere_"+str(sizes[i])+".dat");
		for( int j = 1; j <= sizes[i]; ++j) {
			is >> coord[i][j].x() >> coord[i][j].y() >> coord[i][j].z();
		}
		is.close();
	}

	// {
	// 	vector1<Size> nbr = getclosest(coord[1][1],coord[1]);
	// 	for(int jj = 1; jj <= 7; ++jj) {
	// 		cout << "nbr1 " << jj << " " << nbr[jj] << " " << coord[1][1].distance(coord[1][nbr[jj]]) << endl;
	// 	}
	// 	utility_exit_with_message("aoristn");
	// }

	vector1<int> level(coord[7].size(),7);
	for(int i = 6; i > 0; --i) {
// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif
		for(int j = 1; j <= coord[i].size(); ++j) {
			Real mn = 9e9; int mnk;
			for(int k = 1; k <= coord[7].size(); ++k) {
				Real d2 = coord[7][k].distance_squared( coord[i][j] );
				if( d2 < mn ) {
					mn = d2;
					mnk = k;
				}
			}
			level[mnk] = i;
		}
	}
	// {
	// 	vector1<int> count(7,0);
	// 	for(int i = 1; i <= level.size(); ++i) {
	// 		cerr << level[i] << endl;
	// 		count[level[i]]++;
	// 	}
	// 	cout << "level 1: " << count[1] << endl;
	// 	cout << "level 2: " << count[2]+count[1] << endl;
	// 	cout << "level 3: " << count[3]+count[2]+count[1] << endl;
	// 	cout << "level 4: " << count[4]+count[3]+count[2]+count[1] << endl;
	// 	cout << "level 5: " << count[5]+count[4]+count[3]+count[2]+count[1] << endl;
	// 	cout << "level 6: " << count[6]+count[5]+count[4]+count[3]+count[2]+count[1] << endl;					
	// 	cout << "level 7: " << count[7]+count[6]+count[5]+count[4]+count[3]+count[2]+count[1] << endl;					
	// }

	vector1<vector1<Vec> > lcoord(7);
	for(int l = 1; l <= 7; ++l) {
		for(int i = 1; i <= level.size(); ++i) {
			if(level[i] > l) continue;
			lcoord[l].push_back(coord[7][i]);
		}
		cout << " " << lcoord[l].size();
	}
	cout << endl;
	
	vector1<Vec> hscoord   (40974  ,65535);
	vector1<int> hslevel   (40974  ,65535);
	vector1<int> hsparent  (40974  ,65535);
	vector1<int> hsneighbor(40974*6,65535);
	vector1<int> hschildren(10252*7,65535); // no children for level 7
	vector1<int> lstart(7,0); lstart[1]=0; lstart[2]=12; lstart[3]=44; lstart[4]=166; lstart[5]=648; lstart[6]=2570; lstart[7]=10252;
	int count = 0;
	for(Size l = 1; l <= 7; ++l) {
		cout << l << endl;
		int lcount = 0;
// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif		
		for(Size i = 1; i <= lcoord[l].size(); ++i) {
// #ifdef USE_OPENMP
// #pragma omp critical
// #endif
			count++;
			// level
			hslevel[count] = l;
			// crd
			hscoord[count] = lcoord[l][i];
			// nbr
			vector1<Size> nbr = getclosest(hscoord[count],lcoord[l]);
			// if(l==1) {
			// 	for(int jj = 1; jj <= 7; ++jj)
			// 	cout << "nbr1 " << jj << " " << nbr[jj] << " " << hscoord[count].distance(lcoord[1][jj]) << endl;
			// }
			for(int n = 1; n <= 6; ++n) {
				if( n==6 && nbr[n+1] == 65535)	{
					hsneighbor[6*(count-1)+n] = 65535;					
				} else if(nbr[n+1] == 65535) {
					utility_exit_with_message("getclosest gave 65535 for nbr "+str(n+1));
				} else {
					hsneighbor[6*(count-1)+n] = nbr[n+1]+lstart[l];
				}
			}
			if(l == 7) continue;
			// child / parent
			vector1<Size> child = getclosest(hscoord[count],lcoord[l+1]);
			for(int n = 1; n <= 7; ++n) {
				if( n==7 && child[n] == 65535)	{
					hschildren[7*(count-1)+n] = 65535;					
				} else if(child[n] == 65535) {
					utility_exit_with_message("getclosest gave 65535 for child "+str(n));
				} else {
					hschildren[7*(count-1)+n] = child[n]+lstart[l+1];
					hsparent[child[n]+lstart[l+1]] = count;
				}
			}			
		}
	}
	
	// sanity checks
	vector1<int> nl(7,0),n5(7,0),n6(7,0);
// #ifdef USE_OPENMP
// #pragma omp parallel for schedule(dynamic,1)
// #endif	
	for(Size i = 1; i <= 40974; ++i) {
		int l  = hslevel[i];
		nl[l]++;
		for(int n = 1; n <= 6; ++n) {
			if( n==6 && hsneighbor[6*(i-1)+n]==65535 ) {
				n5[l]++;
			} else if( levelof( hsneighbor[6*(i-1)+n] ) != l ) {
				cout << i << " " << l << " " << n << " " << " " << hsneighbor[6*(i-1)+n] << " " << levelof(hsneighbor[6*(i-1)+n]) << endl;
				utility_exit_with_message("bad neighbor!");
			}
		}
		
		// if( l != 1 && levelof(hsparent[i]) != l-1 ) {
		// 	cout << i << " " << l << " " << hsparent[i] << " " << levelof(hsparent[i]) << endl;
		// 	utility_exit_with_message("bad parent!");
		// }
		if(l == 7) continue;
		for(int n = 1; n <= 7; ++n) {
			if( n==7 && hschildren[7*(i-1)+n]==65535 ) {
				n6[l]++;
			} else if( levelof( hschildren[7*(i-1)+n] ) != l+1 ) {
				cout << i << " " << l << " " << n << " " << hschildren[7*(i-1)+n] << " " << levelof(hschildren[7*(i-1)+n]) << endl;
				utility_exit_with_message("bad child!");
			}
		}
	}	
	for(int i = 1; i <= 7; ++i) {
		cout << "nl " << i << " " << nl[i] << endl;
		if( n5[i] != 12 ) utility_exit_with_message("bad n5! "+str(i)+" " +str(n5[i]));
		if( i != 7 && n6[i] != 12 ) utility_exit_with_message("bad n6! "+str(i)+" "+str(n6[i]));
	}
	cout << "sanity checks pass!" << endl;
	cout << "dumping data file hsphere.dat.gz" << endl;
	
	utility::io::ozstream out("hsphere.dat.gz");
	for(Size i = 1; i <= 40974; ++i) {
		out << hslevel[i] << " " << hsparent[i] << " ";
		out << F(12,10,hscoord[i].x()) << " " << F(12,10,hscoord[i].y()) << " " << F(12,10,hscoord[i].z()) << endl;
		if(hslevel[i]<7) {
			for(Size j = 1; j <= 7; ++j) {
				out << " " << hschildren[7*(i-1)+j];
			}
			out << endl;
		} else {
			out << "65535" << endl;
		}
		for(Size j = 1; j <= 6; ++j) {
			out << " " << hsneighbor[6*(i-1)+j];
		}
		out << endl;
	}
	out.close();

	return 0;
}


