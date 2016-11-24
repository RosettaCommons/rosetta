#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>
#include <core/pose/Pose.hh>
#include <numeric/random/random.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>
#include <protocols/sic_dock/Rose.hh>
#include <sstream>
#include <numeric/xyzTransform.hh>

using std::string;
using utility::vector1;
using ObjexxFCL::format::I;
using ObjexxFCL::format::F;
using ObjexxFCL::format::RJ;
using numeric::min;
using numeric::max;
using std::cout;
using std::cerr;
using std::endl;
using namespace numeric;
using namespace numeric::random;
using std::cout;
using std::endl;
using core::pose::Pose;
using core::pose::PoseOP;
using basic::options::option;
using namespace basic::options::OptionKeys;
using protocols::sic_dock::Rose;

typedef xyzVector<core::Real>    V;
typedef xyzMatrix<core::Real>    M;
typedef xyzTransform<core::Real> X;

core::Real const eps = 0.000001;

int main(int argv, char **argc){

  try {

	devel::init(argv,argc);

	for(core::Size i = 1; i <= 100000; ++i){
		X I,x;
		// std::istringstream iss("1 -9 0\n 0 1 3\n 0 0 1\n 1 2 3\n");
		// iss >> x;

		x   = X( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()), V(gaussian(),gaussian(),gaussian()) );
		X y = X( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()), V(gaussian(),gaussian(),gaussian()) );
		X z = X( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()), V(gaussian(),gaussian(),gaussian()) );
		V u(gaussian(),gaussian(),gaussian());
		V v(gaussian(),gaussian(),gaussian());
		M m( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()) );
		M n( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()) );

		if((   x             ).distance(     ~~x           ) > eps ) utility_exit_with_message("x != ~~x!");
		if((   x*v           ).distance(    (~~x)*v        ) > eps ) utility_exit_with_message("v!= ~~x*v!");
		if((   (m*x)*v       ).distance(   m*(x*v)         ) > eps ) utility_exit_with_message("((m*x)*v != m*(x*v)");
		if((   (x*m)*v       ).distance(   x*(m*v)         ) > eps ) utility_exit_with_message("((x*m)*v != x*(m*v)");
		if((        I        ).distance(   ~x*x             ) > eps ) utility_exit_with_message("x/x != I");
		if((   m*(u-v)+v     ).distance(   X::rot(m,v)*u   ) > eps ) utility_exit_with_message( "m*(u-v)+v != X::rot(m,v)*u");
		if((   x*( y*z)      ).distance(   (x*y)*z         ) > eps ) utility_exit_with_message( "fail");
		if((   x*( y*v)      ).distance(   (x*y)*v         ) > eps ) utility_exit_with_message( "fail");
		// if((  ((~y)*x)*v     ).distance(   (x/y)*v         ) > eps ) utility_exit_with_message( "((~y)*x)*v) != (x/y)*v");
		if((   x*m           ).distance(   X(x.R*m,x.t)    ) > eps ) utility_exit_with_message( "fail");

	}
	cout << "xyzTransform<core::Real> tests pass" << endl;

	if( option[in::file::s]().size() < 2 ) return 0;

	vector1<PoseOP> poses = core::import_pose::poseOPs_from_files(option[in::file::s](), core::import_pose::PDB_file);
	Rose r(poses[1]);
	Rose s(poses[2]);
	// bool dumpcl=true,dumpnc=true;
	for(core::Size i = 1; i <= 10000; ++i){
		X x = X( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()), 30*V(gaussian(),gaussian(),gaussian()) );
		X y = X( rotation_matrix(V(gaussian(),gaussian(),gaussian()),gaussian()), 30*V(gaussian(),gaussian(),gaussian()) );
		Rose const rx = x*r;
		Rose const sy = y*s;
		bool const c1 = (rx).clashes( sy );
		core::Size const c2 = (sy).contacts_naive( rx );
		if( c1 != (0<c2) ) utility_exit_with_message("clash fail");
		// if(c2>0 && c2<10 && dumpcl){
		// 	rx.dump_pdb("clash1_test.pdb");
		// 	sy.dump_pdb("clash2_test.pdb");
		// 	dumpcl=false;
		// }
		// if(!c1 && dumpnc){
		// 	rx.dump_pdb("noclash1_test.pdb");
		// 	sy.dump_pdb("noclash2_test.pdb");
		// 	dumpnc=false;
		// }
	}
	cout << "Rose tests pass" << endl;


	return 0;

  } catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
  }

}
