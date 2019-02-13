#include <utility/excn/Exceptions.hh>
#include <devel/init.hh>

#include <protocols/network/cloud.hh>

#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>

#include <basic/options/option_macros.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>

#include <basic/Tracer.hh>



#include <basic/database/open.hh>


#include <iostream>

#include <future>
#include <chrono>
#include <ctime>
//#include <cstdlib> // atexit
#include <cmath>

using std::string;


static basic::Tracer TR( "cloud_demo" );


char const * SCORE_FILE_DATA = R"_(SEQUENCE:
SCORE: total_score         rms CAPRI_rank       Fnat       I_sc       Irms   Irms_leg    cen_rms dslf_fa13    fa_atr    fa_dun   fa_elec fa_intra_rep fa_intra_sol_xover4              fa_rep              fa_sol hbond_bb_sc hbond_lr_bb    hbond_sc hbond_sr_bb interchain_contact interchain_env interchain_pair interchain_vdw lk_ball_wtd       omega     p_aa_pp pro_close rama_prepro         ref     st_rmsd yhh_planarity description
SCORE:    -295.196      35.124      0.000      0.000      0.000     16.936     17.236     13.195     0.000 -1054.125   275.583  -288.386        2.364              43.218             219.583             647.125     -35.910     -42.367     -14.041     -54.538            -19.000        -15.908          -2.429          0.348     -19.282       9.691     -43.599     3.157       8.266      48.065      13.322         0.000 input_0001_0053
SCORE:    -409.548      19.132      0.000      0.085     -6.809     10.954     10.739     20.539     0.000 -1063.954   257.366  -283.172        2.180              40.290             143.710             641.122     -35.799     -42.367     -16.870     -54.538             10.000        -17.034           0.076          0.035     -23.095       9.691     -43.599     3.157       8.266      48.065      19.084         0.000 input_0001_0067
SCORE:    -380.275      19.892      0.000      0.000     -5.052     11.304     11.034     19.142     0.000 -1064.496   256.821  -294.775        2.190              39.664             171.727             651.620     -35.277     -42.367     -15.757     -54.538              5.000        -16.821          -0.255          0.000     -23.935       9.691     -43.599     6.423       8.266      48.065      19.041         0.000 input_0001_0056
SCORE:    -377.828      25.371      0.000      0.000     -8.759     12.393     11.788     17.563     0.000 -1072.145   261.577  -294.009        2.234              42.053             172.069             654.136     -37.002     -42.367     -16.246     -54.538              1.000        -14.698          -0.255          0.576     -22.438       9.691     -43.599     6.423       8.266      48.065      20.576         0.000 input_0001_0007
SCORE:    -240.275      18.511      0.000      0.000     -3.358     11.433     10.966     12.065     0.000 -1053.480   267.412  -275.647        2.303              41.731             252.094             640.246     -34.889     -42.367     -13.439     -54.538            -18.000        -16.717          -2.098          0.374     -19.758       9.691     -43.599    27.633       8.266      48.065       9.874         0.000 input_0001_0016
SCORE:    -367.760      22.110      0.000      0.021     -4.157     12.245     12.069     18.657     0.000 -1057.905   267.309  -282.149        2.232              40.419             165.872             639.172     -36.127     -42.367     -11.366     -54.538            -20.000        -16.739           0.787          0.724     -23.892       9.691     -43.599     3.157       8.266      48.065      16.658         0.000 input_0001_0070
SCORE:    -406.923      10.155      0.000      0.021     -4.410      6.158      5.646     10.209     0.000 -1068.288   259.297  -291.625        2.198              39.624             145.988             651.636     -35.799     -42.367     -16.209     -54.538              6.000        -15.441          -1.953          0.000     -22.422       9.691     -43.599     3.157       8.266      48.065       9.752         0.000 input_0001_0002
SCORE:    -367.084       8.607      0.000      0.064    -10.199      4.775      4.894      9.681     0.000 -1077.589   261.431  -297.581        2.211              41.520             186.405             659.714     -37.114     -42.367     -17.190     -54.538              2.000        -13.533          -0.035          0.287     -20.833       9.691     -43.599     6.423       8.266      48.065       9.299         0.000 input_0001_0014
SCORE:    -407.970      16.268      0.000      0.043     -8.102      8.323      7.951     16.089     0.000 -1070.114   261.932  -289.366        2.201              39.955             143.038             650.110     -35.424     -42.367     -15.331     -54.538              9.000        -16.917          -0.236          0.000     -23.645       9.691     -43.599     3.157       8.266      48.065      14.879         0.000 input_0001_0071
SCORE:    -406.298      19.081      0.000      0.064     -6.536      7.011      6.717     11.812     0.000 -1072.724   260.530  -290.459        2.219              41.053             150.349             649.991     -37.550     -42.367     -15.866     -54.538              4.000        -14.094          -0.792          0.307     -22.517       9.691     -43.599     3.157       8.266      48.065       8.290         0.000 input_0001_0061
SCORE:    -390.326      14.968      0.000      0.021     -3.161      7.685      7.215     12.398     0.000 -1072.907   258.915  -283.450        2.214              41.180             160.610             646.843     -36.640     -42.367     -15.591     -54.538             -1.000        -13.461           0.482          0.078     -23.441       9.691     -43.599     6.423       8.266      48.065      10.657         0.000 input_0001_0024
SCORE:     321.039      17.820      0.000      0.000     -2.027      9.036      8.539     11.942     0.000 -1076.825   269.498  -284.126        2.237              44.736             855.161             651.670     -36.074     -42.367     -13.323     -54.538            -20.000        -15.138           3.064          0.027     -20.590       9.691     -43.599     3.157       8.266      48.065      10.566         0.000 input_0001_0062
SCORE:    -398.726      14.299      0.000      0.000     -7.702      6.628      6.129     12.244     0.000 -1066.552   256.871  -284.298        2.188              40.314             149.951             643.889     -35.401     -42.367     -14.598     -54.538             -4.000        -14.418          -2.060          0.361     -23.032       9.691     -43.599     6.423       8.266      48.065      13.139         0.000 input_0001_0018
SCORE:    -400.173      17.364      0.000      0.000     -8.643     10.136      9.692     16.323     0.000 -1075.814   259.928  -283.749        2.212              40.515             156.663             649.171     -37.589     -42.367     -15.395     -54.538              3.000        -16.661          -1.025          0.636     -24.792       9.691     -43.599     3.157       8.266      48.065      19.309         0.000 input_0001_0039
SCORE:    -406.368      24.592      0.000      0.043     -7.825     13.373     12.700     23.019     0.000 -1073.531   257.835  -285.455        2.181              40.037             150.145             650.092     -36.596     -42.367     -17.037     -54.538              9.000        -16.133          -2.027          0.000     -22.713       9.691     -43.599     3.157       8.266      48.065      16.453         0.000 input_0001_0005
SCORE:    -401.854      20.491      0.000      0.021     -4.640      9.872      9.762     22.576     0.000 -1077.420   257.985  -286.751        2.187              39.679             153.166             655.687     -35.998     -42.367     -14.592     -54.538              6.000        -16.302           0.027          0.102     -24.473       9.691     -43.599     3.157       8.266      48.065      18.702         0.000 input_0001_0011
SCORE:    -359.359      10.898      0.000      0.085     -5.751      5.352      4.980     11.764     0.000 -1069.573   256.848  -282.495        2.205              39.210             187.024             646.385     -35.600     -42.367     -15.643     -54.538             -3.000        -14.845          -2.411          0.005     -21.753       9.691     -43.599     8.514       8.266      48.065      12.821         0.000 input_0001_0004
SCORE:    -408.340      24.398      0.000      0.000     -8.388     13.184     13.151     18.902     0.000 -1072.901   260.786  -293.834        2.206              40.561             149.667             651.213     -36.243     -42.367     -15.479     -54.538              5.000        -16.475          -2.027          0.587     -24.072       9.691     -43.599     4.238       8.266      48.065      13.676         0.000 input_0001_0065
SCORE:    -405.850       8.811      0.000      0.021     -8.629      5.020      4.743     11.380     0.000 -1075.410   259.954  -285.754        2.219              41.197             149.307             648.798     -35.401     -42.367     -16.565     -54.538              5.000        -14.562          -1.202          0.110     -22.870       9.691     -43.599     3.157       8.266      48.065      12.108         0.000 input_0001_0001
SCORE:    -400.582      19.712      0.000      0.000     -7.304     12.044     11.458     18.460     0.000 -1067.256   256.327  -281.354        2.193              40.209             156.665             640.992     -38.590     -42.367     -14.508     -54.538              1.000        -15.562          -0.258          0.006     -23.933       9.691     -43.599     3.157       8.266      48.065      19.065         0.000 input_0001_0017
SCORE:    -374.680      20.642      0.000      0.021     -8.181     11.404     11.153     17.718     0.000 -1075.919   258.573  -286.221        2.258              41.837             172.019             652.910     -35.478     -42.367     -16.484     -54.538            -13.000        -16.492           0.196          0.782     -22.994       9.691     -43.599     9.300       8.266      48.065      18.815         0.000 input_0001_0008
SCORE:    -405.834      16.534      0.000      0.000    -13.037      9.461      9.438     15.722     0.000 -1086.198   260.191  -290.230        2.214              41.679             158.758             657.384     -39.092     -42.367     -14.592     -54.538              0.000        -16.489          -0.065          0.000     -24.623       9.691     -43.599     3.157       8.266      48.065      17.189         0.000 input_0001_0003
SCORE:     227.980      18.758      0.000      0.021    -13.756     11.197     10.854     15.911     0.000 -1072.626   264.789  -288.684        2.227              40.753             773.936             653.651     -34.121     -42.367     -16.455     -54.538            -15.000        -14.699          -1.351          0.325     -24.166       9.691     -43.599     3.157       8.266      48.065      14.809         0.000 input_0001_0057
SCORE:    -406.069      11.531      0.000      0.170    -12.140      6.268      5.627     13.163     0.000 -1075.674   255.181  -285.105        2.209              40.169             152.762             653.255     -34.831     -42.367     -17.026     -54.538             -7.000        -11.784          -0.745          0.162     -25.683       9.691     -43.599     3.157       8.266      48.065      15.052         0.000 input_0001_0074
SCORE:    -378.017      18.686      0.000      0.043    -11.631     11.361     10.576     18.613     0.000 -1074.661   265.468  -296.571        2.223              42.571             167.484             657.355     -36.769     -42.367     -15.482     -54.538            -16.000        -18.125           0.311          0.060     -21.578       9.691     -43.599     6.423       8.266      48.065      16.521         0.000 input_0001_0028
SCORE:    -404.523       8.902      0.000      0.085     -9.753      4.492      4.332      8.120     0.000 -1071.759   258.127  -289.668        2.207              40.835             155.084             648.313     -37.935     -42.367     -15.972     -54.538             -8.000        -15.767          -0.702          0.000     -22.430       9.691     -43.599     3.157       8.266      48.065       7.163         0.000 input_0001_0015
)_";


using namespace protocols::network;

int main(int argc, char * argv [])
{
	try {
		devel::init(argc, argv);


		//core::import_pose::pose_from_file("../../test/core/io/1QYS.pdb");

		core::pose::Pose pose;
		//core::import_pose::pose_from_file(pose, "../../src/python/PyRosetta/src/test/data/test_in.pdb", core::import_pose::PDB_file);

		//auto pdb_file_name = "/Users/sergey/work/main.ui/source/build/build-qt-debug/0-test.pdb";
		string pdb_file_name;

		auto input_files = basic::options::option[ basic::options::OptionKeys::in::file::s ]();
		if( input_files.empty() ) pdb_file_name = basic::database::full_name("../source/test/core/io/test_in.pdb");
		else pdb_file_name = input_files.front();

		TR << "pdb_file_name:" << pdb_file_name << std::endl;

		core::import_pose::pose_from_file(pose, pdb_file_name, core::import_pose::PDB_file);



		TR << "--------------------------------------------------------------------------------" << std::endl;

		//TR << SCORE_FILE_DATA << std::endl;
		//post_file("docking.score", SCORE_FILE_DATA);
		auto score_lines = utility::split_by_newlines(SCORE_FILE_DATA);


		// auto es_id = NetworkQueue::get_execution_summary_id();
		// TR << "ES_id:" << es_id << std::endl;


		if(true) {
			std::clock_t c_start = std::clock();
			auto t_start = std::chrono::high_resolution_clock::now();

			int old_score_lines_count = 0;

			// const int iterations = 100;
			// for(int i=0; i<iterations; ++i) {

			int i=0;
			while(true) {
				++i;

				int r =  (i / 16) % (pose.size() / 4) + pose.size() / 2;

				pose.set_phi(r, 16.0 * i);
				post_decoy("decoy.pdb", pose);


				int const score_update_period = 24;
				int score_lines_count = score_lines.empty() ? 0 : i/score_update_period % ( score_lines.size() + 1 );
				if( score_lines_count < 2  and  score_lines.size() >= 2 ) score_lines_count = 2;

				if( old_score_lines_count != score_lines_count ) {
					//TR << "SL: " << score_lines_count << std::endl;
					auto sf = utility::join( std::vector<string>(score_lines.begin(), std::next( score_lines.begin(), score_lines_count) ), "\n");

					post_file("docking.score", sf);
					old_score_lines_count = score_lines_count;
				}


				int log_time = 1 + i / 64;
				int n_logs = int( log2(log_time) );
				int log_file_index = log_time % (n_logs + 1);
				post_file("log." + std::to_string(log_file_index), "dummy-" + std::to_string(i) + '\n', true);

				// for(int k = 0; k < log2(i/64); ++k) {
				// 	if( i % int( pow(k+1, 2) ) == 0 ) {

				std::this_thread::sleep_for(std::chrono::milliseconds(100));
			}

			std::clock_t c_end = std::clock();
			auto t_end = std::chrono::high_resolution_clock::now();

			TR << std::fixed << std::setprecision(2) << "CPU time used: "
					  << 1000.0 * (c_end-c_start) / CLOCKS_PER_SEC << " ms\n"
					  << "Wall clock time passed: "
					  << std::chrono::duration<double, std::milli>(t_end-t_start).count()
					  << " ms\n" << std::endl;
		}
		//test();

		//std::cerr << "This app need to be build with extras=serialization! Aborting..." << std::endl;

		//for(auto &r : RESULTS) r.wait();

		std::cerr << "main(): done!" << std::endl;
		return 1;
	} catch (utility::excn::Exception const & e ) {
		TR << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}


// attic
/*

void post_decoy_(core::pose::PoseUP &&pose, std::shared_ptr<AsyncResult> r) // yes, we hold our own future!
{
	// Creating message
	std::ostringstream os;
	pose->dump_pdb(os);

	// Compressing message
	std::ostringstream zmsg;
	zlib_stream::zip_ostream zipper(zmsg, true);
	zipper << os.str();
	zipper.zflush_finalize();

	auto auth = basic::options::option[ basic::options::OptionKeys::auth ]();
	string auth_encoded = "Basic " + utility::io::base64_encode(reinterpret_cast<unsigned char const *>( auth.c_str() ), auth.size());

	httplib::Client cli("localhost", 64082);
	httplib::Headers headers { {"Authorization", auth_encoded}, {"Content-Encoding", "gzip"}, {"User-Agent", "Rosetta/3.7"} };
	auto res = cli.Post("/api.execution/file", headers, zmsg.str(), "application/octet-stream"); // "text/plain"

	std::shared_ptr<AsyncResult> r1 = r;
}

static std::vector< std::future<void> > RESULTS;

void post_decoy_async(core::pose::Pose const &pose)
{
	core::pose::PoseUP p( new core::pose::Pose() );
	p->detached_copy(pose);

	//post_decoy( std::move(p) );

	//results.rs.push_back( std::async(std::launch::async, post_decoy_, std::move(p)) );
	//RESULTS.push_back( std::async(std::launch::async, post_decoy_, std::move(p)) );

	std::shared_ptr<AsyncResult> r = std::make_shared<AsyncResult>();
	r->f = std::async(std::launch::async, post_decoy_, std::move(p), r); // yes, it is magic!
}


void test()
{
	// string user, password;
	// std::tie(user, password) = user_password( basic::options::option[ basic::options::OptionKeys::auth ]() );
	// TR << "user: " << user << "  password: " << password << std::endl;
	auto auth = basic::options::option[ basic::options::OptionKeys::auth ]();
	string auth_encoded = "Basic " + utility::io::base64_encode(reinterpret_cast<unsigned char const *>( auth.c_str() ), auth.size());
	//TR << auth_encoded;

	httplib::Client cli("localhost", 64082);

	{
		//httplib::Headers headers { {"Authorization", auth_encoded} };
		//headers.emplace("Authorization", auth_encoded);

		// auto res = cli.Get("/api/file/379/output/0.docking-pre-pack.score", headers);
		// if(res) {
		// 	TR << "status: " << res->status << std::endl;
		//     TR << res->body << std::endl;
		// }
	}

	{
		httplib::Headers headers { {"Authorization", auth_encoded}, {"Content-Encoding", "gzip"}, {"User-Agent", "Rosetta/3.7"} };

		std::ifstream input("../../test/core/io/1QYS.pdb");

		string payload = (std::stringstream() << input.rdbuf()).str();
		TR << payload << std::endl << "size: " << payload.size() << std::endl;

		//string payload(25600, '#');
		httplib::detail::compress(payload);
		TR << "Compressed size: " << payload.size() << std::endl;


		auto res = cli.Post("/_test", headers, payload, "text/plain");
		if(res) {
			TR << "status: " << res->status << std::endl;
			TR << res->body << std::endl;
		}
	}
}
*/
