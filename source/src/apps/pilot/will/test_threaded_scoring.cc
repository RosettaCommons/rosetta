#include <protocols/init/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <vector>
#include <iostream>

#ifdef USE_OPENMP
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
    int nthreads = 1;
    #ifdef USE_OPENMP
    nthreads = omp_get_max_threads();
    #endif

    std::cout << "lkball_thread_bug test with " << nthreads << " threads" << std::endl;

    std::mt19937 rng((unsigned int)time(0) + 192634);
    std::uniform_real_distribution<> runif(-1.0, 1.0);

    protocols::init::init(argc,argv);

    std::cout << "create separate sfxn for each thread to use" << std::endl;
    std::vector<core::scoring::ScoreFunctionOP> sfxns(nthreads,nullptr);
    for( auto & sfxnop : sfxns ){
        sfxnop = core::scoring::get_score_function();
    }
    std::cout << "created " << sfxns.size() << " sfxns" << std::endl;

    std::cout << "read in poses from -in:file:s" << std::endl;
    utility::vector1<core::pose::PoseOP> poses_in_vector1;
    core::import_pose::read_all_poses(
                        basic::options::option[basic::options::OptionKeys::in::file::s](), poses_in_vector1 );
    runtime_assert_msg(poses_in_vector1.size() > 1, "-in:file:s at least two structures");
    std::cout << "poses_in_vector1.size() " << poses_in_vector1.size() << std::endl;

    int N = 1024;

    std::cout << "preclone " << N << " poses for full list of operations" << std::endl;
    std::vector<core::pose::PoseOP> poses(N,nullptr);
    for(int i = 0; i < N; ++i){
        // std::cout << "making pose clone " << i << std::endl;
        poses[i] = poses_in_vector1[i%(poses_in_vector1.size()-1)+2]->clone();
        int ir = poses[i]->size()/2 + 1;
        // perturb middle phi by +-1 degree
        poses[i]->set_phi(ir, poses[i]->phi(ir) + runif(rng) );
    }


    // // pre-scoring each input pose single-threaded fixes the problem...
    // // this should not effect the actual poses getting scored (already cloned above)
    // // uncomment these lines and things work
    // std::cout << "single threaded score to warm up whatever database files will be needed"
    //           << "without thread issues" << std::endl;
    // for(auto poseop : poses_in_vector1){
    //     sfxns.front()->score(*poseop);
    // }

    // // what about scoring only one... also seems to fix it
    // sfxns.front()->score(*poses.front());

    std::cout << "parallel loop over poses, score each with thread-local sfxn" << std::endl;
    std::vector<double> scores(N,0);
    #ifdef USE_OPENMP
    #pragma omp parallel for schedule(dynamic,1)
    #endif
    for(int i = 0; i < N; ++i){
        // #ifdef USE_OPENMP
        // #pragma omp critical
        // #endif
        // std::cout << "i: " << i << ", thread " << omp_get_thread_num() << std::endl;
        int threadnum = 0;
        #ifdef USE_OPENMP
        threadnum = omp_get_thread_num();
        #endif
        core::scoring::ScoreFunctionOP this_thread_sfxn = sfxns[threadnum];
        core::pose::PoseOP poseop_to_score = poses[i];
        scores[i] = this_thread_sfxn->score(*poseop_to_score);
    }

    std::cout << "scores are:" << std::endl;
    for(int i = 0; i < N; ++i){
        std::cout << "score " << i << " " << scores[i] << std::endl;
    }

    return 0;
}
