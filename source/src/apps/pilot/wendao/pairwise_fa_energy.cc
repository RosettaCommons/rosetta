/// @file
/// @brief


#include <devel/init.hh>

#include <protocols/viewer/viewers.hh>
#include <protocols/moves/Mover.hh>

#include <protocols/simple_moves/PackRotamersMover.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/symmetry/SymmetricScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/PDBInfo.hh>
#include <core/io/pdb/Remarks.hh>
#include <core/id/AtomID.hh>

#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/operation/TaskOperations.hh>
#include <core/conformation/util.hh>


#include <protocols/idealize/IdealizeMover.hh>
#include <protocols/rbsegment_relax/util.hh>

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/types.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/AtomTreeMultifunc.hh>

#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/optimization/symmetry/SymMinimizerMap.hh>
#include <core/optimization/symmetry/SymAtomTreeMultifunc.hh>

#include <core/optimization/CartesianMinimizerMap.hh>
#include <core/optimization/CartesianMultifunc.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/scoring/EnergyGraph.hh>

#include <utility/excn/Exceptions.hh>
#include <protocols/jd2/Job.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <core/kinematics/MoveMap.hh>


#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>


///////////////////////////////////////////////////////////////////////////////

class CustomMover : public protocols::moves::Mover
{
public:
    CustomMover() {}
    void apply( core::pose::Pose &pose)
    {
        using namespace core;
        using namespace core::scoring;
        using namespace core::optimization;
        using namespace core::optimization::symmetry;

        core::scoring::ScoreFunctionOP score_function = core::scoring::get_score_function();

        // repack pose
        core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
        main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
        protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
        pack_mover->task_factory( main_task_factory );
        pack_mover->score_function( score_function );

        // sc min pose
        core::optimization::AtomTreeMinimizer minimizer;
        core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
        options.max_iter(25);
        core::kinematics::MoveMap mm;
        mm.set_bb ( false );
        mm.set_chi ( true );
        minimizer.run( pose, mm, *score_function, options );

        // rescore pose
        (*score_function)(pose);
        protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

        // loop over energy graph
        EnergyGraph const &energy_graph( pose.energies().energy_graph() );
        for (Size i = 1; i <= pose.size(); ++i)
        {
            core::conformation::Residue const &rsd_i( pose.residue(i) );

            for ( utility::graph::Graph::EdgeListConstIter
                    iru  = energy_graph.get_node(i)->const_upper_edge_list_begin(),
                    irue = energy_graph.get_node(i)->const_upper_edge_list_end();
                    iru != irue; ++iru )
            {
                EnergyEdge const *edge( static_cast< EnergyEdge const *> (*iru) );
                Size const e1( edge->get_first_node_ind() );
                Size const e2( edge->get_second_node_ind() );

                Real Erep = (*edge)[core::scoring::fa_rep];
                if (Erep > 10.0)
                {
                    if (pose.residue(e2).aa() > pose.residue(e1).aa())
                    {
                        std::cout << pose.residue(e1).name3() << " "
                                  << pose.residue(e2).name3() << " "
                                  << Erep << " " << job->input_tag() << std::endl;
                    }
                    else
                    {
                        std::cout << pose.residue(e2).name3() << " "
                                  << pose.residue(e1).name3() << " "
                                  << Erep << " " << job->input_tag() << std::endl;
                    }
                }
            }
        }
    }

    virtual std::string get_name() const
    {
        return "CustomMover";
    }
};

///////////////////////////////////////////////////////////////////////////////

void *
my_main( void *)
{
    using namespace protocols::moves;
    using namespace protocols::simple_moves::symmetry;

    SequenceMoverOP seq( new SequenceMover() );
    seq->add_mover( new CustomMover() );

    try
    {
        protocols::jd2::JobDistributor::get_instance()->go( seq );
    }
    catch ( utility::excn::EXCN_Base &excn )
    {
        std::cerr << "Exception: " << std::endl;
        excn.show( std::cerr );
    }

    return 0;
}

///////////////////////////////////////////////////////////////////////////////

int
main( int argc, char *argv [] )
{
    try
    {
        // initialize option and random number system
        devel::init( argc, argv );
        protocols::viewer::viewer_main( my_main );
    }
    catch ( utility::excn::EXCN_Base const &e )
    {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}
