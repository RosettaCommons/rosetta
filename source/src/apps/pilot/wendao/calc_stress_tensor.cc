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

//for tensor
#include <utility>
#include <core/scoring/methods/EnergyMethodOptions.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/AtomVDW.hh>
#include <core/chemical/ResidueType.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <core/types.hh>
#include <numeric/xyz.functions.hh>
#include <core/scoring/MinimizationGraph.hh>
#include <core/optimization/CartesianMinimizerMap.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

using namespace std;
using namespace core::optimization;
using namespace core::scoring;

class Tensor : public utility::pointer::ReferenceCount
{
public:
    Tensor(): tensor_(0), count(0)
    {
        //init
    }

    void check()
    {
        Real p=0;
        p = tensor_.xx()+tensor_.yy()+tensor_.zz();
        cout << count << " " << -1000.0*p/count << endl;
    }

    void sum_to( Real &P, Size &N )
    {
        P += tensor_.xx()+tensor_.yy()+tensor_.zz();
        N += count;
    }

    void add( const Vector &f, const Vector &x, Real rsq )
    {
        Real r = sqrt(rsq);
        numeric::xyzMatrix<Real> tens = numeric::outer_product(f, x);
        tensor_ += tens / (r*r*r); //4/3 pi r^3
        count++;
    }

    numeric::xyzMatrix<Real> tensor_;
    int count;
};

// Map a tensor instance to each pair, global
std::map< std::string, Tensor > storage_tensor;

///////////////////////////////////////////////////////////////////////////////
// Only calculate pair interaction to make sure pair and vdw get balanced
// No one body term considered
class CustomMover : public protocols::moves::Mover
{
public:
    CustomMover():
    atom_vdw_( ScoringManager::get_instance()->get_AtomVDW("centroid_rot") )
    {

    }

    void apply( core::pose::Pose &pose)
    {
        using namespace core;
        using namespace core::scoring;
        using namespace core::optimization;
        using namespace core::optimization::symmetry;

        core::scoring::ScoreFunctionOP score_function = core::scoring::getScoreFunction();
        protocols::simple_moves::SwitchResidueTypeSetMover("centroid_rot").apply(pose);
        AtomVDW const & score_vdw( ScoringManager::get_instance()->get_AtomVDW(
            score_function->energy_method_options().atom_vdw_atom_type_set_name() ) );

        //no need to repack and min(?), just score

        // repack pose
        // core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
        // main_task_factory->push_back( new core::pack::task::operation::RestrictToRepacking );
        // protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;
        // pack_mover->task_factory( main_task_factory );
        // pack_mover->score_function( score_function );

        // sc min pose
        // core::optimization::AtomTreeMinimizer minimizer;
        // core::optimization::MinimizerOptions options( "lbfgs_armijo_nonmonotone", 0.01, true, false, false );
        // options.max_iter(25);
        core::kinematics::MoveMap mm;
        mm.set_bb ( false );
        mm.set_chi ( true );
        // minimizer.run( pose, mm, *score_function, options );

        // rescore pose
        (*score_function)(pose);
        //protocols::jd2::JobOP job( protocols::jd2::JobDistributor::get_instance()->current_job() );

        //setup cart ??
        CartesianMinimizerMap min_map;
        min_map.setup( pose, mm );
        score_function->setup_for_minimizing( pose, min_map );

        assert( pose.energies().minimization_graph() );

        MinimizationGraphCOP mingraph = pose.energies().minimization_graph();

        // how to get pair force ??
        // copy from cartesian_minimize

        /// 2. eval inter-residue derivatives
        for ( utility::graph::Node::EdgeListConstIter edgeit = mingraph->const_edge_list_begin(),
                edgeit_end = mingraph->const_edge_list_end();
                edgeit != edgeit_end; ++edgeit )
        {

            MinimizationEdge const &minedge = static_cast< MinimizationEdge const & > ( (**edgeit) );
            Size const rsd1ind = minedge.get_first_node_ind();
            Size const rsd2ind = minedge.get_second_node_ind();

            conformation::Residue const &rsd1( pose.residue( rsd1ind ));
            conformation::Residue const &rsd2( pose.residue( rsd2ind ));
            ResSingleMinimizationData const &r1_min_data( mingraph->get_minimization_node( rsd1ind )->res_min_data() );
            ResSingleMinimizationData const &r2_min_data( mingraph->get_minimization_node( rsd2ind )->res_min_data() );

            //new interface
            utility::vector1< DerivVectorPair > r1atom_derivs(min_map.atom_derivatives( rsd1ind ));
            utility::vector1< DerivVectorPair > r2atom_derivs(min_map.atom_derivatives( rsd2ind ));
            eval_weighted_atom_derivatives_for_minedge( minedge, rsd1, rsd2,
                    r1_min_data, r2_min_data, pose, score_function->weights(),
                    r1atom_derivs, r2atom_derivs);

            std::string key(rsd1.name3() + "_" + rsd2.name3());
            if (rsd1.aa() > rsd2.aa()) key = rsd2.name3() + "_" + rsd1.name3();

            //cout << key << endl;
            //cout << r1atom_derivs[rsd1.nbr_atom()+1].f2().x() << endl;
            //cout << r2atom_derivs[rsd2.nbr_atom()+1].f2().x() << endl;

            Size const i(rsd1.nbr_atom()+1);
            Size const j(rsd2.nbr_atom()+1);

            //cutoff
            Vector dr = rsd1.xyz(i)-rsd2.xyz(j);
            //if (dr.length_squared()>25.0) continue; //only consider small radii

            Size const i_type( rsd1.atom_type_index(i) );
            Size const j_type( rsd2.atom_type_index(j) );

            //new interface force on CEN_i
            core::Vector f12(r1atom_derivs[i].f2());

            //cal vdw manually, because it's not pairwisealbe
            Real bump_sq = score_vdw(i_type)[j_type];
            Real dis2( dr.length_squared() );
            if ( dis2 < bump_sq ) {
                Real dE_dr_over_r = score_function->weights()[vdw] * 4.0 * (dis2-bump_sq) /bump_sq;
                f12 += dE_dr_over_r * dr;
            }

            //size of CEN_i
            Real r_sq = atom_vdw_(i_type)[j_type];
            //cout << key << ": " << r_sq << ", " << f12.length_squared() << endl;

            if (storage_tensor.find(key) == storage_tensor.end())
            {
                //cout << "New!" << endl;
                Tensor newtensor;
                storage_tensor[key] = newtensor;
            }

            storage_tensor[key].add(f12, dr, r_sq);
            //storage_tensor[key].check();
        }
    }

    virtual std::string get_name() const
    {
        return "CustomMover";
    }

private:
    AtomVDW const & atom_vdw_;
};

///////////////////////////////////////////////////////////////////////////////
void stat_output()
{
    Real P=0.0;
    Size N=0;

    for (std::map< std::string, Tensor >::iterator it=storage_tensor.begin(),
        it_end = storage_tensor.end(); it!=it_end; it++)
    {
        cout << it->first << " ";
        it->second.check();
        it->second.sum_to(P, N);
    }

    cout << "#ALL " << N << " " << -1000.0*P/N << endl;
}


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
    catch (utility::excn::Exception &excn )
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

        stat_output();
    }
    catch (utility::excn::Exception const &e )
    {
        std::cout << "caught exception " << e.msg() << std::endl;
    }
    return 0;
}
