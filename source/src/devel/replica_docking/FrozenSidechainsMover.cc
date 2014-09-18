

#include <devel/replica_docking/FrozenSidechainsMover.hh>
#include <devel/replica_docking/FrozenSidechainsMoverCreator.hh>

#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/import_pose/import_pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/options/option.hh>

#include <protocols/scoring/InterfaceInfo.hh>
#include <protocols/scoring/Interface.hh>
#include <protocols/simple_moves/ReturnSidechainMover.hh>
#include <protocols/jd2/util.hh>

#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

#include <utility/io/ozstream.hh>

#include <basic/Tracer.hh>

#include <numeric/MathVector.hh>

#include <ObjexxFCL/format.hh>

#include <fstream>

using namespace utility::tag;

namespace devel {
namespace replica_docking {

using namespace core;
using namespace protocols;
using namespace core::scoring;
using namespace protocols::scoring;

static thread_local basic::Tracer tr( "devel.replica_docking.FrozenSidechainsMover" );

std::string
FrozenSidechainsMoverCreator::keyname() const
{
  return FrozenSidechainsMoverCreator::mover_name();
}

moves::MoverOP
FrozenSidechainsMoverCreator::create_mover() const {
  return new FrozenSidechainsMover;
}

std::string
FrozenSidechainsMoverCreator::mover_name()
{
  return "FrozenSidechainsMover";
}

///////////////////////////////////////////////////////////////////////

FrozenSidechainsMover::FrozenSidechainsMover()
	: moves::Mover( FrozenSidechainsMoverCreator::mover_name() )
{
}

moves::MoverOP
FrozenSidechainsMover::clone() const {
  return new FrozenSidechainsMover( *this );
}

//FrozenSidechainsMover::~FrozenSidechainsMover(){}

  static bool first( true );
  static utility::vector1<bool> frozen;
  static int ct_dump( 0 );

void
FrozenSidechainsMover::apply(pose::Pose& pose){

  //  allow_copy_chi_( pose.total_residue(), false );
  compute_interface( pose );  // calculate interface
  protocols::scoring::InterfaceInfo const & interfaceinfo( interface_from_pose( pose ) ); // get the interface information
  protocols::scoring::InterfaceCOP interface = interfaceinfo.interface( 1 );   // only work for 1 jump

  // if a residue does not locate in the interface, then move its sidechains back to Crystal state
//   for ( core::Size i = 1; i <= pose.total_residue(); i++ ) {
//     allow_copy_chi_[ i ] = !interface.is_interface( pose.residue( i ) );
//     if ( allow_copy_chi_[i] ) {
//       tr.Debug << i << std::endl;
//     }
//   }

//   if ( recover_sidechains_pose_.is_fullatom() ) {
//     frozen_sidechains_ = new protocols::simple_moves::ReturnSidechainMover( recover_sidechains_pose_, allow_copy_chi_ );
//     tr.Debug << "got input pose, and generated ReturnSidechainMover" << std::endl;
//     frozen_sidechains_->apply( pose );
//   }
  ++ct_dump;
  bool dump( ct_dump % 10000 == 0 );
  if ( first ) {
    // recover_sidechains_pose_ = pose;
    //    recover_sidechains_pose_.dump_pdb(ObjexxFCL::lead_zero_string_of( jd2::current_replica(), 6 ) +"_used_for_freezing.pdb");
    for ( core::Size pos = 1; pos<=pose.total_residue(); pos++ ) {
      frozen.push_back( false );
    }
    first = false;
  }
  for ( core::Size pos = 1; pos <= pose.total_residue(); pos++ ) {
    if ( !interface->is_interface( pos ) ){
      //      tr.Debug << "residue " << pos << " is not interface residue, now move it back to crystal state..." << std::endl;
      if ( dump ) {
	frozen[pos]=true;
      }
      // this was slow
      //      pose.replace_residue( pos, recover_sidechains_pose_.residue( pos ), true );
      //      tr.Debug << "replace_residue used to move sidechains of non-interface residue back to X-ray state" << std::endl;

      //      tr.Debug << "set_chi used to move sidechains of non-interface residue back to X-ray state" << std::endl;
      core::Size nchi = pose.residue_type( pos ).nchi();
      for ( core::Size chi = 1; chi <= nchi; ++chi ) {
	pose.set_chi( chi, pos, recover_sidechains_pose_.chi( chi, pos ) );
      }
    } else {
      //tr.Debug << "residue " << pos << " is interface" << std::endl;
      if ( dump ) {
	frozen[pos]=false;
      }

    }
  }
  if ( tr.Trace.visible() && dump ) {
    std::string numbers=ObjexxFCL::lead_zero_string_of( jd2::current_replica(), 6 ) +"_"+ ObjexxFCL::lead_zero_string_of( ct_dump, 6 );
    std::ofstream dump_pos;
    std::string name("frozen_"+numbers+".txt");
    dump_pos.open( name.c_str() );
    for ( core::Size pos = 1; pos <= pose.total_residue(); pos++ ) {
      if ( ! frozen[pos] ) {
	dump_pos << pos << " ";
      }
    }
    dump_pos << std::endl;
    using namespace ObjexxFCL;
    std::string name2("dump_chis_"+numbers+".txt");
    std::ofstream dumpchi;
    dumpchi.open( name2.c_str() );
    for ( core::Size pos = 1; pos <= pose.total_residue(); pos++ ) {
      dumpchi << format::I( 3, pos) << " " << ( frozen[pos] ? "FROZEN " : "ALLOWED" );
      for ( core::Size chino = 1; chino <= 4; chino++ ) {
	dumpchi << " " << format::RJ( 10, chino <=pose.residue_type(pos).nchi() ? pose.chi( chino, pos ) : 0.0 );
      }
      dumpchi << " | ";
      for ( core::Size chino = 1; chino <= 4; chino++ ) {
	dumpchi << " " << format::RJ( 10, chino <=recover_sidechains_pose_.residue_type(pos).nchi() ? recover_sidechains_pose_.chi( chino, pos ) : 0.0 );
      }
      dumpchi << std::endl;
    }
    pose.dump_pdb( "frozen_"+numbers+".pdb" );
  }
}

std::string
FrozenSidechainsMover::get_name() const {
  return FrozenSidechainsMoverCreator::mover_name();
}

void
FrozenSidechainsMover::parse_my_tag(
  utility::tag::TagCOP tag,
  basic::datacache::DataMap &,
  filters::Filters_map const &,
  moves::Movers_map const &,
  pose::Pose const& pose
) {
  if ( tag->hasOption("recover_sidechains_pdb")) {
    std::string pdbfile = tag->getOption< std::string > ( "recover_sidechains_pdb" );
    core::import_pose::pose_from_pdb( recover_sidechains_pose_, pdbfile );
    tr.Debug << "sidechains of non-interface residues will be frozen back to the sidechains of the given recover_sidechains pdb file" << std::endl;
  } else {
    recover_sidechains_pose_ = pose;
    tr.Debug << "sidechains of the non-interface residues will be frozen back to the sidechain of input pose" << std::endl;
  }
}

void
FrozenSidechainsMover::compute_interface(
	core::pose::Pose & pose
	) const
{
  InterfaceInfo & interface( nonconst_interface_from_pose( pose ) );

  /// initialize the cenlist info:
  /// only if they have not been calculated since the last score
  if ( !interface.calculated() ) {

    // initialize values
    interface.initialize();

    // compute interpolated number of neighbors at various distance cutoffs
    interface.calculate( pose );
  }

  interface.calculated() = true;
}

/// @details Pose must already contain a Interface object or this method will fail
InterfaceInfo const &
FrozenSidechainsMover::interface_from_pose( core::pose::Pose const & pose ) const
{
  //using core::pose::datacache::CacheableDataType::INTERFACE_INFO;
  return *( static_cast< InterfaceInfo const * >(pose.data().get_const_ptr( core::pose::datacache::CacheableDataType::INTERFACE_INFO )() ) );
}

/// @details Either returns a non-const reference to the Interface object that already exists
/// in the pose, or creates a new Interface object, places it in the pose, and then returns
/// a non-const reference to it
InterfaceInfo &
FrozenSidechainsMover::nonconst_interface_from_pose( core::pose::Pose & pose ) const
{
  //using core::pose::datacache::CacheableDataType::INTERFACE_INFO;

  if ( pose.data().has( core::pose::datacache::CacheableDataType::INTERFACE_INFO ) ) {
    return *( static_cast< InterfaceInfo * >(pose.data().get_ptr( core::pose::datacache::CacheableDataType::INTERFACE_INFO )() ) );
  }
  // else
  InterfaceInfoOP interface = new InterfaceInfo;
  pose.data().set( core::pose::datacache::CacheableDataType::INTERFACE_INFO, interface );
  return *interface;
}

}//replica_docking
}//devel
