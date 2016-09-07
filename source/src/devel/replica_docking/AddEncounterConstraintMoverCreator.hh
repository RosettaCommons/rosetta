/// @brief
/// @author Zhe Zhang
#ifndef INCLUDED_devel_replica_docking_AddEncounterConstraintMoverCreator_hh
#define INCLUDED_devel_replica_docking_AddEncounterConstraintMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace devel {
namespace replica_docking {

class AddEncounterConstraintMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	static std::string mover_name();

};

}
}
#endif
