/// @authoer Zhe Zhang
#ifndef INCLUDED_protocols_replica_docking_ModulatedMoverCreator_hh
#define INCLUDED_protocols_replica_docking_ModulatedMoverCreator_hh

#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace replica_docking {

class ModulatedMoverCreator : public protocols::moves::MoverCreator {
public:
	protocols::moves::MoverOP create_mover() const override;
	std::string keyname() const override;
	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

}
}
#endif  // INCLUDED_protocols_replica_docking_ModulatedMoverCreator_hh
