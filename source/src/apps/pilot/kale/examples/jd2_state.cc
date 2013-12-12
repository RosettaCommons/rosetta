#include <core/types.hh>
#include <core/pose/Pose.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <devel/init.hh>
#include <basic/Tracer.hh>

using namespace std;

class StatefulMover : public protocols::moves::Mover {

public:

	StatefulMover() {
		cout << "  StatefulMover()" << endl;
		state_ = 0;
	}

	StatefulMover(core::Size value) {
		cout << "  StatefulMover(value)" << endl;
		cout << state_ << endl;
		state_ = value;
	}

	StatefulMover(StatefulMover const & other) {
		cout << "  StatefulMover(other)" << endl;
		cout << state_ << endl;
		state_ = other.state_;
	}

	StatefulMover& operator = (StatefulMover const & other) {
		cout << "  operator = (other)" << endl;
		cout << state_ << endl;
		state_ = other.state_;
	}

	void apply(core::pose::Pose & pose) {
		cout << "  apply(pose): "   << state_ << endl;
	}

	string get_name() const {
		return "  StatefulMover";
	}

	protocols::moves::MoverOP fresh_instance() const { 
		cout << "  fresh_instance()" << endl;
		return new StatefulMover;
	}

	bool reinitialize_for_each_job() const {
		cout << "  reinitialize_for_each_job()" << endl;
		return false;
	}

	bool reinitialize_for_new_input() const {
		cout << "  reinitialize_for_new_input()" << endl;
		return false;
	}

private:
	core::Size state_;
};

int main(int argc, char* argv[])
{
	devel::init(argc, argv);

	protocols::moves::MoverOP mover = new StatefulMover(5);

	protocols::jd2::JobDistributor::get_instance()->go(mover);

	return 0;
}
