
#include <ui/task/functor.h>

#include <QDebug>

namespace ui {
namespace task {


void Functor::state(State state)
{
	state_ = state;
	if(state == State::_running_) {
		Q_EMIT started();
		Q_EMIT tick();
	}
	if(state == State::_finished_) {
		Q_EMIT finished();
		Q_EMIT final();
	}
}

void Functor::execute()
{
	qDebug() << "Executing Functor<" << name() << ">...";
	state(State::_running_);
	run();
}

bool Functor::is_running() const
{
	return state_ == State::_running_;
}

bool Functor::is_finished() const
{
	return state_ == State::_finished_;
}

// return <current_value, max_value> for execution progress
std::pair<int, int> Functor::progress() const
{
	// if( is_finished() ) return std::make_pair(1, 1);
	// else return std::make_pair(0, 1);
	if( is_running() ) return std::make_pair(0, 1);
	else return std::make_pair(0, 0);
}


void Functor::run()
{
	state(State::_finished_);
}

void Functor::finish()
{
	state(State::_finished_);
}



void FunctorSequence::add_functor(ValueType const &v)
{
	sequence_.push_back(v);
}

void FunctorSequence::run()
{
	if( !sequence_.empty() ) {
		qDebug() << "FunctorSequence<" << name() << ">::run()...";

		for(uint i=0; i<sequence_.size(); ++i) {
			connect(sequence_[i].get(), &Functor::started, this, &FunctorSequence::next);
			if( i < sequence_.size() - 1 ) {
				connect(sequence_[i].get(), &Functor::tick,  this, &Functor::tick);
				connect(sequence_[i].get(), &Functor::final, this, &Functor::tick);
				connect(sequence_[i].get(), &Functor::final, sequence_[i+1].get(), &Functor::execute);
			}
		}
		connect(sequence_.back().get(), &Functor::final, this, &Functor::finish);

		sequence_.front()->execute();
		// connect(this, &FunctorSequence::final,
		// 		[]{
		// 		qDebug() << "test";
		// 		});
	}
	else finish();
}

// return <current_value, max_value> for execution progress
std::pair<int, int> FunctorSequence::progress() const
{
	if( is_running() ) {
		std::pair<int, int> t {0, sequence_.size()};
		for(auto const & f : sequence_) {
			auto r = f->progress();
			t.first  += r.first;
			t.second += r.second;

			// qDebug() << "FunctorSequence::progress(): checking: f->is_running() for " << f->name() << " " << f->is_running();
			// if( f->is_running() ) {
			// 	auto r = f->progress();
			// 	t.first  += r.first;
			// 	t.second += r.second;
			// }
			// else t.first += 1;
		}
		//qDebug() << "FunctorSequence::progress(): " << t;
		return t;
	}
	else if( is_finished() ) return std::make_pair(0, 0);
	else return std::make_pair(0, 0);
}



FunctorNetworkCall::FunctorNetworkCall(QString const &name, CallCallback const &, QObject *parent) : Functor(name, parent)
{
}

void FunctorNetworkCall::run()
{
	connect(&request_, &NetworkCall::finished, this, &Functor::finish);
	if(callback_) callback_(request_);

	//request_.call("http://ui.graylab.jhu.edu/api/task");

	//auto nc = new NetworkCall;
	//nc->call("http://ui.graylab.jhu.edu/api/task");
}


} // namespace task
} // namespace ui
