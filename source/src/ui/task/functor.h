#pragma once

#include <ui/task/functor.fwd.h>

#include <ui/task/util.h>

#include <QObject>

#include <deque>

namespace ui {
namespace task {



///
/// Something that could be executed and emit `done` on completion
///
class Functor : public QObject
{
    Q_OBJECT
public:
	enum class State {
		_none_,     /// have not run yet
		_running_,  /// still running
		_finished_, /// finished
	};

	using QObject::QObject;

	Functor(QString const &name, QObject *parent = Q_NULLPTR) : QObject(parent), name_(name) {}


	/// return current state
	State state() const { return state_; }

	/// return true if Functor is now executing
	bool is_running() const;

	/// return true if Functor is finished
	bool is_finished() const;

	QString name() const { return name_; }

	void execute();

	// return <current_value, max_value> for execution progress
	virtual std::pair<int, int> progress() const;


public Q_SLOTS:
	/// change state to _finished_
	void finish();

Q_SIGNALS:
	/// indicate start of this Functor, emitted only once
	void started();

	/// Might be emitted periodically if Functor have loops. At least one tick() event is guaranteed to happened between start and finish.
	void tick();

	/// indicate end of this Functor execution, emitted only once
	void finished();

	/// last emitted singal, before Functor consider fully done, emitted only once
	void final();

private Q_SLOTS:
	virtual void run();


protected:
	void state(State state);

private:

	QString name_;

	State state_ = State::_none_;
};

/// FunctorSequence is essentially a subroutine-like construct that allow you to call sequence of function one-by-one. It guarantee that the next in line function will be called only after all previous calls has finished.
class FunctorSequence : public Functor
{
    Q_OBJECT

public:
	using ValueType = FunctorSP;

	using Functor::Functor;
	FunctorSequence(QString const & name, std::initializer_list<ValueType> sequence);

	void add_functor(ValueType const &v);

	// return <current_value, max_value> for execution progress
	std::pair<int, int> progress() const override;


Q_SIGNALS:
	/// execute each time when new item starts running
	void next();

public Q_SLOTS:


private Q_SLOTS:
	void run() override;

private Q_SLOTS:

private:
	std::vector<ValueType> sequence_;
};


/// FunctorASyncSequence: same as FunctorSequence but start excution of all functors at once
class FunctorASyncSequence : public Functor
{
    Q_OBJECT

public:
	using ValueType = FunctorSP;

	using Functor::Functor;
	FunctorASyncSequence(QString const & name, std::initializer_list<ValueType> sequence);

	void add_functor(ValueType const &v);

	// return <current_value, max_value> for execution progress
	std::pair<int, int> progress() const override;

private Q_SLOTS:
	void run() override;

	void item_finished();

private:
	std::vector<ValueType> sequence_;
};


/// Functor to execute NetworkCall operation untill success
/// Network call will be initiated though user provided CallCallback
class FunctorNetworkCall : public Functor
{
    Q_OBJECT
public:
	using CallCallback = std::function< void (NetworkCall &) >;

	FunctorNetworkCall(QString const &name, CallCallback const & f = nullptr, QObject *parent = Q_NULLPTR);

	void callback(CallCallback const &callback) { callback_ = callback; }

	QByteArray result() const { return request_.result(); }

	int status_code() const { return request_.status_code(); }

	QJsonDocument result_as_json() const { return request_.result_as_json(); }

private Q_SLOTS:
	void run() override;

private:
	CallCallback callback_;

	NetworkCall request_;
};



} // namespace task
} // namespace ui
