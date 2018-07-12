#pragma once

#include <ui/network/argument.fwd.h>

#include <ui/network/bowman_thread.h>

#include <core/pose/Pose.fwd.hh>

#include <json.hpp>

#include <QLabel>
#include <QPointer>

namespace ui {
namespace network {


/// Argument specification for HAL RPC, will be used as `prototype` when constructing particular Call
/// we use protected inheritance from QObject so derived classes could create slots if they need to
class Argument : public QObject
{
    Q_OBJECT

protected:
	Argument();
	//Argument(json const &);

public:

	~Argument() override;

	//std::string name() const { return name_; }
	//bool optional() const { return optional_; }
	std::string const & description() const { return description_; }


	/// factory method: tries to analyze JSON object and build on of Argument's derived class from it.
	/// return null on failure
	static ArgumentSP create(json const &);


	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	virtual ArgumentSP clone() const = 0;


	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	virtual QWidget * create_input_widget(QWidget * parent = nullptr) const = 0;


	/// Encode argument into JSON object
	virtual json encode() const = 0;


	/// qDebug helper
	virtual operator QString() const = 0;


	// init from JSON, note that this version is always succeeded
	void init(json const &);


protected:
	Argument(Argument const &);

private:
	//bool optional_ = true;

	std::string description_;
};


class BooleanArgument : public Argument
{
    Q_OBJECT

public:

	// try to init from JSON, return true on success and false otherwise
	bool init(json const &);


	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	ArgumentSP clone() const override;


	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	QWidget * create_input_widget(QWidget * parent = nullptr) const override;


	/// Encode argument into JSON object
	json encode() const override;


	/// qDebug helper
	operator QString() const override;


private Q_SLOTS:
	void on_toggled(bool checked);


private:
	bool default_ = false;

	bool value_ = false;
};



class IntegerArgument : public Argument
{
    Q_OBJECT

public:

	// try to init from JSON, return true on success and false otherwise
	bool init(json const &);


	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	ArgumentSP clone() const override;


	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	QWidget * create_input_widget(QWidget * parent = nullptr) const override;


	/// Encode argument into JSON object
	json encode() const override;


	/// qDebug helper
	operator QString() const override;


private Q_SLOTS:
	void on_value_changed(int value);


private:
	int default_ = 0;

	int value_ = 0;

	int min_ = 0;
	int max_ = 0;
};


class FloatArgument : public Argument
{
    Q_OBJECT

public:

	// try to init from JSON, return true on success and false otherwise
	bool init(json const &);

	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	ArgumentSP clone() const override;

	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	QWidget * create_input_widget(QWidget * parent = nullptr) const override;


	/// Encode argument into JSON object
	json encode() const override;


	operator QString() const override;

private Q_SLOTS:
	void on_value_changed(double value);

private:
	double default_ = 0.0;

	double value_ = 0.0;

	double min_ = 0.0;
	double max_ = 0.0;
};


class StringArgument : public Argument
{
    Q_OBJECT

public:

	// try to init from JSON, return true on success and false otherwise
	bool init(json const &);

	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	ArgumentSP clone() const override;

	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	QWidget * create_input_widget(QWidget * parent = nullptr) const override;


	/// Encode argument into JSON object
	json encode() const override;

	operator QString() const override;

	void value(std::string const &v) { value_  = v; }
	std::string value() const { return value_; }

	std::string get_default() const { return default_; } // `default` is now reserverd keyword so we have to use different name

private Q_SLOTS:
	void on_text_changed(QString const &);

private:
	std::string value_, default_;
};


class PathArgument : public StringArgument
{
    Q_OBJECT

public:
	enum class Kind {file, directory};

	PathArgument(Kind kind);


	// try to init from JSON, return true on success and false otherwise
	bool init(json const &);

	/// clone Argument. Note that underlying QObject is not clonned so you will need to re-connect any signal/slots on a new instance manually
	ArgumentSP clone() const override;

	/// Create widget that could be used to edit Argument value. Created widget will be connection to `this` instance so `this->value` could be updated transperently
	QWidget * create_input_widget(QWidget * parent = nullptr) const override;


	/// Encode argument into JSON object
	//json encode() const override;

	operator QString() const override;


private Q_SLOTS:
	void on_button_clicked(bool);

private:
	void update_label_text(QString const &) const;

	mutable QPointer<QLabel>  label_;

	Kind kind_ = Kind::file;
};



} // namespace network
} // namespace ui
