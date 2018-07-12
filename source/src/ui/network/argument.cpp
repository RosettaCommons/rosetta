#include <ui/network/argument.h>

#include <protocols/network/util.hh>

#include <utility/json_utilities.hh>

#include <QSpinBox>
#include <QLineEdit>
#include <QCheckBox>
//#include <QVBoxLayout>
#include <QHBoxLayout>
#include <QPushButton>
#include <QFileDialog>

#include <QDebug>

#include <limits>
#include <algorithm>

using std::string;


using namespace protocols::network;

using utility::extract_value_if_present;

namespace ui {
namespace network {

Argument::Argument() = default;


// Argument::Argument(json const &j)
// {
// 	auto it_optional = j.find(_f_optional_);
// 	if( it_optional != j.end()  and  it_optional->is_boolean() ) optional_ = it_optional.value();
// }


Argument::~Argument() = default;


Argument::Argument(Argument const &other) : QObject()
{
	//name_     = other.name_;
	//optional_ = other.optional_;
	description_ = other.description_;
}


void Argument::init(json const &j)
{
	//extract_value_if_present(j, _f_optional_, optional_);
	extract_value_if_present(j, _f_description_, description_);
}


template <typename A, class ... Types>
std::function< ArgumentSP(json const &) > create_argument(Types ... args)
{
	//qDebug() << "create_argument:" << j.dump(2).c_str();

	auto creator = [args...](json const &j) -> std::shared_ptr<A>
		{
		 std::shared_ptr<A> o = std::make_shared<A>(args...);
		 if( o->init(j) ) return o;
		 else return std::shared_ptr<A>();
		};
	return creator;
}


/// factory method: tries to analyze JSON object and build on of Argument's derived class from it.
/// return null on failure
ArgumentSP Argument::create(json const &j)
{
	//qDebug() << "Argument::create:" << j.dump(2).c_str();

	//auto it_name = j.find(_f_name_);
	auto it_type = j.find(_f_type_);

	if( it_type != j.end() ) {

		std::map<string, std::function< ArgumentSP(json const &) > > creators = {
			{ _t_boolean_, create_argument<BooleanArgument>() },
			{ _t_integer_, create_argument<IntegerArgument>() },
			{ _t_float_,   create_argument<FloatArgument>()   },
			{ _t_string_,  create_argument<StringArgument>()  },

			{ _t_file_,      create_argument<PathArgument>(PathArgument::Kind::file) },
			{ _t_directory_, create_argument<PathArgument>(PathArgument::Kind::directory) },

			{ _t_pose_,    create_argument<StringArgument>() }, /// place holder for Pose-type payload
		};

		auto it = creators.find( it_type.value() );
		if( it != creators.end() ) return it->second(j);
	}

	return ArgumentSP();
}




/// --------------------------------
/// BooleanArgument
/// --------------------------------
bool BooleanArgument::init(json const &j)
{
	//qDebug() << "BooleanArgument::init:" << j.dump(2).c_str();

	Argument::init(j);

	if(  extract_value_if_present(j, _f_default_, default_) ) value_ = default_;

	//qDebug() << "BooleanArgument::init: default=" << default_;

	return true;
}

ArgumentSP BooleanArgument::clone() const
{
	return std::make_shared<BooleanArgument>(*this);
}


void BooleanArgument::on_toggled(bool checked)
{
	//qDebug() << "BooleanArgument::toggled:" << checked;
	value_ = checked;
}


QWidget * BooleanArgument::create_input_widget(QWidget * parent) const
{
	auto w = new QCheckBox(parent);
	w->setChecked(value_);

	connect(w, &QCheckBox::toggled, this, &BooleanArgument::on_toggled);

	return w;
}


/// Encode argument into JSON object
json BooleanArgument::encode() const
{
	return json(value_);
}


BooleanArgument::operator QString() const
{
	auto r = QString("BooleanArgument(%1, default=%2)").arg(value_).arg(default_);
	return r + ")";
}



/// --------------------------------
/// IntegerArgument
/// --------------------------------
bool IntegerArgument::init(json const &j)
{
	//qDebug() << "IntegerArgument::init:" << j.dump(2).c_str();

	Argument::init(j);

	// min_ = j.value(_f_min_, std::numeric_limits<int>::lowest() );
	// max_ = j.value(_f_max_, std::numeric_limits<int>::max() );

	min_ = std::numeric_limits<int>::lowest();
	max_ = std::numeric_limits<int>::max();

	extract_value_if_present(j, _f_min_, min_);
	extract_value_if_present(j, _f_max_, max_);

	if(  extract_value_if_present(j, _f_default_, default_) ) value_ = default_;
	else default_ = value_ = std::min( std::max(0, min_), max_);

	//qDebug() << "IntegerArgument::init: min=" << min_ << " max=" << max_ << " default=" << default_;

	return true;
}

ArgumentSP IntegerArgument::clone() const
{
	return std::make_shared<IntegerArgument>(*this);
}


void IntegerArgument::on_value_changed(int value)
{
	value_ = value;
}


QWidget * IntegerArgument::create_input_widget(QWidget * parent) const
{
	auto w = new QSpinBox(parent);
	w->setRange(min_, max_);
	w->setValue(value_);

	connect(w, static_cast< void (QSpinBox::*)(int) >( &QSpinBox::valueChanged ), this, &IntegerArgument::on_value_changed);

	return w;
}


/// Encode argument into JSON object
json IntegerArgument::encode() const
{
	return json(value_);
}


IntegerArgument::operator QString() const
{
	auto r = QString("IntegerArgument(%1, default=%2").arg(value_).arg(default_);
	if( min_ != std::numeric_limits<int>::lowest() ) r += QString("min=%1, ").arg(min_);
	if( max_ != std::numeric_limits<int>::max() ) r += QString("max=%1, ").arg(max_);
	return r + ")";
}


/// --------------------------------
/// FloatArgument
/// --------------------------------
bool FloatArgument::init(json const &j)
{
	//qDebug() << "FloatArgument::init:" << j.dump(2).c_str();

	Argument::init(j);

	//min_ = j.value(_f_min_, std::numeric_limits<double>::lowest() );
	//max_ = j.value(_f_max_, std::numeric_limits<double>::max() );

	// string s_min = j.value(_f_min_, "");
	// if( not s_min.empty() ) min_ = std::stod(s_min);

	// string s_max = j.value(_f_max_, "");
	// if( not s_max.empty() ) max_ = std::stod(s_max);

	min_ = std::numeric_limits<double>::lowest();
	max_ = std::numeric_limits<double>::max();

	extract_value_if_present(j, _f_min_, min_);
	extract_value_if_present(j, _f_max_, max_);

	if(  extract_value_if_present(j, _f_default_, default_) ) value_ = default_;
	else default_ = value_ = std::min( std::max(0.0, min_), max_);

	// qDebug() << "FloatArgument::init: min=" << min_ << " max=" << max_ << " default=" << default_;

	// qDebug() << "FloatArgument::init: std::max(0.0, min_)=" << std::max<double>(0.0, min_);

	// qDebug() << "FloatArgument::init: std::numeric_limits<int>::min()=" << std::numeric_limits<int>::min();
	// qDebug() << "FloatArgument::init: std::numeric_limits<double>::min()=" << std::numeric_limits<double>::min();
	// qDebug() << "FloatArgument::init: std::numeric_limits<double>::lowest()=" << std::numeric_limits<double>::lowest();

	return true;
}

ArgumentSP FloatArgument::clone() const
{
	return std::make_shared<FloatArgument>(*this);
}

void FloatArgument::on_value_changed(double value)
{
	value_ = value;
}

QWidget * FloatArgument::create_input_widget(QWidget * parent) const
{
	auto w = new QDoubleSpinBox(parent);
	w->setRange(min_, max_);
	w->setDecimals(5);
	w->setValue(value_);
	if( max_ - min_ < 10 ) w->setSingleStep( (max_ - min_) / 100.0 );

	connect(w, static_cast< void (QDoubleSpinBox::*)(double) >( &QDoubleSpinBox::valueChanged ), this, &FloatArgument::on_value_changed);

	//qDebug() << "FloatArgument::create_input_widget:" << value_;

	return w;
}

/// Encode argument into JSON object
json FloatArgument::encode() const
{
	return json(value_);
}

FloatArgument::operator QString() const
{
	auto r = QString("FloatArgument(%1, default=%2").arg(value_).arg(default_);
	if( min_ != std::numeric_limits<double>::lowest() ) r += QString("min=%1, ").arg(min_);
	if( max_ != std::numeric_limits<double>::max() ) r += QString("max=%1, ").arg(max_);
	return r + ")";
}


/// --------------------------------
/// StringArgument
/// --------------------------------
bool StringArgument::init(json const &j)
{
	//qDebug() << "StringArgument::init:" << j.dump(2).c_str();
	Argument::init(j);

	if(  extract_value_if_present(j, _f_default_, default_) ) value_ = default_;

	return true;
}

ArgumentSP StringArgument::clone() const
{
	return std::make_shared<StringArgument>(*this);
}

void StringArgument::on_text_changed(QString const &text)
{
	value_ = text.toStdString();
}

QWidget * StringArgument::create_input_widget(QWidget * parent) const
{
	auto w = new QLineEdit(parent);
	w->setText( QString::fromStdString(value_) );

	connect(w, &QLineEdit::textChanged, this, &StringArgument::on_text_changed);


	return w;
}

/// Encode argument into JSON object
json StringArgument::encode() const
{
	return json(value_);
}

StringArgument::operator QString() const
{
	return QString("StringArgument(%1, default=%2)").arg( QString::fromStdString(value_) ).arg( QString::fromStdString(default_) );
}


/// --------------------------------
/// PathArgument
/// --------------------------------
PathArgument::PathArgument(Kind kind)
{
	kind_ = kind;
}

bool PathArgument::init(json const &j)
{
	//qDebug() << "PathArgument::init:" << j.dump(2).c_str();

	StringArgument::init(j);

	string type;
	if(  extract_value_if_present(j, _f_type_, type) ) kind_ = (type == _t_file_) ? Kind::file : Kind::directory;

	return true;
}

ArgumentSP PathArgument::clone() const
{
	return std::make_shared<PathArgument>(*this);
}


QWidget * PathArgument::create_input_widget(QWidget * parent) const
{
	auto w = new QWidget(parent);
	w->setMinimumWidth(256);

	auto layout = new QHBoxLayout();

	label_ = new QLabel(w);
	auto button = new QPushButton("set", w);

	layout->addWidget(label_);
	layout->addWidget(button);

	w->setLayout(layout);

	update_label_text( QString::fromStdString( value() ) );

	connect(button, &QPushButton::clicked, this, &PathArgument::on_button_clicked);  // QAbstractButton

	return w;
}


void PathArgument::on_button_clicked(bool)
{
	QString file_name;

	if( kind_ == Kind::file ) file_name = QFileDialog::getOpenFileName(nullptr, tr("Open file"), "", tr("All files (*.*)"), Q_NULLPTR/*, QFileDialog::DontUseNativeDialog*/);
	else file_name = QFileDialog::getExistingDirectory(nullptr, tr("Select directory"), "", QFileDialog::ShowDirsOnly | QFileDialog::DontResolveSymlinks);

	if( not file_name.isEmpty() ) {
		value( file_name.toStdString() );
		update_label_text(file_name);
	}
}


void PathArgument::update_label_text(QString const &v) const
{
	if(label_) label_->setText( v.isEmpty() ? "None" : QFileInfo(v).fileName() );
}


PathArgument::operator QString() const
{
	return QString("PathArgument(%1, kind=%2 default=%3)").arg( QString::fromStdString( value() ) ).arg(kind_ == Kind::file ? "file" : "directory").arg( QString::fromStdString( get_default() ) );
}

} // namespace network
} // namespace ui
