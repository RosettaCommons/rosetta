#pragma once

#include <QException>
#include <QDebug>

namespace ui {
namespace util {


class BadFileFormatException : public QException
{
public:
	BadFileFormatException(QString const &msg="") : msg_(msg) {
		qDebug() << "BadFileFormatException: " << msg_;
	}

    void raise() const {
		throw *this;
	}

    BadFileFormatException *clone() const { return new BadFileFormatException(*this); }

private:
	QString msg_;
};

} // namespace util
} // namespace ui
