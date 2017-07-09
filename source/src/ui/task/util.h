#pragma once

#include <ui/util/exception.h>

#include <map>
#include <memory>

#include <QDataStream>

namespace ui {
namespace task {


quint64 const _std_map_QDataStream_magic_number_   = 0xFFF2C04ABB492107;


template <typename K, typename T>
QDataStream &operator<<(QDataStream &out, std::map<K, std::shared_ptr<T> > const&m)
{
	out << _std_map_QDataStream_magic_number_;

	out << (qint64) m.size();

	for(auto const & p : m) out << p.first << *p.second;

	return out;
}


template <typename K, typename T>
QDataStream &operator>>(QDataStream &in, std::map<K, std::shared_ptr<T> > &r)
{
	quint64 magic;
	in >> magic;
	if( magic != _std_map_QDataStream_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _Node_magic_number_: read %1, was expecting %2...").arg(magic).arg(_std_map_QDataStream_magic_number_) );

	using Map = std::map<K, std::shared_ptr<T> >;
	Map m;

	qint64 size;
	in >> size;

	for(int i=0; i<size; ++i) {
		std::pair<K, std::shared_ptr<T>> v;
		v.second = std::make_shared<T>();

		in >> v.first >> *v.second;

		m.insert(  std::move(v) );
	}

	std::swap(m, r);
	return in;
}


} // namespace task
} // namespace ui
