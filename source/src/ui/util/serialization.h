#pragma once

#include <QDataStream>

#include <map>
#include <vector>

namespace ui {
namespace util {

quint64 const _Node_magic_number_                     = 0x1A5CD6BF94D5E57F;
quint64 const _TaskSyncer_NodeStrategy_magic_number_  = 0xFFF2D9FACD201111;
quint64 const _TaskSyncer_TaskStrategy_magic_number_  = 0xFFF2D9FACD201121;
quint64 const _Task_magic_number_                     = 0xFFF2D9FACD279EE1;
quint64 const _Project_magic_number_                  = 0xFFF1D6BF94D5E57F;

quint64 const _std_map_QDataStream_magic_number_      = 0xFFF2C04ABB492107;
quint64 const _std_vector_QDataStream_magic_number_   = 0xFFF2C04ABB492108;

} // namespace util
} // namespace ui

// we define operator<< and operator>> in same namespace as QDataStream

template <typename K, typename T>
QDataStream &operator<<(QDataStream &out, std::map<K, std::shared_ptr<T> > const&m)
{
	out << ui::util::_std_map_QDataStream_magic_number_;

	out << (qint64) m.size();

	for(auto const & p : m) out << p.first << *p.second;

	return out;
}


template <typename K, typename T>
QDataStream &operator>>(QDataStream &in, std::map<K, std::shared_ptr<T> > &r)
{
	quint64 magic;
	in >> magic;
	if( magic != ui::util::_std_map_QDataStream_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _std_map_QDataStream_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_std_map_QDataStream_magic_number_) );

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



template <typename T>
QDataStream &operator<<(QDataStream &out, std::vector<std::shared_ptr<T> > const &v)
{
	out << ui::util::_std_vector_QDataStream_magic_number_;

	out << (qint64) v.size();

	for(auto const & e : v) out << *e;

	return out;
}

template <typename T>
QDataStream &operator>>(QDataStream &in, std::vector<std::shared_ptr<T> > &r)
{
	quint64 magic;
	in >> magic;
	if( magic != ui::util::_std_vector_QDataStream_magic_number_ ) throw ui::util::BadFileFormatException( QString("Invalid _std_vector_QDataStream_magic_number_: read %1, was expecting %2...").arg(magic).arg(ui::util::_std_vector_QDataStream_magic_number_) );

	using Vector = std::vector< std::shared_ptr<T> >;
	Vector v;

	qint64 size;
	in >> size;

	for(int i=0; i<size; ++i) {
		auto e = std::make_shared<T>();
		in >> *e;
		v.push_back( std::move(e) );
	}

	std::swap(v, r);
	return in;
}
