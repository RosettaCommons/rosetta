#pragma once

#include <ui/util/exception.h>
#include <ui/task/project.fwd.h>


#include <map>
#include <memory>

#include <QDataStream>

// NetworkCall extra includes
#include <QJsonDocument>
#include <QNetworkReply>
#include <QNetworkAccessManager>
#include <QPointer>


namespace ui {
namespace task {

using TimeStamp = int64_t;

TimeStamp get_utc_timestamp();

/// calculate appropriate delay in sec until next network retry should be made
int get_retry_interval();


/// Simple wrapper around QNetworkReply to model RPC calls
/// Store results and perform retry if needed
class NetworkCall final : public QObject
{
    Q_OBJECT
public:
	NetworkCall(QObject * parent=nullptr);

	void call(QString const & address, QNetworkAccessManager::Operation operation, QJsonDocument const &post_data=QJsonDocument());
	void call(QString const & address, QNetworkAccessManager::Operation operation = QNetworkAccessManager::Operation::GetOperation, QByteArray const &post_data=QByteArray(), QString const &content_type="application/octet-stream");

	QByteArray result() const;
	QJsonDocument json() const;


Q_SIGNALS:
	void finished();

private Q_SLOTS:
	void network_operation_is_finished();

private:
	void abort_network_operation();

	QString address_;
	QString content_type_;
	QByteArray post_data_;
	QNetworkAccessManager::Operation operation_;

	QPointer<QNetworkReply> reply_;
	QByteArray result_;
};



/// Try to save project. If project file name is not set then ask user for it. If always_ask_for_file_name is true - always ask for file name.
/// return true on success
bool save_project(Project &, bool always_ask_for_file_name);

/// Check if system have enough configuration for Task to be submitted. Right now this include credentials and file-name for project.
bool check_submit_requirements(Project &project);


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
