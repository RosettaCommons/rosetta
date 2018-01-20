#pragma once

#include <ui/util/exception.h>
#include <ui/task/project.fwd.h>

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


QNetworkAccessManager * network_access_manager();

QString server_base_url();

/// deprecated: basic url for Node-base API
QString server_url();

QString task_api_url();



/// Simple wrapper around QNetworkReply to model RPC calls
/// Store results and perform retry if needed
class NetworkCall final : public QObject
{
    Q_OBJECT
public:
	using TerminationCodes = std::vector<int>;

	NetworkCall(QObject * parent=nullptr);
	~NetworkCall();

	void set_termination_codes(TerminationCodes const & termination_codes) { termination_codes_ = termination_codes; }

	using CustomHeader = std::pair<QByteArray, QByteArray>;

	void call(QString const & address, QNetworkAccessManager::Operation operation, QJsonDocument const &post_data);
	void call(QString const & address, QNetworkAccessManager::Operation operation = QNetworkAccessManager::Operation::GetOperation, QByteArray const &post_data=QByteArray(), QString const &content_type="application/octet-stream", std::vector<CustomHeader> const & raw_headers=std::vector<CustomHeader>());

	QByteArray result() const;
	int status_code() const { return status_code_; }

	QJsonDocument result_as_json() const;

Q_SIGNALS:
	void finished();

public Q_SLOTS:

private Q_SLOTS:
	void call();
	void network_operation_is_finished();

private:
	void abort_network_operation();

	QString address_;
	QString content_type_;
	QByteArray post_data_;
	QNetworkAccessManager::Operation operation_;
	std::vector<CustomHeader> raw_headers_;

	QNetworkReply * reply_ = nullptr;

	TerminationCodes termination_codes_;

	int status_code_ = 0;
	QByteArray result_;
};



/// Try to save project. If project file name is not set then ask user for it. If always_ask_for_file_name is true - always ask for file name.
/// return true on success
bool save_project(Project &, bool always_ask_for_file_name);

/// Check if system have enough configuration for Task to be submitted. Right now this include credentials and file-name for project.
bool check_submit_requirements(Project &project);

} // namespace task
} // namespace ui
