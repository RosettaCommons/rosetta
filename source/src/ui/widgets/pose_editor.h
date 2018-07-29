#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include <ui/network/bowman.h>

#include <core/pose/Pose.fwd.hh>


#include <utility/json_utilities.hh>


#include <QMainWindow>

#include <QListWidget>

namespace Ui {
class PoseEditor;
}

namespace ui {
namespace widgets {

class PoseEditor : public QMainWindow
{
    Q_OBJECT

public:
    explicit PoseEditor(QWidget *parent = 0);
    ~PoseEditor();

private Q_SLOTS:
	//void client_connected(std::string const &);
	//void client_disconnected(std::string const &);
	//void specification_received(std::string const &, JSON_CSP const &j);

	void result_received(core::pose::PoseOP const &, JSON_CSP const &);
	void progress_data_received(core::pose::PoseOP const &, JSON_CSP const &);

    void on_action_open_pose_triggered();
    void on_action_save_pose_triggered();
    void on_action_save_pose_as_triggered();

    void on_apply_clicked();
    void on_abort_clicked();
    void on_color_toggled(bool);
	void on_pause_toggled(bool);

	void on_functions_double_clicked(ui::network::FunctionID const &);


private:
	bool save_pose(bool always_ask_for_file_name);

	QString file_name_;

	//std::map<std::string, QListWidget *> back_ends_;
	//int hal_unique_index_ = 1;

    Ui::PoseEditor *ui;
};

} // namespace widgets
} // namespace ui

#endif // MAINWINDOW_H
