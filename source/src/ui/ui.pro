QT += core gui network widgets charts

CONFIG += object_parallel_to_source c++11 no_keywords

TARGET = ui
TEMPLATE = lib

DEFINES += BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS PTR_STD UNUSUAL_ALLOCATOR_DECLARATION

INCLUDEPATH = $$PWD/../../src $$PWD/../../src/platform/macos

QMAKE_CXXFLAGS += \
    -isystem$$PWD/../../external \
    -isystem$$PWD/../../external/include \
    -isystem$$PWD/../../external/boost_1_55_0 \
    -isystem$$PWD/../../external/dbio \
    -isystem$$PWD/../../external/dbio/sqlite3 \
    -isystem$$PWD/../../external/libxml2/include


SOURCES += \
    ui_lib_test.cpp \
    config/config_dialog.cpp \
    config/util.cpp \
    task/job_submit.cpp \
    task/job_view.cpp \
    task/file.cpp \
    task/functor.cpp \
    task/node.cpp \
    task/project.cpp \
    task/project_model.cpp \
    task/project_view.cpp \
    task/task.cpp \
    task/task_syncer.cpp \
    task/task_submit.cpp \
    task/task_view.cpp \
    task/util.cpp \
    ui_protocols/helical_bundle/HelicalBundleDialogueWidget.cpp \
    ui_protocols/helical_bundle/HelixOptionWidget.cpp \
    ui_core/pose_draw/SimplePoseDrawOpenGLWidget.cpp \
    ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.cpp \
    viewers/score_file_view.cpp


HEADERS  += \
    ui_lib_test.h \
    config/config_dialog.h \
    config/util.h \
    task/job_submit.h \
    task/job_view.h \
    task/file.fwd.h \
    task/file.h \
    task/functor.fwd.h \
    task/functor.h \
    task/node.h \
    task/project.fwd.h \
    task/project.h \
    task/project_model.h \
    task/task_syncer.h \
    task/task_submit.h \
    task/project_view.h \
    task/task.fwd.h \
    task/task.h \
    task/util.h \
    task/task_view.h \
    ui_protocols/helical_bundle/HelicalBundleDialogueWidget.h \
    ui_protocols/helical_bundle/HelixOptionWidget.h \
    ui_protocols/helical_bundle/HelicalBundleDialogueWidget.fwd.h \
    ui_protocols/helical_bundle/HelixOptionWidget.fwd.h \
    ui_core/pose_draw/SimplePoseDrawOpenGLWidget.h \
    ui_core/pose_draw/SimplePoseDrawOpenGLWidget.fwd.h \
    ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.h \
    ui_protocols/helical_bundle/HelicalBundlePoseDrawOpenGLWidget.fwd.h \
    util/exception.h \
    util/serialization.h \
    viewers/score_file_view.h

FORMS    += \
    config/config_dialog.ui \
    task/job_submit.ui \
    task/job_view.ui \
    task/task_submit.ui \
    task/project_view.ui \
    task/task_view.ui \
    ui_protocols/helical_bundle/HelicalBundleDialogueWidget.ui \
    ui_protocols/helical_bundle/HelixOptionWidget.ui \
    viewers/score_file_view.ui

LIBS += \
        -L$$OUT_PWD/../rosetta/external         -lexternal \
        -L$$OUT_PWD/../rosetta/cifparse         -lcifparse \
        -L$$OUT_PWD/../rosetta/libxml2          -llibxml2 -lz \
        -L$$OUT_PWD/../rosetta/ObjexxFCL        -lObjexxFCL \
        -L$$OUT_PWD/../rosetta/utility          -lutility \
        -L$$OUT_PWD/../rosetta/numeric          -lnumeric \
        -L$$OUT_PWD/../rosetta/basic            -lbasic \
        -L$$OUT_PWD/../rosetta/core_1           -lcore_1 \
        -L$$OUT_PWD/../rosetta/core_2           -lcore_2 \
        -L$$OUT_PWD/../rosetta/core_3           -lcore_3 \
        -L$$OUT_PWD/../rosetta/core_4           -lcore_4 \
        -L$$OUT_PWD/../rosetta/core_5           -lcore_5 \
        -L$$OUT_PWD/../rosetta/protocols_1      -lprotocols_1 \
        -L$$OUT_PWD/../rosetta/protocols_a_2    -lprotocols_a_2 \
        -L$$OUT_PWD/../rosetta/protocols_b_2    -lprotocols_b_2 \
        -L$$OUT_PWD/../rosetta/protocols_3      -lprotocols_3 \
        -L$$OUT_PWD/../rosetta/protocols_4      -lprotocols_4 \
        -L$$OUT_PWD/../rosetta/protocols_a_5    -lprotocols_a_5 \
        -L$$OUT_PWD/../rosetta/protocols_b_5    -lprotocols_b_5 \
        -L$$OUT_PWD/../rosetta/protocols_c_5    -lprotocols_c_5 \
        -L$$OUT_PWD/../rosetta/protocols_d_5    -lprotocols_d_5 \
        -L$$OUT_PWD/../rosetta/protocols_e_5    -lprotocols_e_5 \
        -L$$OUT_PWD/../rosetta/protocols_f_5    -lprotocols_f_5 \
        -L$$OUT_PWD/../rosetta/protocols_g_5    -lprotocols_g_5 \
        -L$$OUT_PWD/../rosetta/protocols_h_5    -lprotocols_h_5 \
        -L$$OUT_PWD/../rosetta/protocols_a_6    -lprotocols_a_6 \
        -L$$OUT_PWD/../rosetta/protocols_b_6    -lprotocols_b_6 \
        -L$$OUT_PWD/../rosetta/protocols_c_6    -lprotocols_c_6 \
        -L$$OUT_PWD/../rosetta/protocols_d_6    -lprotocols_d_6 \
        -L$$OUT_PWD/../rosetta/protocols_e_6    -lprotocols_e_6 \
        -L$$OUT_PWD/../rosetta/protocols_7      -lprotocols_7 \
        -L$$OUT_PWD/../rosetta/protocols_8      -lprotocols_8

macx {
LIBS += \
        -framework OpenGL \
        -framework GLUT
}

unix:!macx {
LIBS += \
        -lglut \
        -lGLU \
        -lGL
}
