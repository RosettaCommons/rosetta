QT += core gui widgets

CONFIG += object_parallel_to_source c++11

TARGET = workbench
TEMPLATE = app

DEFINES += BOOST_ERROR_CODE_HEADER_ONLY BOOST_SYSTEM_NO_DEPRECATED BOOST_MATH_NO_LONG_DOUBLE_MATH_FUNCTIONS PTR_STD UNUSUAL_ALLOCATOR_DECLARATION

INCLUDEPATH = $$PWD/../../../../src $$PWD/../../../../src/platform/macos

QMAKE_CXXFLAGS += \
    -isystem$$PWD/../../../../external \
    -isystem$$PWD/../../../../external/include \
    -isystem$$PWD/../../../../external/boost_1_55_0 \
    -isystem$$PWD/../../../../external/dbio \
    -isystem$$PWD/../../../../external/dbio/sqlite3 \
    -isystem$$PWD/../../../../external/libxml2/include


SOURCES += \
        main.cpp \
        mainwindow.cpp


HEADERS  += mainwindow.h

FORMS    += mainwindow.ui


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
        -L$$OUT_PWD/../rosetta/protocols_a_4    -lprotocols_a_4 \
        -L$$OUT_PWD/../rosetta/protocols_b_4    -lprotocols_b_4 \
        -L$$OUT_PWD/../rosetta/protocols_c_4    -lprotocols_c_4 \
        -L$$OUT_PWD/../rosetta/protocols_d_4    -lprotocols_d_4 \
        -L$$OUT_PWD/../rosetta/protocols_e_4    -lprotocols_e_4 \
        -L$$OUT_PWD/../rosetta/protocols_f_4    -lprotocols_f_4 \
        -L$$OUT_PWD/../rosetta/protocols_g_4    -lprotocols_g_4 \
        -L$$OUT_PWD/../rosetta/protocols_h_4    -lprotocols_h_4 \
        -L$$OUT_PWD/../rosetta/protocols_a_5    -lprotocols_a_5 \
        -L$$OUT_PWD/../rosetta/protocols_b_5    -lprotocols_b_5 \
        -L$$OUT_PWD/../rosetta/protocols_c_5    -lprotocols_c_5 \
        -L$$OUT_PWD/../rosetta/protocols_d_5    -lprotocols_d_5 \
        -L$$OUT_PWD/../rosetta/protocols_e_5    -lprotocols_e_5 \
        -L$$OUT_PWD/../rosetta/protocols_6      -lprotocols_6 \
        -L$$OUT_PWD/../rosetta/protocols_7      -lprotocols_7 \
        -L$$OUT_PWD/../ui                       -lui
