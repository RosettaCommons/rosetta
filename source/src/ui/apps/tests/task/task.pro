QT += testlib core network widgets
QT -= gui

CONFIG += qt console warn_on depend_includepath testcase
CONFIG -= app_bundle

INCLUDEPATH = $$PWD/../../../../../src $$PWD/../../../../../src/platform/macos


TEMPLATE = app

SOURCES += node_test.cpp

LIBS += \
        -L$$OUT_PWD/../../ui                       -lui \
        -L$$OUT_PWD/../../rosetta/protocols_7      -lprotocols_7 \
        -L$$OUT_PWD/../../rosetta/protocols_6      -lprotocols_6 \
        -L$$OUT_PWD/../../rosetta/protocols_e_5    -lprotocols_e_5 \
        -L$$OUT_PWD/../../rosetta/protocols_d_5    -lprotocols_d_5 \
        -L$$OUT_PWD/../../rosetta/protocols_c_5    -lprotocols_c_5 \
        -L$$OUT_PWD/../../rosetta/protocols_b_5    -lprotocols_b_5 \
        -L$$OUT_PWD/../../rosetta/protocols_a_5    -lprotocols_a_5 \
        -L$$OUT_PWD/../../rosetta/protocols_h_4    -lprotocols_h_4 \
        -L$$OUT_PWD/../../rosetta/protocols_g_4    -lprotocols_g_4 \
        -L$$OUT_PWD/../../rosetta/protocols_f_4    -lprotocols_f_4 \
        -L$$OUT_PWD/../../rosetta/protocols_e_4    -lprotocols_e_4 \
        -L$$OUT_PWD/../../rosetta/protocols_d_4    -lprotocols_d_4 \
        -L$$OUT_PWD/../../rosetta/protocols_c_4    -lprotocols_c_4 \
        -L$$OUT_PWD/../../rosetta/protocols_b_4    -lprotocols_b_4 \
        -L$$OUT_PWD/../../rosetta/protocols_a_4    -lprotocols_a_4 \
        -L$$OUT_PWD/../../rosetta/protocols_3      -lprotocols_3 \
        -L$$OUT_PWD/../../rosetta/protocols_b_2    -lprotocols_b_2 \
        -L$$OUT_PWD/../../rosetta/protocols_a_2    -lprotocols_a_2 \
        -L$$OUT_PWD/../../rosetta/protocols_1      -lprotocols_1 \
        -L$$OUT_PWD/../../rosetta/core_5           -lcore_5 \
        -L$$OUT_PWD/../../rosetta/core_4           -lcore_4 \
        -L$$OUT_PWD/../../rosetta/core_3           -lcore_3 \
        -L$$OUT_PWD/../../rosetta/core_2           -lcore_2 \
        -L$$OUT_PWD/../../rosetta/core_1           -lcore_1 \
        -L$$OUT_PWD/../../rosetta/basic            -lbasic \
        -L$$OUT_PWD/../../rosetta/numeric          -lnumeric \
        -L$$OUT_PWD/../../rosetta/utility          -lutility \
        -L$$OUT_PWD/../../rosetta/ObjexxFCL        -lObjexxFCL \
        -L$$OUT_PWD/../../rosetta/libxml2          -llibxml2 \
        -L$$OUT_PWD/../../rosetta/cifparse         -lcifparse \
        -L$$OUT_PWD/../../rosetta/external         -lexternal \
        -lz
