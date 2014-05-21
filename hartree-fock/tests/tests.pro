#-------------------------------------------------
#
# Project created by QtCreator 2014-05-19T21:04:33
#
#-------------------------------------------------

TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

unix: LIBS += -L$$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/lib/ -larmadillo

INCLUDEPATH += $$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/include
DEPENDPATH += $$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/include

unix: LIBS += -L$$OUT_PWD/../src/ -lhartree-fock

INCLUDEPATH += $$PWD/../src
DEPENDPATH += $$PWD/../src

unix: LIBS += -L$$PWD/../../../../../../usr/local/lib/ -lUnitTest++

INCLUDEPATH += $$PWD/../../../../../../usr/local/include
DEPENDPATH += $$PWD/../../../../../../usr/local/include

unix: PRE_TARGETDEPS += $$PWD/../../../../../../usr/local/lib/libUnitTest++.a
