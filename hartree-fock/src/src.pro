#-------------------------------------------------
#
# Project created by QtCreator 2014-05-19T20:21:28
#
#-------------------------------------------------
CONFIG -= app_bundle
CONFIG-= qt
CONFIG += c++11
TARGET = hartree-fock
TEMPLATE = lib

DEFINES += SRC_LIBRARY

SOURCES += \
    primitive.cpp \
    integrator.cpp \
    hartreefocksolver.cpp \
    electronicsystem.cpp \
    berylliumhf.cpp \
    atomicorbitals.cpp \
    boysfunction.cpp

HEADERS += \
    primitive.h \
    integrator.h \
    hartreefocksolver.h \
    electronicsystem.h \
    berylliumhf.h \
    atomicorbitals.h \
    boysfunction.h

OTHER_FILES +=

unix: LIBS += -L/usr/local/Cellar/armadillo/4.000.4/lib/ -larmadillo
INCLUDEPATH += /usr/local/Cellar/armadillo/4.000.4/include
DEPENDPATH += /usr/local/Cellar/armadillo/4.000.4/include
