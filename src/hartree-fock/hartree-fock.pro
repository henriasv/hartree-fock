TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    hartreefocksolver.cpp \
    electronicsystem.cpp \
    berylliumhf.cpp \
    atomicorbitals.cpp

HEADERS += \
    hartreefocksolver.h \
    electronicsystem.h \
    berylliumhf.h \
    atomicorbitals.h


LIBS += -L/usr/local/Cellar/armadillo/4.000.4/lib -llapack -lblas -larmadillo

INCLUDEPATH += /usr/local/Cellar/armadillo/4.000.4/include
