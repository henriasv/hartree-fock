TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

unix: LIBS += -L/usr/local/Cellar/armadillo/4.000.4/lib/ -larmadillo

INCLUDEPATH += /usr/local/Cellar/armadillo/4.000.4/include
DEPENDPATH += /usr/local/Cellar/armadillo/4.000.4/include

unix: LIBS += -L/usr/local/lib/ -lUnitTest++

INCLUDEPATH += /usr/local/include
DEPENDPATH += /usr/local/include

unix: PRE_TARGETDEPS += /usr/local/lib/libUnitTest++.a
unix: LIBS += -L$$OUT_PWD/../src/ -lhartree-fock

INCLUDEPATH += $$PWD/../src
DEPENDPATH += $$PWD/../src
