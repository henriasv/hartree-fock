TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += -std=c++11

SOURCES += main.cpp


win32:CONFIG(release, debug|release): LIBS += -L$$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/lib/release/ -larmadillo
else:win32:CONFIG(debug, debug|release): LIBS += -L$$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/lib/debug/ -larmadillo
else:unix: LIBS += -L$$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/lib/ -larmadillo

INCLUDEPATH += $$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/include
DEPENDPATH += $$PWD/../../../../../../usr/local/Cellar/armadillo/4.000.4/include

INCLUDEPATH += $$PWD/../src
DEPENDPATH += $$PWD/../src


unix: LIBS += -L$$OUT_PWD/../src/ -lhartree-fock
