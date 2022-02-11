# QMake Build file
QMAKE_CC = gcc
QMAKE_LINK_C = gcc

QMAKE_CXX = g++
QMAKE_LINK = g++


QMAKE_LFLAGS_X86_64 = -arch x86_64
QMAKE_CFLAGS_X86_64 = -arch x86_64
QMAKE_CXXFLAGS_X86_64 = -arch x86_64

QMAKE_CFLAGS_RELEASE += -g
QMAKE_CXXFLAGS_RELEASE += -g -std=gnu++0x 
QMAKE_CFLAGS_DEBUG += -g -Wall -Wextra 
QMAKE_CXXFLAGS_DEBUG += -g -std=gnu++0x -Wall -Wextra 

TEMPLATE = app \
    console
CONFIG += release
CONFIG -= app_bundle
CONFIG -= qt
HEADERS += ellipsoid.hpp 

SOURCES += ellipsoid.cpp main.cpp

TARGET = boundingelipsoid
INCLUDEPATH = 

LIBS += -Lm

PKGCONFIG += 
QT -= gui \
    core
