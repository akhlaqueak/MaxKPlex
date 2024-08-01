# ********     Author: Lijun Chang    ******
# ******** Email: ljchang@outlook.com ******
#
CC=g++ -flto -O3
CFLAGS=-c -I. -std=c++20

debug: CC+= -g
	all

prof: CC+= -pg
	all

all: kPlexS

kPlexS: 
	${CC} Graph.cpp -o kPlexTBs


