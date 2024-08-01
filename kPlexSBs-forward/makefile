# ********     Author: Lijun Chang    ******
# ******** Email: ljchang@outlook.com ******
#
CC=g++ -flto -O3
CFLAGS=-c -I. -std=c++20

debug: CC+= -g
debug: kPlexS

prof: CC+= -pg
prof: kPlexS

all: kPlexS

kPlexS: 
	${CC} Graph.cpp -o kPlexTBs


