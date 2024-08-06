# ********     Author: Lijun Chang    ******
# ******** Email: ljchang@outlook.com ******
#
CC=g++ -flto -O3
CFLAGS=-c -I. -std=c++20

all: kPlexS

debug: CC+= -g
debug: kPlexS

prof: CC= g++ -pg
prof: kPlexS

kPlexS: 
	${CC} Graph.cpp -o kPlexTBs


