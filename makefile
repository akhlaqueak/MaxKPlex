# ********     Author: Lijun Chang    ******
# ******** Email: ljchang@outlook.com ******
#
CC=g++ -flto -O3 -std=c++17 -I.

all: kPlexS

kPlexS: main.cpp
	${CC} main.cpp -o kPlexS

