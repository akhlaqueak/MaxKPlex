CC=g++ -I../common/ -std=c++20 -O3

all: main

debug: CC+= -g
debug: main

run: main
	./kplex -g ..\datasets\as-caida20071105.bin -k 2

main: main.cpp
	git pull
	${CC}  main.cpp -o kplex

local: main.cpp
	${CC}  main.cpp -o kplex
	.\kplex.exe -g ..\datasets\as-caida20071105.bin -k 2
