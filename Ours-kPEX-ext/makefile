all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o Maple main.cpp Graph.cpp -w -std=c++17
clean:
	rm -rf Maple
debug:
	g++ -g -DNDEBUG -o Maple main.cpp Graph.cpp
prof:
	g++ -pg -DNDEBUG -march=native -o Maple-pg main.cpp Graph.cpp -w