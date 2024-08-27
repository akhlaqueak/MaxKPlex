all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o Maple main.cpp Graph.cpp -w
clean:
	rm -rf Maple
debug:
	g++ -g -DNDEBUG -march=native -o Maple main.cpp Graph.cpp -w
prof:
	g++ -pg -DNDEBUG -march=native -o Maple main.cpp Graph.cpp -w