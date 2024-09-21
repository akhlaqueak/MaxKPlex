all : 
	g++ -flto -Ofast -DNDEBUG -march=native -o Maple Graph.cpp -w
clean:
	rm -rf Maple
debug:
	g++ -g -DNDEBUG -o Maple Graph.cpp
prof:
	g++ -pg -DNDEBUG -march=native -o Maple-pg Graph.cpp -w