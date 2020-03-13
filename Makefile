run_popcorn: popcorn.cpp
	g++ -std=c++11 popcorn.cpp -o run_popcorn

clean:
	rm run_popcorn
