all: popcorn_demo test_max_substrings

popcorn_demo: popcorn_demo.cpp popcorn.cpp
	g++ -std=c++11 popcorn_demo.cpp -o popcorn_demo

test_max_substrings: test_max_substrings.cpp popcorn.cpp
	g++ -std=c++11 test_max_substrings.cpp -o test_max_substrings

clean:
	rm popcorn_demo test_max_substrings
