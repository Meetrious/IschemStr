#pragma once
#include <chrono>
#include <iostream>

class Timer {
	std::chrono::time_point<std::chrono::steady_clock> start, end;
public: 
	Timer() { start = std::chrono::high_resolution_clock::now(); }

	void ClickStart() {
		start = std::chrono::high_resolution_clock::now();
	}
	void ClickEnd() {
		end = std::chrono::high_resolution_clock::now();
	}

	std::chrono::duration<float> CountInterval() {
		return end - start;
	}

	~Timer() = default;
};