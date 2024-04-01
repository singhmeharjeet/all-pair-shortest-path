#ifndef UTILS_H
#define UTILS_H

#include <limits.h>

#include <atomic>
#include <condition_variable>
#include <iostream>
#include <mutex>

#include "cxxopts.h"
#include "get_time.h"

struct CustomBarrier {
	int num_of_threads_;
	int current_waiting_;
	int barrier_call_;
	std::mutex my_mutex_;
	std::condition_variable my_cv_;

	CustomBarrier(int t_num_of_threads)
		: num_of_threads_(t_num_of_threads), current_waiting_(0), barrier_call_(0) {}

	void wait() {
		std::unique_lock<std::mutex> u_lock(my_mutex_);
		int c = barrier_call_;
		current_waiting_++;
		if (current_waiting_ == num_of_threads_) {
			current_waiting_ = 0;
			// unlock and send signal to wake up
			barrier_call_++;
			u_lock.unlock();
			my_cv_.notify_all();
			return;
		}
		my_cv_.wait(u_lock, [&] { return (c != barrier_call_); });
		//  Condition has been reached. return
	}
};

#endif
