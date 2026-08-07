/*
 * Copyright 2011, Ben Langmead <langmea@cs.jhu.edu>
 *
 * This file is part of Bowtie 2.
 *
 * Bowtie 2 is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Bowtie 2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Bowtie 2.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef THREADING_H_
#define THREADING_H_

#include <iostream>
#include <thread>
#include <chrono>
#include <atomic>
#include <mutex>

#define MUTEX_T std::mutex

/**
 * Pause inside a spin-wait loop.
 *
 * These exist as named functions on purpose. Every wait loop in the codebase
 * spelled its own pause inline, wrapped in TinyThread++'s platform macros:
 *
 *     while(!done) {
 *     #if defined(_TTHREAD_WIN32_)
 *         Sleep(1);
 *     #elif defined(_TTHREAD_POSIX_)
 *         nanosleep(&ts, NULL);
 *     #endif
 *     }
 *
 * Those macros come from tinythread.h/fast_mutex.h, so the loop body silently
 * empties out in any translation unit that stops including them -- and an
 * empty loop polling a plain, non-atomic flag has no side effects, so an
 * optimizing compiler is free to hoist the load or delete the loop outright.
 * clang at -O3 deletes it: the whole wait becomes a single `ret`. The waiter
 * then runs straight past the worker it was supposed to wait for.
 *
 * Routing every wait through these functions takes the preprocessor out of the
 * picture, so a wait can no longer compile down to nothing. The fence makes
 * that a guarantee rather than a hope: it keeps the surrounding non-atomic
 * flag load from being hoisted out of the loop.
 *
 * This fixes the compile-away hazard, not the underlying data race -- the
 * flags these loops poll are still plain bools written by other threads
 * (blockwise_sa.h's _done, hgfm.h's ThreadParam::done, pat.h's SRA done).
 * Making those atomic is a separate change.
 */
inline void threadYield() {
	std::this_thread::yield();
	std::atomic_thread_fence(std::memory_order_acquire);
}

inline void threadSleepMs(unsigned ms) {
	std::this_thread::sleep_for(std::chrono::milliseconds(ms));
	std::atomic_thread_fence(std::memory_order_acquire);
}

/**
 * Wrap a lock; obtain lock upon construction, release upon destruction.
 */
class ThreadSafe {
public:
    ThreadSafe(MUTEX_T* ptr_mutex, bool locked = true) {
		if(locked) {
		    this->ptr_mutex = ptr_mutex;
		    ptr_mutex->lock();
		}
		else
		    this->ptr_mutex = NULL;
	}

	~ThreadSafe() {
	    if (ptr_mutex != NULL)
	        ptr_mutex->unlock();
	}
    
private:
	MUTEX_T *ptr_mutex;
};

#endif
