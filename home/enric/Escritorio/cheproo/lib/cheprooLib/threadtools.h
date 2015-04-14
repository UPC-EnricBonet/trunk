#ifndef THREADTOOLS_H
#define THREADTOOLS_H


#include <functional>
#ifdef IsMultiThreaded
#include <thread>
#include <iostream>

class threadtools
{
public:
  static void serial_for(int ibegin, int iend, std::function<void (int) > afunc)
  {
      for (int i = ibegin; i < iend; i++)
      {
	  afunc(i);
      }
  };


  static void parallel_for(int ibegin, int iend, std::function<void (int) > afunc)
  {

      unsigned long const length=iend - ibegin;
      if(!length) return;
      
      unsigned long const min_per_thread = 25;
      unsigned long const max_threads    = (length+min_per_thread-1)/min_per_thread;

      unsigned long const hardware_threads=std::thread::hardware_concurrency();

      unsigned long const num_threads=std::min(hardware_threads!=0?hardware_threads:2,max_threads);

      unsigned long const block_size=length/num_threads;  


      std::vector<std::thread> threads;

      int block_start=ibegin;
      
      for(unsigned long i=0;i< num_threads;i++)
      {
	  int block_end = block_start + block_size;

	  threads.push_back(std::thread(threadtools::serial_for, block_start, block_end, afunc));

	  block_start=block_end;
      }    


    // join threads
      for(unsigned long i=0;i<(num_threads);++i)
      {
	  threads[i].join();
      }
    return;
  };
};


#endif
#ifndef IsMultiThreaded
 class threadtools
{
public:
  static void serial_for(int ibegin, int iend, std::function<void (int) > afunc)
  {
      for (int i = ibegin; i < iend; i++)
      {
	    afunc(i);
      }
  };

  // if we do not use multithreading, then the parallel_for is actually serial
  static void parallel_for(int ibegin, int iend, std::function<void (int) > afunc)
  {
      for (int i = ibegin; i < iend; i++)
      {
	    afunc(i);
      }
  };
 };

#endif
#endif
