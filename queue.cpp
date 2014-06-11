// afin Queue class

#include <queue>
#include <mutex>
#include <condition_variable>
#include "queue.hpp"

using namespace std;

// pop front item and return value
template< typename T >
T Queue<T>::pop(){
  unique_lock<mutex> mlock( mtx );
  while( qu.empty() ){
    cv.wait(mlock);
  }

  auto item = qu.front();
  qu.pop();
  return item;
}

// pass item by reference and return front of queue through the reference
template< typename T >
void Queue<T>::pop( T& item ){
  unique_lock<mutex> mlock( mtx );
  while( qu.empty() ){
    cv.wait(mlock);
  }

  item = qu.front();
  qu.pop();
}

// push item onto back of queue
template< typename T >
void Queue<T>::push( const T& item ){
  unique_lock<mutex> mlock( mtx );
  qu.push( item );
  mlock.unlock();
  cv.notify_one();
}
