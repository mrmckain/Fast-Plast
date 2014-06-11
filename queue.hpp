// afin FIFO queue class

#ifndef QUEUE_H
#define QUEUE_H

#include <queue>
#include <mutex>
#include <condition_variable>

using namespace std;

template <typename T>
class Queue{
  private:
    queue<T> qu;
    mutex mtx;
    condition_variable cv;
    
  public:
    // pop front item and return value
    T pop();

    // pass item by reference and return front of queue through the reference
    void pop( T& item );

    // push item onto back of queue
    void push( const T& item );
};

#endif
