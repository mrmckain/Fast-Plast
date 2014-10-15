// $Author: benine $
// $Date$
// $Log$
// Contains the Mismatch class for afin

#include "mismatch.hpp"

// Mismatch Class Code
Mismatch::Mismatch(){
  score = 1.0;
  length = 0;
  index_i = 0;
  index_j = 0;
  end_i = 0;
  end_j = 0;
}

Mismatch::Mismatch( double score, int length, int index_i, int index_j, int end_i, int end_j ): score(score), length(length), index_i(index_i), index_j(index_j), end_i(end_i), end_j(end_j){}

// set mismatch score
void Mismatch::set_score( double score ){
  this->score = score;
}

// set length
void Mismatch::set_length( int length ){
  this->length = length;
}

// set index_i
void Mismatch::set_index_i( int index ){
  index_i = index;
}

// set index_j
void Mismatch::set_index_j( int index ){
  index_j = index;
}

// set index
void Mismatch::set_indices( int index_i, int index_j ){
  this->index_i = index_i;
  this->index_j = index_j;
}

// set end_i
void Mismatch::set_end_i( int end_i ){
  this->end_i = end_i;
}

// set end_j
void Mismatch::set_end_j( int end_j ){
  this->end_j = end_j;
}

// return mismatch score
double Mismatch::get_score(){
  return score;
}

// return length
int Mismatch::get_length(){
  return length;
}

// return index i
int Mismatch::get_index_i(){
  return index_i;
}

// return index j
int Mismatch::get_index_j(){
  return index_j;
}

//return end_i
int Mismatch::get_end_i(){
  return end_i;
}

//return end_j
int Mismatch::get_end_j(){
  return end_j;
}
