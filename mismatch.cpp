// $Author: benine $
// $Date$
// $Log$
// Contains the Mismatch class for afin

#include "mismatch.hpp"

using namespace std;

// Mismatch Class Code
Mismatch::Mismatch(){
  id = -1;
  score = 1.0;
  length = 0;
  orientation = 0;
  match_index_i = 0;
  match_index_j = 0;
  rev = false;
}

Mismatch::Mismatch( int id, double score, int length, int orientation, int match_index_i, int match_index_j ): id(id), score(score), length(length), orientation(orientation), match_index_i(match_index_i), match_index_j(match_index_j){
  rev = false;
}

Mismatch::Mismatch( int id, int orientation, int match_index_i, int match_index_j ): id(id), orientation(orientation), match_index_i(match_index_i), match_index_j(match_index_j){
  score = 1.0;
  length = 0;
  rev = false;
}

// set id
void Mismatch::set_id( int id ){
  this->id = id;
}

// set mismatch score
void Mismatch::set_score( double score ){
  this->score = score;
}

// set length
void Mismatch::set_length( int length ){
  this->length = length;
}

// set orientation
void Mismatch::set_orientation( int orientation ){
  this->orientation = orientation;
}

// set index
void Mismatch::set_indices( int match_index_i, int match_index_j ){
  this->match_index_i = match_index_i;
  this->match_index_j = match_index_j;
}

// set rev
void Mismatch::set_rev( bool rev ){
  this->rev = rev;
}

// return id
int Mismatch::get_id(){
  return id;
}

// return mismatch score
double Mismatch::get_score(){
  return score;
}

// return length
int Mismatch::get_length(){
  return length;
}

// return orientation
int Mismatch::get_orientation(){
  return orientation;
}

// return index i
int Mismatch::get_index_i(){
  return match_index_i;
}

// return index j
int Mismatch::get_index_j(){
  return match_index_j;
}

// return rev
bool Mismatch::get_rev(){
  return rev;
}
