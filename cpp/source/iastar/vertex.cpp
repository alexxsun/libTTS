#include "vertex.h"

#include <math.h>
#include <stdlib.h>

using namespace std;

Vertex::Vertex() {}

Vertex::Vertex(vector<float> coords) { coordinates = coords; }

void Vertex::changeCoordinate(float coord, int index) { coordinates[index] = coord; }

vector<float>& Vertex::getCoordinates() { return coordinates; }

void Vertex::addCoordinate(int i, float f) {
  this->coordinates[i + 3] = f;  // todo: is it OK? the size of coordinates?
}

void Vertex::addLastCoor(float f) { this->coordinates.push_back(f); }

float Vertex::getCoordinate(int index) const { return coordinates[index]; }

int Vertex::getCoboundaryMaxDim() {
  int max = 0;
  for (map<int, vector<int> >::iterator it = partial_coboundaryTop.begin(); it != partial_coboundaryTop.end(); it++)
    if (max < it->first) max = it->first;
  return max;
}

int Vertex::getCoboundaryTopNum() {
  int sum = 0;
  for (map<int, vector<int> >::iterator it = partial_coboundaryTop.begin(); it != partial_coboundaryTop.end(); it++)
    sum += it->second.size();
  return sum;
}

int Vertex::getCoboundaryTopNum(int dim) {
  map<int, vector<int> >::iterator it = partial_coboundaryTop.find(dim);
  if (it != partial_coboundaryTop.end())
    return it->second.size();
  else
    return 0;
}

void Vertex::addPartialCoboundaryTop(int dim, int index) { partial_coboundaryTop[dim].push_back(index); }

map<int, vector<int> >& Vertex::getPartialCoboundaryTopRelations() { return partial_coboundaryTop; }

vector<int>* Vertex::getPartialCoboundaryTop(int dim) {
  map<int, vector<int> >::iterator it = partial_coboundaryTop.find(dim);
  if (it != partial_coboundaryTop.end())
    return &(it->second);
  else
    return NULL;
}

float Vertex::euclideanDistance(Vertex& v) {
  if (this->coordinates.size() > v.getCoordinates().size()) {
    float dist = 0;
    for (uint i = 0; i < coordinates.size(); i++) {
      if (i >= v.getCoordinates().size())
        dist += pow(coordinates[i], 2);
      else
        dist += pow(coordinates[i] - v.getCoordinate(i), 2);
    }

    return sqrt(dist);
  } else {
    float dist = 0;
    for (uint i = 0; i < v.getCoordinates().size(); i++) {
      if (i >= coordinates.size())
        dist += pow(v.getCoordinate(i), 2);
      else
        dist += pow(coordinates[i] - v.getCoordinate(i), 2);
    }

    return sqrt(dist);
  }
}

float Vertex::euclideanDistance_xy(Vertex& v) {
  double del_x = this->getCoordinate(0) - v.getCoordinate(0);
  double del_y = this->getCoordinate(1) - v.getCoordinate(1);
  return sqrt(del_x * del_x + del_y * del_y);
}

void Vertex::middlePoint(Vertex& v) {
  if (coordinates.size() == v.getCoordinates().size()) {
    for (int i = 0; i < coordinates.size(); i++) {
      coordinates[i] = (coordinates[i] + v.getCoordinate(i)) / 2.0;
    }
  }
}

bool Vertex::operator<(const Vertex& v) const {
  for (int i = 0; i < coordinates.size(); i++) {
    if (coordinates[i] != v.getCoordinate(i)) return coordinates[i] < v.getCoordinate(i);
  }

  return true;
}
