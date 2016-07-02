#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

struct Grid {
  float x;
  float y;
};

typedef std::pair<int, Grid> city;
bool sort_x(const city& left, const city& right)
{
  return left.second.x < right.second.x ;
}

bool sort_y(const city& left, const city& right)
{
  return left.second.y < right.second.y ;
}

std::map<int, Grid> read_input(std::string filename) {
  std::ifstream ifs(filename);
  std::string x, y;
  std::map<int, Grid> cities;
  int index = 0;
  getline(ifs, x);
  while (1) {
    getline(ifs, x, ',');
    getline(ifs, y);
    if (!ifs)
      break;
    cities[index] = {std::stof(x), std::stof(y)};
    index++;    
  }
  return cities;
}

float distance(Grid city1, Grid city2) {
  return sqrt(pow((city1.x - city2.x), 2) + pow((city1.y - city2.y), 2));
}

std::vector<std::vector<float>> calculateEachDistance(
  std::map<int, Grid> cities) {
  std::vector<std::vector<float>> dist(
    cities.size(), std::vector<float> (cities.size()));
  for (std::string::size_type i = 0; i < cities.size(); i++) {
     for (std::string::size_type j = 0; j < cities.size(); j++) {
       dist[i][j] = dist[j][i] = distance(cities[i], cities[j]);
     }
  }
  return dist;
}

float slope(Grid city1, Grid city2) {
  return (city1.y - city2.y) / (city1.x - city2.x);
}

std::set<int> makeCircle(
  Grid xmin, Grid xmax,Grid ymin, Grid ymax,
  std::map<int, Grid>* cities) {
  std::set<int> outCircle;
  for (std::string::size_type i = 0; i < cities->size(); i++) {
    Grid city = cities->at(i);
    if ((city.y >  slope(xmin, ymax) * (city.x - xmin.x) + xmin.y)
        || (city.y >  slope(xmax, ymax) * (city.x - xmax.x) + xmax.y)
        || (city.y <  slope(xmin, ymin) * (city.x - xmin.x) + xmin.y) 
        || (city.y <  slope(xmax, ymin) * (city.x - xmin.x) + xmin.y)) {
      outCircle.insert(i);
      cities->erase(i);
    }
  }
  return outCircle;
}

int nearestCity(std::vector<std::vector<float>> dist,
                std::set<int> unvisitedCities, int currentCity) {
  int minCity;
  float min = 1000000;
  for (int i : unvisitedCities) {
    if (dist[currentCity][i] < min) {
      minCity = i;
      min = dist[currentCity][i];
    }
  }
  return minCity;
}

std::vector<int> solve(std::map<int, Grid> cities) {
  int nextCity, currentCity = 0;
  std::set<int> unvisitedCities, visitSoonCities;
  std::map<int, Grid> xSort = cities;
  std::map<int, Grid> ySort = cities;
  std::vector<int> solution{0};
  float sum = 0;
  int before = 0;
  sort(xSort.begin(), xSort.end(), sort_x);
  sort(ySort.begin(), ySort.end(), sort_y)
  for (std::string::size_type i = 1; i < cities.size(); i++) {
    unvisitedCities.insert(i);
  }
  std::vector<std::vector<float>> dist = calculateEachDistance(cities);
  while (unvisitedCities.size() > 0) {
    nextCity = nearestCity(dist, unvisitedCities, currentCity);
    unvisitedCities.erase(nextCity);
    solution.push_back(nextCity);
    currentCity = nextCity;
  }
   for (auto city : solution) {
    sum += dist[before][city];
    before = city;
  }
   sum += dist[0][solution.back()];
   std::cout << sum << std::endl;
  return solution;
}

int main(int argc, char *argv[]) {
  if (argc != 2) {
    std::cout << "Illegal input." << std::endl;
    return -1;
  }
  std::map<int, Grid> cities = read_input(argv[1]);
  std::vector<int> solution = solve(cities);
  printf("index\n");
  for (int city : solution) {
    printf("%d\n", city);
  }
  return 0;
}
