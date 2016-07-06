#include "TSP.h"
#include <iostream>
#include <limits>
#include <assert.h>

using namespace std;

int main(int argc, const char *argv[])
{
  /* 90% mutation probability, 2% mutation probability */
  if (argc < 2) {
    cout << "no input file" << endl;
    return -1;
  }
  TSP *tsp = new TSP(argv[1], 0.9, 0.02);
  size_t generations = 0, generationsWithoutImprovement = 0;
  double bestFitness = -1;
  while(generationsWithoutImprovement < 10)
  {
    std::cout << generations << std::endl;
    tsp->nextPopulation();
    ++generations;
    double newFitness = tsp->getBestFitness();
    if(newFitness > bestFitness)
    {
      bestFitness = newFitness;
      generationsWithoutImprovement = 0;
    }
    else
    {
      ++generationsWithoutImprovement;
    }
  }
  cout << "index" << endl;
  // for (string::size_type i = 0; i < tsp->getBestPathString().size(); i++)
  cout << "\t-Path: " << tsp->getBestPathString() << endl;
  delete tsp;
  return 0;
}
