import sys
import math

from common import print_solution, read_input
from chris import compute
# christofides = __import__("python-christofides")
def distance(city1, city2):
    return math.sqrt((city1[0] - city2[0]) ** 2 + (city1[1] - city2[1]) ** 2)

def make_distanceMatrix(cities):
    N = len(cities)
    
    dist = [[0] * N for i in range(N)]
    for i in range(N):
        for j in range(N):
            dist[i][j] = dist[j][i] = int(distance(cities[i], cities[j]))
    return dist

if __name__ == '__main__':
    assert len(sys.argv) > 1
    dist = make_distanceMatrix(read_input(sys.argv[1]))
    TSP = compute(dist)
    solution = TSP['indexes']
    print_solution(solution)
