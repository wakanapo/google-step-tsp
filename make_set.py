import sys

from common import read_input, print_solution
from optimize import optimize
from opt_2 import opt_2
from greedy import greedy
from or_opt import or_opt
from vortex import vortex
from kruskal import kruskal_greedy
from my_random import solve
from combine import optimize2, read_path

if __name__ == '__main__':
    assert len(sys.argv) > 1
    optimize = optimize(read_input(sys.argv[1]))
    optimize2 = optimize2(read_input(sys.argv[1]), read_path(sys.argv[2]))
    greedy = greedy(read_input(sys.argv[1]))
    opt_2 = opt_2(read_input(sys.argv[1]))
   
    print_solution(or_opt(read_input(sys.argv[1])))
    print_solution(vortex(read_input(sys.argv[1])))
    print_solution(kruskal_greedy(read_input(sys.argv[1])))
    print_solution(greedy)
    print_solution(optimize)
    print_solution(optimize)
    print_solution(optimize)
    print_solution(optimize)
    print_solution(optimize2)
    print_solution(optimize2)
    print_solution(optimize2)
    print_solution(optimize2)
    print_solution(optimize2)
    print_solution(opt_2)
    print_solution(opt_2)
    print_solution(solve(read_input(sys.argv[1])))
    print_solution(solve(read_input(sys.argv[1])))
    print_solution(solve(read_input(sys.argv[1])))
    print_solution(solve(read_input(sys.argv[1])))
    print_solution(solve(read_input(sys.argv[1])))
    print_solution(solve(read_input(sys.argv[1])))
   
