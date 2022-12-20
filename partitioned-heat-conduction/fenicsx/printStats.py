import pstats
from pstats import SortKey
p = pstats.Stats('Dirichlet.prof')
p.strip_dirs().sort_stats(1).print_stats(20)
