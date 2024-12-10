
| Day             | Exercise              | Code type   | Source code | Exercises in jupyter | Pretty plots | Manual updated | Code reviewed |
| --------------- | --------------------- | ----------- | ----------- | -------------------- | ------------ | -------------- | ------------- |
| 1 (StatMech)    | Distributions         | python      | ✅            | ✅                     |              |                |               |
|                 | Random walk (1D & 2D) | numba       | ✅           |                      |              |                |               |
| 2 (Basic MC)    | Pi                    | python      | ✅           | ✅                     | ✅            |                |               |
|                 | Hard Disks            | numba & c++ | ✅           | ✅                     |           |                |               |
|                 | LJ - MC               | c++         | ✅           |                      |              |                |               |
| 3 (Basic MD)    | Muller Brown          | numba       | ✅           |                      |              |                |               |
|                 | LJ - MD               | c++         | ✅           | ✅                    |              |                |               |
|                 | LJ code optimization  | c++         | ✅           | ✅                    |              |                |               |
| 4 (Advanced MC) | Parallel Tempering    | numba       | ✅           |                      | ✅            |                |               |
|                 | $\mu VT$ LJ       | c++         |            |                      |              |                |               |
|                 | NPT LJ                | c++         | ✅           |                      |              |                |               |
| 5 (Advanced MD) | Umbrella sampling     | numba       |             |                      |              |                |               |
|                 | Nosé-Hoover           | c++         |  ✅        |                      |              |                |               |
|                 | Ewald (opt.)          | numba       | ✅           |                      |              |                |               |
Other:
- Manual how to use code



Errors or unknowns:
- LJ steps run 81604378640 for initial step
- LJ movie.pdb loading issue using mdtraj (loading ok in vmd)
