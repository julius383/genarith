
# Genetic Algorithm for Guessing Expressions

Try and guess an expression such as `4 * 5` when given input of a number such
as `20`. I wrote an [article][1] describing this project.

## Installing

1. clone this repo
1. `cd` into directory
1. start julia repl and enter package mode with `]`
1. run `activate .` to activate project
1. install dependencies with `instantiate`
1. in shell run `sqlite3 data.db < db.sql` to create table for saving stats

## Running

Include `genarith.jl` into the repl and run the `ga` function with 
`target = somevalue`. Optionally uncomment the last few lines of `genarith.jl`
to be able to do `julia --project genarith.jl somevalue` to run the program
from the CLI.

[1]: https://julius383.github.io/posts-output/2023-11-27-using-a-genetic-algorithm-to-guess-arithmetic-equations-in-julia/
