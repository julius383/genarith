using Random;
using Term;
using Logging;
using Base: convert;
using ProgressMeter;
using StatsBase: sample, Weights;
using SQLite;
using JSON;

const GENE_SIZE = 4
const EXPRESSION_LENGTH = 9
const CHROMOSOME_SIZE = GENE_SIZE * EXPRESSION_LENGTH
# TODO: add these as default parameters to relevant functions
const CROSSOVER_RATE = Float64(0.5)
const MUTATION_RATE = Float64(0.01)

io = open("log.txt", "w+")
logger = SimpleLogger(io)
global_logger(logger)


db = SQLite.DB("data.db")

Expression = Vector{UInt64}
EncodedExpression = UInt64


expression_chars =
    ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '+', '-', '*', '/']

isnumber(x) = x >= 0 && x <= 9
isoperator(x) = x > 9 && x <= 13

"""
    encode(expression::Expression)::EncodedExpression

Convert expression from Vector{UInt64} to UInt64.

Relies of the `GENE_SIZE` constant.

See also [`decode`](@ref)
"""
function encode(expression::Expression)::EncodedExpression
    res = EncodedExpression(0)
    for i in expression
        res <<= GENE_SIZE
        res |= i
    end
    return res
end

"""
    decode(expression::EncodedExpression)::Expression

Convert expression from UInt64 to Vector{UInt64}.

Relies of the `EXPRESSION_LENGTH` and `GENE_SIZE` constants.

See also [`encode`](@ref)
"""
function decode(expression::EncodedExpression)::Expression
    expr = []
    expr_mask = (1 << GENE_SIZE) - 1
    for _ in (1:EXPRESSION_LENGTH)
        v = expression & expr_mask
        pushfirst!(expr, v)
        expression >>= GENE_SIZE
    end
    return expr
end


Base.convert(::Type{EncodedExpression}, x::Expression) = encode(x)
Base.convert(::Type{Expression}, x::EncodedExpression) = decode(x)


"""
    make(str::String)::Expression

Creates an expression from a string representation.

# Examples
```julia-repl
julia> make("1 + 2")
3-element Vector{UInt64}:
 0x0000000000000001
 0x000000000000000a
 0x0000000000000002
```
"""
function make(str::String)::Expression
    e = []
    str = first(split(str, "="))
    for i in str
        if !(i in expression_chars)
            continue
        else
            push!(e, findfirst(isequal(i), expression_chars) - 1)
        end
    end
    return clean(convert(Expression, e))
end

"""
    isvalid(expression::Expression)::Bool

Check if expression is well-formed.

A well-formed expression follows the general pattern [0-9]([+-*/][0-9])*.
"""
function isvalid(expression::Expression)
    nums = [expression[i] for i in (1:2:length(expression))]
    ops = [expression[i] for i in (2:2:length(expression))]
    return all(isnumber, nums) && all(isoperator, ops)
end

function isvalid(expression::EncodedExpression)
    return isvalid(convert(Expression, expression))
end

"""
    randomexpression(; num_terms = nothing)::Expression

Generates a random expression of variable length.

`num_terms` is the number of terms the expression should have
and must be an odd number less than `EXPRESSION_LENGTH`.
"""
function random_expression(; num_terms = nothing)::Expression
    # TODO: try and include logic to generate expression with a minimum and maximum value
    num_terms = isnothing(num_terms) ? rand(1:2:EXPRESSION_LENGTH) : num_terms
    num_terms % 2 == 0 &&
        num_terms <= EXPRESSION_LENGTH &&
        error("num_terms must be odd number")
    expression = []
    for i in (1:num_terms)
        if i % 2 != 0
            push!(expression, rand(Vector(1:9)))
        else
            push!(expression, rand(Vector(10:13)))
        end
    end
    return expression
end

"""
    clean(expression::Expression)

Turns a malformed expression into one that satisfies `is_valid`.

# Examples

```julia-repl
julia> clean(make("1 - - 0 + 5 /"))
5-element Vector{UInt64}:
 0x0000000000000001
 0x000000000000000b
 0x0000000000000000
 0x000000000000000a
 0x0000000000000005

julia> show_expression(ans)
"1 - 0 + 5 = 6.0"
```

See also [`is_valid`](@ref)
"""
function clean(expression::Expression)::Expression
    if isvalid(expression) || isempty(expression)
        return expression
    end
    expression = copy(expression)
    num = true
    result = []
    for i in expression
        if isnumber(i) && num
            push!(result, i)
            num = !num
        elseif isoperator(i) && !num
            push!(result, i)
            num = !num
        end
    end
    if num
        pop!(result)
    end
    return result
end

"""
    evaluate(expression::Expression)

Computes the value of evaluating the expression in strictly left to right order.

# Examples
```julia-repl
julia> evaluate(make("4 * 9 / 2 + 5 - 3 "))
20.0

julia> evaluate(make("7 - 3 * 6 + 4 / 8"))
3.5
```
"""
function evaluate(expression::Expression)
    expression = clean(expression)
    result = Float64(popfirst!(expression))
    while !isempty(expression)
        o = popfirst!(expression)
        l = popfirst!(expression)
        if '*' == expression_chars[o+1]
            result *= l
        elseif '/' == expression_chars[o+1]
            result /= l
        elseif '+' == expression_chars[o+1]
            result += l
        elseif '-' == expression_chars[o+1]
            result -= l
        end
    end
    return result
end

function evaluate(expression::EncodedExpression)
    return evaluate(convert(Expression, expression))
end


"""
    show_expression(expression::Expression, evaluate::Bool = true)

Returns the human readable format of an expression.

`evaluate` controls whether to also include the result of evaluating
the expression.

# Examples

```julia-repl
julia> showexpression(0x00000007c6a8d3b6)
"7 * 6 + 8 / 3 - 6 = 10.67"
```
"""
function showexpression(expression::Expression, show_result::Bool = true)
    expression = clean(expression)
    if isempty(expression)
        return "empty expression"
    end
    return join([expression_chars[i+1] for i in expression], " ") * (
        show_result ? " = " * string(round(evaluate(expression), digits = 4)) :
        ""
    )
end

function showexpression(expression::EncodedExpression, show_result::Bool = true)
    return showexpression(convert(Expression, expression), show_result)
end

asbitstring(expression::EncodedExpression) =
    string(expression, base = 2, pad = CHROMOSOME_SIZE)

"""
    show_encoding(expression::EncodedExpression)

Returns the bitstring representation of an encoded expression.

Bits are grouped into `GENE_SIZE` chunks for ease of comparisson.

# Examples
```julia-repl
julia> showencoding(0x00000007c6a8d3b6)
"0111 1100 0110 1010 1000 1101 0011 1011 0110"

See also [`printencoding`](@ref)
```
"""
function showencoding(expression::EncodedExpression)
    bitstr = asbitstring(expression)
    n = length(bitstr)
    formatted_bitstr = String[]

    for i = 1:GENE_SIZE:n
        push!(formatted_bitstr, bitstr[i:min(i + (GENE_SIZE - 1), n)])
    end
    return join(formatted_bitstr, ' ')
end

printencoding(expression::EncodedExpression) =
    println(@bold showencoding(expression))

"""
    printdiff(e1::EncodedExpression, e2::EncodedExpression)

Prints the difference between the bits of two encoded expressions.

Red is used when the second is has 0 has 1 in a position and green
for vice versa.

# Examples
```julia-repl
julia> printdiff(0x00000007c6a8d3b6, 0x00000000a6c6b9d5)
0111 1100 0110 1010 1000 1101 0011 1011 0110
0000 1010 0110 1100 0110 1011 1001 1101 0101
```
"""
function print_diff(e1::EncodedExpression, e2::EncodedExpression)
    s1 = showencoding(e1)
    s2 = showencoding(e2)
    r1 = ""
    r2 = ""
    for (a, b) in zip(s1, s2)
        if a == b
            r1 *= a
            r2 *= b
        else
            r1 *= a
            r2 *= b == '0' ? (@red (string(b))) : (@green (string(b)))
        end
    end
    println(@bold r1)
    println(@bold r2)
end


createpopulation(population_size = 100)::Vector{EncodedExpression} = map(
    Base.Fix1(convert, EncodedExpression),
    [random_expression() for _ in (1:population_size)],
)

"""
    fitness_score(expression, target)

Calculates how far away the result of an expression is from the desired target.
"""
function fitness_score(expression, target)
    if isempty(expression) || iszero(expression)
        return 0
    end
    return 1 / abs(target - evaluate(expression))
end

"""
    mutate(expression::EncodedExpression, rate::Float64 = MUTATION_RATE)

Flips a singular bit in an encoded expression according to a given `rate`.
"""
function mutate(expression::EncodedExpression, rate::Float64 = MUTATION_RATE)
    if rand(Float64) <= rate
        expression = copy(expression)
        position = rand(1:CHROMOSOME_SIZE)
        mask = 1 << position
        return xor(expression, mask)
    else
        return expression
    end
end

"""
    crossover(e1::EncodedExpression, e2::EncodedExpression, rate::Float64 = CROSSOVER_RATE)

Swaps the bits in two encoded expression to the right of a random position according to a given rate.

# Examples
```julia-repl
julia> crossover(0x00000007c6a8d3b6, 0x00000000a6c6b9d5)
(0x00000007a6c6b9d5, 0x00000000c6a8d3b6)
```
"""
function crossover(
    e1::EncodedExpression,
    e2::EncodedExpression,
    rate::Float64 = CROSSOVER_RATE,
)
    if rand(Float64) <= rate
        e1 = copy(e1)
        e2 = copy(e2)
        position = rand(1:CHROMOSOME_SIZE-1)
        mask = EncodedExpression(0)
        for i = 0:position
            mask |= (EncodedExpression(1) << i)
        end
        complement = ~mask

        c1_bits = e1 & mask
        c2_bits = e2 & mask

        r1 = e1 & complement
        r2 = e2 & complement

        r1 |= c2_bits
        r2 |= c1_bits
        return (r1, r2)
    else
        return (e1, e2)
    end
end

"""Implements fitness propotionate selection."""
function roulette_wheel(fitness, population, size = 2)
    weights = fitness / sum(fitness)
    return sample(population, Weights(weights), size, replace = false)
end

function printsolution(e::Expression, target)
    @info("found solution for $(target):", showexpression(e, false))
end

function printsolution(e::EncodedExpression, target)
    @info("found solution for $(target):", showexpression(e, false))
end

function savestats(stats)
    @info("saving results to database")
    stats[:population] = JSON.json(map(showexpression, stats[:population]))
    stmt = DBInterface.prepare(
        db,
        """
        INSERT INTO experiment_stats 
                       (target, 
                       population_size, max_generations, generations_run, mutations, crossovers, solution, population)
                       VALUES (:target, :population_size, :max_generations, :generations_run, :mutations, :crossovers, :solution, :population)
        """,
    )
    res = DBInterface.execute(stmt, stats)
    DBInterface.close!(stmt)
    return res
end

function targetexists(target)
    r = DBInterface.execute(
        db,
        "select id from experiment_stats where target = ?",
        [target],
    )
    return !isempty(r)
end


"""
    ga(; target = 20.0, population_size = 1000, max_generations = 500, show_progress = true)

Runs genetic algorithm to try and find expression that evaluates to
`target`.

See also [`mutate`](@ref) [`crossover`](@ref) [`fitness_score`](@ref)"""
function ga(;
    target = 20.0,
    population_size = 1000,
    max_generations = 500,
    show_progress = true,
    save_results = true,
)
    # 4 * 9 / 2 + 5 - 3 = 20 
    # 7 - 3 * 6 + 4 / 8 = 3.5
    @info("searching for solution to $(target)")
    population = createpopulation(population_size)
    fitness = fitness_score.(population, target)
    stats = Dict{Symbol,Union{Nothing,Real,String,Vector{UInt64}}}(
        :target => target,
        :population_size => length(population),
        :max_generations => max_generations,
        :generations_run => 1,
        :mutations => 0,
        :crossovers => 0,
        :solution => nothing,
        :population => nothing,
    )
    prog = ProgressUnknown(
        "Searching for solution to $(target):",
        spinner = true,
        enabled = show_progress,
    )
    if Inf in fitness
        e = population[findfirst(isinf, fitness)]
        ProgressMeter.finish!(prog)
        printsolution(e, target)
        stats[:solution] = showexpression(e, false)
        stats[:population] = population
        save_results && savestats(stats)
        return stats[:solution]
    end

    for g in (1:max_generations)
        for _ in (1:population_size)
            ProgressMeter.next!(
                prog,
                showvalues = [
                    (:generations_run, stats[:generations_run]),
                    (:mutations, stats[:mutations]),
                    (:crossovers, stats[:crossovers]),
                ],
            )
            p1, p2 = roulette_wheel(fitness, population, 2)

            e1, e2 = crossover(p1, p2)
            e1, e2 = mutate(e1), mutate(e2)

            f1, f2 = fitness_score(e1, target), fitness_score(e2, target)
            if isnan(f1) || isnan(f2)
                continue
            end
            stats[:crossovers] += e1 == p1 ? 0 : 1
            stats[:mutations] += e1 == p1 ? 0 : 1
            stats[:mutations] += e2 == p2 ? 0 : 1

            for e in [e1, e2]
                if isinf(fitness_score(e, target))
                    ProgressMeter.finish!(prog)
                    printsolution(e, target)
                    stats[:solution] = showexpression(e, false)
                    stats[:population] = population
                    save_results && savestats(stats)
                    return stats[:solution]
                end
            end
            least_fit = argmin(fitness)
            population[least_fit] = f1 < f2 ? e1 : e2
            fitness[least_fit] = f1 < f2 ? f1 : f2
        end
        stats[:generations_run] = g
    end
    stats[:population] = population
    save_results && savestats(stats)
    @info(
        "unable to find solution to $(target) within $(max_generations) and population of $(population_size)"
    )
end

function run()
    vals = rand(1000:2000, (200, 1))
    for v in vals
        if targetexists(v)
            continue
        end
        ga(target = convert(Float64, v))
        flush(io)
    end
    DBInterface.close!(db)
end

run()

# if !isempty(ARGS)
#     ga(target = parse(Float64, ARGS[1]))
# end
