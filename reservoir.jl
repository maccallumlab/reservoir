import Base.show
import StatsBase: Weights, sample


"""
    Walker(landscape, pos)

Immutable structure to hold the state of a 1-D Walker.
"""
mutable struct Walker
    landscape::Array{Float64, 1}
    pos::Int
end

Base.show(io::IO, x::Walker) = print(io, "Walker with position ", x.pos)


"""
    walk(walker, steps)

Take `steps` with a 1-D walker.
"""
function walk(walker::Walker, steps::Int)
    current = walker.pos
    @assert walker.pos > 0 && walker.pos <= size(walker.landscape, 1)
    deltas = rand(-1:2:1, steps)
    for step = 1:steps
        delta = deltas[step]
        proposed = current + delta
        if proposed < 1 || proposed > size(walker.landscape, 1)
            continue
        end
        e1 = walker.landscape[current]
        e2 = walker.landscape[proposed]
        accepted = false
        if e2 < e1
            accepted = true
        else
            metrop = exp(e1 - e2)
            if rand() < metrop
                accepted = true
            end
        end
        if accepted
            current = proposed
        end
    end
    walker.pos = current
    nothing
end


"""
    ReservoirEntry(index, weight)

Holds the index to a structure with an associated weight.
"""
struct ReservoirEntry
    index::Int
    weight::Float64
end


"""
    Reservoir(targetsize, contents)

A reservoir with target size `targetsize`.
"""
mutable struct Reservoir
    targetsize::Int
    lowwater::Int
    highwater::Int
    maxratio::Float64
    contents::Array{ReservoirEntry, 1}
    totalweight::Float64
    nsamples::Int
end

Base.show(io::IO, r::Reservoir) = print(io, "Reservoir with ", r.targetsize, " entries")

Reservoir(N, low, high, maxratio, contents) = Reservoir(N, low, high, maxratio, contents, 0.0, 0)
Reservoir(N, low, high, maxratio) = Reservoir(N, low, high, maxratio, [])


"""
    insert!(reservoir, entry)

Insert `entry` into `reservoir`.
"""
function insert!(reservoir::Reservoir, entry::ReservoirEntry)
    reservoir.totalweight += entry.weight
    reservoir.nsamples += 1
    meanweight = reservoir.totalweight / reservoir.nsamples
    target = entry.weight / meanweight
    nmin = Int(floor(target))
    if rand() > (target - nmin)
        count = nmin
    else
        count = nmin + 1
    end
    for i=1:count
        push!(reservoir.contents, ReservoirEntry(entry.index, meanweight))
    end
    resample!(reservoir)
end


"""
    poprandom!(reservoir)

Remove a random entry from `reservoir` and return it.
"""
function poprandom!(reservoir::Reservoir)::ReservoirEntry
    weights = Weights([x.weight for x in reservoir.contents])
    index = sample(weights)
    entry = reservoir.contents[index]
    deleteat!(reservoir.contents, index)
    resample!(reservoir)
    entry
end


"""
    resample!(reservoir)

Resample to obtain exactly `reservoir.targetsize` entries.
"""
function resample!(res::Reservoir)
    weights = [e.weight for e in res.contents]
    
    # check if we need to resample
    needsresample = false
    if size(res.contents, 1) < res.lowwater
        needsresample = true
    end
    if size(res.contents, 1) > res.highwater
        needsresample = true
    end
    if maximum(weights) / minimum(weights) > res.maxratio
        needsresample = true
    end
    
    if needsresample
        total_weight = sum(weights)
        new_weight = total_weight / res.targetsize
        sampled_indices = systematicresample(weights, res.targetsize)
        res.contents = [ReservoirEntry(res.contents[i].index, new_weight) for i in sampled_indices]
    end
end


"""
    systematicresample(weights, n)

Return `n` indicies using systematic resampling with `weights`.
"""
function systematicresample(weights::Array{Float64, 1}, n::Int)::Array{Int, 1}
    probs = weights / sum(weights)
    output = zeros(Int, n)
    u = rand() / n
    j = 1
    sumq = probs[j]
    
    for i=1:n
        while sumq < u
            j = j + 1
            sumq = sumq + probs[j]
        end
        output[i] = j
        u = u + 1.0 / n
    end
    output
end


struct Waterfall
    walkers::Array{Walker, 1}
    reservoirs::Array{Reservoir, 1}
    history::Array{Int, 2}
    weighthistory::Array{Float64, 2}

    function Waterfall(walkers, reservoirs, history, weighthistory)
        n = size(walkers, 1)
        @assert size(reservoirs, 1) == n
        @assert size(history, 1) == n
        @assert size(weighthistory, 1) == n
        new(walkers, reservoirs, history, weighthistory)
    end
end

Base.show(io::IO, wf::Waterfall) = print(io, "Waterfall with ", size(wf.walkers, 1), " walkers")

function runwaterfall(walkers::Array{Walker, 1}, steps::Int, walksteps::Int,
    targetsize::Int, lowwater::Int, highwater::Int, maxratio::Float64, usetopres::Bool)
    N = size(walkers, 1)
    
    # The top reservoir is different, because we want all structures
    # to have a weight of 1.0, so we need to add them all at the start
    reservoirs = [Reservoir(targetsize, lowwater, highwater, maxratio) for i in 1:N-1]
    topres = Reservoir(targetsize, lowwater, highwater, maxratio, [ReservoirEntry(1, 1.0) for i=1:targetsize])
    push!(reservoirs, topres)
    
    # Create the history
    history = zeros(Int, N, steps)
    weighthistory = zeros(Float64, N, steps)
    
    # create the waterfall
    wf = Waterfall(walkers, reservoirs, history, weighthistory)
    
    # Setup initial state of history.
    # We pre-allocate the arrays to the right
    # size for speed.
    for i=1:N
        wf.history[i, 1] = wf.walkers[i].pos
        if i==N
            wf.weighthistory[i, 1] = 1.0
        else
            wf.weighthistory[i, 1] = 1e-99
        end
    end
    
    # Setup the initial state of the reservoir.
    # The top reservoir has already been handled
    for i=1:N-1
        insert!(wf.reservoirs[i], ReservoirEntry(1, 1e-99))
    end
    
    #
    # Do all of our steps
    #
    for step=2:steps
        #
        # update the top walker
        #
        if !usetopres
            walk(wf.walkers[N], walksteps)
            w = wf.walkers[N]
            wf.history[N, step] = w.pos
            wf.weighthistory[N, step] = 1.0
            neww = computew(wf.walkers[N], wf.walkers[N-1], w.pos)
            insert!(wf.reservoirs[N-1], ReservoirEntry(step, neww))
        else
            entry = poprandom!(wf.reservoirs[N])
            oldpos = wf.history[N, entry.index]
            oldw = entry.weight
            wf.walkers[N].pos = oldpos
            walk(wf.walkers[N], walksteps)
            w = wf.walkers[N]
            wf.history[N, step] = w.pos
            wf.weighthistory[N, step] = oldw
            neww = computew(wf.walkers[N], wf.walkers[N-1], w.pos)
            insert!(wf.reservoirs[N-1], ReservoirEntry(step, neww))
            insert!(wf.reservoirs[N], ReservoirEntry(step, oldw))
        end
        
        #
        # Update the middle walkers
        #
        for i=N-1:-1:2
            # Get a random structure from the reservoir
            entry = poprandom!(wf.reservoirs[i])
            oldpos = wf.history[i+1, entry.index]
            oldw = entry.weight
            # Update the walker
            wf.walkers[i].pos = oldpos
            walk(wf.walkers[i], walksteps)
            w = wf.walkers[i]
            # Update the history
            wf.history[i, step] = w.pos
            # Compute the weight
            neww = computew(wf.walkers[i], wf.walkers[i-1], w.pos)
            # Put the new structure in the next lowest reservoir
            insert!(wf.reservoirs[i-1], ReservoirEntry(step, neww*oldw))
            wf.weighthistory[i, step] = neww * oldw
        end
        
        #
        # Update the bottom walker
        #
        entry = poprandom!(wf.reservoirs[1])
        oldpos = wf.history[2, entry.index]
        oldw = entry.weight
        # Update the walker
        wf.walkers[1].pos = oldpos
        walk(wf.walkers[1], walksteps)
        w = wf.walkers[1]
        # Update the history
        wf.history[1, step] = w.pos
        wf.weighthistory[1, step] = oldw
    end
    wf
end


function computew(w1::Walker, w2::Walker, x::Int)::Float64
    e1 = w1.landscape[x]
    e2 = w2.landscape[x]
    exp(e1 - e2)
end