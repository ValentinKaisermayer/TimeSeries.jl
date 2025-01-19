"""
    +(ta::TimeArray)

Unary plus operator for `TimeArray`.
"""
function Base.:(+)(ta::TimeArray)
    return .+ta
end

"""
    -(ta::TimeArray)

Unary minus operator for `TimeArray`.
"""
function Base.:(-)(ta::TimeArray)
    return .-ta
end

###### lag, lead ################

"""
    lag(ta::TimeArray{T, N}, n::Int=1; padding::Bool=false, period::Int=0) where {T, N}

Lag the time series by `n` periods. Optionally pad with `NaN` values.

# Arguments
- `ta::TimeArray{T, N}`: The input time series.
- `n::Int=1`: Number of periods to lag.
- `padding::Bool=false`: Whether to pad with `NaN` values.
- `period::Int=0`: Deprecated, use `n` instead.

# Returns
- `TimeArray`: The lagged time series.
"""
function lag(ta::TimeArray{T,N}, n::Int=1; padding::Bool=false, period::Int=0) where {T,N}
    if period != 0
        @warn("the period kwarg is deprecated, use lag(ta::TimeArray, period::Int) instead")
        n = period
    end

    # TODO: apply `unchecked`
    if padding
        paddedvals = [NaN * ones(n, size(ta, 2)); values(ta)[1:end-n, :]]
        ta = TimeArray(timestamp(ta), paddedvals, colnames(ta), meta(ta))
    else
        ta = TimeArray(timestamp(ta)[1+n:end], values(ta)[1:end-n, :], colnames(ta), meta(ta))
    end

    N == 1 && (ta = ta[colnames(ta)[1]])
    return ta
end  # lag

"""
    lead(ta::TimeArray{T, N}, n::Int=1; padding::Bool=false, period::Int=0) where {T, N}

Lead the time series by `n` periods. Optionally pad with `NaN` values.

# Arguments
- `ta::TimeArray{T, N}`: The input time series.
- `n::Int=1`: Number of periods to lead.
- `padding::Bool=false`: Whether to pad with `NaN` values.
- `period::Int=0`: Deprecated, use `n` instead.

# Returns
- `TimeArray`: The lead time series.
"""
function lead(ta::TimeArray{T,N}, n::Int=1; padding::Bool=false, period::Int=0) where {T,N}
    if period != 0
        @warn("the period kwarg is deprecated, use lead(ta::TimeArray, period::Int) instead")
        n = period
    end

    if padding
        paddedvals = [values(ta)[1+n:end, :]; NaN * ones(n, length(colnames(ta)))]
        ta = TimeArray(timestamp(ta), paddedvals, colnames(ta), meta(ta))
    else
        ta = TimeArray(timestamp(ta)[1:end-n], values(ta)[1+n:end, :], colnames(ta), meta(ta))
    end

    N == 1 && (ta = ta[colnames(ta)[1]])
    return ta
end  # lead

###### diff #####################

"""
    diff(ta::TimeArray, n::Int=1; padding::Bool=false, differences::Int=1)

Compute the difference of the time series.

# Arguments
- `ta::TimeArray`: The input time series.
- `n::Int=1`: Number of periods to lag for the difference.
- `padding::Bool=false`: Whether to pad with `NaN` values.
- `differences::Int=1`: Number of differences to compute.

# Returns
- `TimeArray`: The differenced time series.
"""
function Base.diff(ta::TimeArray, n::Int=1; padding::Bool=false, differences::Int=1)
    cols = colnames(ta)
    for d in 1:differences
        ta = ta .- lag(ta, n, padding=padding)
    end
    colnames(ta)[:] = cols
    return ta
end  # diff

###### percentchange ############

"""
    percentchange(ta::TimeArray, returns::Symbol=:simple; padding::Bool=false, method::AbstractString="")

Compute the percent change of the time series.

# Arguments
- `ta::TimeArray`: The input time series.
- `returns::Symbol=:simple`: Type of returns, either `:simple` or `:log`.
- `padding::Bool=false`: Whether to pad with `NaN` values.
- `method::AbstractString=""`: Deprecated, use `returns` instead.

# Returns
- `TimeArray`: The percent change time series.
"""
function percentchange(ta::TimeArray, returns::Symbol=:simple; padding::Bool=false, method::AbstractString="")
    if method != ""
        @warn("the method kwarg is deprecated, use percentchange(ta, :methodname) instead")
        returns = Symbol(method)
    end

    cols = colnames(ta)
    ta = returns == :log ? diff(log.(ta), padding=padding) :
         returns == :simple ? expm1.(percentchange(ta, :log, padding=padding)) :
         throw(ArgumentError("returns must be either :simple or :log"))
    colnames(ta)[:] = cols

    return ta
end  # percentchange

###### moving ###################

"""
    moving(f, ta::TimeArray{T,1}, window::Integer; padding::Bool = false) where {T}

Apply user-defined function `f` to a 1D `TimeArray` with window size `w`.

# Arguments
- `f`: The function to apply.
- `ta::TimeArray{T,1}`: The input time series.
- `window::Integer`: The window size.
- `padding::Bool=false`: Whether to pad with `NaN` values.

# Returns
- `TimeArray`: The resulting time series after applying the function.
"""
function moving(f, ta::TimeArray{T,1}, window::Integer; padding::Bool=false) where {T}
    ts = padding ? timestamp(ta) : @view(timestamp(ta)[window:end])
    A = values(ta)
    vals = map(i -> f(@view A[i:i+(window-1)]), 1:length(ta)-window+1)
    padding && (vals = [fill(NaN, window - 1); vals])
    TimeArray(ta; timestamp=ts, values=vals)
end

"""
    moving(f, ta::TimeArray{T,2}, window::Integer; padding::Bool = false, dims::Integer = 1, colnames::AbstractVector{Symbol} = _colnames(ta)) where {T}

Apply user-defined function `f` to a 2D `TimeArray` with window size `w`.

# Arguments
- `f`: The function to apply.
- `ta::TimeArray{T,2}`: The input time series.
- `window::Integer`: The window size.
- `padding::Bool=false`: Whether to pad with `NaN` values.
- `dims::Integer=1`: The dimension along which to apply the function.
- `colnames::AbstractVector{Symbol}=_colnames(ta)`: Column names for the resulting time series.

# Returns
- `TimeArray`: The resulting time series after applying the function.
"""
function moving(f, ta::TimeArray{T,2}, window::Integer; padding::Bool=false, dims::Integer=1, colnames::AbstractVector{Symbol}=_colnames(ta)) where {T}
    if !(dims ∈ (1, 2))
        throw(ArgumentError("invalid dims $dims"))
    end

    ts = padding ? timestamp(ta) : @view(timestamp(ta)[window:end])
    A = values(ta)

    if dims == 1
        vals = similar(@view(A[window:end, :]))
        for i ∈ 1:size(vals, 1), j ∈ 1:size(vals, 2)
            vals[i, j] = f(@view(A[i:i+(window-1), j]))
        end
    else # case of dims = 2
        vals = mapreduce(i -> f(view(A, i-window+1:i, :)), vcat, window:size(A, 1))
        if size(vals, 2) != length(colnames)
            throw(DimensionMismatch(
                "the output dims should match the length of columns, " *
                "please set the keyword argument `colnames` properly."
            ))
        end
    end

    padding && (vals = [fill(NaN, (window - 1), size(vals, 2)); vals])
    TimeArray(ta; timestamp=ts, values=vals, colnames=colnames)
end

###### upto #####################

"""
    upto(f, ta::TimeArray{T, 1}) where {T}

Apply user-defined function `f` cumulatively to a 1D `TimeArray`.

# Arguments
- `f`: The function to apply.
- `ta::TimeArray{T, 1}`: The input time series.

# Returns
- `TimeArray`: The resulting time series after applying the function.
"""
function upto(f, ta::TimeArray{T,1}) where {T}
    vals = zero(values(ta))
    for i = 1:length(vals)
        vals[i] = f(values(ta)[1:i])
    end
    TimeArray(timestamp(ta), vals, colnames(ta), meta(ta))
end

"""
    upto(f, ta::TimeArray{T, 2}) where {T}

Apply user-defined function `f` cumulatively to a 2D `TimeArray`.

# Arguments
- `f`: The function to apply.
- `ta::TimeArray{T, 2}`: The input time series.

# Returns
- `TimeArray`: The resulting time series after applying the function.
"""
function upto(f, ta::TimeArray{T,2}) where {T}
    vals = zero(values(ta))
    for i = 1:size(vals, 1), j = 1:size(vals, 2)
        vals[i, j] = f(values(ta)[1:i, j])
    end
    TimeArray(timestamp(ta), vals, colnames(ta), meta(ta))
end

###### basecall #################

"""
    basecall(ta::TimeArray, f::Function; cnames=colnames(ta))

Apply a function `f` to the values of the `TimeArray`.

# Arguments
- `ta::TimeArray`: The input time series.
- `f::Function`: The function to apply.
- `cnames`: Column names for the resulting time series.

# Returns
- `TimeArray`: The resulting time series after applying the function.
"""
function basecall(ta::TimeArray, f::Function; cnames=colnames(ta))
    return TimeArray(timestamp(ta), f(values(ta)), cnames, meta(ta))
end

###### uniform observations #####

"""
    uniformspaced(ta::TimeArray)

Check if the time series has uniform spacing.

# Arguments
- `ta::TimeArray`: The input time series.

# Returns
- `Bool`: `true` if the time series is uniformly spaced, `false` otherwise.
"""
function uniformspaced(ta::TimeArray)
    gap1 = timestamp(ta)[2] - timestamp(ta)[1]
    i, n, is_uniform = 2, length(ta), true
    while is_uniform & (i < n)
        is_uniform = gap1 == (timestamp(ta)[i+1] - timestamp(ta)[i])
        i += 1
    end
    return is_uniform
end  # uniformspaced

"""
    uniformspace(ta::TimeArray{T, N}) where {T, N}

Create a uniformly spaced time series.

# Arguments
- `ta::TimeArray{T, N}`: The input time series.

# Returns
- `TimeArray`: The uniformly spaced time series.
"""
function uniformspace(ta::TimeArray{T,N}) where {T,N}
    min_gap = minimum(timestamp(ta)[2:end] - timestamp(ta)[1:end-1])
    newtimestamp = timestamp(ta)[1]:min_gap:timestamp(ta)[end]
    emptyta = TimeArray(collect(newtimestamp), zeros(length(newtimestamp), 0), Symbol[], meta(ta))
    ta = merge(emptyta, ta, method=:left)
    N == 1 && (ta = ta[colnames(ta)[1]])
    return ta
end  # uniformspace

###### dropnan ####################

"""
    dropnan(ta::TimeArray, method::Symbol = :all)

Drop `NaN` values from the time series.

# Arguments
- `ta::TimeArray`: The input time series.
- `method::Symbol=:all`: Method to drop `NaN` values, either `:all` or `:any`.

# Returns
- `TimeArray`: The time series with `NaN` values dropped.
"""
function dropnan(ta::TimeArray, method::Symbol=:all)
    return method == :all ? ta[findall(reshape(values(any(.!isnan.(ta), dims=2)), :))] :
           method == :any ? ta[findall(reshape(values(all(.!isnan.(ta), dims=2)), :))] :
           throw(ArgumentError("dropnan method must be :all or :any"))
end
