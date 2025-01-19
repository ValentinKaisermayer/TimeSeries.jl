###### type definition ##########

abstract type AbstractTimeSeries{T,N,D} end

"""
    TimeArray{T,N,D<:TimeType,A<:AbstractArray{T,N}} <: AbstractTimeSeries{T,N,D}

# Constructors

    TimeArray(timestamp, values[, colnames, meta = nothing])
    TimeArray(ta::TimeArray; timestamp, values, colnames, meta)
    TimeArray(data::NamedTuple; timestamp, meta = nothing)
    TimeArray(table; timestamp::Symbol, timeparser::Callable = identity)

The second constructor yields a new `TimeArray` with the new given fields.
Note that the unchanged fields will be shared, there aren't any copy for the
underlying arrays.

The third constructor builds a `TimeArray` from a `NamedTuple`.

# Arguments

- `timestamp::AbstractVector{<:TimeType}`: a vector of sorted timestamps,

- `timestamp::Symbol`: the column name of the time index from the source table.
  The constructor is used for the Tables.jl package integration.

- `values::AbstractArray`: a data vector or matrix. Its number of rows
  should match the length of `timestamp`.

- `colnames::Vector{Symbol}`: the column names. Its length should match
  the column of `values`.

- `meta::Any`: a user-defined metadata.

- `timeparser::Callable`: a mapping function for converting the source time index.
  For instance, `Dates.unix2datetime` is a common case.

# Examples

    data = (datetime = [DateTime(2018, 11, 21, 12, 0), DateTime(2018, 11, 21, 13, 0)],
            col1 = [10.2, 11.2],
            col2 = [20.2, 21.2],
            col3 = [30.2, 31.2])
    ta = TimeArray(data; timestamp = :datetime, meta = "Example")

"""
struct TimeArray{T,N,D<:TimeType,A<:AbstractArray{T,N}} <: AbstractTimeSeries{T,N,D}
    timestamp::Vector{D}
    values::A
    colnames::Vector{Symbol}
    meta::Any

    function TimeArray{T,N,D,A}(
        timestamp::AbstractVector{D},
        values::A,
        colnames::Vector{Symbol},
        meta::Any;
        unchecked=false) where {T,N,D<:TimeType,A<:AbstractArray{T,N}}
        nrow = size(values, 1)
        ncol = size(values, 2)
        colnames = copy(colnames)

        unchecked && return new(timestamp, values, replace_dupes!(colnames), meta)

        nrow != length(timestamp) && throw(DimensionMismatch("values must match length of timestamp"))
        ncol != length(colnames) && throw(DimensionMismatch("column names must match width of array"))

        issorted(timestamp) && return new(
            timestamp, values, replace_dupes!(colnames), meta)

        timestamp_r = reverse(timestamp)
        issorted(timestamp_r) && return new(
            timestamp_r, reverse(values, dims=1), replace_dupes!(colnames), meta)

        throw(ArgumentError("timestamps must be monotonic"))
    end
end

###### outer constructor ########

TimeArray(d::AbstractVector{D}, v::AbstractArray{T,N},
    c::Vector{Symbol}=gen_colnames(size(v, 2)),
    m::Any=nothing; args...) where {T,N,D<:TimeType} =
    TimeArray{T,N,D,typeof(v)}(d, v, c, m; args...)

TimeArray(d::D, v::AbstractArray{T,N},
    c::Vector{Symbol}=gen_colnames(size(v, 2)),
    m::Any=nothing; args...) where {T,N,D<:TimeType} =
    TimeArray{T,N,D,typeof(v)}([d], v, c, m; args...)

TimeArray(ta::TimeArray;
    timestamp=_timestamp(ta), values=_values(ta),
    colnames=_colnames(ta), meta=_meta(ta), args...) =
    TimeArray(timestamp, values, colnames, meta; args...)

function TimeArray(data::NamedTuple; timestamp::Symbol, meta=nothing, args...)
    columns = (key for key in keys(data) if key != timestamp)
    dat = hcat((data[key] for key in columns)...)
    TimeArray(data[timestamp], dat, collect(columns), meta; args...)
end

###### conversion ###############

function Base.convert(::Type{TimeArray{Float64,N}}, x::TimeArray{Bool,N}) where {N}
    return TimeArray(timestamp(x), Float64.(values(x)), colnames(x), meta(x); unchecked=true)
end

function Base.convert(x::TimeArray{Bool,N}) where {N}
    return convert(TimeArray{Float64,N}, x::TimeArray{Bool,N})
end

###### copy ###############

function Base.copy(ta::TimeArray)
    return TimeArray(timestamp(ta), values(ta), colnames(ta), meta(ta); unchecked=true)
end

###### length ###################

function Base.length(ata::AbstractTimeSeries)
    return length(timestamp(ata))
end

###### size #####################

function Base.size(ta::TimeArray)
    return size(values(ta))
end

function Base.size(ta::TimeArray, dim)
    return size(values(ta), dim)
end

###### ndims #####################

function Base.ndims(ta::AbstractTimeSeries{T,N}) where {T,N}
    return N
end

###### iteration protocol ########

@generated function Base.iterate(ta::AbstractTimeSeries{T,N}, i=1) where {T,N}
    val = (N == 1) ? :(values(ta)[i]) : :(values(ta)[i, :])

    quote
        i > length(ta) && return nothing
        ((timestamp(ta)[i], $val), i + 1)
    end
end

###### equal ####################

"""
    ==(x::TimeArray, y::TimeArray)

If `true`, all fields of `x` and `y` should be equal,
meaning that the two `TimeArray`s have the same values at the same points in time,
the same colnames and the same metadata.

Implies

```julia
x.timestamp == y.timestamp &&
x.values    == y.values    &&
x.colnames  == y.colnames  &&
x.meta      == y.meta
```
"""
function Base.:(==)(x::TimeArray{T,N}, y::TimeArray{S,M}) where {T,S,N,M}
    return false
end

function Base.:(==)(x::TimeArray{T,N}, y::TimeArray{S,N}) where {T,S,N}
    return all(f -> getfield(x, f) == getfield(y, f), fieldnames(TimeArray))
end

function Base.isequal(x::TimeArray{T,N}, y::TimeArray{S,M}) where {T,S,N,M}
    return false
end

function Base.isequal(x::TimeArray{T,N}, y::TimeArray{S,N}) where {T,S,N}
    return all(f -> isequal(getfield(x, f), getfield(y, f)), fieldnames(TimeArray))
end

function Base.hash(x::TimeArray, h::UInt)
    return sum(f -> hash(getfield(x, f), h), fieldnames(TimeArray))
end

###### eltype #####################

Base.eltype(::AbstractTimeSeries{T,1,D}) where {T,D} = Tuple{D,T}
Base.eltype(::AbstractTimeSeries{T,2,D}) where {T,D} = Tuple{D,Vector{T}}

###### show #####################

function Base.summary(io::IO, ta::TimeArray)
    return show(io, ta)
end

function Base.show(io::IO, ta::TimeArray)
    nrow = size(values(ta), 1)
    ncol = size(values(ta), 2)
    print(io, "$(nrow)×$(ncol) $(typeof(ta))")
    if nrow != 0
        print(io, " $(timestamp(ta)[1]) to $(timestamp(ta)[end])")
    else  # e.g. TimeArray(Date[], [])
        return
    end
end

function Base.show(io::IO, ::MIME"text/plain", ta::TimeArray; allrows=!get(io, :limit, false), allcols=!get(io, :limit, false))
    nrow = size(values(ta), 1)
    ncol = size(values(ta), 2)

    show(io, ta) # summary line

    nrow == 0 && return

    println(io)

    if allcols && allrows
        crop = :none
    elseif allcols
        crop = :vertical
    elseif allrows
        crop = :horizontal
    else
        crop = :both
    end

    data = hcat(timestamp(ta), values(ta))
    header = vcat("", string.(colnames(ta)))
    pretty_table(io, data;
        header=header,
        newline_at_end=false,
        reserved_display_lines=2,
        row_label_alignment=:r,
        header_alignment=:l,
        crop=crop,
        vcrop_mode=:middle,
    )
end

###### getindex #################

function Base.getindex(ta::TimeArray)
    throw(BoundsError(typeof(ta), []))
end

@propagate_inbounds function Base.getindex(ta::TimeArray, n::Integer)
    # avoid conversion to column vector
    TimeArray(timestamp(ta)[n], values(ta)[n:n, :], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,1}, n::Integer) where {T}
    TimeArray(timestamp(ta)[n], values(ta)[[n]], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray, r::UnitRange{<:Integer})
    TimeArray(timestamp(ta)[r], values(ta)[r, :], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,1}, r::UnitRange{<:Integer}) where {T}
    TimeArray(timestamp(ta)[r], values(ta)[r], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray, a::AbstractVector{<:Integer})
    TimeArray(timestamp(ta)[a], values(ta)[a, :], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,1}, a::AbstractVector{<:Integer}) where {T}
    TimeArray(timestamp(ta)[a], values(ta)[a], colnames(ta), meta(ta))
end

@propagate_inbounds function Base.getindex(ta::TimeArray, s::Symbol)
    n = findcol(ta, s)
    TimeArray(timestamp(ta), values(ta)[:, n], Symbol[s], meta(ta), unchecked=true)
end

@propagate_inbounds function Base.getindex(ta::TimeArray, ss::Symbol...)
    return getindex(ta, collect(ss))
end

@propagate_inbounds function Base.getindex(ta::TimeArray, ss::Vector{Symbol})
    TimeArray(ta; values=values(ta)[:, map(s -> findcol(ta, s), ss)], colnames=ss)
end

@propagate_inbounds function Base.getindex(ta::TimeArray, rows::Union{AbstractVector{<:Integer},Colon}, cols::AbstractVector{Symbol})
    TimeArray(
        ta;
        timestamp=timestamp(ta)[rows],
        values=values(ta)[rows, map(s -> findcol(ta, s), cols)],
        colnames=cols,
        unchecked=true)
end

@propagate_inbounds function Base.getindex(ta::TimeArray, n::Integer, cols)
    return getindex(ta, [n], cols)
end

@propagate_inbounds function Base.getindex(ta::TimeArray, rows, col::Symbol)
    return getindex(ta, rows, [col])
end

@propagate_inbounds function Base.getindex(ta::TimeArray, n::Integer, col::Symbol)
    return getindex(ta, [n], [col])
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,N,D}, d::D) where {T,N,D}
    idxs = searchsorted(timestamp(ta), d)
    length(idxs) == 1 ? ta[idxs[1]] : nothing
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,N,D}, dates::Vector{D}) where {T,N,D}
    dates = sort(dates)
    idxs, _ = overlap(timestamp(ta), dates)
    ta[idxs]
end

@propagate_inbounds function Base.getindex(ta::TimeArray{T,N,D}, r::StepRange{D}) where {T,N,D}
    return ta[collect(r)]
end

@propagate_inbounds function Base.getindex(ta::TimeArray, k::TimeArray{Bool,1})
    return ta[findwhen(k)]
end

function Base.lastindex(ta::TimeArray, d::Integer=1)
    (d == 1) ? length(timestamp(ta)) :
    (d == 2) ? length(colnames(ta)) :
    1
end

function Base.eachindex(ta::TimeArray)
    return Base.OneTo(length(timestamp(ta)))
end

###### getproperty/propertynames #################

function Base.getproperty(ta::AbstractTimeSeries, c::Symbol)
    return ta[c]
end

function Base.propertynames(ta::TimeArray)
    return colnames(ta)
end

###### element wrappers ###########

"""
    timestamp(ta::TimeArray)

Get the time index of a `TimeArray`.
"""
function timestamp(ta::TimeArray)
    return getfield(ta, :timestamp)
end

"""
    values(ta::TimeArray)

Get the underlying value table of a `TimeArray`.
"""
function values(ta::TimeArray)
    return getfield(ta, :values)
end

"""
    colnames(ta::TimeArray)

Get the column names of a `TimeArray`.

# Examples

```julia-repl
julia> colnames(ohlc)
4-element Array{Symbol,1}:
 :Open
 :High
 :Low
 :Close
```
"""
function colnames(ta::TimeArray)
    return getfield(ta, :colnames)
end

"""
    meta(ta::TimeArray)

Get the user-defined metadata of a `TimeArray`.
"""
function meta(ta::TimeArray)
    return getfield(ta, :meta)
end

# internal use, to avoid name collision
function _timestamp(ta::TimeArray)
    return getfield(ta, :timestamp)
end

function _values(ta::TimeArray)
    return getfield(ta, :values)
end

function _colnames(ta::TimeArray)
    return getfield(ta, :colnames)
end

function _meta(ta::TimeArray)
    return getfield(ta, :meta)
end
