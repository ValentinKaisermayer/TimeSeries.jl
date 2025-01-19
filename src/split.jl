# when ############################

"""
    when(ta::TimeArray, period::Function, t::Integer)

Filter the `TimeArray` `ta` to include only the entries where the `period` function applied to the timestamps equals `t`.
"""
function when(ta::TimeArray, period::Function, t::Integer)
    ta[findall(period.(timestamp(ta)) .== t)]
end

"""
    when(ta::TimeArray, period::Function, t::String)

Filter the `TimeArray` `ta` to include only the entries where the `period` function applied to the timestamps equals `t`.
"""
function when(ta::TimeArray, period::Function, t::String)
    ta[findall(period.(timestamp(ta)) .== t)]
end

# from, to ######################

"""
    from(ta::TimeArray{T,N,D}, d::D) where {T,N,D}

Return a subset of the `TimeArray` `ta` starting from the timestamp `d`.
"""
function from(ta::TimeArray{T,N,D}, d::D) where {T,N,D}
    if length(ta) == 0
        return ta
    elseif d < timestamp(ta)[1]
        return ta
    elseif d > timestamp(ta)[end]
        return ta[1:0]
    else
        return ta[searchsortedfirst(timestamp(ta), d):end]
    end
end

"""
    to(ta::TimeArray{T,N,D}, d::D) where {T,N,D}

Return a subset of the `TimeArray` `ta` up to the timestamp `d`.
"""
function to(ta::TimeArray{T,N,D}, d::D) where {T,N,D}
    if length(ta) == 0
        return ta
    elseif d < timestamp(ta)[1]
        return ta[1:0]
    elseif d > timestamp(ta)[end]
        return ta
    else
        return ta[1:searchsortedlast(timestamp(ta), d)]
    end
end

###### findall ##################

"""
    findall(ta::TimeArray{Bool,1})

Return the indices of `true` values in the `TimeArray` `ta`.
"""
function Base.findall(ta::TimeArray{Bool,1})
    findall(values(ta))
end

"""
    findall(f::Function, ta::TimeArray{T,1}) where {T}

Return the indices of values in the `TimeArray` `ta` for which the function `f` returns `true`.
"""
function Base.findall(f::Function, ta::TimeArray{T,1}) where {T}
    findall(f, values(ta))
end

"""
    findall(f::Function, ta::TimeArray{T,2}) where {T}

Return the indices of rows in the `TimeArray` `ta` for which the function `f` returns `true`.
"""
function Base.findall(f::Function, ta::TimeArray{T,2}) where {T}
    A = values(ta)
    collect(i for i in axes(A, 1) if f(view(A, i, :)))
end

###### findwhen #################

"""
    findwhen(ta::TimeArray{Bool,1})

Return the timestamps of `true` values in the `TimeArray` `ta`.
"""
function findwhen(ta::TimeArray{Bool,1})
    timestamp(ta)[findall(values(ta))]
end

###### head, tail ###########

"""
    head(ta::TimeArray{T,N}, n::Int=6) where {T,N}

Return the first `n` entries of the `TimeArray` `ta`.
"""
@generated function head(ta::TimeArray{T,N}, n::Int=6) where {T,N}
    new_values = (N == 1) ? :(values(ta)[1:n]) : :(values(ta)[1:n, :])

    quote
        new_timestamp = timestamp(ta)[1:n]
        TimeArray(new_timestamp, $new_values, colnames(ta), meta(ta))
    end
end

"""
    tail(ta::TimeArray{T,N}, n::Int=6) where {T,N}

Return the last `n` entries of the `TimeArray` `ta`.
"""
@generated function tail(ta::TimeArray{T,N}, n::Int=6) where {T,N}
    new_values = (N == 1) ? :(values(ta)[start:end]) : :(values(ta)[start:end, :])

    quote
        start = length(ta) - n + 1
        new_timestamp = timestamp(ta)[start:end]
        TimeArray(new_timestamp, $new_values, colnames(ta), meta(ta))
    end
end

###### first, last ###########

"""
    first(ta::TimeArray)

Return the first entry of the `TimeArray` `ta`.
"""
function Base.first(ta::TimeArray)
    head(ta, 1)
end

"""
    last(ta::TimeArray)

Return the last entry of the `TimeArray` `ta`.
"""
function Base.last(ta::TimeArray)
    tail(ta, 1)
end
