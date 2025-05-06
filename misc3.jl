N = 9

function ST_indexing(pos::Int64, N::Int64)

    if mod(N,2) == 1
        ST_index = mod1(2*pos - 1, N)

    # If pos is in an even site, considering ST starts in an odd site:
    elseif pos > Int64(N/2)
        ST_index = mod1(2*pos, N)
    else
        ST_index = mod1(1 + 2*(pos - 1), N)
    end

    return ST_index
end

for i in 1:N
    print(ST_indexing(i, N))
end