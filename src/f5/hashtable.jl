# TAKEN AND ADJUSTED FROM GROEBNER.JL
mutable struct Hashvalue
    hash::MonHash
    divmask::DivMask
end

function copy_hashvalue(x::Hashvalue)
    Hashvalue(x.hash, x.divmask)
end

# Hashtable designed to store monomials
mutable struct MonomialHashtable{N}
    exponents::Vector{Monomial{N}}

    # maps exponent hash to its position in exponents array
    hashtable::Vector{MonIdx}

    # stores hashes, division masks,
    # and other valuable info
    # for each hashtable enrty
    hashdata::Vector{Hashvalue}

    #= Monom divisibility =#
    # divisor map to check divisibility faster
    divmap::Vector{UInt32}
    # bits per div variable
    ndivbits::Int

    size::Int
    # elements added
    load::Int
    #
    offset::Int
end

#------------------------------------------------------------------------------

# Returns the next look-up position in the table 
function nexthashindex(h::MonHash, j::MonHash, mod::MonHash)
    (h + j) & mod + MonHash(1)
end

#------------------------------------------------------------------------------

# initialize and set fields for basis hashtable
function initialize_basis_hash_table(::Val{N}, initial_size::Int=2^16) where N

    # for now at most 32 variables
    if N > 32
        error("At most 32 variables currently supported.")
    end
    
    # not necessary to create `initial_size` exponents
    exponents = Vector{Monomial{N}}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(MonIdx, initial_size)

    # exponents[1:load] cover all stored exponents
    # , also exponents[1] is zeroed by default
    load = 1
    size = initial_size

    # exponents array starts from index offset,
    # We store buffer array at index 1
    offset = 2

    # initialize fast divisibility params
    int32bits = 32
    ndivbits = div(int32bits, N)
    # division mask stores at least 1 bit
    # per each of first ndivvars variables
    ndivbits == 0 && (ndivbits += 1)
    divmap = Vector{DivMask}(undef, N * ndivbits)

    # first stored exponent used as buffer lately
    exponents[1] = monomial(SVector{N}(zeros(Exp, N)))

    MonomialHashtable{N}(exponents, hashtable, hashdata,
                         divmap, ndivbits,
                         size, load, offset)
end

# initialize hashtable for `symbolic_preprocessing` 
# These are of the same purpose and structure as basis hashtable,
# but are more local oriented
function initialize_secondary_hash_table(basis_ht::MonomialHashtable{N}) where {N}
    # 2^6 seems to be the best out of 2^5, 2^6, 2^7
    initial_size = 2^6

    exponents = Vector{Monomial{N}}(undef, initial_size)
    hashdata = Vector{Hashvalue}(undef, initial_size)
    hashtable = zeros(MonIdx, initial_size)

    # preserve division info
    divmap = basis_ht.divmap
    ndivbits = basis_ht.ndivbits

    load = 1
    size = initial_size
    offset = 2

    # first stored exponent used as buffer lately
    exponents[1] = monomial(SVector{N}(zeros(Exp, N)))

    MonomialHashtable(exponents, hashtable, hashdata,
                      divmap, ndivbits,
                      size, load, offset)
end

function select_tablesize(nvars, syssize)
    nvars = ring.nvars
    tablesize = 2^10
    if nvars > 4
        tablesize = 2^14
    end
    if nvars > 7
        tablesize = 2^16
    end

    if syssize < 3
        tablesize = divch(tablesize, 2)
    end
    if syssize < 2
        tablesize = divch(tablesize, 2)
    end

    return tablesize
end

#------------------------------------------------------------------------------

ht_resize_threshold() = 0.4
ht_needs_resize(size, load, added) = (load + added)/size > ht_resize_threshold()

function check_enlarge_hashtable!(ht::MonomialHashtable, added::Integer)
    newsize = ht.size
    while ht_needs_resize(newsize, ht.load, added)
        newsize *= 2
    end
    if newsize != ht.size
        ht.size = newsize
        @assert ispow2(ht.size)

        resize!(ht.hashdata, ht.size)
        resize!(ht.exponents, ht.size)
        ht.hashtable = zeros(Int, ht.size)
        
        mod = MonHash(ht.size - 1)

        for i in ht.offset:ht.load
            # hash for this elem is already computed
            he = ht.hashdata[i].hash
            hidx = he
            @inbounds for j in MonHash(1):MonHash(ht.size)
                hidx = nexthashindex(he, j, mod)
                !iszero(ht.hashtable[hidx]) && continue
                ht.hashtable[hidx] = i
                break
            end
        end
    end
    nothing
end

#------------------------------------------------------------------------------

# if hash collision happened
function ishashcollision(ht::MonomialHashtable, vidx, e, he)
    # if not free and not same hash
    @inbounds if ht.hashdata[vidx].hash != he
        return true
    end
    # if not free and not same monomial
    @inbounds if ht.exponents[vidx].exps != e
        return true
    end
    false
end

function insert_in_hash_table!(ht::MonomialHashtable{N}, e::Monomial{N}) where {N}
    # generate hash
    he = Base.hash(e)

    # find new elem position in the table
    hidx = he
    # power of twoooo
    @assert ispow2(ht.size)
    mod = MonHash(ht.size - 1)
    i = MonHash(1)
    hsize = MonHash(ht.size)

    @inbounds while i < hsize
        hidx = nexthashindex(he, i, mod)

        vidx = ht.hashtable[hidx]

        # if free
        iszero(vidx) && break

        # if not free and not same hash
        if ishashcollision(ht, vidx, e, he)
            i += MonHash(1)
            continue
        end

        # already present in hashtable
        return vidx
    end

    # add its position to hashtable, and insert exponent to that position
    vidx = MonIdx(ht.load + 1)
    @inbounds ht.hashtable[hidx] = vidx
    @inbounds ht.exponents[vidx] = copy(e)
    divm = divmask(e, ht.divmap, ht.ndivbits)
    @inbounds ht.hashdata[vidx] = Hashvalue(he, divm)  

    ht.load += 1

    return vidx
end

#------------------------------------------------------------------------------

#=
    Having `ht` filled with monomials from input polys,
    computes ht.divmap and divmask for each of already stored monomials
=#
function fill_divmask!(ht::MonomialHashtable{N}) where N

    min_exp = Vector{UInt64}(undef, N)
    max_exp = Vector{UInt64}(undef, N)

    e = ht.exponents[ht.offset].exps

    @inbounds for i in 1:N
        min_exp[i] = e[i]
        max_exp[i] = e[i]
    end

   @inbounds for i in ht.offset:ht.load # TODO: offset
       make_dense!(e, ht.exponents[i])
       e = ht.exponents[i].exps
       for j in 1:N
           if e[j] > max_exp[j]
               max_exp[j] = e[j]
               continue
           end
           if e[j] < min_exp[j]
               min_exp[j] = e[j]
           end
       end
   end

    ctr = 1
    steps = UInt32(0)
    @inbounds for i in 1:N
        steps = div(max_exp[i] - min_exp[i], UInt32(ht.ndivbits))
        (iszero(steps)) && (steps += UInt32(1))
        for j in 1:ht.ndivbits
            ht.divmap[ctr] = steps
            steps += UInt32(1)
            ctr += 1
        end
    end
    @inbounds for vidx in ht.offset:ht.load
        m = ht.exponents[vidx]
        divm = divmask(m, ht.divmap, ht.ndivbits)
        hsh = Base.hash(m)
        ht.hashdata[vidx] = Hashvalue(hsh, divm)
    end

    nothing
end

#------------------------------------------------------------------------------

# add monomials from `poly` multiplied by exponent vector `etmp`
# with hash `htmp` to hashtable `symbol_ht`,
# and substitute hashes in row
function insert_multiplied_poly_in_hash_table!(row::Vector{MonIdx},
                                               htmp::MonHash,
                                               etmp::Monomial{N},
                                               poly::Vector{MonIdx},
                                               ht::MonomialHashtable{N},
                                               symbol_ht::MonomialHashtable{N}) where N

    # length of poly to add
    len = length(poly)

    mod = MonHash(symbol_ht.size - 1)

    bexps = ht.exponents
    bdata = ht.hashdata

    sexps = symbol_ht.exponents
    sdata = symbol_ht.hashdata

    l = 1 # hardcoding 1 does not seem nice =(
    @label Letsgo
    @inbounds while l <= len
        # we iterate over all monoms of the given poly,
        # multiplying them by htmp/etmp,
        # and inserting into symbolic hashtable

        # hash is linear, so that
        # hash(e1 + e2) = hash(e1) + hash(e2)
        h = htmp + bdata[poly[l]].hash

        e = bexps[poly[l]]

        lastidx = symbol_ht.load + 1
        enew = sexps[1]
        enew = mul!(enew, etmp, e)

        # insert into hashtable
        k = h

        i = MonHash(1)
        ssize = MonHash(symbol_ht.size)
        @inbounds while i <= ssize
            k = nexthashindex(h, i, mod)
            vidx = symbol_ht.hashtable[k]

            # if index is free
            iszero(vidx) && break

            if ishashcollision(symbol_ht, vidx, enew, h)
                i += MonHash(1)
                continue
            end
            
            # hit
            row[l] = vidx
            l += 1

            @goto Letsgo
        end
        # miss

        # add multiplied exponent to hash table        
        sexps[lastidx] = Monomial(enew.deg, copy(enew.exps))
        symbol_ht.hashtable[k] = lastidx

        divm = divmask(enew, symbol_ht.divmap, symbol_ht.ndivbits)
        sdata[lastidx] = Hashvalue(h, divm)

        row[l] = lastidx
        l += 1
        symbol_ht.load += 1
    end

    row
end

function multiplied_poly_to_matrix_row!(
    symbolic_ht::MonomialHashtable{N}, basis_ht::MonomialHashtable{N},
    htmp::MonHash, etmp::Monomial{N}, poly::Vector{MonIdx}) where {N}

    row = similar(poly)
    check_enlarge_hashtable!(symbolic_ht, length(poly))

    insert_multiplied_poly_in_hash_table!(row, htmp, etmp, poly, basis_ht, symbolic_ht)
end

#------------------------------------------------------------------------------

function insert_in_basis_hash_table_pivots(
    row::Vector{ColIdx},
    ht::MonomialHashtable{N},
    symbol_ht::MonomialHashtable{N},
    col2hash::Vector{MonIdx}) where {N}

    check_enlarge_hashtable!(ht, length(row))

    sdata = symbol_ht.hashdata
    sexps = symbol_ht.exponents

    mod = MonHash(ht.size - 1)
    bdata = ht.hashdata
    bexps = ht.exponents
    bhash = ht.hashtable

    l = 1
    @label Letsgo
    @inbounds while l <= length(row)
        hidx = col2hash[row[l]]

        # symbolic hash
        h = sdata[hidx].hash

        lastidx = ht.load + 1
        bexps[lastidx] = sexps[hidx]
        e = bexps[lastidx]

        k = h
        i = MonHash(1)
        @inbounds while i <= ht.size
            k = nexthashindex(h, i, mod)
            hm = bhash[k]

            iszero(hm) && break

            if ishashcollision(ht, hm, e, h)
                i += MonHash(1)
                continue
            end

            row[l] = hm
            l += 1
            @goto Letsgo
        end

        bhash[k] = pos = lastidx
        row[l] = pos
        l += 1

        bdata[pos] = Hashvalue(h, sdata[hidx].divmask)

        ht.load += 1
    end

    nothing
end
