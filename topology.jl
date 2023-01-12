connectivitydict = Dict(
    2 => 2,
    3 => 2,
    10 => 2,
    11 => 4,
    12 => 2,    #
    13 => 2,    
    14 => 4,    #
    15 => 2,    #
    27 => 2,
    37 => 2,
    48 => 2,
    49 => 2,
    50 => 2,
    52 => 2,    #
    53 => 4,
    54 => 4,    #
    56 => 4,
    58 => 4,    #
    60 => 2,
    66 => 2,
    68 => 2,
    70 => 2,    #
    75 => 2,
    77 => 2,
    81 => 2,    #
    83 => 2,
    84 => 2,
    85 => 4,
    86 => 2,
    87 => 3,
    88 => 2,
    103 => 4,
    124 => 4,
    128 => 6,
    130 => 8,
    147 => 2,   #
    148 => 2,   #
    162 => 2,
    163 => 4,
    164 => 2,   #
    165 => 4,
    166 => 2,   #
    167 => 4,   #
    168 => 2,
    171 => 3,
    172 => 3, 
    175 => 2,
    176 => 4,   #
    184 => 4,
    192 => 4,
    201 => 4,
    203 => 2)   #

# Space groups that have filling-enforced topology
fillingenforced = [13, 48, 49, 50, 68, 86]

degendimdict = Dict(
    79 => 2,
    87 => 2,
    97 => 2,
    98 => 2,
    107 => 2,
    139 => 2,
    209 => 2,
    218 => 6,
    220 => 6,
    222 => 6,
    225 => 2,
    223 => 6,
    230 => 6
)

# Space groups with 6-band degeneracies
degenklabdict = Dict(
    79 => "P",
    87 => "P",
    97 => "P",
    98 => "P",
    218 => "R",
    220 => "H",
    222 => "R",
    223 => "R",
    230 => "H"
)

# Space groups with 6-band degeneracies
degenksymdict = Dict(
    79 => :P,
    87 => :P,
    97 => :P,
    98 => :P,
    218 => :R,
    220 => :H,
    222 => :R,
    223 => :R,
    230 => :H
)

degenchardict = Dict(
    # Not sure how to interpret character tables for sg=218 and 220
    79 => "P₂P₂",
    87 => "P₃P₄",
    97 => "P₃P₄",
    98 => "P₁",
    222 => "R₄",
    223 => "R₄",
    230 => "H₄"
)


function isfilling(sg)
    global fillingenforced
    return sg in fillingenforced
end

function minconnectivity(sg; order=1)
    global connectivitydict
    if haskey(connectivitydict, sg)
        return connectivitydict[sg]
    else
        @warn """This space group is not in the connectivitydict. 
            This may be because either this space group is not topological, 
            or I haven't finished implementing this."""
        return 0
    end
end

function getdegendim(sg)
    if haskey(degendimdict, sg)
        return degendimdict[sg]
    else
        @warn """This space group is not in the degeneracies dict. 
            This may be because either this space group is not topological, 
            or I haven't finished implementing this."""
        return 1
    end
end

function getdegenklab(sg)
    if haskey(degenklabdict, sg)
        return degenklabdict[sg]
    else
        @warn """This space group is not in the degeneracies dict. 
            This may be because either this space group is not topological, 
            or I haven't finished implementing this."""
    end
end

function getdegenksym(sg)
    if haskey(degenksymdict, sg)
        return degenksymdict[sg]
    else
        @warn """This space group is not in the degeneracies dict. 
            This may be because either this space group is not topological, 
            or I haven't finished implementing this."""
    end
end

function getdegenchar(sg)
    if haskey(degenchardict, sg)
        return degenchardict[sg]
    else
        @warn """This space group is not in the degeneracies dict. 
            This may be because either this space group is not topological, 
            or I haven't finished implementing this."""
    end
end


# CharacterTable(get_lgirreps(222, 3)["R"])
# println(bandirsd["R"])    # Print degeneracies
# println(lgirsd["R"])
# println(CharacterTable(lgirsd["R"]))
# println(bands)
# println(nds)