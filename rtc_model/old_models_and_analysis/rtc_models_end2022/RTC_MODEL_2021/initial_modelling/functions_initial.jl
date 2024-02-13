function separate_species!(species, sol)
    solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
    b = solDF[!, :b];
    c = solDF[!, :c];
    d = solDF[!, :d];
    p = solDF[!, :p];
    return b, c, d, p
end

function separate_species2!(species, sol)
    solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
    a = solDF[!, :a];
    b = solDF[!, :b];
    c = solDF[!, :c];
    d = solDF[!, :d];
    e = solDF[!, :e];
    p = solDF[!, :p];
    return a, b, c, d, e, p
end

function separate_species3!(species, sol)
    solDF = DataFrame([[j[i] for j in sol.u] for i=1:length(sol.u[1])], species);
    b = solDF[!, :b];
    c = solDF[!, :c];
    d = solDF[!, :d];
    e = solDF[!, :e];
    p = solDF[!, :p];
    return b, c, d, e, p
end