function var = PolyShape(a, xi, der)

if a == 1
    if der == 0
        var = 0.5 * (1-xi);
    elseif der == 1
        var = -0.5;
    else
    end
elseif a == 2
    if der == 0
        var = 0.5 * (1+xi);
    elseif der == 1
        var = 0.5;
    else
    end
end