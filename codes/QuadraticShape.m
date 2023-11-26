function var = QuadraticShape(a, xi, der)

if a == 1
    if der == 0
        var = 0.5*xi*(xi-1);
    elseif der == 1
        var = xi-0.5;
    else
    end
elseif a == 2
    if der == 0
        var = 1-xi^2;
    elseif der == 1
        var = -2*xi;
    else
    end
elseif a == 3
    if der == 0
        var = 0.5*xi*(xi+1);
    elseif der == 1
        var = xi+0.5;
    else
    end
end