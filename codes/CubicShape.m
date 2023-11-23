function var = CubicShape(a, xi, der)

if a == 1
    if der == 0
        var = -9/16*xi^3+9/16*xi^2+1/16*xi-1/16;
    elseif der == 1
        var = -27/16*xi^2+9/8*xi+1/16;
    else
    end
elseif a == 2
    if der == 0
        var = 27/16*xi^3-9/16*xi^2-27/16*xi+9/16;
    elseif der == 1
        var = 81/16*xi^2-9/8*xi-27/16;
    else
    end
elseif a == 3
    if der == 0
        var = -27/16*xi^3-9/16*xi^2+27/16*xi+9/16;
    elseif der == 1
        var = -81/16*xi^2-9/8*xi+27/16;
    else
    end
elseif a == 4
    if der == 0
        var = 9/16*xi^3+9/16*xi^2-1/16*xi-1/16;
    elseif der == 1
        var =  27/16*xi^2+9/8*xi-1/16;
    else
    end
end