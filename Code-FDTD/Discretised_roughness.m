%Switching to a function with only 3 states [0 1 2] characterized by a
%certain directional periodicity
%Should be generalized for n states in argument if necessary

function [ New, x, y ] = Discretised_roughness(N,rL,h,lx,ly) 
 x = linspace(-rL/2,rL/2,N); 
 y = linspace(-rL/2,rL/2,N);
[X,Y] = meshgrid(x,y);

RD2D = Roughness_Defined_2D(N,rL,h,lx,ly);
    for i = 1:length(RD2D)
        for j = 1:length(RD2D)
            if (RD2D(i,j) > 0.6e-8)
                New(i,j) = 2;
            elseif (RD2D(i,j) > -0.6e-8)
                New(i,j) = 1;
            else
                New(i,j) = 0;
            end
        end
    end      
end