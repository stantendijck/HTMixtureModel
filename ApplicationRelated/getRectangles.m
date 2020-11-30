function A = getRectangles(xgrid, ygrid)

A = zeros(1,4); cntr = 1;
for ix = 1:length(xgrid)-1
    for iy = 1:length(ygrid)-1
        A(cntr,:) = [xgrid(ix),xgrid(ix+1),ygrid(iy),ygrid(iy+1)];
        cntr = cntr + 1;
    end
end


end