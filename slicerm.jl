function slicerm(g, dim, x0, xargs; N=1000, w=0.5, m=10)
    #= multi-D slice sampler, g: log of distribution function, 
    dim number of dimensions
    x0 start point, 
    xargs extra arguments for distribution, N*dim samples to return,
        w "out step" distance, m max number of outward steps.
    
        Note that you only have to worry about log scaling in 
	'vertical' comparisons
        in these Monte Carlo chain calculations.
    =#
    xs = zeros(N*dim,dim)  # array that will be returned
    xs[1,:] = x0  #x0 should be vector of length dim
    x1 = zeros(1,dim)
    L = zeros(1,dim)
    R = zeros(1,dim)
    way = zeros(1,dim)  # which axis to go along in space
    for i in 2:N*dim
        for d in 1:dim
            way = 0*way #clear it
            way[d] = 1.0 #set nonzero in direction we go for slicing
            y0 = g(x0,xargs)  #height of distribution at x0
            y = y0+log(rand()) # height for slice
            #start stepping out
            U = rand()   # between 0 and 1
            L = x0 - w*way*U
            R = L + w*way
            V = rand()
            J = floor(m*V)
            K = (m - 1) - J
            while J > 0 && y < g(L,xargs)
                L = L - w*way
                J = J - 1
            end
            while K > 0 && y < g(R,xargs)
                R = R + w*way
                K = K - 1
            end
            #now should be stepped out beyond distribution at slice level
            # work back in if no value found:
            Lbar, Rbar = L, R
            while true 
                U = rand()
                x1 = Lbar + U * (Rbar - Lbar)  # vector subtraction should be 
						# correct dir
                if y < g(x1,xargs)
                    break # exit while loop
                end
                if x1[d] < x0[d]
                    Lbar = x1
                    else 
                    Rbar = x1
                end
            end
            xs[i,:] = x1 # found an acceptable point, record in chain (a row)
            x0 = x1 # set initial to new point for next round.
        end # for d
    end  #for N
    return xs 
end
    
    
