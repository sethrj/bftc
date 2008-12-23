function oneDTransport()
   % SN transport in one dimension.
   % Author: Seth R. Johnson
   % Licensed under BSD
   
   N = 4;
   I = 8;
   
   numEquations = I*N; %equations from the finite difference scheme
   numUnknowns  = (I+1)*N; %unknown values of the flux

   %%%%% generate problem data
   
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, 1.0);
   sigma_t  = createXsn(I, 1.0);
   sigma_s0 = createXsn(I, 0.1);
   q        = createSource(I, N, 1.0);
   reflect  = 1;
   
   %%%%% create the transport matrices explicitly
   
   transportMatrix = sparse(numUnknowns, numUnknowns);
   scatterMatrix   = sparse(numUnknowns, numUnknowns);
   
   for i = 1:I
      for n = 1:N
         r = (i-1)*N + n;
         transportMatrix(r,r + N) = transportMatrix(r,r + N) ...
            + mu(n);
         transportMatrix(r,r) = transportMatrix(r,r) ...
            - mu(n);
         
         transportMatrix(r,r + N) = transportMatrix(r,r + N) ...
               + sigma_t(i) * delta_x(i) /2;
         transportMatrix(r,r) = transportMatrix(r,r) ...
               + sigma_t(i) * delta_x(i) /2;
         
         scatterTerm = sigma_s0(i) * delta_x(i) / 2;
         
         for m = 1:N
            s = (i-1)*N + m;
            
            scatterMatrix(r,s + N) = scatterMatrix(r,s + N) ...
                  + scatterTerm * w(m) /2;
               
            scatterMatrix(r,s) = scatterMatrix(r,s) ...
                  + scatterTerm * w(m) /2;
         end
      end
   end
   
   boundariesMatrix = createBoundariesMatrix(I, N, reflect);
   
   figure(1)
   spy(transportMatrix)
   
   figure(2)
   spy(scatterMatrix)
   
   figure(3)
   spy(q)
   
   transportMatrix = transportMatrix - scatterMatrix;
   
   transportMatrix = transportMatrix + boundariesMatrix;
   
   figure(4)
   imagesc(transportMatrix)
   axis square
   
   psi = transportMatrix \ q;
   
%   disp(psi)
   disp(fliplr((reshape(psi, N, I+1))'))
   
   %disp(mu)
   %disp(w)
end

function [phin] = discreteToMoment(n, psi, mu, w, I)
   
end

function [bm]   = createBoundariesMatrix(I, N, reflecting)
   numUnknowns = (I+1)*N;
   bm = sparse(numUnknowns, numUnknowns);
   
   if reflecting == 1
      % set reflecting on left side (i = 0) for incidident directions (n > N/2)
      r = I*N + 1;
      for n = (N/2 + 1):N
         bm(r, n)       = 1;
         bm(r, n - N/2) = -1;
         r = r + 1;
      end
      
      % set reflecting on left side (i = 0) for incidident directions (n < N/2)
      for n = 1:N/2
         bm(r, I*N + n)       = 1;
         bm(r, I*N + n + N/2) = -1;
         r = r + 1;
      end
   else
      % set incoming flux at the boundaries to zero
      r = I*N + 1;
      for n = (N/2 + 1):N
         bm(r, n)       = 1;
         r = r + 1;
      end
      for n = 1:N/2
         bm(r, I*N + n)       = 1;
         r = r + 1;
      end
   end
end

function [q] = createSource(I, N, q0)
   %source in the finite difference equations
   q = ones(I*N, 1) * q0;
   
   %zeros that go along with the boundary condition equations
   q = [q; zeros(N, 1)];
end

function [sigma_x] = createXsn(I, xsn)
   sigma_x = ones(1, I) .* xsn;
end

function [deltas] = createGrid(I, width)
   deltas = ones(I) .* width ./ I;
   % deltas = diff(gridPoints)
end

function [mu w] = createQuadrature(N)
   %N is ordinate set: 2, 4, 8, ...
   %Return gauss-legendre quadrature set
   if (N == 2)
      mu = [.577350269189626];
      w  = [1.0];
   elseif (N == 4)
      mu = [.339981043584856 .861136311594053];
      w  = [.652145154862546 .347854845137454];
   elseif (N == 8)
      mu = [.183434642495650 .525532409916329 .796666477413627 .960289856497536];
      w  = [.362683783378363 .313706645877887 .222381034453374 .101228536290376];
   else
      error('Sorry, invalid N for quadrature')
   end
   
   mu = [-fliplr(mu) mu];
   w  = [-fliplr(w) w];
end