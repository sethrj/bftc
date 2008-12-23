function oneDTransport()
   % SN transport in one dimension.
   % Author: Seth R. Johnson
   % Licensed under BSD
   
   I = 128; % number of spatial cells
   N = 8;  % number of discrete ordinates
   
%   runSourceDriven(I, N);
   runEigenvalue(I, N);
end   

function runEigenvalue(I, N)
   
   numEigenvals = 5;
   
   %%%%% generate problem data
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, 50.0);
   sigma_t  = createXsn(I, 1.2);
   sigma_s0 = createXsn(I, 0.4);
   nusigma_f= createXsn(I, 0.8);
   reflect  = [1 1];
   
%   sigma_t(240:250) = 10.0;
   nusigma_f(I/2 + (1:I/32)) = 0.0;
   sigma_t(I/2 + (1:I/16))   = 10.0;
   
   %%%%% create the transport matrices explicitly
    transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w);
    bm = createBoundariesMatrix(I, N, reflect);
    transportMatrix = transportMatrix + bm;
%   transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w);...
%                     + createBoundariesMatrix(I, N, reflect);
   
   fissionMatrix   = createFissionMatrix(delta_x, nusigma_f, mu, w);
                 %    + createBoundariesMatrix(I, N, reflect);

   if (I*N < 150)
      figure(3)
      imagesc(fissionMatrix)
      axis square
   end

   
   %%%%% annnnd.... OMG TRANSPORT:
   operator = @(x) transportMatrix\(fissionMatrix * x);
   options.disp = 0;
   [psi, lambda] = eigs(operator, (I+1)*N, numEigenvals, 'lm', options);
   
   %%%%% fix up the eigenvectors
   
   lambda = diag(lambda);
   
   phi = zeros(I+1, numEigenvals);
   % scale the eigenvectors and put them into the flux
   for i=1:numEigenvals
       phi(:,i) = discreteToMoment(0, psi(:,i), mu, w);
       
       sgn = sign(sum(phi(:,i)));
       if (sgn == 0)
          sgn = sign(phi(1,i));
       end
       phi(:,i) = phi(:,i) * sgn;
       phi(:,i) = phi(:,i)./max(phi(:,i)).*lambda(i) ;
   end
   
   %%%%% display the results
   
   figure(1)
   plot(real(lambda), imag(lambda), 'bx')
   disp(sprintf('dominant eigenvalue: %.4g', lambda(1)))
   
   figure(2)
   clf
   grid = [0 cumsum(delta_x)];
   legendText = {};
   for i = 1:numEigenvals
      plot(grid, phi(:,i)')
      hold all
      legendText{end+1} = sprintf('%d', i); %#ok<AGROW>
   end
   xlabel('x')
   ylabel('\phi')
   legend(legendText)
end

function runSourceDriven(I, N)
   %%%%% generate problem data
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, 6.0);
   
   sigma_t  = createXsn(I, 1.0);
   sigma_s0 = createXsn(I, 0.5);   
   q        = createSource(N, delta_x, 1.0);
   reflect  = [1 1];
   
   sigma_t(12:14) = 5.0;
   
   %%%%% create the transport matrices explicitly
   
   transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w)...
                     + createBoundariesMatrix(I, N, reflect);
   
   figure(4)
   imagesc(log10(abs(transportMatrix)))
   axis square
   
   %%%%% annnnd.... OMG TRANSPORT:
   psi = transportMatrix \ q;
   
   %%%%% display the results
%   disp((reshape(psi, N, I+1))')
   
   figure(5)
   plot([0 cumsum(delta_x)], discreteToMoment(0,psi, mu, w)')
   xlabel('x')
   ylabel('\phi')
end

function [phil] = discreteToMoment(l, psi, mu, w)
   N  = length(mu);
   Ip = length(psi) / N; % Ip = I + 1
   
   phil = zeros(Ip, 1);
   
   if (l == 0)
      r = 1;
      for i = 1:Ip
         for n = 1:N
            phil(i) = phil(i) + w(n) * psi(r);
            r = r + 1;
         end
      end
   else
      error('Higher moments not yet supported.')
   end
end

function [bm]   = createBoundariesMatrix(I, N, reflecting)
   numUnknowns = (I+1)*N;
   bm = sparse(numUnknowns, numUnknowns);
   assert(length(reflecting) == 2, 'reflecting should have left, right')
   
   if reflecting(1) == 1
      % set reflecting on left side (i = 0) for incidident directions (n > N/2)
      r = I*N + 1;
      for n = (N/2 + 1):N
         bm(r, n)       = 1;
         bm(r, N + 1 - n) = -1;
         r = r + 1;
      end
   else
      % set incoming flux at the boundaries to zero
      r = I*N + 1;
      for n = (N/2 + 1):N
         bm(r, n)       = 1;
         r = r + 1;
      end
   end
   
   if reflecting(2) == 1
      % set reflecting on left side (i = 0) for incidident directions (n < N/2)
      for n = 1:N/2
         bm(r, I*N + n)       = 1;
         bm(r, I*N + N + 1 - n) = -1;
         r = r + 1;
      end
   else
      % set incoming flux at the boundaries to zero
      for n = 1:N/2
         bm(r, I*N + n)       = 1;
         r = r + 1;
      end
   end
end

function [transportMatrix] = createTransportMatrix(delta_x, sigma_t, ...
                                                   sigma_s0, mu, w)
   % generate the SN finite difference equations
   N = length(mu);
   I = length(sigma_t);
   numUnknowns  = (I+1)*N;

   transportMatrix = sparse(numUnknowns, numUnknowns);
   
   for i = 1:I
      for n = 1:N
         r = (i-1)*N + n;
         
         %streaming term
         transportMatrix(r,r + N) = transportMatrix(r,r + N) ...
            + mu(n);
         transportMatrix(r,r) = transportMatrix(r,r) ...
            - mu(n);
         
         %absorption term
         transportMatrix(r,r + N) = transportMatrix(r,r + N) ...
               + sigma_t(i) * delta_x(i) /2;
         transportMatrix(r,r) = transportMatrix(r,r) ...
               + sigma_t(i) * delta_x(i) /2;
         
         %scattering term         
         scatterTerm = sigma_s0(i) * delta_x(i) / 2;
         for m = 1:N
            s = (i-1)*N + m;
            
            transportMatrix(r,s + N) = transportMatrix(r,s + N) ...
                  - scatterTerm * w(m) /2;
               
            transportMatrix(r,s)     = transportMatrix(r,s) ...
                  - scatterTerm * w(m) /2;
         end
      end
   end
end

function [fissionMatrix] = createFissionMatrix(delta_x, nusigma_f, mu, w)
   % generate the SN finite difference equations
   N = length(mu);
   I = length(nusigma_f);
   numUnknowns  = (I+1)*N;

   fissionMatrix = sparse(numUnknowns, numUnknowns);
   
   for i = 1:I
      for n = 1:N
         r = (i-1)*N + n;
         
         %fission term         
         scatterTerm = nusigma_f(i) * delta_x(i) / 2;
         for m = 1:N
            s = (i-1)*N + m;
            
            fissionMatrix(r,s + N) = fissionMatrix(r,s + N) ...
                  + scatterTerm * w(m) /2;
               
            fissionMatrix(r,s)     = fissionMatrix(r,s) ...
                  + scatterTerm * w(m) /2;
         end
      end
   end
end

function [q] = createSource(N, delta_x, q0)
   %uniform, isotropic source in the finite difference equations
   I = length(delta_x);
   q = q0 ./ 2 .* delta_x;
   q = repmat(q, N, 1);
   q = reshape(q, I * N, 1);
   
   %zeros that go along with the boundary condition equations
   q = [q; zeros(N, 1)];
end

function [sigma_x] = createXsn(I, xsn)
   sigma_x = ones(1, I) .* xsn;
end

function [deltas] = createGrid(I, width)
   deltas = ones(1, I) .* width ./ I;
   % deltas = diff(gridPoints)
end

function [mu w] = createQuadrature(N)
   %N is ordinate set: 2, 4, 8, ...
   %Return gauss-legendre quadrature set
   if (N == 2)
      mu = [.577350269189626]; %#ok<NBRAK>
      w  = [1.0]; %#ok<NBRAK>
   elseif (N == 4)
      mu = [.339981043584856 .861136311594053];
      w  = [.652145154862546 .347854845137454];
   elseif (N == 8)
      mu = [.183434642495650 .525532409916329 .796666477413627 .960289856497536];
      w  = [.362683783378363 .313706645877887 .222381034453374 .101228536290376];
   elseif (N == 16)
      mu = [0.989400934991650,0.944575023073233,0.865631202387832,...
         0.755404408355003,0.617876244402644,0.458016777657227,...
         0.281603550779259,0.0950125098376370;];
      w = [0.0271524594117540,0.0622535239386480,0.0951585116824930,...
         0.124628971255534,0.149595988816577,0.169156519395003,...
         0.182603415044924,0.189450610455067;];
   elseif (N == 32)
      mu = [0.0483076656877380,0.144471961582796,0.239287362252137,...
         0.331868602282128,0.421351276130635,0.506899908932229,...
         0.587715757240762,0.663044266930215,0.732182118740290,...
         0.794483795967942,0.849367613732570,0.896321155766052,...
         0.934906075937740,0.964762255587506,0.985611511545268,...
         0.997263861849482];
      w = [0.0965400885147260,...
         0.0956387200792750,0.0938443990808050,0.0911738786957640,...
         0.0876520930044040,0.0833119242269470,0.0781938957870700,...
         0.0723457941088490,0.0658222227763620,0.0586840934785360,...
         0.0509980592623760,0.0428358980222270,0.0342738629130210,...
         0.0253920653092620,0.0162743947309060,0.007018610009470];
   else
      error('Sorry, invalid N for quadrature')
   end
   
   assert(length(mu) == length(w));
   assert(abs(sum(w) - 1.0) < 2*eps);
   assert(all(mu > 0.0) && all(mu < 1.0));
   
   mu = [-fliplr(mu) mu];
   w  = [fliplr(w) w];
   
   assert(all(w > 0));
end