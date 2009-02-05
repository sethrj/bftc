function oneDTransportLD()
   % SN transport in one dimension, linear discontinuous.
   % Author: Seth R. Johnson
   % Licensed under BSD
   
   global I N sigma_t  sigma_s0  q reflect  width;
   
   I = 6;  % number of spatial cells
   N = 4;  % number of discrete ordinates
   
   sigma_t  = createXsn(I, 2.0);
   sigma_s0 = createXsn(I, 1.0);
   reflect  = [0 0]; %reflecting bounds on left and right
   width    = 3.0;
   
   q        = createSource(I, 1.0);
   
   % change a piece of the cross sections
   %sigma_t(3) = 10;
   %q(3) = 0;
   %   sigma_t(I/2 - I/32 + (1:I/16))   = 15.0;
   %   q(I/2 - I/32 + (1:I/16))         = 0.0;
   
   
   %%% uncomment any of these three lines to run different codes
   
   %  runTransportMatrixTest()
   %  runIterateMfpEigenvalues()
     runSourceDriven();
   %  runEigenvalue();
end

function runIterateMfpEigenvalues()
   % look at how changing the width of a problem with two reflecting
   % boundaries changes the eigenvalue clustering
   
   global transportMatrix;
   global reflect  width;
   
   %change problem parameters
   reflect  = [1 1];
   
   
   iterate = [1 2 5 10 20]; %width of the problem
   
   shapes = {'bx' 'rx' 'gx' 'kx' 'bo' 'ro' 'go' 'ko'};
   legtext = cell(1, length(iterate));
   
   figure(4)
   clf
   for i = 1:length(iterate)
      width = iterate(i);
      legtext{i} = sprintf('%.2g', width);
      
      runSourceDriven();
      
      [lambda] = eig(full(transportMatrix));
      figure(4)
      plot(real(lambda), imag(lambda), shapes{i})
      hold on
   end
   hold off
   legend(legtext)
end

function runTransportMatrixTest()
   global transportMatrix;
   global I N reflect sigma_t  sigma_s0 width;
      
   %%%%% generate problem data
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, width);
   
   %%%%% create the transport matrices explicitly
   
   transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w)...
   + createBoundariesMatrix(I, N, reflect);
   figure(1)
   advanceSpy(transportMatrix, 4, 0);
   drawmatlabels();
end

function runEigenvalue()
   global I N nusigma_f sigma_t  sigma_s0  reflect width;
   
   numEigenvals = 6;
   
   %%%%% generate problem data
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, width);
   nusigma_f= createXsn(I, 1.5);
   
   %   sigma_t(240:250) = 10.0;
   nusigma_f(I/2 + (1:I/32)) = 0.0;
   sigma_t(I/2 + (1:I/32))   = 10.0;
   
   %%%%% create the transport matrices explicitly
   transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w);
   bm              = createBoundariesMatrix(I, N, reflect);
   transportMatrix = transportMatrix + bm;
   
   fissionMatrix   = createFissionMatrix(delta_x, nusigma_f, mu, w);
   
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
   xlabel('Re \lambda')
   ylabel('Im \lambda')
   
   figure(2)
   clf
   grid = [0 cumsum(delta_x)];
   legendText = cell(1, numEigenvals);
   for e = 1:numEigenvals
      plot(grid, phi(:,e)')
      hold all
      legendText{e} = sprintf('%d', e);
   end
   xlabel('x')
   ylabel('\phi')
   legend(legendText)
end

function runSourceDriven()
   global transportMatrix;
   
   global I N q sigma_t  sigma_s0  reflect width;
   
   
   %%%%% generate problem data
   [mu w]   = createQuadrature(N);
   delta_x  = createGrid(I, width);
   
   %%%%% create the transport matrices explicitly
   
   transportMatrix = createTransportMatrix(delta_x, sigma_t,sigma_s0, mu, w)...
      + createBoundariesMatrix(I, N, reflect);
   
   q = isotropicSourceVector(N, q);
   
   figure(3)
   if (I*N < 150)
      %   imagesc(log10(abs(transportMatrix)))
      %      imagesc(transportMatrix)
      %      axis square
      advanceSpy([transportMatrix q], 3, 0);
      numUnknowns = I*N*2 + N;
      line([1 1]*numUnknowns + 0.5, [0 1]*numUnknowns + 0.5,'Color','k')
   end
   
   %%%%% annnnd.... OMG TRANSPORT:
   %tic
   psi = transportMatrix \ q;
   %toc
   
   %%%%% display the results
   grid = [0 cumsum(delta_x)];
%    figure(1)
%    plot(grid, discreteToMoment(0,psi, mu, w)')
%    xlabel('x')
%    ylabel('\phi')
   
   figure(2)
   plotAngularFlux(grid, psi, mu)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotAngularFlux(x, psi, mu)
   global I N thelabels

   % two x's for psi on the left and right
   x2 = reshape(repmat(x, 2, 1), 1, (I+1)*2);
   x2 = x2(1:end);
   
   psi2 = nan*ones((I+1)*2, N);
   psil = cell((I+1)*2, N);
   a = 1;
   for n = N/2+1:N
      psi2(1,n) = psi(a);
      psil(1,n) = thelabels(a);
      a = a+1;
   end
   
   for i = 1:I %spatial cell
      for n = 1:N %ordinate         
         for g = 1:2 %basis function, 1 = left, 2 = right
            psi2(2*i - 1 + g,n) = psi(a);
            psil(2*i - 1 + g,n) = thelabels(a);
            a = a+1;
         end
      end
   end
   
   for n = 1:N/2
      psi2(end,n) = psi(a);
      psil(end,n) = thelabels(a);
      a = a+1;
   end
   
   clf   
   %    [X Y] = meshgrid(x, mu);
   %    surf(X, Y, psi2');
   %    xlabel('x')
   %    ylabel('\mu')
   %    zlabel('\psi_n')
   %    return
   
   disp([num2cell(x2') psil])
   disp([x2' psi2])
   %if (N < 10)
      legendText = cell(1, N);
      %2-D plot with different colors
      for n = 1:N
         plot(x2, psi2(:,n))
         hold all
         legendText{n} = sprintf('\\mu=%.4f', mu(n));
      end
      xlabel('x')
      ylabel('\psi_n')
      legend(legendText)
%    else
%       %3-D plot
%       for n = 1:N
%          plot3(x, mu(n) * ones(1, I), psi2(:,n))
%          hold on
%       end
%       xlabel('x')
%       ylabel('\mu')
%       zlabel('\psi_n')
%       grid on
%    end
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


function drawmatlabels()
   global thelabels;
   N = length(thelabels);
   
   text(0,0,'n,D,i','HorizontalAlignment','right','VerticalAlignment','middle');
   for r = 1:length(thelabels)
      text(0,r,thelabels{r},'HorizontalAlignment','right','VerticalAlignment','middle');
      text(r,0,thelabels{r},'HorizontalAlignment','left','VerticalAlignment','middle',...
         'Rotation',90);
   end
   grid on
   set(gca, 'YAxisLocation','right')
   set(gca,'XTick',1:N);
   set(gca,'YTick',1:N);
   set(gca,'GridLineStyle',':')
end


function [bm]   = createBoundariesMatrix(I, N, reflecting)
   numUnknowns  = I*N*2 + N; %I cells, N directions, 2 basis funcs
   bm = sparse(numUnknowns, numUnknowns);
   assert(length(reflecting) == 2, 'reflecting should have left, right')
   
   global thelabels
   
   if reflecting(1) == 1
      % set reflecting on left side (i = 1) for incidident directions (n > N/2)
      % sets \Psi_{n,R,0} to \Psi_{n',L,1}
      assert(false)
   else
      % set \Psi_{n,R,0} for \mu_n > 0
      for n = N/2+1:N
         r = n - N/2;
         bm(r, r) = 1;
         thelabels{r} = sprintf('%d,%s,%d',n,'R',0);
      end
   end
   
   if reflecting(2) == 1
      % set reflecting on right side (i = I) for incidident directions (n < N/2)
      % sets \Psi_{n,L,I+1} to \Psi_{n',R,1}
      assert(false)
   else
      % sets \Psi_{n,L,I+1} for \mu_n > 0
      for n = 1:N/2
         r = I*N*2 + n + N/2;
         bm(r, r) = 1; 
         thelabels{r} = sprintf('%d,%s,%d',n,'L',I+1);
      end
   end
end

function [transportMatrix] = createTransportMatrix(delta_x, sigma_t, ...
      sigma_s0, mu, w)
   % generate the SN linear discontinous equations
   N = length(mu);
   I = length(sigma_t);
   numUnknowns  = I*N*2 + N; %I cells, N directions, 2 basis funcs
   
   massmatrix = [1/3 1/6; 1/6 1/3]; %regular 
   %massmatrix = [1/2 0; 0 1/2]; %lumped
   
   transportMatrix = sparse(numUnknowns, numUnknowns);
   
   global thelabels;
   
   thelabels = cell(1, numUnknowns);
   
   lr = 'LR';
   
   for i = 1:I %spatial cell
      for n = 1:N %ordinate
         
         r = ((i-1)*N + (n-1))*2 + N/2; 
         
         % set the indices for \psi_{n,i-1/2} and \psi_{n,i+1/2}
         if (n > N/2) %going right
            right = r + 2;       %our own right
            if (i == 1)
               left  = n - N/2; %boundary on the left
            else
               left  = r - N*2 + 2; %(i-1)'s right
            end
         else
            if (i == I)
               right = numUnknowns; %boundary on the right
            else
               right = r + N*2 + 1; %(i+1)'s left, same angle
            end
            left  = r + 1;       %our own left
         end

         
         if (left < 1)
            disp(sprintf('Left is less than 1: i=%d,n=%d, left=%d',...
               i,n,left));
            left = 1;
         end
         
         if (right > numUnknowns)
            disp(sprintf('Right is greater than %d: i=%d,n=%d, left=%d',...
               numUnknowns,i,n,right));
            left = 1;
         end
         
         for g = 1:2 %basis function, 1 = left, 2 = right
            thelabels{r+g} = sprintf('%d,%s,%d',n,lr(g),i);
            %streaming term
            if (g == 1)
                % \psi_{n,i}
               transportMatrix(r+g,r+1) = transportMatrix(r+g,r+1) ...
                  + mu(n) / delta_x(i) * 0.5;
               transportMatrix(r+g,r+2) = transportMatrix(r+g,r+2) ...
                  + mu(n) / delta_x(i) * 0.5;
               
                % - \psi_{n,i-1/2}
               transportMatrix(r+g,left) = transportMatrix(r+g,left) ...
                 - mu(n) / delta_x(i) ;
            else
                % \psi_{n,i+1/2}
               transportMatrix(r+g,right) = transportMatrix(r+g,right) ...
                  + mu(n) / delta_x(i);
               
                % - \psi_{n,i}
               transportMatrix(r+g,r+1) = transportMatrix(r+g,r+1) ...
                  - mu(n) / delta_x(i) * 0.5;
               transportMatrix(r+g,r+2) = transportMatrix(r+g,r+2) ...
                  - mu(n) / delta_x(i) * 0.5;
            end

            %absorption term
            transportMatrix(r+g,r+1) = transportMatrix(r+g,r+1) ...
              + 1;% + massmatrix(g,1)*sigma_t(i);
            transportMatrix(r+g,r+2) = transportMatrix(r+g,r+2) ...
              + 1;% + massmatrix(g,2)*sigma_t(i);

            %scattering term
            scatterTerm = sigma_s0(i) / 2;
            for m = 1:N
               s = ((i-1)*N + (m-1))*2 + N/2;
               % \Phi_{L,i}
               transportMatrix(r+g,s+1)     = transportMatrix(r+g,s+1) ...
                  - scatterTerm * w(m);
               
               % \Phi_{R,i}
               transportMatrix(r+g,s+2)     = transportMatrix(r+g,s+2) ...
                  - scatterTerm * w(m);
            end
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

function [q] = createSource(I, q0)
   q = ones(1, I) .* q0;
end

function [newq] = isotropicSourceVector(N, q)
   %uniform, isotropic source
   I = length(q);
   q = q ./ 4 ; % 2 for isotropic, 2 for splitting it up into sum of left and right
   newq = repmat(q, N*2, 1);
   newq = reshape(newq, I * N * 2, 1);
   
   %zeros that go along with the boundary condition equations
   newq = [zeros(N/2, 1); newq; zeros(N/2, 1)];
end

function [sigma_x] = createXsn(I, xsn)
   sigma_x = ones(1, I) .* xsn;
end

function [deltas] = createGrid(I, width)
   deltas = ones(1, I) .* width ./ I;
   % deltas = diff(gridPoints)
end

function [mu w] = createQuadrature(N)
   % N is ordinate set: 2, 4, 8, ...
   % Return gauss-legendre quadrature set
   % 
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