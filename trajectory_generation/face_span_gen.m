
function  V = face_span_gen(X, Y, mu)
    addpath('polytopes_2017_10_04_v1.9\')
    U = [-1,  0, -mu,  0,  0,  0;
         +1,  0, -mu,  0,  0,  0;
          0, -1, -mu,  0,  0,  0;
          0,  1, -mu,  0,  0,  0;
          0,  0, - Y, -1,  0,  0;
          0,  0, - Y,  1,  0,  0;
          0,  0, - X,  0, -1,  0;
          0,  0, - X,  0,  1,  0;
         -Y, -X, -(X+Y)*mu,  mu,  mu, -1;
         -Y,  X, -(X+Y)*mu,  mu, -mu, -1;
          Y, -X, -(X+Y)*mu, -mu,  mu, -1;
          Y,  X, -(X+Y)*mu, -mu, -mu, -1;
          Y,  X, -(X+Y)*mu,  mu,  mu,  1;
          Y, -X, -(X+Y)*mu,  mu, -mu,  1;
         -Y,  X, -(X+Y)*mu, -mu,  mu,  1;
         -Y, -X, -(X+Y)*mu, -mu, -mu,  1];

    V = span_of_face(U);
    %V = licols(V1);
end


function [Xsub,idx]=licols(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%
%    [Xsub,idx]=licols(X)
%
%in:
%
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%
%out:
%
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
   if ~nnz(X) %X has no non-zeros and hence no independent columns
       
       Xsub=[]; idx=[];
       return
   end
   if nargin<2, tol=1e-10; end
   
           
     [Q, R, E] = qr(X,0); 
     
     if ~isvector(R)
      diagr = abs(diag(R));
     else
      diagr = abs(R(1));   
     end
     %Rank estimation
     r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
     idx=sort(E(1:r));
     Xsub=X(:,idx);
end
