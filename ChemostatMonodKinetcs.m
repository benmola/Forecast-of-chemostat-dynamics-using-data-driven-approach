% x1 - S substrate 
% x2 - X Boimass (Monod)
% D  - u Dilution rate (control input) 

function f = ChemostatMonodKinetcs( t,x,u )

f = [ ((u*20- u.*x(1,:))- (0.6*x(1,:)./(x(1,:)+3)).*x(2,:));
    ( (0.6*x(1,:)./(x(1,:)+3)).*x(2,:)-(u.*x(2,:))) ];
end

%Model from 

%   dx1 = D*(Sin-S)-k*�(S)*X = u*(Sin-x(1))- k*�(x(1))*x(2)
%   dx2   = (k*�(S)-D)*X  = (k*�(x(1))-u)*x(2)
%   y = S = x1

%   �(Sin)= �max*(Sin/(Sin + K))

%   Parameters �(S)= m*(S/(K+S))m=0.6;K=3; k=1; Sin=20;

%   Dm = �*(1-sqrt())
% Constraints (scaled to [i,i] in the final model)

% x1min =0 ; x1max = 20;
% x2min =0 ; x2max =5 ;
% umin = 0; umax = 0.49;
