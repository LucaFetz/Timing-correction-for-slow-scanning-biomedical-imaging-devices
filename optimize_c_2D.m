function C = optimize_c_2D(H, lambda, measurements)
%minimize inverse problem ||y - Hc|| + R to find c.
%input: forward model H, measurements (y) as a column vector, regulation tuning lambda
%output: optimized C


H2 = LinOpMatrix(H);  %put H to LinOpMatrix format
LS=CostL2([],measurements);  % Least-Sqaures data term
F=LS*H2; %composition of cost and H
R = CostL2(size(H,1)); %regulation term. size(H,1) = size(C)
F2 = F+lambda*R; %how to integrate the regulation term?
%regarder evolcost pour voir nb iterations

GD=OptiGradDsct(F2); %GradientDescent
GD.OutOp=OutputOpti(1,[],40); %prendre df0 avec un pas de 1 et pas dt
GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
GD.maxiter=500;         % max number of iterations
GD.run(measurements); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

C = GD.xopt;

end