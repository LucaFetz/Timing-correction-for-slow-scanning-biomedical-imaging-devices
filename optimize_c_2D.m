function C = optimize_c_2D(H, lambda, measurements, opti_type, regul_type)
%minimize inverse problem ||y - Hc|| + R to find c.
%input: forward model H, measurements (y) as a column vector, regulation tuning lambda
%output: optimized C


H = LinOpMatrix(H);  %put H to LinOpMatrix format
% If true, applyHtH method will save its result. Hence for two consecutive HtH*x with the
% same x, no computation is done for the second call
H.memoizeOpts.applyHtH=1;  

LS=CostL2([],measurements);  % Least-Sqaures data term
F=LS*H; %composition of cost and H
switch regul_type
    case 'L2'
        R = CostL2(size(H,1)); %regulation term. size(H,1) = size(C)
    case 'L1'
        R = CostL1(size(H,1));
    case 'L21'
        R = CostTV(size(H,1));
end

F2 = F+lambda*R; 

% For the CostL2, the precomputation save Hty=H'*y and then the gradient is
% computed using H.HtH(x) - Hty and the apply is computed using the HtH
% method
F2.doPrecomputation=1;
%regarder evolcost pour voir nb iterations

switch opti_type
    case 'GradDsct'
        GD=OptiGradDsct(F2); %GradientDescent
        GD.OutOp=OutputOpti(1,[],40); %prendre df0 avec un pas de 1 et pas dt
        GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
        GD.maxiter=500;         % max number of iterations
        GD.run(measurements); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

        C = GD.xopt;
    case 'ADMM'
        
end

end