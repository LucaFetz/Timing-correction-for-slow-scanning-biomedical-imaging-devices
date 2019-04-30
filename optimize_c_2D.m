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
    case "L2"
        R = CostL2(H.sizein); %regulation term. size(H,1) = size(C)
    case "L1"
        R = CostL1(H.sizein);
    case "TV"
        R = CostTV(H.sizein);
    case "L12"
        R=CostMixNorm21([H.sizein,1],3);
end

F2 = F+lambda*R; 

% For the CostL2, the precomputation save Hty=H'*y and then the gradient is
% computed using H.HtH(x) - Hty and the apply is computed using the HtH
% method
F2.doPrecomputation=1;
%regarder evolcost pour voir nb iterations

switch opti_type
    case "GradDsct"
        GD=OptiGradDsct(F2); %GradientDescent
        GD.OutOp=OutputOpti(1,[],40); %prendre df0 avec un pas de 1 et pas dt
        GD.ItUpOut=2;           % call OutputOpti update every ItUpOut iterations
        GD.maxiter=500;         % max number of iterations
        GD.run(measurements); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

        C = GD.xopt;
        
    case "ConjGrad"
        %Note: ConjGrad is used only with L2 cost
        if(lambda ~= 0)
            A = H.makeHtH() + lambda * LinOPIdentity(H.sizein);
        else
            A = H.makeHtH();
        end
        b = H'*measurements;
        CG=OptiConjGrad(A,b);  
        CG.OutOp=OutputOptiConjGrad(1,dot(measurements(:),measurements(:)),[],40);  
        CG.ItUpOut=2; 
        CG.maxiter=500;                             % max number of iterations
        CG.run(measurements);                                  % run the algorithm 

        C = CG.xopt;
        
    case "ADMM"
        Fn={lambda*R};
        Hn = {H};
        rho_n=[1e-1];
        ADMM=OptiADMM(F,Fn,Hn,rho_n);
        %ADMM.OutOp=OutputOptiSNR(1,[],10);
        % STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower
        % than 1e-4 or when the distance between two successive step is lower than 1e-4
        ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);   
        ADMM.ItUpOut=1;             % call OutputOpti update every ItUpOut iterations
        ADMM.maxiter=200;           % max number of iterations
        ADMM.run(zeros(size(measurements)));   % run the algorithm 
        
        C = ADMM.xopt;
    case "FBS"
        FBS=OptiFBS(F,lambda*R);
        %FBS.OutOp=OutputOptiSNR(1,[],10);
        FBS.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);  
        FBS.ItUpOut=1;          % call OutputOpti update every ItUpOut iterations
        FBS.fista=true;         % activate fista
        FBS.maxiter=200;        % max number of iterations
        FBS.run(zeros(size(measurements)));% run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

        C = FBS.xopt;
end

end