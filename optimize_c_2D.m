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
        R = LinOpGrad(H.sizein,[],'mirror');
        R = CostMixNorm21(R.sizein,numel(R.sizeout))*R;
    case "L12"
        R=CostMixNorm21([H.sizein,1],3);
    case "Tikhonov"
        R = LinOpGrad(H.sizein,[],'mirror');
        R = CostL2(H.sizein);
end
if lambda==0
    F2 = F;
    opti_type = "GradDsct";
else
    F2 = F+lambda*R; 
end
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
        GD.maxiter=2000;         % max number of iterations
        GD.run(measurements); % run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 
        
        %plot(GD.OutOp.evolcost)
        
        C = GD.xopt;
        
    case "ConjGrad"
        %Note: ConjGrad is used only with L2 cost
        if(lambda ~= 0)
            A = H.makeHtH() + lambda * LinOpIdentity(H.sizein);
        else
            A = H.makeHtH();
        end
        b = H'*measurements;
        CG=OptiConjGrad(A,b);  
        CG.OutOp=OutputOptiConjGrad(1,dot(measurements(:),measurements(:)),[],40);  
        CG.ItUpOut=2; 
        CG.maxiter=500;                             % max number of iterations
        CG.run(measurements);                                  % run the algorithm 
        
        %plot(CG.OutOp.evolcost)
        
        C = CG.xopt;
        
    case "ADMM"
        Fn={lambda*R};
        Hn = {LinOpIdentity(H.sizein)};
        rho_n=[1e-1];
        ADMM=OptiADMM(F,Fn,Hn,rho_n);
        %ADMM.OutOp=OutputOptiSNR(1,[],10);
        % STOP when the sum successives C = F*x + Fn{1}*Hn{1}*x is lower
        % than 1e-4 or when the distance between two successive step is lower than 1e-4
        ADMM.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);   
        ADMM.ItUpOut=1;             % call OutputOpti update every ItUpOut iterations
        ADMM.maxiter=200;           % max number of iterations
        ADMM.OutOp=OutputOpti(1,[],40);
        ADMM.run(zeros(size(measurements)));   % run the algorithm 
        

        %plot(ADMM.OutOp.iternum,ADMM.OutOp.evolcost,'LineWidth',1.5);
        C = ADMM.xopt;
    case "FBS"
        FBS=OptiFBS(F,lambda*R);
        %FBS.OutOp=OutputOptiSNR(1,[],10);
        FBS.CvOp=TestCvgCombine(TestCvgCostRelative(1e-4), 'StepRelative',1e-4);  
        FBS.ItUpOut=1;          % call OutputOpti update every ItUpOut iterations
        FBS.fista=true;         % activate fista
        FBS.maxiter=200;        % max number of iterations
        FBS.OutOp=OutputOpti(1,[],40);
        FBS.run(zeros(size(measurements)));% run the algorithm (Note that gam is fixed automatically to 1/F.lip here since F.lip is defined and since we do not have setted gam) 

        %plot(FBS.OutOp.evolcost)
        
        C = FBS.xopt;
end

end