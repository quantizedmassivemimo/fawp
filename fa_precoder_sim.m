% =========================================================================
% -- Simulator for Finite-Alphabet Wiener Filter Precoding (FAWP)
% -------------------------------------------------------------------------
% -- Simulator for the paper:
% -- Oscar Castañeda, Sven Jacobsson, Giuseppe Durisi, Tom Goldstein, and
% -- Christoph Studer, "Finite-Alphabet Wiener Filter Precoding for mmWave
% -- Massive MU-MIMO Systems," Asilomar Conference on Signals, Systems,
% -- and Computers, November 2019, pp. 178-183
% -------------------------------------------------------------------------
% -- (c) 2019 Oscar Castañeda, Christoph Studer, and Sven Jacobsson
% -- e-mail: caoscar@ethz.ch, studer@ethz.ch, and sven.jacobsson@ericsson.com
% =========================================================================

function fa_precoder_sim(varargin)

% -- set up default/custom parameters

if isempty(varargin)

    disp('using default simulation settings and parameters...')

    % set default simulation parameters
    par.runId = 0;              % simulation ID (used to reproduce results)    
    par.U = 16;                 % number of single-antenna users
    par.B = 256;                % number of base-station antennas (B>>U)
    par.mod = '16QAM';          % modulation type: 'BPSK','QPSK','16QAM','64QAM','8PSK'
    par.trials = 1e4;           % number of Monte-Carlo trials (transmissions)
    par.NTPdB_list = ...        % list of normalized transmit power [dB]
        0:2:20;                 % values to be simulated    
    par.betaEst = 'pilot';      % method for estimating the precoding factor: 'perfect', 'pilot'
    par.pilot = 'sqrtEs';       % use square root of par.Es as pilot;
                                % a specific symbol could also be used, but
                                % be careful not to exceed TxPower constraint
    par.numPilots = 1;          % number of slots to use for pilots for betaEst;
                                % used only if par.betaEst='pilot'
    par.precoder = ...          % precoding scheme(s) to be evaluated:
        ...                     % 'WF', 'Pre-FAWP-WF', 'Post-FAWP-WF',
        ...                     % 'Pre-FAWP-FBS', 'Post-FAWP-FBS'
        ...                     % if you are using 'WF' and any other of its
        ...                     % variants, WF must be computed first!
        {'WF','Pre-FAWP-WF','Post-FAWP-WF','Pre-FAWP-FBS','Post-FAWP-FBS'};
    par.save = true;            % save results
    par.plot = true;            % plot results
    par.plotEVM = false;        % plot EVM (true) or BER (false)
    par.rhopower = 1;           % power constraint is rhopower^2
    par.FAWP.levels = 2^1;      % number of levels to use in finite-alphabet
                                % matrices. A value of 2^B means that B bits are
                                % being completely used
    par.FAWPFBS.WFinit = false; % initialize FAWP-FBS with the FAWP-WF solution,
                                % otherwise use MRT solution
    % manually indicate FAWP-FBS parameters for those scenarios that were
    % not explored in the paper; note that parameter tuning is required to
    % obtain the best performance
    % -- Pre-FAWP-FBS
    par.PreFAWPFBS.iters_fixed = 10;
    par.PreFAWPFBS.opt_fixed = 1.10;
    par.PreFAWPFBS.tau_fixed = 2^-8;
    par.PreFAWPFBS.push_fixed = 1.25;
    % -- Post-FAWP-FBS
    par.PostFAWPFBS.iters_fixed = 10;
    par.PostFAWPFBS.opt_fixed = 19;
    par.PostFAWPFBS.tau_fixed = 2^-9;
    par.PostFAWPFBS.push_fixed = 1.25;

else

    disp('use custom simulation settings and parameters...')
    par = varargin{1};       % only argument is par structure

end

% -- initialization

% make sure par.plotEVM only has effect if par.plot is set
par.plotEVM = par.plotEVM && par.plot;

% total number of time slots used with the same channel, for both beta
% estimation and payload transmission
par.dataSlots = 1; % number of time slots for data, only 1 is supported for now!
par.timeSlots = par.dataSlots + strcmp(par.betaEst,'pilot')*par.numPilots;

% make sure at least one pilot is specified
if(strcmp(par.betaEst,'pilot')&&(par.numPilots<1))
    error('If you are training with pilots, there should be at least one pilot')
end

% if Pre-FAWP-FBS is running with par.FAWPFBS.WFinit, make sure that
% Pre-FAWP-WF is also being computed; otherwise, set variables to zero so
% that execution continues
if( ~any(strcmp(par.precoder,'Pre-FAWP-WF')) && ...
        any(strcmp(par.precoder,'Pre-FAWP-FBS')) && par.FAWPFBS.WFinit )
    par.precoder(2:end+1) = par.precoder;
    par.precoder{1} = 'Pre-FAWP-WF';
else
    WFPre_X = 0;
end
% same for Post-FAWP-FBS
if( ~any(strcmp(par.precoder,'Post-FAWP-WF')) && ...
        any(strcmp(par.precoder,'Post-FAWP-FBS')) && par.FAWPFBS.WFinit )
    par.precoder(2:end+1) = par.precoder;
    par.precoder{1} = 'Post-FAWP-WF';
else
    WFPost_V = 0;
end

% if any of the WF versions is being computed, make sure that WF is running
if( ~any(strcmp(par.precoder,'WF')) && ...
        (any(contains(par.precoder,'WF')) ) )
    par.precoder(2:end+1) = par.precoder;
    par.precoder{1} = 'WF';
end

% load parameters for Pre-FAWP-FBS
if(any(strcmp(par.precoder,'Pre-FAWP-FBS')))
    PreFAWPFBSidx = log2(par.FAWP.levels);
    PreFAWPFBSiters = [10 5 5];
    par.PreFAWPFBS.iters = PreFAWPFBSiters(PreFAWPFBSidx);
    PreFAWPFBSfile = ['./params/PreFAWPFBS_Rayleigh_U' int2str(par.U) ...
        '_B' int2str(par.B) '_Levels' int2str(par.FAWP.levels) ...
        '_Iters' int2str(par.PreFAWPFBS.iters) '.mat'];
    if isfile(PreFAWPFBSfile)
        PreFAWPFBSparams = load(PreFAWPFBSfile);
        par.PreFAWPFBS.opt = PreFAWPFBSparams.opt;
        par.PreFAWPFBS.tau = PreFAWPFBSparams.tau;
        par.PreFAWPFBS.push = PreFAWPFBSparams.push;
    else
        warn_msg = [ 'No trained Pre-FAWP-FBS parameters available for ' ...
                     'this system and Pre-FAWP-FBS configuration. Using ' ...
                     'parameters in par.PreFAWPFBS. Note that parameter ' ...
                     'tuning is needed to obtain the best performance.' ];
        warning(warn_msg);
        par.PreFAWPFBS.iters = par.PreFAWPFBS.iters_fixed;
        par.PreFAWPFBS.opt = par.PreFAWPFBS.opt_fixed * ones(par.PreFAWPFBS.iters,1);
        par.PreFAWPFBS.tau = par.PreFAWPFBS.tau_fixed * ones(par.PreFAWPFBS.iters,1);
        par.PreFAWPFBS.push = par.PreFAWPFBS.push_fixed * ones(par.PreFAWPFBS.iters,1);
    end
end

% parameters for Post-FAWP-FBS
if(any(strcmp(par.precoder,'Post-FAWP-FBS')))
    if (par.B==256)&&(par.U==16)&&(par.FAWP.levels<=2^3)
        PostFAWPFBSidx = log2(par.FAWP.levels);
        PostFAWPFBSiters = [10 10 10];
        PostFAWPFBSopt = [19 17 17];
        PostFAWPFBStau = [2^-9 2^-9 2^-9];    
        PostFAWPFBSpush = [1.25 1.05 1.01];
        par.PostFAWPFBS.iters = PostFAWPFBSiters(PostFAWPFBSidx);
        par.PostFAWPFBS.opt = PostFAWPFBSopt(PostFAWPFBSidx)*ones(par.PostFAWPFBS.iters,1);
        par.PostFAWPFBS.tau = PostFAWPFBStau(PostFAWPFBSidx)*ones(par.PostFAWPFBS.iters,1);
        par.PostFAWPFBS.push = PostFAWPFBSpush(PostFAWPFBSidx)*ones(par.PostFAWPFBS.iters,1);
    else
        warn_msg = [ 'No trained Post-FAWP-FBS parameters available for ' ...
                     'this system and Post-FAWP-FBS configuration. Using ' ...
                     'parameters in par.PostFAWPFBS. Note that parameter ' ...
                     'tuning is needed to obtain the best performance.' ];
        warning(warn_msg);
        par.PostFAWPFBS.iters = par.PostFAWPFBS.iters_fixed;
        par.PostFAWPFBS.opt = par.PostFAWPFBS.opt_fixed * ones(par.PostFAWPFBS.iters,1);
        par.PostFAWPFBS.tau = par.PostFAWPFBS.tau_fixed * ones(par.PostFAWPFBS.iters,1);
        par.PostFAWPFBS.push = par.PostFAWPFBS.push_fixed * ones(par.PostFAWPFBS.iters,1);
    end
end

% use runId random seed (enables reproducibility)
rng(par.runId);

% simulation name (used for saving results) if not indicated
if(~isfield(par,'simName'))
    par.simName = ['Rayleigh_U',num2str(par.U),'xB',num2str(par.B),...
        '_',par.mod, '_', num2str(par.trials),'Trials'];
end

% set up Gray-mapped constellation alphabet (according to IEEE 802.11)
switch (par.mod)
    case 'BPSK'
        par.symbols = [ -1 1 ];
    case 'QPSK'
        par.symbols = [ -1-1i,-1+1i,+1-1i,+1+1i ];
    case '16QAM'
        par.symbols = [ -3-3i,-3-1i,-3+3i,-3+1i, ...
            -1-3i,-1-1i,-1+3i,-1+1i, ...
            +3-3i,+3-1i,+3+3i,+3+1i, ...
            +1-3i,+1-1i,+1+3i,+1+1i ];
    case '64QAM'
        par.symbols = [ -7-7i,-7-5i,-7-1i,-7-3i,-7+7i,-7+5i,-7+1i,-7+3i, ...
            -5-7i,-5-5i,-5-1i,-5-3i,-5+7i,-5+5i,-5+1i,-5+3i, ...
            -1-7i,-1-5i,-1-1i,-1-3i,-1+7i,-1+5i,-1+1i,-1+3i, ...
            -3-7i,-3-5i,-3-1i,-3-3i,-3+7i,-3+5i,-3+1i,-3+3i, ...
            +7-7i,+7-5i,+7-1i,+7-3i,+7+7i,+7+5i,+7+1i,+7+3i, ...
            +5-7i,+5-5i,+5-1i,+5-3i,+5+7i,+5+5i,+5+1i,+5+3i, ...
            +1-7i,+1-5i,+1-1i,+1-3i,+1+7i,+1+5i,+1+1i,+1+3i, ...
            +3-7i,+3-5i,+3-1i,+3-3i,+3+7i,+3+5i,+3+1i,+3+3i ];
    case '8PSK'
        par.symbols = [ exp(1i*2*pi/8*0), exp(1i*2*pi/8*1), ...
            exp(1i*2*pi/8*7), exp(1i*2*pi/8*6), ...
            exp(1i*2*pi/8*3), exp(1i*2*pi/8*2), ...
            exp(1i*2*pi/8*4), exp(1i*2*pi/8*5) ];
end

% compute symbol energy
par.Es = mean(abs(par.symbols).^2);

% precompute bit labels
par.bps = log2(length(par.symbols)); % number of bits per symbol
par.bits = de2bi(0:length(par.symbols)-1,par.bps,'left-msb');

% track simulation time
time_elapsed = 0;

% -- start simulation

% - initialize result arrays (detector x normalized transmit power)
% vector error rate
res.VER = zeros(length(par.precoder),length(par.NTPdB_list));
% symbol error rate
res.SER = zeros(length(par.precoder),length(par.NTPdB_list));
% bit error rate
res.BER = zeros(length(par.precoder),length(par.NTPdB_list));
% denominator of EVM
res.VM = zeros(length(par.precoder),length(par.NTPdB_list));
% mean-square error (but note that the channel is not fixed)
res.MSE = zeros(length(par.precoder),length(par.NTPdB_list));
% SINDR
res.SINDR = zeros(length(par.precoder),length(par.NTPdB_list));
% transmit power
res.TxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% receive power
res.RxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% pilot transmit power
res.PilotTxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% pilot receive power
res.PilotRxPower = zeros(length(par.precoder),length(par.NTPdB_list));
% simulation time
res.time = zeros(length(par.precoder),length(par.NTPdB_list));

% compute noise variances to be considered
N0_list = (par.rhopower^2)*(10.^(-par.NTPdB_list/10));

% generate random bit stream (antenna x bit x trial)
bits = randi([0 1],par.U,par.bps,par.trials);

% trials loop
tic
for tt=1:par.trials

    % generate transmit symbol
    if(par.plotEVM)
        % we test with Gaussian numbers for EVM
        s = sqrt(0.5*par.Es)*(randn(par.U,1)+1i*randn(par.U,1));
    else
        idx = bi2de(bits(:,:,tt),'left-msb')+1;
        s = par.symbols(idx).';
    end

    % generate iid Gaussian channel matrix and noise vector
    % For a fixed channel, we will simulate tranmission over several slots
    n = sqrt(0.5)*(randn(par.U,par.timeSlots)+1i*randn(par.U,par.timeSlots));
    H = sqrt(0.5)*(randn(par.U,par.B)+1i*randn(par.U,par.B));

    % normalized transmit power loop
    for kk=1:length(par.NTPdB_list)

        % set noise variance
        N0 = N0_list(kk);

        % algorithm loop
        for dd=1:length(par.precoder)

            % record time used by the precoder
            starttime = toc;

            % precoders
            switch (par.precoder{dd})
                case 'WF'    % Wiener-Filter precoding (infinite precision)
                    [P, WF_Q] = WF(par, H, N0);
                case 'Pre-FAWP-WF'
                    [P, WFPre_X] = PreFAWPWF(par, H, N0, WF_Q);
                case 'Post-FAWP-WF'
                    [P, WFPost_V] = PostFAWPWF(par, H, N0, WF_Q);
                case 'Pre-FAWP-FBS'
                    P = PreFAWPFBS(par, H, N0, WFPre_X);
                case 'Post-FAWP-FBS'
                    P = PostFAWPFBS(par, H, N0, WFPost_V);
                otherwise
                    error('par.precoder not specified')
            end

            % record beamforming simulation time
            res.time(dd,kk) = res.time(dd,kk) + (toc-starttime);

            % pilot time slots loop ---------------------------------------
            if(strcmp(par.betaEst,'pilot'))
                yPilots = zeros(par.U,par.numPilots);
            end
            for pts=1:strcmp(par.betaEst,'pilot')*par.numPilots
                % generate pilot
                switch (par.pilot)
                    case 'sqrtEs'
                        pilotS = sqrt(par.Es);
                    otherwise
                        pilotS = par.pilot;
                end
                pilotx = P*(pilotS*ones(par.U,1));

                % transmit pilot
                Hpx = H*pilotx;
                yPilots(:,pts) = Hpx + sqrt(N0)*n(:,pts);

                % extract transmit and receive power
                res.PilotTxPower(dd,kk) = res.PilotTxPower(dd,kk) + ...
                    mean(sum(abs(pilotx).^2));
                res.PilotRxPower(dd,kk) = res.PilotRxPower(dd,kk) + ...
                    mean(sum(abs(Hpx).^2))/par.U;
            end % pilot time slots loop -----------------------------------

            % estimating beta ---------------------------------------------
            % we incorporate the prior knowledge that beta is positive real
            switch (par.betaEst)
                case 'perfect'
                    beta = real(1./diag(H*P));
                case 'pilot'
                    beta = pilotS./real(mean(yPilots,2));
                    beta = max(0.1,beta); % so that beta is positive
                otherwise
                    error('par.betaEst not specified')
            end % estimating beta -----------------------------------------

            % data transmission time slots loop ---------------------------
            for dts=1:par.dataSlots

                % transmit data over noisy channel
                x = P*s;
                Hx = H*x;
                y = Hx + sqrt(N0)*...
                    n(:,dts+strcmp(par.betaEst,'pilot')*par.numPilots);

                % extract transmit and receive power
                res.TxPower(dd,kk) = res.TxPower(dd,kk) + ...
                    mean(sum(abs(x).^2));
                res.RxPower(dd,kk) = res.RxPower(dd,kk) + ...
                    mean(sum(abs(Hx).^2))/par.U;

                % UEs scale with the precoding factor beta
                shat = beta.*y;

                if(~par.plotEVM)

                    % -- perform UE-side detection
                    [~,idxhat] = min(abs(shat*ones(1,length(par.symbols)) ...
                        -ones(par.U,1)*par.symbols).^2,[],2);
                    bithat = par.bits(idxhat,:);

                    % -- compute error and complexity metrics
                    err = (idx~=idxhat);
                    res.VER(dd,kk) = res.VER(dd,kk) + any(err);
                    res.SER(dd,kk) = res.SER(dd,kk) + sum(err)/par.U;
                    res.BER(dd,kk) = res.BER(dd,kk) + ...
                        sum(sum(bits(:,:,tt)~=bithat))/(par.U*par.bps);

                end
                res.MSE(dd,kk) = res.MSE(dd,kk) + norm(shat - s)^2;
                res.VM(dd,kk) = res.VM(dd,kk) + norm(s)^2;
                res.SINDR(dd,kk) = res.SINDR(dd,kk) + norm(s)^2/norm(shat - s)^2;

            end % data tx time slots loop ---------------------------------

        end % algorithm loop

    end % NTP loop

    % keep track of simulation time
    if toc>10
        time=toc;
        time_elapsed = time_elapsed + time;
        fprintf('estimated remaining simulation time: %3.0f min.\n',...
            time_elapsed*(par.trials/tt-1)/60);
        tic
    end

end % trials loop

% normalize results
res.VER = res.VER/par.trials;
res.SER = res.SER/par.trials;
res.BER = res.BER/par.trials;
res.EVM = sqrt(res.MSE./res.VM).*100;
res.MSE = res.MSE/par.trials;
res.SINDR = res.SINDR/par.trials;
res.TxPower = res.TxPower/par.trials;
res.RxPower = res.RxPower/par.trials;
res.PilotTxPower = res.PilotTxPower/par.trials;
res.PilotRxPower = res.PilotRxPower/par.trials;
res.time = res.time/par.trials;
res.time_elapsed = time_elapsed;

% -- save final results (par and res structures)

if par.save
    [~,~]=mkdir('./results');
    save(['./results/' par.simName '_' num2str(par.runId) ],'par','res');
end

% -- show results (generates fairly nice Matlab plots)

if par.plot

    % - BER results
    style_color = {'black','#D95319','#77AC30','#0072BD','#7E2F8E','#EDB120','#4DBEEE'};
    style_line = {'-','--','-.','--',':','-.',':'};
    style_marker = {'none','square','^','o','diamond','x','*'};
    figure()
    for dd=1:length(par.precoder)
        if(par.plotEVM)
            plot(par.NTPdB_list,res.EVM(dd,:), ...
                 'Color',style_color{dd},'LineStyle',style_line{dd}, ...
                 'Marker',style_marker{dd},'LineWidth',2);
        else
            semilogy(par.NTPdB_list,res.BER(dd,:), ...
                 'Color',style_color{dd},'LineStyle',style_line{dd}, ...
                 'Marker',style_marker{dd},'LineWidth',2);
        end
        if (dd==1)
            hold on
        end
    end
    if(par.plotEVM)
        % Mark minimum EVM (%) value required by 3GPP-5G NR in section
        % 6.5.2.2
        minEVMmods = {'QPSK','16-QAM','64-QAM','256-QAM'};
        minEVM = [17.5, 12.5, 8, 3.5]; % {QPSK, 16QAM, 64QAM, 256QAM}
        for dd=1:length(minEVMmods)
            plot(par.NTPdB_list,minEVM(dd)*ones(size(par.NTPdB_list)),'r--','LineWidth',2);
            text(min(par.NTPdB_list), minEVM(dd)+1, minEVMmods{dd},'Color','red');
        end
    end
    hold off
    grid on
    box on
    xlabel('normalized transmit power [dB]','FontSize',12)
    if(par.plotEVM)
        ylabel('error-vector magnitude (EVM) [%]','FontSize',12);
    else
        ylabel('uncoded bit error rate (BER)','FontSize',12);
    end
    if length(par.NTPdB_list) > 1
        if(par.plotEVM)
            axis([min(par.NTPdB_list) max(par.NTPdB_list) 1 35]);
        else
            axis([min(par.NTPdB_list) max(par.NTPdB_list) 1e-4 1]);
        end
    end
    if(par.plotEVM)
        legend(par.precoder,'FontSize',12,'location','northeast')
    else
        legend(par.precoder,'FontSize',12,'location','southwest')
    end
    plot_title = [int2str(log2(par.FAWP.levels)) '-bit FAWP; beta estimation: ', par.betaEst];
    title(plot_title)
    set(gca,'FontSize',12);

end

end

%% WF:
%  Wiener-Filter precoding (infinite precision)

function [P, Q] = WF(par, H, N0)

% precoding matrix components
k = par.U*N0/(par.rhopower)^2;
Q = (H'*H+k*eye(par.B))\H';

% beamforming factor
betaTx = sqrt(par.Es*trace(Q'*Q)/(par.rhopower)^2);

% precoding matrix
P = Q/betaTx;

end

%% Pre-FAWP-WF:
%  Quantize the WF solution per column to obtain the low-resolution matrix
%  A of a pre-FAWP matrix

function [P, X] = PreFAWPWF(par, H, N0, WF_Q)

% precoding matrix components
k = par.U*N0/(par.rhopower)^2;

% quantize WF per column
% -- normalize WF per column
WF_QR = real(WF_Q);
WF_QI = imag(WF_Q);
max_abs_W = max(abs([WF_QR; WF_QI]));
max_abs_W = ones(par.B,1) * max_abs_W;
WR = WF_QR./max_abs_W;
WI = WF_QI./max_abs_W;
W = WR + 1i*WI;
% -- quantize the normalized WF
X = UniSymQuantiz(W,par.FAWP.levels);

% alpha
alpha_num = diag(X'*H').';
alpha_den = sum(abs(H*X).^2,1)+k*sum(abs(X).^2,1);
alpha = alpha_num./alpha_den;

% finite-alphabet matrix
Q = X*diag(alpha);

% beamforming factor
betaTx = sqrt(par.Es*trace(Q'*Q)/(par.rhopower)^2);

% precoding matrix
P = Q/betaTx;

end

%% Post-FAWP-WF:
%  Quantize the WF solution per row to obtain the low-resolution matrix Z
%  of a post-FAWP matrix

function [P, V] = PostFAWPWF(par, H, N0, WF_Q)

% precoding matrix components
k = par.U*N0/(par.rhopower)^2;

% quantize WF per row
% -- normalize WF per row
WF_QR = real(WF_Q);
WF_QI = imag(WF_Q);
max_abs_W = max(abs([WF_QR WF_QI]),[],2);
max_abs_W = max_abs_W * ones(1,par.U);
WR = WF_QR./max_abs_W;
WI = WF_QI./max_abs_W;
W = WR + 1i*WI;
% -- quantize the normalized WF
V = UniSymQuantiz(W,par.FAWP.levels);

% zeta
zeta_num = diag(H'*V').';
zeta_den = sum(abs(H'*V').^2,1)+k*sum(abs(V').^2,1);
zeta = zeta_num./zeta_den;

% finite-alphabet matrix
Q = diag(zeta)*V;

% beamforming factor
betaTx = sqrt(par.Es*trace(Q'*Q)/(par.rhopower)^2);

% precoding matrix
P = Q/betaTx;

end

%% Pre-FAWP-FBS:
%  Use FBS on a per-column basis to obtain the low-resolution matrix A of a
%  pre-FAWP matrix

function P = PreFAWPFBS(par, H, N0, WFPre_X)

% parameters ------------------
iters = par.PreFAWPFBS.iters;
tau = par.PreFAWPFBS.tau;
opt = par.PreFAWPFBS.opt;
push = par.PreFAWPFBS.push;
%-------------------------------

k = par.U*N0/(par.rhopower)^2;

if(~par.FAWPFBS.WFinit)
    X = H'; % MRT initializer
else
    % Pre-FAWP-WF: We need to scale it so that it is within -1 and 1
    % and so that it gets quantized correctly if no iterations are done
    maxlevelWF = (par.FAWP.levels-1)/(1+mod(par.FAWP.levels,2)); % absolute value of max level in WF
    maxlevelFAWP = 1-1/par.FAWP.levels; % absolute value of max level in FAWP-FBS
    X = (WFPre_X./maxlevelWF).*maxlevelFAWP;
end

for ii=1:iters
    HX = H*X;
    gradF = H'*(HX-opt(ii)*diag(diag(HX)));
    Z = X-tau(ii)*gradF;
    Z = push(ii)*Z;
    X = min(max(real(Z),-1),1) + 1i*min(max(imag(Z),-1),1);
end

% quantize X per column, with 1 as max value
XQ = UniSymQuantiz(X,par.FAWP.levels);

% alpha
alpha_num = diag(XQ'*H').';
alpha_den = sum(abs(H*XQ).^2,1)+k*sum(abs(XQ).^2,1);
alpha = alpha_num./alpha_den;

% finite-alphabet matrix
Q = XQ*diag(alpha);

% beamforming factor
betaTx = sqrt(par.Es*trace(Q'*Q)/(par.rhopower)^2);

% precoding matrix
P = Q/betaTx;

end

%% Post-FAWP-FBS:
%  Use FBS on a per-row basis to obtain the low-resolution matrix Z of a
%  post-FAWP matrix

function P = PostFAWPFBS(par, H, N0, WFPost_V)

% parameters ------------------
iters = par.PostFAWPFBS.iters;
tau = par.PostFAWPFBS.tau;
opt = par.PostFAWPFBS.opt;
push = par.PostFAWPFBS.push;
%-------------------------------

k = par.U*N0/(par.rhopower)^2;

if(~par.FAWPFBS.WFinit)
    VH = H; % MRT initializer
else
    % Post-FAWP-WF: We need to scale it so that it is within -1 and 1
    % and so that it gets quantized correctly if no iterations are done
    maxlevelWF = (par.FAWP.levels-1)/(1+mod(par.FAWP.levels,2)); % absolute value of max level in WF
    maxlevelFAWP = 1-1/par.FAWP.levels; % absolute value of max level in FAWP-FBS
    VH = ((WFPost_V')./maxlevelWF).*maxlevelFAWP;
end

for ii=1:iters
    HVH = H'*VH;
    gradF = H*(HVH-opt(ii)*diag(diag(HVH)));
    Z = VH-tau(ii)*gradF;
    Z = push(ii)*Z;
    VH = min(max(real(Z),-1),1) + 1i*min(max(imag(Z),-1),1);
end
V = VH';

% quantize V
VQ = UniSymQuantiz(V,par.FAWP.levels);

% zeta
zeta_num = diag(H'*VQ').';
zeta_den = sum(abs(H'*VQ').^2,1)+k*sum(abs(VQ').^2,1);
zeta = zeta_num./zeta_den;

% finite-alphabet matrix
Q = diag(zeta)*VQ;

% beamforming factor
betaTx = sqrt(par.Es*trace(Q'*Q)/(par.rhopower)^2);

% precoding matrix
P = Q/betaTx;

end

%% Uniform Symmetric Quantizer between -1 and +1

function [XQ] = UniSymQuantiz(X,numLevels)

XR = real(X);
XI = imag(X);
% - first, saturate the extreme bins
scale = 1-1/numLevels;
XR = XR/scale;
XI = XI/scale;
XR(abs(XR)>1) = sign(XR(abs(XR)>1));
XI(abs(XI)>1) = sign(XI(abs(XI)>1));
% - map from [-1,+1] range to [0,1] to quantize with rounding function    
XR = round(0.5*(numLevels-1)*(XR+1));
XI = round(0.5*(numLevels-1)*(XI+1));
% - finish quantization by returning from [0,1] range to [-1,+1]
XR = (2/(numLevels-1))*XR-1;
XI = (2/(numLevels-1))*XI-1;
XQ = XR + 1i*XI;

end
