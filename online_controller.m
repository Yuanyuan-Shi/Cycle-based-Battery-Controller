close all;
clear
clc
load('PJM_Reg_Signal_2013_06-201405.mat');

%% prepare regulation signal
reg = RegD_signal(:,2:2:43201); %load regulation signal, resolution=4s
num = 21600;
days = 365;
reg = reshape(reg', num*days, 1);
totalNum = length(reg);
hour = 2;
inputDim = hour*900; % 1 hour
sampleNum = floor((totalNum-inputDim)/inputDim)+1;
dataR = zeros(sampleNum,inputDim);
for i = 1: sampleNum
    dataR(i, :) = reg((i-1)*inputDim+1:i*inputDim);
end

%% parameters for the simulation
r = dataR(819,(1:900));
ts = 4/3600;
delta = 0;
[r, opts] = init_opts(r, ts, delta);

% Inputs and Parameter initialization
%   opts: data structure for parameters
%   opts.ts - control time step
%   opts.N - dispatch time interval
%   opts.lambda_p - regulation mismatch penalty price
%   opts.lambda_e - per MWh market energy price
%   opts.yita_ch - battery charging efficiency
%   opts.yita_dc - battery discharging efficiency
%   opts.soc_init - battery starting SoC level
%   opts.soc_max//opt.soc_min - battery maximum and minimum SoC level
%   opts.B_P - battery power rating 
%   opts.B_E - battery energy rating
%   Exponential battery degradation model: phi(u) = k1*(cycle depth)^k2
%   opts.k1 - cycle aging stress function linear coefficient
%   opts.k2 - cycle aging stress function nonlinear coefficient
%   r - regulatin signal
%   b - battery power output solution point, initialized as the same as the 
%   all zeros. It can be decomposed as b = b_dc - b_ch
%   soc - battery soc level
%% Online rainflow control
b = zeros(1,opts.N);
sig = zeros(1,opts.N);
eco = zeros(1,opts.N); 
limit_soc = zeros(1,opts.N);
soc = zeros(1,opts.N+1); %SoC level at the beginning of each control time step
d2_a = zeros(1,opts.N);
c2_a = zeros(1,opts.N);
d_cur = 0;

flag = 0;
for i = 1:opts.N
  % online rainflow
  if(r(i)>0) %discharging signal
    if(i==1)
      soc(1) = opts.soc_init;
      b(1) = min([opts.dc_max*opts.B_E/(opts.ts/opts.yita_dc),r(i),(opts.soc_init-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc)]);
      b(1) = max([0,b(1)]);
      d2_a(1) = 0;
      c2_a(1) = 1;
      sig(i) = r(i);
      eco(i) = opts.dc_max*opts.B_E/(opts.ts/opts.yita_dc);
      limit_soc(i) = (opts.soc_init-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc);
    else
      soc(i) = (soc(i-1)*opts.B_E-(b(i-1)/opts.yita_dc)*opts.ts)/opts.B_E;
      if(sign(r(i))*sign(r(i-1))>0||sign(r(i))*sign(r(i-1))==0)%same cycle as the previouss
          if(flag == 1)
            b(i) = 0;
            d_cur = d2_a(i-1);
          else
            ss = soc(1:i);
            ext=sig2ext(ss);
            a=rainflow(ext,1);
            d_cur = 2*a(1,end); %identify current charging half cycle depth
            b(i) = min([(opts.dc_max-d_cur)*opts.B_E/(opts.ts/opts.yita_dc), r(i),(soc(i)-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc)]);
            b(i) = max([0,b(i)]);
            if(b(i)==(soc(i)-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc)) flag=1;end
            if(b(i)== (opts.dc_max-d_cur)*opts.B_E/(opts.ts/opts.yita_dc)) flag=1;end
            if(b(i)==0) flag = 1;end
          end
      else %if regulation signal changes direction
          ss = soc(1:i);
          d_cur = 0; %a new cycle
          flag = 0;
          b(i) = min([(opts.dc_max-d_cur)*opts.B_E/(opts.ts/opts.yita_dc),r(i),(soc(i)-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc)]);
          if(b(i)==(soc(i)-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc)) flag=1;end
          if(b(i)==(opts.dc_max-d_cur)*opts.B_E/(opts.ts/opts.yita_dc)) flag=1;end
      end
      d2_a(i) = d_cur;
      c2_a(i) = 1;
      sig(i) = r(i);
      eco(i) =(opts.dc_max-d_cur)*opts.B_E/(opts.ts/opts.yita_dc);
      limit_soc(i) = (soc(i)-opts.soc_min)*opts.B_E/(opts.ts/opts.yita_dc);
    end
  
  else %charging cycle: r(t)<0 charging cycle
      if(i==1)
        soc(1) = opts.soc_init;
        b(1) = -min([opts.ch_max*opts.B_E/(opts.ts*opts.yita_ch),-r(i),(opts.soc_max-opts.soc_init)*opts.B_E/(opts.ts*opts.yita_ch)]);
        d2_a(1) = 0;
        c2_a(1) = -1;
        sig(i) = r(i);
        eco(i) = -opts.ch_max*opts.B_E/(opts.ts*opts.yita_ch);
        limit_soc(i) = -(opts.soc_max-opts.soc_init)*opts.B_E/(opts.ts*opts.yita_ch);
      else
          soc(i) = (soc(i-1)*opts.B_E-b(i-1)*opts.yita_ch*opts.ts)/opts.B_E;
          if(sign(r(i))*sign(r(i-1))>0||sign(r(i))*sign(r(i-1))==0)
              if(flag == 1)%Hit the SoC/Eco boundary in one cycle
                b(i) = 0;
                d_cur = d2_a(i-1);
              else
                ss = soc(1:i);
                ext=sig2ext(ss);
                a=rainflow(ext,1);
                d_cur = 2*a(1,end);%identify current charging half cycle depth
                b(i) = min([(opts.ch_max-d_cur)*opts.B_E/(opts.ts*opts.yita_ch), -r(i),(opts.soc_max-soc(i))*opts.B_E/(opts.ts*opts.yita_ch)]);
                b(i) = max([b(i),0]);
                b(i) = -b(i);
                if(-b(i)== max([(opts.ch_max-d_cur)*opts.B_E/(opts.ts*opts.yita_ch),0])) flag=1; end
                if(-b(i)==(opts.soc_max-soc(i))*opts.B_E/(opts.ts*opts.yita_ch)) flag=1; end
                if(b(i)==0) flag = 1;end
              end
          else %if regulation signal changes direction
              ss = soc(1:i);
              d_cur = 0; %a new cycle
              flag = 0;
              b(i) = -min([(opts.ch_max-d_cur)*opts.B_E/(opts.ts*opts.yita_ch),-r(i),(opts.soc_max-soc(i))*opts.B_E/(opts.ts*opts.yita_ch)]);
              if(-b(i)==(opts.ch_max-d_cur)*opts.B_E/(opts.ts*opts.yita_ch)) flag=1; end
              if(-b(i)==(opts.soc_max-soc(i))*opts.B_E/(opts.ts*opts.yita_ch)) flag=1; end
          end
          d2_a(i) = d_cur;
          c2_a(i) = -1;
          sig(i) = r(i);
          eco(i) = -(opts.ch_max-d_cur)*opts.B_E/(opts.ts*opts.yita_ch);
          limit_soc(i) = -(opts.soc_max-soc(i))*opts.B_E/(opts.ts*opts.yita_ch);
      end
  end
end


%% Calculate cost
b_dc(1,:) = (b(1,:)+abs(b(1,:)))/2; %discharge is positive
b_ch(1,:) = b_dc(1,:)-b(1,:); %charging component;
soc1 = cal_soc(b, opts.soc_init, opts.yita_ch, opts.yita_dc, opts.B_E, opts.ts);
C = cal_cost(b_ch, b_dc, r', soc1, opts);
%%
figure
subplot(4,1,1)
hold all
plot(b,'-or');
plot(sig,'-b');
% plot(eco,'--g');
legend('bat','sig');
hold off
title(sprintf('Battery control Output, Cost: %.1f',C.value));
subplot(4,1,2)
plot(soc(1:opts.N))
title('SoC')
subplot(4,1,3)
stem(d2_a.*c2_a)
title('Cycle Depth');
subplot(4,1,4)
stem(c2_a)
title('Charging or discharging Cycle');



