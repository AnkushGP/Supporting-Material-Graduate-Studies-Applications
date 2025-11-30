
function [pll_opt_jit_intN,pll_opt_ugb1_intN,pll_opt_ugb2_intN,phase_margin_1,phase_margin_2] = NoiseAnalysis_PLL_IntN(single_loop,fig_num,zero_fac,pole_fac,zero_num,pole_num,fref_fac,fvco_fac,ro_usecase)
%% single_loop = 1 : Single loop PLL noise analysis;
%% single_loop = 0 : Cascaded PLL noise analysis


freq = logspace(3,9,100000);
w_rad = 2*pi.*freq;
s = 1j*w_rad;

fref = 20e6*fref_fac;
Tref = 1/fref;

fvco = 2.4e9*fvco_fac;
fvco1 = fref*20;
fvco2 = fref*120;
Tvco = 1/fvco;
Tvco1 = 1/fvco1;
Tvco2 = 1/fvco2;

Ndiv1 = fvco1/fref;
Ndiv2 = fvco2/fvco1;
Ndiv = fvco/fref;

ref_psd = -160 + 3;

Kref = (10^(ref_psd/10))*((fref/20e6)^2)*(1-ro_usecase);
Kvco = (2)*((fvco/2.4e9)^2); 
Kvco1 = 2*((fvco1/2.4e9)^2);
Kvco2 = 2*((fvco2/2.4e9)^2);


Kpole = pole_fac;
Kzero = zero_fac;
L=zero_num; %2
M=pole_num; %2

switch single_loop
    case 1
        x=1;
        for fugb = logspace(log10(fref/1000),log10(fref/5),100)   
            Wugb = 2*pi*fugb;

            LG_1 = 0;
            LG_2 = 1;
            for L_idx = 1:1:L
                LG_1 = (LG_1) + (Wugb^(L_idx)).*(1./s.^(L_idx)).*(1/Kzero^(L_idx-1));
            end
            
            for M_idx = 1:1:M
                LG_2 = (LG_2).*(1./(1 + (s./(Wugb*(Kpole^(M_idx))))));    
            end
            LG_tf = LG_1.*LG_2;
            LG_tf_dB = abs(log(abs(LG_tf)));
            LG_tf_pm_idx = find(LG_tf_dB == min(LG_tf_dB));
            phase_margin_1 = 180 + angle(LG_tf(LG_tf_pm_idx))*(180/pi);
            phase_margin_2 = "NA"
            
            pll_highPassTF = ((abs(1./(1+LG_tf))).^2);
            pll_lowPassTF = ((abs((Ndiv*LG_tf)./(1 + (LG_tf)))).^2);
            
            vco_noise_psdTotal = ((Kvco./(freq.^2)).*pll_highPassTF);

            ref_noise_psdTotal = (Kref)*(pll_lowPassTF);
            
            
            pll_noise_psdTotal = vco_noise_psdTotal + ref_noise_psdTotal;

            pll_powerTotal(x) = trapz(freq(1:length(freq)),pll_noise_psdTotal(1:length(freq)));
            ref_powerTotal(x) = trapz(freq(1:length(freq)),ref_noise_psdTotal(1:length(freq)));
            vco_powerTotal(x) = trapz(freq(1:length(freq)),vco_noise_psdTotal(1:length(freq)));

            pll_jitterTotal(x) = sqrt(pll_powerTotal(x)).*(Tvco/(2*pi));
            
            x = x+1;
        end
        
        fugb = logspace(log10(fref/1000),log10(fref/5),100);
        pll_MinPowerTotal_idx = find(pll_powerTotal == min(pll_powerTotal));
        pll_MinJitterTotal_idx = find(pll_jitterTotal == min(pll_jitterTotal));
        % find gives us the index in the vector 
        ref_cont_atMinJitter = ref_powerTotal(pll_MinPowerTotal_idx)/pll_powerTotal(pll_MinPowerTotal_idx);
        vco_cont_atMinJitter = vco_powerTotal(pll_MinPowerTotal_idx)/pll_powerTotal(pll_MinPowerTotal_idx);
        
        pll_opt_jit_intN = min(pll_jitterTotal);
        pll_opt_ugb1_intN = fugb(pll_MinJitterTotal_idx);
        pll_opt_ugb2_intN = "NA";
        
        figure(fig_num)
        loglog(fugb,1e15*pll_jitterTotal);
        hold on
        xlabel('$fu (Hz)$','interpreter','latex')
        ylabel('$Jrms (fs)$','interpreter','latex')        
    case 0
        x=1;
        for fugb2 = logspace(log10(fvco1/1000),log10(fvco1/5),100)
            Wugb2 = 2*pi*fugb2;

            LG2_1 = 0;
            LG2_2 = 1;
            for L_idx = 1:1:L
                LG2_1 = (LG2_1) + (Wugb2^(L_idx)).*(1./s.^(L_idx)).*(1/Kzero^(L_idx-1));
            end
            
            for M_idx = 1:1:M
                LG2_2 = (LG2_2).*(1./(1 + (s./(Wugb2*(Kpole^(M_idx))))));    
            end
            LG2_tf = LG2_1.*LG2_2;
            LG2_tf_dB = abs(log(abs(LG2_tf)));
            LG2_tf_pm_idx = find(LG2_tf_dB == min(LG2_tf_dB));
            phase_margin_2 = 180 + angle(LG2_tf(LG2_tf_pm_idx))*(180/pi);
            y=1;
            for fugb1 = logspace(log10(fref/1000),log10(fref/5),100)
                Wugb1 = 2*pi*fugb1;
                LG1_1 = 0;
                LG1_2 = 1;
                for L_idx = 1:1:L
                    LG1_1 = (LG1_1) + (Wugb1^(L_idx)).*(1./s.^(L_idx)).*(1/Kzero^(L_idx-1));
                end
                for M_idx = 1:1:M
                    LG1_2 = (LG1_2).*(1./(1 + (s./(Wugb1*(Kpole^(M_idx))))));    
                end
                LG1_tf = LG1_1.*LG1_2;
                LG1_tf_dB = abs(log(abs(LG1_tf)));
                LG1_tf_pm_idx = find(LG1_tf_dB == min(LG1_tf_dB));
                phase_margin_1 = 180 + angle(LG1_tf(LG1_tf_pm_idx))*(180/pi);
                
                pll_highPassTF_loop1 = ((abs(1./(1+LG1_tf))).^2);
                pll_lowPassTF_loop1 = ((abs((Ndiv1*LG1_tf)./(1 + (LG1_tf)))).^2);
                
                pll_highPassTF_loop2 = ((abs(1./(1+LG2_tf))).^2);
                pll_lowPassTF_loop2 = ((abs((Ndiv2*LG2_tf)./(1 + (LG2_tf)))).^2);
                
            
                vco_noise_psdTotal_loop1 = ((Kvco1./(freq.^2)).*pll_highPassTF_loop1);
                vco_noise_psdTotal_loop2 = ((Kvco2./(freq.^2)).*pll_highPassTF_loop2);
                
                ref_noise_psdTotal_loop1 = (Kref)*(pll_lowPassTF_loop1);

                pll_noise_psdTotal_loop1 = vco_noise_psdTotal_loop1 + ref_noise_psdTotal_loop1;
    
                pll_noise_psdTotal = (pll_noise_psdTotal_loop1.*(pll_lowPassTF_loop2)) + vco_noise_psdTotal_loop2;
                
                pll_powerTotal(x,y) = trapz(freq(1:length(freq)),pll_noise_psdTotal(1:length(freq)));
                
                pll_jitterTotal(x,y) = sqrt(pll_powerTotal(x,y)).*(Tvco2/(2*pi));
                y=y+1;
                
            end
            x=x+1;
        end
        
       pll_opt_jit_intN = min(pll_jitterTotal(:));
       [pll_MinJitterTotal_Xidx,pll_MinJitterTotal_Yidx] = find(pll_jitterTotal == pll_opt_jit_intN);
       
       fugb2 = logspace(log10(fvco1/1000),log10(fvco1/5),100);
       fugb1 = logspace(log10(fref/1000),log10(fref/5),100);
       figure(fig_num)
       mesh(fugb2,fugb1,1e15*pll_jitterTotal);
       xlabel('$f_{ugb2}(Hz)$','interpreter','latex')
       ylabel('$f_{ugb1}(Hz)$','interpreter','latex')
       zlabel('$Jitter(fs)$','interpreter','latex')
        
       pll_opt_ugb1_intN = fugb1(pll_MinJitterTotal_Yidx);
       pll_opt_ugb2_intN = fugb2(pll_MinJitterTotal_Xidx);
       
end

end
