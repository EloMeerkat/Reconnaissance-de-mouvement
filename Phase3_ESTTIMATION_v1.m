%==========================================================================
% Programme de teste
% F. Nougarou, Juillet 2017
%==========================================================================
clear all; close all; clc; system('TASKKILL /IM MyoWinform.exe >nul');clc;

%==========================================================================
load data_APP_Coef_LDA.mat
affiche =1;
%==========================================================================
val_touch = []; %pour être en mesure de faire un break lors du programme
mem_aff = 6000;
elec = 8; %8 elec par bracelet
nb_don = nb_myo*elec+1; %+1 pour le tick donnant la fréquence 1/(ticks*10^-7)
Nb_elec = nb_myo*elec;
if nb_myo == 1
    Fs = 200;
elseif nb_myo == 2
    Fs = 150;
elseif nb_myo == 3
    Fs = 100;
end
%-- Ouverture du programme ------------------------------------------------
dim = 0;
Data_EMG_mat = [];

ind = 0;
system(['MyoWinform.exe ' num2str(nb_myo) ' &']);
fileID1 = fopen('emg.dat','rb');
n = 0;
dim = 0;
while 1
    fseek(fileID1,0,'eof');
    n_dim = ftell(fileID1);
    if n_dim >= (dim + nb_don*8*20) %lit à chaque fois qu'il y a 20 valeurs pour chaque type de données
        ind = ind + 1;
        fseek(fileID1,dim,'bof');
        Data_reformat = fread(fileID1,[nb_don 20],'double');
        dim = dim + (nb_don*8*20);
        Data_file_reformat = [];
        Data_file_reformat = Data_reformat(2:end,:)';
        for elec = 1:Nb_elec
            %-- Critères calcul -------------------------------------------
            data_wind0 = []; data_wind0 = Data_file_reformat(:,elec);
            %-- Filtrage --------------------------------------------------
            data_wind = [];
            interval_frq = [20];
            Data_temp = []; Data_temp = data_wind0;
            ord = 2; Wn = interval_frq/(Fs/2);ftype = 'high';[b,a] = butter(ord,Wn,ftype);
            data_wind = filtfilt(b,a,Data_temp);
            Data_file_reformat_F(:,elec) = data_wind;
        end
        Data_EMG_mat = [Data_EMG_mat;Data_file_reformat_F];
        
        for elec = 1:Nb_elec
            Data_EMG_all(:,elec) = Data_file_reformat_F(:,elec);
            data_wind = Data_EMG_all(:,elec);
            %-- critère RMS -
            RMS(ind,elec) = norm(data_wind)/sqrt(length(data_wind));
            %-- critère MAV -
            MAV(ind,elec) = mean(abs(data_wind));
            %-- critère ZC -
            pos_PZ = []; pos_PZ = find((data_wind(1:end-1).*data_wind(2:end))<0);
            ZF(ind,elec) = length(pos_PZ);
            %-- Wavelength -
            WL(ind,elec) = sum(abs(data_wind(1:end-1,1)-data_wind(2:end,1)));
            
            SQRT_VAR(ind,elec) = sqrt(var(data_wind));
            VARi(ind,elec) = var(data_wind);
            LOG_VAR(ind,elec) = log(var(data_wind));
            
            %-- SSC -
            data_wind_pos = []; data_wind_pos = data_wind-min(data_wind);
            change_slope = []; change_slope = (data_wind_pos(2:end)-data_wind_pos(1:end-1));
            pos_CS = []; pos_CS = find((change_slope(1:end-1).*change_slope(2:end))<0);
            SSC(ind,elec) = length(pos_CS);
        end
        
        %-- Classification ------------------------------------------------
        Input_class = [];
        Input_class = [RMS(ind,:) MAV(ind,:) WL(ind,:) ZF(ind,:) SSC(ind,:)];        
        Input_LDA_app_U_norm_all(ind,:) = (Input_class-val_norm_mean(1,:))./val_norm_std;
         
        
        L_A_perSamp_Myo = ([1 Input_LDA_app_U_norm_all(ind,:)] * W_LDA_Myo');
        P_A_perSamp_Myo = zeros(1,Nb_mov);
        P_A_perSamp_Myo = exp(L_A_perSamp_Myo) ./ repmat(sum(exp(L_A_perSamp_Myo),2),[1 Nb_mov]);
        pos_max = []; val_max = []; [val_max pos_max] = max(P_A_perSamp_Myo');
        Output_LDA_Myo_app(ind,1) = (pos_max-1)';
        Output_LDA_Myo_app_dec2(ind,1) = (pos_max-1)';
        
        if ind>1
            if Output_LDA_Myo_app_dec2(ind,1) == 1
                if val_max<0.99999
                    Output_LDA_Myo_app_dec2(ind,1) = Output_LDA_Myo_app_dec2(ind-1,1);
                end
            else
                if val_max<0.99;
                    Output_LDA_Myo_app_dec2(ind,1) = Output_LDA_Myo_app_dec2(ind-1,1);
                end
            end
        end
        
        if affiche == 1
            h=figure(10);
            val_touch = [];
            val_touch = double(get(h,'CurrentCharacter')); % compare the values to the list
            if size(Output_LDA_Myo_app_dec2)<=mem_aff
                plot(Output_LDA_Myo_app_dec2(1:ind,1),'r');hold on;
                ylim([-0.1 Nb_mov+.1]);grid on;
            else
                plot(Output_LDA_Myo_app_dec2(end-mem_aff2:end,1),'r');hold on;
                hold off;
                ylim([-0.1 Nb_mov+.1]);grid on;
                xlim([0 mem_aff])
            end
            if val_touch == 27
                system('TASKKILL /IM MyoWinform.exe >nul');
                %                 close all
                break
            end
        end
    end
end
system('TASKKILL /IM MyoWinform.exe >nul');clc;
disp('== Merci! ==============================================')
