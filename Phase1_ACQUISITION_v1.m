%==========================================================================
% Programme d'acquistion de données
%  - acquision séparées pour les APP mouvement par mouvement avec phase
%    de de force modérée et force maximale
%
% E. Vannier, décembre 2017
%==========================================================================
clear all; close all; clc; system('TASKKILL /IM MyoWinfor m.exe >nul');clc;

%======================================================== ==================
Mat_mov = [1:5];
Nb_mov  = length(Mat_mov);
Nb_ample = 100;
Nb_bloc = 20;
Nb_iter_mov_app = (Nb_ample)*ones(Nb_mov,1);
mem_aff = 200;
affiche = 1;
val_version = '29_novembre';
%==========================================================================
TH_Max = [];
val_touch = []; %pour être en mesure de faire un break lors du programme
nb_myo = 1; %définit le nombre de  bracelet EMG's(1,2 ou 3)
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
% system(['MyoWinform.exe ' num2str(nb_myo) ' &']);
% fileID1 = fopen('emg.dat','rb');

dim = 0;
ind = 0;
target_charact_all = [];
Data_EMG_mat = [];

tmp_mov = 0;

for b = 2:2
    disp(sprintf('==> Acquisition des données sur le bras %d ! ===========================',b));
    for m = 16:Nb_bloc
        tmp_mov=0;
        for n = Mat_mov
            tmp_mov = tmp_mov + 1;
            disp(sprintf('==> Mouv. %d == Press. to GO! ===========================',n));
            target_charact = [];
            Features_Bras =[];
            Data_EMG_mat = [];
            ind = 0;
            
            if m < 10
                name_file = sprintf('Bloc_0%d_Mouv%d_bras%d.mat',m,n,b);
            else
                name_file = sprintf('Bloc_%d_Mouv%d_bras%d.mat',m,n,b);
                
            end
            
            pause
            system(['MyoWinform.exe ' num2str(nb_myo) ' &']);
            fileID1 = fopen('emg.dat','rb');
            dim = 0;
            
            RMS = []; MAV= [];ZF = [];WL = []; SSC = [];
            while 1
                fseek(fileID1,0,'eof');
                n_dim = ftell(fileID1);
                if n_dim >= (dim + nb_don*8*20) %lit à chaque fois qu'il y a 20 valeurs pour chaque type de données
                    ind = ind + 1;
                    fseek(fileID1,dim,'bof');
                    Data_reformat = fread(fileID1,[nb_don 20],'double');
                    dim = dim + (nb_don*8*20);
                    Data_file_reformat = [];
                    
                    %--> ici
                    Data_file_reformat = Data_reformat(2:end,:)';
                    %-- filtrage --------------------------------------------------
                    for elec = 1:Nb_elec
                        %-- Critères calcul -------------------------------------------
                        data_wind0 = []; data_wind0 = Data_file_reformat(:,elec);
                        %-- Filtrage --------------------------------------------------
                        data_wind = [];
                        interval_frq = [20]; 
                        Data_temp = []; Data_temp = data_wind0;
                        ord = 2; Wn = interval_frq/(Fs/2);ftype = 'high';[b1,a] = butter(ord,Wn,ftype);
                        data_wind = filtfilt(b1,a,Data_temp);
                        Data_file_reformat_F(:,elec) = data_wind;
                    end
                    Data_EMG_mat = [Data_EMG_mat;Data_file_reformat_F];
                    
                    %-- Calcul des critères discriminants -------------------------
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
                        %-- SSC -
                        data_wind_pos = []; data_wind_pos = data_wind-min(data_wind);
                        change_slope = []; change_slope = (data_wind_pos(2:end)-data_wind_pos(1:end-1));
                        pos_CS = []; pos_CS = find((change_slope(1:end-1).*change_slope(2:end))<0);
                        SSC(ind,elec) = length(pos_CS);
                    end
                    
                    RMS_med(ind,1) = median(RMS(ind,:));
                    
                    if affiche == 1
                        h = figure(1);
                        
                        val_touch = double(get(h,'CurrentCharacter'));
                        target_charact(ind,1) = 1;
                        
                        if size(Data_EMG_mat)<=mem_aff
                            subplot(2,1,1)
                            plot(RMS_med(1:end,1),'r');hold off
                            subplot(2,1,2)
                            plot(mean(Data_EMG_mat(:,:)'));hold off;
                            ylim([-129 129]);
                        else
                            subplot(2,1,1)
                            plot(RMS_med(2:end,1),'r');hold off
                            subplot(2,1,2)
                            plot(mean(Data_EMG_mat(end-mem_aff:end,:)'));hold off;
                            ylim([-129 129]);
                            xlim([0 mem_aff]);
                        end
                        val_touch = double(get(h,'CurrentCharacter'));
                        
                        if val_touch == 27
                            tmp_mov = tmp_mov - 1;
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            close all;
                            break
                        elseif ind == Nb_iter_mov_app(tmp_mov,1)
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            system('TASKKILL /IM MyoWinform.exe >nul');clc;
                            disp('Sauvegarde');
                            
                            Features_Bras = [RMS MAV WL ZF SSC];
                            Target_Bras = (n-1)*ones(length(target_charact),1);
                            save(name_file,'Features_Bras','Target_Bras');
                            break
                            
                        end
                        val_touch = [];
                        pause(0.005)
                    end
                end
            end
        system('TASKKILL /IM MyoWinform.exe >nul');clc;
        end
    end
end


disp('== Merci! ==========================================================')
%%=========================================================================