%==========================================================================
% Programme d'apprentissage des méthodes de classification e
%
% F. Nougarou, octobre 2017
%==========================================================================
clear all; close all; clc;
%==========================================================================
Mat_mov = [1:5];

Nb_mov = length(Mat_mov);
Nb_mov_u = Nb_mov;
val_version = '29_novembre';
% 'target_charact','RMS','MAV','ZF','WL','SSC','LOG_VAR','VARi','SQRT_VAR','mov','Nb_ample_Moder','Nb_ample_Max','nb_block');
%==========================================================================
%-- Extraction des données APP
Target_U = [];
Input_U_class = [];
for bras = 1:2
    for mov = Mat_mov
        name_load_file = [];
        name_load_file = sprintf('APP_B%d_M%d_%s.mat',bras,mov,val_version);
        load(name_load_file)
        nb_block=2;
        block=1;
        Target_U = [Target_U;
            (mov-1)*ones(length(target_charact(1:(length(target_charact)/nb_block),:)),1)];

        Input_U_class = [Input_U_class;
            RMS(1:(length(RMS)/nb_block),:) MAV(1:(length(MAV)/nb_block),:) WL(1:(length(WL)/nb_block),:) ZF(1:(length(ZF)/nb_block),:) SSC(1:(length(SSC)/nb_block),:)]; 

    end
end

%-- Normalisation des Matrices --------------------------------------------
pos_Mov_ref = []; pos_Mov_ref = find(Target_U(1:(length(RMS)/nb_block),:) == 0);

val_norm_std = []; val_norm_std = std(Input_U_class(pos_Mov_ref,:));
val_norm_mean = []; val_norm_mean = mean(Input_U_class(pos_Mov_ref,:));


Input_U_class_norm = [];
for nn = 1:(size(Input_U_class,1)/(2*nb_block))
    Input_U_class_norm(nn,:) = (Input_U_class(nn,:)-val_norm_mean(1,:))./val_norm_std;
end


%== Classification ========================================================
Ind_EnsClass = Mat_mov;
Nb_mov = length(Ind_EnsClass);
pos_class = []; pos_class = find(Target_U(1:(length(Target_U)/(2*nb_block)),:) <=Ind_EnsClass(end)-1);

Input_U_class_norm_U = Input_U_class_norm(pos_class,:);
Target_UU = Target_U(pos_class,1);
 
W_LDA_Myo = [];
W_LDA_Myo = LDA(Input_U_class_norm_U,Target_UU);


Input_U_class_norm = [];
for nn = 1:(size(Input_U_class,1))
    Input_U_class_norm(nn,:) = (Input_U_class(nn,:)-val_norm_mean(1,:))./val_norm_std;
end
pos_class = []; pos_class = find(Target_U <=Ind_EnsClass(end)-1);

Input_U_class_norm_U = Input_U_class_norm(pos_class,:);
Target_UU = Target_U(pos_class,1);

fact_div =  nanmax(nanmax(abs(W_LDA_Myo)));
L_A_Myo = [];
L_A_Myo = ([ones((length(Target_UU)),1) Input_U_class_norm_U] * W_LDA_Myo');
P_A_Myo = [];
P_A_Myo = exp(L_A_Myo) ./ repmat(sum(exp(L_A_Myo),2),[1 Nb_mov]);
pos_max = []; val_max = [];
[val_max pos_max] = max(P_A_Myo(:,:)');
outputs_A_dec = [];
outputs_A_dec = (pos_max-1)';

pour_A_Myo = zeros(Nb_mov,1);
for m = 1:Nb_mov
    posi = []; posi = find(Target_UU == m-1);
    pos_val = []; pos_val = find(outputs_A_dec(posi) == m-1);
    pour_A_Myo(m,1) = length(pos_val)/length(posi)*100;
end
pour_A_Myo
mean(pour_A_Myo)

%-- Sauvegarde des données ------------------------------------------------
save data_APP_Coef_LDA.mat W_LDA_Myo val_norm_mean val_norm_std Nb_mov Mat_mov nb_myo;
%%=========================================================================

