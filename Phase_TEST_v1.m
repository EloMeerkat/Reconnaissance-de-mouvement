
%-- Phase TEST ------------------------------------------------------------

%-- Boucle appel et construction ------------------------------------------
Target_all_bras1 = [];
Target_all_bras2 = [];
Features_all_bras1 = [];
Features_all_bras2 = [];
Nb_bloc = 20;
Nb_mov = 5;

for b = 1:2
    for m = 1:Nb_bloc
        Features_Bras_bloc = [];
        Target_Bras_bloc = [];
        for n = 1:Nb_mov
            if m < 10
                name_file = sprintf('Bloc_0%d_Mouv%d_bras%d.mat',m,n,b);
            else
                name_file = sprintf('Bloc_%d_Mouv%d_bras%d.mat',m,n,b);
            end
            load(name_file,'Target_Bras','Features_Bras');
            Features_Bras_bloc = [Features_Bras_bloc Features_Bras'];
            Target_Bras_bloc = [Target_Bras_bloc Target_Bras'];                      
        end
        if b == 1
            Features_all_bras1(:,:,m) = Features_Bras_bloc;
            Target_all_bras1(:,m) = Target_Bras_bloc';
        else
             Features_all_bras2(:,:,m) = Features_Bras_bloc;
             Target_all_bras2(:,m) = Target_Bras_bloc';
        end
    end
end


Bloc_app = 1:10;

for n = Bloc_app
    %-- construction data APP Pour le bras 1 ------------------------------
    Features_all_bras1_APP = [];
    Target_all_bras1_APP = [];
    for m = 1:n
        Features_all_bras1_APP = [Features_all_bras1_APP Features_all_bras1(:,:,m)];
        Target_all_bras1_APP   = [Target_all_bras1_APP Target_all_bras1(:,m)'];
    end
    nb_col = size(Features_all_bras1_APP,2);
    new_order = randperm(nb_col);
    Features_all_bras1_APP = Features_all_bras1_APP(:,new_order);
    Target_all_bras1_APP   = Target_all_bras1_APP(1,new_order)-1;
    
    %-- construction data DET Pour le bras 1 ------------------------------
    Features_all_bras1_DET = [];
    Target_all_bras1_DET = [];
    for m = n+1:Nb_bloc
        Features_all_bras1_DET = [Features_all_bras1_DET Features_all_bras1(:,:,m)];
        Target_all_bras1_DET   = [Target_all_bras1_DET Target_all_bras1(:,m)'];
    end    
    nb_col = size(Features_all_bras1_DET,2);
    new_order = randperm(nb_col);
    Features_all_bras1_DET = Features_all_bras1_DET(:,new_order);
    Target_all_bras1_DET   = Target_all_bras1_DET(1,new_order);
    
    %-- LDA pour le bras 1 ------------------------------------------------
    B1_W_LDA_Myo = [];
    B1_W_LDA_Myo = LDA(Features_all_bras1_APP', Target_all_bras1_APP');
    
     %-- construction data APP Pour le bras 2 ------------------------------
    Features_all_bras2_APP = [];
    Target_all_bras2_APP = [];
    for m = 1:n
        Features_all_bras2_APP = [Features_all_bras2_APP Features_all_bras2(:,:,m)];
        Target_all_bras2_APP   = [Target_all_bras2_APP Target_all_bras2(:,m)'];
    end
    nb_col = size(Features_all_bras2_APP,2);
    new_order = randperm(nb_col);
    Features_all_bras2_APP = Features_all_bras2_APP(:,new_order);
    Target_all_bras2_APP   = Target_all_bras2_APP(1,new_order)-1;
    
    %-- construction data DET Pour le bras 2 ------------------------------
    Features_all_bras2_DET = [];
    Target_all_bras2_DET = [];
    for m = n+1:Nb_bloc
        Features_all_bras2_DET = [Features_all_bras2_DET Features_all_bras2(:,:,m)];
        Target_all_bras2_DET   = [Target_all_bras2_DET Target_all_bras2(:,m)'];
    end    
    nb_col = size(Features_all_bras2_DET,2);
    new_order = randperm(nb_col);
    Features_all_bras2_DET = Features_all_bras2_DET(:,new_order);
    Target_all_bras2_DET   = Target_all_bras2_DET(1,new_order);
    
    %-- LDA pour le bras 2 ------------------------------------------------
    B2_W_LDA_Myo = [];
    B2_W_LDA_Myo = LDA(Features_all_bras2_APP', Target_all_bras2_APP');
    
    B11_L_A_Myo = [];
    B11_L_A_Myo = ([ones((length(Target_all_bras1_DET')),1) Features_all_bras1_DET'] * B1_W_LDA_Myo');
    B11_P_A_Myo = [];
    B11_P_A_Myo = exp(B11_L_A_Myo) ./ repmat(sum(exp(B11_L_A_Myo),2),[1 Nb_mov]);
    pos_max = []; val_max = [];
    [val_max pos_max] = max(B11_P_A_Myo(:,:)');
    B11_outputs_A_dec = [];
    B11_outputs_A_dec = (pos_max-1)';
    
    B12_L_A_Myo = [];
    B12_L_A_Myo = ([ones((length(Target_all_bras2_DET')),1) Features_all_bras2_DET'] * B1_W_LDA_Myo');
    B12_P_A_Myo = [];
    B12_P_A_Myo = exp(B12_L_A_Myo) ./ repmat(sum(exp(B12_L_A_Myo),2),[1 Nb_mov]);
    pos_max = []; val_max = [];
    [val_max pos_max] = max(B12_P_A_Myo(:,:)');
    B12_outputs_A_dec = [];
    B12_outputs_A_dec = (pos_max-1)';
    
    B21_L_A_Myo = [];
    B21_L_A_Myo = ([ones((length(Target_all_bras1_DET')),1) Features_all_bras1_DET'] * B2_W_LDA_Myo');
    B21_P_A_Myo = [];
    B21_P_A_Myo = exp(B21_L_A_Myo) ./ repmat(sum(exp(B21_L_A_Myo),2),[1 Nb_mov]);
    pos_max = []; val_max = [];
    [val_max pos_max] = max(B21_P_A_Myo(:,:)');
    B21_outputs_A_dec = [];
    B21_outputs_A_dec = (pos_max-1)';
    
    B22_L_A_Myo = [];
    B22_L_A_Myo = ([ones((length(Target_all_bras2_DET')),1) Features_all_bras2_DET'] * B2_W_LDA_Myo');
    B22_P_A_Myo = [];
    B22_P_A_Myo = exp(B22_L_A_Myo) ./ repmat(sum(exp(B22_L_A_Myo),2),[1 Nb_mov]);
    pos_max = []; val_max = [];
    [val_max pos_max] = max(B22_P_A_Myo(:,:)');
    B22_outputs_A_dec = [];
    B22_outputs_A_dec = (pos_max-1)';

    B11_pour_A_Myo = zeros(Nb_mov,1);
    B12_pour_A_Myo = zeros(Nb_mov,1);
    B21_pour_A_Myo = zeros(Nb_mov,1);
    B22_pour_A_Myo = zeros(Nb_mov,1);
    for m = 1:Nb_mov
        B11_posi = 0; B11_posi = find(Target_all_bras1_DET == m-1);
        B11_pos_val = []; B11_pos_val = find(B11_outputs_A_dec(B11_posi) == m-1);
        B11_pour_A_Myo(m,1) =( 1-length(B11_pos_val)/length(B11_posi))*100;
        
        B12_posi = 0; B12_posi = find(Target_all_bras2_DET == m-1);
        B12_pos_val = []; B12_pos_val = find(B12_outputs_A_dec(B12_posi) == m-1);
        B12_pour_A_Myo(m,1) =( 1-length(B12_pos_val)/length(B12_posi))*100;
        
        B21_posi = 0; B21_posi = find(Target_all_bras1_DET == m-1);
        B21_pos_val = []; B21_pos_val = find(B21_outputs_A_dec(B21_posi) == m-1);
        B21_pour_A_Myo(m,1) =( 1-length(B21_pos_val)/length(B21_posi))*100;
        
        B22_posi = 0; B22_posi = find(Target_all_bras2_DET == m-1);
        B22_pos_val = []; B22_pos_val = find(B22_outputs_A_dec(B22_posi) == m-1);
        B22_pour_A_Myo(m,1) =( 1-length(B22_pos_val)/length(B22_posi))*100;
    end
    
    B11_mean_pour_A_Myo(n,1)= mean(B11_pour_A_Myo);
    B12_mean_pour_A_Myo(n,1)= mean(B12_pour_A_Myo);
    B21_mean_pour_A_Myo(n,1)= mean(B21_pour_A_Myo);
    B22_mean_pour_A_Myo(n,1)= mean(B22_pour_A_Myo);
    
end

subplot(2,2,1)
plot(B11_mean_pour_A_Myo)
axis([1 10 0 100])
title('Apprentissage du bras droit et test pour le bras droit')
ylabel('Pourcentage d erreur')
xlabel('Nombre de blocs utilisé pour l apprentissage /20')
subplot(2,2,2)
plot(B12_mean_pour_A_Myo)
axis([1 10 0 100])
title('Apprentissage du bras droit et test pour le bras gauche')
ylabel('Pourcentage d erreur')
xlabel('Nombre de blocs utilisé pour l apprentissage /20')
subplot(2,2,3)
plot(B22_mean_pour_A_Myo)
axis([1 10 0 100])
title('Apprentissage du bras gauche et test pour le bras gauche')
ylabel('Pourcentage d erreur')
xlabel('Nombre de blocs utilisé pour l apprentissage /20')
subplot(2,2,4)
plot(B21_mean_pour_A_Myo)
axis([1 10 0 100])
title('Apprentissage du bras gauche et test pour le bras droit')
ylabel('Pourcentage d erreur')
xlabel('Nombre de blocs utilisé pour l apprentissage /20')

