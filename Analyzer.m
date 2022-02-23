clc;
%clear;
close all;
dbstop if error
% Reading Data

% % Way 1 : First I imported data from the excell file manually by matlab excel reader named 'Fe10Ni10Si5Gparun2m' variable, then :
% data = table2array(Fe10Ni10Si5Gparun2m(:,1));
% data = transpose(data);

% Way 2 :
file_name = 'Example file MatLab (june 6).xlsx'; %'Fe10Ni10Si_5Gpa_run2m.xlsx';
data = importdata(file_name);
data = transpose(data.data);

if strcmp(file_name,'Example file MatLab (june 6).xlsx') || strcmp(file_name,'Example file MatLab.xlsx')
    data(690:758)=[];
end

a = data;
b = [data(2:end) 0];
c = [data(3:end) 0 0];
d = [data(4:end) 0 0 0];
ab = a.*b;
ac = a.*c;
ad = a.*d;
negative_index = ab<0 & ac<0 & ad<0;
figure(1);plot(data,'.b');
title('Step1 : Temporary result');
xlabel('Index')
ylabel('Data')
grid on;grid minor;

negative_index_number = [find(negative_index)];
aa = a(negative_index_number);

neg_diff = [];
k = 1;
interval_neg2=[];
for i=1:(length(negative_index_number)-1)
        interval_neg = [negative_index_number(i)+1:negative_index_number(i+1)];
        if sum(data(interval_neg)<0)>2 && mode(data(interval_neg))<-1e-4
            temp_data_neg = data(interval_neg);
            temp_data_neg = temp_data_neg(temp_data_neg<0);
            interval_neg = interval_neg(temp_data_neg<0);
            figure(1);
            hold on;
            plot(interval_neg,temp_data_neg,'.r');
            
            temp_data_before = data(interval_neg(1)-10:interval_neg(1)-1);
            temp_data_after = data(interval_neg(end)+1:interval_neg(end)+10);
            if k==44
                1;
            end
            n=1.4;
            n2 = max(ceil(n*length(interval_neg)),10);
            if abs(abs(abs(median(temp_data_before))-abs(median(temp_data_neg)))- abs(abs(median(temp_data_after))-abs(median(temp_data_neg))))<1.3e-4
                interval_neg2 = [interval_neg2 interval_neg(1) interval_neg(end)];
            elseif (abs(abs(median(temp_data_before))-abs(median(temp_data_neg))) < abs(abs(median(temp_data_after))-abs(median(temp_data_neg))))
                interval_pos = [interval_neg(1)-n2:interval_neg(1)-1];
                temp_data_pos = data(interval_pos);
                
                V_ind = [min([interval_pos , interval_neg]):max([interval_pos , interval_neg])];
                VP_ind = V_ind(data(V_ind)>=0);
                isempty(VP_ind)
                VP = data(VP_ind);
                VN_ind = V_ind(data(V_ind)<0);
                VN = data(VN_ind);
                dif_VN = diff(VN);
                T_a_ind = [VN_ind(end)+1:VN_ind(end)+10];
                T_a = data(T_a_ind);
                [~,ind1] = max(abs(diff(VP(1:end-3))));
                if ~isempty(ind1)
                    ind = VP_ind(1)+ind1;
                    VP_ind = [ind:VP_ind(end)];
                    VP = data(VP_ind);
                    T_b_ind = ind-10:ind-1;
                    T_b = data(T_b_ind);
                    V_ind = [VP_ind(1):VN_ind(end)];
                    final.VP{k} = [VP;VP_ind];
                    final.VN{k} = [VN;VN_ind];
                    final.Tb{k} = [T_b;T_b_ind];
                    final.Ta{k} = [T_a;T_a_ind];
                    k = k+1;
                    plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                    pause(1e-5);
                end
            else
                interval_pos = [interval_neg(end)+1:interval_neg(end)+n2];
                temp_data_pos = data(interval_pos);
                
                V_ind = [min([interval_pos , interval_neg]):max([interval_pos , interval_neg])];
                VP_ind = V_ind(data(V_ind)>=0);
                isempty(VP_ind)
                VP = data(VP_ind);
                VN_ind = V_ind(data(V_ind)<0);
                VN = data(VN_ind);
                dif_VN = diff(VN);
                if ~isempty(VN_ind)
                    T_b_ind = [VN_ind(1)-10:VN_ind(1)-1];
                    T_b = data(T_b_ind);
                    [~,ind1] = max(abs(diff(VP(3:end))));
                    if ~isempty(ind1)
                        ind = VP_ind(1)+ind1+2;
                        VP_ind = [VP_ind(1):ind-1];
                        VP = data(VP_ind);
                        T_a_ind = ind:ind+9;
                        T_a = data(T_a_ind);
                        V_ind = [VN_ind(1):VP_ind(end)];
                        final.VP{k} = [VP;VP_ind];
                        final.VN{k} = [VN;VN_ind];
                        final.Tb{k} = [T_b;T_b_ind];
                        final.Ta{k} = [T_a;T_a_ind];
                        k = k+1;
                        plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                        pause(1e-5);
                    end
                end
            end            
        end
end
for jj=1:length(interval_neg2)-1
   interval_neg = interval_neg2(jj):interval_neg2(jj+1);
   n2 = max(ceil(n*length(interval_neg)),10);
   if length(interval_neg)>1 && mode(data(interval_neg))<-1e-4
        temp_data_neg = data(interval_neg);
        temp_data_neg = temp_data_neg(temp_data_neg<0);
        interval_neg = interval_neg(temp_data_neg<0);
        figure(1);
        hold on;
        plot(interval_neg,temp_data_neg,'.r')

        temp_data_before = data(interval_neg(1)-n2:interval_neg(1)-1);
        temp_data_after = data(interval_neg(end)+1:interval_neg(end)+n2);
        if max(abs(diff(temp_data_before)))>max(abs(diff(temp_data_after)))
            interval_pos = [interval_neg(1)-n2:interval_neg(1)-1];
            temp_data_pos = data(interval_pos);
            
            V_ind = [min([interval_pos , interval_neg]):max([interval_pos , interval_neg])];
            VP_ind = V_ind(data(V_ind)>=0);
            isempty(VP_ind)
            VP = data(VP_ind);
            VN_ind = V_ind(data(V_ind)<0);
            VN = data(VN_ind);
            dif_VN = diff(VN);
            T_a_ind = [VN_ind(end)+1:VN_ind(end)+10];
            T_a = data(T_a_ind);
            [~,ind1] = max(abs(diff(VP(1:end-3))));
            if ~isempty(ind1)
                ind = VP_ind(1)+ind1;
                VP_ind = [ind:VP_ind(end)];
                VP = data(VP_ind);
                T_b_ind = ind-10:ind-1;
                T_b = data(T_b_ind);
                V_ind = [VP_ind(1):VN_ind(end)];
                final.VP{k} = [VP;VP_ind];
                final.VN{k} = [VN;VN_ind];
                final.Tb{k} = [T_b;T_b_ind];
                final.Ta{k} = [T_a;T_a_ind];
                k = k+1;
                plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                pause(1e-5);
            end
        else
            interval_pos = [interval_neg(end)+1:interval_neg(end)+n2];
            temp_data_pos = data(interval_pos);
            
            V_ind = [min([interval_pos , interval_neg]):max([interval_pos , interval_neg])];
            VP_ind = V_ind(data(V_ind)>=0);
            isempty(VP_ind)
            VP = data(VP_ind);
            VN_ind = V_ind(data(V_ind)<0);
            VN = data(VN_ind);
            dif_VN = diff(VN);
            if ~isempty(VN_ind)
                T_b_ind = [VN_ind(1)-10:VN_ind(1)-1];
                T_b = data(T_b_ind);
                [~,ind1] = max(abs(diff(VP(3:end))));
                if ~isempty(ind1)
                    ind = VP_ind(1)+ind1+2;
                    VP_ind = [VP_ind(1):ind-1];
                    VP = data(VP_ind);
                    T_a_ind = ind:ind+9;
                    T_a = data(T_a_ind);
                    V_ind = [VN_ind(1):VP_ind(end)];
                    final.VP{k} = [VP;VP_ind];
                    final.VN{k} = [VN;VN_ind];
                    final.Tb{k} = [T_b;T_b_ind];
                    final.Ta{k} = [T_a;T_a_ind];
                    k = k+1;
                    plot(V_ind,data(V_ind),'*k',T_b_ind,T_b,'.r',T_a_ind,T_a,'.g');
                    pause(1e-5);
                end
            end
        end
   end
end

figure(2);
plot(data,'.b'); hold on;
title('Step2 : Final Result')
xlabel('Index')
ylabel('Data')
grid on;grid minor;

VP_dif_max = [];
VN_dif_max = [];
Ta_dif_max = [];
Tb_dif_max = [];
for kk = 1:length(final.VP)
     if kk==1
         Ta_ind = final.Ta{kk}(2,:);
         VP_b_ind = final.VP{kk+1}(2,:);
         VN_b_ind = final.VN{kk+1}(2,:);
         for i=1:length(Ta_ind)
             item = Ta_ind(i);
             if any(VP_b_ind==item) || any(VN_b_ind==item)
                 Ta_ind(i) = -1;
             end
         end
         Ta = [final.Ta{kk}(1,Ta_ind~=-1);Ta_ind(Ta_ind~=-1)];
         final.Ta{kk} = Ta;
     elseif kk==length(final.VP)
         Tb_ind = final.Tb{kk}(2,:);
         VP_a_ind = final.VP{kk-1}(2,:);
         VN_a_ind = final.VN{kk-1}(2,:);
         for i=1:length(Ta_ind)
             item = Tb_ind(i);
             if any(VP_a_ind==item) || any(VN_a_ind==item)
                 Tb_ind(i) = -1;
             end
         end
         Tb = [final.Tb{kk}(1,Tb_ind~=-1);Tb_ind(Tb_ind~=-1)];
         final.Tb{kk} = Tb;
     else
        Ta_ind = final.Ta{kk}(2,:);
        Tb_ind = final.Tb{kk}(2,:);
        VP_a_ind = final.VP{kk-1}(2,:);
        VP_b_ind = final.VP{kk+1}(2,:);
        VN_a_ind = final.VN{kk-1}(2,:);
        VN_b_ind = final.VN{kk+1}(2,:);
        for i=1:length(Ta_ind)
           item = Ta_ind(i);
           if any(VP_b_ind==item) || any(VP_a_ind==item) || any(VN_b_ind==item) || any(VN_a_ind==item)
               Ta_ind(i) = -1;
           end
           item = Tb_ind(i);
           if any(VP_b_ind==item) || any(VP_a_ind==item) || any(VN_b_ind==item) || any(VN_a_ind==item)
               Tb_ind(i) = -1;
           end
        end
        Ta = [final.Ta{kk}(1,Ta_ind~=-1);Ta_ind(Ta_ind~=-1)];
        final.Ta{kk} = Ta;
        Tb = [final.Tb{kk}(1,Tb_ind~=-1);Tb_ind(Tb_ind~=-1)];
        final.Tb{kk} = Tb;
     end
     VP = final.VP{kk}(1,:);
     VN = final.VN{kk}(1,:);
     Ta = final.Ta{kk}(1,:);
     Tb = final.Tb{kk}(1,:);
     VP_ind = final.VP{kk}(2,:);
     VN_ind = final.VN{kk}(2,:);
     Ta_ind = final.Ta{kk}(2,:);
     Tb_ind = final.Tb{kk}(2,:);

     n3 = 10;
     n4 = 4;
     th_min = 1.85e-4;1.39e-4; % default = 0;
     th_max = 1000;%0.04; %default = 1000;
     VP_dif_max = [VP_dif_max max(abs(diff(VP)))];
     VN_dif_max = [VN_dif_max max(abs(diff(VN)))];
     Ta_dif_max = [Ta_dif_max max(abs(diff(VP)))];
     Tb_dif_max = [Tb_dif_max max(abs(diff(VP)))];

%      VP_dif = mean(abs(diff(VP)));
     VP_dif = sort(abs(diff(VP)));
     if length(VP_dif)>n4
         VP_dif = VP_dif(end-n4);
         VP_dif1 = [abs(diff(VP)) max(abs(diff(VP)))];
         VP_dif2 = [max(abs(diff(VP))) abs(diff(VP))];
         VP_del_ind = (VP_dif1>n3*VP_dif & VP_dif2>n3*VP_dif) | (abs(VP)<th_min) | (abs(VP)>th_max);
         VP = VP(~VP_del_ind);
         VP_ind = VP_ind(~VP_del_ind);
     end
%      VN_dif = mean(abs(diff(VN)));
     VN_dif = sort(abs(diff(VN)));
     if length(VN_dif)>n4
         VN_dif = VN_dif(end-n4);
         VN_dif1 = [abs(diff(VN)) max(abs(diff(VN)))];
         VN_dif2 = [max(abs(diff(VN))) abs(diff(VN))];
         VN_del_ind = (VN_dif1>n3*VN_dif & VN_dif2>n3*VN_dif) | (abs(VN)<th_min) | (abs(VN)>th_max);
         VN = VN(~VN_del_ind);
         VN_ind = VN_ind(~VN_del_ind);
     end
     
%      Ta_dif = mean(abs(diff(Ta)));
     Ta_dif = sort(abs(diff(Ta)));
     if length(Ta_dif)>n4
         Ta_dif = Ta_dif(end-n4);
         Ta_dif1 = [abs(diff(Ta)) max(abs(diff(Ta)))];
         Ta_dif2 = [max(abs(diff(Ta))) abs(diff(Ta))];
         Ta_del_ind = Ta_dif1>n3*Ta_dif & Ta_dif2>n3*Ta_dif;
         Ta = Ta(~Ta_del_ind);
         Ta_ind = Ta_ind(~Ta_del_ind);
     end
%      Tb_dif = mean(abs(diff(Tb)));
     Tb_dif = sort(abs(diff(Tb)));
     if length(Tb_dif)>n4
         Tb_dif = Tb_dif(end-n4);
         Tb_dif1 = [abs(diff(Tb)) max(abs(diff(Tb)))];
         Tb_dif2 = [max(abs(diff(Tb))) abs(diff(Tb))];
         Tb_del_ind = (Tb_dif1>n3*Tb_dif & Tb_dif2>n3*Tb_dif);
         Tb = Tb(~Tb_del_ind);
         Tb_ind = Tb_ind(~Tb_del_ind);
     end
     
     final.VP{kk} = [VP;VP_ind];
     final.VN{kk} = [VN;VN_ind];
     final.Ta{kk} = [Ta;Ta_ind];
     final.Tb{kk} = [Tb;Tb_ind];
     final_avg.VP(kk) = mean(VP);
     final_avg.VN(kk) = mean(VN);
     final_avg.Ta(kk) = mean(Ta);
     final_avg.Tb(kk) = mean(Tb);
     plot(final.VP{kk}(2,:),final.VP{kk}(1,:),'*k',...
         final.VN{kk}(2,:),final.VN{kk}(1,:),'*r',...
         final.Ta{kk}(2,:),final.Ta{kk}(1,:),'.r',...
         final.Tb{kk}(2,:),final.Tb{kk}(1,:),'.g');
     pause(1e-5);
end