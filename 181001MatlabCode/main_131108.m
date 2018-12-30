% wavelet + Higuchi's fractal dimension + sliding window
% Lu Yang 
% 2013.9.23
% 2013.9.25:
%   1. use linkage and plist function to calculate the distance and draw
%   tree
%   2. add comments
% 2013.9.26:
%   1. test other distance measure
%   2. enlarge the number of protein sequences
% 2013.10.27:
%   1. add 20 ND5 protein sequence
%   2. make all the fd sequences same length
% 2013.11.8:
%   1. add adcc features
clear;clc;

% for level_i=1:5
% read data and set parameters
% raw_data = textread('ND5.txt','%s');
%raw_data = textread('ND5_20.txt','%s');
raw_data = textread('ND5_new.txt','%s');
n = length(raw_data)/2;
% use_level = 2;
for use_level = 2:2
% window_width = 11;
for window_width = 17:17


% separate name and sequence from the file
for i=1:n
    protein_name(i) = raw_data((i-1)*2+1);
    protein_sequence(i) = raw_data(i*2);
end
number_map = [1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,13,14,15,16,17,0,18,19,0,20,0];
h_map = [1.8 2.5 -3.5 -3.5 2.8 -0.4 -3.2 4.5 -3.9 3.8 1.9 -3.5 -1.6 -3.5 -4.5 -0.8 -0.7 4.2 -0.9 -1.3];
p_map_1 = [2.34 1.71 2.09 2.19 1.83 2.34 1.82 2.36 2.18 2.36 2.28 2.02 1.99 2.17 2.17 2.21 2.63 2.32 2.38 2.20];
p_map_2 = [9.69 10.78 9.82 9.67 9.13 9.6 9.17 9.68 8.95 9.6 9.21 8.8 10.6 9.13 9.04 9.15 10.43 9.62 9.39 9.11];
p_map_3 = [6 5.07 2.77 3.22 5.48 5.97 7.59 6.02 9.74 5.98 5.74 5.41 6.3 5.65 10.76 5.68 5.6 5.96 5.89 5.66];
w_map = [71.08 103.15 115.09 129.12 147.18 57.05 137.14 113.16 128.18 113.16 131.2 114.11 97.12 128.13 156.19 87.08 101.11 99.13 186.22 163.18];

% mix_map = mapstd(h_map)+mapstd(p_map_1)+mapstd(p_map_2)+...
%     mapstd(p_map_3)+mapstd(w_map);

h_code = [1,3,0,0,1,3,0,1,0,1,1,0,3,0,0,2,2,1,2,2];
h_map_s = mapstd(h_map);
p_map_1s = mapstd(p_map_1);
p_map_2s = mapstd(p_map_2);
p_map_3s = mapstd(p_map_3);
w_map_s = mapstd(w_map);

for j=1:n
    clc;
    disp(use_level);
    disp(window_width);
    disp(j);
    % map each residue to the number and then corresponding feature
    protein_number = number_map(protein_sequence{j}-'A'+1);
    w_feature = w_map(protein_number);
    h_feature = h_map(protein_number);
%     h_feature = h_code(protein_number);
    p_feature = [p_map_1(protein_number),...
        p_map_2(protein_number)];
%     p_feature = p_map_1(protein_number) + ...
%             p_map_2(protein_number);
    
    % set used feature as the initial signal
%     ac = [h_feature,p_feature];
    ac = p_feature;

    for i=1:use_level
        % select a mother wavelet and use dwt to calculate ac and dc 
        [ac,dc] = dwt(ac,'haar');
        ac_length = length(ac);
        % slide the window, calculate the fractal dimension for each window
        for k=1:ac_length-window_width+1
            fd(k) = hfd(ac(k:k+window_width-1),min(floor(window_width/2),5));
        end
%         if i==use_level && j <= 10
%             subplot(4,n/2,j);plot(ac);
%             subplot(4,n/2,j+10);plot(fd);
%         end
%         axis tight
%         if i==use_level && j > 10
%             subplot(4,n/2,j+10);plot(ac);
%             subplot(4,n/2,j+20);plot(fd);
%         end
%         axis tight
        % record current fractal dimension and ac
        fd_record{j,i} = fd;
        ac_record{j,i} = ac;
        fd_dim_record(j,i) = length(fd);
        ac_dim_record(j,i) = length(ac);
        clear fd
    end
end

% the 8th protein is longer than others, delete the last bit of it
for i=1:use_level
    fd_min_length = min(fd_dim_record(:,i));
    ac_min_length = min(ac_dim_record(:,i));
    for j=1:n
        if fd_dim_record(j,i)>fd_min_length
            fd_record{j,i}(fd_min_length+1:end) = [];
        end
        if ac_dim_record(j,i)>ac_min_length
            ac_record{j,i}(ac_min_length+1:end) = [];
        end
    end
%     ac_record{8,i}(end) = [];
end

% connect all the fractal dimension of a specified level to a matrix 
fd2 = cat(1,fd_record{:,use_level})';
% figure;
% calculate the distance and draw tree
% h=figure('visible','off');
figure;
dendrogram(linkage(pdist(fd2','cosine')*100,'single'),'orientation','left','label',protein_name);
% file_name = strcat('result_131108/nonstd/p1&p2_d_',num2str(use_level),'_',num2str(window_width),'.jpg');
% saveas(h,file_name,'jpg');
end
end
load chirp
sound(y,Fs)