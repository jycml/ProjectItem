% wavelet + Higuchi's fractal dimension + sliding window
%  Yang 
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
clear;clc;

% for level_i=1:5
% read data and set parameters
raw_data = textread('ND5.txt','%s');
%raw_data = textread('ND5_20.txt','%s');
n = length(raw_data)/2;
max_level = 5;
use_level = 3;
window_width = 11;
%for window_width = 5:30


% separate name and sequence from the file
for i=1:n
    protein_name(i) = raw_data((i-1)*2+1);
    protein_sequence(i) = raw_data(i*2);
end
number_map = [1,0,2,3,4,5,6,7,8,0,9,10,11,12,0,13,14,15,16,17,0,18,19,0,20,0];
h_map = [1.8 2.5 -3.5 -3.5 2.8 -0.4 -3.2 4.5 -3.9 3.8 1.9 -3.5 -1.6 -3.5 -4.5 -0.8 -0.7 4.2 -0.9 -1.3];
%h_code = [1,3,0,0,1,3,0,1,0,1,1,0,3,0,0,2,2,1,2,2];

for j=1:n
    % map each residue to the number and then corresponding feature
    protein_number = number_map(protein_sequence{j}-'A'+1);
    h_feature = h_map(protein_number);
%     h_feature = h_code(protein_number);
    % set used feature as the initial signal
    ac = h_feature;

    for i=1:max_level
        % select a mother wavelet and use dwt to calculate ac and dc 
        [ac,dc] = dwt(ac,'haar');
%         subplot(max_level,n*2,(i-1)*n*2+j);plot(ac);
        ac_length = length(ac);
        % slide the window, calculate the fractal dimension for each window
        for k=1:ac_length-window_width+1
            fd(k) = hfd(ac(k:k+window_width-1),min(floor(window_width/2),5));
        end
%         subplot(max_level,n*2,(i-1)*n*2+j+n);plot(fd);
        % record current fractal dimension and ac
        fd_record{j,i} = fd;
        ac_record{j,i} = ac;
        fd_dim_record(j,i) = length(fd);
        ac_dim_record(j,i) = length(ac);
        clear fd
    end
end

% the 8th protein is longer than others, delete the last bit of it
for i=1:max_level
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
figure;
dendrogram(linkage(pdist(fd2','cosine'),'single'),'orientation','left','label',protein_name)
%end
