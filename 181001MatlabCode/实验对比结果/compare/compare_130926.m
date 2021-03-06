% compare different methods for protein sequences similarity
% Lu Yang 
% 2013.9.26:
%   1. compare our proposed method, Su' model, Zhang' model, Yao' model and
%   results generated by MEGA
clear; clc;

load result_130926/our
load result_130926/su
load result_130926/zhang
load result_130926/yao
load result_130926/mega
load result_130926/protein_name

figure;dendrogram(linkage(pdist(our','cosine'),'single'),'orientation','left','label',protein_name)
xlabel('Our model')
figure;dendrogram(linkage(pdist(su','cosine'),'single'),'orientation','left','label',protein_name)
xlabel('Su''s model')
figure;dendrogram(linkage(zhang,'single'),'orientation','left','label',protein_name)
xlabel('Zhang''s model')
figure;dendrogram(linkage(yao,'single'),'orientation','left','label',protein_name)
xlabel('Yao''s model')
figure;dendrogram(linkage(mega,'single'),'orientation','left','label',protein_name)
xlabel('MEGA')