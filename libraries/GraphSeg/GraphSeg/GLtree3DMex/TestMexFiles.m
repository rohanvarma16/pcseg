%Script to compile and test GLTree files


clc
clear all
close all


%% Compile
fprintf('COMPILING:\n')
 
mex BuildGLTree.cpp
fprintf('\tBuildGLTree : mex succesfully completed.\n') 

mex KNNSearch.cpp
fprintf('\tKNNSearch : mex succesfully completed.\n') 

mex DeleteGLTree.cpp
fprintf('\tDeleteGLTree : mex succesfully completed.\n\n') 



%% Runs a few test

N=10000;%number of reference points
p=rand(N,3);


Nq=10;%number of query points
qp=rand(Nq,3);

k=4;%k neighbours for k Search

fprintf('RANDOM POINTS GENERATED\n\n')


fprintf('BUILDING THE DATA STRUCTURE:\n')
tic
ptrtree=BuildGLTree(p);
fprintf('\tGLTree built in %4.4f s\n\treturned pointer %4.0f:\n\n',toc,ptrtree);


% 
% 
fprintf('START K NEAREST NEIGHBOR SEARCH:\n')
tic
[KNNG,distances]=KNNSearch(p,p,ptrtree,k);
fprintf('\t  NNsearch of%4.0f neighbors in %4.0f reference points and %4.0f query points took: %4.4f s\n\n',k,N,Nq,toc);


fprintf('DELETING THE TREE\n\n')
DeleteGLTree(ptrtree);
 

fprintf('TEST SUCCESFULLY COMPLETED !!!\n\n')

%plot the NNG


figure(1)
title([num2str(k),' Neighbours'],'fontsize',14);
axis equal
hold on
plot3(p(:,1),p(:,2),p(:,3),'g.')
p1x=p(KNNG(:,k),1);
p1y=p(KNNG(:,k),2);
p1z=p(KNNG(:,k),3);
    for j=1:k-1%just the firs point(itself);
        
    p2x=p(KNNG(:,j),1);
    p2y=p(KNNG(:,j),2);
    p2z=p(KNNG(:,j),3);
    
    
    plot3([p1x,p2x]',[p1y,p2y]',[p1z,p2z]','r-')
    end
