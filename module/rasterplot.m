function [] = rasterplot(res,param)
% A function to show scatterplots and save them in folder output

ne = param.ne;
ni = param.ni;
spike = res.spike;
num_E = sum(spike(1,1:ne));
num_I = sum(spike(1,ne+1:ne+ni));
coor_E = zeros(num_E,2);
coor_I = zeros(num_I,2);
index_E = 1;
index_I = 1;
for i=1:ne
num_Ei = spike(1,i);
coor_E(index_E:index_E+num_Ei-1,1) = i;
coor_E(index_E:index_E+num_Ei-1,2) = spike(2:1+num_Ei,i)*1000;
index_E = index_E + num_Ei;
end

for i=(ne+1):(ne+ni)
num_Ii = spike(1,i);
coor_I(index_I:index_I+num_Ii-1,1) = i;
coor_I(index_I:index_I+num_Ii-1,2) = spike(2:1+num_Ii,i)*1000;
index_I = index_I + num_Ii;
end

scatter(coor_E(:,2), coor_E(:,1),10,'.','r');
hold on
scatter(coor_I(:,2), coor_I(:,1),10,'.','b');
ylabel('Index');
set(gca,'fontsize',11);
end


