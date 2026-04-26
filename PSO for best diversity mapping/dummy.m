format compact
clear
clc
clf reset

% ----------

load("maps_total.mat")

field_names = string(fieldnames(maps_total));
for n=1:numel(field_names)
    maps_total.(field_names(n)).factor;
end

input_vector = rand(1,9);

% percent_dist = [
% 0.9 %TiO2
% 0.8 %MgO
% 0.5 %FeO
% 0.3 %CaO
% 0.3 %AlO3
% 0.7 %SiO2
% 0.4 %H2
% 0.5 %Th
% 0.9 %gradient 
% ];
percent_dist = input_vector.';

factor_list = [
"TiO2"
"MgO"
"FeO"
"CaO"
"Al2O3"
"SiO2"
"H2"
"Th"
"gradient"
];

standard_size = [500,1000];

percent_map = single(zeros(standard_size));

for n=1:numel(field_names)

    map_tmp = maps_total.(field_names(n)).map;
    map_factor_tmp = maps_total.(field_names(n)).factor;
    map_bounds = [min(map_tmp(:)),max(map_tmp(:))];

    map_norm = map_tmp./max(map_tmp(:));

    ind_f = find(map_factor_tmp==factor_list);
    percent_limit = percent_dist(ind_f);

    percent_region = single(map_norm > percent_limit);
    percent_map = percent_map + percent_region;
end

full_criteria_map = percent_map == single(numel(field_names)-1);

regions = bwconncomp( imgaussfilt(single(full_criteria_map),0.5)>0 ).NumObjects;

sum(full_criteria_map(:))

hold on
grid on
axis tight equal
imagesc(flipud( imgaussfilt(single(full_criteria_map),0.5)>0 ))