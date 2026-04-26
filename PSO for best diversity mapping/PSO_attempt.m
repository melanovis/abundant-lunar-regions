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

%% PSO
n_DOF = 9;
var_size = [1, n_DOF];
var_min = 0;
var_max = 1;

max_iters = 100;
population = 8*20;
phi_1 = 2.05;
phi_2 = 2.05;
phi = phi_1+phi_2;
kappa = 1;
chi = 2*kappa/abs(2-phi-sqrt(phi^2-4*phi));
w = chi;
w_damp = 0.75;
w_original = 1;

c1 = chi*phi_1;
c2 = chi*phi_2;
% c1 = 1;
% c2 = 1;

max_velocity = (var_max-var_min)*0.3;
min_velocity = -max_velocity;

empty_particle.position = [];
empty_particle.velocity = [];
empty_particle.fitness = [];
empty_particle.best.position = [];
empty_particle.best.fitness = [];

particle_init_series = [];
particle_init_block = repmat(empty_particle, population, 1);
particle = repmat(empty_particle, population, 1);

checkgood_threshold = round(interp1([0,1],[20,40],rand()));
init_checkgood_quantity = 0;
checkgood_series = [];

if ~gcp().Connected %start up the cores
    delete(gcp('nocreate'));
    parpool('local',8);
end

parfor n=1:population
    particle(n).position = unifrnd(var_min,var_max,var_size);
    
    input_vector = particle(n).position;
    [fitness, ~, ~] = asses_criteria(input_vector,factor_list,maps_total);

    particle(n).fitness = fitness;
    particle(n).velocity = zeros(var_size); 
    particle(n).best.position = particle(n).position;
    particle(n).best.fitness = particle(n).fitness;
end

global_best.fitness = -inf;

for n=1:population
    if particle(n).best.fitness > global_best.fitness
        global_best.position = particle(n).best.position;
        global_best.fitness = particle(n).best.fitness;
    end
    particle(n).velocity = rand(var_size).*max_velocity;
end

fprintf("------------------\n")

iter = 1;
ag_countdown = 10;
steps_forward = 0;
global_best_map = single(zeros(standard_size));
global_best_dist = [];

while iter <= max_iters

    parfor n=1:population
        
        particle(n).velocity = w*particle(n).velocity ...
            + c1*rand(var_size).*(particle(n).best.position - particle(n).position) ...
            + c2*rand(var_size).*(global_best.position - particle(n).position);
        
        particle(n).velocity = max(particle(n).velocity, min_velocity);
        particle(n).velocity = min(particle(n).velocity, max_velocity);
    
        particle(n).position = particle(n).position + particle(n).velocity;
        
        particle(n).position = max(particle(n).position, var_min);
        particle(n).position = min(particle(n).position, var_max);
    
        input_vector = particle(n).position;
        [fitness, ~, ~] = asses_criteria(input_vector,factor_list,maps_total);

        particle(n).fitness = fitness;
    
        %update personal best
        if particle(n).fitness > particle(n).best.fitness
            particle(n).best.position = particle(n).position;
            particle(n).best.fitness = particle(n).fitness;
        end
       
    end

    update_plot = false;
    for n=1:population
        if particle(n).best.fitness > global_best.fitness
            global_best = particle(n).best;
            input_vector = particle(n).best.position;
            [~, output_vector, region_map] = asses_criteria(input_vector,factor_list,maps_total);
            global_best_map = region_map;
            global_best_dist = output_vector;
            steps_forward = steps_forward+1;
            update_plot = true;
        end
    end

    bestfitnesss(iter) = global_best.fitness;
    
    if iter > 20
        if std(bestfitnesss(end-3:end)) < 5e-3
            ag_countdown = ag_countdown-1;
        else
            ag_countdown = 10;
        end
        if ag_countdown <= 0
            fitness_score = [];
            for n=1:population
                fitness_score(n) = particle(n).fitness;
            end
            [~,ind_sort] = sort(fitness_score);
    
            if rand() < 0.5
                ind_ag = ind_sort(1:round(length(ind_sort)*0.5));
            else
                ind_ag = ind_sort;
            end

            for n=1:length(ind_ag)
                particle(ind_ag(n)).position = rand(1,n_DOF)*var_max;
                particle(ind_ag(n)).velocity = rand(1,n_DOF)*max_velocity;
            end
    
            w = w_original;
            ag_countdown = 10;
            fprintf("\n agitating.\n")
        end
    end

    if update_plot
        scatter(nan,nan)
        hold on
        grid on
        axis tight equal
        colormap(gray)
        imagesc(flipud(global_best_map));
        drawnow()
        hold off
    end
    
    w = w * w_damp;

    fprintf("\n completed %i, fitness: %3.2f, steps foward: %i.\n", iter, global_best.fitness,steps_forward)
    global_best_dist
    fprintf("\n")

    iter = iter+1;
end


function [fitness,output_vector, region_map] = asses_criteria(input_vector,factor_list,maps_total)

    percent_dist = interp1([0,1],[1e-3,0.99],input_vector.');

    percent_dist = round(percent_dist*100,2)/100;

    field_names = string(fieldnames(maps_total));
    
    standard_size = [500,1000];
    percent_map = single(zeros(standard_size));
    
    for n=1:numel(field_names)
    
        map_tmp = maps_total.(field_names(n)).map;
        map_factor_tmp = maps_total.(field_names(n)).factor;
    
        map_tmp = map_tmp-min(map_tmp(:));
    
        map_norm = map_tmp./max(map_tmp(:));
        
        ind_f = find(map_factor_tmp==factor_list);
        percent_limit = percent_dist(ind_f);
    
        percent_region = single(map_norm > percent_limit);
        percent_map = percent_map + percent_region;
    end
    
    full_criteria_map = percent_map == single(numel(field_names)-1);

    pixels_filled = sum(full_criteria_map(:));
    
    regions = bwconncomp( imgaussfilt(single(full_criteria_map),2)>0 ).NumObjects;
    
    output_vector = percent_dist.'*100;

    %fitness = mean( percent_dist(1:6) ) + mean( percent_dist(7:8) )*0.3; %ignoring gradient
    fitness = mean(percent_dist([1:7,8:9])); 
    if percent_dist(9) < 0.75 || percent_dist(1) < 0.75
        fitness = 0;
    end

    max_regions_desired = 20;
    min_regions_desired = 5;
    min_pixels_filled = 10;

    if regions > max_regions_desired
        fitness = fitness * exp( -(regions-max_regions_desired)/100 );
    end
    
    if regions < min_regions_desired
        fitness = 0;
    end

    if pixels_filled < min_pixels_filled
        fitness = 0;
    end

    region_map = full_criteria_map;

    fprintf("-")
end