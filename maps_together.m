format compact
clear
clc
clf reset

% ----------

load("craters_formatted.mat")
load("maps_total.mat")

field_names = string(fieldnames(maps_total));
for n=1:numel(field_names)
    maps_total.(field_names(n)).factor;
end

ref_list = [
"comprehensive mapping of lunar surface "+newline+"chemistry by adding chang'e-5 samples with "+newline+"deep learning, Yang et.al 2023"
"global lunar FeO mapping via wavelet "+newline+"autoencoder feature learning from M3 "+newline+"hyperspectral data, Fernandez-Diaz et.al 2025"
"global hydrogen abundances on the lunar "+newline+"surface, Lawerence et.al 2022"
"lunar heat flow: global predictions and "+newline+"reduced heat flux, Siegler et.al 2022"
"jaanga moon heightmaps, Jaanga 2026"
];

explenation = [
"This plot uses nonlinear optimization to compare " + newline + ...
"overlapping regions of lunar abundance maps from " + newline + ...
"different authors to find the regions with the greatest " + newline + ...
"simultaneous abundances which would potentially be " + newline + ...
"of most interest to industrial sites. The labelled regions " + newline + ...
"(which only span a roughly 5.45 by 5.45km region " + newline + ...
"each) possess the abundance percentiles in the minerals " + newline + ...
"shown below. Eg Th being in the 74th percentile " + newline + ...
"means it has better Th abundance than 74% of the " + newline + ...
"rest of the lunar surface."
];
explenation = "How does this plot work?"+ newline + explenation;
explenation = explenation + newline + "Please note this map is approximate and should not "+newline+"be taken too seriously, precise abundances "+newline+"are disputed due to reliance on remote sensing.";


factor_percentile = [
"TiO2", "99"
"MgO", "74.4"
"FeO", "38.31"
"CaO", "35.41"
"Al2O3", "45.19"
"SiO2", "75.83"
"H2", "67.88"
"Th", "74.26"
"gradient", "88.31"
];

region_labels = [
"Kepler abundantia prima"
"Kepler abundantia secunda"
"Kepler abundantia tertia"
"T Mayer abundantia prima"
"T Mayer abundantia secunda"
"Plato abundantia"
"Fra Mauro abundantia"
"Lalande abundantia"
];

label_offsets = [
3,-4
8,15
15,0
10,5
10,0
5,-5
7,2
3,-2.5
];

lunar_map = imread("lunar_map.png");

standard_size = [500,1000];
lunar_map = imresize(lunar_map,standard_size);

percent_map = single(zeros(standard_size));

for n=1:numel(field_names)

    map_tmp = maps_total.(field_names(n)).map;
    map_org = map_tmp;
    map_factor_tmp = maps_total.(field_names(n)).factor;
    
    map_tmp = map_tmp-min(map_tmp(:));

    map_norm = map_tmp./max(map_tmp(:));

    ind_f = find(map_factor_tmp==factor_percentile(:,1));
    percent_limit = single(str2double(factor_percentile(ind_f,2)))/100;

    map_bounds(ind_f,:) = [min(map_org(:)),max(map_org(:))];

    percent_region = single(map_norm > percent_limit);
    percent_map = percent_map + percent_region;
end

full_criteria_map = percent_map == single(numel(field_names)-1);

cell_obj = bwconncomp(full_criteria_map);
for n=1:numel(cell_obj.PixelIdxList)
    cell_tmp = cell_obj.PixelIdxList(n);
    cell_list = single(cell_tmp{:,:});
    row = [];
    col = [];
    for m=1:numel(cell_list)
        [row(m),col(m)] = ind2sub(standard_size, cell_list(m));
    end
    abundance_regions(n,:) = single( [mean(col),standard_size(1) - mean(row)+1, numel(cell_list)] );
end

for n=1:height(abundance_regions)
    %fprintf("-------\n")
    lat_spec = interp1([1,height(lunar_map)],[-90,90],abundance_regions(n,2));
    long_spec = interp1([1,width(lunar_map)],[-180,180],abundance_regions(n,1));
    % n
    % long_spec
    % lat_spec
    abundance_regions(n,1:2) = [long_spec,lat_spec];
end


scatter(nan,nan,"k")
hold on
for n=1:2
    scatter(nan,nan,"k")
end
grid on
axis tight equal

xlim([-180,180])
ylim([-90,90])

surfacemap = flipud(lunar_map);
wm = image(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

set(gca,'ydir','normal')
uistack(wm,'down')

ax = gca;
ax.FontSize = 20;

for n=1:height(abundance_regions)
    plot([0,label_offsets(n,1)]+abundance_regions(n,1) , [0,label_offsets(n,2)]+abundance_regions(n,2) , "w")
end

scatter(abundance_regions(:,1), abundance_regions(:,2), 15 ,"filled",markerfacecolor=[255,85,0]./255,MarkerEdgeColor=[0,0,0])
text(abundance_regions(:,1)+label_offsets(:,1), abundance_regions(:,2)+label_offsets(:,2), " " + string(region_labels) + " (" + string(round(abundance_regions(:,1),1))+"°, "+string(round(abundance_regions(:,2),1)) +"°)" ,fontsize = 11, Color=[1,1,1], Interpreter="latex")

for n=1:height(crater_labels)
    if crater_latlong(n,2) < 170 && crater_diameter(n) > 120
        scatter(crater_latlong(n,2),crater_latlong(n,1),"rx")
        text(crater_latlong(n,2),crater_latlong(n,1)," "+string(crater_labels(n)),fontsize = 7,Color=[1,1,1])
    end
end

%scatter(-22.95,5.44,"r","filled")

legend_strings = [
"$TiO_2$", "wrt$\%$"
"$MgO$", "wrt$\%$"
"$FeO$", "wrt$\%$"
"$CaO$", "wrt$\%$"
"$Al_2O_3$", "wrt$\%$"
"$SiO_2$", "wrt$\%$"
"$H_2$", "ppm"
"$Th$", "ppm"
"Terrain $\nabla$ (flatness)", "$\%$"
];
for n=1:height(legend_strings)
    legend_strings(n) = legend_strings(n,1) + " - " + factor_percentile(n,2) + "$\%$";
end
for n=1:height(legend_strings)-1
    tolerance = round(abs(map_bounds(n,2)-map_bounds(n,1))*0.12,2);
    legend_strings(n) = legend_strings(n,1) + " (" + round( interp1([0,100],[map_bounds(n,:)], str2double(factor_percentile(n,2)) ), 2) + "$\pm$" + tolerance + " " + legend_strings(n,2)+")";
end

legend_str = "";
for n=1:height(legend_strings)
    legend_str = legend_str + "$\cdot$ " + legend_strings(n) + newline;
end

legend_str = newline + " " + newline + "Factor percentiles ($\pm$12$\%$ bounds)" + newline + legend_str;

ref_str = newline + ""+ newline + "Refrences:" + newline;
for n=1:height(ref_list)
    ref_str = ref_str + "$\cdot$ " + ref_list(n) + newline;
end

legend_obj = legend([explenation],[legend_str],[ref_str],"TextColor","w", Interpreter="latex", FontSize=13.2,location="northeastoutside",Box="off");

title("Prime lunar real estate for industrialization", Interpreter="latex", FontSize=20)
ax = gca;
ax.TitleHorizontalAlignment = 'left';

set(gcf, 'Color', [0,0,0])
set(ax, 'Color', [0,0,0])

ax.XColor = 'w';
ax.YColor = 'w';
ax.ZColor = 'w';
ax.Title.Color  = 'w';
ax.XLabel.Color = 'w';
ax.YLabel.Color = 'w';

ax.GridColor = [1 1 1];
ax.MinorGridColor = [1 1 1];
ax.GridAlpha = 0.2;

set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

im_raw = getframe(gcf);
imwrite(im_raw.cdata, "map.png");
