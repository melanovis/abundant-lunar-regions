format compact
clear
clc
clf reset

% ----------

load("craters_formatted.mat")

mat_names = [
"yang_results.mat"
"diaz.mat"
"lawerence.mat"
"jaanga.mat"
"siegler.mat"
];

maps_total = struct();
ind_t = 1;
for n=1:numel(mat_names)
    
    load(mat_names(n))
    field_names = string(fieldnames(map_struct));

    for m=1:numel(field_names)
        maps_total = setfield(maps_total,"index_"+string(ind_t),"map",single(map_struct.(field_names(m)).map));
        maps_total = setfield(maps_total,"index_"+string(ind_t),"label",map_struct.(field_names(m)).label);
        maps_total = setfield(maps_total,"index_"+string(ind_t),"factor",map_struct.(field_names(m)).factor);
        ind_t = ind_t+1;

        map_tmp = single(map_struct.(field_names(m)).map);
        label_tmp = map_struct.(field_names(m)).label;

        scatter(nan,nan,"w","filled")
        hold on
        grid on
        axis tight equal

        xlim([-180,180])
        ylim([-90,90])

        if n==4
            map_tmp = 1-map_tmp;
        end

        surfacemap = flipud(map_tmp);
        wm = imagesc(surfacemap,'xdata',[-180,180],'ydata',[-90 90]);

        cmap = interp1([0,0.2,0.4,0.6,0.8,1], [[0 0 0]; [0.259 0.039 0.408]; [0.584 0.149 0.404]; [0.867 0.318 0.227]; [0.98 0.647 0.039]; [0.98 1 0.643]], linspace(0, 1, 1e3));
        colormap(cmap)

        h = colorbar;
        h.Color=[1,1,1];
        clim([min(map_tmp(:)),max(map_tmp(:))])

        ax = gca;
        ax.FontSize = 20;

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

        for n=1:height(crater_labels)
            if crater_latlong(n,2) < 170 && crater_latlong(n,1) < 80 && crater_diameter(n) > 100
                scatter(crater_latlong(n,2),crater_latlong(n,1),"rx")
                text(crater_latlong(n,2),crater_latlong(n,1)," "+string(crater_labels(n)),fontsize = 8,Color=[1,1,1])
            end
        end

        legend(label_tmp, 'TextColor', 'w', FontSize=15, location = "northwest", Interpreter="latex")
        legend boxoff
        
        set(findall(gcf,'-property','FontSize'), 'FontName', 'Times')

        drawnow()
        hold off

        im_raw = getframe(gcf);
        imwrite(im_raw.cdata, "map_"+string(ind_t-1)+".png");


        %pause
    end
end
save("maps_total.mat","maps_total")