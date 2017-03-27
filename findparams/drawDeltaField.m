pic_id = 0;
base_folder = '../1_single_bubble_wdlap/';
saved_folder = './saved/';
for pic_id = 100
    rhoField = load([base_folder 'rho' int2str(pic_id)]);
    vxField = load([base_folder 'vx' int2str(pic_id)]);
    %vxWatch = vxField;
    vyField = load([base_folder 'vy' int2str(pic_id)]);

    wdname = [base_folder 'wd' int2str(pic_id)];
    if (exist(wdname, 'file') == 2)
        wdfield = load(wdname);
        figure;
        gca = pcolor(wdfield);
        title('Wd');
        set(gca, 'LineStyle','none');
        colorbar;
    end

    lapname = [base_folder 'la' int2str(pic_id)];
    if (exist(lapname, 'file') == 2)
        lapfield = load(lapname);
        figure;
        gca = pcolor(lapfield);
        title('Lap');
        set(gca, 'LineStyle','none');
        colorbar;
    end
    
    % Interpolation.
    for i = 1:(size(vxField, 2)-1)
        vxField(:,i) = (vxField(:,i) + vxField(:,i+1)) ./ 2;
    end
    vxField = vxField(1:end, 1:end-1);
    for i = 1:(size(vyField, 1)-1)
        vyField(i,:) = (vyField(i,:) + vyField(i+1,:)) ./ 2;
    end
    vyField = vyField(1:end-1, 1:end);

    figure;
    quiver(vxField, vyField);
    title('速度场');
    saveas(gcf,[saved_folder 'vel' int2str(pic_id)],'bmp');

    figure;
    gca = pcolor(rhoField);
    caxis([0,0.7]);
    set(gca, 'LineStyle','none');
    %quiver(vxWatch, zeros(size(vxWatch)));
    %quiver(vxField, vyField);
    title('密度场');
    colorbar;
    saveas(gcf,[saved_folder 'pic' int2str(pic_id)],'bmp');

    % close all;
end
